# =============================================================================
# nhn.R — National Hydro Network download and processing
# =============================================================================
# Both streams and lakes are read from the same GDB work units. GDBs are
# downloaded once to NHN_RAW_DIR and reused across all projects permanently.
#
# Cache key: the sorted set of WSCSSDA codes intersecting the AOI (wscssda_key
# from build_aoi()). This ensures that sites spanning multiple HydroBasins
# polygons but covered by the same NHN work units share a single download.


# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

#' List NHN GDB zip filenames for a WSCSSDA code via FTP directory listing
#'
#' @param wscssda 4-character code, lowercase (e.g. "10ca")
#' @param ftp_base Base FTP URL
#' @return Character vector of matching zip filenames
nhn_list_ftp <- function(wscssda, ftp_base) {
  wscmda  <- substr(wscssda, 1, 2)
  dir_url <- paste(ftp_base, wscmda, "", sep = "/")

  raw <- tryCatch(
    RCurl::getURL(dir_url, ftp.use.epsv = FALSE, dirlistonly = TRUE,
                  .encoding = "UTF-8"),
    error = function(e) {
      warning("FTP listing failed for WSCMDA '", wscmda, "': ", e$message)
      ""
    }
  )

  all_files <- strsplit(raw, "\r?\n")[[1]]
  pattern   <- paste0("^nhn_rhn_", wscssda, ".*_gdb_en\\.zip$")
  matched   <- all_files[grepl(pattern, all_files, ignore.case = TRUE)]

  if (length(matched) == 0)
    warning("No NHN GDB files found on FTP for WSCSSDA '", wscssda, "'")

  matched
}


#' Ensure a NHN GDB is present locally, downloading from FTP if needed
#'
#' @param wscssda  4-character WSCSSDA code
#' @param ftp_base Base FTP URL
#' @return Character vector of local .gdb directory paths for this WSCSSDA
nhn_ensure_gdbs <- function(wscssda, ftp_base) {
  wscmda <- substr(wscssda, 1, 2)

  # Check for already-extracted GDB directories before hitting FTP
  all_dirs  <- list.dirs(NHN_RAW_DIR, full.names = FALSE, recursive = FALSE)
  cached    <- all_dirs[grepl(paste0("^nhn_rhn_", wscssda, ".*_gdb_en$"),
                              all_dirs, ignore.case = TRUE)]

  zip_stems <- if (length(cached) > 0) {
    message("  Using cached GDB(s) for ", wscssda, ": ",
            paste(cached, collapse = ", "))
    cached
  } else {
    zips <- nhn_list_ftp(wscssda, ftp_base)
    if (length(zips) == 0) return(character(0))
    sub("\\.zip$", "", zips)
  }

  # Download and extract any stems not yet on disk
  for (stem in zip_stems) {
    gdb_dir <- file.path(NHN_RAW_DIR, stem)
    if (dir.exists(gdb_dir)) next

    zip_file <- paste0(stem, ".zip")
    zip_dest <- file.path(NHN_RAW_DIR, zip_file)
    url      <- paste(ftp_base, wscmda, zip_file, sep = "/")

    message("  Downloading ", zip_file, " (may take several minutes)...")
    success <- nhn_download_zip(url, zip_dest)
    if (!success) next

    unzip(zip_dest, exdir = gdb_dir)
    file.remove(zip_dest)
  }

  # Return paths to .gdb directories (ESRI GDB is a directory, not a file)
  lapply(zip_stems, function(stem) {
    gdb_dir <- file.path(NHN_RAW_DIR, stem)
    gdbs    <- list.dirs(gdb_dir, full.names = TRUE, recursive = TRUE)
    gdbs    <- gdbs[grepl("\\.gdb$", gdbs, ignore.case = TRUE)]
    if (length(gdbs) == 0) {
      warning("No .gdb directory found in: ", gdb_dir)
      return(NULL)
    }
    gdbs[1]
  }) |>
    Filter(Negate(is.null), x = _) |>
    unlist()
}


#' Download a zip file with retries
#'
#' @param url      Remote URL
#' @param destfile Local destination path
#' @param retries  Number of attempts (default 3)
#' @return Logical — TRUE on success
nhn_download_zip <- function(url, destfile, retries = 3L) {
  orig_timeout <- getOption("timeout")
  on.exit(options(timeout = orig_timeout))
  options(timeout = 600L)

  for (attempt in seq_len(retries)) {
    if (attempt > 1) {
      message("  Retry ", attempt, "/", retries)
      Sys.sleep(5)
    }

    result <- tryCatch(
      {
        download.file(url, destfile, mode = "wb", quiet = FALSE)
        0L
      },
      warning = function(w) {
        if (grepl("downloaded length", conditionMessage(w))) {
          if (file.exists(destfile)) file.remove(destfile)
          return(1L)
        }
        invokeRestart("muffleWarning")
        0L
      },
      error = function(e) {
        warning("Download error: ", e$message)
        if (file.exists(destfile)) file.remove(destfile)
        1L
      }
    )

    ok <- result == 0L && file.exists(destfile) && file.size(destfile) > 0
    if (!ok) next

    # Verify zip integrity before returning success
    zip_ok <- tryCatch(nrow(unzip(destfile, list = TRUE)) > 0, error = function(e) FALSE)
    if (!zip_ok) {
      warning("Downloaded zip appears corrupt: ", basename(destfile))
      file.remove(destfile)
      next
    }

    return(TRUE)
  }

  warning("Failed to download after ", retries, " attempts: ", basename(destfile),
          "\nTry downloading manually and placing in: ", NHN_RAW_DIR)
  FALSE
}


#' Read one layer from a NHN GDB, clip to AOI, and return geometry-only sf
#'
#' @param gdb_path Path to the .gdb directory
#' @param layer    Layer name to read
#' @param aoi      sf polygon to clip to
#' @param keep_types Character vector of geometry types to keep after clipping
#' @return sf object or NULL if no features remain
nhn_read_layer <- function(gdb_path, layer, aoi, keep_types) {
  tryCatch(
    {
      raw <- st_read(gdb_path, layer = layer, quiet = TRUE) |>
        st_zm() |>
        st_transform(st_crs(aoi)) |>
        st_intersection(st_union(aoi))

      raw <- raw[st_geometry_type(raw) %in% keep_types, ]
      if (nrow(raw) == 0) return(NULL)

      # Use st_geometry() — NHN GDBs name their column "SHAPE" not "geometry"
      st_sf(geometry = st_geometry(raw))
    },
    error = function(e) {
      warning("Could not read layer '", layer, "' from ", basename(gdb_path),
              ": ", e$message)
      NULL
    }
  )
}


#' Merge a list of sf objects, casting all to a uniform geometry type
#'
#' @param sf_list   List of sf objects
#' @param cast_type Target geometry type (e.g. "MULTILINESTRING")
#' @return Single merged sf object
nhn_merge <- function(sf_list, cast_type) {
  sf_list <- Filter(Negate(is.null), sf_list)
  if (length(sf_list) == 0) return(NULL)

  do.call(rbind, lapply(sf_list, function(x) {
    tryCatch(
      st_cast(x, cast_type),
      error = function(e) {
        warning("Could not cast to ", cast_type, ": ", e$message)
        x
      }
    )
  }))
}


# -----------------------------------------------------------------------------
# Public functions
# -----------------------------------------------------------------------------

#' Download and merge NHN stream lines for an AOI
#'
#' Reads NHN_HN_PrimaryDirectedNLFlow_1 (or NHN_HN_NLFLOW_1 as fallback)
#' from all GDB work units intersecting the AOI. GDBs are downloaded to
#' NHN_RAW_DIR if not already present.
#'
#' @param aoi sf polygon defining the area of interest (EPSG:3979)
#' @return sf MULTILINESTRING object clipped to AOI, or NULL if none found
download_nhn_streams <- function(aoi) {
  ftp_base <- "ftp://ftp.maps.canada.ca/pub/nrcan_rncan/vector/geobase_nhn_rhn/gdb_en"
  wscssda_codes <- nhn_wscssda_for_aoi(aoi)
  fs::dir_create(NHN_RAW_DIR)

  stream_list <- list()

  for (wscssda in wscssda_codes) {
    gdb_paths <- nhn_ensure_gdbs(wscssda, ftp_base)

    for (gdb_path in gdb_paths) {
      layers <- tryCatch(st_layers(gdb_path)$name, error = function(e) NULL)
      if (is.null(layers)) next

      # Prefer primary directed network; fall back to full network
      candidates <- c("NHN_HN_PrimaryDirectedNLFlow_1", "NHN_HN_NLFLOW_1")
      layer      <- intersect(candidates, layers)[1]

      if (is.na(layer)) {
        warning("No flow layer found in ", basename(gdb_path),
                "\nAvailable: ", paste(layers, collapse = ", "))
        next
      }

      result <- nhn_read_layer(gdb_path, layer, aoi,
                               keep_types = c("LINESTRING", "MULTILINESTRING"))
      if (!is.null(result)) {
        message("  ", nrow(result), " stream features from ", basename(gdb_path))
        stream_list[[gdb_path]] <- result
      }
    }
  }

  merged <- nhn_merge(stream_list, "MULTILINESTRING")
  if (is.null(merged))
    warning("No stream features found for this AOI — proceeding without burn-in")
  merged
}


#' Download and merge NHN lake polygons for an AOI
#'
#' Reads NHN_HD_WATERBODY_2 from all GDB work units intersecting the AOI.
#' GDBs must already be present (run download_nhn_streams() first, or ensure
#' burn_streams = TRUE for this site so GDBs are downloaded).
#'
#' @param aoi sf polygon defining the area of interest (EPSG:3979)
#' @return sf MULTIPOLYGON object clipped to AOI, or NULL if none found
download_nhn_lakes <- function(aoi) {
  wscssda_codes <- nhn_wscssda_for_aoi(aoi)
  lake_list     <- list()

  for (wscssda in wscssda_codes) {
    # Lakes reuse GDBs already downloaded for streams — no FTP call needed
    all_dirs <- list.dirs(NHN_RAW_DIR, full.names = FALSE, recursive = FALSE)
    stems    <- all_dirs[grepl(paste0("^nhn_rhn_", wscssda, ".*_gdb_en$"),
                               all_dirs, ignore.case = TRUE)]

    if (length(stems) == 0) {
      message("  No cached GDBs for ", wscssda,
              " — skipping lakes (run streams first)")
      next
    }

    for (stem in stems) {
      gdb_path <- list.dirs(file.path(NHN_RAW_DIR, stem),
                            full.names = TRUE, recursive = TRUE)
      gdb_path <- gdb_path[grepl("\\.gdb$", gdb_path, ignore.case = TRUE)][1]
      if (is.na(gdb_path)) next

      layers <- tryCatch(st_layers(gdb_path)$name, error = function(e) NULL)
      if (is.null(layers) || !"NHN_HD_WATERBODY_2" %in% layers) next

      result <- nhn_read_layer(gdb_path, "NHN_HD_WATERBODY_2", aoi,
                               keep_types = c("POLYGON", "MULTIPOLYGON"))
      if (!is.null(result)) {
        message("  ", nrow(result), " lake features from ", basename(gdb_path))
        lake_list[[gdb_path]] <- result
      }
    }
  }

  nhn_merge(lake_list, "MULTIPOLYGON")
}


#' Return WSCSSDA codes for all NHN work units intersecting an AOI
#'
#' @param aoi sf polygon (any CRS)
#' @return Character vector of lowercase WSCSSDA codes
nhn_wscssda_for_aoi <- function(aoi) {
  nhn_index  <- st_read(NHN_INDEX_PATH, quiet = TRUE) |>
    st_transform(st_crs(aoi))
  intersects <- st_intersects(nhn_index, st_union(aoi), sparse = FALSE)[, 1]
  codes      <- unique(tolower(trimws(nhn_index$WSCSSDA[intersects])))

  if (length(codes) == 0)
    stop("No NHN work units intersect the AOI.")

  message("  NHN work units: ", paste(codes, collapse = ", "))
  codes
}
