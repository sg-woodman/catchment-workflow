# 03_prepare_streams_burn.R
# ---------------------------------------------------------------------------
# Resolves config$streams_burn for the raw-dem terrain tier only (burning
# has to happen before flow direction is derived — see
# 00_resolve_config.R's tier-vs-burn validation). Called from
# 02_prepare_terrain.R::prepare_raw_dem_tier().
#
#   source = "nhn_auto"  — download_nhn() below, ported from the predecessor
#     project /Users/sam/Documents/cfs/catchment_delineation/code/functions.R
#     (confirmed working there: FTP directory listing per NHN work unit,
#     local-cache-first, retry + corrupt-zip handling, layer-name fallback
#     chain). Adapted to this repo's cw_inform/cw_warn/cw_abort/ensure_dir/
#     cache_exists conventions — logic unchanged.
#   source = "supplied"  — reads config$streams_burn$path directly (any
#     vector format sf::st_read() handles), clips to the group AOI.
#   source = "none"      — never reaches this file (00_resolve_config.R
#     only allows streams_burn to matter for the dem tier, and
#     02_prepare_terrain.R only calls resolve_streams_burn() when
#     source != "none").
#
# Requires workflow/R/stream/burn_streams.R already sourced for
# burn_streams_into_dem() (reused unmodified — it only needs flowlines +
# a dem path + an output path, no CRS assumptions beyond "same CRS as the
# dem", which is guaranteed here since both are resolved against
# config$working_crs).
#
# download_nhn() requires the RCurl package for FTP directory listing —
# only when source = "nhn_auto" is actually used; not part of this
# project's general required-package list.
#
# Dependencies: sf, fs, glue, cli (via utils.R); RCurl (nhn_auto only)
# ---------------------------------------------------------------------------

#' Resolve streams_burn for one group and burn into its DEM
#'
#' @param config          Resolved config from resolve_engine_config()
#' @param aoi             sfc polygon. Group AOI, in config$working_crs
#' @param dem_path        Character. Path to the group's dem.tif
#' @param dem_burned_path Character. Output path for the burned DEM
#' @param group_id        Character. Group identifier (log messages)
#' @return Invisibly, dem_burned_path if burning succeeded, else NULL
resolve_streams_burn <- function(config, aoi, dem_path, dem_burned_path, group_id) {
  if (!exists("burn_streams_into_dem", mode = "function")) {
    cw_abort(paste(
      "resolve_streams_burn() requires workflow/R/stream/burn_streams.R",
      "to be sourced first (defines burn_streams_into_dem())."
    ))
  }

  aoi_sf <- sf::st_as_sf(aoi)

  flowlines <- switch(
    config$streams_burn$source,
    nhn_auto = download_nhn(
      aoi            = aoi_sf,
      nhn_index_path = config$nhn_index_path,
      nhn_raw_dir    = config$nhn_raw_dir
    ),
    supplied = {
      cw_inform(glue::glue("Group '{group_id}': reading supplied burn-in streams..."))
      raw <- sf::st_read(config$streams_burn$path, quiet = TRUE) |>
        sf::st_transform(sf::st_crs(aoi_sf))
      clipped <- tryCatch(
        sf::st_intersection(raw, aoi_sf),
        error = function(e) {
          cw_warn(glue::glue(
            "Group '{group_id}': error clipping supplied streams — {e$message}"
          ))
          raw[0, ]
        }
      )
      if (nrow(clipped) == 0) NULL else clipped
    }
  )

  if (is.null(flowlines) || nrow(flowlines) == 0) {
    cw_warn(glue::glue(
      "Group '{group_id}': no burn-in streams available — dem.tif will be ",
      "used as-is for breaching."
    ))
    return(invisible(NULL))
  }

  # Persist to the group cache, mirroring stream/burn_streams.R's
  # prepare_nhn_layers() (old pipeline) — required so that
  # engine/04_delineate_site.R's clip_flowlines_to_catchment() (guarded on
  # cache_exists(flowlines_path)) can write each site's own streams.gpkg.
  # Found missing entirely: this function used to compute flowlines
  # in-memory and hand them straight to burn_streams_into_dem() without
  # ever writing them to disk, so every engine-based raw-dem project
  # silently never produced a per-site streams.gpkg — no warning,
  # since 04_delineate_site.R's guard just skips the call when the path
  # doesn't exist. dem_path's directory is the group's cache_dir (the
  # caller always passes grp_cache/dem.tif — see 02_prepare_terrain.R's
  # prepare_raw_dem_tier()). Deliberately flowlines only, not a
  # waterbodies.gpkg companion — the old pipeline's waterbodies.gpkg was
  # unused downstream by anything the engine relies on (lake-bisection
  # checking does its own on-the-fly NHN fetch, not a persisted group
  # file), and fetching it unconditionally here would add real runtime to
  # every group's terrain prep for no consumer.
  flowlines_path <- fs::path(fs::path_dir(dem_path), "flowlines.gpkg")
  if (!cache_exists(flowlines_path)) {
    sf::st_write(flowlines, flowlines_path, delete_dsn = TRUE, quiet = TRUE)
    cw_inform(glue::glue(
      "Group '{group_id}': flowlines.gpkg written ({nrow(flowlines)} feature(s))."
    ))
  }

  burn_streams_into_dem(
    flowlines       = flowlines,
    dem_path        = dem_path,
    dem_burned_path = dem_burned_path,
    group_id        = group_id
  )
}

# -- NHN auto-download (ported from catchment_delineation/code/functions.R) --

#' List NHN GDB zip filenames available on the FTP server for one WSCSSDA
#'
#' @param wscssda  4-character WSCSSDA base code, lowercase (e.g. "01aa")
#' @param ftp_base Base FTP URL for the gdb_en directory
#' @return Character vector of matching zip filenames, or character(0)
list_nhn_ftp_files <- function(wscssda, ftp_base) {
  wscmda <- substr(wscssda, 1, 2)
  dir_url <- paste(ftp_base, wscmda, "", sep = "/") # trailing slash needed

  filenames <- tryCatch(
    {
      if (!requireNamespace("RCurl", quietly = TRUE)) {
        cw_abort(paste(
          "RCurl package is required for NHN FTP directory listing.",
          "Install with: install.packages('RCurl')"
        ))
      }
      raw <- RCurl::getURL(
        dir_url,
        ftp.use.epsv = FALSE, dirlistonly = TRUE, .encoding = "UTF-8"
      )
      strsplit(raw, "\r?\n")[[1]]
    },
    error = function(e) {
      cw_warn(glue::glue("Could not list FTP directory for WSCMDA '{wscmda}': {e$message}"))
      character(0)
    }
  )

  pattern <- paste0("^nhn_rhn_", wscssda, ".*_gdb_en\\.zip$")
  matched <- filenames[grepl(pattern, filenames, ignore.case = TRUE)]

  if (length(matched) == 0) {
    cw_warn(glue::glue(
      "No NHN GDB files found on FTP for WSCSSDA '{wscssda}' ",
      "(directory: {dir_url}, pattern: {pattern})"
    ))
  }

  matched
}

#' Download NHN stream lines covering an AOI from the GDB source
#'
#' Queries the NRCan FTP server for National Hydro Network (NHN) File
#' Geodatabases, extracts the primary directed flow layer, clips to the
#' AOI, and returns merged stream lines as an sf object. Checks a local
#' cache directory first — a WSCSSDA whose GDB is already extracted in
#' nhn_raw_dir is used directly, no FTP listing needed.
#'
#' The NHN index shapefile must already exist locally (nhn_index_path).
#' Download once from:
#'   https://ftp.maps.canada.ca/pub/nrcan_rncan/vector/geobase_nhn_rhn/
#'     index/NHN_INDEX_WORKUNIT_LIMIT_2.zip
#'
#' @param aoi             sf polygon defining the area of interest (any CRS)
#' @param nhn_index_path  Full path to the NHN index shapefile
#' @param nhn_raw_dir     Directory for raw NHN GDB downloads (shared data,
#'                        kept permanently alongside other shared datasets)
#' @return sf object of stream lines clipped to AOI (CRS matches aoi),
#'         or NULL if no streams are found
download_nhn <- function(aoi, nhn_index_path, nhn_raw_dir) {
  # Use ftp:// protocol — the https:// mirror blocks directory listing
  ftp_base <- "ftp://ftp.maps.canada.ca/pub/nrcan_rncan/vector/geobase_nhn_rhn/gdb_en"

  if (!fs::file_exists(nhn_index_path)) {
    cw_abort(glue::glue(
      "NHN index shapefile not found: {nhn_index_path}\n",
      "Download once from: https://ftp.maps.canada.ca/pub/nrcan_rncan/vector/",
      "geobase_nhn_rhn/index/NHN_INDEX_WORKUNIT_LIMIT_2.zip\n",
      "Then set nhn_index_path in run_config."
    ))
  }

  cw_inform("Reading NHN work unit index...")
  nhn_index <- sf::st_read(nhn_index_path, quiet = TRUE) |>
    sf::st_transform(sf::st_crs(aoi))

  aoi_units <- nhn_index[
    sf::st_intersects(nhn_index, sf::st_union(aoi), sparse = FALSE)[, 1],
  ]

  if (nrow(aoi_units) == 0) {
    cw_warn("No NHN work units intersect the AOI. Proceeding without stream burn-in.")
    return(NULL)
  }

  wscssda_codes <- unique(tolower(trimws(aoi_units$WSCSSDA)))
  cw_inform(glue::glue("NHN WSCSSDA codes required: {paste(wscssda_codes, collapse = ', ')}"))

  ensure_dir(nhn_raw_dir)

  stream_list <- list()

  for (wscssda in wscssda_codes) {
    wscmda <- substr(wscssda, 1, 2)

    # -- Check local cache first -----------------------------------------
    cached_stems <- list.dirs(nhn_raw_dir, full.names = FALSE, recursive = FALSE)
    cached_stems <- cached_stems[grepl(
      paste0("^nhn_rhn_", wscssda, ".*_gdb_en$"), cached_stems, ignore.case = TRUE
    )]

    if (length(cached_stems) > 0) {
      cw_inform(glue::glue(
        "Using {length(cached_stems)} cached GDB(s) for WSCSSDA {wscssda}: ",
        "{paste(cached_stems, collapse = ', ')}"
      ))
      zip_files <- paste0(cached_stems, ".zip")
    } else {
      cw_inform(glue::glue("Listing FTP directory for WSCSSDA {wscssda}..."))
      zip_files <- list_nhn_ftp_files(wscssda, ftp_base)
      if (length(zip_files) == 0) {
        cw_warn(glue::glue("No files found on FTP for WSCSSDA {wscssda} — skipping."))
        next
      }
      cw_inform(glue::glue(
        "Found {length(zip_files)} GDB file(s) to download: {paste(zip_files, collapse = ', ')}"
      ))
    }

    for (zip_file in zip_files) {
      zip_stem <- sub("\\.zip$", "", zip_file)
      gdb_dir  <- file.path(nhn_raw_dir, zip_stem)
      url      <- paste(ftp_base, wscmda, zip_file, sep = "/")

      if (!dir.exists(gdb_dir)) {
        zip_dest <- file.path(nhn_raw_dir, zip_file)
        max_tries <- 3
        success <- FALSE
        orig_timeout <- getOption("timeout")
        options(timeout = 600) # 10 minutes — large files over FTP need this
        on.exit(options(timeout = orig_timeout), add = TRUE)

        for (attempt in seq_len(max_tries)) {
          if (attempt > 1) {
            cw_inform(glue::glue("Retry {attempt}/{max_tries}: {zip_file}"))
            Sys.sleep(5)
          } else {
            cw_inform(glue::glue("Downloading: {zip_file} (may take several minutes)"))
          }

          dl_result <- tryCatch(
            {
              download.file(url = url, destfile = zip_dest, mode = "wb", quiet = FALSE)
              0L
            },
            warning = function(w) {
              if (grepl("downloaded length", conditionMessage(w), fixed = TRUE)) {
                cw_warn(glue::glue("Incomplete download (truncated): {zip_file}"))
                if (file.exists(zip_dest)) file.remove(zip_dest)
                return(1L)
              }
              invokeRestart("muffleWarning")
              0L
            },
            error = function(e) {
              cw_warn(glue::glue("Download error (attempt {attempt}): {e$message}"))
              if (file.exists(zip_dest)) file.remove(zip_dest)
              1L
            }
          )

          if (dl_result == 0L && file.exists(zip_dest) && file.size(zip_dest) > 0) {
            success <- TRUE
            break
          }
        }

        options(timeout = orig_timeout)
        on.exit(NULL)

        if (!success) {
          cw_warn(glue::glue(
            "Failed to download after {max_tries} attempts: {zip_file} (url: {url}). ",
            "Try downloading manually and placing in: {nhn_raw_dir}"
          ))
          if (file.exists(zip_dest)) file.remove(zip_dest)
          next
        }

        zip_ok <- tryCatch(
          {
            nms <- unzip(zip_dest, list = TRUE)
            nrow(nms) > 0
          },
          error = function(e) FALSE
        )
        if (!zip_ok) {
          cw_warn(glue::glue("Downloaded file appears corrupt (bad zip): {zip_file} — removing."))
          file.remove(zip_dest)
          next
        }

        unzip(zip_dest, exdir = gdb_dir)
        file.remove(zip_dest)
      }

      gdb_path <- list.dirs(gdb_dir, full.names = TRUE, recursive = TRUE)
      gdb_path <- gdb_path[grepl("\\.gdb$", gdb_path, ignore.case = TRUE)][1]

      if (is.na(gdb_path)) {
        cw_warn(glue::glue("No .gdb directory found in: {gdb_dir}"))
        next
      }
      cw_inform(glue::glue("Found GDB: {basename(gdb_path)}"))

      gdb_layers <- tryCatch(
        sf::st_layers(gdb_path)$name,
        error = function(e) {
          cw_warn(glue::glue(
            "Cannot read GDB (GDAL driver incompatibility?): {basename(gdb_path)} — ",
            "{e$message}. Skipping {zip_stem}. If needed, convert manually: ",
            "ogr2ogr -f GPKG out.gpkg {gdb_path}"
          ))
          NULL
        }
      )
      if (is.null(gdb_layers)) next

      candidate_names <- c("NHN_HN_PrimaryDirectedNLFlow_1", "NHN_HN_NLFLOW_1")
      flow_layer <- intersect(candidate_names, gdb_layers)[1]

      if (is.na(flow_layer)) {
        flow_layer <- gdb_layers[grepl(
          "(?i)primarydirected|nlflow|directed.*flow|flow.*directed",
          gdb_layers, perl = TRUE
        )][1]
      }
      if (is.na(flow_layer)) {
        cw_warn(glue::glue(
          "No flow layer found in: {basename(gdb_path)} ",
          "(available: {paste(gdb_layers, collapse = ', ')})"
        ))
        next
      }

      cw_inform(glue::glue("Reading '{flow_layer}' from {zip_stem}..."))

      streams_unit <- tryCatch(
        {
          raw <- sf::st_read(gdb_path, layer = flow_layer, quiet = TRUE) |>
            sf::st_zm() |>
            sf::st_transform(sf::st_crs(aoi)) |>
            sf::st_intersection(sf::st_union(aoi))

          raw <- raw[sf::st_geometry_type(raw) %in% c("LINESTRING", "MULTILINESTRING"), ]

          if (nrow(raw) == 0) NULL else sf::st_sf(geometry = sf::st_geometry(raw))
        },
        error = function(e) {
          cw_warn(glue::glue("Failed to read layer from {zip_stem}: {e$message}"))
          NULL
        }
      )

      if (!is.null(streams_unit) && nrow(streams_unit) > 0) {
        cw_inform(glue::glue("Loaded {nrow(streams_unit)} features from {zip_stem}"))
        stream_list[[zip_stem]] <- streams_unit
      }
    }
  }

  stream_list <- Filter(Negate(is.null), stream_list)
  if (length(stream_list) == 0) {
    cw_warn("No stream features found in any NHN work unit. Proceeding without stream burn-in.")
    return(NULL)
  }

  # Cast to MULTILINESTRING before merging — mixing LINESTRING and
  # MULTILINESTRING across GDBs produces a corrupt ShapeType downstream.
  streams <- do.call(
    rbind,
    lapply(stream_list, function(x) {
      geom <- sf::st_geometry(x)
      tryCatch(
        sf::st_cast(sf::st_sf(geometry = geom), "MULTILINESTRING"),
        error = function(e) {
          cw_warn(glue::glue("Could not cast geometry to MULTILINESTRING: {e$message}"))
          sf::st_sf(geometry = geom)
        }
      )
    })
  )

  cw_inform(glue::glue("NHN streams loaded: {nrow(streams)} features total"))
  streams
}
