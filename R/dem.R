# =============================================================================
# dem.R — DEM loading and hydrological conditioning
# =============================================================================
# Each function writes one product to a named path inside cache_dir and
# returns that path. targets tracks the path + content hash via tar_file(),
# so expensive WBT steps are never re-run unless their inputs change.
#
# cache_dir is keyed on wscssda_key (the sorted set of NHN work unit codes
# intersecting the AOI). All sites sharing the same work units share one
# cache_dir and therefore one set of DEM products — including sites that
# span multiple HydroBasins polygons.
#
# Dependency chain (enforced by function signatures):
#
#   make_cache_dir(wscssda_key)
#       └─► load_mrdem(aoi, cache_dir)            → dem_raw.tif
#               └─► burn_streams(dem_file,         → dem_burned.tif
#                                streams, cache_dir)
#                       └─► breach_dem(dem_file,   → dem_breached.tif
#                                      cache_dir)
#                               └─► flow_direction(dem_file,   → fdr.tif
#                                                  cache_dir)
#                                       └─► flow_accumulation(fdr_file, → fac.tif
#                                                              cache_dir)
#
# Intermediate WBT bridge files (e.g. writing a SpatRaster to disk before
# passing to a WBT function) still use tempfile — these are transient and
# not worth persisting.


#' Build and return the cache directory path for a wscssda_key
#'
#' Creates the directory if it does not exist. Called as the first target in
#' the DEM branch so all downstream targets have a stable path to write to.
#'
#' @param wscssda_key Sorted, underscore-joined WSCSSDA codes (from build_aoi())
#' @return Character path to the cache directory
make_cache_dir <- function(wscssda_key) {
  path <- file.path(CACHE_DIR, wscssda_key)
  fs::dir_create(path)
  path
}


#' Crop the MRDEM Virtual Raster to an AOI and write to cache
#'
#' @param aoi       sf polygon (any CRS)
#' @param cache_dir Path returned by make_cache_dir()
#' @return Path to dem_raw.tif
load_mrdem <- function(aoi, cache_dir) {
  out <- file.path(cache_dir, "dem_raw.tif")
  if (file.exists(out)) {
    message("  Using cached dem_raw.tif")
    return(out)
  }

  if (!file.exists(MRDEM_PATH))
    stop("MRDEM VRT not found: ", MRDEM_PATH)

  mrdem    <- rast(MRDEM_PATH)
  aoi_proj <- project(vect(aoi), crs(mrdem))

  if (!relate(ext(aoi_proj), ext(mrdem), relation = "intersects"))
    stop("AOI does not overlap MRDEM extent. Check pour point coordinates.")

  message("  Cropping MRDEM to AOI...")
  dem <- crop(mrdem, aoi_proj)
  message("  DEM: ", nrow(dem), " rows x ", ncol(dem), " cols | CRS: ",
          crs(dem, describe = TRUE)$code)

  writeRaster(dem, out, overwrite = FALSE)
  out
}


#' Burn NHN stream lines into a DEM and write to cache
#'
#' If streams is NULL (burn_streams = FALSE for all sites sharing this
#' cache_dir), copies dem_raw.tif to dem_burned.tif unchanged so that
#' downstream targets can always depend on dem_burned.tif existing.
#'
#' After burning, cells below the valid DEM minimum are clamped to prevent
#' WBT's nodata sentinel (~-32768) from leaking into the breaching step as
#' extreme negative elevations.
#'
#' @param dem_file  Path to dem_raw.tif (from load_mrdem())
#' @param streams   sf MULTILINESTRING or NULL
#' @param cache_dir Path returned by make_cache_dir()
#' @return Path to dem_burned.tif
burn_streams <- function(dem_file, streams, cache_dir) {
  out <- file.path(cache_dir, "dem_burned.tif")
  if (file.exists(out)) {
    message("  Using cached dem_burned.tif")
    return(out)
  }

  if (is.null(streams)) {
    message("  No streams — copying dem_raw.tif to dem_burned.tif")
    file.copy(dem_file, out)
    return(out)
  }

  dem <- rast(dem_file)

  # Streams written to a tempfile — transient WBT bridge, not worth caching.
  # Cast to MULTILINESTRING to avoid WBT ShapeType panic on mixed geometries.
  streams_tmp <- tempfile(fileext = ".shp")
  st_write(
    st_cast(
      st_transform(streams, paste0("EPSG:", crs(dem, describe = TRUE)$code)),
      "MULTILINESTRING"
    ),
    streams_tmp,
    quiet = TRUE
  )

  message("  Burning streams into DEM (", nrow(streams), " features)...")
  wbt_fill_burn(dem = dem_file, streams = streams_tmp, output = out)

  # Clamp sentinel values — WBT writes its nodata fill value into cells
  # outside the stream vector extent when the DEM has no nodata set
  burned    <- rast(out)
  valid_min <- min(values(dem), na.rm = TRUE)
  n_bad     <- sum(values(burned) < valid_min, na.rm = TRUE)
  if (n_bad > 0) {
    message("  Clamping ", n_bad, " sub-minimum cells in burned DEM")
    burned[burned < valid_min] <- valid_min
    writeRaster(burned, out, overwrite = TRUE)
  }

  out
}


#' Breach depressions in a DEM and write to cache
#'
#' @param dem_file  Path to dem_burned.tif (from burn_streams())
#' @param cache_dir Path returned by make_cache_dir()
#' @return Path to dem_breached.tif
breach_dem <- function(dem_file, cache_dir) {
  out <- file.path(cache_dir, "dem_breached.tif")
  if (file.exists(out)) {
    message("  Using cached dem_breached.tif")
    return(out)
  }

  message("  Breaching depressions...")
  wbt_breach_depressions_least_cost(
    dem    = dem_file,
    output = out,
    dist   = WBT_BREACH_DIST,
    fill   = FALSE
  )
  out
}


#' Compute D8 flow direction and write to cache
#'
#' @param dem_file  Path to dem_breached.tif (from breach_dem())
#' @param cache_dir Path returned by make_cache_dir()
#' @return Path to fdr.tif
flow_direction <- function(dem_file, cache_dir) {
  out <- file.path(cache_dir, "fdr.tif")
  if (file.exists(out)) {
    message("  Using cached fdr.tif")
    return(out)
  }

  message("  Computing D8 flow direction...")
  wbt_d8_pointer(dem = dem_file, output = out)
  out
}


#' Compute D8 flow accumulation and write to cache
#'
#' Takes the flow direction raster as input — not the DEM. The signature
#' enforces this: passing dem_breached here is a type error.
#'
#' @param fdr_file  Path to fdr.tif (from flow_direction())
#' @param cache_dir Path returned by make_cache_dir()
#' @return Path to fac.tif
flow_accumulation <- function(fdr_file, cache_dir) {
  out <- file.path(cache_dir, "fac.tif")
  if (file.exists(out)) {
    message("  Using cached fac.tif")
    return(out)
  }

  message("  Computing flow accumulation...")
  wbt_d8_flow_accumulation(input = fdr_file, output = out, out_type = "cells")
  out
}
