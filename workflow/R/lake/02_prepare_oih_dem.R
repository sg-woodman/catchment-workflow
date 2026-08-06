# 02_prepare_oih_dem.R
# ---------------------------------------------------------------------------
# Caches the OIH Enforced DEM and derives hydrological products from the
# OIH Enhanced Flow Direction for the project study area.
#
# Unlike the stream workflow (R/stream/02–04), this module does NOT run
# WhiteboxTools depression breaching — the OIH Enhanced Flow Direction is
# already hydrologically conditioned. The D8 pointer is recoded from OIH
# encoding to WhiteboxTools encoding, then flow accumulation and stream
# extraction are run using WhiteboxTools.
#
# OIH → WhiteboxTools D8 encoding reclassification (one-step rotation):
#   OIH 128 (NE) → WBT 1     OIH 1  (E)  → WBT 2
#   OIH 2   (SE) → WBT 4     OIH 4  (S)  → WBT 8
#   OIH 8   (SW) → WBT 16    OIH 16 (W)  → WBT 32
#   OIH 32  (NW) → WBT 64    OIH 64 (N)  → WBT 128
#
# All outputs use file-existence caching — steps are skipped if the output
# already exists in cache_dir.
#
# Outputs (written to cache_dir/):
#   dem.tif        : OIH Enforced DEM (copied locally for consistent pathing)
#   d8_pntr.tif    : D8 flow pointer (recoded to WhiteboxTools encoding)
#   flow_accum.tif : D8 flow accumulation (cells)
#   streams.tif    : Binary stream network (1 = stream cell)
#
# Dependencies: terra, whitebox, fs, cli (via utils.R)
# ---------------------------------------------------------------------------

#' Prepare OIH hydrological products for a lake project
#'
#' Caches the OIH DEM and derives flow pointer, flow accumulation, and stream
#' network. Skips any step whose output already exists in cache_dir.
#'
#' @param cache_dir          Character. Project cache directory
#' @param oih_dem_path       Character. Path to OIH Enforced DEM (.tif)
#' @param oih_flow_dir_path  Character. Path to OIH Enhanced Flow Direction (.tif)
#' @param stream_threshold   Integer. Flow accumulation threshold for stream
#'                           extraction (cells). Default 100.
#'
#' @return invisibly TRUE. Called for side effects.
prepare_oih_products <- function(
  cache_dir,
  oih_dem_path,
  oih_flow_dir_path,
  stream_threshold = 100
) {
  ensure_dir(cache_dir)

  dem_path    <- fs::path(cache_dir, "dem.tif")
  pntr_path   <- fs::path(cache_dir, "d8_pntr.tif")
  accum_path  <- fs::path(cache_dir, "flow_accum.tif")
  stream_path <- fs::path(cache_dir, "streams.tif")

  if (!fs::file_exists(oih_dem_path)) {
    cw_abort(glue::glue("OIH DEM not found at: {oih_dem_path}"))
  }
  if (!fs::file_exists(oih_flow_dir_path)) {
    cw_abort(glue::glue("OIH Enhanced Flow Direction not found at: {oih_flow_dir_path}"))
  }

  # -- Step 1: Cache OIH Enforced DEM -----------------------------------------
  if (cache_exists(dem_path)) {
    cw_inform("dem.tif found in cache — skipping")
  } else {
    cw_inform("Caching OIH Enforced DEM...")
    dem <- terra::rast(oih_dem_path)
    terra::writeRaster(
      dem, dem_path,
      overwrite = TRUE,
      gdal = c("COMPRESS=LZW", "BIGTIFF=IF_SAFER")
    )
    cw_inform(glue::glue(
      "dem.tif cached: {terra::nrow(dem)} x {terra::ncol(dem)}, ",
      "res={round(terra::res(dem)[1])} m, ",
      "CRS=EPSG:{terra::crs(dem, describe = TRUE)$code}"
    ))
  }

  # -- Step 2: Recode OIH Enhanced Flow Direction to WhiteboxTools encoding ----
  if (cache_exists(pntr_path)) {
    cw_inform("d8_pntr.tif found in cache — skipping")
  } else {
    cw_inform("Recoding OIH Enhanced Flow Direction to WhiteboxTools encoding...")

    oih_recode_matrix <- matrix(
      c(128, 1, 1, 2, 2, 4, 4, 8, 8, 16, 16, 32, 32, 64, 64, 128),
      ncol = 2, byrow = TRUE
    )

    oih_pntr <- terra::rast(oih_flow_dir_path)
    pntr_recoded <- terra::classify(oih_pntr, oih_recode_matrix, others = NA)

    terra::writeRaster(
      pntr_recoded, pntr_path,
      overwrite = TRUE,
      datatype = "INT1U",
      gdal = c("COMPRESS=LZW")
    )

    vals <- sort(unique(terra::values(pntr_recoded), na.rm = TRUE))
    cw_inform(glue::glue(
      "d8_pntr.tif cached — unique values: {paste(vals, collapse = ', ')}"
    ))
  }

  # -- Step 3: Flow accumulation -----------------------------------------------
  if (cache_exists(accum_path)) {
    cw_inform("flow_accum.tif found in cache — skipping")
  } else {
    cw_inform("Computing D8 flow accumulation...")
    whitebox::wbt_d8_flow_accumulation(
      input    = fs::path_abs(pntr_path),
      output   = fs::path_abs(accum_path),
      out_type = "cells",
      pntr     = TRUE
    )
    cw_inform("flow_accum.tif cached")
  }

  # -- Step 4: Stream extraction -----------------------------------------------
  if (cache_exists(stream_path)) {
    cw_inform("streams.tif found in cache — skipping")
  } else {
    cw_inform(glue::glue(
      "Extracting stream network (threshold = {stream_threshold} cells)..."
    ))
    whitebox::wbt_extract_streams(
      flow_accum = fs::path_abs(accum_path),
      output     = fs::path_abs(stream_path),
      threshold  = stream_threshold
    )
    cw_inform("streams.tif cached")
  }

  cw_inform("OIH product preparation complete.")
  invisible(TRUE)
}

#' Load hydrological rasters from a lake project cache directory
#'
#' Returns a named list of SpatRasters suitable for use in reclip_outputs().
#' Rasters are loaded lazily (not read into memory).
#'
#' Names in the returned list follow the stream workflow convention so that
#' reclip_outputs() and downstream modules can handle both methods uniformly:
#'   dem           <- dem.tif
#'   flow_pointer  <- d8_pntr.tif
#'   flow_accum    <- flow_accum.tif
#'   streams       <- streams.tif
#'
#' @param cache_dir Character. Project cache directory
#' @return Named list of SpatRasters
load_lake_cache_rasters <- function(cache_dir) {
  name_map <- c(
    dem          = "dem.tif",
    flow_pointer = "d8_pntr.tif",
    flow_accum   = "flow_accum.tif",
    streams      = "streams.tif"
  )

  rasters <- purrr::imap(name_map, function(filename, raster_name) {
    path <- fs::path(cache_dir, filename)
    if (!cache_exists(path)) {
      cw_warn(glue::glue(
        "Lake cache raster not found: {path} — run prepare_oih_products() first"
      ))
      return(NULL)
    }
    terra::rast(path)
  })

  rasters[!vapply(rasters, is.null, logical(1))]
}
