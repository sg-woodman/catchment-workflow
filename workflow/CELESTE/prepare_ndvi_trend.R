# prepare_ndvi_trend.R
# =============================================================================
# Prepares the NDVI trend LOI raster (Theil-Sen slope + Mann-Kendall p-value,
# 1984-2025) for the CELESTE stream hydroweight stage (run_celeste.R, Stage
# 9). Kept as its own script for the same reason as prepare_harvest_regen.R
# — it's a distinct, slow (minutes per group), one-time data-prep step, not
# part of the interactive delineation stage sequence.
#
# One raster per group, written to
# cache/CELESTE/hydroweight_loi/ndvi_trend/<group_id>.tif, 2 bands
# ("slope", "p_value") — matches sens_slope_trend()'s return shape exactly,
# so no reshaping is needed before handing it to the hydroweight stage as a
# path_lazy LOI (06_hydroweight_attributes.R already loops any multi-layer
# continuous LOI one band at a time and uses each band's own name() as the
# output column prefix — see run_loi_attributes_stream_multilayer()).
#
# Operates on the CLEANED tiles from prepare_ndvi.R's clean_ndvi_tiles()
# (0 -> NA masked-pixel fix + /10000 rescale to true NDVI), not the raw
# data/ndvi/*.tif exports directly — see that script's header for why both
# are necessary (raw 0 is Earth Engine's masked-pixel export fill value,
# not a genuine NDVI reading; treating it as a real annual observation
# would corrupt the Theil-Sen slope with false "crashes to zero" in
# cloudy/gap years instead of correctly excluding them via min_obs).
#
# Design / performance notes (see CLAUDE.md's "known quirks" section):
#
#   - sens_slope_trend() cost is ~0.5 ms per VALID pixel, essentially
#     independent of the raster's total footprint or tile count (confirmed
#     directly: COC/MOR/TUR — single, fully-valid tiles — and NIP/NBE —
#     multi-tile, mostly-NA mosaics — all landed within ~15% of the same
#     time-per-valid-cell). Theil-Sen's combn(n_years, 2) pairwise
#     comparison is what actually costs time, not raster geometry. This is
#     why this script mosaics ONLY each group's own source NDVI tiles
#     (regex-matched by group_id in the filename) rather than cropping the
#     full multi-province mosaic VRT (tiles span three widely-separated
#     provinces — Ontario: KEN/NIP/TUR; New Brunswick: NBE/COC; PEI: MOR —
#     NOT transcontinental; confirmed by reprojecting each tile's extent to
#     lon/lat) to the group's AOI bounding box — the latter was tested on
#     COC (a small group) and took 752.9s against an 80.7M-cell, 99%-NA
#     crop, vs. 248.5s for the same group's own 508K-cell tile. Comparable
#     valid-pixel count, ~3x the wall time, purely from iterating over
#     padding cells that were never going to have data — a HydroBasins-
#     derived AOI (plus buffer) is typically much larger than the tight
#     footprint Earth Engine actually exported around each group's sites.
#   - terra::app(cores > 1) was tested (NIP, cores = 8) and made this
#     dramatically SLOWER, not faster: 2896.8s vs. 306.5s single-threaded —
#     a ~9.5x regression, not just a failure to parallelize. Most likely
#     cause: terra::app()'s parallel path uses a PSOCK cluster (not fork),
#     so each of the 8 workers re-loads terra/trend and re-opens the VRT's
#     source tiles independently, and every chunk's data + result is
#     serialized over a socket — overhead that swamps the actual per-pixel
#     work at this raster size. Do not enable cores > 1 here without
#     re-investigating that regression first.
#   - Total measured single-threaded cost across all 6 CELESTE groups was
#     ~29 min (COC 248.5s, MOR 162.0s, TUR 226.3s, NIP 306.5s, NBE 673.7s,
#     KEN 123.9s) — a one-time cost, cached to disk per group exactly like
#     harvest_regen, not repeated on subsequent run_celeste.R sources.
#
# Usage (from run_celeste.R, after sourcing workflow/raster_attributes.R and
# workflow/CELESTE/prepare_ndvi.R):
#   source(here("workflow/CELESTE/prepare_ndvi_trend.R"))
#   prepare_ndvi_trend_rasters(group_manifest, cache_dir = cache_dir)
#
# Dependencies: terra, trend (via sens_slope_trend()), fs, glue; sens_slope_
#   trend() from workflow/raster_attributes.R and clean_ndvi_tiles() from
#   workflow/CELESTE/prepare_ndvi.R must be sourced first.
# =============================================================================

#' Rasterize the NDVI trend (Theil-Sen slope + Mann-Kendall p-value) LOI,
#' once per group, from that group's own NDVI source tiles
#'
#' Skips any group whose output raster already exists (cache_exists()) —
#' safe to call on every run_celeste.R source.
#'
#' @param group_manifest sf tibble from build_group_manifest() — only
#'   group_id is used (to know which groups to process); no AOI cropping is
#'   done here (see design note above — masking to crop_to happens
#'   implicitly downstream when the hydroweight stage crops the LOI to each
#'   site's catchment).
#' @param ndvi_dir  Character. Directory of raw per-region NDVI tiles, one
#'   file per HydroBasins-ish export region, filename containing the
#'   group_id (case-insensitive, optionally prefixed "Celeste_", optionally
#'   suffixed with a tile number — e.g. "Landsat_NDVI_1984_2025_Celeste_
#'   Nip3.tif"). Default here::here("data/ndvi").
#' @param cache_dir Character. Project cache root — rasters are written to
#'   cache_dir/hydroweight_loi/ndvi_trend/<group_id>.tif.
#' @param years     Integer vector, same length/order as each tile's band
#'   count — the calendar year each band represents. Default 1984:2025
#'   (CELESTE's NDVI tiles' native range; passed through to sens_slope_
#'   trend()'s `x` so the slope is per-year, not per-band-index).
#' @param min_obs   Passed to sens_slope_trend() — minimum non-NA years
#'   required to fit a trend for a given pixel. Default 4.
#'
#' @return Character vector of the (existing-or-newly-written) per-group
#'   raster paths, invisibly. A group with no matching NDVI tiles is
#'   skipped with a warning (same pattern as harvest_regen's MOR gap — the
#'   hydroweight stage's "all NA after crop/mask" check excludes it
#'   downstream rather than erroring).
prepare_ndvi_trend_rasters <- function(
  group_manifest,
  ndvi_dir = here::here("data/ndvi"),
  cache_dir,
  years = 1984:2025,
  min_obs = 4
) {
  out_dir <- fs::path(cache_dir, "hydroweight_loi", "ndvi_trend")
  files <- clean_ndvi_tiles(ndvi_dir = ndvi_dir, cache_dir = cache_dir)
  written <- character(0)

  for (grp in unique(group_manifest$group_id)) {
    group_raster_path <- fs::path(out_dir, paste0(grp, ".tif"))
    written <- c(written, group_raster_path)
    if (cache_exists(group_raster_path)) next

    matched <- match_group_tiles(files, grp)
    if (length(matched) == 0) {
      cw_warn(glue::glue(
        "No NDVI tiles matched for group '{grp}'; skipping NDVI trend raster ",
        "(will read as NA downstream, same as harvest_regen's MOR gap)."
      ))
      next
    }

    fs::dir_create(out_dir)

    # Mosaic only this group's own tiles, not a crop of the full
    # multi-province mosaic VRT — see design note above for why this
    # matters for runtime, not just tidiness.
    mini_vrt <- terra::vrt(
      matched,
      filename = tempfile(fileext = ".vrt"),
      options = c("-vrtnodata", "nan"),
      overwrite = TRUE
    )

    cw_inform(glue::glue(
      "Computing NDVI trend for group '{grp}' ({length(matched)} tile(s), ",
      "{terra::ncell(mini_vrt)} cells) — this can take several minutes."
    ))

    trend <- sens_slope_trend(mini_vrt, x = years, min_obs = min_obs)

    terra::writeRaster(trend, group_raster_path, overwrite = TRUE)
  }

  invisible(written)
}
