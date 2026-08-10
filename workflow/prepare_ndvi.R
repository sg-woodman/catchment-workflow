# prepare_ndvi.R
# =============================================================================
# Prepares the CELESTE NDVI (continuous) LOI raster for the stream
# hydroweight stage (workflow/run_celeste.R, Stage 9). Kept as its own
# script — same reasoning as workflow/prepare_harvest_regen.R: this is a
# workflow-specific data-prep concern (source data, how it's tiled, what
# needs mosaicking/reprojecting), and how it's handled will genuinely
# differ between projects (a continuous time-series mosaic here; a
# rasterized-per-group categorical layer with temporal precedence for
# harvest/regen; something else entirely for a future project) — none of
# that belongs baked into the standardized Stage 1-9 runner or the generic
# hydroweight module.
#
# Source: data/ndvi/*.tif — 12 regional tiles (COC, KEN1-4, MOR, NIP1-5,
# TUR), 1984-2025 annual composites, all already EPSG:3979/30m/42-band-
# aligned (verified with check_tile_consistency() — no reprojection or band
# alignment needed here; a future batch that isn't this consistent would
# need that handled in this script before mosaicking, not downstream).
#
# Mosaicked into one VRT via build_mosaic_vrt(), which sets vrtnodata so
# areas outside every tile read as NA, not 0 — without that fix a gap reads
# as valid-looking zeros indistinguishable from real data. Returned for use
# via loi_layers[[i]]$path_lazy rather than $path because the 12 tiles'
# combined bounding box is transcontinental (CELESTE sites span NWT,
# Saskatchewan, Ontario, Quebec, Newfoundland) even though actual coverage
# is small clustered patches — caching that whole extent as one reprojected
# file would be enormous for no benefit.
#
# KNOWN GAP: the NBE group (20 sites) has no NDVI tile at all — confirmed
# empirically (every NBE site's LOI crop is 100% NA after the vrtnodata
# fix). Those sites are silently excluded downstream by the hydroweight
# module's existing "all NA after crop/mask" check, not an error. Add an
# NBE tile to data/ndvi/ once available and rerun.
#
# Usage (from run_celeste.R, after sourcing workflow/raster_attributes.R):
#   source(here("workflow/prepare_ndvi.R"))
#   ndvi_vrt_path <- prepare_ndvi_vrt(cache_dir = cache_dir)
#
# Dependencies: terra, fs (via utils.R); build_mosaic_vrt() from
#   workflow/raster_attributes.R must be sourced first.
# =============================================================================

#' Mosaic the NDVI tiles into one VRT for use as a hydroweight path_lazy LOI
#'
#' Always rebuilt (not cache_exists()-gated like
#' prepare_harvest_regen_rasters()'s per-group rasters) — a VRT is a
#' lightweight XML pointer, not materialized data, so rebuilding is cheap,
#' and skipping would mean a newly-added tile in ndvi_dir (this has already
#' happened once: two NBE tiles added after the initial NDVI batch) doesn't
#' get picked up until something manually deletes the stale VRT.
#'
#' @param ndvi_dir  Character. Directory of source NDVI tiles. Default
#'   here("data/ndvi").
#' @param cache_dir Character. Project cache root — VRT is written to
#'   cache_dir/hydroweight_loi/ndvi_mosaic.vrt.
#'
#' @return Character path to the mosaic VRT, invisibly.
prepare_ndvi_vrt <- function(ndvi_dir = here::here("data/ndvi"), cache_dir) {
  vrt_path <- fs::path(cache_dir, "hydroweight_loi", "ndvi_mosaic.vrt")
  fs::dir_create(fs::path_dir(vrt_path))

  build_mosaic_vrt(
    files = list.files(ndvi_dir, pattern = "[.]tif$", full.names = TRUE),
    vrt_path = vrt_path
  )

  invisible(vrt_path)
}
