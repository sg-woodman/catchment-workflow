# prepare_harvest_regen.R
# =============================================================================
# Prepares the harvest/regen disturbance LOI raster for the CAM stream
# hydroweight stage (workflow/CAM/run_cam_streams.R, Stage 7). Adapted from
# workflow/CELESTE/prepare_harvest_regen.R's prepare_harvest_regen_rasters()
# / rasterize_competing_classes() approach, simplified for two real
# structural differences from CELESTE:
#
#   1. Single source, single group. CAM's AOI is entirely inside Ontario
#      MNR tenure — the Irving/NB source CELESTE also uses (NBE/COC) has no
#      relevance here, so this script only ever reads
#      shared_data/raw/harvest/ontario_harvest.gdb. And since
#      run_cam_streams.R's grouping$strategy = "whole_domain", there is
#      exactly one group ("CAM_streams", = config$project_id) rather than
#      CELESTE's six HydroBasins-derived groups.
#   2. CAM's group_manifest has a FLAT cache_dir (`group_manifest$cache_dir`
#      equals `config$cache_dir` directly — no per-group subfolder; see
#      workflow/R/engine/01_build_group_manifest.R's "whole_domain"
#      strategy). The per-group DEM template is therefore read straight
#      from `cache_dir/dem_breached.tif`, NOT
#      `cache_dir/<group_id>/dem_breached.tif` the way CELESTE's script
#      reconstructs it — that CELESTE-style path would be wrong here.
#
# Real (not bounding-box) coverage of ontario_harvest.gdb over CAM's AOI
# was confirmed directly before writing this script — CLAUDE.md documents
# that a bbox-only check gave false positives for this exact dataset in
# CELESTE's case, so it wasn't skipped here either. Feature counts inside
# the CAM AOI (+1km buffer), after the same geometry cleaning
# rasterize_competing_classes() itself does (st_cast/st_make_valid — the
# raw layers include curved MULTISURFACE geometries that break a naive
# st_filter()):
#   Harvest_CC02: 38,999 / 204,475      Harvest_CC17: 26,613 / 155,472
#   Regen_Seed:    3,096 /  17,833      Regen_Natural: 32,423 / 136,860
#   Regen_Plant:  26,506 / 122,899
# All five real, substantial coverage — no gap analogous to CELESTE's MOR.
#
# year_range 2002-2024 (Ontario's native range for these layers, same as
# CELESTE's Ontario source) — no NB source here, so no cross-source
# comparability concern to reconcile.
#
# Usage (from run_cam_streams.R, after sourcing workflow/raster_attributes.R):
#   source(here("workflow/CAM/prepare_harvest_regen.R"))
#   prepare_cam_harvest_regen_rasters(group_manifest, cache_dir = cache_dir)
#
# Dependencies: terra, sf, dplyr, fs (via utils.R); rasterize_competing_
#   classes() from workflow/raster_attributes.R must be sourced first.
# =============================================================================

#' Class levels for the harvest/regen categorical layer — same scheme
#' CELESTE's harvest_regen LOI uses, so the two projects' output columns
#' are directly comparable if ever analyzed together.
harvest_regen_levels <- data.frame(ID = 0:2, Class = c("other", "harvest", "regen"))

#' Rasterize the harvest/regen disturbance layer for CAM's one whole-domain
#' group, from the Ontario MNR harvest gdb
#'
#' Skips if the output raster already exists (cache_exists()) — safe to
#' call on every run_cam_streams.R source.
#'
#' @param group_manifest sf tibble from build_engine_group_manifest() — used
#'   to resolve the group's AOI (for spatially filtering source layers) and
#'   cache directory (for the DEM template). Expected to have exactly one
#'   row (whole_domain grouping).
#' @param cache_dir      Character. Project cache root — the raster is
#'   written to cache_dir/hydroweight_loi/harvest_regen/<group_id>.tif, and
#'   the DEM template is read from cache_dir/dem_breached.tif directly (flat
#'   layout — see header).
#' @param gdb_path       Character. Path to the Ontario MNR harvest gdb.
#' @param buckets        Named list passed to rasterize_competing_classes().
#'   Default: CC-only harvest scope (Harvest_CC02 + Harvest_CC17 = the full
#'   2002-2024 clear-cut record, two non-overlapping vintage exports of the
#'   same product), matching CELESTE's Ontario harvest scope.
#' @param year_range     Numeric length-2. Passed to
#'   rasterize_competing_classes(). Default c(2002, 2024) — Ontario's native
#'   range for these layers.
#'
#' @return The (existing-or-newly-written) raster path, invisibly.
prepare_cam_harvest_regen_rasters <- function(
  group_manifest,
  cache_dir,
  gdb_path = "~/Documents/cfs/shared_data/raw/harvest/ontario_harvest.gdb",
  buckets = list(
    harvest = c("Harvest_CC02", "Harvest_CC17"),
    regen   = c("Regen_Seed", "Regen_Natural", "Regen_Plant")
  ),
  year_range = c(2002, 2024)
) {
  if (nrow(group_manifest) != 1) {
    cw_abort(paste(
      "prepare_cam_harvest_regen_rasters() expects a single-row",
      "group_manifest (grouping$strategy = 'whole_domain'); got",
      nrow(group_manifest), "rows."
    ))
  }

  group_id <- group_manifest$group_id[1]
  group_raster_path <- fs::path(
    cache_dir, "hydroweight_loi", "harvest_regen", paste0(group_id, ".tif")
  )
  if (cache_exists(group_raster_path)) {
    return(invisible(group_raster_path))
  }

  fs::dir_create(fs::path_dir(group_raster_path))
  gdb_path <- path.expand(gdb_path)
  # Flat cache_dir (see header) — NOT cache_dir/<group_id>/dem_breached.tif.
  grp_template <- terra::rast(fs::path(cache_dir, "dem_breached.tif"))

  hr <- rasterize_competing_classes(
    buckets = buckets,
    gdb_path = gdb_path,
    template = grp_template,
    year_range = year_range,
    crop_to = group_manifest
  )
  terra::writeRaster(
    hr,
    group_raster_path,
    overwrite = TRUE,
    datatype = "INT1U",
    gdal = "PHOTOMETRIC=MINISBLACK" # avoid a real GDAL quirk: an
    # exactly-3-band Byte GeoTIFF gets auto-tagged as RGB, silently
    # renaming bands to red/green/blue on every subsequent read
  )

  invisible(group_raster_path)
}
