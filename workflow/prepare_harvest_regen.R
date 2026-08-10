# prepare_harvest_regen.R
# =============================================================================
# Prepares the Ontario harvest/regen disturbance LOI raster for the CELESTE
# stream hydroweight stage (workflow/run_celeste.R, Stage 9). Kept as its
# own script rather than inline in the runner — it's a distinct, sometimes
# slow (~1 minute+ per group, dominated by terra::rasterize() over each
# group's full DEM extent) one-time data-prep step, not part of the
# interactive delineation stage sequence.
#
# Source: shared_data/raw/harvest/ontario_harvest.gdb (Ontario only). Only
# NIP, TUR, and KEN CELESTE sites fall inside Ontario's harvest-tracked
# extent (confirmed by real spatial-filter feature counts against the gdb,
# not just a bounding-box check) — COC/MOR/NBE sites get no coverage and
# are silently excluded downstream by the hydroweight module's existing
# "all NA after crop/mask" check, same pattern as the NDVI/NBE gap.
#
# Rasterized once PER GROUP, not once for all of Ontario — like the NDVI
# VRT, the combined extent across groups is far larger than actual site
# coverage. See workflow/raster_attributes.R's rasterize_competing_classes()
# for the rasterization + temporal-precedence logic itself; this script is
# just the CELESTE/Ontario-specific configuration and per-group driver loop.
#
# Two competing classes (harvest, regen), most-recent-AR_YEAR-wins where
# they'd otherwise overlap. Ties — which is EVERY within-year overlap,
# since a per-year band only ever includes same-year features from both
# classes by construction — go to whichever bucket is named LAST in
# `buckets`; regen is last by default, so regen wins ties, matching
# "harvested then regenerated -> regen" as the more ecologically current
# state.
#
# CC-only by default (Harvest_CC02 + Harvest_CC17 = the full 2002-2024
# clear-cut record, two non-overlapping vintage exports of the same
# product). To include seed-tree/selection and shelterwood harvest too,
# pass a `buckets` list with Harvest_SE02/17 and Harvest_SH02/17 added to
# the harvest vector — no other code changes needed.
#
# Usage (from run_celeste.R, after sourcing workflow/raster_attributes.R):
#   source(here("workflow/prepare_harvest_regen.R"))
#   prepare_harvest_regen_rasters(group_manifest, cache_dir = cache_dir)
#
# Dependencies: terra, fs (via utils.R); rasterize_competing_classes() from
#   workflow/raster_attributes.R must be sourced first.
# =============================================================================

#' Class levels for the harvest/regen categorical layer
harvest_regen_levels <- data.frame(ID = 0:2, Class = c("other", "harvest", "regen"))

#' Rasterize the Ontario harvest/regen disturbance layer, once per group
#'
#' Skips any group whose output raster already exists (cache_exists()) —
#' safe to call on every run_celeste.R source.
#'
#' @param group_manifest sf tibble from build_group_manifest() — used to
#'   resolve each group's AOI (for spatially filtering the gdb) and cache
#'   directory (for the group's DEM template).
#' @param groups         Character vector of group_ids to prepare. Must fall
#'   within Ontario's harvest-tracked extent. Default c("NIP", "TUR", "KEN")
#'   — confirmed via real spatial-filter feature counts to be the only
#'   CELESTE groups inside it.
#' @param gdb_path       Character. Path to ontario_harvest.gdb.
#' @param cache_dir      Character. Project cache root — rasters are written
#'   to cache_dir/hydroweight_loi/harvest_regen/<group_id>.tif, and each
#'   group's DEM template is read from cache_dir/<group_id>/dem_breached.tif.
#' @param buckets        Named list passed to rasterize_competing_classes().
#'   Default is CC-only, per confirmed scope (see file header). Order
#'   matters: the LAST-named bucket wins ties.
#'
#' @return Character vector of the (existing-or-newly-written) per-group
#'   raster paths, invisibly.
prepare_harvest_regen_rasters <- function(
  group_manifest,
  groups = c("NIP", "TUR", "KEN"),
  gdb_path = "~/Documents/cfs/shared_data/raw/harvest/ontario_harvest.gdb",
  cache_dir,
  buckets = list(
    harvest = c("Harvest_CC02", "Harvest_CC17"),
    regen   = c("Regen_Seed", "Regen_Natural", "Regen_Plant")
  )
) {
  gdb_path <- path.expand(gdb_path)
  written <- character(0)

  for (grp in groups) {
    group_raster_path <- fs::path(
      cache_dir, "hydroweight_loi", "harvest_regen", paste0(grp, ".tif")
    )
    written <- c(written, group_raster_path)
    if (cache_exists(group_raster_path)) next

    fs::dir_create(fs::path_dir(group_raster_path))
    grp_aoi <- group_manifest[group_manifest$group_id == grp, ]
    grp_template <- terra::rast(fs::path(cache_dir, grp, "dem_breached.tif"))

    hr <- rasterize_competing_classes(
      buckets = buckets,
      gdb_path = gdb_path,
      template = grp_template,
      crop_to = grp_aoi
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
  }

  invisible(written)
}
