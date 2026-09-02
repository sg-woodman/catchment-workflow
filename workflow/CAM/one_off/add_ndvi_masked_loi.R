# add_ndvi_masked_loi.R
# =============================================================================
# One-off: adds the "ndvi_masked"/"ndvi_trend_masked" LOIs (see
# workflow/CAM/README.md's "NDVI lake masking" note and
# prepare_ndvi_masked.R's header) alongside the existing "ndvi"/
# "ndvi_trend" LOIs, and regenerates the genuine unmasked "ndvi"/
# "ndvi_trend" baseline — without re-running any catchment-producing stage
# (1-6). Same discipline as workflow/CELESTE/one_off/add_canlcc_loi.R /
# add_ndvi_masked_loi.R: Stages 1-6 are untouched, so no catchment.gpkg,
# catchment_clipped.gpkg, all_catchments*.gpkg, or catchment_metrics.csv
# is read for writing or modified.
#
# SUPERSEDES an earlier (same-session) attempt that masked "ndvi"/
# "ndvi_trend" IN PLACE (inside clean_cam_ndvi_tile()), overwriting their
# cache — reverted for the same reason as CELESTE's equivalent revert (see
# workflow/CELESTE/one_off/add_ndvi_masked_loi.R's header): it makes a
# real problem look like a masking bug when there's no unmasked baseline
# left to compare against. This script instead:
#   1. Regenerates the genuine UNMASKED "ndvi"/"ndvi_trend" baseline
#      (prepare_cam_ndvi_site_rasters()/prepare_cam_ndvi_trend_site_
#      rasters(), now reverted to their original unmasked-only behavior)
#      — the one unavoidable recompute, since the pre-masking unmasked
#      trend was overwritten in place and never separately cached.
#   2. Prepares "ndvi_masked"/"ndvi_trend_masked" (prepare_ndvi_masked.R)
#      — cheap here, since the masked files from the earlier in-place
#      attempt were preserved (moved to cache_dir/hydroweight_loi/
#      ndvi_masked_clean/ and .../ndvi_trend_src_masked/ before any code
#      changed) rather than discarded, so cache_exists() finds them
#      immediately with no re-mask or sens_slope_trend() re-run.
#   3. Recomputes hydroweight for all 4 LOIs, replaces every existing
#      "ndvi_*"/"ndvi_trend_*" column (currently the STALE, masked-in-
#      place values from the reverted attempt) with the fresh set — the
#      genuinely unmasked ones plus the new "ndvi_masked_*"/
#      "ndvi_trend_masked_*" ones — by (site, version).
#   4. Regenerates output_dir/tidy/.
#
# POST-RUN NOTE: the full-batch STEP 3 call above left 17 scattered site
# x version rows with no ndvi-family data at all, including several with
# no plausible "lake-dominated catchment" story (Bell, Marina,
# Smoothwater, SUD12) — a terra temp-file/GC-race class issue already
# documented elsewhere in this project for large multi-LOI runs, not a
# real per-site problem. Re-running calculate_hydroweight_attributes_
# stream() for just the affected site_ids resolved every case; see
# workflow/CAM/README.md's "NDVI lake masking" section for the full
# writeup. Check the final CAM_streams_hydroweight.csv for any all-NA
# ndvi-family row after running this script, and retry that subset the
# same way if any turn up.
#
# Usage: source from an R session in the project root, after
# output/CAM/stream_delineation/CAM_streams_hydroweight.csv already
# exists.
# =============================================================================

library(sf)
library(terra)
library(whitebox)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(tibble)
library(fs)
library(cli)
library(glue)
library(here)

source(here("workflow/R/utils.R"))
source(here("workflow/R/stream/burn_streams.R"))
source(here("workflow/R/stream/delineate_sites.R"))
source(here("workflow/R/stream/hydroweight_attributes.R"))
source(here("workflow/R/engine/00_resolve_config.R"))
source(here("workflow/R/engine/01_build_group_manifest.R"))
source(here("workflow/gee_utils.R")) # prepare_polygons_for_gee()/group_polygons_by_adjacency() — used by build_cam_ndvi_site_map()
source(here("workflow/raster_attributes.R")) # sens_slope_trend(), mask_out_waterbodies()
source(here("workflow/CAM/prepare_ndvi.R"))
source(here("workflow/CAM/prepare_ndvi_trend.R"))
source(here("workflow/CAM/prepare_ndvi_masked.R"))
source(here("workflow/CAM/tidy_outputs.R"))

# =============================================================================
# STAGE 1 (config + group manifest only — no terrain/delineation calls)
# =============================================================================

PROJECT_ID  <- "CAM"
DELINEATION <- "stream_delineation"

output_dir <- here("output", PROJECT_ID, DELINEATION)
cache_dir  <- here("cache", PROJECT_ID, DELINEATION)

OIH_DEM_PATH      <- "/Users/sam/Downloads/IntegratedHydrologyNE/EnforcedDEM.tif"
OIH_FLOW_DIR_PATH <- "/Users/sam/Downloads/IntegratedHydrologyNE/EnhancedFlowDirection.tif"

oih_recode_matrix <- matrix(
  c(128, 1, 1, 2, 2, 4, 4, 8, 8, 16, 16, 32, 32, 64, 64, 128),
  ncol = 2, byrow = TRUE
)

STREAM_THRESHOLD <- 100

EXCLUDED_SITE_IDS <- c("SUD11")

sites_raw <- readr::read_csv(here("data/cam_stream_sites_raw.csv"), show_col_types = FALSE)

sites <- sites_raw |>
  dplyr::mutate(
    site_id = stream_id |>
      gsub("\\s+", "_", x = _) |>
      gsub("[^A-Za-z0-9_-]", "", x = _),
    site_name = stream_id
  ) |>
  dplyr::filter(!site_id %in% EXCLUDED_SITE_IDS) |>
  dplyr::select(site_id, site_name, lon, lat)

run_config <- list(
  project_id = "CAM_streams",
  output_dir = output_dir,
  cache_dir  = cache_dir,
  sites      = sites,

  dem = list(path = OIH_DEM_PATH),
  flow_direction = list(path = OIH_FLOW_DIR_PATH, recode = oih_recode_matrix),
  flow_pointer = NULL,
  flow_accum   = NULL,
  crs = NULL,
  stream_threshold = STREAM_THRESHOLD,

  streams_burn = list(source = "none"),

  lake_polygons = NULL,
  lake_buffer_m = 30,

  grouping = list(strategy = "whole_domain"),

  loi_layers = NULL
)

config <- resolve_engine_config(run_config)

gm <- build_engine_group_manifest(config)
sites <- gm$sites
group_manifest <- gm$group_manifest

cw_inform(glue::glue(
  "Config + group manifest resolved (Stage 1 only): {nrow(sites)} site(s), ",
  "{nrow(group_manifest)} group(s). No terrain/delineation stage touched."
))

# =============================================================================
# STEP 1 — regenerate the genuine UNMASKED baseline (now reverted code)
# =============================================================================

cw_inform("Regenerating ndvi (unmasked baseline)...")
prepare_cam_ndvi_site_rasters(sites, output_dir = output_dir, cache_dir = cache_dir)

cw_inform("Regenerating ndvi_trend (unmasked baseline) — this can take ~45 min...")
prepare_cam_ndvi_trend_site_rasters(sites, output_dir = output_dir, cache_dir = cache_dir)

# =============================================================================
# STEP 2 — prepare ndvi_masked/ndvi_trend_masked (cheap — already-computed
# masked files from the reverted attempt were preserved, moved into place
# before any code changed).
# =============================================================================

cw_inform("Preparing ndvi_masked (reusing preserved masked files if present)...")
prepare_cam_ndvi_masked_site_rasters(sites, output_dir = output_dir, cache_dir = cache_dir)

cw_inform("Preparing ndvi_trend_masked (reusing preserved masked files if present)...")
prepare_cam_ndvi_trend_masked_site_rasters(sites, output_dir = output_dir, cache_dir = cache_dir)

# =============================================================================
# STEP 3 — recompute hydroweight for all 4 ndvi-family LOIs
# =============================================================================

ndvi_family_layers <- list(
  list(
    path_template = fs::path(cache_dir, "hydroweight_loi", "ndvi", "{site_id}.tif"),
    name = "ndvi", type = "continuous"
  ),
  list(
    path_template = fs::path(cache_dir, "hydroweight_loi", "ndvi_trend", "{site_id}.tif"),
    name = "ndvi_trend", type = "continuous",
    stats = c("distwtd_mean", "distwtd_sd", "mean", "sd", "median", "min", "max")
  ),
  list(
    path_template = fs::path(cache_dir, "hydroweight_loi", "ndvi_masked", "{site_id}.tif"),
    name = "ndvi_masked", type = "continuous"
  ),
  list(
    path_template = fs::path(cache_dir, "hydroweight_loi", "ndvi_trend_masked", "{site_id}.tif"),
    name = "ndvi_trend_masked", type = "continuous",
    stats = c("distwtd_mean", "distwtd_sd", "mean", "sd", "median", "min", "max")
  )
)

ndvi_results <- calculate_hydroweight_attributes_stream(
  sites = sites,
  group_manifest = group_manifest,
  output_dir = output_dir,
  cache_dir = cache_dir,
  loi_layers = ndvi_family_layers,
  catchment_versions = c("unclipped", "clipped"),
  raster_crs = config$working_crs # EPSG:3161 — not the function's default EPSG:3979
)

cw_inform(glue::glue("ndvi family recomputed for {nrow(ndvi_results)} site x version row(s)."))

# =============================================================================
# STEP 4 — REPLACE every ndvi_*/ndvi_trend_* column in
# CAM_streams_hydroweight.csv (the existing ones are stale — masked-in-
# place values from the reverted attempt) with the fresh set (genuinely
# unmasked + the two new masked LOIs)
# =============================================================================

hydroweight_path <- fs::path(output_dir, "CAM_streams_hydroweight.csv")
if (!fs::file_exists(hydroweight_path)) {
  cw_abort(glue::glue("Not found: {hydroweight_path} — run run_cam_streams.R first."))
}

hw_existing <- readr::read_csv(hydroweight_path, show_col_types = FALSE)

old_ndvi_cols <- grep("^ndvi_", names(hw_existing), value = TRUE)
cw_inform(glue::glue(
  "Dropping {length(old_ndvi_cols)} stale ndvi-family column(s) before rejoining."
))
hw_stripped <- dplyr::select(hw_existing, -dplyr::all_of(old_ndvi_cols))

new_ndvi_cols <- setdiff(names(ndvi_results), c("site", "version"))
hw_merged <- dplyr::left_join(hw_stripped, ndvi_results, by = c("site", "version"))

n_before <- nrow(hw_existing)
n_after  <- nrow(hw_merged)
if (n_after != n_before) {
  cw_abort(glue::glue(
    "Row count changed after join ({n_before} -> {n_after}) — the join key ",
    "(site, version) didn't line up 1:1. Aborting without writing."
  ))
}

n_missing <- sum(is.na(hw_merged[[new_ndvi_cols[1]]]))
if (n_missing > 0) {
  cw_warn(glue::glue(
    "{n_missing} existing row(s) got no ndvi match (site/version present in ",
    "the CSV but not recomputed) — check for a group/site mismatch before ",
    "trusting the tidy ndvi tables fully."
  ))
}

readr::write_csv(hw_merged, hydroweight_path)
cw_inform(glue::glue(
  "ndvi family replaced in {hydroweight_path} ({n_before} row(s) unchanged, ",
  "{length(old_ndvi_cols)} old column(s) dropped, {length(new_ndvi_cols)} new ",
  "column(s) added)."
))

# =============================================================================
# STEP 5 — regenerate tidy tables
# =============================================================================

tidy_cam_outputs(output_dir = output_dir)
