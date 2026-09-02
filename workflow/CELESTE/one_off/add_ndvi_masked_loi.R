# add_ndvi_masked_loi.R
# =============================================================================
# One-off: adds the "ndvi_masked"/"ndvi_trend_masked" LOIs (see
# workflow/CELESTE/README.md's "NDVI lake masking" note and
# prepare_ndvi_masked.R's header) alongside the existing "ndvi"/
# "ndvi_trend" LOIs, and regenerates the genuine unmasked "ndvi"/
# "ndvi_trend" baseline — without re-running any catchment-producing stage
# (1-6). Same discipline as add_canlcc_loi.R: Stages 1-6 are untouched, so
# no catchment.gpkg, catchment_clipped.gpkg, all_catchments*.gpkg, or
# catchment_metrics.csv is read for writing or modified.
#
# SUPERSEDES an earlier (same-session) attempt that masked "ndvi"/
# "ndvi_trend" IN PLACE, overwriting their cache — reverted because it
# made a real, unrelated problem (MOR-CRANE/MOR-KEN returning no data at
# all, for a hydroweight()-level reason having nothing to do with
# masking) look like a masking bug. This script instead:
#   1. Regenerates the genuine UNMASKED "ndvi"/"ndvi_trend" baseline
#      (prepare_ndvi_per_group_rasters()/prepare_ndvi_trend_rasters(),
#      now reverted to their original unmasked-only behavior) — the one
#      unavoidable recompute, since the pre-masking unmasked trend was
#      overwritten in place and never separately cached.
#   2. Prepares "ndvi_masked"/"ndvi_trend_masked" (prepare_ndvi_masked.R)
#      — cheap here, since the masked group rasters from the earlier
#      in-place attempt were preserved (moved to cache_dir/hydroweight_
#      loi/ndvi_masked/ and .../ndvi_trend_masked/ before any code
#      changed) rather than discarded, so cache_exists() finds them
#      immediately with no NHN re-fetch or sens_slope_trend() re-run.
#   3. Recomputes hydroweight for all 4 LOIs, replaces every existing
#      "ndvi_*"/"ndvi_trend_*" column (currently the STALE, masked-in-
#      place values from the reverted attempt) with the fresh set — the
#      genuinely unmasked ones plus the new "ndvi_masked_*"/
#      "ndvi_trend_masked_*" ones — by (site, version).
#   4. Regenerates output_dir/tidy/.
#
# POST-RUN NOTE: the full-batch STEP 3 call above left 5 scattered site x
# version rows with no ndvi-family data at all (not the same sites twice
# across CELESTE/CAM, and not reproducible in isolation) — a terra temp-
# file/GC-race class issue already documented elsewhere in this project
# for large multi-LOI runs, not a real per-site problem. Re-running
# calculate_hydroweight_attributes_stream() for just the affected
# site_ids resolved every case; see workflow/CELESTE/README.md's "NDVI
# lake masking" section for the full writeup. Check the final
# CELESTE_hydroweight.csv for any all-NA ndvi-family row after running
# this script, and retry that subset the same way if any turn up.
#
# Usage: source from an R session in the project root, after
# output/CELESTE/CELESTE_hydroweight.csv already exists.
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
source(here("workflow/R/stream/group_sites.R"))
source(here("workflow/R/stream/burn_streams.R")) # find_nhn_sheets()/read_nhn_from_gdb() — needed by fetch_nhn_lakes_for_aoi()
source(here("workflow/R/stream/hydroweight_attributes.R"))
source(here("workflow/R/engine/00_resolve_config.R"))
source(here("workflow/R/engine/01_build_group_manifest.R"))
source(here("workflow/R/engine/prepare_lake_conditioning.R")) # fetch_nhn_lakes_for_aoi()
source(here("workflow/raster_attributes.R")) # mask_out_waterbodies()
source(here("workflow/CELESTE/prepare_ndvi.R"))
source(here("workflow/CELESTE/prepare_ndvi_trend.R"))
source(here("workflow/CELESTE/prepare_ndvi_masked.R"))
source(here("workflow/CELESTE/tidy_outputs.R"))

# =============================================================================
# STAGE 1 (config + group manifest only — no terrain/delineation calls)
# =============================================================================

PROJECT_ID <- "CELESTE"
output_dir <- here("output", PROJECT_ID)
cache_dir  <- here("cache", PROJECT_ID)

mrdem_vrt       <- "~/Documents/cfs/shared_data/raw/dem/mrdem-30-dtm.vrt"
hydrobasins_dir <- "/Users/sam/Documents/cfs/shared_data/raw/hydro/watersheds/HydroBasins"
nhn_dir         <- "/Users/sam/Documents/cfs/shared_data/raw/hydro/networks/NHN/gdb"
nhn_index       <- "/Users/sam/Documents/cfs/shared_data/raw/hydro/networks/NHN/NHN_INDEX_WORKUNIT_LIMIT_2/NHN_INDEX_22_INDEX_WORKUNIT_LIMIT_2.shp"

sites_sf <- st_read(here("data/celeste_milli_sites_clean_corrected.gpkg"), quiet = TRUE)
coords <- st_coordinates(sites_sf)
sites <- sites_sf |>
  st_drop_geometry() |>
  as_tibble() |>
  mutate(lon = coords[, "X"], lat = coords[, "Y"]) |>
  select(site_id, site_name, lon, lat, group_id, burn_streams, aoi_buffer_m) |>
  mutate(burn_streams = dplyr::if_else(group_id == "COC", FALSE, burn_streams))

run_config <- list(
  project_id = PROJECT_ID,
  output_dir = output_dir,
  cache_dir  = cache_dir,
  sites      = sites,

  dem = list(path = mrdem_vrt),
  flow_direction = NULL,
  flow_pointer   = NULL,
  crs = NULL,
  stream_threshold = 1000,

  streams_burn   = list(source = "nhn_auto"),
  nhn_index_path = nhn_index,
  nhn_raw_dir    = nhn_dir,

  lake_conditioning = list(source = "nhn_auto"),

  lake_polygons = NULL,
  lake_buffer_m = 30,

  grouping = list(
    strategy = "hydrobasins",
    hydrobasins_dir = hydrobasins_dir,
    default_buffer_m = 1000
  ),

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
prepare_ndvi_per_group_rasters(group_manifest, cache_dir = cache_dir)

cw_inform("Regenerating ndvi_trend (unmasked baseline) — this can take ~30 min...")
prepare_ndvi_trend_rasters(group_manifest, cache_dir = cache_dir)

# =============================================================================
# STEP 2 — prepare ndvi_masked/ndvi_trend_masked (cheap — already-computed
# masked group rasters from the reverted attempt were preserved, moved
# into place before any code changed).
# =============================================================================

cw_inform("Preparing ndvi_masked (reusing preserved masked group rasters if present)...")
prepare_ndvi_masked_group_rasters(
  group_manifest, cache_dir = cache_dir,
  nhn_index_path = nhn_index, nhn_raw_dir = nhn_dir
)

cw_inform("Preparing ndvi_trend_masked (reusing preserved masked group rasters if present)...")
prepare_ndvi_trend_masked_rasters(
  group_manifest, cache_dir = cache_dir,
  nhn_index_path = nhn_index, nhn_raw_dir = nhn_dir
)

# =============================================================================
# STEP 3 — recompute hydroweight for all 4 ndvi-family LOIs
# =============================================================================

ndvi_family_layers <- list(
  list(
    path_lazy = fs::path(cache_dir, "hydroweight_loi", "ndvi", "{group_id}.tif"),
    name = "ndvi", type = "continuous",
    stats = c("distwtd_mean", "distwtd_sd", "mean", "sd", "median", "min", "max")
  ),
  list(
    path_lazy = fs::path(cache_dir, "hydroweight_loi", "ndvi_trend", "{group_id}.tif"),
    name = "ndvi_trend", type = "continuous",
    stats = c("distwtd_mean", "distwtd_sd", "mean", "sd", "median", "min", "max")
  ),
  list(
    path_lazy = fs::path(cache_dir, "hydroweight_loi", "ndvi_masked", "{group_id}.tif"),
    name = "ndvi_masked", type = "continuous",
    stats = c("distwtd_mean", "distwtd_sd", "mean", "sd", "median", "min", "max")
  ),
  list(
    path_lazy = fs::path(cache_dir, "hydroweight_loi", "ndvi_trend_masked", "{group_id}.tif"),
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
  catchment_versions = c("unclipped", "clipped")
)

cw_inform(glue::glue("ndvi family recomputed for {nrow(ndvi_results)} site x version row(s)."))

# =============================================================================
# STEP 4 — REPLACE every ndvi_*/ndvi_trend_* column in CELESTE_hydroweight.csv
# (the existing ones are stale — masked-in-place values from the reverted
# attempt) with the fresh set (genuinely unmasked + the two new masked LOIs)
# =============================================================================

hydroweight_path <- fs::path(output_dir, "CELESTE_hydroweight.csv")
if (!fs::file_exists(hydroweight_path)) {
  cw_abort(glue::glue("Not found: {hydroweight_path} — run run_celeste.R first."))
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

tidy_celeste_outputs(output_dir = output_dir)
