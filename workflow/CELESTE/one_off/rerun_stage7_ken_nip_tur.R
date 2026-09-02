# rerun_stage7_ken_nip_tur.R
# =============================================================================
# One-off script: reruns Stage 7 (LOI prep + hydroweighting) for KEN/NIP/TUR
# only, after the 2026-08-29 upfront lake-conditioning terrain fix
# (commit ac5f908) reconditioned those 3 groups' DEMs and redelineated their
# sites but never got a Stage 7 rerun — confirmed via file timestamps:
# CELESTE_engine_hydroweight.csv (Aug 27 08:18) predates all 3 groups'
# dem_lakes_flattened.tif (Aug 29 ~14:00/15:19/15:56).
#
# Does NOT touch COC/MOR/NBE — merges KEN/NIP/TUR rows into the existing
# combined CSV via merge_rows_into_csv() rather than overwriting it wholesale
# (same pattern as engine/99_rerun_sites's rerun_engine_sites()).
#
# IMPORTANT cache-staleness gotcha this script exists to handle: every LOI
# (ndvi, ndvi_trend, harvest_regen) caches a per-SITE cropped/reprojected
# raster (cache_dir/hydroweight_loi/<loi>/<site_id>_reprojected.tif),
# cropped to that site's catchment geometry at the time it was built
# (resolve_site_loi_raster()'s crop_to = catchment_sf). Those per-site
# caches for KEN/NIP/TUR sites were built 2026-08-26, against the OLD
# (lake-bisected) catchment geometry — stale regardless of whether the LOI
# source itself is DEM-dependent. harvest_regen's GROUP-level raster
# (cache_dir/hydroweight_loi/harvest_regen/{KEN,NIP,TUR}.tif) is also stale
# since it rasterizes onto a template read from THIS run's own
# dem_breached.tif, which changed. ndvi/ndvi_trend's group-level rasters
# (then production_cache_dir, now just cache_dir after the 2026-08-31 rename) are NOT DEM-dependent and are left alone. All
# stale per-site + group caches for the 3 groups are deleted below before
# recomputing.
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
source(here("workflow/R/stream/burn_streams.R"))
source(here("workflow/R/stream/delineate_sites.R"))
source(here("workflow/R/stream/hydroweight_attributes.R"))
source(here("workflow/R/remove_upstream.R"))
source(here("workflow/R/reclip_outputs.R"))
source(here("workflow/R/catchment_metrics.R"))
source(here("workflow/R/engine/00_resolve_config.R"))
source(here("workflow/R/engine/01_build_group_manifest.R"))
source(here("workflow/R/engine/02_prepare_terrain.R"))
source(here("workflow/R/engine/03_prepare_streams_burn.R"))
source(here("workflow/R/engine/prepare_lake_conditioning.R"))
source(here("workflow/R/engine/04_delineate_site.R"))
source(here("workflow/R/engine/99_rerun_sites")) # drop_redundant_clipped_rows(), merge_rows_into_csv()
source(here("workflow/raster_attributes.R"))
source(here("workflow/CELESTE/prepare_ndvi.R"))
source(here("workflow/CELESTE/prepare_ndvi_trend.R"))
source(here("workflow/CELESTE/prepare_harvest_regen.R"))
source(here("workflow/CELESTE/tidy_outputs.R"))

# -- Config / sites / group_manifest (Stage 1, identical to run_celeste.R) ---

PROJECT_ID <- "CELESTE" # was "CELESTE_engine" when this script ran (2026-08-31);
# updated to match the output/cache consolidation done the same day — see
# run_celeste.R's PROJECT_ID comment. Left in the repo as a record of the
# KEN/NIP/TUR Stage 7 rerun; not expected to be re-run again.
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
  project_id = PROJECT_ID, output_dir = output_dir, cache_dir = cache_dir, sites = sites,
  dem = list(path = mrdem_vrt), flow_direction = NULL, flow_pointer = NULL,
  crs = NULL, stream_threshold = 1000,
  streams_burn = list(source = "nhn_auto"), nhn_index_path = nhn_index, nhn_raw_dir = nhn_dir,
  lake_polygons = NULL, lake_buffer_m = 30,
  grouping = list(strategy = "hydrobasins", hydrobasins_dir = hydrobasins_dir, default_buffer_m = 1000),
  loi_layers = NULL
)
config <- resolve_engine_config(run_config)

gm <- build_engine_group_manifest(config)
sites_all <- gm$sites
group_manifest_all <- gm$group_manifest

cw_inform(glue::glue("Full project: {nrow(sites_all)} sites, {nrow(group_manifest_all)} groups."))

# -- Stage 4 rerun (full sites — cheap, needed for fresh n_erased used by ----
# -- drop_redundant_clipped_rows; also self-verifies Stage 4/5 currency) -----

upstream_results <- remove_upstream_catchments(sites_all, output_dir)

# -- Scope down to KEN/NIP/TUR for Stage 7 -----------------------------------

target_groups <- c("KEN", "NIP", "TUR")
sites <- dplyr::filter(sites_all, group_id %in% target_groups)
group_manifest <- dplyr::filter(group_manifest_all, group_id %in% target_groups)

cw_inform(glue::glue(
  "Scoped to {nrow(sites)} site(s) across {nrow(group_manifest)} group(s): ",
  "{paste(target_groups, collapse = ', ')}."
))

# -- Invalidate stale per-site + group-level LOI caches for these sites -----

for (loi in c("ndvi", "ndvi_trend", "harvest_regen")) {
  loi_cache_dir <- fs::path(cache_dir, "hydroweight_loi", loi)
  if (!fs::dir_exists(loi_cache_dir)) next
  stale <- fs::path(loi_cache_dir, paste0(sites$site_id, "_reprojected.tif"))
  stale <- stale[fs::file_exists(stale)]
  if (length(stale) > 0) {
    cw_inform(glue::glue("Deleting {length(stale)} stale per-site '{loi}' cache file(s)..."))
    fs::file_delete(stale)
  }
}

harvest_regen_group_rasters <- fs::path(
  cache_dir, "hydroweight_loi", "harvest_regen", paste0(target_groups, ".tif")
)
harvest_regen_group_rasters <- harvest_regen_group_rasters[fs::file_exists(harvest_regen_group_rasters)]
if (length(harvest_regen_group_rasters) > 0) {
  cw_inform(glue::glue(
    "Deleting {length(harvest_regen_group_rasters)} stale group-level harvest_regen raster(s)..."
  ))
  fs::file_delete(harvest_regen_group_rasters)
}

# -- Stage 7 — LOI prep + hydroweighting, scoped to KEN/NIP/TUR -------------

cw_inform("Preparing ndvi (reusing production cache — no DEM dependency)...")
prepare_ndvi_per_group_rasters(group_manifest, cache_dir = cache_dir)

cw_inform("Preparing ndvi_trend (reusing production cache — no DEM dependency)...")
prepare_ndvi_trend_rasters(group_manifest, cache_dir = cache_dir)

cw_inform("Preparing harvest_regen (FRESH — depends on this run's own dem_breached.tif)...")
prepare_harvest_regen_rasters(group_manifest, cache_dir = cache_dir)

loi_layers <- list(
  list(
    path_lazy = fs::path(cache_dir, "hydroweight_loi", "harvest_regen", "{group_id}.tif"),
    name = "harvest_regen", type = "categorical", class_levels = harvest_regen_levels
  ),
  list(
    path_lazy = fs::path(cache_dir, "hydroweight_loi", "ndvi", "{group_id}.tif"),
    name = "ndvi", type = "continuous",
    stats = c("distwtd_mean", "distwtd_sd", "mean", "sd", "median", "min", "max")
  ),
  list(
    path_lazy = fs::path(cache_dir, "hydroweight_loi", "ndvi_trend", "{group_id}.tif"),
    name = "ndvi_trend", type = "continuous",
    stats = c("distwtd_mean", "distwtd_sd", "mean", "sd", "median", "min", "max")
  )
)

hw_results_new <- calculate_hydroweight_attributes_stream(
  sites = sites,
  group_manifest = group_manifest,
  output_dir = output_dir,
  cache_dir = cache_dir,
  loi_layers = loi_layers,
  catchment_versions = c("unclipped", "clipped")
)

cw_inform(glue::glue("Stage 7 produced {nrow(hw_results_new)} row(s) for KEN/NIP/TUR."))

# -- Merge into the existing combined CSV, never overwrite wholesale --------

hw_csv_path <- fs::path(output_dir, "CELESTE_hydroweight.csv")
merge_rows_into_csv(
  new_rows = hw_results_new,
  csv_path = hw_csv_path,
  key_cols = c("site", "version"),
  filter_fn = function(df) drop_redundant_clipped_rows(df, upstream_results, site_col = "site")
)

# -- Stage 8 — regenerate tidy long-format tables from the now-current -----
# -- wide CSVs (cheap, always a full rebuild from source — see the file's --
# -- own header) --------------------------------------------------------

tidy_celeste_outputs(output_dir = output_dir)

cw_inform("Done: Stage 7 (+ Stage 8 tidy tables) rerun for KEN/NIP/TUR.")
