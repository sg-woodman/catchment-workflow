# run_celeste_engine_hydroweight.R
# =============================================================================
# LOI prep + hydroweighting for the engine-based CELESTE run
# (output/CELESTE_engine / cache/CELESTE_engine — see run_celeste_engine.R).
# Reuses CELESTE's existing prepare_ndvi_per_group_rasters()/
# prepare_ndvi_trend_rasters()/prepare_harvest_regen_rasters() (workflow/
# CELESTE/prepare_*.R) completely unmodified — no engine-specific LOI code
# needed, since all three already take group_manifest + cache_dir as plain
# arguments, decoupled from which pipeline (stream vs. engine) produced the
# group_manifest.
#
# CACHE REUSE, PER-LOI (deliberate, not uniform):
#   - "ndvi" and "ndvi_trend": their per-group raster depends ONLY on
#     group_id (matched against data/ndvi/*.tif filenames) — NOT on the
#     group's DEM/breach at all (confirmed by reading both functions: the
#     only per-group state they touch is group_id itself, cache_dir is
#     just where to write). REUSED DIRECTLY from cache/CELESTE (the
#     production run) — cache_dir = here("cache/CELESTE") is passed for
#     these two, so cache_exists() finds every group's raster already
#     computed (confirmed present for all 6 groups) and skips instantly.
#     This is not an approximation — the correct result is identical
#     regardless of which pipeline delineated the sites, since neither
#     function reads anything terrain-derived.
#   - "harvest_regen": genuinely DEM-dependent — prepare_harvest_regen_
#     rasters() rasterizes onto a template read from
#     cache_dir/<group_id>/dem_breached.tif. Computed FRESH here, pointed
#     at cache/CELESTE_engine (this run's own breach grid), since the
#     engine's dem_breached.tif is a real, independently-computed raster,
#     not guaranteed to share an identical grid with production's (already
#     confirmed materially different for MOR, whose burn-in coverage
#     roughly doubled after the site-coordinate correction). MOR is
#     expected to stay excluded regardless (pre-existing, documented gap —
#     zero real coverage from either harvest/regen source, confirmed
#     against actual gdb features, not a caching artifact).
#
# STAGE ORDER: remove_upstream_catchments() AND reclip_outputs() must both
# run first. remove_upstream_catchments() produces catchment_clipped.gpkg,
# needed for the "clipped" hydroweight version at all. reclip_outputs() is
# ALSO required, not optional — confirmed the hard way (2026-08-26): the
# "clipped" version's stream-based weighting schemes (iEucS/iFLS/HAiFLS)
# silently come out missing (not erroring — process_hw_site_stream()'s
# target_S resolution falls back to dropping stream-based schemes when
# streams_clipped.tif isn't found, per resolve_hw_streams_path_stream())
# unless reclip_outputs() has already run to produce that file. An earlier
# version of this script skipped reclip_outputs() on the reasoning that
# "hydroweighting reads catchment polygons + group-level rasters directly,
# it doesn't need per-site re-clipped rasters" — true for the DEM/flow_
# accum rasters, but NOT true for the streams raster specifically, which
# has separate clipped/unclipped file variants
# (resolve_hw_streams_path_stream() looks for streams_clipped.tif when
# version == "clipped", streams.tif otherwise) — only reclip_outputs()
# produces the former. calculate_catchment_metrics() (production's Stage
# 8) is still not run here — hydroweighting doesn't need it.
#
# Usage: source after run_celeste_engine.R (and its COC-burn-in fix) have
# already been run — this script assumes output/CELESTE_engine already has
# a catchment.gpkg for every site.
# =============================================================================

library(sf); library(terra); library(whitebox); library(dplyr); library(tidyr)
library(purrr); library(readr); library(tibble); library(fs); library(cli); library(glue); library(here)

source(here("workflow/R/utils.R"))
source(here("workflow/R/stream/01_group_sites.R"))
source(here("workflow/R/stream/03_burn_streams.R"))
source(here("workflow/R/stream/05_delineate_sites.R"))
source(here("workflow/R/stream/06_hydroweight_attributes.R")) # calculate_hydroweight_attributes_stream()
source(here("workflow/R/06_remove_upstream.R"))
source(here("workflow/R/07_reclip_outputs.R")) # reclip_outputs() — see STAGE ORDER note above for why this is required, not optional
source(here("workflow/R/08_catchment_metrics.R")) # calculate_catchment_metrics() — standard stage in every engine project (workflow/templates/run_engine_template.R's own Stage 6), missing here only because this script started as an ad-hoc validation script, not the production sequence
source(here("workflow/R/engine/00_resolve_config.R"))
source(here("workflow/R/engine/01_build_group_manifest.R"))
source(here("workflow/R/engine/02_prepare_terrain.R"))
source(here("workflow/R/engine/03_prepare_streams_burn.R"))
source(here("workflow/R/engine/04_delineate_site.R"))
source(here("workflow/R/engine/99_rerun_sites")) # drop_redundant_clipped_rows()
source(here("workflow/raster_attributes.R")) # sens_slope_trend(), rasterize_competing_classes()
source(here("workflow/CELESTE/prepare_ndvi.R"))
source(here("workflow/CELESTE/prepare_ndvi_trend.R"))
source(here("workflow/CELESTE/prepare_harvest_regen.R"))

# =============================================================================
# Rebuild config/group_manifest exactly as run_celeste_engine.R does
# (corrected site coordinates, COC no-burn-in override)
# =============================================================================

PROJECT_ID <- "CELESTE_engine"
output_dir <- here("output", PROJECT_ID)
cache_dir  <- here("cache", PROJECT_ID)
production_cache_dir <- here("cache/CELESTE") # for ndvi/ndvi_trend reuse only — see file header

mrdem_vrt       <- "~/Documents/cfs/shared_data/raw/dem/mrdem-30-dtm.vrt"
hydrobasins_dir <- "/Users/sam/Documents/cfs/shared_data/raw/hydro/watersheds/HydroBasins"
nhn_dir         <- "/Users/sam/Documents/cfs/shared_data/raw/hydro/networks/NHN/gdb"
nhn_index       <- "/Users/sam/Documents/cfs/shared_data/raw/hydro/networks/NHN/NHN_INDEX_WORKUNIT_LIMIT_2/NHN_INDEX_22_INDEX_WORKUNIT_LIMIT_2.shp"

sites_sf <- st_read(here("data/celeste_milli_sites_clean_corrected.gpkg"), quiet = TRUE) # includes the 6 manually-corrected pour points
coords <- st_coordinates(sites_sf)
sites <- sites_sf |>
  st_drop_geometry() |> as_tibble() |>
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
sites <- gm$sites
group_manifest <- gm$group_manifest

# =============================================================================
# Prerequisite: remove-upstream (produces catchment_clipped.gpkg per site)
# =============================================================================

upstream_results <- remove_upstream_catchments(sites, output_dir)
print(upstream_results)

# =============================================================================
# Prerequisite: reclip_outputs (produces streams_clipped.tif + other
# *_clipped.tif per site — see STAGE ORDER note above)
# =============================================================================

reclip_results <- reclip_outputs(sites = sites, output_dir = output_dir, group_manifest = group_manifest)
print(table(reclip_results$status))

# =============================================================================
# Catchment morphometric metrics — standard stage for every catchment
# delineation workflow in this project (CAM already has this; CELESTE_engine
# didn't, only because this script predates that being made explicit)
# =============================================================================

metrics <- calculate_catchment_metrics(sites = sites, output_dir = output_dir)
metrics <- drop_redundant_clipped_rows(metrics, upstream_results, site_col = "site_id")
ref_table <- build_metrics_reference_table()
write_metrics_outputs(metrics = metrics, ref_table = ref_table, output_dir = output_dir)
print(metrics)

# =============================================================================
# LOI prep
# =============================================================================

cw_inform("Preparing ndvi (reusing production cache — no DEM dependency)...")
prepare_ndvi_per_group_rasters(group_manifest, cache_dir = production_cache_dir)

cw_inform("Preparing ndvi_trend (reusing production cache — no DEM dependency)...")
prepare_ndvi_trend_rasters(group_manifest, cache_dir = production_cache_dir)

cw_inform("Preparing harvest_regen (FRESH — depends on this run's own dem_breached.tif)...")
prepare_harvest_regen_rasters(group_manifest, cache_dir = cache_dir)

loi_layers <- list(
  list(
    path_lazy = here("cache/CELESTE_engine", "hydroweight_loi", "harvest_regen", "{group_id}.tif"),
    name = "harvest_regen", type = "categorical", class_levels = harvest_regen_levels
  ),
  list(
    path_lazy = here("cache/CELESTE", "hydroweight_loi", "ndvi", "{group_id}.tif"), # production cache — see header
    name = "ndvi", type = "continuous",
    stats = c("distwtd_mean", "distwtd_sd", "mean", "sd", "median", "min", "max")
  ),
  list(
    path_lazy = here("cache/CELESTE", "hydroweight_loi", "ndvi_trend", "{group_id}.tif"), # production cache — see header
    name = "ndvi_trend", type = "continuous",
    stats = c("distwtd_mean", "distwtd_sd", "mean", "sd", "median", "min", "max")
  )
)

# =============================================================================
# Hydroweighting
# =============================================================================

hw_results <- calculate_hydroweight_attributes_stream(
  sites = sites,
  group_manifest = group_manifest,
  output_dir = output_dir,
  cache_dir = cache_dir,
  loi_layers = loi_layers,
  catchment_versions = c("unclipped", "clipped")
)

hw_results <- drop_redundant_clipped_rows(hw_results, upstream_results, site_col = "site")

write_csv(hw_results, fs::path(output_dir, "CELESTE_engine_hydroweight.csv"))
print(hw_results)

cw_inform(glue::glue(
  "\nHydroweighting complete: {nrow(hw_results)} row(s) written to ",
  "{output_dir}/CELESTE_engine_hydroweight.csv."
))
