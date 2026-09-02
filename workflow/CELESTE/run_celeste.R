# run_celeste.R
# =============================================================================
# Top-level runner for CELESTE (all 6 HydroBasins groups — COC, KEN, MOR,
# NBE, NIP, TUR — 132 stream/river sites), on the modular, input-driven
# delineation engine (workflow/R/engine/) — the standard approach for every
# project in this repo now, CELESTE included. Supersedes an earlier,
# non-modular run_celeste.R that ran directly on workflow/R/stream/ (retired
# 2026-08 — see git history if you need it; workflow/R/stream/'s remaining
# files — group_sites.R, burn_streams.R, delineate_sites.R,
# hydroweight_attributes.R — are genuinely shared building blocks the
# engine reuses directly, not old-pipeline-specific, and stay in place).
#
# See workflow/CELESTE/README.md for the full rationale behind every
# CELESTE-specific decision below (site coordinates, COC's no-burn-in
# override, stage ordering, LOI cache-reuse rules, lake-bisection
# handling) and for how to reproduce this run from scratch.
#
# WORKFLOW STAGES:
#   Stage 1 — Resolve config, build group manifest (HydroBasins level-6
#              grouping, delegates to stream/group_sites.R unmodified)
#   Stage 2 — Prepare terrain (crop MRDEM, burn NHN streams in, breach,
#              upfront lake conditioning, D8 pointer, flow accumulation,
#              extract streams) — once per group, the expensive part
#   Stage 3 — Delineate catchments (point pour point, Jenson snap)
#   Stage 4 — Remove upstream nested catchments
#   Stage 5 — Re-clip rasters/flowlines to clipped catchments (REQUIRED
#              before hydroweighting — see README)
#   Stage 6 — Catchment morphometric metrics
#   Stage 7 — LOI prep (ndvi, ndvi_trend, harvest_regen) + distance-weighted
#              hydroweighting
#   Stage 8 — Tidy (long-format) plotting-ready tables
#
# Usage: source from an R session in the project root. Genuinely long —
# 6 groups (MOR's DEM crop is ~5x the smallest group's), ~90 min total on
# a cold cache (terrain prep + harvest_regen rasterization + hydroweighting
# dominate; ndvi/ndvi_trend reuse is near-instant).
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
source(here("workflow/R/stream/group_sites.R")) # build_group_manifest() — reused unmodified by the engine's "hydrobasins" strategy
source(here("workflow/R/stream/burn_streams.R")) # burn_streams_into_dem() — reused unmodified by resolve_streams_burn()
source(here("workflow/R/stream/delineate_sites.R")) # snap_pour_point(), delineate_watershed(), etc. — reused unmodified by engine/04
source(here("workflow/R/stream/hydroweight_attributes.R")) # calculate_hydroweight_attributes_stream() — shared with CAM, reused unmodified
source(here("workflow/R/remove_upstream.R"))
source(here("workflow/R/reclip_outputs.R")) # reclip_outputs() — required before Stage 7, see STAGE ORDER note above
source(here("workflow/R/catchment_metrics.R"))
source(here("workflow/R/engine/00_resolve_config.R"))
source(here("workflow/R/engine/01_build_group_manifest.R"))
source(here("workflow/R/engine/02_prepare_terrain.R"))
source(here("workflow/R/engine/03_prepare_streams_burn.R"))
source(here("workflow/R/engine/prepare_lake_conditioning.R")) # resolve_lake_conditioning() -- opt-in via config$lake_conditioning, off by default
source(here("workflow/R/engine/04_delineate_site.R"))
source(here("workflow/R/engine/99_rerun_sites")) # drop_redundant_clipped_rows(); also rerun_engine_site_watershed()/rerun_engine_sites() for post-hoc pour-point corrections
source(here("workflow/raster_attributes.R")) # sens_slope_trend(), rasterize_competing_classes(), mask_out_waterbodies() — used by the CELESTE prepare_*.R scripts below
source(here("workflow/CELESTE/prepare_ndvi.R"))
source(here("workflow/CELESTE/prepare_ndvi_trend.R"))
source(here("workflow/CELESTE/prepare_ndvi_masked.R")) # lake-masked "ndvi_masked"/"ndvi_trend_masked" LOIs, alongside (not in place of) the above
source(here("workflow/CELESTE/prepare_harvest_regen.R"))
source(here("workflow/CELESTE/tidy_outputs.R")) # tidy_celeste_outputs() — Stage 8

# =============================================================================
# CONFIGURATION
# =============================================================================

PROJECT_ID <- "CELESTE" # output/cache directory name — see README.md's
# "migration provenance" note for the retired pre-engine pipeline this
# consolidated with.
output_dir <- here("output", PROJECT_ID)
cache_dir  <- here("cache", PROJECT_ID)
fs::dir_create(output_dir, recurse = TRUE)
fs::dir_create(cache_dir, recurse = TRUE)

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
  mutate(burn_streams = dplyr::if_else(group_id == "COC", FALSE, burn_streams)) # see README.md's "COC runs without stream burn-in" note

cw_inform(glue::glue(
  "CELESTE engine run: {nrow(sites)} site(s) across ",
  "{length(unique(sites$group_id))} group(s): ",
  "{paste(sort(unique(sites$group_id)), collapse = ', ')}."
))

run_config <- list(
  project_id = PROJECT_ID,
  output_dir = output_dir,
  cache_dir  = cache_dir,
  sites      = sites,

  dem = list(path = mrdem_vrt),
  flow_direction = NULL,
  flow_pointer   = NULL,
  crs = NULL, # NULL = MRDEM's own native CRS (EPSG:3979)
  stream_threshold = 1000,

  streams_burn   = list(source = "nhn_auto"), # per group, from that group's own AOI — see README.md's "Stream burn-in" note
  nhn_index_path = nhn_index,
  nhn_raw_dir    = nhn_dir,

  # Opt-in upfront lake conditioning (workflow/R/engine/prepare_lake_
  # conditioning.R) — see README.md's "Upfront lake conditioning" note for
  # the verification results and why fix_lake_bisection.R remains as a
  # fallback for any residual.
  lake_conditioning = list(source = "nhn_auto"),

  lake_polygons = NULL,
  lake_buffer_m = 30, # unused in point pour-point mode

  grouping = list(
    strategy = "hydrobasins",
    hydrobasins_dir = hydrobasins_dir,
    default_buffer_m = 1000
  ),

  loi_layers = NULL # set in Stage 7 below
)

config <- resolve_engine_config(run_config)

# =============================================================================
# STAGE 1 — Group manifest
# =============================================================================

gm <- build_engine_group_manifest(config)
sites <- gm$sites
group_manifest <- gm$group_manifest
print(group_manifest)

# =============================================================================
# STAGE 2 — Terrain prep (crop, burn, breach, D8, accumulation, streams)
# =============================================================================

prepare_engine_terrain(config, group_manifest)

burn_status <- purrr::map_dfr(group_manifest$group_id, function(grp) {
  burned_path <- fs::path(cache_dir, grp, "dem_burned.tif")
  breached_path <- fs::path(cache_dir, grp, "dem_breached.tif")
  tibble(
    group_id = grp,
    dem_burned_written = fs::file_exists(burned_path),
    dem_breached_written = fs::file_exists(breached_path)
  )
})
cw_inform("\n--- Burn-in status per group (Stage 2 summary) ---")
print(burn_status)
readr::write_csv(burn_status, fs::path(output_dir, "burn_in_status.csv"))

# =============================================================================
# STAGE 3 — Delineate catchments
# =============================================================================

results <- delineate_engine_catchments(
  config = config, sites = sites, group_manifest = group_manifest,
  snap_dist = 200, min_cells = 10
)
print(results)

flagged <- dplyr::filter(results, flagged)
if (nrow(flagged) > 0) {
  cw_warn(glue::glue("\n{nrow(flagged)} site(s) flagged — review pour points:"))
  print(flagged[, c("site_id", "catchment_cells", "catchment_km2", "flag_reason")])
}

all_catchments <- purrr::map(sites$site_id, function(sid) {
  p <- fs::path(site_output_dir(output_dir, sid), "catchment.gpkg")
  if (!cache_exists(p)) {
    return(NULL)
  }
  sf::st_read(p, quiet = TRUE)
}) |>
  purrr::compact() |>
  dplyr::bind_rows()

sf::st_write(all_catchments, fs::path(output_dir, "all_catchments.gpkg"), delete_dsn = TRUE, quiet = TRUE)

cw_inform(glue::glue(
  "\nStage 3 complete. {nrow(results)} site(s) processed, ",
  "{nrow(all_catchments)} catchment(s) written to {output_dir}/all_catchments.gpkg."
))

# =============================================================================
# STAGE 4 — Remove upstream nested catchments
# =============================================================================

upstream_results <- remove_upstream_catchments(sites, output_dir)
print(upstream_results)

# =============================================================================
# STAGE 5 — Re-clip rasters and flowlines to clipped catchments
# =============================================================================
# REQUIRED before Stage 7 — see STAGE ORDER note above.

reclip_results <- reclip_outputs(sites = sites, output_dir = output_dir, group_manifest = group_manifest)
print(table(reclip_results$status))

# =============================================================================
# STAGE 6 — Catchment morphometric metrics
# =============================================================================

metrics <- calculate_catchment_metrics(sites = sites, output_dir = output_dir)
metrics <- drop_redundant_clipped_rows(metrics, upstream_results, site_col = "site_id")
ref_table <- build_metrics_reference_table()
write_metrics_outputs(metrics = metrics, ref_table = ref_table, output_dir = output_dir)
print(metrics)

# =============================================================================
# STAGE 7 — LOI prep + distance-weighted hydroweighting
# =============================================================================
# See README.md's "LOI layers" section for why ndvi/ndvi_trend reuse
# across reruns but harvest_regen is always recomputed fresh.

cw_inform("Preparing ndvi (cached under cache_dir — no DEM dependency)...")
prepare_ndvi_per_group_rasters(group_manifest, cache_dir = cache_dir)

cw_inform("Preparing ndvi_trend (cached under cache_dir — no DEM dependency)...")
prepare_ndvi_trend_rasters(group_manifest, cache_dir = cache_dir)

cw_inform("Preparing ndvi_masked (lake-masked, alongside ndvi — see prepare_ndvi_masked.R's header)...")
prepare_ndvi_masked_group_rasters(
  group_manifest, cache_dir = cache_dir,
  nhn_index_path = nhn_index, nhn_raw_dir = nhn_dir
)

cw_inform("Preparing ndvi_trend_masked (lake-masked, alongside ndvi_trend)...")
prepare_ndvi_trend_masked_rasters(
  group_manifest, cache_dir = cache_dir,
  nhn_index_path = nhn_index, nhn_raw_dir = nhn_dir
)

cw_inform("Preparing harvest_regen (FRESH — depends on this run's own dem_breached.tif)...")
prepare_harvest_regen_rasters(group_manifest, cache_dir = cache_dir)

# canlcc: same national land-cover product + class scheme as both CAM
# projects (CAN_LLC_2020.tif) — declared path_lazy, same as CAM streams'
# canlcc entry, so each site's LOI is cropped straight from the source
# once per site rather than reprojecting/caching the whole ~2.9 GB raster
# up front. No project-specific prep script needed (a single-year,
# ready-to-use raster — nothing to mosaic or rasterize). CAN_LLC_2020.tif
# is already in EPSG:3979, matching CELESTE's own working CRS exactly, so
# this needs no reprojection at all (unlike CAM's EPSG:3161 runs).
CANLCC_PATH <- "/Users/sam/Documents/cfs/shared_data/raw/landcover/CAN_LLC_2020.tif"

canlcc_levels <- data.frame(
  stringsAsFactors = FALSE,
  ID = c(1L, 2L, 5L, 6L, 8L, 10L, 11L, 12L, 13L, 14L, 15L, 16L, 17L, 18L, 19L),
  Class = c(
    "needleleaf_forest", "taiga_needleleaf_forest", "broadleaf_deciduous_forest",
    "mixed_forest", "temperate_shrubland", "temperate_grassland",
    "shrubland_lichen_moss", "grassland_lichen_moss", "barren_lichen_moss",
    "wetland", "cropland", "barren", "urban", "water", "snow_ice"
  )
)

loi_layers <- list(
  list(path_lazy = CANLCC_PATH, name = "canlcc", type = "categorical", class_levels = canlcc_levels),
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

hw_results <- calculate_hydroweight_attributes_stream(
  sites = sites,
  group_manifest = group_manifest,
  output_dir = output_dir,
  cache_dir = cache_dir,
  loi_layers = loi_layers,
  catchment_versions = c("unclipped", "clipped")
)

hw_results <- drop_redundant_clipped_rows(hw_results, upstream_results, site_col = "site")

write_csv(hw_results, fs::path(output_dir, "CELESTE_hydroweight.csv"))
print(hw_results)

cw_inform(glue::glue(
  "\nHydroweighting complete: {nrow(hw_results)} row(s) written to ",
  "{output_dir}/CELESTE_hydroweight.csv."
))

# =============================================================================
# STAGE 8 — Tidy (long-format) plotting-ready tables
# =============================================================================
# See workflow/CELESTE/tidy_outputs.R's header for the exact output columns/
# shapes and per-LOI year-handling rules.

tidy_celeste_outputs(output_dir = output_dir)
