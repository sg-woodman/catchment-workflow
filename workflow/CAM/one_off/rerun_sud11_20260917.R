# rerun_sud11_20260917.R
# =============================================================================
# One-off: adds SUD11 to the CAM streams project after its source longitude
# was corrected (lon = -51.15433 -> -81.15433, matching neighbours SUD12/
# VER01 at ~-81.0) in both data/cam_stream_sites_raw.csv and
# data/cam_stream_sites_raw.gpkg. SUD11 was never delineated before (it was
# listed in EXCLUDED_SITE_IDS from initial delivery onward), so this is a
# resnap_site_ids rerun, not an edited-snap correction — reset_site() is a
# safe no-op on a site with no existing output dir.
#
# Reconstructs only what rerun_engine_sites() actually needs (Stage 1
# config+group_manifest, Stage 7's LOI-layers setup) — Stage 2 terrain is
# untouched (shared whole_domain group, unaffected by adding a site).
# rerun_engine_sites() itself handles Stage 3 delineation (resnap), then
# cascades remove-upstream (full sites list)/reclip/metrics/hydroweight
# scoped to SUD11 (+ any cascaded neighbor), merging into the existing
# combined CSVs rather than overwriting them wholesale.
#
# EXCLUDED_SITE_IDS in run_cam_streams.R has been updated to character(0)
# to match — see workflow/CAM/README.md's "Site list and the SUD11
# exclusion" note.
# =============================================================================

library(sf); library(terra); library(whitebox); library(dplyr); library(tidyr)
library(purrr); library(readr); library(tibble); library(fs); library(cli)
library(glue); library(here)

source(here("workflow/R/utils.R"))
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
source(here("workflow/R/engine/04_delineate_site.R"))
source(here("code/reset_workflow.R"))
source(here("workflow/R/engine/99_rerun_sites"))
source(here("workflow/gee_utils.R"))
source(here("workflow/raster_attributes.R"))
source(here("workflow/CAM/prepare_ndvi.R"))
source(here("workflow/CAM/prepare_ndvi_trend.R"))
source(here("workflow/CAM/prepare_ndvi_masked.R"))
source(here("workflow/CAM/prepare_harvest_regen.R"))

# -- Config (identical to run_cam_streams.R) ---------------------------------

PROJECT_ID  <- "CAM"
DELINEATION <- "stream_delineation"
output_dir <- here("output", PROJECT_ID, DELINEATION)
cache_dir  <- here("cache", PROJECT_ID, DELINEATION)

OIH_DEM_PATH       <- "/Users/sam/Downloads/IntegratedHydrologyNE/EnforcedDEM.tif"
OIH_FLOW_DIR_PATH  <- "/Users/sam/Downloads/IntegratedHydrologyNE/EnhancedFlowDirection.tif"
oih_recode_matrix <- matrix(
  c(128, 1, 1, 2, 2, 4, 4, 8, 8, 16, 16, 32, 32, 64, 64, 128),
  ncol = 2, byrow = TRUE
)
STREAM_THRESHOLD <- 100
SNAP_DIST_M      <- 200
MIN_CELLS        <- 10

EXCLUDED_SITE_IDS <- character(0) # SUD11's coordinate is now corrected — see run_cam_streams.R
sites_raw <- readr::read_csv(here("data/cam_stream_sites_raw.csv"), show_col_types = FALSE)
sites <- sites_raw |>
  dplyr::mutate(
    site_id = stream_id |> gsub("\\s+", "_", x = _) |> gsub("[^A-Za-z0-9_-]", "", x = _),
    site_name = stream_id
  ) |>
  dplyr::filter(!site_id %in% EXCLUDED_SITE_IDS) |>
  dplyr::select(site_id, site_name, lon, lat)

stopifnot("SUD11" %in% sites$site_id)
stopifnot(dplyr::filter(sites, site_id == "SUD11")$lon == -81.15433)

run_config <- list(
  project_id = "CAM_streams", output_dir = output_dir, cache_dir = cache_dir, sites = sites,
  dem = list(path = OIH_DEM_PATH),
  flow_direction = list(path = OIH_FLOW_DIR_PATH, recode = oih_recode_matrix),
  flow_pointer = NULL, flow_accum = NULL, crs = NULL, stream_threshold = STREAM_THRESHOLD,
  streams_burn = list(source = "none"),
  lake_polygons = NULL, lake_buffer_m = 30,
  grouping = list(strategy = "whole_domain"),
  loi_layers = NULL
)
config <- resolve_engine_config(run_config)

gm <- build_engine_group_manifest(config)
sites <- gm$sites
group_manifest <- gm$group_manifest

cw_inform(glue::glue("Reconstructed state: {nrow(sites)} site(s), including SUD11. Terrain NOT rerun — shared whole_domain group, unaffected by adding a site."))

# -- Stage 7 LOI-layers setup (identical to run_cam_streams.R) ---------------

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

# Each of these is cache-aware per site/group — with 43 sites already
# prepared, this only does real work for SUD11.
prepare_cam_ndvi_site_rasters(sites, output_dir = output_dir, cache_dir = cache_dir)
prepare_cam_ndvi_trend_site_rasters(sites, output_dir = output_dir, cache_dir = cache_dir)
prepare_cam_ndvi_masked_site_rasters(sites, output_dir = output_dir, cache_dir = cache_dir)
prepare_cam_ndvi_trend_masked_site_rasters(sites, output_dir = output_dir, cache_dir = cache_dir)
prepare_cam_harvest_regen_rasters(group_manifest, cache_dir = cache_dir)

loi_layers <- list(
  list(path_lazy = CANLCC_PATH, name = "canlcc", type = "categorical", class_levels = canlcc_levels),
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
  ),
  list(
    path_lazy = fs::path(cache_dir, "hydroweight_loi", "harvest_regen", "{group_id}.tif"),
    name = "harvest_regen", type = "categorical", class_levels = harvest_regen_levels
  )
)

# =============================================================================
# Add SUD11 — resnap from its now-corrected raw coordinate
# =============================================================================

cw_inform("\n===== Adding SUD11 (corrected longitude, never previously delineated) =====")
rerun_engine_sites(
  resnap_site_ids = "SUD11",
  sites = sites, group_manifest = group_manifest, config = config,
  output_dir = output_dir, cache_dir = cache_dir, loi_layers = loi_layers,
  snap_dist = SNAP_DIST_M, min_cells = MIN_CELLS
)

# =============================================================================
# Stage 8 — refresh plotting-ready long tables to include SUD11
# =============================================================================

source(here("workflow/CAM/tidy_outputs.R"))
tidy_cam_outputs(output_dir = output_dir)

cw_inform("\n===== Done. SUD11 added to CAM streams outputs. =====")
