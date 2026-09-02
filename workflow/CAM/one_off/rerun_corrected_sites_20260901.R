# rerun_corrected_sites_20260901.R
# =============================================================================
# One-off: applies Sam's manual review of the 2026-08-31 clean CAM streams
# rerun's 7 reverted lake-bisection sites (see [[cam-streams-clean-rerun]]
# memory once written):
#   - NCMN, SUD17, SUD200: leave as-is (Sam's call after manual QGIS review)
#   - Daisy, Wolf: pour_point_snapped.shp manually corrected by Sam in QGIS
#     (confirmed via mtime — both edited 2026-09-01, well after the clean
#     rerun's own delineation) -> rerun via edited_snap_site_ids, no lake
#     correction involved
#   - SUD102, SUD103, Tilton: genuine lake bisection -> rerun via
#     workflow/CAM/fix_lake_bisection.R's correct_lake_bisected_sites()
#
# Reconstructs only what rerun_engine_sites()/correct_lake_bisected_sites()
# actually need (Stage 1/1b config+group_manifest, Stage 7's LOI-layers
# setup) — skips re-running Stage 2-6 wholesale, since terrain/delineation/
# metrics for every OTHER site are already correct and cached; both rerun
# mechanisms cascade remove-upstream/reclip/metrics/hydroweight scoped to
# the affected sites (plus any cascaded neighbor) and MERGE into the
# existing combined CSVs on their own.
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

EXCLUDED_SITE_IDS <- c("SUD11")
sites_raw <- readr::read_csv(here("data/cam_stream_sites_raw.csv"), show_col_types = FALSE)
sites <- sites_raw |>
  dplyr::mutate(
    site_id = stream_id |> gsub("\\s+", "_", x = _) |> gsub("[^A-Za-z0-9_-]", "", x = _),
    site_name = stream_id
  ) |>
  dplyr::filter(!site_id %in% EXCLUDED_SITE_IDS) |>
  dplyr::select(site_id, site_name, lon, lat)

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

cw_inform(glue::glue("Reconstructed state: {nrow(sites)} site(s). Terrain/delineation NOT rerun — already cached from the 2026-08-31 clean rerun."))

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

prepare_cam_ndvi_site_rasters(sites, output_dir = output_dir, cache_dir = cache_dir)
prepare_cam_ndvi_trend_site_rasters(sites, output_dir = output_dir, cache_dir = cache_dir)
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
    path_lazy = fs::path(cache_dir, "hydroweight_loi", "harvest_regen", "{group_id}.tif"),
    name = "harvest_regen", type = "categorical", class_levels = harvest_regen_levels
  )
)

# =============================================================================
# 1. Daisy, Wolf — manually-edited snapped pour point, no lake correction
# =============================================================================

cw_inform("\n===== Rerunning Daisy, Wolf (manually-edited snap) =====")
rerun_engine_sites(
  edited_snap_site_ids = c("Daisy", "Wolf"),
  sites = sites, group_manifest = group_manifest, config = config,
  output_dir = output_dir, cache_dir = cache_dir, loi_layers = loi_layers
)

# =============================================================================
# 2. SUD102, SUD103, Tilton — genuine lake bisection, reactive fix
# =============================================================================

cw_inform("\n===== Correcting SUD102, SUD103, Tilton (lake bisection) =====")
source(here("workflow/CAM/fix_lake_bisection.R")) # runs validate_catchment_lake_intersections(), defines correct_lake_bisected_sites()

correct_lake_bisected_sites(
  site_ids = c("SUD102", "SUD103", "Tilton"),
  sites = sites, group_manifest = group_manifest, config = config,
  output_dir = output_dir, cache_dir = cache_dir, loi_layers = loi_layers
)

cw_inform("\n===== Done. NCMN/SUD17/SUD200 intentionally left as-is. =====")
