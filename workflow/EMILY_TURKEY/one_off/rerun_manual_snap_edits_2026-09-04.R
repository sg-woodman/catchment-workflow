# rerun_manual_snap_edits_2026-09-04.R
# =============================================================================
# One-off record: reruns 3 EMILY_TURKEY sites whose pour_point_snapped.shp
# was manually edited in QGIS (mtimes confirmed well after the original
# snap: HARM027, GAR012, NORSHOR002 — all edited ~14:43-14:46, originally
# snapped ~10:27-10:56). Uses rerun_engine_sites()'s edited_snap_site_ids
# path (workflow/R/engine/99_rerun_sites) — reads the edited .shp directly,
# no re-snap — then cascades remove-upstream (full sites list),
# reclip/metrics/hydroweight (affected sites + any cascaded neighbor),
# merging into the existing combined CSVs rather than overwriting them.
#
# Rebuilds config/sites/group_manifest exactly as run_emily_turkey.R's
# Stage 1 does (all cache hits — terrain/most sites already delineated, so
# this is fast) rather than sourcing the whole runner, to avoid redoing
# Stage 2-7 work for the other 174 unaffected sites.
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

# =============================================================================
# Rebuild config/sites/group_manifest — identical to run_emily_turkey.R's
# Stage 1 (same paths, same site derivation, same manual_groups config)
# =============================================================================

PROJECT_ID <- "EMILY_TURKEY"
output_dir <- here("output", PROJECT_ID)
cache_dir  <- here("cache", PROJECT_ID)

OIH_NE_DEM_PATH      <- "/Users/sam/Downloads/IntegratedHydrologyNE/EnforcedDEM.tif"
OIH_NE_FLOW_DIR_PATH <- "/Users/sam/Downloads/IntegratedHydrologyNE/EnhancedFlowDirection.tif"
OIH_NC_DEM_PATH      <- "/Users/sam/Downloads/IntegratedHydrologyNC/EnforcedDEM.tif"
OIH_NC_FLOW_DIR_PATH <- "/Users/sam/Downloads/IntegratedHydrologyNC/EnhancedFlowDirection.tif"

oih_recode_matrix <- matrix(
  c(128, 1, 1, 2, 2, 4, 4, 8, 8, 16, 16, 32, 32, 64, 64, 128),
  ncol = 2, byrow = TRUE
)

sites_raw <- readr::read_csv(here("data/emily_turkey_sites_raw.csv"), show_col_types = FALSE)
sites <- sites_raw |>
  dplyr::mutate(
    site_id = sample_site |>
      gsub("'", "", x = _) |>
      gsub("\\s+", "_", x = _) |>
      gsub("[^A-Za-z0-9_-]", "", x = _),
    site_id = dplyr::if_else(
      sample_site == "HARM030'" & sample_code == 149,
      paste0(site_id, "_", sample_code),
      site_id
    ),
    site_name = sample_site
  ) |>
  dplyr::select(site_id, site_name, lon, lat, location_project, date_collected)
stopifnot(!anyDuplicated(sites$site_id))

run_config <- list(
  project_id = "EMILY_TURKEY",
  output_dir = output_dir,
  cache_dir  = cache_dir,
  sites      = sites,
  dem = NULL, flow_direction = NULL, flow_pointer = NULL, flow_accum = NULL,
  crs = NULL, stream_threshold = 100,
  streams_burn = list(source = "none"),
  lake_polygons = NULL, lake_buffer_m = 30,
  grouping = list(
    strategy = "manual_groups",
    groups = list(
      list(group_id = "NE", dem = list(path = OIH_NE_DEM_PATH),
           flow_direction = list(path = OIH_NE_FLOW_DIR_PATH, recode = oih_recode_matrix)),
      list(group_id = "NC", dem = list(path = OIH_NC_DEM_PATH),
           flow_direction = list(path = OIH_NC_FLOW_DIR_PATH, recode = oih_recode_matrix))
    )
  ),
  loi_layers = NULL
)

config <- resolve_engine_config(run_config)
gm <- build_engine_group_manifest(config)
sites <- gm$sites
group_manifest <- gm$group_manifest

# =============================================================================
# Rerun the 3 manually-edited sites
# =============================================================================

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
  list(path_lazy = CANLCC_PATH, name = "canlcc", type = "categorical", class_levels = canlcc_levels)
)

affected <- rerun_engine_sites(
  edited_snap_site_ids = c("HARM027", "GAR012", "NORSHOR002"),
  sites = sites, group_manifest = group_manifest, config = config,
  output_dir = output_dir, cache_dir = cache_dir, loi_layers = loi_layers
)
print(affected)
