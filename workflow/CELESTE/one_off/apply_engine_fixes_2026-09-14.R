# apply_engine_fixes_2026-09-14.R
# =============================================================================
# One-off: applies two shared-engine bug fixes (found and fixed the same
# day while investigating CAM streams, then checked directly against
# CELESTE's own output rather than assumed) to CELESTE via
# rerun_engine_sites().
#
# Found by two read-only diagnostic passes over every CELESTE site's
# output before touching anything:
#   1. 13 sites had an invalid (self-intersecting-ring) catchment.gpkg —
#      MOR-CRANE, NBE4/5/6, NBI1/4/5, WF_30K_INTLO, WK6DW/6UP,
#      WK_100K_INTHI, WK_30K_INTLO, WK_50K_INTLO. None had actually lost
#      their pour point (0/132 containment failures — CELESTE dodged that
#      part), but the geometry was still malformed, risking subtly wrong
#      area_km2/perim_km/etc. in catchment_metrics.csv (no validation
#      there) and only band-aided, not fixed, at hydroweight-read-time.
#   2. 19 of CELESTE's then-68 "clipped" rows were boundary-touch false
#      positives — a "nested" neighbor that only touched the focal
#      catchment at a shared flow divide, with real interior overlap
#      confirmed ~0, still counted by the old st_intersects()-only nested
#      check. Confirmed directly by symmetric-difference: those 19 sites'
#      catchment_clipped.gpkg was geometrically byte-identical to
#      catchment.gpkg.
#
# resnap_site_ids = the 13 invalid-geometry sites. Re-delineating them
# regenerates catchment.gpkg via the now-fixed watershed_to_polygon()
# (terra::as.polygons(dissolve=TRUE), commit 6dc39cc) — cheap: reuses the
# already-cached group terrain; only the polygonization step actually
# changes, not the underlying watershed raster.
#
# rerun_engine_sites() ALSO reruns remove_upstream_catchments() for the
# FULL 132-site list unconditionally (see its own docstring) — this picks
# up the boundary-touch nested-overlap fix (commit 4010fb6) project-wide,
# not just for the 13 resnapped sites, and its merge_rows_into_csv()
# filter_fn drops every now-redundant "clipped" row from the FULL
# catchment_metrics.csv/CELESTE_hydroweight.csv (all 19 found, not just
# any overlapping with the 13).
#
# VERIFIED (re-ran both diagnostic passes after, not just trusted the
# rerun log): invalid catchment.gpkg 13 -> 0, pour-point containment
# stayed 132/132, "clipped" rows 68 -> 49 sites and the surviving 49 are
# an exact match to the sites independently confirmed as genuine
# erasures beforehand. rerun_engine_sites() also auto-detected 2
# cascaded neighbors (NBE1, WF_100K_INTHI — sites whose own clipped
# catchment depended on one of the 13 that changed shape) and reran
# those too. Stage 8 tidy tables regenerated to match.
#
# Confirmed NOT affected by the other two fixes committed the same day:
# the reclip_outputs.R always-NA bug (CELESTE's only no-burn-in group,
# COC, has 2 sites with no nested neighbor, so no clipped row was ever
# written for it) and the hydroweight all-NA-column crash (CELESTE's
# sites table does carry an eligible all-NA column, aoi_buffer_m, but
# every one of 132 sites was already present in CELESTE_hydroweight.csv
# before this script ran, so it never actually crashed in production —
# the site_id-only fix now protects it regardless, going forward).

library(sf); library(terra); library(whitebox); library(dplyr); library(tidyr)
library(purrr); library(readr); library(tibble); library(fs); library(cli)
library(glue); library(here)

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
source(here("code/reset_workflow.R"))
source(here("workflow/R/engine/99_rerun_sites"))
source(here("workflow/raster_attributes.R"))
source(here("workflow/CELESTE/prepare_ndvi.R"))
source(here("workflow/CELESTE/prepare_ndvi_trend.R"))
source(here("workflow/CELESTE/prepare_ndvi_masked.R"))
source(here("workflow/CELESTE/prepare_harvest_regen.R")) # defines harvest_regen_levels, used in loi_layers below
source(here("workflow/CELESTE/tidy_outputs.R"))

# -- Config (identical to run_celeste.R) -------------------------------------

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
  project_id = PROJECT_ID, output_dir = output_dir, cache_dir = cache_dir, sites = sites,
  dem = list(path = mrdem_vrt), flow_direction = NULL, flow_pointer = NULL,
  crs = NULL, stream_threshold = 1000,
  streams_burn = list(source = "nhn_auto"), nhn_index_path = nhn_index, nhn_raw_dir = nhn_dir,
  lake_conditioning = list(source = "nhn_auto"),
  lake_polygons = NULL, lake_buffer_m = 30,
  grouping = list(strategy = "hydrobasins", hydrobasins_dir = hydrobasins_dir, default_buffer_m = 1000),
  loi_layers = NULL
)
config <- resolve_engine_config(run_config)

gm <- build_engine_group_manifest(config)
sites <- gm$sites
group_manifest <- gm$group_manifest

cw_inform(glue::glue("Reconstructed state: {nrow(sites)} site(s). Terrain/delineation cache reused as-is."))

# -- loi_layers (identical to run_celeste.R Stage 7) --------------------------

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

# -- The fix -------------------------------------------------------------------

invalid_geom_sites <- c(
  "MOR-CRANE", "NBE4", "NBE5", "NBE6", "NBI1", "NBI4", "NBI5",
  "WF_30K_INTLO", "WK6DW", "WK6UP", "WK_100K_INTHI", "WK_30K_INTLO", "WK_50K_INTLO"
)

cw_inform("\n===== rerun_engine_sites(): regenerating 13 invalid-geometry catchments + full remove-upstream rerun =====")
affected <- rerun_engine_sites(
  resnap_site_ids = invalid_geom_sites,
  sites = sites, group_manifest = group_manifest, config = config,
  output_dir = output_dir, cache_dir = cache_dir, loi_layers = loi_layers
)
cw_inform(glue::glue("Affected sites (target + cascaded neighbors): {paste(affected, collapse = ', ')}"))

# -- Stage 8: tidy outputs (full regenerate, cheap) ---------------------------
cw_inform("\n===== Stage 8: tidy_celeste_outputs() =====")
tidy_celeste_outputs(output_dir = output_dir)

cw_inform("\nDone.")
