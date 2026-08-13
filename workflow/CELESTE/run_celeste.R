# run_celeste.R
# =============================================================================
# Top-level runner for the CELESTE stream-site catchment delineation
# workflow, on the shared workflow/ layout (see workflow/R/stream/ and the
# modules shared with the lake pipeline: utils.R, 06_remove_upstream.R,
# 07_reclip_outputs.R, 08_catchment_metrics.R).
#
# Supersedes code/delineate_catchments.R + R/*.R for this project — verified
# via workflow/verify_stream_migration.R to reproduce the same catchments
# (see that script's header for the full methodology and the two documented,
# investigated exceptions).
#
# USAGE:
#   1. Verify the paths in the CONFIGURATION section below
#   2. Source this file or run sections interactively in R
#
# WORKFLOW STAGES:
#   Stage 1 — Validate sites and build group manifest
#   Stage 2 — Prepare DEMs (crop MRDEM to group AOIs)
#   Stage 3 — Prepare NHN layers and burn streams into DEM
#   Stage 4 — Run WhiteboxTools conditioning (breach, pointer, accumulation,
#              stream extraction, hillshade)
#   Stage 5 — Delineate catchments for all sites
#   Stage 6 — Remove upstream nested catchments (clipping)
#   Stage 7 — Re-clip rasters and flowlines to clipped catchments
#   Stage 8 — Catchment morphometric metrics (clipped + unclipped)
#   Stage 9 — Distance-weighted catchment attributes (hydroweight), for both
#              catchment versions, using the pour point as the O-scheme
#              target (see workflow/R/stream/06_hydroweight_attributes.R)
#
# CACHING:
#   Group-level rasters (DEM, breached DEM, flow products) are cached in
#   cache/CELESTE/. Re-running skips any step whose output already exists.
#   Use code/reset_workflow.R to clear cache selectively or fully.
#
# CORRECTING POUR POINTS:
#   If a catchment looks wrong, inspect pour_point_snapped.shp and
#   streams_tmp.tif in QGIS for that site. Edit pour_point_snapped.shp
#   to move the point to the correct stream cell, then:
#     source(here("workflow/R/stream/99_rerun_sites"))
#     rerun_site_watershed("site_id", sites, group_manifest, output_dir)
#   This reads the edited snapped pour point directly, without rerunning
#   snapping or the group-level stages.
# =============================================================================

# -- Packages ------------------------------------------------------------------
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

# -- Source modules -------------------------------------------------------------
# Core, standardized modules (shared across every project on this workflow/
# layout — stream, lake, and any future project) live under workflow/R/ and
# workflow/raster_attributes.R. CELESTE-specific data prep (this runner and
# its LOI prepare_*.R scripts) lives under workflow/CELESTE/ instead, kept
# physically separate so another project (e.g. CAM lakes) can have its own
# prepare_*.R scripts without the two ever colliding or needing to agree on
# project-specific details.
source(here("workflow/R/utils.R"))
source(here("workflow/R/stream/01_group_sites.R"))
source(here("workflow/R/stream/02_prepare_dem.R"))
source(here("workflow/R/stream/03_burn_streams.R"))
source(here("workflow/R/stream/04_run_whitebox.R"))
source(here("workflow/R/stream/05_delineate_sites.R"))
source(here("workflow/R/stream/99_rerun_sites"))
source(here("workflow/R/stream/06_hydroweight_attributes.R"))
source(here("workflow/R/06_remove_upstream.R"))
source(here("workflow/R/07_reclip_outputs.R"))
source(here("workflow/R/08_catchment_metrics.R"))
source(here("code/reset_workflow.R"))
source(here("workflow/raster_attributes.R")) # reclassify_categorical(), sens_slope_trend(), rasterize_competing_classes() — used ad hoc to prepare LOI rasters below, not part of Stages 1-9
source(here("workflow/CELESTE/prepare_ndvi.R")) # prepare_ndvi_vrt() — NDVI (continuous) LOI prep, see Stage 9
source(here("workflow/CELESTE/prepare_ndvi_trend.R")) # prepare_ndvi_trend_rasters() — NDVI Sen's-slope trend (continuous) LOI prep, see Stage 9
source(here("workflow/CELESTE/prepare_harvest_regen.R")) # prepare_harvest_regen_rasters() — Ontario/NB harvest/regen (categorical) LOI prep, see Stage 9

# =============================================================================
# CONFIGURATION
# =============================================================================

PROJECT_ID <- "CELESTE"

output_dir <- here("output", PROJECT_ID)
cache_dir <- here("cache", PROJECT_ID)
dir_create(output_dir)
dir_create(cache_dir)

# Sites — data/celeste_milli_sites_clean.gpkg (132 sites, 6 groups: COC, KEN,
# MOR, NBE, NIP, TUR), built once from the raw site export. To rebuild it
# from an updated source CSV, see the one-off snippet at the bottom of this
# file rather than editing this block.
sites_sf <- st_read(here("data/celeste_milli_sites_clean.gpkg"), quiet = TRUE)
coords <- st_coordinates(sites_sf)
sites_csv <- sites_sf |>
  st_drop_geometry() |>
  as_tibble() |>
  mutate(lon = coords[, "X"], lat = coords[, "Y"]) |>
  select(site_id, site_name, lon, lat, group_id, burn_streams, aoi_buffer_m)

# Path to MRDEM .vrt file
mrdem_vrt <- "~/Documents/cfs/shared_data/raw/dem/mrdem-30-dtm.vrt"

# Path to NHN root directory (containing nhn_rhn_*_gdb_en subfolders)
nhn_dir <- "/Users/sam/Documents/cfs/shared_data/raw/hydro/networks/NHN/gdb"

# Path to NHN index shapefile
nhn_index <- "/Users/sam/Documents/cfs/shared_data/raw/hydro/networks/NHN/NHN_INDEX_WORKUNIT_LIMIT_2/NHN_INDEX_22_INDEX_WORKUNIT_LIMIT_2.shp"

# Path to HydroBasins root directory (containing 'north_america' and 'arctic')
hydrobasins_dir <- "/Users/sam/Documents/cfs/shared_data/raw/hydro/watersheds/HydroBasins"

# Global parameters (can be overridden per group via sites tibble columns)
default_buffer_m <- 1000 # buffer around HydroBasins polygon (metres)
target_crs <- 3979 # output CRS — matches MRDEM native CRS
snap_dist <- 200 # pour point snap distance (metres)
default_stream_threshold <- 1000 # flow accumulation threshold for streams
max_dist <- 10 # max breach path length (cells)
min_cells <- 10 # minimum catchment size before flagging

# =============================================================================
# STAGE 1 — Validate sites and build group manifest
# =============================================================================

check_packages()

sites <- validate_sites_tibble(sites_csv)

group_manifest <- build_group_manifest(
  sites = sites,
  output_dir = output_dir,
  cache_dir = cache_dir,
  hydrobasins_dir = hydrobasins_dir,
  default_buffer_m = default_buffer_m
)

write_group_aois(group_manifest, cache_dir)

message(glue(
  "\nGroup manifest built. Inspect {cache_dir}/group_aois.gpkg in QGIS to ",
  "verify AOI extents before proceeding to Stage 2.\n"
))

# =============================================================================
# STAGE 2 — Prepare DEMs
# =============================================================================

prepare_dem(
  group_manifest = group_manifest,
  mrdem_vrt = mrdem_vrt,
  target_crs = target_crs
)

# =============================================================================
# STAGE 3 — NHN layers and stream burning
# =============================================================================

prepare_nhn_layers(
  group_manifest = group_manifest,
  nhn_dir = nhn_dir,
  nhn_index = nhn_index
)

# =============================================================================
# STAGE 4 — WhiteboxTools hydrological conditioning
# =============================================================================

run_whitebox(
  group_manifest = group_manifest,
  max_dist = max_dist,
  flat_increment = NULL,
  fill = TRUE,
  default_stream_threshold = default_stream_threshold
)

wb_check <- verify_whitebox_outputs(group_manifest)
print(wb_check)

# =============================================================================
# STAGE 5 — Delineate catchments for all sites
# =============================================================================

results <- delineate_catchments(
  sites = sites,
  group_manifest = group_manifest,
  output_dir = output_dir,
  snap_dist = snap_dist,
  min_cells = min_cells
)
print(results)

flagged <- dplyr::filter(results, flagged)
if (nrow(flagged) > 0) {
  message(glue("\n{nrow(flagged)} site(s) flagged — review pour points:"))
  print(flagged[, c(
    "site_id",
    "catchment_cells",
    "catchment_km2",
    "flag_reason"
  )])
}

# Combine all catchment polygons into a single file for QA in QGIS
all_catchments <- purrr::map(sites$site_id, function(sid) {
  p <- fs::path(site_output_dir(output_dir, sid), "catchment.gpkg")
  if (!cache_exists(p)) {
    return(NULL)
  }
  sf::st_read(p, quiet = TRUE)
}) |>
  purrr::compact() |>
  dplyr::bind_rows()

sf::st_write(
  all_catchments,
  fs::path(output_dir, "all_catchments.gpkg"),
  delete_dsn = TRUE,
  quiet = TRUE
)

message(glue(
  "\nStage 5 complete. {nrow(results)} site(s) processed. ",
  "Combined catchments saved to {output_dir}/all_catchments.gpkg"
))

# =============================================================================
# STAGE 6 — Remove upstream nested catchments
# =============================================================================

remove_upstream_catchments(sites, output_dir)

# =============================================================================
# STAGE 7 — Re-clip rasters and flowlines to clipped catchments
# =============================================================================

reclip_outputs(
  sites = sites,
  group_manifest = group_manifest,
  output_dir = output_dir
)

# =============================================================================
# STAGE 8 (optional) — Catchment morphometric metrics
# =============================================================================
# Comment out if metrics are not needed for this project.

metrics <- calculate_catchment_metrics(sites = sites, output_dir = output_dir)
ref_table <- build_metrics_reference_table()
write_metrics_outputs(
  metrics = metrics,
  ref_table = ref_table,
  output_dir = output_dir
)
print(metrics)

# =============================================================================
# STAGE 9 (optional) — Distance-weighted catchment attributes (hydroweight)
# =============================================================================
# Comment out if hydroweight attributes are not needed for this project.
#
# LOI (layer of interest) rasters — three ways to declare a layer,
# corresponding to the three fields loi_layers[[i]] can set (exactly one):
#   $path          — a single raster already covering the full site extent
#   $path_template — one pre-clipped raster per site, "{site_id}" placeholder
#                    (only needs to cover the UNCLIPPED catchment; the
#                    clipped-version pass crops it down further)
#   $path_lazy     — one shared source too large to reproject/cache whole
#                    (e.g. a mosaic VRT over scattered regional tiles) —
#                    cropped to each site's catchment before ever being
#                    reprojected or written to disk. See
#                    workflow/R/stream/06_hydroweight_attributes.R's
#                    resolve_site_loi_raster() for how this differs from
#                    $path.
#
# Each LOI's actual prep (source data, mosaicking/reprojection/
# rasterization, temporal precedence, etc.) lives in its own workflow-
# specific script, not here — that logic is exactly the part expected to
# differ between projects, so it stays out of both this standardized
# runner and the generic hydroweight module. This block just calls each
# one and gets back a path to feed into loi_layers below.
#
#   NDVI (continuous)                         -> workflow/CELESTE/prepare_ndvi.R
#   NDVI trend (continuous, slope + p_value)  -> workflow/CELESTE/prepare_ndvi_trend.R
#   Harvest/regen (categorical, Ontario + NB) -> workflow/CELESTE/prepare_harvest_regen.R
#
# To add another LOI, write a similar prepare_*.R script under
# workflow/CELESTE/ (reusing workflow/raster_attributes.R's generic helpers
# where they fit) and call it here the same way.

prepare_ndvi_per_group_rasters(group_manifest, cache_dir = cache_dir)
prepare_ndvi_trend_rasters(group_manifest, cache_dir = cache_dir) # slow — see that script's header
prepare_harvest_regen_rasters(group_manifest, cache_dir = cache_dir)

loi_layers <- list(
  list(
    path_lazy = here(
      "cache", PROJECT_ID, "hydroweight_loi", "harvest_regen", "{group_id}.tif"
    ),
    name = "harvest_regen",
    type = "categorical",
    class_levels = harvest_regen_levels
  ),
  list(
    path_lazy = here(
      "cache", PROJECT_ID, "hydroweight_loi", "ndvi", "{group_id}.tif"
    ),
    name = "ndvi",
    type = "continuous",
    # sum/cell_count/NA_cell_count are meaningless for this analysis (per
    # user request) — excluded from computation entirely, not just hidden
    # from output (see run_loi_attributes_stream_multilayer_continuous()'s
    # numeric_stats param).
    stats = c("distwtd_mean", "distwtd_sd", "mean", "sd", "median", "min", "max")
  ),
  list(
    path_lazy = here(
      "cache", PROJECT_ID, "hydroweight_loi", "ndvi_trend", "{group_id}.tif"
    ),
    name = "ndvi_trend",
    type = "continuous",
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

write_csv(
  hw_results,
  here("output", PROJECT_ID, paste0(PROJECT_ID, "_hydroweight.csv"))
)

print(hw_results)

# =============================================================================
# One-off: rebuilding data/celeste_milli_sites_clean.gpkg from a new export
# =============================================================================
# Run manually, not as part of the normal Stage 1-8 sequence above, when a
# fresh site list is provided.
#
# sites_raw <- read_csv("~/Downloads/celeste_site_locations_cleaned.csv") |>
#   distinct()
# sites_sf <- st_as_sf(sites_raw, coords = c("lon", "lat"), crs = 4326)
# plot(sites_sf)
# st_write(
#   sites_sf,
#   here("data/celeste_milli_sites_clean.gpkg"),
#   delete_dsn = TRUE
# )
