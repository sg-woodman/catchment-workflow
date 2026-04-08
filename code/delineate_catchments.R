# delineate_catchments.R
# =============================================================================
# Top-level runner for the catchment delineation workflow.
#
# USAGE:
#   1. Define your sites using sites_template.R (or sites_template.csv)
#   2. Fill in the paths in the CONFIGURATION section below
#   3. Source this file or run sections interactively
#
# WORKFLOW STAGES:
#   Stage 1 — Validate sites and build group manifest
#   Stage 2 — Prepare DEMs (crop MRDEM to group AOIs)
#   Stage 3 — Prepare NHN layers and burn streams into DEM
#   Stage 4 — Run WhiteboxTools conditioning (breach, pointer, accumulation,
#              stream extraction, hillshade)
#   Stage 5 — Delineate catchments for all sites
#
# CACHING:
#   Group-level rasters (DEM, breached DEM, flow products) are cached in
#   cache/<group_id>/. Re-running skips any step whose output already exists.
#   Use reset_workflow.R to clear cache selectively or fully.
#
# CORRECTING POUR POINTS:
#   If a catchment looks wrong, inspect pour_point_snapped.shp and
#   streams_tmp.tif in QGIS for that site. Edit pour_point_snapped.shp
#   to move the point to the correct stream cell, then:
#     source("reset_workflow.R")
#     reset_site("site_id")   # clears only that site's output folder
#   Re-run Stage 5 — the edited snapped pour point will be detected and
#   the workflow will resume from the watershed step.
# =============================================================================

# -- Packages ----------------------------------------------------------------
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

# -- Source all modules ------------------------------------------------------
source("R/utils.R")
source("R/01_group_sites.R")
source("R/02_prepare_dem.R")
source("R/03_burn_streams.R")
source("R/04_run_whitebox.R")
source("R/05_delineate_sites.R")
source("reset_workflow.R")

# =============================================================================
# CONFIGURATION — fill in before running
# =============================================================================

# Sites — defined as a tibble (recommended) or path to CSV
# Option A: tibble (edit sites_template.R and source it)
source("sites_template.R") # loads `sites` tibble into environment
# Option B: CSV file
# sites_csv <- "sites_template.csv"

# Root directories
output_dir <- "output"
cache_dir <- "cache"

# Path to MRDEM .vrt file
mrdem_vrt <- "/path/to/mrdem-30-dtm.vrt"

# Path to NHN root directory (containing nhn_rhn_*_gdb_en subfolders)
nhn_dir <- "/Users/sam/Documents/cfs/shared_data/raw/hydro/networks/NHN/gdb"

# Path to NHN index shapefile
nhn_index <- "/Users/sam/Documents/cfs/shared_data/raw/hydro/networks/NHN/NHN_INDEX_WORKUNIT_LIMIT_2/NHN_INDEX_WORKUNIT_LIMIT_2.shp"

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

# Validate sites — use validate_sites_tibble() if using tibble approach,
# or validate_sites() if using CSV
sites <- validate_sites_tibble(sites)
# sites <- validate_sites(sites_csv)  # CSV alternative

group_manifest <- build_group_manifest(
  sites = sites,
  output_dir = output_dir,
  cache_dir = cache_dir,
  hydrobasins_dir = hydrobasins_dir,
  default_buffer_m = default_buffer_m
)

# Write group AOIs for visual inspection in QGIS before proceeding
write_group_aois(group_manifest, cache_dir)

message(glue(
  "\nGroup manifest built. Inspect cache/group_aois.gpkg in QGIS to verify ",
  "AOI extents before proceeding to Stage 2.\n"
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

# Quick verification — check all group outputs exist and flow accum looks sane
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

# Combine all catchment polygons into a single file for QA in QGIS
all_catchments <- purrr::map(sites$site_id, function(sid) {
  path <- fs::path(site_output_dir(output_dir, sid), "catchment.gpkg")
  if (!cache_exists(path)) {
    return(NULL)
  }
  sf::st_read(path, quiet = TRUE)
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
  "\nWorkflow complete. {nrow(results)} site(s) processed. ",
  "Combined catchments saved to {output_dir}/all_catchments.gpkg"
))

# Report any flagged sites
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

# =============================================================================
# STAGE 6 (optional) — Catchment morphometric metrics
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
