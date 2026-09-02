# add_canlcc_loi.R
# =============================================================================
# One-off backfill: adds the newly-introduced "canlcc" LOI (see
# workflow/CELESTE/README.md's "canlcc" note) to an already-complete
# CELESTE_hydroweight.csv, without re-running any catchment-producing
# stage (1-6) — Stages 1-6 are untouched, so no catchment.gpkg,
# catchment_clipped.gpkg, all_catchments*.gpkg, or catchment_metrics.csv
# is read for writing or modified.
#
# Approach: rebuilds only the in-memory objects calculate_hydroweight_
# attributes_stream() needs (sites/group_manifest/config, via Stage 1 —
# pure AOI/config resolution, no raster writes), computes hydroweight
# attributes for JUST the "canlcc" LOI (identical distance-weighting
# result whether canlcc is computed alone or alongside every other LOI —
# hydroweight() itself doesn't depend on which LOI is attached afterward),
# then LEFT JOINs the new canlcc columns onto the existing
# CELESTE_hydroweight.csv by (site, version) — the join naturally
# restricts to whatever rows already survived drop_redundant_clipped_
# rows() in the original run, so that filtering doesn't need to be
# reproduced here. Finally reruns tidy_celeste_outputs() (a pure CSV
# reshape) to regenerate output_dir/tidy/, including the new
# canlcc_long.csv.
#
# Usage: source from an R session in the project root, after
# output/CELESTE/CELESTE_hydroweight.csv already exists from a prior
# run_celeste.R run.
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
source(here("workflow/R/stream/hydroweight_attributes.R"))
source(here("workflow/R/engine/00_resolve_config.R"))
source(here("workflow/R/engine/01_build_group_manifest.R"))
source(here("workflow/CELESTE/tidy_outputs.R"))

# =============================================================================
# STAGE 1 (config + group manifest only — no terrain/delineation calls)
# =============================================================================
# Identical config to run_celeste.R's own Stage 1, so sites/group_manifest
# resolve to exactly the same values already on disk.

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
# canlcc-only hydroweight (matches run_celeste.R's canlcc entry exactly)
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

canlcc_layers <- list(
  list(path_lazy = CANLCC_PATH, name = "canlcc", type = "categorical", class_levels = canlcc_levels)
)

canlcc_results <- calculate_hydroweight_attributes_stream(
  sites = sites,
  group_manifest = group_manifest,
  output_dir = output_dir,
  cache_dir = cache_dir,
  loi_layers = canlcc_layers,
  catchment_versions = c("unclipped", "clipped")
)

cw_inform(glue::glue("canlcc computed for {nrow(canlcc_results)} site x version row(s)."))

# =============================================================================
# Merge into the existing CELESTE_hydroweight.csv (join by site, version —
# NOT merge_rows_into_csv(), since we're adding COLUMNS to existing rows,
# not replacing/adding ROWS)
# =============================================================================

hydroweight_path <- fs::path(output_dir, "CELESTE_hydroweight.csv")
if (!fs::file_exists(hydroweight_path)) {
  cw_abort(glue::glue("Not found: {hydroweight_path} — run run_celeste.R first."))
}

hw_existing <- readr::read_csv(hydroweight_path, show_col_types = FALSE)

new_canlcc_cols <- setdiff(names(canlcc_results), c("site", "version"))
already_present <- intersect(new_canlcc_cols, names(hw_existing))
if (length(already_present) > 0) {
  cw_abort(glue::glue(
    "canlcc column(s) already present in {hydroweight_path} — refusing to ",
    "overwrite: {paste(already_present, collapse = ', ')}"
  ))
}

hw_merged <- dplyr::left_join(hw_existing, canlcc_results, by = c("site", "version"))

n_before <- nrow(hw_existing)
n_after  <- nrow(hw_merged)
if (n_after != n_before) {
  cw_abort(glue::glue(
    "Row count changed after join ({n_before} -> {n_after}) — the join key ",
    "(site, version) didn't line up 1:1. Aborting without writing."
  ))
}

n_missing_canlcc <- sum(is.na(hw_merged[[paste0("canlcc_", canlcc_levels$Class[1], "_lumped_prop")]]))
if (n_missing_canlcc > 0) {
  cw_warn(glue::glue(
    "{n_missing_canlcc} existing row(s) got no canlcc match (site/version ",
    "present in the CSV but not in canlcc_results) — check for a group/site ",
    "mismatch before trusting canlcc_long.csv fully."
  ))
}

readr::write_csv(hw_merged, hydroweight_path)
cw_inform(glue::glue(
  "canlcc merged into {hydroweight_path} ({n_before} row(s) unchanged, ",
  "{length(new_canlcc_cols)} new column(s) added)."
))

# =============================================================================
# Regenerate tidy tables (pure CSV reshape — no catchment/raster I/O)
# =============================================================================

tidy_celeste_outputs(output_dir = output_dir)
