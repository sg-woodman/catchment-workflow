# run_celeste_engine.R
# =============================================================================
# Full CELESTE run (all 6 HydroBasins groups, 132 sites) on the modular
# engine (workflow/R/engine/), using the CORRECT, current site coordinates
# (data/celeste_milli_sites_clean.gpkg — confirmed by the user 2026-08-26 to
# match what they were actually sent, correct site count and coordinates).
#
# WHY THIS EXISTS, SEPARATE FROM run_celeste.R:
#   The existing production run (output/CELESTE/, from workflow/R/stream/)
#   was built from DIFFERENT, apparently-outdated site coordinates for 112
#   of 132 sites (up to ~100 m off, confirmed by comparing every site's
#   current source coordinate against what's actually stored in each site's
#   own output/CELESTE/<site>/pour_point.gpkg — see
#   workflow/CELESTE/test_engine_parity.R's TUR-only precursor to this
#   script for the single-group version of that finding). This script is a
#   from-scratch delineation of every CELESTE site using the CURRENT,
#   correct coordinates, on the engine (not workflow/R/stream/, per the
#   user's explicit ask — this doubles as the engine's first real multi-
#   group run, not just a single-group parity smoke test).
#
#   output/CELESTE and cache/CELESTE (the production run) are DELIBERATELY
#   UNTOUCHED by this script — separate project_id, separate output_dir/
#   cache_dir — so the two can be compared side by side. Do not point this
#   script at the production directories.
#
# BURN-IN: uses streams_burn = "nhn_auto", NOT "supplied" (unlike
# test_engine_parity.R's TUR-only precursor, which reused a single cached
# flowlines.gpkg for one group). Reason: with corrected site coordinates,
# each group's AOI needs to be freely re-derived (a shifted site could,
# in principle, change which HydroBasins polygon/NHN work units apply),
# and "supplied" only accepts one fixed vector for the whole run,
# unsuited to CELESTE's 6 independent regional groups. download_nhn()
# (engine/03_prepare_streams_burn.R) checks nhn_raw_dir for an
# already-downloaded GDB matching the required WSCSSDA code BEFORE ever
# hitting the FTP server — pointed at the same local NHN directory
# run_celeste.R already uses, this resolves entirely from local cache, no
# network needed (verified: every group's local GDB is already present
# from the original run).
#
# WHAT TO WATCH FOR (per the user's explicit ask): whether burn-in
# actually applies for every group, or whether any group's (corrected)
# AOI now falls outside NHN coverage and silently proceeds unburned.
# Checked before this run: in the CURRENT production cache (cache/CELESTE/),
# all 6 groups have non-empty flowlines.gpkg (23K-113K features) and a
# real dem_burned.tif, and every site's burn_streams flag is TRUE — no
# group was already falling back to "no burn" as of the last production
# run. This script logs and summarizes burn-in status per group again
# under the corrected coordinates, since an AOI shift COULD in principle
# change NHN coverage even if it didn't for TUR specifically (verified
# unchanged in test_engine_parity.R).
#
# Usage: source from an R session in the project root. This is a genuinely
# long run — 6 groups, the largest (MOR) DEM crop is ~5x TUR's (the group
# test_engine_parity.R used, which took ~3 min for terrain prep alone) —
# expect this to take a good while longer than that single-group test.
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
source(here("workflow/R/stream/01_group_sites.R")) # build_group_manifest() — reused unmodified by the engine's "hydrobasins" strategy
source(here("workflow/R/stream/03_burn_streams.R")) # burn_streams_into_dem() — reused unmodified by resolve_streams_burn()
source(here("workflow/R/stream/05_delineate_sites.R")) # snap_pour_point(), delineate_watershed(), etc. — reused unmodified by engine/04
source(here("workflow/R/engine/00_resolve_config.R"))
source(here("workflow/R/engine/01_build_group_manifest.R"))
source(here("workflow/R/engine/02_prepare_terrain.R"))
source(here("workflow/R/engine/03_prepare_streams_burn.R"))
source(here("workflow/R/engine/04_delineate_site.R"))

# =============================================================================
# CONFIGURATION — mirrors run_celeste.R's own values (see file header for
# what's deliberately different: nhn_auto instead of a per-group supplied path)
# =============================================================================

PROJECT_ID <- "CELESTE_engine"
output_dir <- here("output", PROJECT_ID)
cache_dir  <- here("cache", PROJECT_ID)
fs::dir_create(output_dir, recurse = TRUE)
fs::dir_create(cache_dir, recurse = TRUE)

mrdem_vrt       <- "~/Documents/cfs/shared_data/raw/dem/mrdem-30-dtm.vrt"
hydrobasins_dir <- "/Users/sam/Documents/cfs/shared_data/raw/hydro/watersheds/HydroBasins"
nhn_dir         <- "/Users/sam/Documents/cfs/shared_data/raw/hydro/networks/NHN/gdb"
nhn_index       <- "/Users/sam/Documents/cfs/shared_data/raw/hydro/networks/NHN/NHN_INDEX_WORKUNIT_LIMIT_2/NHN_INDEX_22_INDEX_WORKUNIT_LIMIT_2.shp"

sites_sf <- st_read(here("data/celeste_milli_sites_clean.gpkg"), quiet = TRUE)
coords <- st_coordinates(sites_sf)
sites <- sites_sf |>
  st_drop_geometry() |>
  as_tibble() |>
  mutate(lon = coords[, "X"], lat = coords[, "Y"]) |>
  select(site_id, site_name, lon, lat, group_id, burn_streams, aoi_buffer_m)

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
  crs = NULL, # NULL = MRDEM's own native CRS (EPSG:3979) — matches run_celeste.R's target_crs = 3979
  stream_threshold = 1000, # matches run_celeste.R's default_stream_threshold

  streams_burn   = list(source = "nhn_auto"), # see file header for why, not "supplied"
  nhn_index_path = nhn_index,
  nhn_raw_dir    = nhn_dir,

  lake_polygons = NULL,
  lake_buffer_m = 30, # unused in point pour-point mode

  grouping = list(
    strategy = "hydrobasins",
    hydrobasins_dir = hydrobasins_dir,
    default_buffer_m = 1000 # matches run_celeste.R's default_buffer_m
  ),

  loi_layers = NULL
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
# Runs once per group — the expensive, one-time part. Logs (and this script
# summarizes below) whether burn-in actually applied per group.

prepare_engine_terrain(config, group_manifest)

burn_status <- purrr::map_dfr(group_manifest$group_id, function(grp) {
  fl_path <- fs::path(cache_dir, grp, "flowlines.gpkg")
  burned_path <- fs::path(cache_dir, grp, "dem_burned.tif")
  breached_path <- fs::path(cache_dir, grp, "dem_breached.tif")
  n_fl <- if (fs::file_exists(fl_path)) nrow(sf::st_read(fl_path, quiet = TRUE)) else NA_integer_
  tibble(
    group_id = grp,
    n_flowline_features = n_fl,
    dem_burned_written = fs::file_exists(burned_path),
    dem_breached_written = fs::file_exists(breached_path)
  )
})
cw_inform("\n--- Burn-in status per group (STAGE 2 summary) ---")
print(burn_status)
readr::write_csv(burn_status, fs::path(output_dir, "burn_in_status.csv"))
no_burn <- dplyr::filter(burn_status, is.na(n_flowline_features) | n_flowline_features == 0)
if (nrow(no_burn) > 0) {
  cw_warn(glue::glue(
    "\n{nrow(no_burn)} group(s) have NO burn-in streams (empty or missing ",
    "flowlines): {paste(no_burn$group_id, collapse = ', ')} — these groups' ",
    "dem_breached.tif was derived from the UNBURNED dem.tif. Compare against ",
    "cache/CELESTE/<group>/flowlines.gpkg (the production run) for the same groups."
  ))
} else {
  cw_inform("All groups have non-empty burn-in streams — consistent with the production run.")
}

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
  "{nrow(all_catchments)} catchment(s) written to ",
  "{output_dir}/all_catchments.gpkg. output/CELESTE and cache/CELESTE ",
  "(the production run) were not touched by this script."
))
