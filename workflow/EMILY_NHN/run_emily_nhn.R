# run_emily_nhn.R
# =============================================================================
# Top-level runner for EMILY_NHN — the SAME 177 sites as EMILY_TURKEY
# (workflow/EMILY_TURKEY/run_emily_turkey.R), delineated with CELESTE's
# approach instead of CAM's: MRDEM (raw DEM, needs breach) + NHN auto-
# download stream burn-in + upfront NHN lake conditioning + HydroBasins
# grouping, EPSG:3979 — see workflow/CELESTE/run_celeste.R, the direct
# model for every config choice below. Point pour points (Jenson snap),
# same as EMILY_TURKEY.
#
# WHY A SEPARATE PROJECT rather than re-running EMILY_TURKEY differently:
# same sites, deliberately different terrain/conditioning method, so the
# two projects' outputs (output/EMILY_TURKEY/ vs output/EMILY_NHN/) can be
# compared directly — that's the whole point of this run.
#
# GROUPING: unlike CELESTE (group_id pre-assigned as semantic project
# labels — "COC"/"NIP"/etc — when its sites gpkg was first built),
# EMILY_NHN's sites have no such pre-existing structure. group_id is
# derived here by workflow/EMILY_NHN/assign_hydrobasins_groups.R —
# spatially assigning each site to the HydroBasins level-6 polygon it
# falls in (group_id = "HYBAS_<id>"). Verified against the real data
# before use: 177/177 sites matched exactly one polygon each, producing 5
# groups (sizes 2/4/7/23/141) — see that file's header for why adjacent
# polygons are NOT merged into fewer groups (they're all mutually
# touching; merging would collapse everything into one 177-site group).
#
# BURN-IN: streams_burn$source = "nhn_auto" for every site, no per-group
# override (unlike CELESTE's COC exception — not applicable here; the
# user explicitly asked for burn-in on for all EMILY_NHN sites).
#
# WORKFLOW STAGES: identical shape to run_celeste.R, minus Stage 8 (no
# tidy_outputs.R equivalent written for this project) and Stage 7 pared
# down to CANLCC only (no per-project ndvi/ndvi_trend/harvest_regen prep
# exists for this AOI — same scope EMILY_TURKEY's Stage 7 used):
#   Stage 1 — Resolve config, build group manifest (HydroBasins grouping,
#              group_id from assign_hydrobasins_group_ids() above)
#   Stage 2 — Prepare terrain (crop MRDEM, burn NHN streams in, upfront
#              NHN lake conditioning, breach, D8 pointer, flow
#              accumulation, extract streams) — once per group
#   Stage 3 — Delineate catchments (point pour point, Jenson snap)
#   Stage 4 — Remove upstream nested catchments
#   Stage 5 — Re-clip rasters/flowlines to clipped catchments (REQUIRED
#              before Stage 7)
#   Stage 6 — Catchment morphometric metrics
#   Stage 7 — CANLCC hydroweight (same path/class table as EMILY_TURKEY's
#              Stage 7 and CELESTE's canlcc LOI)
#
# CACHING: cache/EMILY_NHN/HYBAS_<id>/ (one subfolder per HydroBasins
# group). Re-running skips any step whose output already exists.
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
source(here("workflow/R/stream/group_sites.R")) # build_group_manifest(), resolve_hydrobasins_region() — reused unmodified by the engine's "hydrobasins" strategy
source(here("workflow/R/stream/burn_streams.R")) # burn_streams_into_dem() — reused unmodified by resolve_streams_burn()
source(here("workflow/R/stream/delineate_sites.R")) # snap_pour_point(), delineate_watershed(), etc. — reused unmodified by engine/04
source(here("workflow/R/stream/hydroweight_attributes.R")) # calculate_hydroweight_attributes_stream() — shared with CAM/CELESTE, reused unmodified
source(here("workflow/R/remove_upstream.R"))
source(here("workflow/R/reclip_outputs.R")) # reclip_outputs() — required before Stage 7
source(here("workflow/R/catchment_metrics.R"))
source(here("workflow/R/engine/00_resolve_config.R"))
source(here("workflow/R/engine/01_build_group_manifest.R"))
source(here("workflow/R/engine/02_prepare_terrain.R"))
source(here("workflow/R/engine/03_prepare_streams_burn.R"))
source(here("workflow/R/engine/prepare_lake_conditioning.R")) # resolve_lake_conditioning() -- opt-in via config$lake_conditioning
source(here("workflow/R/engine/04_delineate_site.R"))
source(here("workflow/R/engine/99_rerun_sites")) # drop_redundant_clipped_rows(); rerun_engine_site_watershed()/rerun_engine_sites() for post-hoc corrections
source(here("workflow/EMILY_NHN/assign_hydrobasins_groups.R")) # assign_hydrobasins_group_ids() — see header above

# =============================================================================
# CONFIGURATION
# =============================================================================

PROJECT_ID <- "EMILY_NHN"
output_dir <- here("output", PROJECT_ID)
cache_dir  <- here("cache", PROJECT_ID)
fs::dir_create(output_dir, recurse = TRUE)
fs::dir_create(cache_dir, recurse = TRUE)

mrdem_vrt       <- "~/Documents/cfs/shared_data/raw/dem/mrdem-30-dtm.vrt"
hydrobasins_dir <- "/Users/sam/Documents/cfs/shared_data/raw/hydro/watersheds/HydroBasins"
nhn_dir         <- "/Users/sam/Documents/cfs/shared_data/raw/hydro/networks/NHN/gdb"
nhn_index       <- "/Users/sam/Documents/cfs/shared_data/raw/hydro/networks/NHN/NHN_INDEX_WORKUNIT_LIMIT_2/NHN_INDEX_22_INDEX_WORKUNIT_LIMIT_2.shp"

# =============================================================================
# SITE DEFINITIONS
# =============================================================================
# Same source data and site_id derivation as EMILY_TURKEY (see that
# project's README for the HARM030 disambiguation rationale) — reusing
# data/emily_turkey_sites_raw.csv directly rather than duplicating it,
# since it's the same underlying site list either way.

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
  dplyr::select(site_id, site_name, lon, lat)

stopifnot(!anyDuplicated(sites$site_id))

# group_id: derived from HydroBasins level-6 membership, not carried over
# from EMILY_TURKEY (which never assigned one — "whole_domain"/
# "manual_groups" don't need it). See assign_hydrobasins_groups.R's header.
sites <- assign_hydrobasins_group_ids(sites, hydrobasins_dir = hydrobasins_dir, hybas_level = 6)

# burn_streams: TRUE for every site, per explicit instruction — no CELESTE-
# style per-group override here.
sites$burn_streams <- TRUE

# aoi_buffer_m: NA for every site -> grouping$default_buffer_m (1000 m,
# same as CELESTE) applies uniformly. Column must exist (get_groups()/
# validate_sites_impl() in workflow/R/utils.R reference it directly) even
# though no group overrides it here.
sites$aoi_buffer_m <- NA_real_

cw_inform(glue::glue(
  "EMILY_NHN engine run: {nrow(sites)} site(s) across ",
  "{length(unique(sites$group_id))} group(s): ",
  "{paste(sort(unique(sites$group_id)), collapse = ', ')}."
))

# =============================================================================
# STAGE 1 — Resolve config, build group manifest
# =============================================================================

run_config <- list(
  project_id = PROJECT_ID,
  output_dir = output_dir,
  cache_dir  = cache_dir,
  sites      = sites,

  dem = list(path = mrdem_vrt),
  flow_direction = NULL,
  flow_pointer   = NULL,
  crs = NULL, # NULL = MRDEM's own native CRS (EPSG:3979), same as CELESTE
  stream_threshold = 1000, # same default as CELESTE (same terrain source)

  streams_burn   = list(source = "nhn_auto"), # ALWAYS on for EMILY_NHN, per instruction — no per-group override
  nhn_index_path = nhn_index,
  nhn_raw_dir    = nhn_dir,

  # Upfront lake conditioning (workflow/R/engine/prepare_lake_conditioning.R)
  # — same as CELESTE's current config, prevents lake bisection instead of
  # patching it after the fact. See CELESTE's README for the verification
  # this was based on.
  lake_conditioning = list(source = "nhn_auto"),

  lake_polygons = NULL, # point pour point mode
  lake_buffer_m = 30,   # unused in point pour-point mode

  grouping = list(
    strategy = "hydrobasins",
    hydrobasins_dir = hydrobasins_dir,
    default_buffer_m = 1000
  ),

  loi_layers = NULL # set in Stage 7 below
)

config <- resolve_engine_config(run_config)

gm <- build_engine_group_manifest(config)
sites <- gm$sites
group_manifest <- gm$group_manifest
print(group_manifest)

# =============================================================================
# STAGE 2 — Prepare terrain (crop, burn, breach, D8, accumulation, streams)
# =============================================================================

prepare_engine_terrain(config, group_manifest)

burn_status <- purrr::map_dfr(group_manifest$group_id, function(grp) {
  burned_path   <- fs::path(cache_dir, grp, "dem_burned.tif")
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
# REQUIRED before Stage 7.

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
# STAGE 7 — CANLCC hydroweight
# =============================================================================
# Same source raster + class lookup table as EMILY_TURKEY/CAM/CELESTE.
# CANLCC is already EPSG:3979, matching this project's own working CRS
# exactly (MRDEM's native CRS) — no reprojection needed, same as CELESTE.

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

hw_results <- calculate_hydroweight_attributes_stream(
  sites = sites,
  group_manifest = group_manifest,
  output_dir = output_dir,
  cache_dir = cache_dir,
  loi_layers = loi_layers,
  catchment_versions = c("unclipped", "clipped"),
  raster_crs = config$working_crs # EPSG:3979 — matches the function's own default here, passed explicitly for consistency with EMILY_TURKEY/CAM's runners
)

hw_results <- drop_redundant_clipped_rows(hw_results, upstream_results, site_col = "site")

write_csv(hw_results, fs::path(output_dir, paste0(PROJECT_ID, "_hydroweight.csv")))
print(hw_results)

cw_inform(glue::glue(
  "\nHydroweighting complete: {nrow(hw_results)} row(s) written to ",
  "{output_dir}/{PROJECT_ID}_hydroweight.csv."
))
