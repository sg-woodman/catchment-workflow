# run_cam_streams.R
# =============================================================================
# Top-level runner for the CAM stream-site catchment delineation workflow,
# on the modular, input-driven engine (workflow/R/engine/) rather than
# either of the two existing project-specific pipelines.
#
# WHY THIS ISN'T ON workflow/R/stream/ (CELESTE's pipeline) OR
# workflow/R/lake/ (CAM lakes' pipeline):
#   These 39 sites are stream sampling locations near (but not matching)
#   the 45 CAM lake sites, plus a separate set of MOE long-term monitoring
#   gauges — delineated as POINT pour points (Jenson snap), like CELESTE,
#   but against the SAME OIH DEM/flow-direction CAM lakes uses (EPSG:3161),
#   not MRDEM (EPSG:3979) — a combination neither existing pipeline
#   supports directly. See workflow/R/engine/ for the general mechanism;
#   see /Users/sam/.claude/plans/i-need-to-run-golden-nova.md for the full
#   design writeup.
#
# WORKFLOW STAGES:
#   Stage 1 — Resolve config, build group manifest (single "whole_domain"
#              group here — OIH is one project-wide, already-conditioned
#              raster set, same as CAM lakes; no HydroBasins grouping)
#   Stage 2 — Prepare terrain (recode OIH flow direction to WhiteboxTools
#              encoding, derive flow accumulation + streams; no breaching —
#              OIH is already hydrologically conditioned)
#   Stage 3 — Delineate catchments (point pour point, Jenson snap)
#   Stage 4 — Remove upstream nested catchments
#   Stage 5 — Re-clip rasters/flowlines to clipped catchments
#   Stage 6 — Catchment morphometric metrics (optional)
#   Stage 7 — Distance-weighted catchment attributes / hydroweight
#              (optional) — reuses workflow/R/stream/06_hydroweight_
#              attributes.R's calculate_hydroweight_attributes_stream()
#              completely unmodified, just passing raster_crs =
#              config$working_crs (EPSG:3161) instead of its default
#              EPSG:3979. This works because the engine's group_manifest
#              has the exact same shape CELESTE's does, and
#              02_prepare_terrain.R always materializes the canonical
#              dem_breached.tif/flow_accum.tif filenames that function
#              expects, regardless of terrain source.
#
# CACHING:
#   cache/CAM/stream_delineation/ (one flat "group" — no per-region
#   subfolders, since grouping$strategy = "whole_domain"). Re-running
#   skips any step whose output already exists.
#
# RERUNNING AFTER A CORRECTION (workflow/R/engine/99_rerun_sites):
#   Two scenarios, both handled by rerun_engine_sites() in one call —
#   redoes delineation for just the named site(s), then remove-upstream
#   (always the full sites list — cheap, and a changed site's shape can
#   change which OTHER sites it's nested within), reclip, metrics, and
#   hydroweight for just the affected sites (including any cascaded
#   neighbor whose clipped catchment also changed), MERGED into the
#   existing catchment_metrics.csv / *_hydroweight.csv rather than
#   overwritten wholesale:
#
#   1. Raw coordinate was wrong, fixed in data/cam_stream_sites_raw.csv:
#        sites_raw <- readr::read_csv(here("data/cam_stream_sites_raw.csv"))
#        # ...re-derive `sites` the same way as above...
#        rerun_engine_sites(
#          resnap_site_ids = "Bell", sites = sites, group_manifest = group_manifest,
#          config = config, output_dir = output_dir, cache_dir = cache_dir,
#          loi_layers = loi_layers # or NULL to skip re-running hydroweight
#        )
#
#   2. Snapped pour point was wrong, manually edited in QGIS (edit
#      output/<site_id>/pour_point_snapped.shp directly, then):
#        rerun_engine_sites(
#          edited_snap_site_ids = "Wolf", sites = sites, group_manifest = group_manifest,
#          config = config, output_dir = output_dir, cache_dir = cache_dir,
#          loi_layers = loi_layers
#        )
#
#   Both in one call (each site handled the correct way):
#        rerun_engine_sites(
#          resnap_site_ids = "Bell", edited_snap_site_ids = "Wolf",
#          sites = sites, group_manifest = group_manifest, config = config,
#          output_dir = output_dir, cache_dir = cache_dir, loi_layers = loi_layers
#        )
#
# DATA QUALITY FLAG:
#   Source xlsx row 30, SUD11, has lon = -51.2 — 2500+ km off from its
#   neighbors (SUD12/VER01, both ~-81.0), almost certainly a typo. Excluded
#   below via EXCLUDED_SITE_IDS; verify against source before re-including.
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

# -- Source modules --------------------------------------------------------
source(here("workflow/R/utils.R"))
source(here("workflow/R/stream/03_burn_streams.R")) # burn_streams_into_dem() — reused by engine 03 (not exercised: streams_burn source = "none" here)
source(here("workflow/R/stream/05_delineate_sites.R")) # load_group_rasters(), snap_pour_point(), delineate_watershed(), clip_rasters_to_catchment(), clip_flowlines_to_catchment() — reused unmodified by engine 04
source(here("workflow/R/stream/06_hydroweight_attributes.R")) # calculate_hydroweight_attributes_stream() — reused unmodified, raster_crs overridden below
source(here("workflow/R/06_remove_upstream.R"))
source(here("workflow/R/07_reclip_outputs.R"))
source(here("workflow/R/08_catchment_metrics.R"))
source(here("workflow/R/engine/00_resolve_config.R"))
source(here("workflow/R/engine/01_build_group_manifest.R"))
source(here("workflow/R/engine/02_prepare_terrain.R"))
source(here("workflow/R/engine/03_prepare_streams_burn.R"))
source(here("workflow/R/engine/04_delineate_site.R"))
source(here("code/reset_workflow.R")) # reset_site() — used by rerun_engine_sites()'s resnap path
source(here("workflow/R/engine/99_rerun_sites")) # rerun_engine_site_watershed(), rerun_engine_sites() — see "RERUNNING AFTER A CORRECTION" above
source(here("workflow/gee_utils.R")) # prepare_polygons_for_gee()/group_polygons_by_adjacency() — used by prepare_ndvi.R to re-derive the site->NDVI-file mapping
source(here("workflow/raster_attributes.R")) # sens_slope_trend() — used by prepare_ndvi_trend.R
source(here("workflow/CAM/prepare_ndvi.R")) # build_cam_ndvi_site_map(), clean_cam_ndvi_tile(), prepare_cam_ndvi_site_rasters() — CAM stream-site NDVI LOI prep (see Stage 7)
source(here("workflow/CAM/prepare_ndvi_trend.R")) # prepare_cam_ndvi_trend_site_rasters() — Sen's-slope NDVI trend LOI prep
source(here("workflow/CAM/prepare_harvest_regen.R")) # prepare_cam_harvest_regen_rasters() — Ontario MNR harvest/regen LOI prep
# If a future CAM-streams-like project needs HydroBasins-based grouping
# (grouping$strategy = "hydrobasins" — a large national-mosaic DEM worth
# cropping per region, unlike OIH's single project-wide raster set), also
# source: workflow/R/stream/01_group_sites.R

# =============================================================================
# CONFIGURATION — verify paths before running
# =============================================================================

PROJECT_ID  <- "CAM"
DELINEATION <- "stream_delineation"

output_dir <- here("output", PROJECT_ID, DELINEATION)
cache_dir  <- here("cache", PROJECT_ID, DELINEATION)
fs::dir_create(output_dir, recurse = TRUE)
fs::dir_create(cache_dir, recurse = TRUE)

# OIH Enforced DEM and Enhanced Flow Direction — same physical files
# run_cam_lakes.R uses (kept untouched by this file).
OIH_DEM_PATH       <- "/Users/sam/Downloads/IntegratedHydrologyNE/EnforcedDEM.tif"
OIH_FLOW_DIR_PATH  <- "/Users/sam/Downloads/IntegratedHydrologyNE/EnhancedFlowDirection.tif"

# OIH -> WhiteboxTools D8 encoding reclassification (one-step rotation),
# same table as workflow/R/lake/02_prepare_oih_dem.R — generic in the
# engine, applied here as CAM's own config value:
#   OIH 128 (NE) -> WBT 1     OIH 1  (E)  -> WBT 2
#   OIH 2   (SE) -> WBT 4     OIH 4  (S)  -> WBT 8
#   OIH 8   (SW) -> WBT 16    OIH 16 (W)  -> WBT 32
#   OIH 32  (NW) -> WBT 64    OIH 64 (N)  -> WBT 128
oih_recode_matrix <- matrix(
  c(128, 1, 1, 2, 2, 4, 4, 8, 8, 16, 16, 32, 32, 64, 64, 128),
  ncol = 2, byrow = TRUE
)

STREAM_THRESHOLD <- 100 # flow accum threshold for stream extraction (cells) — same default as run_cam_lakes.R
SNAP_DIST_M      <- 200 # point-mode pour point snap distance (m)
MIN_CELLS        <- 10  # minimum catchment size before flagging

# =============================================================================
# SITE DEFINITIONS
# =============================================================================
# Converted from /Users/sam/Downloads/Stream Coordinates Summer and Fall -
# Cameron Lefebvre.xlsx, Sheet1, via data/cam_stream_sites_raw.csv (39
# sites: 25 "Summer"/CRADLES-lakes-adjacent + 14 "Fall"/MOE long-term
# monitoring gauges; blank separator row already dropped in the CSV).

EXCLUDED_SITE_IDS <- c("SUD11") # lon = -51.2 in source — 2500+ km off from
# neighbouring SUD12/VER01 (both ~-81.0) — almost certainly a typo in the
# source workbook. Verify against Cameron Lefebvre's original data before
# re-including.

sites_raw <- readr::read_csv(here("data/cam_stream_sites_raw.csv"), show_col_types = FALSE)

sites <- sites_raw |>
  dplyr::mutate(
    site_id = stream_id |>
      gsub("\\s+", "_", x = _) |>
      gsub("[^A-Za-z0-9_-]", "", x = _),
    site_name = stream_id
  ) |>
  dplyr::filter(!site_id %in% EXCLUDED_SITE_IDS) |>
  dplyr::select(site_id, site_name, lon, lat)

# =============================================================================
# STAGE 1 — Resolve config, build group manifest
# =============================================================================

run_config <- list(
  project_id = "CAM_streams", # internal label / group_id (whole_domain = one group)
  output_dir = output_dir,
  cache_dir  = cache_dir,
  sites      = sites,

  dem = list(path = OIH_DEM_PATH), # elevation surface for per-site clipping/output + hydroweight — decoupled from conditioning (flow_direction supersedes it below)
  flow_direction = list(path = OIH_FLOW_DIR_PATH, recode = oih_recode_matrix), # already hydrologically conditioned — no breaching
  flow_pointer = NULL,
  flow_accum   = NULL,
  crs = NULL, # NULL = match OIH's own native CRS (EPSG:3161) — already projected, metre-based, suitable for delineation, so no override needed here
  stream_threshold = STREAM_THRESHOLD,

  streams_burn = list(source = "none"), # OIH flow direction is already conditioned — burning would have to happen before conditioning, not after

  lake_polygons = NULL, # point pour point mode
  lake_buffer_m = 30,

  grouping = list(strategy = "whole_domain"), # OIH covers the whole domain as one project-wide raster set, same as CAM lakes today

  loi_layers = NULL # set in Stage 7 below
)

config <- resolve_engine_config(run_config)

gm <- build_engine_group_manifest(config)
sites <- gm$sites
group_manifest <- gm$group_manifest

print(group_manifest)

message(glue(
  "\nStage 1 complete. {nrow(sites)} site(s) in {nrow(group_manifest)} group(s). ",
  "Proceed to Stage 2.\n"
))

# =============================================================================
# STAGE 2 — Prepare terrain
# =============================================================================
# Caches dem.tif / dem_breached.tif (alias of dem.tif — no real breaching,
# OIH is pre-conditioned) / flow_pointer.tif / flow_accum.tif / streams.tif /
# hillshade.tif in cache/CAM/stream_delineation/.

prepare_engine_terrain(config, group_manifest)

# Spot-check: load streams.tif and dem.tif in QGIS before proceeding. If
# streams look wrong, adjust STREAM_THRESHOLD and re-run Stage 2 (delete
# cache/CAM/stream_delineation/streams.tif first).

# =============================================================================
# STAGE 3 — Delineate catchments (point pour point)
# =============================================================================

results <- delineate_engine_catchments(
  config = config, sites = sites, group_manifest = group_manifest,
  snap_dist = SNAP_DIST_M, min_cells = MIN_CELLS
)
print(results)

flagged <- dplyr::filter(results, flagged)
if (nrow(flagged) > 0) {
  message(glue("\n{nrow(flagged)} site(s) flagged — review pour points:"))
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

message(glue(
  "\nStage 3 complete. {nrow(results)} site(s) processed. ",
  "Combined catchments saved to {output_dir}/all_catchments.gpkg"
))

# =============================================================================
# STAGE 4 — Remove upstream nested catchments
# =============================================================================
# upstream_results is kept — Stages 6-7 use its n_erased column to drop
# "clipped" rows that are pure duplicates of "unclipped" (no nested
# upstream catchment was actually erased) from the output CSVs.

upstream_results <- remove_upstream_catchments(sites, output_dir)

# =============================================================================
# STAGE 5 — Re-clip rasters and flowlines to clipped catchments
# =============================================================================

reclip_outputs(sites = sites, output_dir = output_dir, group_manifest = group_manifest)

# =============================================================================
# STAGE 6 (optional) — Catchment morphometric metrics
# =============================================================================

metrics <- calculate_catchment_metrics(sites = sites, output_dir = output_dir)
# Drop "clipped" rows that duplicate "unclipped" (no nested upstream
# catchment actually erased for that site) — see Stage 4's comment.
metrics <- drop_redundant_clipped_rows(metrics, upstream_results, site_col = "site_id")
ref_table <- build_metrics_reference_table()
write_metrics_outputs(metrics = metrics, ref_table = ref_table, output_dir = output_dir)
print(metrics)

# =============================================================================
# STAGE 7 (optional) — Distance-weighted catchment attributes (hydroweight)
# =============================================================================
# Reuses CAM's existing project-wide LOI rasters (from run_cam_lakes.R), but
# as `path_lazy` rather than `path` entries — deliberately, not just for
# consistency. With `path`, calculate_hydroweight_attributes_stream()
# prepares ONE shared SpatRaster per LOI up front and reuses that SAME
# object across every site x version job (76 here: 38 sites x 2 versions);
# each job's crop()/mask() then operates on the FULL project-wide extent
# every time before narrowing down to one catchment. Confirmed in practice
# to eventually crash with a terra "[readStart] file does not exist"
# error partway through a run (a temp-file GC race — the shared object
# survives many iterations of heavy repeated crop/mask churn on the full
# extent, which is exactly the failure mode that class of terra bug needs).
# `path_lazy` avoids this structurally, not just by luck: each site's LOI
# is pre-cropped down to just its own catchment ONCE, cached to a real
# permanent per-site file (hydroweight_loi/<loi>/<site_id>_reprojected.tif),
# *before* run_loi_attributes_stream() ever touches it — so there's far
# less data movement per job, and no long-lived shared object at all.
# CELESTE already uses `path_lazy` for its own LOIs without issue.
#
# NDVI (the "ndvi" entry below): SUPERSEDED 2026-08-26. Originally wired to
# /Users/sam/Downloads/NDVI_CAM_sites_BAP.tif, a raster built for the LAKE-
# site locations — that source genuinely doesn't cover SUD22 (~8 km outside
# its extent), so those columns read NA under it. But that's a bad reason
# to leave SUD22 uncovered: workflow/CAM/prepare_ndvi.R's
# prepare_cam_ndvi_site_rasters() already builds a properly stream-
# catchment-matched CONTINUOUS NDVI layer from the SAME per-catchment/group
# Google Earth Engine exports (data/ndvi/CAM/*.tif) that "ndvi_trend" below
# uses — confirmed directly that SUD22 maps to
# Landsat_NDVI_1984_2025_cam_sud17.tif (shared with SUD17) with full
# coverage. "ndvi" now uses that source instead of the lake-site raster —
# real coverage for every site "ndvi_trend" covers, SUD22 included.
#
# "ndvi_trend" is the SEPARATE Sen's-slope trend computed from the same
# export data (workflow/CAM/prepare_ndvi_trend.R, mirroring
# workflow/CELESTE/prepare_ndvi_trend.R) — kept alongside "ndvi" since the
# two answer different questions (a static mean/distribution vs. a
# directional trend), not because "ndvi" needed a different source.

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

# One-time (per distinct source file), cached to disk thereafter (a cheap
# clean-and-symlink, not the trend computation below) — see
# prepare_ndvi.R's header.
prepare_cam_ndvi_site_rasters(sites, output_dir = output_dir, cache_dir = cache_dir)

# One-time (per distinct source file), cached to disk thereafter — see
# prepare_ndvi_trend.R's header for the ~43 min total estimate on a fully
# cold cache (dominated by the largest merged-group files: SUD12+VER01,
# SUD17+SUD22, the NCMN cluster, ILD02).
prepare_cam_ndvi_trend_site_rasters(sites, output_dir = output_dir, cache_dir = cache_dir)

# One-time (single whole-domain group), cached to disk thereafter. Ontario
# MNR harvest/regen — real (not bbox-only) coverage over CAM's AOI
# confirmed before writing this; see prepare_harvest_regen.R's header.
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
    stats = c("distwtd_mean", "distwtd_sd", "mean", "sd", "median", "min", "max") # sum/cell_count meaningless for a trend, same exclusions CELESTE's ndvi_trend LOI uses
  ),
  list(
    path_lazy = fs::path(cache_dir, "hydroweight_loi", "harvest_regen", "{group_id}.tif"),
    name = "harvest_regen", type = "categorical", class_levels = harvest_regen_levels
  )
)

hw_results <- calculate_hydroweight_attributes_stream(
  sites = sites,
  group_manifest = group_manifest,
  output_dir = output_dir,
  cache_dir = cache_dir,
  loi_layers = loi_layers,
  catchment_versions = c("unclipped", "clipped"),
  raster_crs = config$working_crs # EPSG:3161 — not the function's default EPSG:3979
)

# Drop "clipped" rows that duplicate "unclipped" — see Stage 4's comment.
hw_results <- drop_redundant_clipped_rows(hw_results, upstream_results, site_col = "site")

write_csv(hw_results, here("output", PROJECT_ID, DELINEATION, paste0(PROJECT_ID, "_streams_hydroweight.csv")))
print(hw_results)

# =============================================================================
# STAGE 8 (optional) — Plotting-ready long tables from Stages 6-7's outputs
# =============================================================================
# Reshapes catchment_metrics.csv and CAM_streams_hydroweight.csv (both wide,
# one row per site x version) into 5 small, purpose-shaped long tables —
# one per LOI (canlcc/ndvi/ndvi_trend/harvest_regen) plus one for the
# morphometric metrics — rather than one wide row per site or one giant
# generic melt. The original wide CSVs are untouched and remain the right
# choice for statistical modelling (one row per site x version, one column
# per covariate); these are specifically for plotting (ggplot2 time series,
# composition, significance-filter plots, ...). See workflow/CAM/
# tidy_outputs.R's header for the exact output columns/shapes and the
# per-LOI year-handling rules.

source(here("workflow/CAM/tidy_outputs.R"))
tidy_cam_outputs(output_dir = output_dir)
