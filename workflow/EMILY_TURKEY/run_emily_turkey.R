# run_emily_turkey.R
# =============================================================================
# Top-level runner for the EMILY_TURKEY stream-site catchment delineation
# workflow, on the modular, input-driven engine (workflow/R/engine/).
# Point pour points, OIH pre-conditioned terrain, EPSG:3161 — same shape as
# CAM streams (workflow/CAM/run_cam_streams.R) EXCEPT grouping: this
# project's 177 sites split across TWO disjoint OIH terrain tiles with
# zero overlap (24 sites only inside "IntegratedHydrologyNE", 153 only
# inside "IntegratedHydrologyNC" — confirmed directly, not assumed; CAM's
# sites happen to sit entirely inside the NE tile, which is why CAM never
# needed this). Uses grouping$strategy = "manual_groups" — a group per OIH
# tile, each with its own terrain source, site->group assignment
# auto-detected by which tile's raster actually has data at each
# coordinate. See workflow/EMILY_TURKEY/README.md for full rationale
# (source data, the HARM030 site_id disambiguation, and the discovery of
# the two-tile split).
#
# WORKFLOW STAGES:
#   Stage 1 — Resolve config, build group manifest (2 groups: "NE" and
#              "NC", one per OIH tile — grouping$strategy = "manual_groups")
#   Stage 2 — Prepare terrain per group (recode OIH flow direction to
#              WhiteboxTools encoding, derive flow accumulation + streams;
#              no breaching — OIH is already hydrologically conditioned)
#   Stage 3 — Delineate catchments (point pour point, Jenson snap)
#   Stage 4 — Remove upstream nested catchments
#   Stage 5 — Re-clip rasters/flowlines to clipped catchments
#   Stage 6 — Catchment morphometric metrics
#   Stage 7 — Distance-weighted catchment attributes (hydroweight) — CANLCC
#              land cover only for now; no NDVI/harvest LOI prep exists yet
#              for this project's AOI.
#
# CACHING: cache/EMILY_TURKEY/NE/ and cache/EMILY_TURKEY/NC/ (one
# subfolder per manual group — no DELINEATION subfolder above that; this
# project has only one delineation approach, per workflow/templates/
# run_engine_template.R's own guidance). Re-running skips any step whose
# output already exists.
#
# RERUNNING AFTER A CORRECTION (workflow/R/engine/99_rerun_sites):
# See run_cam_streams.R's header for the full two-scenario walkthrough
# (raw-coordinate fix vs. manually-edited snap fix) — rerun_engine_sites()
# works identically here; it doesn't need to know which manual group a
# site belongs to, since that's re-derived from group_manifest either way.
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
source(here("workflow/R/stream/burn_streams.R")) # burn_streams_into_dem() — reused by engine 03 (not exercised: streams_burn source = "none" here)
source(here("workflow/R/stream/delineate_sites.R")) # load_group_rasters(), snap_pour_point(), delineate_watershed(), clip_rasters_to_catchment(), clip_flowlines_to_catchment() — reused unmodified by engine 04
source(here("workflow/R/stream/hydroweight_attributes.R")) # calculate_hydroweight_attributes_stream() — reused unmodified, raster_crs overridden below
source(here("workflow/R/remove_upstream.R"))
source(here("workflow/R/reclip_outputs.R"))
source(here("workflow/R/catchment_metrics.R"))
source(here("workflow/R/engine/00_resolve_config.R"))
source(here("workflow/R/engine/01_build_group_manifest.R"))
source(here("workflow/R/engine/02_prepare_terrain.R"))
source(here("workflow/R/engine/03_prepare_streams_burn.R"))
source(here("workflow/R/engine/04_delineate_site.R"))
source(here("code/reset_workflow.R")) # reset_site() — used by rerun_engine_sites()'s resnap path
source(here("workflow/R/engine/99_rerun_sites")) # rerun_engine_site_watershed(), rerun_engine_sites() — see run_cam_streams.R's header

# =============================================================================
# CONFIGURATION — verify paths before running
# =============================================================================

PROJECT_ID <- "EMILY_TURKEY"

output_dir <- here("output", PROJECT_ID)
cache_dir  <- here("cache", PROJECT_ID)
fs::dir_create(output_dir, recurse = TRUE)
fs::dir_create(cache_dir, recurse = TRUE)

# OIH Enforced DEM and Enhanced Flow Direction — same physical files
# run_cam_lakes.R / run_cam_streams.R use (kept untouched by this file).
OIH_NE_DEM_PATH      <- "/Users/sam/Downloads/IntegratedHydrologyNE/EnforcedDEM.tif"
OIH_NE_FLOW_DIR_PATH <- "/Users/sam/Downloads/IntegratedHydrologyNE/EnhancedFlowDirection.tif"
OIH_NC_DEM_PATH      <- "/Users/sam/Downloads/IntegratedHydrologyNC/EnforcedDEM.tif"
OIH_NC_FLOW_DIR_PATH <- "/Users/sam/Downloads/IntegratedHydrologyNC/EnhancedFlowDirection.tif"

# OIH -> WhiteboxTools D8 encoding reclassification (one-step rotation),
# same table as run_cam_streams.R / workflow/R/lake/02_prepare_oih_dem.R —
# generic in the engine, applied here as this project's own config value:
#   OIH 128 (NE) -> WBT 1     OIH 1  (E)  -> WBT 2
#   OIH 2   (SE) -> WBT 4     OIH 4  (S)  -> WBT 8
#   OIH 8   (SW) -> WBT 16    OIH 16 (W)  -> WBT 32
#   OIH 32  (NW) -> WBT 64    OIH 64 (N)  -> WBT 128
oih_recode_matrix <- matrix(
  c(128, 1, 1, 2, 2, 4, 4, 8, 8, 16, 16, 32, 32, 64, 64, 128),
  ncol = 2, byrow = TRUE
)

STREAM_THRESHOLD <- 100 # flow accum threshold for stream extraction (cells) — same default as CAM streams (same OIH terrain)
SNAP_DIST_M      <- 200 # point-mode pour point snap distance (m)
MIN_CELLS        <- 10  # minimum catchment size before flagging

# =============================================================================
# SITE DEFINITIONS
# =============================================================================
# Converted from /Users/sam/Downloads/DataCodes.csv, via
# data/emily_turkey_sites_raw.csv (177 sites: N/W decimal-degree WGS84
# coordinates, "Sample Site" as the raw site identifier, kept unmodified in
# the raw CSV so the sanitizing/disambiguation below stays visible here
# rather than being hidden in a pre-cleaned file).

sites_raw <- readr::read_csv(here("data/emily_turkey_sites_raw.csv"), show_col_types = FALSE)

sites <- sites_raw |>
  dplyr::mutate(
    site_id = sample_site |>
      gsub("'", "", x = _) |>
      gsub("\\s+", "_", x = _) |>
      gsub("[^A-Za-z0-9_-]", "", x = _),
    # "HARM030'" (sample_code 149, 46.889617/-84.29985, 29-07-2019) collapses
    # onto "HARM030" (sample_code 148, 46.81364/-84.36151, 26-07-2019) after
    # the sanitize above — these are two genuinely distinct sites (different
    # coordinates and collection dates), not a duplicate entry. Disambiguate
    # by appending the Project Sample Code to the second row only.
    site_id = dplyr::if_else(
      sample_site == "HARM030'" & sample_code == 149,
      paste0(site_id, "_", sample_code),
      site_id
    ),
    site_name = sample_site
  ) |>
  dplyr::select(site_id, site_name, lon, lat, location_project, date_collected)

stopifnot(!anyDuplicated(sites$site_id)) # defensive: resolve_engine_config() would also catch this, but fails with a clearer message here, at the point of derivation

# =============================================================================
# STAGE 1 — Resolve config, build group manifest
# =============================================================================

run_config <- list(
  project_id = "EMILY_TURKEY", # internal label — group_ids come from grouping$groups below, not this
  output_dir = output_dir,
  cache_dir  = cache_dir,
  sites      = sites,

  # NOTE: top-level dem/flow_direction/flow_pointer stay NULL for
  # "manual_groups" — each group below supplies its own terrain instead of
  # one shared globally (see resolve_engine_config()'s manual_groups guard).
  dem = NULL,
  flow_direction = NULL,
  flow_pointer = NULL,
  flow_accum   = NULL,
  crs = NULL, # NULL = match the first manual group's ("NE") native CRS (EPSG:3161) — both tiles are already EPSG:3161, so no reprojection happens for either
  stream_threshold = STREAM_THRESHOLD,

  streams_burn = list(source = "none"), # OIH flow direction is already conditioned — burning would have to happen before conditioning, not after

  lake_polygons = NULL, # point pour point mode
  lake_buffer_m = 30,

  # "manual_groups": this project's 177 sites split across two disjoint
  # OIH terrain tiles with zero overlap (confirmed directly — see
  # README.md). Each group supplies its own dem/flow_direction; site ->
  # group assignment is auto-detected by which tile's raster actually has
  # data at each site's coordinate (build_engine_group_manifest() aborts
  # loudly if any site matches zero or more than one group).
  grouping = list(
    strategy = "manual_groups",
    groups = list(
      list(
        group_id = "NE",
        dem = list(path = OIH_NE_DEM_PATH),
        flow_direction = list(path = OIH_NE_FLOW_DIR_PATH, recode = oih_recode_matrix)
      ),
      list(
        group_id = "NC",
        dem = list(path = OIH_NC_DEM_PATH),
        flow_direction = list(path = OIH_NC_FLOW_DIR_PATH, recode = oih_recode_matrix)
      )
    )
  ),

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
# hillshade.tif in cache/EMILY_TURKEY/.

prepare_engine_terrain(config, group_manifest)

# Spot-check: load streams.tif and dem.tif in QGIS before proceeding. If
# streams look wrong, adjust STREAM_THRESHOLD and re-run Stage 2 (delete
# cache/EMILY_TURKEY/streams.tif first).

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
# STAGE 6 — Catchment morphometric metrics
# =============================================================================

metrics <- calculate_catchment_metrics(sites = sites, output_dir = output_dir)
# Drop "clipped" rows that duplicate "unclipped" (no nested upstream
# catchment actually erased for that site) — see Stage 4's comment.
metrics <- drop_redundant_clipped_rows(metrics, upstream_results, site_col = "site_id")
ref_table <- build_metrics_reference_table()
write_metrics_outputs(metrics = metrics, ref_table = ref_table, output_dir = output_dir)
print(metrics)

# =============================================================================
# STAGE 7 — Distance-weighted catchment attributes (hydroweight)
# =============================================================================
# CANLCC land cover only, for now — no NDVI/harvest-regen LOI prep exists
# yet for this project's AOI. Same source raster + class lookup table as
# CAM streams/CAM lakes (workflow/CAM/run_cam_streams.R).

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
  raster_crs = config$working_crs # EPSG:3161 — not the function's default EPSG:3979
)

# Drop "clipped" rows that duplicate "unclipped" — see Stage 4's comment.
hw_results <- drop_redundant_clipped_rows(hw_results, upstream_results, site_col = "site")

write_csv(hw_results, here("output", PROJECT_ID, paste0(PROJECT_ID, "_hydroweight.csv")))
print(hw_results)
