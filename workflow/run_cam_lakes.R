# run_cam_lakes.R
# =============================================================================
# Top-level runner for the CAM lake catchment delineation workflow, on the
# modular engine (workflow/R/engine/) — the standard approach for every
# project in this repo now. Supersedes an earlier, non-modular run on
# workflow/R/lake/'s own pipeline (retired — see git history), after
# validation: delineated catchments compared directly against that
# production run, 45/45 sites, IoU = 1.0000 and 0% area difference on
# every site (exact geometric match — OIH needs no breach step, so this
# migration avoided the MRDEM-breach nondeterminism CELESTE's own
# migration had to work through), and hydroweight canlcc values matched
# to within ~0.1% (normal raster-resampling noise).
#
# For correcting a single bad lake catchment: edit
# output/<site_id>/lake_pourpoint.tif in QGIS, then use
# rerun_engine_lake_site_watershed() (workflow/R/engine/99_rerun_sites) —
# same tool the old pipeline's rerun_watershed_lake() provided, ported to
# the engine before this pipeline was retired specifically so that
# capability wasn't lost.
#
# WHY THIS WORKS WITH SO LITTLE NEW CODE: the engine's lake-polygon
# delineation path already existed (workflow/R/engine/04_delineate_site.R's
# delineate_engine_lake_site(), adapted line-for-line from the old
# lake/03_delineate_lakes.R's delineate_single_lake() when the engine was
# first built) but had never actually been run before this. Terrain prep
# needed no new code either — workflow/CAM/run_cam_streams.R already runs
# OIH's Enhanced Flow Direction through engine/02_prepare_terrain.R's
# recode+accumulate+extract-streams path, identical physical data this
# project uses. So this runner was mostly wiring together already-built,
# already-proven pieces:
#
#   Stage 1 (match_lake_polygons) and Stage 4/6/7 (remove_upstream_lake_
#   catchments, calculate_catchment_metrics, calculate_hydroweight_
#   attributes) are workflow/R/lake/*.R functions, REUSED COMPLETELY
#   UNMODIFIED — none of them have any real dependency on which pipeline
#   produced their inputs, confirmed by reading each one's actual file
#   dependencies (not assumed):
#     - match_lake_polygons(): pure spatial/name matching, no cache_dir
#       dependency at all.
#     - calculate_hydroweight_attributes(): defaults to cache_dir/
#       dem_breached.tif + cache_dir/flow_accum.tif — BOTH already
#       materialized by engine/02_prepare_terrain.R (dem_breached.tif via
#       alias_dem_breached(), since OIH needs no real breaching — same
#       mechanism run_cam_streams.R already relies on). No d8_pntr.tif
#       dependency found in this file.
#     - remove_upstream_lake_catchments(): the ONE real gap — hardcodes
#       cache_dir/d8_pntr.tif (the old pipeline's own naming), but the
#       engine's terrain prep writes flow_pointer.tif instead. Fixed with
#       a cheap symlink alias right after Stage 2 — same pattern
#       engine/02_prepare_terrain.R's own alias_dem_breached() already
#       establishes for an analogous naming gap — rather than modifying
#       the shared function itself.
#
# WORKFLOW STAGES:
#   Stage 1 — Match lake sites to OHN waterbody polygons
#   Stage 2 — Prepare terrain (OIH DEM + Enhanced Flow Direction -> recode,
#              accumulate, extract streams) via the engine
#   Stage 3 — Delineate lake catchments (buffered lake polygon as pour
#              point) via the engine's lake-mode path
#   Stage 4 — Remove upstream lake catchments
#   Stage 5 — Re-clip rasters to clipped catchments
#   Stage 6 — Catchment morphometric metrics
#   Stage 7 — Distance-weighted catchment attributes (hydroweight)
#
# Usage: source from an R session in the project root.
# =============================================================================

library(sf)
library(terra)
library(whitebox)
library(hydroweight)
library(dplyr)
library(purrr)
library(readr)
library(tibble)
library(fs)
library(cli)
library(glue)
library(here)

source(here("workflow/R/utils.R"))
source(here("workflow/R/lake/01_match_lake_polygons.R"))
source(here("workflow/R/lake/04_remove_upstream_lakes.R"))
source(here("workflow/R/lake/05_hydroweight_attributes.R"))
source(here("workflow/R/07_reclip_outputs.R"))
source(here("workflow/R/08_catchment_metrics.R"))
source(here("workflow/R/stream/01_group_sites.R")) # unused directly (whole_domain grouping), but 01_build_group_manifest.R checks for it defensively
source(here("workflow/R/stream/05_delineate_sites.R")) # load_group_rasters() — used by delineate_engine_catchments()
source(here("workflow/R/engine/00_resolve_config.R"))
source(here("workflow/R/engine/01_build_group_manifest.R"))
source(here("workflow/R/engine/02_prepare_terrain.R"))
source(here("workflow/R/engine/03_prepare_streams_burn.R"))
source(here("workflow/R/engine/04_delineate_site.R"))

# =============================================================================
# CONFIGURATION
# =============================================================================

PROJECT_ID <- "CAM_lakes_engine"
output_dir <- here("output/CAM/lake_delineation_engine")
cache_dir  <- here("cache", PROJECT_ID)
fs::dir_create(output_dir, recurse = TRUE)
fs::dir_create(cache_dir, recurse = TRUE)

OIH_LAKES_PATH    <- "/Users/sam/Documents/cfs/shared_data/raw/hydro/waterbodies/ON_waterbodies/ohn_waterbodies_valid.gpkg"
OIH_DEM_PATH      <- "/Users/sam/Downloads/IntegratedHydrologyNE/EnforcedDEM.tif"
OIH_FLOW_DIR_PATH <- "/Users/sam/Downloads/IntegratedHydrologyNE/EnhancedFlowDirection.tif"

# OIH -> WhiteboxTools D8 encoding reclassification (one-step rotation) —
# same table run_cam_streams.R uses for the same physical data.
oih_recode_matrix <- matrix(
  c(128, 1, 1, 2, 2, 4, 4, 8, 8, 16, 16, 32, 32, 64, 64, 128),
  ncol = 2, byrow = TRUE
)

BUFFER_M       <- 10   # shoreline snap buffer for spatial lake matching (m)
NAME_DIST_MAX  <- 0.15 # max Jaro-Winkler distance for name matching
NAME_SEARCH_KM <- 40   # name match search radius (km)

LAKE_BUFFER_M        <- 30    # buffer applied to lake polygon pour point (m)
STREAM_THRESHOLD     <- 100   # flow accum threshold for stream extraction (cells)
MIN_UPSTREAM_AREA_M2 <- 10000 # minimum upstream lake area to include (m2, 1 ha)

# =============================================================================
# SITE DEFINITIONS — identical to workflow/run_cam_lakes.R (production);
# this is project data, not pipeline-specific, so it's unchanged here.
# =============================================================================

cam_sites_raw <- tibble::tribble(
  ~lake_name           , ~lat            , ~lon              ,
  "Chief"              , 46.362683       , -81.006352        ,
  "Arcand"             , 46.823374       , -80.292216        ,
  "Cucumber"           , 46.833425       , -80.316343        ,
  "Turtleshell"        , 46.886178       , -80.260631        ,
  "Upper Shakwa"       , 46.804235       , -81.969215        ,
  "Tee"                , 46.780054       , -81.913564        ,
  "Rushbrook"          , 46.745399       , -81.923317        ,
  "Nepahwin"           , 46.460605420467 , -80.9893512429419 ,
  "Daisy"              , 46.443509       , -80.900269        ,
  "George"             , 46.014566       , -81.406996        ,
  "Acid"               , 46.033666       , -81.444221        ,
  "Lumsden"            , 46.022442       , -81.428004        ,
  "Teardrop"           , 46.045789       , -81.414736        ,
  "Ruth-Roy"           , 46.09469        , -81.242431        ,
  "Wolf"               , 46.848869       , -80.643334        ,
  "Dewdney"            , 46.86242        , -80.635807        ,
  "Manitou"            , 46.825398       , -80.274342        ,
  "Bowland"            , 47.096559       , -80.846233        ,
  "Sam Martin"         , 46.879153       , -80.808045        ,
  "Bell"               , 46.123672       , -81.240847        ,
  "Tilton"             , 46.354884       , -81.068328        ,
  "Smoothwater"        , 47.42526        , -80.690605        ,
  "Marina"             , 47.394273       , -80.663699        ,
  "Whirligig"          , 47.380333       , -80.636416        ,
  "Aurora Whitepine"   , 47.394552       , -80.647887        ,
  "Low"                , 46.100279       , -81.566133        ,
  "Helen"              , 46.106352       , -81.561851        ,
  "Chiniguchi"         , 46.89681        , -80.67749         ,
  "Sunnywater"         , 47.398657       , -80.629523        ,
  "Lulu"               , 47.398399       , -80.748682        ,
  "Jerry"              , 47.362905       , -80.658052        ,
  "McCulloch"          , 47.317233       , -80.702609        ,
  "Marjorie"           , 46.923602       , -80.624567        ,
  "Laundrie"           , 47.142382       , -80.851195        ,
  "Paradise"           , 46.97395        , -80.759125        ,
  "Wawiashkashi"       , 46.780758       , -80.353666        ,
  "Baby"               , 46.464036       , -80.860249        ,
  "Clearwater"         , 46.378460       , -81.046010        ,
  "Lohi"               , 46.393304       , -81.036809        ,
  "Hannah"             , 46.442684       , -81.030408        ,
  "Middle"             , 46.443665       , -81.024360        ,
  "Swan"               , 46.365670       , -81.067828        ,
  "Sans Chambre"       , 46.722986       , -81.135476        ,
  "Little Whitepine"   , 47.378210       , -80.641877        ,
  "Whitepine (Mcleod)" , 47.281703       , -80.826524
)

sites <- cam_sites_raw |>
  dplyr::mutate(
    site_id = lake_name |>
      gsub("\\s+", "_", x = _) |>
      gsub("[^A-Za-z0-9_-]", "", x = _),
    site_name = lake_name
  ) |>
  dplyr::select(site_id, site_name, lake_name, lon, lat)

manual_id_lookup <- tibble::tibble(
  OGF_ID = c(
    115461047L, 120009986L, 117510443L, 117510487L,
    120010058L, 115460709L, 120010149L, 119924928L
  ),
  lake_name = c(
    "Marina", "Tilton", "Middle", "Hannah",
    "Lohi", "Jerry", "Daisy", "Chiniguchi"
  )
)

excluded_ogf_ids <- c(120010146L, 119924308L)

# =============================================================================
# STAGE 1 — Match lake sites to OHN polygons
# =============================================================================

polys_all <- match_lake_polygons(
  sites = sites,
  lakes_path = OIH_LAKES_PATH,
  buffer_m = BUFFER_M,
  name_dist_max = NAME_DIST_MAX,
  name_search_km = NAME_SEARCH_KM,
  manual_id_lookup = manual_id_lookup,
  excluded_ogf_ids = excluded_ogf_ids
)

cw_inform(glue::glue("Stage 1 complete. {nrow(polys_all)} lake polygon(s) matched."))

# =============================================================================
# CONFIG — resolve the engine run
# =============================================================================

run_config <- list(
  project_id = PROJECT_ID,
  output_dir = output_dir,
  cache_dir  = cache_dir,
  sites      = sites,

  dem = list(path = OIH_DEM_PATH),
  flow_direction = list(path = OIH_FLOW_DIR_PATH, recode = oih_recode_matrix),
  flow_pointer = NULL,
  crs = NULL, # NULL = OIH's own native CRS (EPSG:3161)
  stream_threshold = STREAM_THRESHOLD,

  streams_burn = list(source = "none"), # OIH is already hydrologically conditioned

  lake_polygons = polys_all,
  lake_buffer_m = LAKE_BUFFER_M,

  grouping = list(strategy = "whole_domain"), # OIH covers the whole domain as one project-wide raster set, same as CAM streams

  loi_layers = NULL # set in Stage 7 below
)

config <- resolve_engine_config(run_config)

# =============================================================================
# STAGE 1b — Group manifest
# =============================================================================

gm <- build_engine_group_manifest(config)
sites <- gm$sites
group_manifest <- gm$group_manifest
print(group_manifest)

# =============================================================================
# STAGE 2 — Terrain prep (recode OIH flow direction, accumulate, extract streams)
# =============================================================================

prepare_engine_terrain(config, group_manifest)

# Fix for remove_upstream_lake_catchments() (workflow/R/lake/
# 04_remove_upstream_lakes.R), which hardcodes cache_dir/d8_pntr.tif — the
# engine writes flow_pointer.tif instead. Symlink, not a copy — cheap, and
# terra::rast() reads through a symlink identically to a real file
# (confirmed elsewhere in this project). Same pattern engine/02_prepare_
# terrain.R's own alias_dem_breached() already uses for an analogous gap.
d8_pntr_alias <- fs::path(cache_dir, "d8_pntr.tif")
if (!fs::file_exists(d8_pntr_alias) && !fs::link_exists(d8_pntr_alias)) {
  fs::link_create(fs::path(cache_dir, "flow_pointer.tif"), d8_pntr_alias, symbolic = TRUE)
  cw_inform("Created d8_pntr.tif alias -> flow_pointer.tif (for remove_upstream_lake_catchments()).")
}

# =============================================================================
# STAGE 3 — Delineate lake catchments
# =============================================================================

results <- delineate_engine_catchments(
  config = config, sites = sites, group_manifest = group_manifest,
  min_cells = 10
)
print(results)

flagged <- dplyr::filter(results, flagged)
if (nrow(flagged) > 0) {
  cw_warn(glue::glue("\n{nrow(flagged)} site(s) flagged — review catchments:"))
  print(flagged[, c("site_id", "catchment_km2", "flag_reason")])
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
# STAGE 4 — Remove upstream lake catchments
# =============================================================================

results_stage4 <- remove_upstream_lake_catchments(
  sites = sites,
  polys_all = polys_all,
  lakes_path = OIH_LAKES_PATH,
  cache_dir = cache_dir,
  output_dir = output_dir,
  lake_buffer_m = LAKE_BUFFER_M,
  min_upstream_area_m2 = MIN_UPSTREAM_AREA_M2,
  exclude_waterbody_types = c("Pond")
)
print(results_stage4)

# =============================================================================
# STAGE 5 — Re-clip rasters to clipped catchments
# =============================================================================

results_stage5 <- reclip_outputs(
  sites = sites,
  output_dir = output_dir,
  cache_dir = cache_dir # lake path — no group_manifest needed
)
print(results_stage5)

# =============================================================================
# STAGE 6 — Morphometric metrics
# =============================================================================

metrics <- calculate_catchment_metrics(
  sites = sites,
  output_dir = output_dir,
  lake_polys = polys_all,
  streams_path = NULL # auto-detect per-site streams_clipped.tif
)
ref_table <- build_metrics_reference_table()
write_metrics_outputs(metrics = metrics, ref_table = ref_table, output_dir = output_dir)
print(metrics)

# =============================================================================
# STAGE 7 — Distance-weighted catchment attributes (hydroweight)
# =============================================================================
# Same LOIs as production run_cam_lakes.R: canlcc (land cover) and ndvi.

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

hw_results <- calculate_hydroweight_attributes(
  sites = sites,
  output_dir = output_dir,
  cache_dir = cache_dir,
  loi_layers = list(
    list(path = CANLCC_PATH, name = "canlcc", type = "categorical", class_levels = canlcc_levels),
    list(path = "/Users/sam/Downloads/NDVI_CAM_sites_BAP.tif", name = "ndvi", type = "continuous")
  ),
  lake_polys = polys_all,
  streams_path = NULL # auto-detect per-site streams_clipped.tif / streams.tif
)

write_csv(hw_results, fs::path(output_dir, "CAM_lakes_engine_hydroweight.csv"))
print(hw_results)

cw_inform(glue::glue(
  "\nHydroweighting complete: {nrow(hw_results)} row(s) written to ",
  "{output_dir}/CAM_lakes_engine_hydroweight.csv."
))
