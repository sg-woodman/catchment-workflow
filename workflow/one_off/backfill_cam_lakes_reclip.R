# backfill_cam_lakes_reclip.R
# =============================================================================
# One-off backfill: CAM lakes' Stage 5 (reclip_outputs()) has never actually
# clipped a single raster for any site — its lake branch called
# load_lake_cache_rasters(cache_dir), a function that was never defined
# anywhere in the repo. Every call threw "could not find function", caught
# by reclip_site()'s own tryCatch and downgraded to a per-site warning, so
# run_cam_lakes.R "completed successfully" while silently producing zero
# *_clipped.tif rasters for all 45 sites, every run. Fixed in
# workflow/R/reclip_outputs.R (load_lake_cache_rasters() now defined). Found
# 2026-08-31 while checking whether the CAM output/cache cleanup could
# repoint catchment_delineation_guide.qmd's dem_clipped.tif figure at the
# engine's own output — see [[celeste-output-cache-cleanup]] memory.
#
# Since the bug was in Stage 5 only, Stage 1-4 outputs (sites, polys_all,
# catchment.gpkg, catchment_clipped.gpkg, all_catchments_clipped.gpkg) are
# already correct and cached — this script reconstructs the in-memory state
# Stage 5-7 need (Stage 1 + 1b, cheap: spatial lake matching + config/group
# manifest resolution, no DEM work) and reruns Stage 5 (now fixed), Stage 6
# (metrics — the "clipped" version rows were reading a dem_clipped.tif that
# never existed), and Stage 7 (hydroweight — the "clipped" version was
# silently missing its stream-based weighting schemes for the same reason)
# fresh. Mirrors run_cam_lakes.R's own Stage 1/1b/5/6/7 exactly; skips
# Stage 2/3/4 since their outputs are unaffected and already on disk.
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
source(here("workflow/R/lake/match_lake_polygons.R"))
source(here("workflow/R/lake/remove_upstream_lakes.R"))
source(here("workflow/R/lake/hydroweight_attributes.R"))
source(here("workflow/R/reclip_outputs.R")) # load_lake_cache_rasters() — the fix
source(here("workflow/R/catchment_metrics.R"))
source(here("workflow/R/stream/group_sites.R"))
source(here("workflow/R/stream/delineate_sites.R"))
source(here("workflow/R/engine/00_resolve_config.R"))
source(here("workflow/R/engine/01_build_group_manifest.R"))
source(here("workflow/R/engine/02_prepare_terrain.R"))
source(here("workflow/R/engine/03_prepare_streams_burn.R"))
source(here("workflow/R/engine/04_delineate_site.R"))

# =============================================================================
# CONFIGURATION — identical to run_cam_lakes.R
# =============================================================================

PROJECT_ID <- "CAM_lakes_engine"
output_dir <- here("output/CAM/lake_delineation_engine")
cache_dir  <- here("cache", PROJECT_ID)

OIH_LAKES_PATH    <- "/Users/sam/Documents/cfs/shared_data/raw/hydro/waterbodies/ON_waterbodies/ohn_waterbodies_valid.gpkg"
OIH_DEM_PATH      <- "/Users/sam/Downloads/IntegratedHydrologyNE/EnforcedDEM.tif"
OIH_FLOW_DIR_PATH <- "/Users/sam/Downloads/IntegratedHydrologyNE/EnhancedFlowDirection.tif"

oih_recode_matrix <- matrix(
  c(128, 1, 1, 2, 2, 4, 4, 8, 8, 16, 16, 32, 32, 64, 64, 128),
  ncol = 2, byrow = TRUE
)

BUFFER_M       <- 10
NAME_DIST_MAX  <- 0.15
NAME_SEARCH_KM <- 40

LAKE_BUFFER_M        <- 30
STREAM_THRESHOLD     <- 100
MIN_UPSTREAM_AREA_M2 <- 10000

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
# CONFIG — resolve the engine run (identical to run_cam_lakes.R)
# =============================================================================

run_config <- list(
  project_id = PROJECT_ID, output_dir = output_dir, cache_dir = cache_dir, sites = sites,
  dem = list(path = OIH_DEM_PATH),
  flow_direction = list(path = OIH_FLOW_DIR_PATH, recode = oih_recode_matrix),
  flow_pointer = NULL, crs = NULL, stream_threshold = STREAM_THRESHOLD,
  streams_burn = list(source = "none"),
  lake_polygons = polys_all, lake_buffer_m = LAKE_BUFFER_M,
  grouping = list(strategy = "whole_domain"),
  loi_layers = NULL
)
config <- resolve_engine_config(run_config)

# =============================================================================
# STAGE 1b — Group manifest
# =============================================================================

gm <- build_engine_group_manifest(config)
sites <- gm$sites
group_manifest <- gm$group_manifest

cw_inform(glue::glue("Reconstructed state: {nrow(sites)} site(s). Skipping Stage 2-4 (unaffected, already on disk)."))

# =============================================================================
# STAGE 5 — Re-clip rasters to clipped catchments (THE FIX)
# =============================================================================

results_stage5 <- reclip_outputs(
  sites = sites,
  output_dir = output_dir,
  cache_dir = cache_dir
)
print(results_stage5)
print(table(results_stage5$status))

# =============================================================================
# STAGE 6 — Morphometric metrics
# =============================================================================

metrics <- calculate_catchment_metrics(
  sites = sites,
  output_dir = output_dir,
  lake_polys = polys_all,
  streams_path = NULL
)
ref_table <- build_metrics_reference_table()
write_metrics_outputs(metrics = metrics, ref_table = ref_table, output_dir = output_dir)
print(metrics)

# =============================================================================
# STAGE 7 — Distance-weighted catchment attributes (hydroweight)
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

hw_results <- calculate_hydroweight_attributes(
  sites = sites,
  output_dir = output_dir,
  cache_dir = cache_dir,
  loi_layers = list(
    list(path = CANLCC_PATH, name = "canlcc", type = "categorical", class_levels = canlcc_levels),
    list(path = "/Users/sam/Downloads/NDVI_CAM_sites_BAP.tif", name = "ndvi", type = "continuous")
  ),
  lake_polys = polys_all,
  streams_path = NULL
)

write_csv(hw_results, fs::path(output_dir, "CAM_lakes_engine_hydroweight.csv"))
print(hw_results)

cw_inform(glue::glue(
  "\nBackfill complete: {nrow(hw_results)} hydroweight row(s) written to ",
  "{output_dir}/CAM_lakes_engine_hydroweight.csv."
))
