# run_cam_lakes.R
# =============================================================================
# Top-level runner for the CAM lake catchment delineation workflow.
#
# USAGE:
#   1. Verify the paths in the CONFIGURATION section below
#   2. Source this file or run sections interactively in R
#
# WORKFLOW STAGES:
#   Stage 1 — Match lake sites to OHN waterbody polygons
#   Stage 2 — Prepare OIH hydrological products (DEM, flow pointer,
#              flow accumulation, stream network)
#   Stage 3 — Delineate lake catchments
#   Stage 4 — Remove upstream lake catchments (clipping)
#   Stage 5 — Re-clip rasters and flowlines to clipped catchments
#   Stage 6 — Compute morphometric metrics
#   Stage 7 — Distance-weighted catchment attributes (hydroweight)
#
# CACHING:
#   Project-level rasters (DEM, D8 pointer, flow products) are cached in
#   cache/CAM/. Re-running skips any step whose output already exists.
#   Per-site outputs are in output/CAM/lake_delineation/<site_id>/.
#
# CORRECTING A CATCHMENT:
#   If a catchment is wrong, the lake pour point raster may need editing.
#   Each site has a lake_pourpoint.tif in its output folder. Edit it in
#   QGIS (add or extend cells to capture the correct lake boundary), then:
#     rerun_watershed_lake("site_id", cache_dir, output_dir)
#   This re-runs watershed delineation from the edited pour point without
#   regenerating the raster from the lake polygon.
#
# RESETTING:
#   Delete output/CAM/lake_delineation/<site_id>/ to force a single site
#   to rerun from Stage 3.
#   Delete cache/CAM/ contents to force Stages 2–6 to rerun fully.
# =============================================================================

# -- Packages -----------------------------------------------------------------
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

# -- Source modules -----------------------------------------------------------
source(here("workflow/R/utils.R"))
source(here("workflow/R/lake/01_match_lake_polygons.R"))
source(here("workflow/R/lake/02_prepare_oih_dem.R"))
source(here("workflow/R/lake/03_delineate_lakes.R"))
source(here("workflow/R/lake/04_remove_upstream_lakes.R"))
source(here("workflow/R/lake/05_hydroweight_attributes.R"))
source(here("workflow/R/07_reclip_outputs.R"))
source(here("workflow/R/08_catchment_metrics.R"))

# =============================================================================
# CONFIGURATION — verify paths before running
# =============================================================================

PROJECT_ID <- "CAM"
DELINEATION <- "lake_delineation"

output_dir <- here("output", PROJECT_ID, DELINEATION)
cache_dir <- here("cache", PROJECT_ID)

fs::dir_create(output_dir, recurse = TRUE)
fs::dir_create(cache_dir, recurse = TRUE)

# OHN waterbody polygon layer
OIH_LAKES_PATH <- "/Users/sam/Documents/cfs/shared_data/raw/hydro/waterbodies/ON_waterbodies/ohn_waterbodies_valid.gpkg"

# OIH Enforced DEM and Enhanced Flow Direction
OIH_DEM_PATH <- "/Users/sam/Downloads/IntegratedHydrologyNE/EnforcedDEM.tif"
OIH_FLOW_DIR_PATH <- "/Users/sam/Downloads/IntegratedHydrologyNE/EnhancedFlowDirection.tif"

# Lake matching parameters
BUFFER_M <- 10 # shoreline snap buffer for spatial matching (m)
NAME_DIST_MAX <- 0.15 # max Jaro-Winkler distance for name matching
NAME_SEARCH_KM <- 40 # name match search radius (km)

# Delineation parameters
LAKE_BUFFER_M <- 30 # buffer applied to lake polygon pour point (m)
STREAM_THRESHOLD <- 100 # flow accum threshold for stream extraction (cells)
MIN_UPSTREAM_AREA_M2 <- 10000 # minimum upstream lake area to include (m², 1 ha)

# =============================================================================
# SITE DEFINITIONS
# =============================================================================
# site_id must be alphanumeric/underscore/hyphen (used as the output folder name)
# lake_name is matched against OFFICIAL_N in the OHN waterbody dataset

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

# Derive filesystem-safe site_id from lake_name:
#   spaces → underscores, remove characters outside [A-Za-z0-9_-]
sites <- cam_sites_raw |>
  dplyr::mutate(
    site_id = lake_name |>
      gsub("\\s+", "_", x = _) |>
      gsub("[^A-Za-z0-9_-]", "", x = _),
    site_name = lake_name
  ) |>
  dplyr::select(site_id, site_name, lake_name, lon, lat)

# Manual OGF_ID overrides — sites that spatial and name matching miss
# (verified against OHN data; add new entries here if needed)
manual_id_lookup <- tibble::tibble(
  OGF_ID = c(
    115461047L,
    120009986L,
    117510443L,
    117510487L,
    120010058L,
    115460709L,
    120010149L,
    119924928L
  ),
  lake_name = c(
    "Marina",
    "Tilton",
    "Middle",
    "Hannah",
    "Lohi",
    "Jerry",
    "Daisy",
    "Chiniguchi"
  )
)

# OGF_IDs to exclude from spatial and name matching (known wrong polygons)
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

# Save matched polygons for inspection in QGIS before proceeding
terra::writeVector(
  polys_all,
  here("data", paste0(PROJECT_ID, "_lake_polygons.gpkg")),
  overwrite = TRUE
)

message(glue(
  "\nStage 1 complete. Inspect data/{PROJECT_ID}_lake_polygons.gpkg in QGIS ",
  "to verify polygon matches before proceeding to Stage 2.\n"
))

# =============================================================================
# STAGE 2 — Prepare OIH hydrological products
# =============================================================================
# Caches dem.tif, d8_pntr.tif, flow_accum.tif, streams.tif in cache/CAM/.
# Each step is skipped if the output already exists. Adjust STREAM_THRESHOLD
# if the stream network looks over- or under-extracted in QGIS.

prepare_oih_products(
  cache_dir = cache_dir,
  oih_dem_path = OIH_DEM_PATH,
  oih_flow_dir_path = OIH_FLOW_DIR_PATH,
  stream_threshold = STREAM_THRESHOLD
)

# Spot-check: load streams and DEM in QGIS (cache/CAM/streams.tif,
# cache/CAM/dem.tif) before proceeding. If streams look wrong, adjust
# STREAM_THRESHOLD and re-run Stage 2 (delete cache/CAM/streams.tif first).

# =============================================================================
# STAGE 3 — Delineate lake catchments
# =============================================================================
# Each site gets a folder output/CAM/lake_delineation/<site_id>/ containing:
#   catchment.gpkg      — catchment polygon
#   watershed.tif       — binary watershed raster
#   lake_pourpoint.tif  — buffered lake raster used as pour point (editable)
#   dem.tif             — OIH DEM clipped to catchment
#   flow_pointer.tif    — D8 pointer clipped to catchment
#   flow_accum.tif      — flow accumulation clipped to catchment
#   streams.tif         — stream network clipped to catchment

results_stage3 <- delineate_lake_catchments(
  sites = sites,
  polys_all = polys_all,
  cache_dir = cache_dir,
  output_dir = output_dir,
  lake_buffer_m = LAKE_BUFFER_M
)

print(results_stage3)

# Combine all catchments for QA in QGIS
all_catchments <- purrr::map(sites$site_id, function(sid) {
  path <- fs::path(site_output_dir(output_dir, sid), "catchment.gpkg")
  if (!cache_exists(path)) {
    return(NULL)
  }
  sf::st_read(path, quiet = TRUE) |> dplyr::mutate(site_id = sid)
}) |>
  purrr::compact() |>
  dplyr::bind_rows()

sf::st_write(
  all_catchments,
  here("output", PROJECT_ID, DELINEATION, "all_catchments.gpkg"),
  delete_dsn = TRUE,
  quiet = TRUE
)

# Check flagged sites before proceeding
flagged <- dplyr::filter(results_stage3, flagged)
if (nrow(flagged) > 0) {
  message(glue(
    "\n{nrow(flagged)} site(s) flagged — inspect catchments in QGIS:"
  ))
  print(flagged[, c("site_id", "catchment_km2", "flag_reason")])
}

# --- Correcting a catchment ---
# If a catchment does not fully contain the lake, edit lake_pourpoint.tif
# in QGIS (paint additional cells around the lake boundary), then:
#   rerun_watershed_lake("site_id", cache_dir, output_dir)

# =============================================================================
# STAGE 4 — Remove upstream lake catchments
# =============================================================================
# Finds lakes upstream of each target site within its catchment, delineates
# their catchments, and erases them from the target catchment polygon.
# Results saved as catchment_clipped.gpkg per site.

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

# Inspect output/CAM/lake_delineation/all_catchments_clipped.gpkg in QGIS
# to verify upstream clipping looks correct before proceeding.

# =============================================================================
# STAGE 5 — Re-clip rasters to clipped catchments
# =============================================================================
# Reads catchment_clipped.gpkg for each site and clips all group rasters to
# that boundary. Writes <name>_clipped.tif files alongside the originals.

results_stage5 <- reclip_outputs(
  sites = sites,
  output_dir = output_dir,
  cache_dir = cache_dir # lake path — no group_manifest needed
)

print(results_stage5)

# =============================================================================
# STAGE 6 — Morphometric metrics
# =============================================================================
# Computes geometry, areal, and relief metrics for both the unclipped and
# clipped catchments. Writes catchment_metrics.csv and a reference table.

metrics <- calculate_catchment_metrics(
  sites = sites,
  output_dir = output_dir,
  lake_polys = polys_all, # enables lake metric category
  streams_path = NULL # auto-detect per-site streams_clipped.tif
)
ref_table <- build_metrics_reference_table()
write_metrics_outputs(
  metrics = metrics,
  ref_table = ref_table,
  output_dir = here("output", PROJECT_ID, DELINEATION)
)

print(metrics)

# =============================================================================
# STAGE 7 — Distance-weighted catchment attributes (hydroweight)
# =============================================================================
# Computes inverse distance-weighted land cover and/or continuous variable
# attributes for each CAM lake catchment.
#
# loi_layers is a list of raster descriptors. Each element must have:
#   path         — path to the raster file (single or multi-layer)
#   name         — short label used as column prefix in the output table
#   type         — "categorical" (proportions) or "continuous" (statistics)
#   class_levels — data.frame(ID, Class) for readable column names (categorical)
#   stats        — character vector of stats to keep; NULL = keep all (continuous)
#
# streams_path = NULL auto-detects per-site stream rasters. Set NA to disable
# stream-based weighting schemes entirely.

CANLCC_PATH <- "/Users/sam/Documents/cfs/shared_data/raw/landcover/CAN_LLC_2020.tif"

canlcc_levels <- data.frame(
  stringsAsFactors = FALSE,
  ID = c(1L, 2L, 5L, 6L, 8L, 10L, 11L, 12L, 13L, 14L, 15L, 16L, 17L, 18L, 19L),
  Class = c(
    "needleleaf_forest",
    "taiga_needleleaf_forest",
    "broadleaf_deciduous_forest",
    "mixed_forest",
    "temperate_shrubland",
    "temperate_grassland",
    "shrubland_lichen_moss",
    "grassland_lichen_moss",
    "barren_lichen_moss",
    "wetland",
    "cropland",
    "barren",
    "urban",
    "water",
    "snow_ice"
  )
)

hw_results <- calculate_hydroweight_attributes(
  sites = sites,
  output_dir = output_dir,
  cache_dir = cache_dir,
  loi_layers = list(
    list(
      path = CANLCC_PATH,
      name = "canlcc",
      type = "categorical",
      class_levels = canlcc_levels
    ),
    list(
      path = "/Users/sam/Downloads/NDVI_CAM_sites_BAP.tif",
      name = "ndvi",
      type = "continuous"
    )
    # Add further LOI layers here, e.g.:
    # list(path = "/path/to/ndvi_stack.tif", name = "ndvi", type = "continuous")
  ),
  lake_polys = polys_all,
  streams_path = NULL # auto-detect per-site streams_clipped.tif / streams.tif
)
view(hw_results)

hw_results <- hw_results |>
  select(site, contains("canlcc"))

write_csv(
  hw_results,
  here(
    "output",
    PROJECT_ID,
    DELINEATION,
    paste0(PROJECT_ID, "_hydroweight.csv")
  )
)

print(hw_results)
