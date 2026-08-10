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
source(here("workflow/raster_attributes.R")) # reclassify_categorical(), sens_slope_trend() — used ad hoc to prepare LOI rasters below, not part of Stages 1-9

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
# NDVI (data/ndvi/*.tif): 12 regional tiles (COC, KEN1-4, MOR, NIP1-5, TUR),
# 1984-2025 annual composites, all already EPSG:3979/30m/42-band-aligned
# (verified with check_tile_consistency() — no reprojection or band
# alignment needed). Mosaicked into one VRT via build_mosaic_vrt(), which
# sets vrtnodata so areas outside every tile read as NA, not 0 — without
# that fix a gap reads as valid-looking zeros indistinguishable from real
# data. Used here via $path_lazy rather than $path because the 12 tiles'
# combined bounding box is transcontinental (CELESTE sites span NWT,
# Saskatchewan, Ontario, Quebec, Newfoundland) even though actual coverage
# is small clustered patches — caching that whole extent as one file would
# be enormous for no benefit.
#
# KNOWN GAP: the NBE group (20 sites) has no NDVI tile at all — confirmed
# empirically (every NBE site's LOI crop is 100% NA after the vrtnodata
# fix). Those sites are silently excluded from hw_results below (the
# existing "all NA after crop/mask" check), not an error. Add an NBE tile
# to data/ndvi/ once available and rerun.
#
# To reclassify a raw categorical raster or derive a Sen's-slope trend
# raster from a time-series stack before adding it below, see
# workflow/raster_attributes.R (reclassify_categorical(), sens_slope_trend()).

ndvi_vrt_path <- here("cache", PROJECT_ID, "hydroweight_loi", "ndvi_mosaic.vrt")
dir_create(fs::path_dir(ndvi_vrt_path))
build_mosaic_vrt(
  files = list.files(here("data/ndvi"), pattern = "[.]tif$", full.names = TRUE),
  vrt_path = ndvi_vrt_path
)

# Harvest/regen disturbance (Ontario only: shared_data/raw/harvest/
# ontario_harvest.gdb). Only NIP, TUR, and KEN sites fall inside Ontario's
# harvest-tracked extent (confirmed by real spatial-filter feature counts,
# not just a bounding-box check) — COC/MOR/NBE sites get no coverage and
# are silently excluded downstream (the same "all NA after crop/mask"
# path as the NDVI/NBE gap).
#
# Rasterized once per group (not once for all of Ontario — like the NDVI
# VRT, the combined extent across groups is far larger than actual site
# coverage) via rasterize_competing_classes(): one band per AR_YEAR present
# in the group's data, plus one "combined" (all-years) band. Two competing
# classes, most-recent-AR_YEAR-wins where they'd otherwise overlap (ties —
# which is EVERY within-year overlap, since a per-year band only ever
# includes same-year features from both classes by construction — go to
# whichever bucket is named LAST in `buckets` below; regen is last, so
# regen wins ties, matching "harvested then regenerated -> regen" as the
# more ecologically current state).
#
# CC-only for now (Harvest_CC02 + Harvest_CC17 = the full 2002-2024 clear-
# cut record, two non-overlapping vintage exports of the same product).
# To include seed-tree/selection and shelterwood harvest too, add
# Harvest_SE02/17 and Harvest_SH02/17 to the harvest bucket below — no
# other code changes needed.
harvest_regen_groups <- c("NIP", "TUR", "KEN")
ontario_harvest_gdb <- "~/Documents/cfs/shared_data/raw/harvest/ontario_harvest.gdb"

for (grp in harvest_regen_groups) {
  group_raster_path <- here(
    "cache", PROJECT_ID, "hydroweight_loi", "harvest_regen", paste0(grp, ".tif")
  )
  if (cache_exists(group_raster_path)) next

  dir_create(fs::path_dir(group_raster_path))
  grp_aoi <- group_manifest[group_manifest$group_id == grp, ]
  grp_template <- rast(here("cache", PROJECT_ID, grp, "dem_breached.tif"))

  hr <- rasterize_competing_classes(
    buckets = list(
      harvest = c("Harvest_CC02", "Harvest_CC17"),
      regen   = c("Regen_Seed", "Regen_Natural", "Regen_Plant")
    ),
    gdb_path = path.expand(ontario_harvest_gdb),
    template = grp_template,
    crop_to  = grp_aoi
  )
  writeRaster(
    hr,
    group_raster_path,
    overwrite = TRUE,
    datatype = "INT1U",
    gdal = "PHOTOMETRIC=MINISBLACK" # avoid a real GDAL quirk: an
    # exactly-3-band Byte GeoTIFF gets auto-tagged as RGB, silently
    # renaming bands to red/green/blue on every subsequent read
  )
}

harvest_regen_levels <- data.frame(ID = 0:2, Class = c("other", "harvest", "regen"))

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
    path_lazy = ndvi_vrt_path,
    name = "ndvi",
    type = "continuous"
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
