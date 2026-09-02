# run_<PROJECT>.R  [TEMPLATE — copy this file, fill in the TODOs, rename]
# =============================================================================
# Top-level runner template for a new project on the modular, input-driven
# delineation engine (workflow/R/engine/). Delineation runs off whatever
# inputs YOU supply — only a sites table + a DEM (or a pre-conditioned
# flow direction / flow pointer) is required. Everything else (breach,
# flow pointer, flow accumulation, burn-in streams, grouping) is either
# taken from what you provide, or generated with a documented default.
#
# HOW TO USE THIS TEMPLATE:
#   1. Copy to workflow/<PROJECT>/run_<project>.R
#   2. Fill in every TODO in the CONFIGURATION section
#   3. Decide: point pour points (leave lake_polygons = NULL) or lake
#      polygons (supply lake_polygons + sites$lake_name)?
#   4. Decide: does your terrain source need HydroBasins-based grouping
#      (a large national mosaic/VRT, like MRDEM — crop per-region) or is
#      it already one project-wide raster (like OIH — grouping$strategy =
#      "whole_domain", the default)?
#   5. Source this file or run sections interactively in R
#
# See an existing project's own runner under workflow/<PROJECT>/ for a
# filled-in example (e.g. a pre-conditioned flow direction source, point
# pour points, whole_domain grouping).
#
# WORKFLOW STAGES:
#   Stage 1 — Resolve config, build group manifest
#   Stage 2 — Prepare terrain (conditional: breach+pointer from a raw DEM,
#              OR recode+use an already-conditioned flow direction/pointer
#              if supplied)
#   Stage 3 — Delineate catchments (point pour point, or lake polygon if
#              lake_polygons is supplied)
#   Stage 4 — Remove upstream nested catchments
#   Stage 5 — Re-clip rasters/flowlines to clipped catchments
#   Stage 6 (optional) — Catchment morphometric metrics
#   Stage 7 (optional) — Distance-weighted catchment attributes (hydroweight)
#
# CACHING: every output is checked with cache_exists() before writing —
# re-running skips completed steps. Delete the relevant cache/output file
# to force a rerun of that step.
#
# RERUNNING AFTER A CORRECTION (workflow/R/engine/99_rerun_sites):
# rerun_engine_sites() handles two scenarios in one call — a raw
# coordinate that was wrong and has been fixed in your sites source
# (resnap_site_ids — fully re-delineates from the new lon/lat), and a
# snapped pour point that was wrong and has been manually corrected by
# editing output/<site_id>/pour_point_snapped.shp in QGIS directly
# (edited_snap_site_ids — reruns from that edited point, no re-snap).
# Either way it also reruns remove-upstream (always the full sites list —
# cheap, and a changed site's shape can change which OTHER sites it's
# nested within), reclip, metrics, and hydroweight for just the affected
# sites, MERGED into the existing combined outputs rather than overwritten
# wholesale. See an existing project's own runner header for a filled-in
# example.
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
source(here("workflow/R/stream/burn_streams.R")) # burn_streams_into_dem() — needed if streams_burn$source != "none"
source(here("workflow/R/stream/delineate_sites.R")) # reused building blocks for point-mode delineation
source(here("workflow/R/stream/hydroweight_attributes.R")) # calculate_hydroweight_attributes_stream() — used in Stage 7
source(here("workflow/R/remove_upstream.R"))
source(here("workflow/R/reclip_outputs.R"))
source(here("workflow/R/catchment_metrics.R"))
source(here("workflow/R/engine/00_resolve_config.R"))
source(here("workflow/R/engine/01_build_group_manifest.R"))
source(here("workflow/R/engine/02_prepare_terrain.R"))
source(here("workflow/R/engine/03_prepare_streams_burn.R"))
source(here("workflow/R/engine/04_delineate_site.R"))
source(here("code/reset_workflow.R")) # reset_site() — used by rerun_engine_sites()'s resnap path
source(here("workflow/R/engine/99_rerun_sites")) # rerun_engine_site_watershed(), rerun_engine_sites() — see "RERUNNING AFTER A CORRECTION" above
# TODO: if grouping$strategy = "hydrobasins" below, also source:
# source(here("workflow/R/stream/group_sites.R"))

# =============================================================================
# CONFIGURATION — fill in every TODO before running
# =============================================================================

PROJECT_ID  <- "TODO" # e.g. "MYPROJECT"
DELINEATION <- "TODO" # e.g. "stream_delineation" — only needed if this
# project has more than one delineation approach sharing PROJECT_ID
# (otherwise just fold this into output_dir/cache_dir directly, no
# DELINEATION subfolder needed — see an existing project whose runner
# uses only one delineation approach)

output_dir <- here("output", PROJECT_ID, DELINEATION)
cache_dir  <- here("cache", PROJECT_ID, DELINEATION)
fs::dir_create(output_dir, recurse = TRUE)
fs::dir_create(cache_dir, recurse = TRUE)

# -- Terrain: supply AT LEAST ONE of dem$path / flow_pointer$path /
# flow_direction$path. Highest-conditioning tier wins (flow_pointer >
# flow_direction > dem) — see workflow/R/engine/00_resolve_config.R.
# `dem` may ALSO be supplied alongside flow_pointer/flow_direction purely
# to provide an elevation surface for per-site clipping/hydroweight, even
# when it isn't what's used to derive flow direction.

# TODO: raw DEM needing depression breaching + D8 pointer generation, OR
# NULL if you're supplying flow_pointer/flow_direction directly instead.
dem <- list(path = NULL) # e.g. list(path = "~/path/to/mrdem-30-dtm.vrt")

# TODO: an EXTERNAL, already-conditioned flow direction raster needing
# recode to WhiteboxTools' encoding (1/2/4/8/16/32/64/128 clockwise from
# E), OR NULL. `recode` is a 2-column (from, to) matrix passed to
# terra::classify(); NULL if the source is already WBT-encoded.
flow_direction <- NULL # e.g. list(path = "~/path/to/flow_direction.tif", recode = my_recode_matrix)

# TODO: an already-WBT-encoded D8 pointer, OR NULL (skips breach + recode
# entirely if supplied).
flow_pointer <- NULL

# TODO: a pre-computed flow accumulation raster, OR NULL (derived from
# flow_pointer via wbt_d8_flow_accumulation() if absent).
flow_accum <- NULL

# TODO: leave NULL to match your terrain source's own native CRS (the
# default — no reprojection happens). Set to an "EPSG:####" string to
# force every terrain raster to be reprojected to it instead — e.g. if
# your DEM's native CRS is geographic (lon/lat) or otherwise unsuitable
# for watershed delineation. Either way, resolve_engine_config() warns if
# the resulting working CRS is geographic or not metre-based — WhiteboxTools'
# flow algorithms (and this workflow's own distance/area logic) assume a
# projected, metre-based CRS.
crs <- NULL

STREAM_THRESHOLD <- 1000 # flow accum threshold (cells) for stream extraction — only used when streams.tif needs generating
SNAP_DIST_M      <- 200  # point-mode pour point snap distance (m)
MIN_CELLS        <- 10   # minimum catchment size before flagging

# -- Burn-in streams: only applies when starting from a raw dem (no
# flow_pointer/flow_direction supplied) — burning has to happen before
# flow direction is derived, not after.
#   "nhn_auto"  — download NHN automatically for the AOI (requires
#                 nhn_index_path + nhn_raw_dir below, and the RCurl package)
#   "supplied"  — read streams_burn$path directly (any vector format
#                 sf::st_read() handles)
#   "none"      — skip burning (e.g. terrain is already conditioned)
streams_burn <- list(source = "none") # TODO
nhn_index_path <- NULL # TODO if source = "nhn_auto" — path to NHN_INDEX_WORKUNIT_LIMIT_2 shapefile
nhn_raw_dir    <- NULL # TODO if source = "nhn_auto" — shared dir for downloaded/cached NHN GDBs

# -- Pour point geometry: leave lake_polygons = NULL for point pour points
# (Jenson snap). Supply a SpatVector/sf with a `matched_lake` column (see
# workflow/R/lake/match_lake_polygons.R's output shape) and give sites
# a `lake_name` column to switch to lake-polygon pour points instead.
lake_polygons <- NULL # TODO
lake_buffer_m <- 30

# -- Grouping: "whole_domain" (default) treats all sites as one group —
# right if your terrain source already covers the whole site set as one
# raster (e.g. OIH). Use "hydrobasins" only if your terrain source is a
# large national mosaic/VRT worth cropping per region (e.g. MRDEM) — in
# that case your sites table must already carry a group_id column, and you
# need workflow/R/stream/group_sites.R sourced (see above).
grouping <- list(strategy = "whole_domain") # TODO if needed: "hydrobasins", hydrobasins_dir = ..., hybas_level = 6, default_buffer_m = 1000

# =============================================================================
# SITE DEFINITIONS
# =============================================================================
# TODO: build a sites tibble with at least: site_id (alphanumeric/
# underscore/hyphen only — used as the output folder name), site_name,
# lon, lat (WGS84 decimal degrees). Add lake_name if lake_polygons is set
# above. Add group_id yourself only if grouping$strategy = "hydrobasins".
#
# Either read from a file:
#   sites <- readr::read_csv(here("data/<project>_sites.csv"))
# ...or build inline:
#   sites <- tibble::tribble(
#     ~site_id, ~site_name, ~lon,       ~lat,
#     "site01", "Site One", -80.123,    46.456,
#   )

sites <- NULL # TODO

# =============================================================================
# STAGE 1 — Resolve config, build group manifest
# =============================================================================

run_config <- list(
  project_id = PROJECT_ID,
  output_dir = output_dir,
  cache_dir  = cache_dir,
  sites      = sites,
  dem = dem,
  flow_direction = flow_direction,
  flow_pointer   = flow_pointer,
  flow_accum     = flow_accum,
  crs = crs,
  stream_threshold = STREAM_THRESHOLD,
  streams_burn = streams_burn,
  nhn_index_path = nhn_index_path,
  nhn_raw_dir    = nhn_raw_dir,
  lake_polygons = lake_polygons,
  lake_buffer_m = lake_buffer_m,
  grouping = grouping,
  loi_layers = NULL # set in Stage 7 below
)

config <- resolve_engine_config(run_config)

gm <- build_engine_group_manifest(config)
sites <- gm$sites
group_manifest <- gm$group_manifest
print(group_manifest)

# =============================================================================
# STAGE 2 — Prepare terrain
# =============================================================================

prepare_engine_terrain(config, group_manifest)

# =============================================================================
# STAGE 3 — Delineate catchments
# =============================================================================

results <- delineate_engine_catchments(
  config = config, sites = sites, group_manifest = group_manifest,
  snap_dist = SNAP_DIST_M, min_cells = MIN_CELLS
)
print(results)

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
# LOI (layer of interest) rasters — three ways to declare a layer,
# corresponding to the three fields loi_layers[[i]] can set (exactly one):
#   $path          — a single raster already covering the full site extent
#   $path_template — one pre-clipped raster per site, "{site_id}" placeholder
#   $path_lazy     — one shared source too large to reproject/cache whole,
#                    optionally "{group_id}"/"{site_id}" templated
#
# PREFER $path_lazy over $path even for a single flat file with no
# placeholders — it still works fine as a literal path (see
# resolve_site_loi_raster() in workflow/R/stream/06_hydroweight_
# attributes.R). $path prepares ONE shared SpatRaster up front and reuses
# that SAME object across every site x version job, each job's crop()/
# mask() operating on the FULL raster extent every time before narrowing
# down — confirmed in practice, on a real engine-based project, to
# eventually crash with a terra "[readStart] file does not exist" error
# partway through a many-job run (a temp-file GC race from heavy repeated
# crop/mask churn on a long-lived shared object). $path_lazy avoids this
# structurally: each
# site's LOI is pre-cropped down to just its own catchment ONCE, cached to
# a real per-site file, before any further processing touches it — much
# less data movement per job, no long-lived shared object. Reach for
# plain $path only for something already small enough that a full-extent
# crop/mask per job is genuinely cheap.
#
# Each LOI's actual prep (source data, mosaicking/reprojection/
# rasterization, etc.) belongs in its own workflow/<PROJECT>/prepare_*.R
# script, not here — see an existing project's own prepare_*.R scripts for
# the established pattern, and workflow/raster_attributes.R for reusable
# generic helpers (reclassify_categorical(), sens_slope_trend(),
# build_mosaic_vrt(), rasterize_competing_classes()).

# TODO: write workflow/<PROJECT>/prepare_*.R scripts and call them here,
# then build loi_layers, e.g.:
# loi_layers <- list(
#   list(path = "TODO", name = "TODO", type = "categorical", class_levels = TODO),
#   list(path = "TODO", name = "TODO", type = "continuous")
# )
#
# hw_results <- calculate_hydroweight_attributes_stream(
#   sites = sites, group_manifest = group_manifest, output_dir = output_dir,
#   cache_dir = cache_dir, loi_layers = loi_layers,
#   catchment_versions = c("unclipped", "clipped"),
#   raster_crs = config$working_crs # not hardcoded — matches your DEM's own CRS
# )
# # Drop "clipped" rows that duplicate "unclipped" — see Stage 4's comment.
# hw_results <- drop_redundant_clipped_rows(hw_results, upstream_results, site_col = "site")
# write_csv(hw_results, here("output", PROJECT_ID, paste0(PROJECT_ID, "_hydroweight.csv")))
