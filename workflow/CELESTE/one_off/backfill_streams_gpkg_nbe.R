# backfill_streams_gpkg_nbe.R
# =============================================================================
# One-off backfill: NBE was delineated before the engine/03_prepare_streams_
# burn.R fix (this branch) that makes resolve_streams_burn() persist
# flowlines.gpkg to the group cache. Since resolve_streams_burn() is only
# invoked when dem_burned.tif doesn't exist yet (02_prepare_terrain.R's
# cache_exists() guard), the fix alone does NOT retroactively backfill
# already-delineated groups — this script does that for NBE specifically
# (the group catchment_delineation_guide.qmd's worked example depends on).
#
# Writes:
#   cache/CELESTE/NBE/flowlines.gpkg (written as cache/CELESTE_engine/NBE/ at the time; renamed since)    (re-resolved via nhn_auto,
#                                                same source/logic the fix
#                                                would have used originally)
#   cache/CELESTE/NBE/waterbodies.gpkg (ditto)  (one-off, for the guide's
#                                                display figure only — see
#                                                03_prepare_streams_burn.R's
#                                                comment on why this isn't
#                                                part of the generic fix)
#   output/CELESTE/<site>/streams.gpkg (ditto)  (per NBE site, 20 sites)
# =============================================================================

library(sf)
library(dplyr)
library(purrr)
library(glue)
library(fs)
library(here)

source(here("workflow/R/utils.R"))
source(here("workflow/R/stream/group_sites.R"))
source(here("workflow/R/stream/burn_streams.R")) # find_nhn_sheets(), read_nhn_from_gdb()
source(here("workflow/R/stream/delineate_sites.R")) # clip_flowlines_to_catchment()
source(here("workflow/R/engine/00_resolve_config.R"))
source(here("workflow/R/engine/01_build_group_manifest.R"))
source(here("workflow/R/engine/03_prepare_streams_burn.R")) # download_nhn()

PROJECT_ID <- "CELESTE" # was "CELESTE_engine" when this script ran (2026-08-31), before the rename
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
  project_id = PROJECT_ID, output_dir = output_dir, cache_dir = cache_dir, sites = sites,
  dem = list(path = mrdem_vrt), flow_direction = NULL, flow_pointer = NULL,
  crs = NULL, stream_threshold = 1000,
  streams_burn = list(source = "nhn_auto"), nhn_index_path = nhn_index, nhn_raw_dir = nhn_dir,
  lake_polygons = NULL, lake_buffer_m = 30,
  grouping = list(strategy = "hydrobasins", hydrobasins_dir = hydrobasins_dir, default_buffer_m = 1000),
  loi_layers = NULL
)
config <- resolve_engine_config(run_config)

gm <- build_engine_group_manifest(config)
sites_all <- gm$sites
group_manifest_all <- gm$group_manifest

nbe_manifest <- dplyr::filter(group_manifest_all, group_id == "NBE")
nbe_sites <- dplyr::filter(sites_all, group_id == "NBE")
grp_cache <- nbe_manifest$cache_dir[1]
aoi_sf <- sf::st_as_sf(nbe_manifest$aoi)

cw_inform(glue::glue("NBE: {nrow(nbe_sites)} site(s), grp_cache = {grp_cache}"))

# -- flowlines.gpkg (same source resolve_streams_burn() would have used) ----

flowlines_path <- fs::path(grp_cache, "flowlines.gpkg")
if (!cache_exists(flowlines_path)) {
  flowlines <- download_nhn(
    aoi = aoi_sf, nhn_index_path = config$nhn_index_path, nhn_raw_dir = config$nhn_raw_dir
  )
  if (is.null(flowlines) || nrow(flowlines) == 0) {
    cw_abort("NBE: no flowlines resolved — cannot backfill.")
  }
  sf::st_write(flowlines, flowlines_path, delete_dsn = TRUE, quiet = TRUE)
  cw_inform(glue::glue("NBE: flowlines.gpkg written ({nrow(flowlines)} feature(s))."))
} else {
  cw_inform("NBE: flowlines.gpkg already present, skipping.")
}

# -- waterbodies.gpkg (one-off, guide display figure only) ------------------

waterbodies_path <- fs::path(grp_cache, "waterbodies.gpkg")
if (!cache_exists(waterbodies_path)) {
  nhn_idx <- sf::st_read(config$nhn_index_path, quiet = TRUE)
  aoi_idx_crs <- sf::st_transform(aoi_sf, sf::st_crs(nhn_idx)) |> sf::st_buffer(0.05)
  sheets <- find_nhn_sheets(aoi_idx_crs, nhn_idx, config$nhn_raw_dir, "NBE")
  waterbodies <- purrr::map(sheets, function(gdb) {
    read_nhn_from_gdb(gdb, "NHN_HD_WATERBODY_2", "NBE")
  }) |>
    purrr::compact() |>
    dplyr::bind_rows() |>
    sf::st_transform(sf::st_crs(aoi_sf)) |>
    sf::st_make_valid()
  waterbodies <- waterbodies[!sf::st_is_empty(waterbodies), ]
  waterbodies <- sf::st_filter(waterbodies, aoi_sf)
  sf::st_write(waterbodies, waterbodies_path, delete_dsn = TRUE, quiet = TRUE)
  cw_inform(glue::glue("NBE: waterbodies.gpkg written ({nrow(waterbodies)} feature(s))."))
} else {
  cw_inform("NBE: waterbodies.gpkg already present, skipping.")
}

# -- per-site streams.gpkg backfill ------------------------------------------

results <- purrr::map_chr(nbe_sites$site_id, function(sid) {
  site_dir <- site_output_dir(output_dir, sid)
  catchment_path <- fs::path(site_dir, "catchment.gpkg")
  if (!cache_exists(catchment_path)) {
    cw_warn(glue::glue("Site '{sid}': no catchment.gpkg — skipping."))
    return("skipped (no catchment)")
  }
  catchment_sf <- sf::st_read(catchment_path, quiet = TRUE)
  clip_flowlines_to_catchment(
    catchment_sf = catchment_sf, flowlines_path = flowlines_path,
    site_dir = site_dir, site_id = sid
  )
  "done"
})

print(table(results))
cw_inform("NBE streams.gpkg backfill complete.")
