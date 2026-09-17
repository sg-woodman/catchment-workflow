# fix_ncmn_sliver_20260917.R
# =============================================================================
# One-off: patches NCMN's catchment to absorb a small sliver (~1.8 ha, 4
# fragments) that SUD103's clipped-catchment erasure (SUD103 minus NCMN u
# SUD101 u SUD200) leaves behind as a disconnected fragment, tripping
# remove_upstream.R's fragmentation integrity check and falling back to
# SUD103's full unclipped catchment (94 km2) instead of the correct ~49 km2
# clipped result.
#
# Root cause: NCMN was NOT included in the 2026-09-01 lake-bisection fix
# (workflow/CAM/fix_lake_bisection.R) -- it still uses its original,
# uncorrected D8 pointer. SUD103 (and SUD102, Tilton) WERE redelineated
# against a lake-flattened, corrected D8 pointer for that fix. NCMN and
# SUD103 therefore trace their shared boundary from two different flow
# fields, leaving a small strip that SUD103's (corrected) trace claims but
# NCMN's (uncorrected) trace doesn't -- exactly the "neighbor used a
# different D8 pointer" scenario remove_upstream.R's own integrity-check
# comment warns about.
#
# Sam's call (2026-09-17, after reviewing in QGIS): this sliver genuinely
# belongs to NCMN's catchment and NCMN's boundary is the one that's wrong,
# not SUD103's. Fix: union the sliver directly into NCMN's catchment.gpkg
# (a geometric patch, not a full redelineation with the corrected pointer --
# deliberately scoped to just this boundary, not hardcoded into the shared
# engine). This also resolves SUD103's fragmentation failure as a natural
# side effect, since erasing the now-larger NCMN cleanly covers the area
# that used to poke out as a disconnected sliver.
#
# Cascades exactly like rerun_engine_sites() would (remove-upstream on the
# full list, reclip/metrics/hydroweight merged for NCMN + any site whose
# erased_site_ids changes as a result -- expected: SUD103) but doesn't use
# rerun_engine_sites() itself since there's no re-snap/re-delineation step
# here, just a direct geometry patch to an already-delineated site.
# =============================================================================

suppressPackageStartupMessages({
  library(sf); library(terra); library(whitebox); library(dplyr); library(tidyr)
  library(purrr); library(readr); library(tibble); library(fs); library(cli); library(glue); library(here)
})

source(here("workflow/R/utils.R"))
source(here("workflow/R/stream/burn_streams.R"))
source(here("workflow/R/stream/delineate_sites.R"))
source(here("workflow/R/stream/hydroweight_attributes.R"))
source(here("workflow/R/remove_upstream.R"))
source(here("workflow/R/reclip_outputs.R"))
source(here("workflow/R/catchment_metrics.R"))
source(here("workflow/R/engine/00_resolve_config.R"))
source(here("workflow/R/engine/01_build_group_manifest.R"))
source(here("workflow/R/engine/02_prepare_terrain.R"))
source(here("workflow/R/engine/03_prepare_streams_burn.R"))
source(here("workflow/R/engine/04_delineate_site.R"))
source(here("code/reset_workflow.R"))
source(here("workflow/R/engine/99_rerun_sites")) # combine_all_catchments(), drop_redundant_clipped_rows(), merge_rows_into_csv()
source(here("workflow/gee_utils.R"))
source(here("workflow/raster_attributes.R"))
source(here("workflow/CAM/prepare_ndvi.R"))
source(here("workflow/CAM/prepare_ndvi_trend.R"))
source(here("workflow/CAM/prepare_ndvi_masked.R"))
source(here("workflow/CAM/prepare_harvest_regen.R"))

PROJECT_ID  <- "CAM"
DELINEATION <- "stream_delineation"
output_dir <- here("output", PROJECT_ID, DELINEATION)
cache_dir  <- here("cache", PROJECT_ID, DELINEATION)
OIH_DEM_PATH       <- "/Users/sam/Downloads/IntegratedHydrologyNE/EnforcedDEM.tif"
OIH_FLOW_DIR_PATH  <- "/Users/sam/Downloads/IntegratedHydrologyNE/EnhancedFlowDirection.tif"
oih_recode_matrix <- matrix(c(128,1,1,2,2,4,4,8,8,16,16,32,32,64,64,128), ncol = 2, byrow = TRUE)

sites_raw <- readr::read_csv(here("data/cam_stream_sites_raw.csv"), show_col_types = FALSE)
sites <- sites_raw |>
  dplyr::mutate(
    site_id = stream_id |> gsub("\\s+", "_", x = _) |> gsub("[^A-Za-z0-9_-]", "", x = _),
    site_name = stream_id
  ) |>
  dplyr::select(site_id, site_name, lon, lat)

run_config <- list(
  project_id = "CAM_streams", output_dir = output_dir, cache_dir = cache_dir, sites = sites,
  dem = list(path = OIH_DEM_PATH),
  flow_direction = list(path = OIH_FLOW_DIR_PATH, recode = oih_recode_matrix),
  flow_pointer = NULL, flow_accum = NULL, crs = NULL, stream_threshold = 100,
  streams_burn = list(source = "none"),
  lake_polygons = NULL, lake_buffer_m = 30,
  grouping = list(strategy = "whole_domain"),
  loi_layers = NULL
)
config <- resolve_engine_config(run_config)
gm <- build_engine_group_manifest(config)
sites <- gm$sites
group_manifest <- gm$group_manifest

# =============================================================================
# Step 1 -- Recompute the sliver (same math remove_upstream.R itself uses),
# and split it into the genuine NCMN-adjacent fragment vs. 3 single-pixel
# (900 m2 = one 30m DEM cell) artifacts sitting against SUD102's own
# boundary -- confirmed via distance checks (see conversation): only one
# fragment actually touches NCMN (distance 0m); the other three are 1.4-13
# km from NCMN but distance 0m from SUD102, and are exactly one grid cell
# each -- D8 boundary noise near SUD102's outlet, not real drainage area
# belonging to any neighbor. Per Sam's call (2026-09-17): absorb the real
# fragment into NCMN, erase the 3 pixel artifacts directly from SUD103
# (same treatment as remove_upstream.R's own 1 m^2 boundary-touch filter,
# just applied post-hoc here since these were never a "nested site").
# =============================================================================

sf::sf_use_s2(FALSE)

sud103 <- sf::st_read(fs::path(site_output_dir(output_dir, "SUD103"), "catchment.gpkg"), quiet = TRUE) |> sf::st_make_valid()
ncmn   <- sf::st_read(fs::path(site_output_dir(output_dir, "NCMN"),   "catchment.gpkg"), quiet = TRUE) |> sf::st_make_valid()
sud101 <- sf::st_read(fs::path(site_output_dir(output_dir, "SUD101"), "catchment.gpkg"), quiet = TRUE) |> sf::st_make_valid()
sud200 <- sf::st_read(fs::path(site_output_dir(output_dir, "SUD200"), "catchment.gpkg"), quiet = TRUE) |> sf::st_make_valid()

erase_mask <- sf::st_union(sf::st_union(sf::st_union(ncmn, sud101), sud200))
diff_geom <- sf::st_make_valid(sf::st_difference(sud103, erase_mask))
parts <- sf::st_cast(sf::st_union(diff_geom), "POLYGON")
part_areas <- as.numeric(sf::st_area(parts))
main_idx <- which.max(part_areas)
sliver_parts <- parts[-main_idx]
sliver_area_ha <- sum(as.numeric(sf::st_area(sliver_parts))) / 1e4

cw_inform(glue::glue(
  "Sliver: {length(sliver_parts)} fragment(s), {round(sliver_area_ha, 3)} ha total."
))
stopifnot(length(sliver_parts) == 4) # known shape of this specific fix -- abort if the underlying data changed

sliver_sfc <- sf::st_sfc(sliver_parts, crs = sf::st_crs(ncmn))
dist_to_ncmn <- as.numeric(sf::st_distance(sliver_sfc, ncmn))
ncmn_sliver_idx <- which(dist_to_ncmn < 1) # touching NCMN
pixel_sliver_idx <- setdiff(seq_along(sliver_parts), ncmn_sliver_idx)

stopifnot(length(ncmn_sliver_idx) == 1, length(pixel_sliver_idx) == 3)
cw_inform(glue::glue(
  "-> 1 fragment ({round(as.numeric(sf::st_area(sliver_sfc[ncmn_sliver_idx]))/1e4, 3)} ha) touches NCMN -- absorbing into NCMN.\n",
  "-> {length(pixel_sliver_idx)} single-pixel fragment(s) ",
  "({round(sum(as.numeric(sf::st_area(sliver_sfc[pixel_sliver_idx])))/1e4, 3)} ha total) touch SUD102, not NCMN -- ",
  "erasing directly from SUD103 instead."
))

ncmn_sliver_sfc <- sliver_sfc[ncmn_sliver_idx]
pixel_slivers_sfc <- sliver_sfc[pixel_sliver_idx]

# =============================================================================
# Step 2 -- Patch NCMN's catchment.gpkg (genuine adjacent fragment only)
# =============================================================================

ncmn_patched_geom <- sf::st_union(sf::st_union(ncmn), sf::st_union(ncmn_sliver_sfc)) |> sf::st_make_valid()
ncmn_patched_parts <- sf::st_cast(ncmn_patched_geom, "POLYGON", warn = FALSE)
if (length(ncmn_patched_parts) != 1) {
  cw_abort(glue::glue(
    "Patched NCMN catchment is not a single connected polygon ({length(ncmn_patched_parts)} parts) -- ",
    "aborting, needs manual review before applying."
  ))
}

ncmn_patched <- ncmn |>
  dplyr::slice(1) |>
  sf::st_set_geometry(ncmn_patched_geom) |>
  dplyr::mutate(
    n_cells = NA_integer_, # stale after a geometric patch -- not used downstream (catchment_metrics.R recomputes area from geometry directly)
    area_km2 = round(as.numeric(sf::st_area(ncmn_patched_geom)) / 1e6, 4)
  )

cw_inform(glue::glue(
  "NCMN catchment: {round(as.numeric(sf::st_area(ncmn))/1e6, 4)} -> ",
  "{ncmn_patched$area_km2[1]} km2 after patch."
))

ncmn_dir <- site_output_dir(output_dir, "NCMN")
sf::st_write(ncmn_patched, fs::path(ncmn_dir, "catchment.gpkg"), delete_dsn = TRUE, quiet = TRUE)

# -- Force-overwrite NCMN's own per-site unclipped rasters to match the new
# catchment boundary (clip_rasters_to_catchment()'s cache_exists() check
# would otherwise skip these since they already exist from the original
# delineation) -- same force-overwrite pattern rerun_engine_site_watershed()
# uses.
ncmn_site <- dplyr::filter(sites, site_id == "NCMN")
ncmn_grp <- ncmn_site$group_id[1]
ncmn_grp_manifest <- dplyr::filter(group_manifest, group_id == ncmn_grp)
ncmn_grp_cache <- ncmn_grp_manifest$cache_dir[1]
ncmn_group_rasters <- load_group_rasters(ncmn_grp_cache, ncmn_grp)

raster_files <- fs::path(ncmn_dir, paste0(names(ncmn_group_rasters), ".tif"))
fs::file_delete(raster_files[fs::file_exists(raster_files)])
clip_rasters_to_catchment(
  catchment_sf = ncmn_patched, group_rasters = ncmn_group_rasters,
  site_dir = ncmn_dir, site_id = "NCMN"
)

# =============================================================================
# Step 3 -- Rebuild combined catchments, rerun downstream cascade
# =============================================================================

combine_all_catchments(sites, output_dir)

upstream_results <- remove_upstream_catchments(sites, output_dir)
cw_inform("\nSUD103 / NCMN / SUD102 rows in upstream_results (standard mechanism):")
print(dplyr::filter(upstream_results, site_id %in% c("SUD103", "NCMN", "SUD102")))

# -- Manually override SUD103's clipped catchment ---------------------------
# The standard remove_upstream_catchments() mechanism above can only erase
# OTHER SITES' own catchment polygons -- it has no way to also erase the 3
# single-pixel artifacts (they aren't a "site"), so it still hits the
# fragmentation integrity check for SUD103 and falls back to the full
# unclipped catchment, same as before the NCMN patch. Construct the
# correct clipped geometry directly: erase NCMN (now patched) + SUD101 +
# SUD200 + the 3 pixel artifacts from SUD103's unclipped catchment, write
# it in place of whatever the standard mechanism just wrote, and patch
# upstream_results' SUD103 row to match so the CSV merge below treats it
# as a genuine (not redundant) clipped row.
ncmn_patched_native <- sf::st_read(fs::path(site_output_dir(output_dir, "NCMN"), "catchment.gpkg"), quiet = TRUE) |> sf::st_make_valid()
final_erase_mask <- sf::st_union(sf::st_union(sf::st_union(ncmn_patched_native, sud101), sud200), sf::st_union(pixel_slivers_sfc))
sud103_clipped_geom <- sf::st_make_valid(sf::st_difference(sud103, final_erase_mask))
sud103_clipped_parts <- sf::st_cast(sf::st_union(sud103_clipped_geom), "POLYGON", warn = FALSE)
if (length(sud103_clipped_parts) != 1) {
  cw_abort(glue::glue(
    "Manually-corrected SUD103 clipped catchment is still not a single ",
    "connected polygon ({length(sud103_clipped_parts)} parts) -- aborting."
  ))
}

sud103_area_before <- round(as.numeric(sf::st_area(sud103)) / 1e6, 4)
sud103_area_after <- round(as.numeric(sf::st_area(sud103_clipped_geom)) / 1e6, 4)
cw_inform(glue::glue("SUD103 catchment_clipped.gpkg (manually corrected): {sud103_area_before} -> {sud103_area_after} km2"))

sud103_clipped_sf <- sf::st_sf(
  site_id = "SUD103",
  area_km2_before = sud103_area_before,
  area_km2_after = sud103_area_after,
  n_erased = 3L,
  geometry = sf::st_geometry(sud103_clipped_geom)
) |> sf::st_set_crs(sf::st_crs(sud103))

sf::st_write(
  sud103_clipped_sf,
  fs::path(site_output_dir(output_dir, "SUD103"), "catchment_clipped.gpkg"),
  delete_dsn = TRUE, quiet = TRUE
)

# reclip_outputs() (below) defaults to reading the COMBINED
# all_catchments_clipped.gpkg, not each site's own catchment_clipped.gpkg
# -- remove_upstream_catchments() above already wrote that combined file
# using SUD103's PRE-correction (fallback) geometry, so it must be rebuilt
# from the current (corrected) per-site files before reclipping, or
# SUD103's clipped rasters get cropped to the wrong footprint again. (Found
# the hard way: an earlier run of this fix skipped this step and left
# dem_clipped.tif etc. at the stale 92 km2 footprint despite
# catchment_clipped.gpkg itself being correct -- caught only because Sam
# asked "does this mean all output files have been updated?".)
combine_clipped_catchments(sites, output_dir)

upstream_results <- upstream_results |>
  dplyr::mutate(
    status = dplyr::if_else(site_id == "SUD103", "success (manually corrected -- pixel-artifact erase)", status),
    n_erased = dplyr::if_else(site_id == "SUD103", 3L, n_erased),
    erased_site_ids = dplyr::if_else(site_id == "SUD103", "NCMN, SUD101, SUD200", erased_site_ids),
    area_km2_before = dplyr::if_else(site_id == "SUD103", sud103_area_before, area_km2_before),
    area_km2_after = dplyr::if_else(site_id == "SUD103", sud103_area_after, area_km2_after)
  )

target_ids <- c("NCMN", "SUD103") # SUD103 added manually -- its erased_site_ids only reflects
# NCMN/SUD101/SUD200 (real sites), so the standard cascade-detection below
# (which scans erased_site_ids strings) would not have caught it on its own.
cascaded <- upstream_results |>
  dplyr::filter(!is.na(erased_site_ids)) |>
  dplyr::filter(purrr::map_lgl(
    strsplit(erased_site_ids, ",\\s*"),
    function(ids) any(trimws(ids) %in% target_ids)
  )) |>
  dplyr::pull(site_id)
cascaded <- setdiff(cascaded, target_ids)
affected <- union(target_ids, cascaded)
cw_inform(glue::glue("Affected sites (NCMN + SUD103 + cascaded neighbors): {paste(affected, collapse = ', ')}"))

affected_sites <- dplyr::filter(sites, site_id %in% affected)

reclip_outputs(sites = affected_sites, output_dir = output_dir, group_manifest = group_manifest)

new_metrics <- calculate_catchment_metrics(sites = affected_sites, output_dir = output_dir)
merge_rows_into_csv(
  new_metrics, fs::path(output_dir, "catchment_metrics.csv"),
  key_cols = c("site_id", "version"),
  filter_fn = function(df) drop_redundant_clipped_rows(df, upstream_results, "site_id")
)

CANLCC_PATH <- "/Users/sam/Documents/cfs/shared_data/raw/landcover/CAN_LLC_2020.tif"
canlcc_levels <- data.frame(
  stringsAsFactors = FALSE,
  ID = c(1L, 2L, 5L, 6L, 8L, 10L, 11L, 12L, 13L, 14L, 15L, 16L, 17L, 18L, 19L),
  Class = c("needleleaf_forest","taiga_needleleaf_forest","broadleaf_deciduous_forest",
    "mixed_forest","temperate_shrubland","temperate_grassland","shrubland_lichen_moss",
    "grassland_lichen_moss","barren_lichen_moss","wetland","cropland","barren","urban","water","snow_ice")
)
loi_layers <- list(
  list(path_lazy = CANLCC_PATH, name = "canlcc", type = "categorical", class_levels = canlcc_levels),
  list(path_template = fs::path(cache_dir, "hydroweight_loi", "ndvi", "{site_id}.tif"), name = "ndvi", type = "continuous"),
  list(path_template = fs::path(cache_dir, "hydroweight_loi", "ndvi_trend", "{site_id}.tif"), name = "ndvi_trend", type = "continuous",
       stats = c("distwtd_mean","distwtd_sd","mean","sd","median","min","max")),
  list(path_template = fs::path(cache_dir, "hydroweight_loi", "ndvi_masked", "{site_id}.tif"), name = "ndvi_masked", type = "continuous"),
  list(path_template = fs::path(cache_dir, "hydroweight_loi", "ndvi_trend_masked", "{site_id}.tif"), name = "ndvi_trend_masked", type = "continuous",
       stats = c("distwtd_mean","distwtd_sd","mean","sd","median","min","max")),
  list(path_lazy = fs::path(cache_dir, "hydroweight_loi", "harvest_regen", "{group_id}.tif"), name = "harvest_regen", type = "categorical", class_levels = harvest_regen_levels)
)

new_hw <- calculate_hydroweight_attributes_stream(
  sites = affected_sites, group_manifest = group_manifest,
  output_dir = output_dir, cache_dir = cache_dir, loi_layers = loi_layers,
  catchment_versions = c("unclipped", "clipped"), raster_crs = config$working_crs
)
merge_rows_into_csv(
  new_hw, fs::path(output_dir, paste0(config$project_id, "_hydroweight.csv")),
  key_cols = c("site", "version"),
  filter_fn = function(df) drop_redundant_clipped_rows(df, upstream_results, "site")
)

source(here("workflow/CAM/tidy_outputs.R"))
tidy_cam_outputs(output_dir = output_dir)

cw_inform("\n===== Done. NCMN patched, SUD103 sliver resolved. =====")
