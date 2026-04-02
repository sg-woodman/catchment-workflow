# =============================================================================
# run_delineation.R — Manual catchment delineation workflow
# =============================================================================
# Runs the delineation pipeline explicitly without targets. All functions
# are defined in R/ and called directly here in dependency order.
#
# WORKFLOW:
#   1. Add sites to data/sites.csv
#   2. Run sections 1–6 to build DEM products for all unique spatial groups
#   3. Run section 7 to delineate catchments
#   4. If a catchment is wrong, inspect cache_dir/{wscssda_key}/fac.tif and
#      {site_name}_pour_point_snapped.gpkg in QGIS, then set lon_snapped and
#      lat_snapped for that site in sites.csv and re-run section 7 only
#   5. Run section 8 to check for boundary-touching catchments
#   6. Run section 9 to save outputs
#
# CACHING:
#   DEM products (dem_raw.tif through fac.tif) are written to named files in
#   CACHE_DIR/{wscssda_key}/ by the functions in R/dem.R. Re-running this
#   script skips any step whose output file already exists — so adding a new
#   site that shares an existing spatial group only triggers delineation,
#   not the expensive WBT conditioning steps.
#
#   To force re-processing a spatial group (e.g. after changing burn depth),
#   delete the relevant files from CACHE_DIR/{wscssda_key}/ manually.
#
# EXTERNAL CATCHMENTS:
#   Sites with a path in the catchment_file column of sites.csv skip
#   delineation entirely. The provided .gpkg is loaded directly and flows
#   through to the output merge.

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

library(here)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(sf)
library(terra)
library(whitebox)
library(fs)

# Load all functions and constants from R/
lapply(list.files(here("R"), full.names = TRUE), source)

whitebox::wbt_init()
sf::sf_use_s2(FALSE)


# =============================================================================
# 1. Sites
# =============================================================================
# Single source of truth — edit data/sites.csv to add or modify sites.
# Required columns: site_name, lon, lat, stream_threshold, burn_streams
# Optional columns: lon_snapped, lat_snapped (manual pour point correction)
#                   catchment_file (path to existing catchment .gpkg)

sites <- read_csv(here("data/sites.csv"), show_col_types = FALSE)

# Separate sites that need delineation from those with existing catchments
sites_delineate <- sites |> filter(is.na(catchment_file))
sites_external <- sites |> filter(!is.na(catchment_file))

message("Sites to delineate: ", nrow(sites_delineate))
message("External catchments: ", nrow(sites_external))

test_site <- sites |>
  filter(site_name == "PHPP05")


# =============================================================================
# 2. AOI — one per site requiring delineation
# =============================================================================
# build_aoi() returns a one-row sf data frame per site with:
#   - geometry    : buffered HydroBasins polygon (EPSG:3979)
#   - hybas_id    : HydroBasins ID
#   - hybas_level : Pfafstetter level used
#   - wscssda_key : cache key — sorted NHN work unit codes for this AOI

message("\n--- Building AOIs ---")

aoi_list <- test_site |>
  rowwise() |>
  mutate(
    aoi = list(
      build_aoi(lon, lat, NHN_INDEX_PATH)
    )
  ) |>
  ungroup()


# =============================================================================
# 3. Unique spatial groups — one set of DEM products per wscssda_key
# =============================================================================
# Many sites share NHN work units and therefore share a single set of DEM
# products. Deduplicating here ensures the expensive WBT steps run once per
# spatial group, not once per site.

# unnest() flattens the aoi sf object into the data frame — the geometry
# becomes a bare sfc column accessible as st_sf(geometry = geometry) inside
# rowwise() mutate calls in sections 5 and 6.
aoi_groups <- aoi_list |>
  unnest(aoi) |>
  distinct(wscssda_key, .keep_all = TRUE)

message("\nUnique spatial groups (wscssda_key): ", nrow(aoi_groups))
message(paste0("  ", aoi_groups$wscssda_key, collapse = "\n"))


# =============================================================================
# 4. Cache directories — one per wscssda_key
# =============================================================================
# make_cache_dir() creates CACHE_DIR/{wscssda_key}/ if it does not exist
# and returns the path. All DEM products for this group are written here.

aoi_groups <- aoi_groups |>
  rowwise() |>
  mutate(cache_dir = make_cache_dir(wscssda_key)) |>
  ungroup()


# =============================================================================
# 5. NHN layers — streams and lakes per spatial group
# =============================================================================
# GDBs are downloaded to NHN_RAW_DIR permanently (shared across projects).
# download_nhn_streams() is only called when burn_streams = TRUE for at
# least one site in this group.
# download_nhn_lakes() reuses already-downloaded GDBs — no FTP call needed.

message("\n--- Downloading NHN layers ---")

aoi_groups <- aoi_groups |>
  rowwise() |>
  mutate(
    nhn_streams = list(
      if (isTRUE(burn_streams)) {
        download_nhn_streams(st_sf(geometry = geometry))
      } else {
        message(
          "  Skipping stream download for ",
          wscssda_key,
          " (burn_streams = FALSE)"
        )
        NULL
      }
    ),
    nhn_lakes = list(
      download_nhn_lakes(st_sf(geometry = geometry))
    )
  ) |>
  ungroup()


# =============================================================================
# 6. DEM conditioning — per spatial group
# =============================================================================
# Each function writes one named .tif to cache_dir and returns its path.
# If the output file already exists, the function skips processing and
# returns the existing path — re-running this section is always safe.
#
# Inspect any of these files in QGIS at any point:
#   cache_dir/dem_raw.tif      — raw MRDEM crop
#   cache_dir/dem_burned.tif   — stream-burned DEM
#   cache_dir/dem_breached.tif — depression-breached DEM
#   cache_dir/fdr.tif          — D8 flow direction
#   cache_dir/fac.tif          — D8 flow accumulation

message("\n--- DEM conditioning ---")

aoi_groups <- aoi_groups |>
  rowwise() |>
  mutate(
    dem_raw = load_mrdem(st_sf(geometry = geometry), cache_dir),
    # Save NHN layers here — dem_raw.tif now exists so streams can be clipped
    # to the exact raster extent, preventing WBT sentinel values at edges
    nhn_paths = list(save_nhn_layers(
      streams = nhn_streams[[1]],
      lakes = nhn_lakes[[1]],
      dem_file = dem_raw,
      cache_dir = cache_dir
    )),
    dem_burned = burn_streams(dem_raw, nhn_paths[["streams_shp"]], cache_dir),
    dem_breached = breach_dem(dem_burned, cache_dir),
    fdr = flow_direction(dem_breached, cache_dir),
    fac = flow_accumulation(fdr, cache_dir)
  ) |>
  ungroup()

message("\nDEM conditioning complete. Cache locations:")
walk2(
  aoi_groups$wscssda_key,
  aoi_groups$cache_dir,
  ~ message("  ", .x, " → ", .y)
)


# =============================================================================
# 7. Catchment delineation — per site
# =============================================================================
# Each site is matched to its spatial group via wscssda_key to get the
# correct fdr and fac file paths.
#
# Pour point handling:
#   - If lon_snapped / lat_snapped are set in sites.csv, those coordinates
#     are used directly (auto-snapping skipped)
#   - Otherwise the nominal pour point is snapped to the nearest stream cell
#
# The snapped pour point is written to cache_dir/{site_name}_pour_point_snapped.gpkg
# for inspection in QGIS alongside fac.tif.
#
# To correct a bad pour point:
#   1. Open fac.tif and {site_name}_pour_point_snapped.gpkg in QGIS
#   2. Identify the correct stream pixel
#   3. Set lon_snapped / lat_snapped for this site in data/sites.csv
#   4. Re-run this section only (sections 1–6 will skip cached products)

message("\n--- Delineating catchments ---")

catchments_delineated <- sites_delineate |>
  rowwise() |>
  mutate(
    catchment = list({
      # Capture site_name as a local variable before filtering aoi_groups.
      # Inside rowwise(), filter(site_name == site_name) is a tautology —
      # dplyr compares the column to itself rather than to the current row value.
      # Assigning to a local name avoids this ambiguity.
      this_site <- site_name
      grp <- aoi_groups |>
        filter(site_name == this_site) |>
        slice(1)

      delineate_catchment(
        site = pick(everything()),
        fdr_file = grp$fdr,
        fac_file = grp$fac,
        cache_dir = grp$cache_dir
      )
    })
  ) |>
  ungroup()

# Load external catchments
catchments_external <- sites_external |>
  rowwise() |>
  mutate(
    catchment = list(
      st_read(catchment_file, quiet = TRUE) |>
        mutate(
          site_name = site_name,
          pour_point_src = "external",
          area_km2 = as.numeric(st_area(geometry)) / 1e6
        )
    )
  ) |>
  ungroup()

# Combine
catchments <- bind_rows(catchments_delineated, catchments_external)


# =============================================================================
# 8. Boundary check
# =============================================================================
# Warns if any delineated catchment touches the AOI boundary, which indicates
# the upstream area may be truncated. External catchments are not checked.
#
# To fix a truncated catchment: increase hybas_level for that site in
# sites.csv (e.g. add a hybas_level column with value 5 or 4), then
# delete that site's cache directory and re-run from section 2.

message("\n--- Checking catchment boundaries ---")

catchments <- catchments |>
  rowwise() |>
  mutate(boundary_ok = {
    if (!is.na(catchment_file)) {
      TRUE # external catchments not checked
    } else {
      this_site <- site_name
      site_aoi <- aoi_groups |>
        filter(site_name == this_site) |>
        slice(1) |>
        st_sf(geometry = geometry)
      !catchment_touches_boundary(catchment[[1]], site_aoi)
    }
  }) |>
  ungroup()

flagged <- catchments |> filter(!boundary_ok) |> pull(site_name)

if (length(flagged) > 0) {
  warning(
    length(flagged),
    " site(s) have catchments that touch the AOI boundary ",
    "and may be truncated:\n",
    paste0("  ", flagged, collapse = "\n"),
    "\n",
    "Increase hybas_level for these sites in data/sites.csv and re-run."
  )
} else {
  message("All catchments are fully contained within their AOI.")
}


# =============================================================================
# 9. Save outputs
# =============================================================================

message("\n--- Saving outputs ---")

fs::dir_create(OUTPUT_DIR)

# One .gpkg per site
catchments |>
  rowwise() |>
  mutate(out_file = {
    path <- file.path(OUTPUT_DIR, paste0(site_name, "_catchment.gpkg"))
    st_write(catchment[[1]], path, delete_dsn = TRUE, quiet = TRUE)
    message("  Saved: ", path)
    path
  }) |>
  ungroup()

# Merged file — all sites in one .gpkg
all_catchments <- catchments |>
  pull(catchment) |>
  bind_rows()

st_write(
  all_catchments,
  file.path(OUTPUT_DIR, "catchments_all.gpkg"),
  delete_dsn = TRUE,
  quiet = TRUE
)
message("  Saved: ", file.path(OUTPUT_DIR, "catchments_all.gpkg"))

# Lakes — one per spatial group
aoi_groups |>
  filter(!map_lgl(nhn_lakes, is.null)) |>
  rowwise() |>
  mutate(out_file = {
    path <- file.path(OUTPUT_DIR, paste0(wscssda_key, "_lakes.gpkg"))
    st_write(nhn_lakes[[1]], path, delete_dsn = TRUE, quiet = TRUE)
    message("  Saved: ", path)
    path
  }) |>
  ungroup()

message("\nDone. All outputs saved to: ", OUTPUT_DIR)
