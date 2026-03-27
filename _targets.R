# =============================================================================
# _targets.R — Catchment delineation pipeline
# =============================================================================
# Run with: targets::tar_make()
# Inspect:  targets::tar_visnetwork()
# Read:     targets::tar_read(catchments_all)
#
# ADDING A SITE: add one row to data/sites.csv — no code changes needed.
#
# PIPELINE OVERVIEW:
#
#   sites.csv
#       │
#       ├─► [per site]       build_aoi()        → aoi sf + wscssda_key
#       │
#       ├─► [per wscssda_key] make_cache_dir()  → cache_dir/
#       │       ├─► load_mrdem()                → dem_raw.tif
#       │       ├─► download_nhn_streams()      → streams sf
#       │       ├─► download_nhn_lakes()        → lakes sf
#       │       ├─► burn_streams()              → dem_burned.tif
#       │       ├─► breach_dem()                → dem_breached.tif
#       │       ├─► flow_direction()            → fdr.tif
#       │       └─► flow_accumulation()         → fac.tif
#       │
#       └─► [per site]       delineate_catchment() → catchment polygon
#
# CACHING STRATEGY:
#   DEM products (dem_raw through fac) are written to named files in
#   CACHE_DIR/{wscssda_key}/ and tracked by targets via tar_file(). This
#   means they persist across sessions, are inspectable in QGIS, and are
#   shared across projects that point to the same CACHE_DIR.
#
#   NHN sf objects (streams, lakes) and catchment polygons are stored in
#   the targets object store as RDS — fast to regenerate if needed.
#
#   WBT bridge files (transient inputs written before a WBT call) still
#   use tempfile() — they are not products worth persisting.

library(targets)
library(tarchetypes)

# Load all functions from R/
lapply(list.files(here::here("R"), full.names = TRUE), source)

whitebox::wbt_init()
sf::sf_use_s2(FALSE)

tar_option_set(
  packages = c("dplyr", "fs", "sf", "terra", "whitebox")
)


# =============================================================================
# Pipeline
# =============================================================================

list(

  # ---------------------------------------------------------------------------
  # 1. Sites — single source of truth
  # ---------------------------------------------------------------------------

  tar_target(
    sites,
    readr::read_csv(here::here("data/sites.csv"), show_col_types = FALSE)
  ),


  # ---------------------------------------------------------------------------
  # 2. AOI — one per site
  #
  # Returns a one-row sf data frame per site including wscssda_key, which
  # is used as the grouping key for all shared DEM and NHN products.
  # ---------------------------------------------------------------------------

  tar_target(
    aoi,
    build_aoi(
      lon           = sites$lon,
      lat           = sites$lat,
      nhn_index_path = NHN_INDEX_PATH
    ),
    pattern = map(sites)
  ),


  # ---------------------------------------------------------------------------
  # 3. Unique groups — one branch per wscssda_key
  #
  # Deduplicates sites that share NHN work units so the expensive DEM and NHN
  # steps run exactly once per spatial group, not once per site.
  # ---------------------------------------------------------------------------

  tar_target(
    aoi_groups,
    dplyr::bind_cols(sites, aoi) |>
      dplyr::distinct(wscssda_key, .keep_all = TRUE),
    pattern = NULL  # runs once after all aoi branches complete
  ),


  # ---------------------------------------------------------------------------
  # 4. Cache directories — one per wscssda_key
  # ---------------------------------------------------------------------------

  tar_target(
    cache_dir,
    make_cache_dir(aoi_groups$wscssda_key),
    pattern  = map(aoi_groups),
    format   = "file"
  ),


  # ---------------------------------------------------------------------------
  # 5. NHN layers — one set per wscssda_key
  #
  # GDBs are downloaded to NHN_RAW_DIR (permanent shared storage).
  # Streams and lakes are stored as RDS in the targets object store.
  # ---------------------------------------------------------------------------

  tar_target(
    nhn_streams,
    if (isTRUE(aoi_groups$burn_streams)) {
      download_nhn_streams(aoi_groups$aoi)
    } else {
      NULL
    },
    pattern = map(aoi_groups)
  ),

  tar_target(
    nhn_lakes,
    download_nhn_lakes(aoi_groups$aoi),
    pattern = map(aoi_groups)
  ),


  # ---------------------------------------------------------------------------
  # 6. DEM conditioning — all products written to cache_dir as named .tif files
  #
  # tar_file() tracks each output by path + content hash. If the file exists
  # and is unchanged, the target is skipped. This makes the pipeline safe to
  # resume across sessions and machines that share the same CACHE_DIR.
  #
  # The dependency chain is enforced by function signatures — each function
  # accepts only the output of its predecessor, making wrong-input errors
  # structurally impossible.
  # ---------------------------------------------------------------------------

  tar_file(
    dem_raw,
    load_mrdem(aoi_groups$aoi, cache_dir),
    pattern = map(aoi_groups, cache_dir)
  ),

  tar_file(
    dem_burned,
    burn_streams(dem_raw, nhn_streams, cache_dir),
    pattern = map(dem_raw, nhn_streams, cache_dir)
  ),

  tar_file(
    dem_breached,
    breach_dem(dem_burned, cache_dir),
    pattern = map(dem_burned, cache_dir)
  ),

  tar_file(
    fdr,
    flow_direction(dem_breached, cache_dir),
    pattern = map(dem_breached, cache_dir)
  ),

  tar_file(
    fac,
    flow_accumulation(fdr, cache_dir),
    pattern = map(fdr, cache_dir)
  ),


  # ---------------------------------------------------------------------------
  # 7. Catchment delineation — one polygon per site
  #
  # Each site is matched to its fdr/fac via wscssda_key. Stream extraction,
  # pour point snapping, and watershed delineation are site-specific and fast;
  # they use tempfiles internally (see delineate.R).
  # ---------------------------------------------------------------------------

  tar_target(
    catchment,
    {
      # Match this site to its wscssda_key group
      site_key <- aoi$wscssda_key
      idx      <- which(aoi_groups$wscssda_key == site_key)

      delineate_catchment(
        site      = sites,
        fdr_file  = fdr[idx],
        fac_file  = fac[idx],
        cache_dir = cache_dir[idx]   # for writing QC pour point to QGIS-inspectable location
      )
    },
    pattern = map(sites, aoi)
  ),


  # ---------------------------------------------------------------------------
  # 8. Boundary check — warn if any catchment may be truncated
  # ---------------------------------------------------------------------------

  tar_target(
    boundary_check,
    {
      truncated <- purrr::map2_lgl(catchment, aoi, catchment_touches_boundary)
      flagged   <- sites$site_name[truncated]
      if (length(flagged) > 0) {
        warning(
          length(flagged), " site(s) have catchments touching the AOI boundary ",
          "and may be truncated: ", paste(flagged, collapse = ", "), "\n",
          "Increase hybas_level for these sites in data/sites.csv."
        )
      }
      data.frame(site_name = sites$site_name, boundary_ok = !truncated)
    },
    pattern = NULL  # runs once after all catchment branches complete
  ),


  # ---------------------------------------------------------------------------
  # 9. Output — one .gpkg per site + merged file
  # ---------------------------------------------------------------------------

  tar_file(
    catchment_file,
    {
      out <- file.path(OUTPUT_DIR, paste0(sites$site_name, "_catchment.gpkg"))
      fs::dir_create(OUTPUT_DIR)
      sf::st_write(catchment, out, delete_dsn = TRUE, quiet = TRUE)
      out
    },
    pattern = map(sites, catchment)
  ),

  tar_file(
    lakes_file,
    {
      if (is.null(nhn_lakes)) return(character(0))
      # Lakes saved per wscssda_key group, not per site
      out <- file.path(OUTPUT_DIR, paste0(aoi_groups$wscssda_key, "_lakes.gpkg"))
      sf::st_write(nhn_lakes, out, delete_dsn = TRUE, quiet = TRUE)
      out
    },
    pattern = map(aoi_groups, nhn_lakes)
  ),

  tar_target(
    catchments_all,
    {
      merged <- dplyr::bind_rows(catchment)
      out    <- file.path(OUTPUT_DIR, "catchments_all.gpkg")
      sf::st_write(merged, out, delete_dsn = TRUE, quiet = TRUE)
      merged
    },
    pattern = NULL
  )

)
