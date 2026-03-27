# =============================================================================
# delineate.R — Catchment delineation from a snapped pour point
# =============================================================================
# Accepts file paths to fdr.tif and fac.tif (from dem.R).
#
# Pour point handling has two paths:
#   - Automatic: snap the nominal pour point to the nearest stream cell using
#     the Jenson method. The snapped location is written to cache_dir for
#     inspection in QGIS alongside fac.tif and the stream raster.
#   - Manual override: if lon_snapped / lat_snapped are provided in sites.csv,
#     auto-snapping is skipped entirely and the verified coordinates are used
#     directly. Set these after inspecting the auto-snapped result in QGIS.
#
# All intermediate files (stream raster, pour point shapefiles, watershed
# raster, catchment shapefile) use tempfile() — they are site-specific,
# fast to produce, and targets re-runs the whole function as a unit if any
# input changes, so there is no benefit to persisting them between sessions.


#' Snap a pour point to the nearest stream cell
#'
#' Extracts the stream network at the site's threshold, snaps the pour point
#' to the nearest stream pixel using the Jenson method, and writes the snapped
#' location to cache_dir as a .gpkg for QGIS inspection.
#'
#' @param site      One-row tibble: site_name, lon, lat, stream_threshold
#' @param fdr_file  Path to fdr.tif
#' @param fac_file  Path to fac.tif
#' @param cache_dir Path to the wscssda_key cache directory
#' @return sf point in the flow direction raster CRS
snap_pour_point <- function(site, fdr_file, fac_file, cache_dir) {
  fdr <- rast(fdr_file)
  crs_wkt <- paste0("EPSG:", crs(fdr, describe = TRUE)$code)

  pour_pt <- st_sfc(st_point(c(site$lon, site$lat)), crs = 4326) |>
    st_transform(crs_wkt) |>
    st_sf()

  # Extract stream network at this site's threshold
  streams_ras <- tempfile(fileext = ".tif")
  wbt_extract_streams(
    flow_accum = fac_file,
    output     = streams_ras,
    threshold  = site$stream_threshold
  )

  # Snap to nearest stream cell
  pour_shp    <- tempfile(fileext = ".shp")
  snapped_shp <- tempfile(fileext = ".shp")

  st_write(pour_pt, pour_shp, quiet = TRUE)
  wbt_jenson_snap_pour_points(
    pour_pts  = pour_shp,
    streams   = streams_ras,
    output    = snapped_shp,
    snap_dist = WBT_SNAP_DIST_M
  )

  if (!file.exists(snapped_shp))
    stop("Pour point snapping failed for site ", site$site_name)

  snapped <- st_read(snapped_shp, quiet = TRUE)
  if (nrow(snapped) == 0)
    stop("Snapped pour point is empty for site ", site$site_name,
         ". Try increasing WBT_SNAP_DIST_M (currently ", WBT_SNAP_DIST_M, " m).")

  # Write snapped point to cache_dir for QGIS inspection alongside fac.tif.
  # If the result looks wrong: identify the correct pixel in QGIS, then set
  # lon_snapped / lat_snapped for this site in data/sites.csv.
  st_write(
    snapped,
    file.path(cache_dir, paste0(site$site_name, "_pour_point_snapped.gpkg")),
    delete_dsn = TRUE,
    quiet      = TRUE
  )

  snapped
}


#' Delineate a catchment polygon for a single site
#'
#' @param site      One-row tibble: site_name, lon, lat, stream_threshold,
#'                  burn_streams, lon_snapped, lat_snapped
#' @param fdr_file  Path to fdr.tif (from flow_direction())
#' @param fac_file  Path to fac.tif (from flow_accumulation())
#' @param cache_dir Path to the wscssda_key cache directory (for QC outputs)
#' @return One-row sf data frame in CRS_OUTPUT with site metadata and area_km2
delineate_catchment <- function(site, fdr_file, fac_file, cache_dir) {
  message("\n", strrep("=", 60))
  message("  Site: ", site$site_name,
          " | lon=", round(site$lon, 5),
          " | lat=", round(site$lat, 5))
  message(strrep("=", 60))

  fdr     <- rast(fdr_file)
  crs_wkt <- paste0("EPSG:", crs(fdr, describe = TRUE)$code)

  # ---- Pour point: manual override or auto-snap ----
  has_manual <- !is.na(site$lon_snapped) && !is.na(site$lat_snapped)

  if (has_manual) {
    message("  Using manually verified pour point (",
            round(site$lon_snapped, 6), ", ", round(site$lat_snapped, 6), ")")
    snapped <- st_sfc(st_point(c(site$lon_snapped, site$lat_snapped)), crs = 4326) |>
      st_transform(crs_wkt) |>
      st_sf()
  } else {
    message("  Auto-snapping pour point (snap distance = ", WBT_SNAP_DIST_M, " m)...")
    snapped <- snap_pour_point(site, fdr_file, fac_file, cache_dir)
  }

  # Write snapped point to a tempfile for wbt_watershed input.
  # (The persistent QC copy is written inside snap_pour_point() above.)
  snapped_shp <- tempfile(fileext = ".shp")
  st_write(snapped, snapped_shp, quiet = TRUE)

  # ---- Delineate watershed ----
  message("  Delineating watershed...")
  catchment_ras <- tempfile(fileext = ".tif")
  wbt_watershed(d8_pntr = fdr_file, pour_pts = snapped_shp, output = catchment_ras)

  ras <- rast(catchment_ras)
  if (all(is.na(values(ras))))
    stop("Watershed produced all-NA raster for site ", site$site_name,
         ". Check CRS alignment between flow direction and pour point.")

  message("  Watershed cells: ", sum(!is.na(values(ras))))

  # ---- Convert raster to polygon ----
  catchment_shp <- tempfile(fileext = ".shp")
  wbt_raster_to_vector_polygons(input = catchment_ras, output = catchment_shp)

  catchment <- st_read(catchment_shp, quiet = TRUE) |>
    filter(VALUE == 1) |>
    st_union() |>
    st_sf(geometry = _) |>
    st_transform(CRS_OUTPUT) |>
    mutate(
      site_name      = site$site_name,
      lon            = site$lon,
      lat            = site$lat,
      pour_point_src = if (has_manual) "manual" else "auto",
      area_km2       = as.numeric(st_area(geometry)) / 1e6
    )

  message("  Catchment area: ", round(catchment$area_km2, 2), " km²")
  catchment
}
