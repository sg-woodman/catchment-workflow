# =============================================================================
# aoi.R — HydroBasins area of interest
# =============================================================================


#' Select which HydroBasins regional file(s) to search for a given latitude
#'
#' Returns one or both region codes depending on proximity to the 60°N boundary
#' between the North America and Arctic regional files.
#'
#' @param lat Decimal degrees latitude (WGS84)
#' @return Character vector: "na", "ar", or c("na", "ar")
hybas_regions <- function(lat) {
  near <- abs(lat - HYBAS_BOUNDARY_LAT) <= HYBAS_BOUNDARY_ZONE
  if (near)             return(c("na", "ar"))
  if (lat >= HYBAS_BOUNDARY_LAT) return("ar")
  return("na")
}


#' Find the HydroBasins polygon containing a pour point
#'
#' Searches one or both regional files and returns the containing basin polygon,
#' or NULL if not found (allows the caller to try the next Pfafstetter level).
#'
#' @param lon,lat Pour point coordinates (decimal degrees, WGS84)
#' @param level   Pfafstetter level (integer)
#' @return sf polygon in EPSG:4326, or NULL
hybas_containing <- function(lon, lat, level) {
  regions <- hybas_regions(lat)
  pt      <- st_sfc(st_point(c(lon, lat)), crs = 4326)

  for (region in regions) {
    path <- file.path(
      HYDROBASINS_DIRS[[region]],
      paste0("hybas_", region, "_lev", sprintf("%02d", level), "_v1c.shp")
    )
    if (!file.exists(path)) next

    basin <- st_read(path, quiet = TRUE) |>
      st_transform(4326) |>
      (\(b) b[st_contains(b, pt, sparse = FALSE)[, 1], ])()

    if (nrow(basin) > 0) return(basin)
  }
  NULL
}


#' Build a buffered AOI for a pour point, expanding until the catchment fits
#'
#' Iterates through Pfafstetter levels (coarse to fine) until a basin is found.
#' The boundary check in delineate.R determines whether this AOI was large
#' enough — if not, the pipeline re-runs at the next level automatically via
#' the iterative targets branching in _targets.R.
#'
#' Returns an sf data frame with columns:
#'   - geometry   : buffered basin polygon (EPSG:3979)
#'   - hybas_id   : HydroBasins ID (character)
#'   - hybas_level: Pfafstetter level used
#'   - wscssda_key: sorted, collapsed WSCSSDA codes intersecting the AOI
#'                  (used as the NHN and DEM cache key)
#'
#' @param lon,lat          Pour point coordinates (decimal degrees, WGS84)
#' @param nhn_index_path   Path to the NHN work unit index shapefile
#' @param level            Pfafstetter level to use (default: first of HYBAS_LEVELS)
#' @return sf data frame (one row)
build_aoi <- function(lon, lat, nhn_index_path, level = HYBAS_LEVELS[1]) {
  basin <- hybas_containing(lon, lat, level)

  if (is.null(basin)) {
    stop(
      "Pour point (lon=", lon, ", lat=", lat, ") not found in any HydroBasins ",
      "file at level ", level, ".\n",
      "Check coordinates and that regional .shp files are present in ",
      "HYDROBASINS_DIRS."
    )
  }

  aoi_geom <- basin |>
    st_transform(CRS_PROJECTED) |>
    st_buffer(HYBAS_BUFFER_M) |>
    st_union()

  # Identify which NHN work units intersect this AOI — their sorted codes
  # become the cache key for all derived products (DEM, flow, NHN layers).
  # Keying on WSCSSDA rather than hybas_id ensures that catchments spanning
  # multiple HydroBasins polygons still map to a single continuous raster.
  nhn_index  <- st_read(nhn_index_path, quiet = TRUE) |>
    st_transform(CRS_PROJECTED)
  intersects <- st_intersects(nhn_index, aoi_geom, sparse = FALSE)[, 1]
  wscssda    <- sort(unique(tolower(trimws(nhn_index$WSCSSDA[intersects]))))

  if (length(wscssda) == 0) {
    stop("No NHN work units intersect the AOI for (lon=", lon, ", lat=", lat, ").")
  }

  st_sf(
    hybas_id    = as.character(basin$HYBAS_ID[1]),
    hybas_level = level,
    wscssda_key = paste(wscssda, collapse = "_"),
    geometry    = aoi_geom,
    crs         = CRS_PROJECTED
  )
}


#' Check whether a catchment polygon touches the AOI boundary
#'
#' A catchment that reaches the edge of the DEM is artificially truncated —
#' the AOI was too small to contain the full upstream area. This check insets
#' the AOI by BOUNDARY_MARGIN_M and tests whether the catchment is fully
#' interior to that inner polygon.
#'
#' @param catchment sf polygon of the delineated catchment
#' @param aoi       sf polygon of the AOI used for delineation
#' @return Logical — TRUE if the catchment touches the boundary (AOI too small)
catchment_touches_boundary <- function(catchment, aoi) {
  aoi_inner <- st_buffer(st_union(aoi), -BOUNDARY_MARGIN_M)
  !st_within(
    st_union(catchment),
    aoi_inner,
    sparse = FALSE
  )[1, 1]
}
