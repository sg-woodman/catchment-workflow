# delineate_sites.R
# ---------------------------------------------------------------------------
# Reusable building blocks for delineating one site's catchment from
# group-level Whitebox outputs. No standalone entry point lives in this
# file — engine/04_delineate_site.R's delineate_engine_point_site() is the
# active orchestrator, composing these same helpers in the same order (the
# file's own former orchestrators, delineate_catchments()/delineate_site(),
# were removed once their only caller — the completed, retired
# verify_stream_migration.R — was retired; see git history).
#
# Steps (per site), performed by that orchestrator using the helpers below:
#   1. Write pour point as a single-feature .shp file
#   2. Snap pour point to nearest stream using wbt_jenson_snap_pour_points()
#   3. Delineate watershed using wbt_watershed()
#   4. Convert watershed raster to polygon
#   5. Clip all group rasters to catchment extent and mask to catchment
#   6. Write all outputs to site output directory
#   7. Flag suspiciously small catchments
#
# All Whitebox calls use absolute paths and .shp inputs.
#
# Inputs:
#   sites          : Validated sites tibble from validate_sites()
#   group_manifest : sf tibble from build_group_manifest()
#   output_dir     : Root output directory
#   snap_dist      : Global pour point snap distance in metres. Default 200.
#   min_cells      : Minimum catchment size in cells before flagging as
#                    suspicious. Default 10.
#
# Outputs (per site, written to output/<site_id>/):
#   pour_point.shp        : Original pour point (from sites CSV/tibble)
#   pour_point_snapped.shp: Snapped pour point — edit this to correct
#                           misaligned pour points and rerun delineation
#   streams_tmp.tif       : Streams raster cropped to site area, used
#                           for snapping — retained for manual inspection
#   watershed.tif         : Binary watershed raster (1 = in catchment)
#   catchment.gpkg        : Catchment polygon in EPSG:3979
#   pour_point.gpkg       : Snapped pour point in EPSG:3979
#   dem.tif               : DEM clipped and masked to catchment
#   dem_breached.tif      : Breached DEM clipped and masked to catchment
#   flow_pointer.tif      : Flow pointer clipped and masked to catchment
#   flow_accum.tif        : Flow accumulation clipped and masked to catchment
#   streams.tif           : Streams raster clipped and masked to catchment
#   hillshade.tif         : Hillshade clipped and masked to catchment
#   streams.gpkg          : NHN flowlines clipped to catchment polygon
#
# Dependencies: sf, terra, whitebox, dplyr, purrr, fs, cli (via utils.R)
# ---------------------------------------------------------------------------

# -- Group raster loading ----------------------------------------------------

#' Load all group-level rasters into a named list
#'
#' Loads the four rasters needed for delineation and clipping. Called once
#' per group so rasters are not reloaded for every site.
#'
#' @param grp_cache Character. Path to group cache directory
#' @param grp       Character. Group identifier (for log messages)
#' @return Named list of SpatRaster objects: dem, dem_breached, flow_pointer,
#'   flow_accum
load_group_rasters <- function(grp_cache, grp) {
  raster_files <- list(
    dem = fs::path(grp_cache, "dem.tif"),
    dem_breached = fs::path(grp_cache, "dem_breached.tif"),
    flow_pointer = fs::path(grp_cache, "flow_pointer.tif"),
    flow_accum = fs::path(grp_cache, "flow_accum.tif"),
    streams = fs::path(grp_cache, "streams.tif"),
    hillshade = fs::path(grp_cache, "hillshade.tif")
  )

  # Verify all required rasters exist before loading
  missing <- purrr::keep(raster_files, function(p) !cache_exists(p))
  if (length(missing) > 0) {
    cw_abort(glue::glue(
      "Group '{grp}': required rasters missing from cache: ",
      "{paste(names(missing), collapse = ', ')}. ",
      "Run prepare_dem() and run_whitebox() first."
    ))
  }

  purrr::map(raster_files, terra::rast)
}

# -- Pour point helpers ------------------------------------------------------
#
# write_pour_point_shp() (EPSG:3979-hardcoded) used to live here, alongside
# the engine's CRS-dynamic write_pour_point_shp_dynamic() (workflow/R/
# engine/04_delineate_site.R) — kept deliberately separate at the time so
# the engine wouldn't need to touch this "reused unmodified" file. Retired
# 2026-09 once confirmed to have zero live callers (every project migrated
# onto the engine, which only ever called its own _dynamic version) — a
# duplicate that could silently diverge from the version every pipeline
# actually runs is a worse risk than the "unmodified" convention it was
# preserving. See git history if you need the removed function; use
# workflow/R/engine/04_delineate_site.R's write_pour_point_shp() (renamed
# from write_pour_point_shp_dynamic() in the same cleanup) instead.

#' Snap pour point to nearest stream using wbt_jenson_snap_pour_points
#'
#' Uses the streams raster (from wbt_extract_streams) to snap the pour point
#' to the nearest stream cell within snap_dist metres.
#'
#' @param pour_point_shp Character. Path to pour point .shp
#' @param streams        SpatRaster. Binary streams raster from group cache
#' @param site_dir       Character. Site output directory path
#' @param site_id        Character. Site identifier (for log messages)
#' @param snap_dist      Numeric. Snap distance in metres
#' @return Character. Path to snapped pour point .shp
snap_pour_point <- function(
  pour_point_shp,
  streams,
  site_dir,
  site_id,
  snap_dist
) {
  snapped_shp <- fs::path(site_dir, "pour_point_snapped.shp")
  streams_path <- fs::path(site_dir, "streams_tmp.tif")

  # Write streams raster to site dir for WhiteboxTools and inspection
  terra::writeRaster(
    streams,
    filename = streams_path,
    overwrite = TRUE,
    gdal = c("COMPRESS=LZW")
  )

  cw_inform(glue::glue(
    "Site '{site_id}': snapping pour point ({snap_dist} m)..."
  ))

  whitebox::wbt_jenson_snap_pour_points(
    pour_pts = normalizePath(pour_point_shp, mustWork = TRUE),
    streams = normalizePath(streams_path, mustWork = TRUE),
    output = normalizePath(snapped_shp, mustWork = FALSE),
    snap_dist = snap_dist
  )

  if (!fs::file_exists(snapped_shp)) {
    cw_abort(glue::glue(
      "Site '{site_id}': wbt_jenson_snap_pour_points() did not produce ",
      "output. Check that the pour point falls within the group AOI."
    ))
  }

  snapped_shp
}

#' Delineate watershed from snapped pour point
#'
#' @param snapped_shp  Character. Path to snapped pour point .shp
#' @param flow_pointer SpatRaster. D8 flow pointer raster
#' @param site_dir     Character. Site output directory path
#' @param site_id      Character. Site identifier (for log messages)
#' @return Character. Path to watershed raster .tif
delineate_watershed <- function(
  snapped_shp,
  flow_pointer,
  site_dir,
  site_id
) {
  watershed_tif <- fs::path(site_dir, "watershed.tif")
  flow_pointer_path <- fs::path(site_dir, "flow_pointer_tmp.tif")

  # Write flow pointer to site dir for WhiteboxTools
  terra::writeRaster(
    flow_pointer,
    filename = flow_pointer_path,
    overwrite = TRUE,
    gdal = c("COMPRESS=LZW")
  )

  cw_inform(glue::glue("Site '{site_id}': delineating watershed..."))

  whitebox::wbt_watershed(
    d8_pntr = normalizePath(flow_pointer_path, mustWork = TRUE),
    pour_pts = normalizePath(snapped_shp, mustWork = TRUE),
    output = normalizePath(watershed_tif, mustWork = FALSE)
  )

  if (!fs::file_exists(watershed_tif)) {
    cw_abort(glue::glue(
      "Site '{site_id}': wbt_watershed() did not produce output. ",
      "Check that the snapped pour point falls within the flow pointer extent."
    ))
  }

  watershed_tif
}

# watershed_to_polygon() (EPSG:3979-hardcoded) used to live here, alongside
# the engine's CRS-dynamic watershed_to_polygon_dynamic() (workflow/R/
# engine/04_delineate_site.R) — same "kept separate to avoid touching this
# reused-unmodified file" reasoning as write_pour_point_shp() above, and
# retired for the same reason: zero live callers, confirmed directly, and
# a real bug already happened because of it — a fix applied here (adding
# sf::st_make_valid() to guard against a self-touching raster-to-polygon
# artifact) silently did nothing for any real project, because every
# project runs on the engine's own copy, not this one. See git history if
# you need the removed function; use workflow/R/engine/04_delineate_site.R's
# watershed_to_polygon() (renamed from watershed_to_polygon_dynamic() in
# the same cleanup, which is where the actual st_make_valid() fix lives)
# instead.

# -- Flowlines clipping ------------------------------------------------------

#' Clip NHN flowlines from group cache to the catchment polygon
#'
#' Reads the group-level flowlines.gpkg, clips to the catchment boundary,
#' and writes the result as streams.gpkg in the site output directory.
#' If no flowlines intersect the catchment (e.g. ungauged headwater sites),
#' an empty GeoPackage is written and a warning is issued.
#'
#' @param catchment_sf   sf polygon. Catchment boundary in EPSG:3979
#' @param flowlines_path Character. Path to group flowlines.gpkg
#' @param site_dir       Character. Site output directory path
#' @param site_id        Character. Site identifier (for log messages)
#' @return Invisibly NULL. Called for side effects.
clip_flowlines_to_catchment <- function(
  catchment_sf,
  flowlines_path,
  site_dir,
  site_id
) {
  out_path <- fs::path(site_dir, "streams.gpkg")

  if (cache_exists(out_path)) {
    cw_inform(glue::glue("Site '{site_id}': streams.gpkg found, skipping."))
    return(invisible(NULL))
  }

  # If no flowlines exist at group level (e.g. no NHN coverage), write empty
  if (!cache_exists(flowlines_path)) {
    cw_warn(glue::glue(
      "Site '{site_id}': group flowlines.gpkg not found at {flowlines_path}. ",
      "Writing empty streams.gpkg."
    ))
    sf::st_sf(geometry = sf::st_sfc(crs = 3979)) |>
      sf::st_write(out_path, delete_dsn = TRUE, quiet = TRUE)
    return(invisible(NULL))
  }

  flowlines <- sf::st_read(flowlines_path, quiet = TRUE)

  if (nrow(flowlines) == 0) {
    cw_warn(glue::glue(
      "Site '{site_id}': group flowlines.gpkg is empty. ",
      "Writing empty streams.gpkg."
    ))
    sf::st_write(flowlines, out_path, delete_dsn = TRUE, quiet = TRUE)
    return(invisible(NULL))
  }

  # Clip flowlines to catchment polygon
  clipped <- tryCatch(
    sf::st_intersection(flowlines, sf::st_union(catchment_sf)),
    error = function(e) {
      cw_warn(glue::glue(
        "Site '{site_id}': error clipping flowlines — {e$message}. ",
        "Writing empty streams.gpkg."
      ))
      flowlines[0, ]
    }
  )

  # Keep only line geometry types — intersection can return points where
  # stream lines touch the catchment boundary
  clipped <- clipped[
    sf::st_geometry_type(clipped) %in%
      c("LINESTRING", "MULTILINESTRING"),
    ,
    drop = FALSE
  ]

  sf::st_write(clipped, out_path, delete_dsn = TRUE, quiet = TRUE)

  cw_inform(glue::glue(
    "Site '{site_id}': streams.gpkg written ({nrow(clipped)} features)."
  ))

  invisible(NULL)
}

# -- Raster clipping ---------------------------------------------------------

#' Clip and mask all group rasters to the catchment extent
#'
#' Crops each group raster to the catchment bounding box then masks to the
#' catchment polygon. Written to the site output directory.
#'
#' @param catchment_sf  sf polygon. Catchment boundary in EPSG:3979
#' @param group_rasters Named list of SpatRaster from load_group_rasters()
#' @param site_dir      Character. Site output directory path
#' @param site_id       Character. Site identifier (for log messages)
#' @return Invisibly NULL. Called for side effects.
clip_rasters_to_catchment <- function(
  catchment_sf,
  group_rasters,
  site_dir,
  site_id
) {
  catchment_vect <- terra::vect(catchment_sf)

  purrr::iwalk(group_rasters, function(rast, name) {
    out_path <- fs::path(site_dir, paste0(name, ".tif"))

    if (cache_exists(out_path)) {
      cw_inform(glue::glue("Site '{site_id}': {name}.tif found, skipping."))
      return(invisible(NULL))
    }

    clipped <- rast |>
      terra::crop(catchment_vect, snap = "out") |>
      terra::mask(catchment_vect)

    terra::writeRaster(
      clipped,
      filename = out_path,
      overwrite = TRUE,
      gdal = c("COMPRESS=LZW", "BIGTIFF=IF_SAFER")
    )

    cw_inform(glue::glue("Site '{site_id}': {name}.tif written."))
  })

  invisible(NULL)
}
