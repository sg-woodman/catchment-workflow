# 05_delineate_sites.R
# ---------------------------------------------------------------------------
# Delineates catchments for each site using the group-level Whitebox outputs.
#
# Steps (per site):
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

#' Delineate catchments for all sites
#'
#' Iterates over all sites, grouping them by group_id to load group-level
#' rasters once per group rather than once per site. Skips sites where
#' catchment.gpkg already exists in the output directory.
#'
#' @param sites          Validated sites tibble from validate_sites()
#' @param group_manifest sf tibble from build_group_manifest()
#' @param output_dir     Character. Root output directory
#' @param snap_dist      Numeric. Pour point snap distance in metres.
#'                       Default 200.
#' @param min_cells      Integer. Catchments smaller than this many cells are
#'                       flagged as suspiciously small. Default 10.
#'
#' @return A tibble summarising delineation results for all sites, with
#'   columns: site_id, status, catchment_cells, catchment_km2, flagged,
#'   flag_reason
delineate_catchments <- function(
    sites,
    group_manifest,
    output_dir,
    snap_dist = 200,
    min_cells = 10
) {

  # Disable S2 geometry — geographic coordinate operations in sf default to
  # S2 spherical geometry which can produce unexpected results when working
  # with the MRDEM's native geographic CRS (EPSG:4617). Planar geometry is
  # more appropriate for the spatial operations used here.
  sf::sf_use_s2(FALSE)

  cw_inform(glue::glue(
    "Delineating catchments for {nrow(sites)} sites across ",
    "{dplyr::n_distinct(sites$group_id)} groups..."
  ))

  # Process sites group by group so group rasters are loaded once per group
  results <- purrr::map(unique(sites$group_id), function(grp) {

    grp_sites    <- dplyr::filter(sites, group_id == grp)
    grp_manifest <- dplyr::filter(group_manifest, group_id == grp)
    grp_cache    <- grp_manifest$cache_dir

    cw_inform(glue::glue(
      "Group '{grp}': processing {nrow(grp_sites)} site(s)..."
    ))

    # Load group rasters once — reused for all sites in the group
    group_rasters <- load_group_rasters(grp_cache, grp)

    # Process each site in the group
    purrr::map(seq_len(nrow(grp_sites)), function(j) {

      site <- grp_sites[j, ]

      delineate_site(
        site          = site,
        grp_cache     = grp_cache,
        group_rasters = group_rasters,
        output_dir    = output_dir,
        snap_dist     = snap_dist,
        min_cells     = min_cells
      )
    }) |>
      dplyr::bind_rows()

  }) |>
    dplyr::bind_rows()

  # Report flagged sites
  flagged <- dplyr::filter(results, flagged)
  if (nrow(flagged) > 0) {
    cw_warn(glue::glue(
      "{nrow(flagged)} site(s) flagged — review pour point locations:\n",
      "{paste(flagged$site_id, ':', flagged$flag_reason, collapse = '\n')}"
    ))
  } else {
    cw_inform("All catchments passed size check.")
  }

  results
}

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
    dem          = fs::path(grp_cache, "dem.tif"),
    dem_breached = fs::path(grp_cache, "dem_breached.tif"),
    flow_pointer = fs::path(grp_cache, "flow_pointer.tif"),
    flow_accum   = fs::path(grp_cache, "flow_accum.tif"),
    streams      = fs::path(grp_cache, "streams.tif"),
    hillshade    = fs::path(grp_cache, "hillshade.tif")
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

# -- Site delineation --------------------------------------------------------

#' Delineate a single site catchment
#'
#' Internal function. Handles the full per-site workflow: pour point writing,
#' snapping, watershed delineation, polygon conversion, raster clipping, and
#' output writing.
#'
#' @param site          Single-row tibble from sites
#' @param grp_cache     Character. Group cache directory path
#' @param group_rasters Named list of SpatRaster from load_group_rasters()
#' @param output_dir    Character. Root output directory
#' @param snap_dist     Numeric. Snap distance in metres
#' @param min_cells     Integer. Minimum cells threshold for flagging
#'
#' @return Single-row tibble with delineation result for this site
delineate_site <- function(
    site,
    grp_cache,
    group_rasters,
    output_dir,
    snap_dist,
    min_cells
) {

  sid      <- site$site_id
  site_dir <- site_output_dir(output_dir, sid)

  catchment_path <- fs::path(site_dir, "catchment.gpkg")

  # Skip if catchment already exists
  if (cache_exists(catchment_path)) {
    cw_inform(glue::glue("Site '{sid}': catchment.gpkg found, skipping."))

    # Return a result row from existing outputs
    catchment <- sf::st_read(catchment_path, quiet = TRUE)
    cells     <- catchment$n_cells[1]
    res_m     <- terra::res(group_rasters$flow_accum)[1]
    km2       <- round((cells * res_m^2) / 1e6, 4)

    return(tibble::tibble(
      site_id         = sid,
      status          = "skipped (cached)",
      catchment_cells = cells,
      catchment_km2   = km2,
      flagged         = km2 < (min_cells * res_m^2 / 1e6),
      flag_reason     = if (km2 < (min_cells * res_m^2 / 1e6))
        "catchment smaller than min_cells threshold" else NA_character_
    ))
  }

  cw_inform(glue::glue("Site '{sid}': delineating catchment..."))

  tryCatch({

    # -- Step 1: Write pour point as .shp ----------------------------------
    pour_point_shp <- write_pour_point_shp(site, site_dir)

    # -- Step 2: Snap pour point to stream ---------------------------------
    snapped_shp <- snap_pour_point(
      pour_point_shp = pour_point_shp,
      streams        = group_rasters$streams,
      site_dir       = site_dir,
      site_id        = sid,
      snap_dist      = snap_dist
    )

    # -- Step 3: Delineate watershed ---------------------------------------
    watershed_tif <- delineate_watershed(
      snapped_shp  = snapped_shp,
      flow_pointer = group_rasters$flow_pointer,
      site_dir     = site_dir,
      site_id      = sid
    )

    # -- Step 4: Convert watershed raster to polygon -----------------------
    catchment_sf <- watershed_to_polygon(
      watershed_tif = watershed_tif,
      site_dir      = site_dir,
      site_id       = sid
    )

    # Count cells for size check
    watershed_rast <- terra::rast(watershed_tif)
    n_cells        <- sum(terra::values(watershed_rast) == 1, na.rm = TRUE)
    res_m          <- terra::res(group_rasters$flow_accum)[1]
    km2            <- round((n_cells * res_m^2) / 1e6, 4)

    # Flag suspiciously small catchments
    flagged     <- n_cells < min_cells
    flag_reason <- if (flagged) {
      glue::glue(
        "catchment is only {n_cells} cells ({km2} km2) — ",
        "pour point may have snapped to wrong location"
      )
    } else {
      NA_character_
    }

    if (flagged) {
      cw_warn(glue::glue("Site '{sid}': {flag_reason}"))
    } else {
      cw_inform(glue::glue(
        "Site '{sid}': catchment = {n_cells} cells ({km2} km2)."
      ))
    }

    # Add metadata to catchment polygon
    catchment_sf <- catchment_sf |>
      dplyr::mutate(
        site_id  = sid,
        n_cells  = n_cells,
        area_km2 = km2,
        flagged  = flagged
      )

    # -- Step 5: Write snapped pour point to output ------------------------
    # Read snapped pour point — already in EPSG:3979 since all processing
    # runs in the MRDEM native CRS
    snapped_sf <- sf::st_read(snapped_shp, quiet = TRUE) |>
      sf::st_transform(3979) |>
      dplyr::mutate(site_id = sid)

    sf::st_write(
      snapped_sf,
      fs::path(site_dir, "pour_point.gpkg"),
      delete_dsn = TRUE,
      quiet      = TRUE
    )

    # -- Step 6: Write catchment polygon -----------------------------------
    sf::st_write(
      catchment_sf,
      catchment_path,
      delete_dsn = TRUE,
      quiet      = TRUE
    )

    # -- Step 7: Clip and mask group rasters to catchment -----------------
    clip_rasters_to_catchment(
      catchment_sf  = catchment_sf,
      group_rasters = group_rasters,
      site_dir      = site_dir,
      site_id       = sid
    )

    # -- Step 8: Clip NHN flowlines to catchment ---------------------------
    clip_flowlines_to_catchment(
      catchment_sf     = catchment_sf,
      flowlines_path   = fs::path(grp_cache, "flowlines.gpkg"),
      site_dir         = site_dir,
      site_id          = sid
    )

    tibble::tibble(
      site_id         = sid,
      status          = "success",
      catchment_cells = n_cells,
      catchment_km2   = km2,
      flagged         = flagged,
      flag_reason     = flag_reason
    )

  }, error = function(e) {
    cw_warn(glue::glue("Site '{sid}': delineation failed — {e$message}"))

    tibble::tibble(
      site_id         = sid,
      status          = paste("failed:", e$message),
      catchment_cells = NA_integer_,
      catchment_km2   = NA_real_,
      flagged         = TRUE,
      flag_reason     = e$message
    )
  })
}

# -- Pour point helpers ------------------------------------------------------

#' Write a single site's pour point as a .shp file
#'
#' WhiteboxTools requires pour points as shapefiles. The point is written
#' in EPSG:3979 to match the group rasters.
#'
#' @param site    Single-row tibble from sites (with lon, lat in WGS84)
#' @param tmp_dir Character. Path to temporary directory for this site
#' @return Character. Path to the written .shp file
write_pour_point_shp <- function(site, tmp_dir) {

  pour_point_shp <- fs::path(tmp_dir, "pour_point.shp")

  # Write pour point in EPSG:3979 to match the MRDEM native CRS used for
  # all Whitebox processing steps.
  sf::st_as_sf(
    site,
    coords = c("lon", "lat"),
    crs    = 4326
  ) |>
    sf::st_transform(3979) |>
    sf::st_write(
      pour_point_shp,
      delete_dsn = TRUE,
      quiet      = TRUE
    )

  pour_point_shp
}

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

  snapped_shp  <- fs::path(site_dir, "pour_point_snapped.shp")
  streams_path <- fs::path(site_dir, "streams_tmp.tif")

  # Write streams raster to site dir for WhiteboxTools and inspection
  terra::writeRaster(
    streams,
    filename  = streams_path,
    overwrite = TRUE,
    gdal      = c("COMPRESS=LZW")
  )

  cw_inform(glue::glue("Site '{site_id}': snapping pour point ({snap_dist} m)..."))

  whitebox::wbt_jenson_snap_pour_points(
    pour_pts  = normalizePath(pour_point_shp, mustWork = TRUE),
    streams   = normalizePath(streams_path,   mustWork = TRUE),
    output    = normalizePath(snapped_shp,    mustWork = FALSE),
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

  watershed_tif     <- fs::path(site_dir, "watershed.tif")
  flow_pointer_path <- fs::path(site_dir, "flow_pointer_tmp.tif")

  # Write flow pointer to site dir for WhiteboxTools
  terra::writeRaster(
    flow_pointer,
    filename  = flow_pointer_path,
    overwrite = TRUE,
    gdal      = c("COMPRESS=LZW")
  )

  cw_inform(glue::glue("Site '{site_id}': delineating watershed..."))

  whitebox::wbt_watershed(
    d8_pntr  = normalizePath(flow_pointer_path, mustWork = TRUE),
    pour_pts = normalizePath(snapped_shp,       mustWork = TRUE),
    output   = normalizePath(watershed_tif,     mustWork = FALSE)
  )

  if (!fs::file_exists(watershed_tif)) {
    cw_abort(glue::glue(
      "Site '{site_id}': wbt_watershed() did not produce output. ",
      "Check that the snapped pour point falls within the flow pointer extent."
    ))
  }

  watershed_tif
}

#' Convert watershed raster to a catchment polygon
#'
#' Uses wbt_raster_to_vector_polygons() — the same WhiteboxTools engine used
#' for all other processing steps — to convert the binary watershed raster to
#' a polygon. This produces cleaner boundaries than terra::as.polygons() as it
#' operates natively on the Whitebox raster grid.
#'
#' The polygon is filtered to VALUE == 1 (catchment cells) and dissolved to
#' a single feature in EPSG:3979, matching the MRDEM native CRS used throughout
#' all Whitebox processing.
#'
#' @param watershed_tif Character. Path to watershed raster
#' @param site_dir      Character. Site output directory (for tmp .shp output)
#' @param site_id       Character. Site identifier (for log messages)
#' @return sf polygon object in EPSG:3979
watershed_to_polygon <- function(watershed_tif, site_dir, site_id) {

  # wbt_raster_to_vector_polygons() always writes a .shp (not .gpkg)
  catchment_shp <- fs::path(site_dir, "catchment_tmp.shp")

  whitebox::wbt_raster_to_vector_polygons(
    input  = normalizePath(watershed_tif,  mustWork = TRUE),
    output = normalizePath(catchment_shp,  mustWork = FALSE)
  )

  if (!fs::file_exists(catchment_shp)) {
    cw_abort(glue::glue(
      "Site '{site_id}': wbt_raster_to_vector_polygons() did not produce ",
      "output. The watershed raster may be all-NoData — check that the ",
      "snapped pour point falls within the flow pointer extent."
    ))
  }

  catchment_sf <- sf::st_read(catchment_shp, quiet = TRUE) |>
    # VALUE == 1 are catchment cells; 0 is background
    dplyr::filter(VALUE == 1) |>
    sf::st_union() |>
    sf::st_as_sf() |>
    dplyr::rename(geometry = x) |>
    # Ensure output is in EPSG:3979 — matches MRDEM native CRS
    sf::st_transform(3979)

  # Clean up temporary shapefile components
  fs::dir_ls(site_dir, glob = "catchment_tmp.*") |>
    fs::file_delete()

  if (nrow(catchment_sf) == 0 || sf::st_is_empty(catchment_sf$geometry[1])) {
    cw_abort(glue::glue(
      "Site '{site_id}': catchment polygon is empty after filtering VALUE == 1. ",
      "The pour point may be outside the flow pointer extent."
    ))
  }

  catchment_sf
}

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
      c("LINESTRING", "MULTILINESTRING"), ,
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
      filename  = out_path,
      overwrite = TRUE,
      gdal      = c("COMPRESS=LZW", "BIGTIFF=IF_SAFER")
    )

    cw_inform(glue::glue("Site '{site_id}': {name}.tif written."))
  })

  invisible(NULL)
}
