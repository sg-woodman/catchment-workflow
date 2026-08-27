# 04_delineate_site.R
# ---------------------------------------------------------------------------
# One delineation entry point, branching on whether config$lake_polygons is
# supplied (lake-polygon pour point) or NULL (point pour point / Jenson
# snap). Both modes reuse existing, unmodified building blocks wherever
# they're already CRS-generic:
#
#   Point mode reuses, verbatim, from workflow/R/stream/delineate_sites.R
#   (must be sourced first): load_group_rasters(), snap_pour_point(),
#   delineate_watershed(), clip_rasters_to_catchment(),
#   clip_flowlines_to_catchment(). Only write_pour_point_shp() and
#   watershed_to_polygon() hardcode EPSG:3979 there — this file defines
#   CRS-dynamic replacements (write_pour_point_shp_dynamic(),
#   watershed_to_polygon_dynamic()) rather than touching the originals.
#
#   Lake mode adapts workflow/R/lake/03_delineate_lakes.R's
#   delineate_single_lake() logic — that file is already CRS-dynamic
#   (terra::crs(d8_template), no hardcoded EPSG) so no CRS change is
#   needed there; the only real adaptation is the D8 pointer filename
#   (flow_pointer.tif here, matching the canonical names
#   02_prepare_terrain.R always produces, vs. d8_pntr.tif in the lake
#   pipeline) and reading from a group's cache_dir (from group_manifest)
#   instead of a single flat project cache_dir. Preserves the lake
#   pipeline's documented no-trim-before-wbt_watershed() rule verbatim
#   (pour-point raster must share exact extent/resolution with the D8
#   pointer; trim only AFTER delineation).
#
# Dependencies: sf, terra, whitebox, dplyr, purrr, fs, glue, cli (via utils.R)
# ---------------------------------------------------------------------------

#' Delineate catchments for all sites in an engine run
#'
#' @param config         Resolved config from resolve_engine_config()
#' @param sites          Validated sites tibble (group_id populated) from
#'   build_engine_group_manifest()
#' @param group_manifest sf tibble from build_engine_group_manifest()
#' @param snap_dist      Numeric. Point-mode pour point snap distance (m).
#' @param min_cells      Integer. Catchments smaller than this are flagged.
#' @return A tibble summarising delineation results, one row per site
delineate_engine_catchments <- function(
  config, sites, group_manifest, snap_dist = 200, min_cells = 10
) {
  if (!exists("load_group_rasters", mode = "function")) {
    cw_abort(paste(
      "delineate_engine_catchments() requires",
      "workflow/R/stream/delineate_sites.R to be sourced first."
    ))
  }

  is_lake_mode <- !is.null(config$lake_polygons)
  sf::sf_use_s2(FALSE)

  cw_inform(glue::glue(
    "Delineating {nrow(sites)} site(s) across {dplyr::n_distinct(sites$group_id)} ",
    "group(s) [{if (is_lake_mode) 'lake' else 'point'} mode]..."
  ))

  results <- purrr::map(unique(sites$group_id), function(grp) {
    grp_sites    <- dplyr::filter(sites, group_id == grp)
    grp_manifest <- dplyr::filter(group_manifest, group_id == grp)
    grp_cache    <- grp_manifest$cache_dir[1]

    group_rasters <- load_group_rasters(grp_cache, grp)

    purrr::map(seq_len(nrow(grp_sites)), function(j) {
      site <- grp_sites[j, ]
      if (is_lake_mode) {
        delineate_engine_lake_site(
          site = site, grp_cache = grp_cache, config = config,
          output_dir = config$output_dir, min_cells = min_cells
        )
      } else {
        delineate_engine_point_site(
          site = site, grp_cache = grp_cache, group_rasters = group_rasters,
          config = config, output_dir = config$output_dir,
          snap_dist = snap_dist, min_cells = min_cells
        )
      }
    }) |>
      dplyr::bind_rows()
  }) |>
    dplyr::bind_rows()

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

# -- Point mode ----------------------------------------------------------------

#' Delineate one site's catchment, point-pour-point (Jenson snap) mode
delineate_engine_point_site <- function(
  site, grp_cache, group_rasters, config, output_dir, snap_dist, min_cells
) {
  sid <- site$site_id
  site_dir <- site_output_dir(output_dir, sid)
  catchment_path <- fs::path(site_dir, "catchment.gpkg")

  if (cache_exists(catchment_path)) {
    cw_inform(glue::glue("Site '{sid}': catchment.gpkg found, skipping."))
    catchment <- sf::st_read(catchment_path, quiet = TRUE)
    cells <- catchment$n_cells[1]
    res_m <- terra::res(group_rasters$flow_accum)[1]
    km2 <- round((cells * res_m^2) / 1e6, 4)
    return(tibble::tibble(
      site_id = sid, status = "skipped (cached)",
      catchment_cells = cells, catchment_km2 = km2,
      flagged = km2 < (min_cells * res_m^2 / 1e6),
      flag_reason = if (km2 < (min_cells * res_m^2 / 1e6)) {
        "catchment smaller than min_cells threshold"
      } else {
        NA_character_
      }
    ))
  }

  cw_inform(glue::glue("Site '{sid}': delineating catchment..."))

  tryCatch(
    {
      pour_point_shp <- write_pour_point_shp_dynamic(site, site_dir, config$working_crs)

      snapped_shp <- snap_pour_point(
        pour_point_shp = pour_point_shp, streams = group_rasters$streams,
        site_dir = site_dir, site_id = sid, snap_dist = snap_dist
      )

      watershed_tif <- delineate_watershed(
        snapped_shp = snapped_shp, flow_pointer = group_rasters$flow_pointer,
        site_dir = site_dir, site_id = sid
      )

      catchment_sf <- watershed_to_polygon_dynamic(
        watershed_tif = watershed_tif, site_dir = site_dir, site_id = sid,
        working_crs = config$working_crs
      )

      watershed_rast <- terra::rast(watershed_tif)
      n_cells <- sum(terra::values(watershed_rast) == 1, na.rm = TRUE)
      res_m <- terra::res(group_rasters$flow_accum)[1]
      km2 <- round((n_cells * res_m^2) / 1e6, 4)

      flagged <- n_cells < min_cells
      flag_reason <- if (flagged) {
        glue::glue("catchment is only {n_cells} cells ({km2} km2) — pour point may have snapped to wrong location")
      } else {
        NA_character_
      }
      if (flagged) {
        cw_warn(glue::glue("Site '{sid}': {flag_reason}"))
      } else {
        cw_inform(glue::glue("Site '{sid}': catchment = {n_cells} cells ({km2} km2)."))
      }

      catchment_sf <- catchment_sf |>
        dplyr::mutate(site_id = sid, n_cells = n_cells, area_km2 = km2, flagged = flagged)

      snapped_sf <- sf::st_read(snapped_shp, quiet = TRUE) |>
        sf::st_transform(config$working_crs) |>
        dplyr::mutate(site_id = sid)
      sf::st_write(snapped_sf, fs::path(site_dir, "pour_point.gpkg"), delete_dsn = TRUE, quiet = TRUE)

      sf::st_write(catchment_sf, catchment_path, delete_dsn = TRUE, quiet = TRUE)

      clip_rasters_to_catchment(
        catchment_sf = catchment_sf, group_rasters = group_rasters,
        site_dir = site_dir, site_id = sid
      )

      # Only clip NHN flowlines if the group actually has burned-in
      # flowlines cached (e.g. streams_burn$source != "none") — whole_domain/
      # OIH-style runs with no burn-in have nothing to clip here.
      flowlines_path <- fs::path(grp_cache, "flowlines.gpkg")
      if (cache_exists(flowlines_path)) {
        clip_flowlines_to_catchment(
          catchment_sf = catchment_sf, flowlines_path = flowlines_path,
          site_dir = site_dir, site_id = sid
        )
      }

      tibble::tibble(
        site_id = sid, status = "success", catchment_cells = n_cells,
        catchment_km2 = km2, flagged = flagged, flag_reason = flag_reason
      )
    },
    error = function(e) {
      cw_warn(glue::glue("Site '{sid}': delineation failed — {e$message}"))
      tibble::tibble(
        site_id = sid, status = paste("failed:", e$message),
        catchment_cells = NA_integer_, catchment_km2 = NA_real_,
        flagged = TRUE, flag_reason = e$message
      )
    }
  )
}

#' Write a single site's pour point as a .shp file, in the run's resolved
#' working CRS — CRS-dynamic variant of workflow/R/stream/
#' 05_delineate_sites.R's write_pour_point_shp() (which hardcodes EPSG:3979
#' to match MRDEM; this variant matches whatever DEM the run was resolved
#' against instead).
write_pour_point_shp_dynamic <- function(site, tmp_dir, working_crs) {
  pour_point_shp <- fs::path(tmp_dir, "pour_point.shp")
  sf::st_as_sf(site, coords = c("lon", "lat"), crs = 4326) |>
    sf::st_transform(working_crs) |>
    sf::st_write(pour_point_shp, delete_dsn = TRUE, quiet = TRUE)
  pour_point_shp
}

#' Convert a watershed raster to a catchment polygon, in the run's resolved
#' working CRS — CRS-dynamic variant of workflow/R/stream/
#' 05_delineate_sites.R's watershed_to_polygon().
watershed_to_polygon_dynamic <- function(watershed_tif, site_dir, site_id, working_crs) {
  catchment_shp <- fs::path(site_dir, "catchment_tmp.shp")

  whitebox::wbt_raster_to_vector_polygons(
    input = normalizePath(watershed_tif, mustWork = TRUE),
    output = normalizePath(catchment_shp, mustWork = FALSE)
  )

  if (!fs::file_exists(catchment_shp)) {
    cw_abort(glue::glue(
      "Site '{site_id}': wbt_raster_to_vector_polygons() did not produce ",
      "output. The watershed raster may be all-NoData."
    ))
  }

  catchment_sf <- sf::st_read(catchment_shp, quiet = TRUE) |>
    dplyr::filter(VALUE == 1) |>
    sf::st_union() |>
    sf::st_as_sf() |>
    dplyr::rename(geometry = x) |>
    sf::st_transform(working_crs)

  fs::dir_ls(site_dir, glob = "catchment_tmp.*") |> fs::file_delete()

  if (nrow(catchment_sf) == 0 || sf::st_is_empty(catchment_sf$geometry[1])) {
    cw_abort(glue::glue(
      "Site '{site_id}': catchment polygon is empty after filtering VALUE == 1."
    ))
  }

  catchment_sf
}

# -- Lake mode -------------------------------------------------------------

#' Delineate one site's catchment, lake-polygon pour point mode
#'
#' Adapted from workflow/R/lake/03_delineate_lakes.R's
#' delineate_single_lake() — same buffer-rasterize-watershed logic and the
#' same documented no-trim-before-wbt_watershed() rule, generalized to read
#' from group_manifest's per-group cache_dir (flow_pointer.tif /
#' flow_accum.tif / streams.tif / dem.tif — the canonical names
#' 02_prepare_terrain.R always produces) instead of a flat project
#' cache_dir with d8_pntr.tif.
delineate_engine_lake_site <- function(site, grp_cache, config, output_dir, min_cells) {
  sid <- site$site_id
  lake_name <- site$lake_name
  site_dir <- site_output_dir(output_dir, sid)
  ensure_dir(site_dir)

  catchment_path <- fs::path(site_dir, "catchment.gpkg")
  if (cache_exists(catchment_path)) {
    cw_inform(glue::glue("Site '{sid}': catchment.gpkg found — skipping"))
    return(tibble::tibble(
      site_id = sid, status = "skipped (cached)",
      catchment_cells = NA_integer_, catchment_km2 = NA_real_,
      flagged = FALSE, flag_reason = NA_character_
    ))
  }

  lake_polys <- config$lake_polygons
  lake_poly <- lake_polys[lake_polys$matched_lake == lake_name, ]
  if (nrow(lake_poly) == 0) {
    cw_warn(glue::glue("Site '{sid}': no polygon matched for lake_name = '{lake_name}' — skipping"))
    return(tibble::tibble(
      site_id = sid, status = "skipped (no polygon)",
      catchment_cells = NA_integer_, catchment_km2 = NA_real_,
      flagged = TRUE, flag_reason = glue::glue("no polygon matched for lake '{lake_name}'")
    ))
  }
  if (nrow(lake_poly) > 1) {
    lake_poly <- lake_poly[1, ]
    cw_warn(glue::glue("Site '{sid}': multiple polygons matched — using first only"))
  }

  tryCatch(
    {
      pointer_path <- fs::path(grp_cache, "flow_pointer.tif")
      accum_path   <- fs::path(grp_cache, "flow_accum.tif")
      streams_path <- fs::path(grp_cache, "streams.tif")
      dem_path     <- fs::path(grp_cache, "dem.tif")
      d8_template  <- terra::rast(pointer_path)

      watershed_path <- fs::path(site_dir, "watershed.tif")
      pourpoint_path <- fs::path(site_dir, "lake_pourpoint.tif")

      # CRITICAL: use d8_template (full extent, no trim) as the raster
      # template — wbt_watershed() requires the pour point raster to share
      # extent and resolution with the D8 pointer; any trimming causes a
      # mismatch. Trim only AFTER delineation.
      lake_proj     <- terra::project(lake_poly, terra::crs(d8_template))
      lake_buffered <- terra::buffer(lake_proj, width = config$lake_buffer_m)
      lake_rast     <- terra::rasterize(lake_buffered, d8_template, field = 1, background = NA)
      terra::writeRaster(lake_rast, pourpoint_path, overwrite = TRUE)

      whitebox::wbt_watershed(
        d8_pntr  = fs::path_abs(pointer_path),
        pour_pts = fs::path_abs(pourpoint_path),
        output   = fs::path_abs(watershed_path)
      )

      watershed_rast <- terra::rast(watershed_path) |> terra::trim()
      terra::writeRaster(watershed_rast, watershed_path, overwrite = TRUE)

      catchment_sf <- watershed_rast |>
        terra::as.polygons() |>
        terra::project(terra::crs(d8_template)) |>
        sf::st_as_sf() |>
        dplyr::mutate(site_id = sid)

      if (nrow(catchment_sf) == 0 || sf::st_is_empty(sf::st_union(catchment_sf))) {
        cw_abort(glue::glue("Site '{sid}': watershed polygon is empty"))
      }

      sf::st_write(catchment_sf, catchment_path, delete_dsn = TRUE, quiet = TRUE)

      catchment_vect <- terra::vect(catchment_sf)
      rasters_to_clip <- list(
        dem          = terra::rast(dem_path),
        flow_pointer = terra::rast(pointer_path),
        flow_accum   = terra::rast(accum_path),
        streams      = terra::rast(streams_path)
      )
      purrr::iwalk(rasters_to_clip, function(r, name) {
        out <- r |> terra::crop(catchment_vect, snap = "out") |> terra::mask(catchment_vect)
        terra::writeRaster(
          out, fs::path(site_dir, paste0(name, ".tif")),
          overwrite = TRUE, gdal = c("COMPRESS=LZW", "BIGTIFF=IF_SAFER")
        )
      })

      n_cells <- sum(terra::values(watershed_rast) == 1, na.rm = TRUE)
      area_km2 <- round(sum(as.numeric(sf::st_area(catchment_sf))) / 1e6, 4)
      flagged <- n_cells < min_cells
      flag_reason <- if (flagged) glue::glue("catchment only {n_cells} cells ({area_km2} km2)") else NA_character_

      if (flagged) {
        cw_warn(glue::glue("Site '{sid}': FLAGGED — {flag_reason}"))
      } else {
        cw_inform(glue::glue("Site '{sid}': {n_cells} cells, {area_km2} km2"))
      }

      tibble::tibble(
        site_id = sid, status = "success", catchment_cells = n_cells,
        catchment_km2 = area_km2, flagged = flagged, flag_reason = flag_reason
      )
    },
    error = function(e) {
      cw_warn(glue::glue("Site '{sid}': delineation failed — {e$message}"))
      tibble::tibble(
        site_id = sid, status = paste("failed:", e$message),
        catchment_cells = NA_integer_, catchment_km2 = NA_real_,
        flagged = TRUE, flag_reason = e$message
      )
    }
  )
}
