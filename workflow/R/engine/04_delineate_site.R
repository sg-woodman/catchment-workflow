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
#   clip_flowlines_to_catchment(). write_pour_point_shp() and
#   watershed_to_polygon() below are THIS file's own CRS-dynamic versions,
#   the ONLY implementation of either now — they were named
#   write_pour_point_shp_dynamic()/watershed_to_polygon_dynamic() until
#   2026-09, coexisting with EPSG:3979-hardcoded originals of the plain
#   base name in stream/delineate_sites.R (kept separate at the time to
#   avoid touching that "reused unmodified" file). Renamed after those
#   originals were confirmed to have zero live callers and were removed —
#   every project runs on the engine, which only ever called the _dynamic
#   versions, so the "non-dynamic" originals were pure dead code that
#   could (and did — see git history for the incident) silently diverge
#   from what any real pipeline actually runs, because a fix applied to
#   them had no effect on anything.
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
      pour_point_shp <- write_pour_point_shp(site, site_dir, config$working_crs)

      snapped_shp <- snap_pour_point(
        pour_point_shp = pour_point_shp, streams = group_rasters$streams,
        site_dir = site_dir, site_id = sid, snap_dist = snap_dist
      )

      watershed_tif <- delineate_watershed(
        snapped_shp = snapped_shp, flow_pointer = group_rasters$flow_pointer,
        site_dir = site_dir, site_id = sid
      )

      catchment_sf <- watershed_to_polygon(
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
      # pre-conditioned-terrain runs with no burn-in have nothing to clip
      # here. That's a silent, expected skip for a group with
      # burn_streams = FALSE by design — but if burn_streams IS TRUE for this
      # site's group and flowlines.gpkg still isn't there, that's not
      # "nothing to clip," it's Stage 2's streams-burn step having failed
      # to produce it, and streams.gpkg silently never gets written with
      # no trace at all (found 2026-08-31 — engine/03_prepare_streams_
      # burn.R used to compute flowlines in-memory and hand them straight
      # to burn_streams_into_dem() without ever persisting flowlines.gpkg,
      # so this branch always took the silent-skip path for every engine
      # project; fixed there, but this guard had no warning of its own to
      # catch a *future* instance of the same shape). Warn once per such
      # site rather than staying silent.
      flowlines_path <- fs::path(grp_cache, "flowlines.gpkg")
      if (cache_exists(flowlines_path)) {
        clip_flowlines_to_catchment(
          catchment_sf = catchment_sf, flowlines_path = flowlines_path,
          site_dir = site_dir, site_id = sid
        )
      } else if (isTRUE(site$burn_streams[1])) {
        cw_warn(glue::glue(
          "Site '{sid}': group '{site$group_id[1]}' has burn_streams = TRUE ",
          "but {flowlines_path} doesn't exist — streams.gpkg will NOT be ",
          "written for this site. Check Stage 2's streams-burn step for ",
          "this group."
        ))
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
#' working CRS. The sole implementation (see this file's header) — matches
#' whatever DEM the run was resolved against, rather than hardcoding
#' EPSG:3979.
write_pour_point_shp <- function(site, tmp_dir, working_crs) {
  pour_point_shp <- fs::path(tmp_dir, "pour_point.shp")
  sf::st_as_sf(site, coords = c("lon", "lat"), crs = 4326) |>
    sf::st_transform(working_crs) |>
    sf::st_write(pour_point_shp, delete_dsn = TRUE, quiet = TRUE)
  pour_point_shp
}

#' Convert a watershed raster to a catchment polygon, in the run's resolved
#' working CRS. The sole implementation (see this file's header) — matches
#' whatever DEM the run was resolved against, rather than hardcoding
#' EPSG:3979.
#'
#' Uses terra::as.polygons(dissolve = TRUE) directly on the trimmed
#' watershed raster, NOT whitebox::wbt_raster_to_vector_polygons() +
#' sf::st_union() (the previous approach, retired here — see git history).
#' The two are not equivalent for a watershed whose D8 trace connects into
#' the pour-point cell only DIAGONALLY (a valid 8-connected flow path):
#' wbt_raster_to_vector_polygons() groups such a cell into the SAME polygon
#' feature as the rest of the catchment, producing a single self-touching
#' ("bowtie") ring at the corner pinch — invalid per OGC rules. The old
#' code's sf::st_make_valid() call resolves that invalidity, but (confirmed
#' directly, this GEOS version only supports geos_method = "valid_structure",
#' no alternative) by silently DROPPING the smaller lobe entirely rather
#' than preserving it as a separate polygon part, instead of returning a
#' multi-part but valid geometry. Confirmed on real data: EMILY_TURKEY's
#' WABO7161/WABO7164 each lost exactly their own pour-point cell this way —
#' watershed.tif had 28 valid cells, the resulting catchment.gpkg polygon
#' covered only 27 (24300 m2 vs the correct 25200 m2), and the missing cell
#' was the pour point's own — leaving it ~25m outside its own catchment
#' polygon. Silent at delineation time; only surfaced later as an opaque
#' hydroweight::hydroweight() crash ("[writeRaster] there are no cell
#' values") when the DEM crop aligned to that undersized polygon didn't
#' cover the pour point at all. terra::as.polygons(dissolve = TRUE) groups
#' cells by the same 8-connectivity rule but returns a proper (valid)
#' MULTIPOLYGON preserving every connected lobe — confirmed directly on
#' WABO7161's actual watershed.tif: 28/28 cells recovered, pour point
#' correctly contained. terra::trim() first crops to the raster's non-NA
#' extent so as.polygons() only ever processes the small per-site window,
#' not the full group-extent raster watershed.tif is written at.
watershed_to_polygon <- function(watershed_tif, site_dir, site_id, working_crs) {
  watershed_rast <- terra::rast(watershed_tif)
  watershed_rast[watershed_rast != 1] <- NA
  watershed_rast <- terra::trim(watershed_rast)

  if (is.null(watershed_rast) || terra::ncell(watershed_rast) == 0) {
    cw_abort(glue::glue(
      "Site '{site_id}': watershed raster is all-NoData. ",
      "wbt_watershed() may have failed, or the pour point fell outside ",
      "the flow pointer extent."
    ))
  }

  catchment_sf <- terra::as.polygons(watershed_rast, dissolve = TRUE, values = FALSE) |>
    sf::st_as_sf() |>
    sf::st_transform(working_crs) |>
    # Defensive, not expected to fire — as.polygons() already returns a
    # valid MULTIPOLYGON for the diagonal-lobe case above — but cheap
    # insurance against any other GEOS-invalidity source downstream code
    # (hydroweight::hydroweight()'s own clip_region handling in particular)
    # would otherwise trip on.
    sf::st_make_valid()

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
