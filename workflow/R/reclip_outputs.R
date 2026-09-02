# reclip_outputs.R
# ---------------------------------------------------------------------------
# Re-clips flowlines and rasters using the clipped catchment polygons produced
# by remove_upstream.R (stream) or remove_upstream_lakes.R (lake).
#
# Rather than overwriting the original outputs, this writes new files with a
# "_clipped" suffix so both versions remain available for comparison.
#
# Works for BOTH stream and lake projects:
#
#   Stream projects: pass group_manifest (from group_sites.R). Source
#     rasters are loaded from each site's group cache directory. Flowlines
#     (streams_clipped.gpkg) are also generated.
#
#   Lake projects:  pass cache_dir (project-level flat cache — the retired
#     lake pipeline's own 02_prepare_oih_dem.R originally, engine/
#     02_prepare_terrain.R with "whole_domain" grouping now; same file
#     names either way). Source rasters (dem, flow_pointer, flow_accum,
#     streams) are loaded from there via load_lake_cache_rasters(). No
#     flowlines step for lake projects.
#
# Outputs (per site, written to output_dir/<site_id>/):
#   *_clipped.tif          : dem, flow_pointer, flow_accum, streams
#                            (+ dem_breached, hillshade for stream projects)
#   streams_clipped.gpkg   : NHN flowlines clipped to clipped catchment
#                            (stream projects only)
#
# Run after remove_upstream.R (stream) or remove_upstream_lakes.R (lake).
#
# Dependencies: sf, terra, dplyr, purrr, fs, cli (via utils.R)
# ---------------------------------------------------------------------------

#' Re-clip flowlines and rasters using clipped catchments
#'
#' @param sites          Validated sites tibble (needs site_id; group_id needed
#'                       for stream projects)
#' @param output_dir     Character. Root output directory
#' @param catchments     One of:
#'   - NULL (default): reads output_dir/all_catchments_clipped.gpkg
#'   - Character path to a GeoPackage (must have site_id column)
#'   - sf object with site_id column
#' @param group_manifest sf tibble from build_group_manifest(). Provide for
#'                       stream projects. NULL for lake projects.
#' @param cache_dir      Character. Project-level cache directory. Provide for
#'                       lake projects (used in place of group_manifest).
#'                       NULL for stream projects.
#'
#' @return A tibble with columns: site_id, status
reclip_outputs <- function(
  sites,
  output_dir,
  catchments     = NULL,
  group_manifest = NULL,
  cache_dir      = NULL
) {
  if (is.null(group_manifest) && is.null(cache_dir)) {
    cw_abort(
      "Provide either group_manifest (stream projects) or cache_dir (lake projects)."
    )
  }

  sf::sf_use_s2(FALSE)
  catchment_sf <- resolve_catchments_input(catchments, output_dir)

  cw_inform(glue::glue(
    "Re-clipping outputs for {nrow(sites)} site(s)..."
  ))

  results <- purrr::map(sites$site_id, function(sid) {
    reclip_site(
      site_id        = sid,
      catchment_sf   = catchment_sf,
      sites          = sites,
      group_manifest = group_manifest,
      cache_dir      = cache_dir,
      output_dir     = output_dir
    )
  }) |>
    dplyr::bind_rows()

  # "failed: ..." (an error reclip_site()'s own tryCatch actually caught)
  # is distinct from "skipped (...)" (a legitimate, expected non-error
  # outcome — no matching geometry, site not in this group, etc.) — only
  # the former can indicate a structural bug.
  failed <- dplyr::filter(results, grepl("^failed:", status))

  if (nrow(failed) > 0) {
    distinct_errors <- unique(failed$status)

    # A structural bug (missing function, missing shared dependency) fails
    # every site the SAME way, not a per-site data issue — that's exactly
    # how a lake project's reclip stage once silently no-op'd for every
    # site in a real run (a referenced helper function was undefined;
    # see the relevant project's README for the fuller trail): every site
    # produced an identical "failed: could not find function..." row, each
    # individually just a cw_warn(), with nothing checking results$status
    # before later stages ran anyway. Abort loudly HERE instead — once,
    # immediately, with the actual error — rather than let it ride
    # silently to the end of a multi-hour run disguised as routine
    # per-site warnings.
    is_systemic <- length(distinct_errors) == 1 || nrow(failed) / nrow(results) > 0.5
    if (is_systemic) {
      cw_abort(glue::glue(
        "{nrow(failed)} of {nrow(results)} site(s) failed re-clipping, ",
        "{if (length(distinct_errors) == 1) 'all with the SAME error' else 'mostly with the same error'} ",
        "— this looks like a structural bug (missing function/dependency), ",
        "not per-site data issues. Not continuing to Stage 6/7 on top of a ",
        "reclip that didn't actually run.\n",
        "First error: {distinct_errors[1]}\n",
        "Affected site(s): {paste(head(failed$site_id, 10), collapse = ', ')}",
        "{if (nrow(failed) > 10) glue::glue(' (+{nrow(failed) - 10} more)') else ''}"
      ))
    }

    cw_warn(glue::glue(
      "{nrow(failed)} site(s) failed:\n",
      "{paste(failed$site_id, ':', failed$status, collapse = '\n')}"
    ))
  } else {
    cw_inform("All sites re-clipped successfully.")
  }

  results
}

# -- Input resolution ----------------------------------------------------------

#' Resolve the catchments argument into an sf object with a site_id column
#'
#' @param catchments Argument as passed to reclip_outputs()
#' @param output_dir Character. Root output directory
#' @return sf polygon object with site_id column, transformed to EPSG:3979
resolve_catchments_input <- function(catchments, output_dir) {
  if (is.null(catchments)) {
    default_path <- fs::path(output_dir, "all_catchments_clipped.gpkg")
    if (!cache_exists(default_path)) {
      cw_abort(glue::glue(
        "catchments = NULL but {default_path} not found. ",
        "Run upstream removal first, or supply a path/sf object directly."
      ))
    }
    cw_inform(glue::glue("Reading combined catchments from {default_path}"))
    catchment_sf <- sf::st_read(default_path, quiet = TRUE)
  } else if (is.character(catchments)) {
    if (!cache_exists(catchments)) {
      cw_abort(glue::glue("Catchments file not found: {catchments}"))
    }
    cw_inform(glue::glue("Reading catchments from {catchments}"))
    catchment_sf <- sf::st_read(catchments, quiet = TRUE)
  } else if (inherits(catchments, "sf")) {
    catchment_sf <- catchments
  } else {
    cw_abort("catchments must be NULL, a character path, or an sf object.")
  }

  if (!"site_id" %in% names(catchment_sf)) {
    cw_abort("Catchments input has no site_id column.")
  }

  # Transform to a common CRS for intersection ops — EPSG:3979 for stream
  # projects; lake projects (EPSG:3161) are handled by terra ops internally,
  # but we need a consistent CRS here for sf-based flowline clipping.
  sf::st_transform(catchment_sf, 3979)
}

# -- Lake project raster loading -----------------------------------------------

#' Load a lake project's project-level cached rasters for re-clipping
#'
#' Was called by reclip_site()'s lake branch but was, at one point, never
#' defined anywhere in the repo — every call raised "could not find
#' function 'load_lake_cache_rasters'", caught by reclip_site()'s own
#' tryCatch and downgraded to a per-site warning rather than an abort, so
#' reclip_outputs() "succeeded" while silently clipping zero rasters for
#' every lake project site, every run (see the relevant project's README
#' for the fuller trail of how this was found and fixed). Exactly 4
#' rasters, not 6 — matches this file's own header comment ("dem,
#' flow_pointer, flow_accum, streams... + dem_breached, hillshade for
#' stream projects" — i.e. NOT for lake projects), and matches what the
#' retired pre-engine lake pipeline actually produced (dem_clipped.tif,
#' flow_pointer_clipped.tif, flow_accum_clipped.tif, streams_clipped.tif;
#' never dem_breached_clipped.tif or hillshade_clipped.tif).
#'
#' @param cache_dir Character. Project-level (flat, not per-group) cache
#'   directory for a lake project.
#' @return Named list of 4 SpatRasters: dem, flow_pointer, flow_accum,
#'   streams.
load_lake_cache_rasters <- function(cache_dir) {
  raster_files <- list(
    dem = fs::path(cache_dir, "dem.tif"),
    flow_pointer = fs::path(cache_dir, "flow_pointer.tif"),
    flow_accum = fs::path(cache_dir, "flow_accum.tif"),
    streams = fs::path(cache_dir, "streams.tif")
  )

  missing <- purrr::keep(raster_files, function(p) !cache_exists(p))
  if (length(missing) > 0) {
    cw_abort(glue::glue(
      "Lake project cache missing required raster(s): ",
      "{paste(names(missing), collapse = ', ')} (looked in {cache_dir}). ",
      "Run terrain prep first."
    ))
  }

  purrr::map(raster_files, terra::rast)
}

# -- Single-site re-clip -------------------------------------------------------

#' Re-clip flowlines and rasters for a single site
#'
#' Dispatches to the stream path (group_manifest) or lake path (cache_dir)
#' based on which argument is provided.
#'
#' @param site_id        Character. Site identifier
#' @param catchment_sf   sf polygon object with site_id column
#' @param sites          Validated sites tibble
#' @param group_manifest sf tibble or NULL
#' @param cache_dir      Character path or NULL
#' @param output_dir     Character. Root output directory
#'
#' @return Single-row tibble with columns: site_id, status
reclip_site <- function(
  site_id,
  catchment_sf,
  sites,
  group_manifest,
  cache_dir,
  output_dir
) {
  site_dir <- site_output_dir(output_dir, site_id)

  focal <- dplyr::filter(catchment_sf, site_id == !!site_id)
  if (nrow(focal) == 0) {
    cw_warn(glue::glue(
      "Site '{site_id}': no matching geometry in catchments input — skipping."
    ))
    return(tibble::tibble(site_id = site_id, status = "skipped (no geometry)"))
  }
  if (nrow(focal) > 1) {
    focal <- focal |>
      sf::st_union() |>
      sf::st_as_sf() |>
      dplyr::rename(geometry = x) |>
      dplyr::mutate(site_id = site_id)
  }

  tryCatch(
    {
      if (!is.null(group_manifest)) {
        # -- Stream path: load rasters from group cache ----------------------
        site <- dplyr::filter(sites, site_id == !!site_id)
        if (nrow(site) == 0) {
          return(tibble::tibble(
            site_id = site_id, status = "skipped (not in sites)"
          ))
        }
        grp <- site$group_id[1]
        grp_manifest <- dplyr::filter(group_manifest, group_id == !!grp)
        if (nrow(grp_manifest) == 0) {
          return(tibble::tibble(
            site_id = site_id, status = "skipped (group not found)"
          ))
        }
        grp_cache <- grp_manifest$cache_dir
        group_rasters <- load_group_rasters(grp_cache, grp)

        clip_flowlines_to_catchment_clipped(
          catchment_sf   = focal,
          flowlines_path = fs::path(grp_cache, "flowlines.gpkg"),
          site_dir       = site_dir,
          site_id        = site_id
        )
      } else {
        # -- Lake path: load rasters from project-level cache ----------------
        group_rasters <- load_lake_cache_rasters(cache_dir)
        # No flowlines step for lake projects
      }

      # Raster clipping is identical for both paths
      clip_rasters_to_catchment_clipped(
        catchment_sf  = focal,
        group_rasters = group_rasters,
        site_dir      = site_dir,
        site_id       = site_id
      )

      tibble::tibble(site_id = site_id, status = "success")
    },
    error = function(e) {
      cw_warn(glue::glue("Site '{site_id}': re-clip failed — {e$message}"))
      tibble::tibble(site_id = site_id, status = paste("failed:", e$message))
    }
  )
}

# -- Flowlines clipping (clipped catchment) ------------------------------------

#' Clip NHN flowlines to a clipped catchment polygon
#'
#' Stream projects only. Writes streams_clipped.gpkg.
#'
#' @param catchment_sf   sf polygon. Clipped catchment boundary
#' @param flowlines_path Character. Path to group flowlines.gpkg
#' @param site_dir       Character. Site output directory path
#' @param site_id        Character. Site identifier (for log messages)
#' @return invisibly NULL. Called for side effects.
clip_flowlines_to_catchment_clipped <- function(
  catchment_sf,
  flowlines_path,
  site_dir,
  site_id
) {
  out_path <- fs::path(site_dir, "streams_clipped.gpkg")

  if (!cache_exists(flowlines_path)) {
    cw_warn(glue::glue(
      "Site '{site_id}': flowlines.gpkg not found at {flowlines_path}. ",
      "Writing empty streams_clipped.gpkg."
    ))
    sf::st_sf(geometry = sf::st_sfc(crs = 3979)) |>
      sf::st_write(out_path, delete_dsn = TRUE, quiet = TRUE)
    return(invisible(NULL))
  }

  flowlines <- sf::st_read(flowlines_path, quiet = TRUE)
  if (nrow(flowlines) == 0) {
    sf::st_write(flowlines, out_path, delete_dsn = TRUE, quiet = TRUE)
    return(invisible(NULL))
  }

  clipped <- tryCatch(
    sf::st_intersection(flowlines, sf::st_union(catchment_sf)),
    error = function(e) {
      cw_warn(glue::glue(
        "Site '{site_id}': error clipping flowlines — {e$message}. ",
        "Writing empty streams_clipped.gpkg."
      ))
      flowlines[0, ]
    }
  )

  clipped <- clipped[
    sf::st_geometry_type(clipped) %in% c("LINESTRING", "MULTILINESTRING"),
    ,
    drop = FALSE
  ]

  sf::st_write(clipped, out_path, delete_dsn = TRUE, quiet = TRUE)
  cw_inform(glue::glue(
    "Site '{site_id}': streams_clipped.gpkg written ({nrow(clipped)} features)."
  ))
  invisible(NULL)
}

# -- Raster clipping (clipped catchment) ---------------------------------------

#' Clip and mask rasters to a clipped catchment polygon
#'
#' Writes <name>_clipped.tif for each raster in group_rasters. Always
#' overwrites — does not skip if outputs exist.
#'
#' @param catchment_sf  sf polygon. Clipped catchment boundary
#' @param group_rasters Named list of SpatRasters (from load_group_rasters() or
#'                      load_lake_cache_rasters())
#' @param site_dir      Character. Site output directory path
#' @param site_id       Character. Site identifier (for log messages)
#' @return invisibly NULL. Called for side effects.
clip_rasters_to_catchment_clipped <- function(
  catchment_sf,
  group_rasters,
  site_dir,
  site_id
) {
  catchment_vect <- terra::vect(catchment_sf)
  if (length(group_rasters) > 0) {
    catchment_vect <- terra::project(catchment_vect, terra::crs(group_rasters[[1]]))
  }

  purrr::iwalk(group_rasters, function(r, name) {
    out_path <- fs::path(site_dir, paste0(name, "_clipped.tif"))

    clipped <- r |>
      terra::crop(catchment_vect, snap = "out") |>
      terra::mask(catchment_vect)

    terra::writeRaster(
      clipped,
      filename = out_path,
      overwrite = TRUE,
      gdal = c("COMPRESS=LZW", "BIGTIFF=IF_SAFER")
    )
    cw_inform(glue::glue("Site '{site_id}': {name}_clipped.tif written."))
  })

  invisible(NULL)
}
