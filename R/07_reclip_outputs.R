# 07_reclip_outputs.R
# ---------------------------------------------------------------------------
# Re-clips flowlines and rasters using the clipped catchment polygons produced
# by 06_remove_upstream.R. Rather than overwriting the original outputs from
# 05_delineate_sites.R (streams.gpkg, dem.tif, flow_accum.tif, etc.), this
# writes new files with a "_clipped" suffix so both versions remain available
# for comparison.
#
# Accepts either:
#   - The combined all_catchments_clipped.gpkg (default — reads per-site
#     geometry from the site_id column)
#   - A per-site catchment_clipped.gpkg path
#   - An sf object directly
#
# Outputs (per site, written to output/<site_id>/):
#   streams_clipped.gpkg : NHN flowlines clipped to catchment_clipped
#   dem_clipped.tif, dem_breached_clipped.tif, flow_pointer_clipped.tif,
#   flow_accum_clipped.tif, streams_clipped.tif, hillshade_clipped.tif
#     : group rasters clipped and masked to catchment_clipped
#
# Run after 06_remove_upstream.R.
#
# Dependencies: sf, terra, dplyr, purrr, fs (via utils.R)
# ---------------------------------------------------------------------------

#' Re-clip flowlines and rasters using clipped catchments
#'
#' For each site, loads its clipped catchment geometry (from the combined
#' file, a per-site file, or a supplied sf object), then regenerates
#' streams_clipped.gpkg and all group raster *_clipped.tif outputs.
#'
#' @param sites          Validated sites tibble from validate_sites()
#' @param group_manifest sf tibble from build_group_manifest()
#' @param output_dir     Character. Root output directory
#' @param catchments     One of:
#'   - NULL (default): reads output/all_catchments_clipped.gpkg and matches
#'     rows by site_id
#'   - Character path to a GeoPackage: read and matched by site_id if a
#'     site_id column exists, otherwise treated as a single-site file
#'   - sf object: used directly, matched by site_id if present
#'
#' @return A tibble summarising re-clip results for all sites, with columns:
#'   site_id, status
reclip_outputs <- function(
  sites,
  group_manifest,
  output_dir,
  catchments = NULL
) {
  sf::sf_use_s2(FALSE)

  catchment_sf <- resolve_catchments_input(catchments, output_dir)

  cw_inform(glue::glue(
    "Re-clipping outputs for {nrow(sites)} site(s)..."
  ))

  results <- purrr::map(sites$site_id, function(sid) {
    reclip_site(
      site_id = sid,
      catchment_sf = catchment_sf,
      sites = sites,
      group_manifest = group_manifest,
      output_dir = output_dir
    )
  }) |>
    dplyr::bind_rows()

  failed <- dplyr::filter(results, status != "success")
  if (nrow(failed) > 0) {
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
#' @param catchments Argument as passed to reclip_outputs() — NULL, a
#'   character path, or an sf object
#' @param output_dir Character. Root output directory — used as the default
#'   location for all_catchments_clipped.gpkg when catchments is NULL
#'
#' @return sf polygon object with a site_id column
resolve_catchments_input <- function(catchments, output_dir) {
  if (is.null(catchments)) {
    default_path <- fs::path(output_dir, "all_catchments_clipped.gpkg")

    if (!cache_exists(default_path)) {
      cw_abort(glue::glue(
        "catchments = NULL but {default_path} not found. ",
        "Run 06_remove_upstream.R first, or supply a path/sf object directly."
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
    cw_abort(
      "catchments must be NULL, a character path, or an sf object."
    )
  }

  if (!"site_id" %in% names(catchment_sf)) {
    cw_abort(
      "Catchments input has no site_id column — cannot match to sites."
    )
  }

  sf::st_transform(catchment_sf, 3979)
}

# -- Single-site re-clip --------------------------------------------------------

#' Re-clip flowlines and rasters for a single site
#'
#' @param site_id        Character. Site identifier
#' @param catchment_sf   sf polygon object with a site_id column, from
#'   resolve_catchments_input()
#' @param sites          Validated sites tibble from validate_sites()
#' @param group_manifest sf tibble from build_group_manifest()
#' @param output_dir     Character. Root output directory
#'
#' @return Single-row tibble with columns: site_id, status
reclip_site <- function(
  site_id,
  catchment_sf,
  sites,
  group_manifest,
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

  # Identify the site's group to load the correct group rasters and flowlines
  site <- dplyr::filter(sites, site_id == !!site_id)
  if (nrow(site) == 0) {
    cw_warn(glue::glue(
      "Site '{site_id}': not found in sites tibble — skipping."
    ))
    return(tibble::tibble(site_id = site_id, status = "skipped (not in sites)"))
  }

  grp <- site$group_id[1]
  grp_manifest <- dplyr::filter(group_manifest, group_id == !!grp)

  if (nrow(grp_manifest) == 0) {
    cw_warn(glue::glue(
      "Site '{site_id}': group '{grp}' not found in group_manifest — skipping."
    ))
    return(tibble::tibble(
      site_id = site_id,
      status = "skipped (group not found)"
    ))
  }

  grp_cache <- grp_manifest$cache_dir

  tryCatch(
    {
      group_rasters <- load_group_rasters(grp_cache, grp)

      clip_flowlines_to_catchment_clipped(
        catchment_sf = focal,
        flowlines_path = fs::path(grp_cache, "flowlines.gpkg"),
        site_dir = site_dir,
        site_id = site_id
      )

      clip_rasters_to_catchment_clipped(
        catchment_sf = focal,
        group_rasters = group_rasters,
        site_dir = site_dir,
        site_id = site_id
      )

      tibble::tibble(site_id = site_id, status = "success")
    },
    error = function(e) {
      cw_warn(glue::glue("Site '{site_id}': re-clip failed — {e$message}"))
      tibble::tibble(site_id = site_id, status = paste("failed:", e$message))
    }
  )
}

# -- Flowlines clipping (clipped catchment) -------------------------------------

#' Clip NHN flowlines to a clipped catchment polygon
#'
#' Same logic as clip_flowlines_to_catchment() from 05_delineate_sites.R, but
#' writes to streams_clipped.gpkg rather than streams.gpkg so the original
#' (unclipped-catchment) flowlines remain available for comparison.
#'
#' @param catchment_sf   sf polygon. Clipped catchment boundary in EPSG:3979
#' @param flowlines_path Character. Path to group flowlines.gpkg
#' @param site_dir       Character. Site output directory path
#' @param site_id        Character. Site identifier (for log messages)
#' @return Invisibly NULL. Called for side effects.
clip_flowlines_to_catchment_clipped <- function(
  catchment_sf,
  flowlines_path,
  site_dir,
  site_id
) {
  out_path <- fs::path(site_dir, "streams_clipped.gpkg")

  if (!cache_exists(flowlines_path)) {
    cw_warn(glue::glue(
      "Site '{site_id}': group flowlines.gpkg not found at {flowlines_path}. ",
      "Writing empty streams_clipped.gpkg."
    ))
    sf::st_sf(geometry = sf::st_sfc(crs = 3979)) |>
      sf::st_write(out_path, delete_dsn = TRUE, quiet = TRUE)
    return(invisible(NULL))
  }

  flowlines <- sf::st_read(flowlines_path, quiet = TRUE)

  if (nrow(flowlines) == 0) {
    cw_warn(glue::glue(
      "Site '{site_id}': group flowlines.gpkg is empty. ",
      "Writing empty streams_clipped.gpkg."
    ))
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
    "Site '{site_id}': streams_clipped.gpkg written ({nrow(clipped)} features)."
  ))

  invisible(NULL)
}

# -- Raster clipping (clipped catchment) -----------------------------------------

#' Clip and mask group rasters to a clipped catchment polygon
#'
#' Same logic as clip_rasters_to_catchment() from 05_delineate_sites.R, but
#' writes to "<name>_clipped.tif" rather than "<name>.tif" so the original
#' (unclipped-catchment) rasters remain available for comparison. Always
#' overwrites — does not skip if the *_clipped.tif already exists, since this
#' step is intended to be rerun whenever the clipped catchment changes.
#'
#' @param catchment_sf  sf polygon. Clipped catchment boundary in EPSG:3979
#' @param group_rasters Named list of SpatRaster from load_group_rasters()
#' @param site_dir      Character. Site output directory path
#' @param site_id       Character. Site identifier (for log messages)
#' @return Invisibly NULL. Called for side effects.
clip_rasters_to_catchment_clipped <- function(
  catchment_sf,
  group_rasters,
  site_dir,
  site_id
) {
  catchment_vect <- terra::vect(catchment_sf)

  purrr::iwalk(group_rasters, function(rast, name) {
    out_path <- fs::path(site_dir, paste0(name, "_clipped.tif"))

    clipped <- rast |>
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

# Using the combined output from 06 (default)
results <- reclip_outputs(sites, group_manifest, output_dir)

# Using a specific per-site clipped catchment
results <- reclip_outputs(
  sites,
  group_manifest,
  output_dir,
  catchments = "output/site_123/catchment_clipped.gpkg"
)
