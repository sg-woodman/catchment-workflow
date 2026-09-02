# 01_build_group_manifest.R
# ---------------------------------------------------------------------------
# Builds the group_manifest for an engine run, per config$grouping$strategy.
# Every run produces the SAME manifest shape every reused stream-project
# module already relies on (group_id, aoi, cache_dir, n_sites,
# burn_streams) regardless of strategy — downstream modules (terrain prep,
# streams-burn, delineation, hydroweight) never special-case which strategy
# produced the manifest they're reading.
#
# Strategies:
#   "whole_domain" — the trivial case: one flat project-wide raster set,
#     no AOI cropping at all. One row: group_id = project_id, aoi = the
#     full extent of whichever terrain input was resolved in
#     00_resolve_config.R, cache_dir = config$cache_dir itself (no
#     per-group subfolder — matches the flat cache_dir convention a
#     whole-domain project uses throughout).
#   "hydrobasins" — delegates directly to the existing, unmodified
#     workflow/R/stream/group_sites.R::build_group_manifest() /
#     workflow/R/utils.R::build_group_aoi() (HydroBasins level-6 union per
#     user-assigned group_id). Requires sites already carry a group_id
#     column and requires workflow/R/utils.R + workflow/R/stream/
#     group_sites.R to be sourced by the caller.
#
# Dependencies: sf, dplyr, purrr, fs, glue, cli (via utils.R)
# ---------------------------------------------------------------------------

#' Build the group_manifest for an engine run
#'
#' Also assigns `group_id` and `burn_streams` onto config$sites where the
#' strategy implies them (whole_domain), then runs the same
#' validate_sites_tibble() check every reused stream-module function
#' expects, so config$sites is fully interchangeable with what those
#' functions already receive from a hydrobasins-grouped project.
#'
#' @param config Resolved config from resolve_engine_config()
#' @return A list with `sites` (validated tibble, group_id/burn_streams
#'   populated) and `group_manifest` (sf tibble: group_id, burn_streams,
#'   buffer_m, aoi, cache_dir, n_sites)
build_engine_group_manifest <- function(config) {
  ensure_dir(config$output_dir)
  ensure_dir(config$cache_dir)

  burn_streams_flag <- config$streams_burn$source != "none"

  if (config$grouping$strategy == "whole_domain") {
    sites <- config$sites
    if (!"group_id" %in% names(sites)) {
      sites$group_id <- config$project_id
    }
    if (!"burn_streams" %in% names(sites)) {
      sites$burn_streams <- burn_streams_flag
    }
    sites <- validate_sites_tibble(sites)

    cw_inform(glue::glue(
      "Grouping strategy 'whole_domain': {nrow(sites)} site(s) in one group ",
      "('{config$project_id}'), cache_dir = {config$cache_dir}."
    ))

    aoi <- whole_domain_aoi(config)

    group_manifest <- tibble::tibble(
      group_id     = config$project_id,
      burn_streams = burn_streams_flag,
      buffer_m     = NA_real_,
      aoi          = aoi,
      cache_dir    = config$cache_dir,
      n_sites      = nrow(sites)
    ) |>
      sf::st_as_sf(sf_column_name = "aoi", crs = config$working_crs)

    purrr::walk(sites$site_id, function(sid) {
      ensure_dir(site_output_dir(config$output_dir, sid))
    })

    return(list(sites = sites, group_manifest = group_manifest))
  }

  # -- "hydrobasins" strategy — delegate to the existing, unmodified
  # HydroBasins grouping machinery. Requires workflow/R/stream/group_sites.R
  # already sourced (defines build_group_manifest()).
  if (!exists("build_group_manifest", mode = "function")) {
    cw_abort(paste(
      "grouping$strategy = 'hydrobasins' requires",
      "workflow/R/stream/group_sites.R to be sourced first",
      "(defines build_group_manifest())."
    ))
  }

  # build_group_manifest()/build_group_aoi() (workflow/R/utils.R) hardcode
  # EPSG:3979 for the AOI they construct — not generalized to an arbitrary
  # working_crs. Fine when config$working_crs resolved to 3979 anyway (the
  # normal case: hydrobasins grouping is meant for a national terrain
  # mosaic natively in that CRS), but a mismatch here would silently tag
  # every group's AOI with the wrong CRS.
  if (!identical(config$working_crs, "EPSG:3979")) {
    cw_warn(glue::glue(
      "grouping$strategy = 'hydrobasins' builds group AOIs in EPSG:3979 ",
      "(workflow/R/utils.R::build_group_aoi(), not yet generalized), but ",
      "working_crs resolved to {config$working_crs} — group AOIs will be in ",
      "EPSG:3979 regardless. This combination hasn't been exercised; verify ",
      "cache/<group_id>/dem.tif etc. end up in the CRS you expect."
    ))
  }

  sites <- validate_sites_tibble(config$sites)

  group_manifest <- build_group_manifest(
    sites             = sites,
    output_dir        = config$output_dir,
    cache_dir         = config$cache_dir,
    hydrobasins_dir   = config$grouping$hydrobasins_dir,
    default_buffer_m  = config$grouping$default_buffer_m %||% 1000,
    hybas_level       = config$grouping$hybas_level %||% 6
  )

  list(sites = sites, group_manifest = group_manifest)
}

#' Build the trivial whole-domain AOI: the full extent of whichever terrain
#' input was resolved, as an sfc polygon in the working CRS
#'
#' @param config Resolved config from resolve_engine_config()
#' @return An sfc polygon (rectangle, the raster's own extent)
whole_domain_aoi <- function(config) {
  src_path <- switch(
    config$terrain_tier,
    flow_pointer   = config$flow_pointer[["path"]],
    flow_direction = config$flow_direction[["path"]],
    dem            = config$dem[["path"]]
  )
  r <- terra::rast(src_path)
  ext_poly <- terra::as.polygons(terra::ext(r), crs = terra::crs(r))
  sf::st_geometry(sf::st_as_sf(ext_poly)) |>
    sf::st_transform(config$working_crs)
}
