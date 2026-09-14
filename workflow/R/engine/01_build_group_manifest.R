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
#   "manual_groups" — for a project whose sites split across disjoint
#     terrain products with no single raster covering them all (e.g. two
#     non-overlapping regional DEM tiles). One row per config$grouping$groups
#     entry, aoi = that group's OWN terrain source's full extent (same
#     "whole_domain" logic, just once per group instead of once overall),
#     cache_dir = cache_dir/<group_id> (per-group subfolder — there IS more
#     than one group here, unlike whole_domain's single flat one).
#     group_id is NOT a manual sites column (that's "hydrobasins"'
#     convention) — it's auto-detected by testing each site's coordinate
#     against every group's own terrain raster and assigning it to
#     whichever one actually has data there, aborting loudly if a site
#     matches zero or more than one group. This is deliberate: a project
#     needing this strategy exists specifically because no raster covers
#     every site, so silently guessing (or requiring 100+ rows of manual
#     bookkeeping that goes stale the moment a site is added) is worse
#     than an explicit, actionable error.
#
# Dependencies: sf, terra, dplyr, purrr, fs, glue, cli (via utils.R)
# ---------------------------------------------------------------------------

#' Build the group_manifest for an engine run
#'
#' Also assigns `group_id` and `burn_streams` onto config$sites where the
#' strategy implies them (whole_domain, manual_groups), then runs the same
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

  if (config$grouping$strategy == "manual_groups") {
    groups_cfg <- config$grouping$groups

    sites <- config$sites
    if ("group_id" %in% names(sites)) {
      cw_abort(paste(
        "grouping$strategy = 'manual_groups' auto-detects group_id from",
        "terrain coverage — remove the group_id column from sites",
        "(or use grouping$strategy = 'hydrobasins' to assign it manually)."
      ))
    }
    sites$group_id <- detect_site_group_ids(sites, groups_cfg, config$working_crs)
    if (!"burn_streams" %in% names(sites)) {
      sites$burn_streams <- burn_streams_flag
    }
    sites <- validate_sites_tibble(sites)

    group_manifest <- purrr::map(groups_cfg, function(g) {
      n_sites_g <- sum(sites$group_id == g$group_id)
      if (n_sites_g == 0) {
        cw_warn(glue::glue(
          "manual_groups: group '{g$group_id}' matched 0 sites — check its ",
          "terrain source actually covers this project's AOI."
        ))
      }
      aoi <- terrain_source_extent_aoi(
        dem = g$dem, flow_direction = g$flow_direction,
        flow_pointer = g$flow_pointer, terrain_tier = g$terrain_tier,
        working_crs = config$working_crs
      )
      tibble::tibble(
        group_id     = g$group_id,
        burn_streams = burn_streams_flag,
        buffer_m     = NA_real_,
        aoi          = aoi,
        cache_dir    = fs::path(config$cache_dir, g$group_id),
        n_sites      = n_sites_g
      )
    }) |>
      dplyr::bind_rows() |>
      sf::st_as_sf(sf_column_name = "aoi", crs = config$working_crs)

    purrr::walk(group_manifest$cache_dir, ensure_dir)
    purrr::walk(sites$site_id, function(sid) {
      ensure_dir(site_output_dir(config$output_dir, sid))
    })

    cw_inform(glue::glue(
      "Grouping strategy 'manual_groups': {nrow(sites)} site(s) across ",
      "{nrow(group_manifest)} group(s) (auto-detected by terrain coverage)."
    ))

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

#' Build the full-extent AOI (rectangle, working CRS) of whichever terrain
#' input drives a given tier — shared by "whole_domain" (config's single
#' terrain source) and "manual_groups" (each group's own independent
#' terrain source).
#'
#' @return An sfc polygon (rectangle, the raster's own extent)
terrain_source_extent_aoi <- function(dem, flow_direction, flow_pointer, terrain_tier, working_crs) {
  src_path <- switch(
    terrain_tier,
    flow_pointer   = flow_pointer[["path"]],
    flow_direction = flow_direction[["path"]],
    dem            = dem[["path"]]
  )
  r <- terra::rast(src_path)
  ext_poly <- terra::as.polygons(terra::ext(r), crs = terra::crs(r))
  sf::st_geometry(sf::st_as_sf(ext_poly)) |>
    sf::st_transform(working_crs)
}

#' Build the trivial whole-domain AOI: the full extent of whichever
#' terrain input was resolved, as an sfc polygon in the working CRS
#'
#' @param config Resolved config from resolve_engine_config()
#' @return An sfc polygon (rectangle, the raster's own extent)
whole_domain_aoi <- function(config) {
  terrain_source_extent_aoi(
    dem = config$dem, flow_direction = config$flow_direction,
    flow_pointer = config$flow_pointer, terrain_tier = config$terrain_tier,
    working_crs = config$working_crs
  )
}

#' Auto-detect which manual_groups group each site belongs to, by testing
#' coverage against each group's own terrain source (whichever raster
#' drives that group's resolved terrain_tier — see 00_resolve_config.R) —
#' a site belongs to the one group whose raster actually has data at that
#' coordinate. Aborts loudly listing any site matched to zero or more than
#' one group, rather than silently guessing: a project needing
#' "manual_groups" exists specifically because no single raster covers
#' every site, so an unmatched/ambiguous site is a real, actionable
#' configuration problem, not a one-off row to drop quietly.
#'
#' @param sites  tibble with site_id, lon, lat (WGS84 decimal degrees)
#' @param groups config$grouping$groups (each entry already carries a
#'   resolved terrain_tier — see 00_resolve_config.R)
#' @param working_crs "EPSG:####" string
#' @return Character vector of group_id, one per row of `sites`, in order
detect_site_group_ids <- function(sites, groups, working_crs) {
  pts <- sf::st_as_sf(sites, coords = c("lon", "lat"), crs = 4326) |>
    sf::st_transform(working_crs) |>
    terra::vect()

  coverage <- purrr::map(groups, function(g) {
    src_path <- switch(
      g$terrain_tier,
      flow_pointer   = g$flow_pointer[["path"]],
      flow_direction = g$flow_direction[["path"]],
      dem            = g$dem[["path"]]
    )
    r <- terra::rast(src_path)
    !is.na(terra::extract(r, pts)[[2]])
  })
  names(coverage) <- purrr::map_chr(groups, "group_id")
  coverage_mat <- as.data.frame(coverage)
  n_matches <- rowSums(coverage_mat)

  unmatched <- sites$site_id[n_matches == 0]
  if (length(unmatched) > 0) {
    cw_abort(glue::glue(
      "manual_groups: {length(unmatched)} site(s) matched no group's ",
      "terrain source (outside every group's coverage): ",
      "{paste(unmatched, collapse = ', ')}"
    ))
  }
  ambiguous <- sites$site_id[n_matches > 1]
  if (length(ambiguous) > 0) {
    cw_abort(glue::glue(
      "manual_groups: {length(ambiguous)} site(s) matched more than one ",
      "group's terrain source (overlapping coverage — group assignment is ",
      "ambiguous): {paste(ambiguous, collapse = ', ')}"
    ))
  }

  purrr::map_chr(seq_len(nrow(coverage_mat)), function(i) {
    names(coverage_mat)[which(as.logical(coverage_mat[i, ]))]
  })
}
