# 00_resolve_config.R
# ---------------------------------------------------------------------------
# Validates a run_config list and resolves the working CRS from whichever
# terrain input was actually supplied, instead of a hardcoded per-project
# constant. This is the mechanical entry point for the "run off whatever
# inputs are provided" design: every downstream engine module reads its
# behavior off config$*  presence/absence rather than a project-specific
# code path.
#
# run_config shape (see workflow/templates/run_engine_template.R for a full
# annotated skeleton):
#   project_id, output_dir, cache_dir, sites
#   dem            = list(path = ...)                      | NULL
#   flow_direction = list(path = ..., recode = matrix|NULL) | NULL
#   flow_pointer   = list(path = ...)                       | NULL
#   flow_accum     = list(path = ...)                       | NULL
#   crs            = "EPSG:####" | NULL   (NULL = match the terrain tier's
#                    own native CRS, the default; supplying a value forces
#                    every terrain raster to be reprojected to it instead —
#                    see 02_prepare_terrain.R, which already reprojects
#                    whenever a source's CRS differs from working_crs)
#   stream_threshold = integer
#   streams_burn   = list(source = "nhn_auto"|"supplied"|"none", path = ...)
#   lake_conditioning = list(source = "nhn_auto"|"supplied"|"none", path = ...,
#                    min_area_ha = 1, exclude_types = c("Watercourse"))
#                    (raw dem tier only, opt-in, default "none" — flattens
#                    known lakes into the DEM before D8 derivation so every
#                    site in the group is delineated from one lake-aware
#                    flow field, instead of reactively correcting individual
#                    bisected catchments after the fact post-hoc — see
#                    workflow/R/engine/prepare_lake_conditioning.R)
#   nhn_index_path, nhn_raw_dir
#   lake_polygons  = sf/SpatVector | NULL   (NULL = point pour point mode)
#   lake_buffer_m
#   grouping       = list(strategy = "whole_domain"|"hydrobasins"|"manual_groups", ...)
#                    "manual_groups" is for a project whose sites split
#                    across disjoint terrain products with no single raster
#                    covering them all (e.g. two non-overlapping regional
#                    DEM tiles) — top-level dem/flow_direction/flow_pointer
#                    must be NULL, and grouping$groups supplies one
#                    independent terrain source per group instead:
#                      groups = list(
#                        list(group_id = "NE", dem = list(path=...), flow_direction = list(path=..., recode=...)),
#                        list(group_id = "NC", dem = list(path=...), flow_direction = list(path=..., recode=...))
#                      )
#                    Site -> group assignment is auto-detected by terrain
#                    coverage in 01_build_group_manifest.R, not a manual
#                    sites$group_id column (see that file for why).
#   loi_layers
#
# Dependencies: sf, terra, fs, glue, cli (via utils.R)
# ---------------------------------------------------------------------------

# Null-coalescing operator — same definition
# workflow/R/stream/hydroweight_attributes.R uses locally. Defined here (not added to the
# shared utils.R) so the engine tree stays self-contained and doesn't touch
# any file outside workflow/R/engine/.
`%||%` <- function(x, y) if (!is.null(x)) x else y

#' Resolve which terrain tier is usable from a set of terrain inputs
#' (flow_pointer > flow_direction > dem, highest-conditioning wins).
#' Factored out so both the single global resolution path (whole_domain/
#' hydrobasins share one terrain source) and the per-group resolution path
#' ("manual_groups" — each group supplies its own) apply the exact same
#' precedence rule from one place.
#'
#' @param dem,flow_direction,flow_pointer Each NULL or list(path = ...)
#' @return "flow_pointer"/"flow_direction"/"dem", or NULL if none supplied
resolve_terrain_tier_from_inputs <- function(dem, flow_direction, flow_pointer) {
  has_pointer   <- !is.null(flow_pointer[["path"]])
  has_direction <- !is.null(flow_direction[["path"]])
  has_dem       <- !is.null(dem[["path"]])
  if (!has_pointer && !has_direction && !has_dem) {
    return(NULL)
  }
  if (has_pointer) "flow_pointer" else if (has_direction) "flow_direction" else "dem"
}

#' Resolve a raster path's native CRS as an "EPSG:####" string. Factored
#' out so both the single-terrain-source path and manual_groups' first-
#' group default read a CRS the same way.
resolve_native_crs <- function(path, terrain_tier) {
  if (!fs::file_exists(path)) {
    cw_abort(glue::glue(
      "Terrain input for tier '{terrain_tier}' not found: {path}"
    ))
  }
  crs_desc <- terra::crs(terra::rast(path), describe = TRUE)
  if (is.na(crs_desc$code)) {
    cw_abort(glue::glue(
      "Could not resolve an EPSG code from {path} — ",
      "the raster's CRS may be undefined or unrecognized."
    ))
  }
  paste0(crs_desc$authority, ":", crs_desc$code)
}

#' Validate a run_config and resolve its working CRS
#'
#' Checks that exactly one usable terrain-conditioning tier is present
#' (flow_pointer > flow_direction > dem, highest already-conditioned wins),
#' that streams_burn is only configured when it can actually apply (raw dem
#' tier only — burning into an already-conditioned flow direction/pointer
#' makes no sense, since burning has to happen before flow direction is
#' derived, not after), and that grouping.strategy is supported. Adds
#' `working_crs` (an "EPSG:####" string) to the returned config, read from
#' whichever terrain input is highest-tier — never hardcoded.
#'
#' @param config A run_config list (see file header)
#' @return The config list with `working_crs` added, and defaults filled in
#'   for optional fields left unset
resolve_engine_config <- function(config) {
  required_top <- c("project_id", "output_dir", "cache_dir", "sites")
  missing_top <- setdiff(required_top, names(config))
  if (length(missing_top) > 0) {
    cw_abort(glue::glue(
      "run_config is missing required field(s): {paste(missing_top, collapse = ', ')}"
    ))
  }

  # -- Grouping strategy — resolved early. "manual_groups" changes how
  # terrain inputs and working CRS are resolved below (each group supplies
  # its own terrain source instead of one shared globally), so the rest of
  # this function branches on it from here on.
  strategy <- config$grouping[["strategy"]] %||% "whole_domain"
  if (!strategy %in% c("whole_domain", "hydrobasins", "manual_groups")) {
    cw_abort(glue::glue(
      "grouping$strategy must be 'whole_domain', 'hydrobasins', or 'manual_groups' — got '{strategy}'."
    ))
  }
  if (strategy == "hydrobasins" && is.null(config$grouping[["hydrobasins_dir"]])) {
    cw_abort("grouping$strategy = 'hydrobasins' requires grouping$hydrobasins_dir.")
  }
  config$grouping$strategy <- strategy

  if (strategy == "manual_groups") {
    # -- "manual_groups": each group supplies its own independent terrain
    # source (dem/flow_direction/flow_pointer/flow_accum) instead of one
    # shared across every group like "whole_domain"/"hydrobasins" — for a
    # project whose sites split across disjoint terrain products (e.g. two
    # non-overlapping regional DEM tiles) with no single raster covering
    # them all. Top-level dem/flow_direction/flow_pointer must stay unset
    # to avoid an ambiguous dual configuration.
    if (!is.null(config$dem) || !is.null(config$flow_direction) || !is.null(config$flow_pointer)) {
      cw_abort(paste(
        "grouping$strategy = 'manual_groups' supplies terrain per group via",
        "grouping$groups — top-level dem/flow_direction/flow_pointer must",
        "be NULL/absent to avoid an ambiguous dual configuration."
      ))
    }
    groups <- config$grouping[["groups"]]
    if (is.null(groups) || length(groups) == 0) {
      cw_abort("grouping$strategy = 'manual_groups' requires a non-empty grouping$groups list.")
    }
    group_ids <- purrr::map_chr(groups, ~ .x[["group_id"]] %||% NA_character_)
    if (anyNA(group_ids) || any(group_ids == "")) {
      cw_abort("Every entry in grouping$groups must have a non-empty group_id.")
    }
    if (anyDuplicated(group_ids) > 0) {
      cw_abort(glue::glue(
        "Duplicate group_id(s) in grouping$groups: ",
        "{paste(unique(group_ids[duplicated(group_ids)]), collapse = ', ')}"
      ))
    }

    # Resolve each group's own terrain tier up front, same precedence rule
    # as the single-source path, and stash it back onto the group entry —
    # 02_prepare_terrain.R reads grouping$groups[[i]]$terrain_tier directly.
    groups <- purrr::map(groups, function(g) {
      tier <- resolve_terrain_tier_from_inputs(g$dem, g$flow_direction, g$flow_pointer)
      if (is.null(tier)) {
        cw_abort(glue::glue(
          "Group '{g$group_id}' (manual_groups) must supply at least one of ",
          "dem$path, flow_direction$path, or flow_pointer$path."
        ))
      }
      g$terrain_tier <- tier
      g
    })
    config$grouping$groups <- groups

    terrain_tier <- "manual_groups" # sentinel — real per-group tiers live in grouping$groups[[i]]$terrain_tier
  } else {
    # -- Terrain tier: exactly one of flow_pointer / flow_direction / dem
    # must be usable as the highest-available-conditioning source. dem may
    # ALSO be supplied alongside flow_pointer/flow_direction purely to
    # provide an elevation surface for per-site clipping/output/hydroweight
    # — that's not a second "tier", so only flow_pointer/flow_direction
    # compete for tier selection; dem is always allowed to coexist.
    terrain_tier <- resolve_terrain_tier_from_inputs(
      config$dem, config$flow_direction, config$flow_pointer
    )
    if (is.null(terrain_tier)) {
      cw_abort(paste(
        "run_config must supply at least one terrain input:",
        "flow_pointer$path, flow_direction$path, or dem$path."
      ))
    }
  }

  # -- streams_burn only applies to the raw-dem tier — burning has to
  # happen BEFORE flow direction is derived (breach uses the burned DEM as
  # its input), so it's meaningless once a pre-conditioned flow_pointer/
  # flow_direction is supplied directly. Warn (not abort) rather than
  # silently ignoring a config the user may have copy-pasted from another
  # project without adjusting. Skipped for "manual_groups" — each group
  # there dispatches on its OWN resolved tier in 02_prepare_terrain.R, so a
  # single blanket check here would be wrong the moment groups have
  # different tiers (e.g. one pre-conditioned, one raw dem).
  burn_source <- config$streams_burn[["source"]] %||% "none"
  if (!identical(terrain_tier, "manual_groups") && terrain_tier != "dem" && burn_source != "none") {
    cw_warn(glue::glue(
      "streams_burn$source = '{burn_source}' has no effect — terrain tier ",
      "is '{terrain_tier}', which is already conditioned. Burning only ",
      "applies when starting from a raw dem$path. Treating as 'none'."
    ))
    burn_source <- "none"
  }
  if (!burn_source %in% c("nhn_auto", "supplied", "none")) {
    cw_abort(glue::glue(
      "streams_burn$source must be one of 'nhn_auto', 'supplied', 'none' — got '{burn_source}'."
    ))
  }
  if (burn_source == "supplied" && is.null(config$streams_burn[["path"]])) {
    cw_abort("streams_burn$source = 'supplied' requires streams_burn$path.")
  }
  if (burn_source == "nhn_auto" &&
    (is.null(config$nhn_index_path) || is.null(config$nhn_raw_dir))) {
    cw_abort(paste(
      "streams_burn$source = 'nhn_auto' requires nhn_index_path and",
      "nhn_raw_dir to be set in run_config."
    ))
  }
  config$streams_burn$source <- burn_source

  # -- lake_conditioning: same restriction as streams_burn (raw dem tier
  # only — flattening has to happen before D8 derivation) and the same
  # source vocabulary, deliberately mirrored so the two options read the
  # same way in a run_config. Opt-in, default "none" — a project supplying
  # an already-conditioned flow_direction/flow_pointer should have zero
  # behavior change; that pre-conditioned surface is trusted as-is, which
  # is the whole point of supplying one. Same "manual_groups" skip as
  # streams_burn above, for the same reason.
  lake_source <- config$lake_conditioning[["source"]] %||% "none"
  if (!identical(terrain_tier, "manual_groups") && terrain_tier != "dem" && lake_source != "none") {
    cw_warn(glue::glue(
      "lake_conditioning$source = '{lake_source}' has no effect — terrain ",
      "tier is '{terrain_tier}', which is already conditioned. Lake ",
      "conditioning only applies when starting from a raw dem$path. ",
      "Treating as 'none'."
    ))
    lake_source <- "none"
  }
  if (!lake_source %in% c("nhn_auto", "supplied", "none")) {
    cw_abort(glue::glue(
      "lake_conditioning$source must be one of 'nhn_auto', 'supplied', 'none' — got '{lake_source}'."
    ))
  }
  if (lake_source == "supplied" && is.null(config$lake_conditioning[["path"]])) {
    cw_abort("lake_conditioning$source = 'supplied' requires lake_conditioning$path.")
  }
  if (lake_source == "nhn_auto" &&
    (is.null(config$nhn_index_path) || is.null(config$nhn_raw_dir))) {
    cw_abort(paste(
      "lake_conditioning$source = 'nhn_auto' requires nhn_index_path and",
      "nhn_raw_dir to be set in run_config (same fields streams_burn's",
      "nhn_auto uses)."
    ))
  }
  config$lake_conditioning$source <- lake_source
  config$lake_conditioning$min_area_ha <- config$lake_conditioning[["min_area_ha"]] %||% 1
  config$lake_conditioning$exclude_types <- config$lake_conditioning[["exclude_types"]] %||% c("Watercourse")
  # Lake search is scoped to a buffer around the group's own SITES, not its
  # full terrain-conditioning AOI — confirmed directly that the latter is
  # wildly disproportionate for a HydroBasins group (58,928 sq km group AOI
  # vs. 2,663 sq km for the sites' own bare bounding box, for one real
  # group). 20 km default is generous relative to observed catchment sizes
  # (largest confirmed so far ~16 sq km).
  config$lake_conditioning$site_buffer_m <- config$lake_conditioning[["site_buffer_m"]] %||% 20000

  # -- Resolve the terrain tier's own native CRS (always needed — used as
  # the default working_crs, and always logged even when overridden, so a
  # reprojection is visible rather than silent). For "manual_groups" there
  # is no single terrain source to read; the FIRST group's own native CRS
  # is used as the default instead — 02_prepare_terrain.R already
  # reprojects any source whose CRS doesn't match working_crs, so a later
  # group in a different CRS is handled the same way any other mismatch is,
  # not a special case.
  if (identical(terrain_tier, "manual_groups")) {
    first_group <- config$grouping$groups[[1]]
    crs_source_path <- switch(
      first_group$terrain_tier,
      flow_pointer   = first_group$flow_pointer[["path"]],
      flow_direction = first_group$flow_direction[["path"]],
      dem            = first_group$dem[["path"]]
    )
  } else {
    crs_source_path <- switch(
      terrain_tier,
      flow_pointer   = config$flow_pointer[["path"]],
      flow_direction = config$flow_direction[["path"]],
      dem            = config$dem[["path"]]
    )
  }
  native_crs <- resolve_native_crs(crs_source_path, terrain_tier)

  # -- Working CRS: config$crs, if supplied, overrides the terrain tier's
  # native CRS — 02_prepare_terrain.R reprojects every terrain raster to
  # it (it already reprojects whenever a source's CRS differs from
  # working_crs, so an override here needs no separate handling there).
  # Default (config$crs unset) matches the terrain source exactly, so no
  # reprojection happens unless the caller explicitly asks for one.
  if (!is.null(config[["crs"]])) {
    working_crs <- config[["crs"]]
    if (is.na(sf::st_crs(working_crs))) {
      cw_abort(glue::glue("config$crs '{working_crs}' could not be parsed as a valid CRS."))
    }
    if (!identical(working_crs, native_crs)) {
      cw_inform(glue::glue(
        "config$crs = {working_crs} overrides the terrain tier's native CRS ",
        "({native_crs}, from {fs::path_file(crs_source_path)}) — every terrain ",
        "raster will be reprojected."
      ))
    }
  } else {
    working_crs <- native_crs
    if (identical(terrain_tier, "manual_groups") && length(config$grouping$groups) > 1) {
      cw_inform(glue::glue(
        "config$crs not set — working CRS defaults to manual group ",
        "'{config$grouping$groups[[1]]$group_id}''s native CRS ({native_crs}). ",
        "Any other group whose terrain source uses a different CRS will be ",
        "reprojected to match in 02_prepare_terrain.R (same as any other ",
        "CRS mismatch)."
      ))
    }
  }

  check_crs_suitability(working_crs)

  cw_inform(glue::glue(
    "Config resolved: project '{config$project_id}', terrain tier = ",
    "'{terrain_tier}', working CRS = {working_crs}, grouping = '{strategy}'."
  ))

  config$working_crs   <- working_crs
  config$native_crs    <- native_crs
  config$terrain_tier  <- terrain_tier
  config$stream_threshold <- config$stream_threshold %||% 1000L
  config$lake_buffer_m    <- config$lake_buffer_m %||% 30

  config
}

#' Warn if a working CRS looks unsuitable for watershed delineation
#'
#' WhiteboxTools' D8 flow-direction/accumulation algorithms (and this
#' workflow's own distance/area logic — breach distance in cells, buffer
#' widths in metres, stream thresholds by cell count, output areas)
#' assume a projected CRS with uniform, metre-based square cells. Checked
#' regardless of whether working_crs came from the terrain source's own
#' native CRS (the default) or an explicit config$crs override — a bad
#' choice is a bad choice either way. Warns rather than aborts: the
#' caller may have a deliberate reason (e.g. testing), and this workflow
#' has no way to know for certain what's "ideal" for a study area it
#' doesn't recognize.
#'
#' @param working_crs Character. "EPSG:####" (or any sf::st_crs()-parseable
#'   string)
#' @return invisibly TRUE. Called for side effects (cw_warn()).
check_crs_suitability <- function(working_crs) {
  crs_obj <- sf::st_crs(working_crs)

  if (sf::st_is_longlat(crs_obj)) {
    cw_warn(glue::glue(
      "Working CRS ({working_crs}) is geographic (degrees), not projected. ",
      "WhiteboxTools' flow-direction/accumulation algorithms assume a ",
      "projected CRS with uniform, metre-based cell size — a geographic ",
      "CRS will produce distorted or incorrect watershed delineation. ",
      "Set config$crs to an appropriate projected CRS for your study area."
    ))
    return(invisible(TRUE))
  }

  units <- crs_obj$units_gdal
  if (!is.null(units) && !tolower(units) %in% c("metre", "meter", "m")) {
    cw_warn(glue::glue(
      "Working CRS ({working_crs}) uses linear units '{units}', not metres. ",
      "Distance/area parameters throughout this workflow (breach distance, ",
      "buffer widths, stream thresholds by cell count, output areas) assume ",
      "metres — results will be silently wrong in another unit. Set ",
      "config$crs to a metre-based projected CRS."
    ))
  }

  invisible(TRUE)
}
