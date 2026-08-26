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
#   nhn_index_path, nhn_raw_dir
#   lake_polygons  = sf/SpatVector | NULL   (NULL = point pour point mode)
#   lake_buffer_m
#   grouping       = list(strategy = "whole_domain"|"hydrobasins", ...)
#   loi_layers
#
# Dependencies: sf, terra, fs, glue, cli (via utils.R)
# ---------------------------------------------------------------------------

# Null-coalescing operator — same definition workflow/R/stream/
# 06_hydroweight_attributes.R uses locally. Defined here (not added to the
# shared utils.R) so the engine tree stays self-contained and doesn't touch
# any file outside workflow/R/engine/.
`%||%` <- function(x, y) if (!is.null(x)) x else y

#' Validate a run_config and resolve its working CRS
#'
#' Checks that exactly one usable terrain-conditioning tier is present
#' (flow_pointer > flow_direction > dem, highest already-conditioned wins),
# that streams_burn is only configured when it can actually apply (raw dem
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

  # -- Terrain tier: exactly one of flow_pointer / flow_direction / dem must
  # be usable as the highest-available-conditioning source. dem may ALSO be
  # supplied alongside flow_pointer/flow_direction purely to provide an
  # elevation surface for per-site clipping/output/hydroweight — that's not
  # a second "tier", so only flow_pointer/flow_direction compete for tier
  # selection; dem is always allowed to coexist.
  has_pointer   <- !is.null(config$flow_pointer[["path"]])
  has_direction <- !is.null(config$flow_direction[["path"]])
  has_dem       <- !is.null(config$dem[["path"]])

  if (!has_pointer && !has_direction && !has_dem) {
    cw_abort(paste(
      "run_config must supply at least one terrain input:",
      "flow_pointer$path, flow_direction$path, or dem$path."
    ))
  }

  terrain_tier <- if (has_pointer) {
    "flow_pointer"
  } else if (has_direction) {
    "flow_direction"
  } else {
    "dem"
  }

  # -- streams_burn only applies to the raw-dem tier — burning has to
  # happen BEFORE flow direction is derived (breach uses the burned DEM as
  # its input), so it's meaningless once a pre-conditioned flow_pointer/
  # flow_direction is supplied directly. Warn (not abort) rather than
  # silently ignoring a config the user may have copy-pasted from another
  # project without adjusting.
  burn_source <- config$streams_burn[["source"]] %||% "none"
  if (terrain_tier != "dem" && burn_source != "none") {
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

  # -- Grouping strategy
  strategy <- config$grouping[["strategy"]] %||% "whole_domain"
  if (!strategy %in% c("whole_domain", "hydrobasins")) {
    cw_abort(glue::glue(
      "grouping$strategy must be 'whole_domain' or 'hydrobasins' — got '{strategy}'."
    ))
  }
  if (strategy == "hydrobasins" && is.null(config$grouping[["hydrobasins_dir"]])) {
    cw_abort("grouping$strategy = 'hydrobasins' requires grouping$hydrobasins_dir.")
  }
  config$grouping$strategy <- strategy

  # -- Resolve the terrain tier's own native CRS (always needed — used as
  # the default working_crs, and always logged even when overridden, so a
  # reprojection is visible rather than silent).
  crs_source_path <- switch(
    terrain_tier,
    flow_pointer   = config$flow_pointer[["path"]],
    flow_direction = config$flow_direction[["path"]],
    dem            = config$dem[["path"]]
  )
  if (!fs::file_exists(crs_source_path)) {
    cw_abort(glue::glue(
      "Terrain input for tier '{terrain_tier}' not found: {crs_source_path}"
    ))
  }
  crs_desc <- terra::crs(terra::rast(crs_source_path), describe = TRUE)
  if (is.na(crs_desc$code)) {
    cw_abort(glue::glue(
      "Could not resolve an EPSG code from {crs_source_path} — ",
      "the raster's CRS may be undefined or unrecognized."
    ))
  }
  native_crs <- paste0(crs_desc$authority, ":", crs_desc$code)

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
