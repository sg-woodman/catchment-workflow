# 04_run_whitebox.R
# ---------------------------------------------------------------------------
# Runs the WhiteboxTools hydrological conditioning steps for each group.
# Uses the burned DEM if burn_streams = TRUE for the group, otherwise uses
# the raw DEM.
#
# Steps (in order):
#   1. wbt_breach_depressions_least_cost() — removes topographic sinks by
#      breaching rather than filling, which better preserves natural drainage
#   2. wbt_d8_pointer()                   — D8 flow direction raster
#   3. wbt_d8_flow_accumulation()         — flow accumulation raster
#
# All steps use absolute paths for WhiteboxTools compatibility.
# Outputs are cached — steps are skipped if output files already exist.
#
# Inputs:
#   group_manifest : sf tibble from build_group_manifest()
#   max_dist       : Maximum breach path distance in cells. Default 500.
#   flat_increment : Elevation increment for flat areas. Default NULL
#                    (WhiteboxTools determines automatically).
#
# Outputs (per group, written to cache/<group_id>/):
#   dem_breached.tif  : Depression-breached DEM
#   flow_pointer.tif  : D8 flow direction raster
#   flow_accum.tif    : D8 flow accumulation raster
#   streams.tif       : Binary stream network raster (1 = stream cell)
#   hillshade.tif     : Hillshade derived from breached DEM
#
# Dependencies: whitebox, fs, purrr, glue, cli (via utils.R)
# ---------------------------------------------------------------------------

#' Run WhiteboxTools hydrological conditioning for all groups
#'
#' Iterates over all groups in the manifest. For each group, determines the
#' correct input DEM (burned or raw), then runs breach depressions, flow
#' pointer, and flow accumulation. Skips any group where all three outputs
#' already exist in cache.
#'
#' @param group_manifest       sf tibble from build_group_manifest()
#' @param max_dist             Integer. Maximum breach path length in cells
#'                             passed to wbt_breach_depressions_least_cost().
#'                             Default 10 (= 300 m at 30 m resolution).
#'                             Increase for very flat terrain only.
#' @param flat_increment       Numeric or NULL. Small elevation increment added
#'                             to flat areas to enforce drainage. NULL lets
#'                             WhiteboxTools determine automatically. Default NULL.
#' @param fill                 Logical. Whether to fill remaining depressions
#'                             after breaching. Should be TRUE for flat or
#'                             northern landscapes with thermokarst terrain.
#'                             Default TRUE.
#' @param default_stream_threshold Integer. Flow accumulation threshold in cells
#'                             used to define streams. Overridden per group via
#'                             stream_threshold column in sites CSV. Default 1000.
#'
#' @return The input group_manifest invisibly. Called for side effects.
run_whitebox <- function(
    group_manifest,
    max_dist                 = 10,
    flat_increment           = NULL,
    fill                     = TRUE,
    default_stream_threshold = 1000
) {

  purrr::walk(seq_len(nrow(group_manifest)), function(i) {

    grp          <- group_manifest$group_id[i]
    grp_cache    <- group_manifest$cache_dir[i]
    burn_streams <- group_manifest$burn_streams[i]

    # Define output paths
    breached_path   <- fs::path(grp_cache, "dem_breached.tif")
    pointer_path    <- fs::path(grp_cache, "flow_pointer.tif")
    accum_path      <- fs::path(grp_cache, "flow_accum.tif")
    streams_path    <- fs::path(grp_cache, "streams.tif")
    hillshade_path  <- fs::path(grp_cache, "hillshade.tif")

    # Resolve stream threshold: group override if present, else global default.
    # Guard against stream_threshold column being absent from the manifest
    # (e.g. if built from a sites object validated before this column was added).
    grp_threshold <- if ("stream_threshold" %in% names(group_manifest)) {
      group_manifest$stream_threshold[i]
    } else {
      NA_real_
    }
    stream_threshold <- if (length(grp_threshold) == 1 && !is.na(grp_threshold)) {
      as.integer(grp_threshold)
    } else {
      as.integer(default_stream_threshold)
    }

    # Skip group if all outputs already exist
    if (cache_exists(breached_path)  &&
        cache_exists(pointer_path)   &&
        cache_exists(accum_path)     &&
        cache_exists(streams_path)   &&
        cache_exists(hillshade_path)) {
      cw_inform(glue::glue(
        "Group '{grp}': Whitebox outputs found in cache, skipping."
      ))
      return(invisible(NULL))
    }

    # Select input DEM: use burned DEM if burn_streams = TRUE and it exists,
    # otherwise fall back to raw DEM with a warning
    dem_burned_path <- fs::path(grp_cache, "dem_burned.tif")
    dem_raw_path    <- fs::path(grp_cache, "dem.tif")

    if (burn_streams && cache_exists(dem_burned_path)) {
      input_dem <- dem_burned_path
      cw_inform(glue::glue(
        "Group '{grp}': using burned DEM as input."
      ))
    } else {
      if (burn_streams && !cache_exists(dem_burned_path)) {
        cw_warn(glue::glue(
          "Group '{grp}': burn_streams = TRUE but dem_burned.tif not found. ",
          "Falling back to raw DEM — run prepare_nhn_layers() first."
        ))
      }
      input_dem <- dem_raw_path
      cw_inform(glue::glue(
        "Group '{grp}': using raw DEM as input."
      ))
    }

    if (!cache_exists(input_dem)) {
      cw_abort(glue::glue(
        "Group '{grp}': input DEM not found at {input_dem}. ",
        "Run prepare_dem() before run_whitebox()."
      ))
    }

    # Run the three Whitebox steps in sequence
    run_whitebox_group(
      grp              = grp,
      input_dem        = input_dem,
      breached_path    = breached_path,
      pointer_path     = pointer_path,
      accum_path       = accum_path,
      streams_path     = streams_path,
      hillshade_path   = hillshade_path,
      max_dist         = max_dist,
      flat_increment   = flat_increment,
      fill             = fill,
      stream_threshold = stream_threshold
    )
  })

  invisible(group_manifest)
}

#' Run all Whitebox steps for a single group
#'
#' Internal function. Runs breach depressions, D8 pointer, and D8 flow
#' accumulation in sequence. Each step checks whether its output already
#' exists before running, allowing partial reruns if a step failed midway.
#'
#' @param grp              Character. Group identifier (for log messages)
#' @param input_dem        Character. Path to input DEM (burned or raw)
#' @param breached_path    Character. Output path for breached DEM
#' @param pointer_path     Character. Output path for flow pointer
#' @param accum_path       Character. Output path for flow accumulation
#' @param streams_path    Character. Output path for streams raster
#' @param hillshade_path  Character. Output path for hillshade raster
#' @param max_dist        Integer. Maximum breach path length in cells
#' @param flat_increment  Numeric or NULL. Flat area elevation increment
#' @param fill            Logical. Fill remaining depressions after breaching
#' @param stream_threshold Integer. Flow accumulation threshold in cells
#'
#' @return Invisibly returns hillshade_path
run_whitebox_group <- function(
    grp,
    input_dem,
    breached_path,
    pointer_path,
    accum_path,
    streams_path,
    hillshade_path,
    max_dist,
    flat_increment,
    fill,
    stream_threshold
) {

  # -- Step 1: Breach depressions ------------------------------------------
  if (!cache_exists(breached_path)) {

    cw_inform(glue::glue("Group '{grp}': breaching depressions..."))

    # wbt_breach_depressions_least_cost() is preferred over wbt_fill_depressions
    # as it modifies the DEM minimally, preserving natural topography while
    # ensuring all cells drain to an outlet.
    # dist = 10 cells (300 m at 30 m resolution) is conservative — large
    # values create artificial long drainage paths that alter catchment shape.
    # fill = TRUE fills any remaining depressions after breaching. This is
    # necessary for flat northern landscapes with thermokarst lakes and
    # other genuine closed depressions that breaching alone cannot resolve.
    # Without fill = TRUE, interior pit cells cause incorrect flow accumulation.
    whitebox::wbt_breach_depressions_least_cost(
      dem      = normalizePath(input_dem,     mustWork = TRUE),
      output   = normalizePath(breached_path, mustWork = FALSE),
      dist     = max_dist,
      fill     = fill
    )

    if (!cache_exists(breached_path)) {
      cw_abort(glue::glue(
        "Group '{grp}': wbt_breach_depressions_least_cost() did not produce ",
        "output at {breached_path}. Check WhiteboxTools log above."
      ))
    }

    cw_inform(glue::glue("Group '{grp}': dem_breached.tif written."))

  } else {
    cw_inform(glue::glue("Group '{grp}': dem_breached.tif found in cache, skipping."))
  }

  # -- Step 2: D8 flow pointer ---------------------------------------------
  if (!cache_exists(pointer_path)) {

    cw_inform(glue::glue("Group '{grp}': computing D8 flow pointer..."))

    whitebox::wbt_d8_pointer(
      dem    = normalizePath(breached_path, mustWork = TRUE),
      output = normalizePath(pointer_path,  mustWork = FALSE)
    )

    if (!cache_exists(pointer_path)) {
      cw_abort(glue::glue(
        "Group '{grp}': wbt_d8_pointer() did not produce output at ",
        "{pointer_path}. Check WhiteboxTools log above."
      ))
    }

    cw_inform(glue::glue("Group '{grp}': flow_pointer.tif written."))

  } else {
    cw_inform(glue::glue("Group '{grp}': flow_pointer.tif found in cache, skipping."))
  }

  # -- Step 3: D8 flow accumulation ----------------------------------------
  if (!cache_exists(accum_path)) {

    cw_inform(glue::glue("Group '{grp}': computing D8 flow accumulation..."))

    # out_type = "cells" returns raw upstream cell counts, which is more
    # useful for pour point snapping and catchment delineation.
    # pntr = TRUE tells WhiteboxTools the input is a flow pointer raster
    # rather than a DEM — without this the accumulation will be nonsensical.
    whitebox::wbt_d8_flow_accumulation(
      input    = normalizePath(pointer_path, mustWork = TRUE),
      output   = normalizePath(accum_path,   mustWork = FALSE),
      out_type = "cells",
      pntr     = TRUE
    )

    if (!cache_exists(accum_path)) {
      cw_abort(glue::glue(
        "Group '{grp}': wbt_d8_flow_accumulation() did not produce output ",
        "at {accum_path}. Check WhiteboxTools log above."
      ))
    }

    cw_inform(glue::glue("Group '{grp}': flow_accum.tif written."))

  } else {
    cw_inform(glue::glue("Group '{grp}': flow_accum.tif found in cache, skipping."))
  }

  # -- Step 4: Extract streams ---------------------------------------------
  if (!cache_exists(streams_path)) {

    cw_inform(glue::glue(
      "Group '{grp}': extracting streams ",
      "(threshold = {stream_threshold} cells)..."
    ))

    # wbt_extract_streams() produces a binary raster where cells with flow
    # accumulation >= threshold are classified as stream cells (1), others
    # are NoData. This is used for pour point snapping in 05_delineate_sites.R.
    whitebox::wbt_extract_streams(
      flow_accum = normalizePath(accum_path,   mustWork = TRUE),
      output     = normalizePath(streams_path, mustWork = FALSE),
      threshold  = stream_threshold
    )

    if (!cache_exists(streams_path)) {
      cw_abort(glue::glue(
        "Group '{grp}': wbt_extract_streams() did not produce output at ",
        "{streams_path}. Check WhiteboxTools log above."
      ))
    }

    cw_inform(glue::glue("Group '{grp}': streams.tif written."))

  } else {
    cw_inform(glue::glue("Group '{grp}': streams.tif found in cache, skipping."))
  }

  # -- Step 5: Hillshade ---------------------------------------------------
  if (!cache_exists(hillshade_path)) {

    cw_inform(glue::glue("Group '{grp}': computing hillshade..."))

    # wbt_hillshade() is computed from the breached DEM rather than the raw
    # DEM so that the shading reflects the hydrologically conditioned surface.
    # Default azimuth (315°) and altitude (30°) produce natural-looking results.
    whitebox::wbt_hillshade(
      dem    = normalizePath(breached_path,   mustWork = TRUE),
      output = normalizePath(hillshade_path,  mustWork = FALSE)
    )

    if (!cache_exists(hillshade_path)) {
      cw_abort(glue::glue(
        "Group '{grp}': wbt_hillshade() did not produce output at ",
        "{hillshade_path}. Check WhiteboxTools log above."
      ))
    }

    cw_inform(glue::glue("Group '{grp}': hillshade.tif written."))

  } else {
    cw_inform(glue::glue("Group '{grp}': hillshade.tif found in cache, skipping."))
  }

  invisible(hillshade_path)
}

#' Verify Whitebox outputs for all groups
#'
#' Checks that all three Whitebox outputs exist and reports basic raster
#' properties. Useful for confirming the conditioning steps completed
#' correctly before running site-level delineation.
#'
#' @param group_manifest sf tibble from build_group_manifest()
#' @return A tibble with one row per group summarising output properties
verify_whitebox_outputs <- function(group_manifest) {

  purrr::map(seq_len(nrow(group_manifest)), function(i) {

    grp       <- group_manifest$group_id[i]
    grp_cache <- group_manifest$cache_dir[i]

    files <- c(
      dem_breached = fs::path(grp_cache, "dem_breached.tif"),
      flow_pointer = fs::path(grp_cache, "flow_pointer.tif"),
      flow_accum   = fs::path(grp_cache, "flow_accum.tif"),
      streams      = fs::path(grp_cache, "streams.tif"),
      hillshade    = fs::path(grp_cache, "hillshade.tif")
    )

    purrr::imap(files, function(path, name) {
      if (!cache_exists(path)) {
        tibble::tibble(
          group_id  = grp,
          file      = name,
          exists    = FALSE,
          size_mb   = NA_real_,
          min_val   = NA_real_,
          max_val   = NA_real_
        )
      } else {
        r <- terra::rast(path)
        v <- terra::values(r, na.rm = TRUE)
        tibble::tibble(
          group_id  = grp,
          file      = name,
          exists    = TRUE,
          size_mb   = round(fs::file_size(path) / 1e6, 1),
          min_val   = round(min(v), 2),
          max_val   = round(max(v), 2)
        )
      }
    }) |>
      dplyr::bind_rows()

  }) |>
    dplyr::bind_rows()
}
