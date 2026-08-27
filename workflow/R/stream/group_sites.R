# 01_group_sites.R
# ---------------------------------------------------------------------------
# Groups sites into processing units, builds group AOIs, and sets up the
# cache and output directory structure.
#
# Inputs:
#   sites          : Validated sites tibble from validate_sites()
#   output_dir     : Root directory for per-site output folders
#   cache_dir      : Root directory for group-level cached rasters
#   default_buffer_m : Global AOI buffer in metres (overridden per group via
#                    aoi_buffer_m column in sites CSV)
#
# Outputs:
#   A tibble (group manifest) with one row per group containing:
#     group_id, burn_streams, aoi (sfc geometry), cache_dir, n_sites
#   Directory structure created on disk:
#     cache/<group_id>/       — for group-level raster derivatives
#     output/<site_id>/       — for per-site outputs
#
# Dependencies: sf, dplyr, fs, cli (via utils.R)
# ---------------------------------------------------------------------------

#' Build the group manifest and set up all directories
#'
#' This is the main entry point for this module. It validates that each group
#' has an internally consistent set of group-level parameters, constructs the
#' AOI for each group using HydroBasins level 6 polygons, creates all required
#' directories, and returns a group manifest tibble used by downstream modules.
#'
#' @param sites            Validated sites tibble from validate_sites()
#' @param output_dir       Root directory where per-site output folders will
#'                         be created
#' @param cache_dir        Root directory where group-level raster cache
#'                         folders will be created
#' @param hydrobasins_dir  Path to HydroBasins root directory containing
#'                         'north_america' and 'arctic' subfolders
#' @param default_buffer_m Buffer in metres added on top of the HydroBasins
#'                         polygon(s). Used for any group that does not
#'                         specify aoi_buffer_m in the CSV. Default 1000 m.
#' @param hybas_level      HydroBasins level to use for AOI definition.
#'                         Default 6.
#'
#' @return A tibble with one row per group and columns:
#'   group_id, burn_streams, buffer_m, aoi, cache_dir, n_sites
build_group_manifest <- function(
  sites,
  output_dir,
  cache_dir,
  hydrobasins_dir,
  default_buffer_m = 1000,
  hybas_level = 6
) {
  cw_inform("Building group manifest and setting up directories...")

  # Convert sites to sf for spatial operations
  sites_sf <- sites_to_sf(sites)

  # Build one row per group
  groups <- get_groups(sites)

  group_manifest <- purrr::map(seq_len(nrow(groups)), function(i) {
    grp <- groups$group_id[i]

    # Resolve buffer: use group override if present, else global default
    grp_buffer <- groups$aoi_buffer_m[i]
    buffer_m <- if (!is.na(grp_buffer)) grp_buffer else default_buffer_m

    cw_inform(glue::glue(
      "Group '{grp}': buffer = {buffer_m} m",
      "{if (!is.na(grp_buffer)) ' (group override)' else ' (global default)'}"
    ))

    # Build AOI polygon for this group in EPSG:3979 using HydroBasins
    aoi <- build_group_aoi(
      sites_sf = sites_sf,
      group = grp,
      hydrobasins_dir = hydrobasins_dir,
      default_buffer_m = default_buffer_m,
      hybas_level = hybas_level
    )

    # Count sites in this group
    n_sites <- sum(sites$group_id == grp)

    # Resolve directory paths for this group
    grp_cache <- group_cache_dir(cache_dir, grp)

    tibble::tibble(
      group_id = grp,
      burn_streams = groups$burn_streams[i],
      buffer_m = buffer_m,
      aoi = aoi, # sfc geometry column
      cache_dir = grp_cache,
      n_sites = n_sites
    )
  }) |>
    dplyr::bind_rows() |>
    sf::st_as_sf(sf_column_name = "aoi", crs = 3979)

  # Report group summary
  cw_inform(glue::glue(
    "Found {nrow(group_manifest)} group(s) covering {nrow(sites)} site(s):"
  ))
  purrr::walk(seq_len(nrow(group_manifest)), function(i) {
    grp <- group_manifest[i, ]
    cw_inform(glue::glue(
      "  {grp$group_id}: {grp$n_sites} site(s), ",
      "burn_streams = {grp$burn_streams}, ",
      "buffer = {grp$buffer_m} m"
    ))
  })

  # Create group cache directories
  cw_inform("Creating cache directories...")
  purrr::walk(group_manifest$cache_dir, ensure_dir)

  # Create per-site output directories
  cw_inform("Creating output directories...")
  purrr::walk(sites$site_id, function(sid) {
    ensure_dir(site_output_dir(output_dir, sid))
  })

  cw_inform("Directory setup complete.")

  group_manifest
}

#' Write the group manifest AOIs to disk as a GeoPackage for inspection
#'
#' Useful for visually verifying that group AOIs cover all sites correctly
#' before running the full workflow. Writes to cache_dir/group_aois.gpkg.
#'
#' @param group_manifest Group manifest tibble from build_group_manifest()
#' @param cache_dir      Root cache directory
write_group_aois <- function(group_manifest, cache_dir) {
  out_path <- fs::path(cache_dir, "group_aois.gpkg")

  group_manifest |>
    dplyr::select(group_id, burn_streams, buffer_m, n_sites) |>
    sf::st_write(out_path, delete_dsn = TRUE, quiet = TRUE)

  cw_inform(glue::glue("Group AOIs written to: {out_path}"))

  invisible(out_path)
}

#' Check which groups still need processing
#'
#' Compares expected cache outputs against what exists on disk. Returns a
#' filtered version of group_manifest containing only groups where one or more
#' expected cache files are missing. Used by downstream modules to skip groups
#' that have already been fully processed.
#'
#' Expected cache files checked here are the final Whitebox outputs (pointer
#' and flow accumulation), since these are the last steps in the group
#' processing chain. If both exist, all upstream steps are assumed complete.
#'
#' @param group_manifest Group manifest tibble from build_group_manifest()
#' @return A filtered group manifest containing only groups needing processing
groups_needing_processing <- function(group_manifest) {
  needs_processing <- purrr::map_lgl(
    seq_len(nrow(group_manifest)),
    function(i) {
      grp_cache <- group_manifest$cache_dir[i]

      # These are the terminal outputs of the group processing chain
      pointer_path <- fs::path(grp_cache, "flow_pointer.tif")
      accum_path <- fs::path(grp_cache, "flow_accum.tif")

      # Group needs processing if either terminal file is missing
      !cache_exists(pointer_path) || !cache_exists(accum_path)
    }
  )

  n_skip <- sum(!needs_processing)
  n_run <- sum(needs_processing)

  if (n_skip > 0) {
    cw_inform(glue::glue(
      "Skipping {n_skip} already-processed group(s): ",
      "{paste(group_manifest$group_id[!needs_processing], collapse = ', ')}"
    ))
  }
  if (n_run > 0) {
    cw_inform(glue::glue(
      "Processing {n_run} group(s): ",
      "{paste(group_manifest$group_id[needs_processing], collapse = ', ')}"
    ))
  }

  group_manifest[needs_processing, ]
}
