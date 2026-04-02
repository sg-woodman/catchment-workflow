# utils.R
# ---------------------------------------------------------------------------
# Shared utility functions for the catchment delineation workflow.
# Sourced by delineate_catchments.R and all module scripts.
#
# Dependencies: sf, terra, dplyr, tidyr, fs, cli, readr
# ---------------------------------------------------------------------------

# -- Logging -----------------------------------------------------------------

#' Emit an informational message with a timestamp
#' @param msg Character string to print
cw_inform <- function(msg) {
  cli::cli_inform(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg))
}

#' Emit a warning message with a timestamp
#' @param msg Character string to print
cw_warn <- function(msg) {
  cli::cli_warn(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg))
}

#' Emit an error and stop execution
#' @param msg Character string to print
cw_abort <- function(msg) {
  cli::cli_abort(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg))
}

# -- Input validation --------------------------------------------------------

#' Validate and read the sites CSV file
#'
#' Reads a CSV file (comment lines beginning with "#" are ignored), then
#' runs all validation checks via validate_sites_impl(). The optional
#' aoi_buffer_m column overrides the global buffer default at the group level.
#' Returns a validated tibble.
#'
#' @param sites_csv Path to the sites CSV file
#' @return A validated tibble with columns: site_id, site_name, lon, lat,
#'   group_id, burn_streams, aoi_buffer_m (NA where not specified)
validate_sites <- function(sites_csv) {
  # Read CSV, skipping comment lines beginning with "#".
  # aoi_buffer_m is optional — col_double() yields NA for blank cells.
  sites <- readr::read_csv(
    sites_csv,
    comment = "#",
    col_types = readr::cols(
      site_id = readr::col_character(),
      site_name = readr::col_character(),
      lon = readr::col_double(),
      lat = readr::col_double(),
      group_id = readr::col_character(),
      burn_streams = readr::col_logical(),
      aoi_buffer_m = readr::col_double() # NA if blank
    ),
    show_col_types = FALSE
  )

  # If aoi_buffer_m column is absent entirely, add it as all-NA
  if (!"aoi_buffer_m" %in% names(sites)) {
    sites <- dplyr::mutate(sites, aoi_buffer_m = NA_real_)
  }

  validate_sites_impl(sites)
}

#' Validate a sites tibble provided directly (as an alternative to a CSV)
#'
#' Accepts a tibble instead of a CSV path, coerces column types, and runs
#' the same validation checks as validate_sites(). Use this when building
#' your site list programmatically in R.
#'
#' @param sites A tibble with columns: site_id, site_name, lon, lat, group_id,
#'   burn_streams, and optionally aoi_buffer_m
#' @return A validated tibble with aoi_buffer_m added as NA if absent
validate_sites_tibble <- function(sites) {
  # Ensure aoi_buffer_m exists, adding it as all-NA if absent
  if (!"aoi_buffer_m" %in% names(sites)) {
    sites <- dplyr::mutate(sites, aoi_buffer_m = NA_real_)
  }

  # Coerce column types to match what validate_sites() produces from CSV
  sites <- sites |>
    dplyr::mutate(
      site_id = as.character(site_id),
      site_name = as.character(site_name),
      lon = as.double(lon),
      lat = as.double(lat),
      group_id = as.character(group_id),
      burn_streams = as.logical(burn_streams),
      aoi_buffer_m = as.double(aoi_buffer_m)
    )

  validate_sites_impl(sites)
}

#' Internal validation logic shared by validate_sites() and
#' validate_sites_tibble()
#'
#' Operates purely on a tibble — no file I/O. Not intended to be called
#' directly; use validate_sites() or validate_sites_tibble() instead.
#'
#' @param sites A tibble with expected columns already present and typed
#' @return The validated tibble, unchanged if all checks pass
validate_sites_impl <- function(sites) {
  # Check required columns (aoi_buffer_m is optional, excluded here)
  required_cols <- c(
    "site_id",
    "site_name",
    "lon",
    "lat",
    "group_id",
    "burn_streams"
  )
  missing_cols <- setdiff(required_cols, names(sites))
  if (length(missing_cols) > 0) {
    cw_abort(paste(
      "sites input is missing required columns:",
      paste(missing_cols, collapse = ", ")
    ))
  }

  # Check for missing values in any required column
  has_na <- sites |>
    dplyr::filter(dplyr::if_any(dplyr::all_of(required_cols), is.na)) |>
    dplyr::pull(site_id)
  if (length(has_na) > 0) {
    cw_abort(paste(
      "Missing values found in required columns for sites:",
      paste(has_na, collapse = ", ")
    ))
  }

  # Check for duplicate site_ids
  dupes <- sites |>
    dplyr::filter(duplicated(site_id)) |>
    dplyr::pull(site_id)
  if (length(dupes) > 0) {
    cw_abort(paste(
      "Duplicate site_id values found:",
      paste(dupes, collapse = ", ")
    ))
  }

  # Check site_id and group_id contain only alphanumeric characters,
  # underscores, or hyphens (no spaces or special characters)
  invalid_site_ids <- sites |>
    dplyr::filter(!grepl("^[A-Za-z0-9_\\-]+$", site_id)) |>
    dplyr::pull(site_id)
  if (length(invalid_site_ids) > 0) {
    cw_abort(paste(
      "site_id values must contain only letters, numbers, underscores, or hyphens:",
      paste(invalid_site_ids, collapse = ", ")
    ))
  }

  invalid_group_ids <- sites |>
    dplyr::filter(!grepl("^[A-Za-z0-9_\\-]+$", group_id)) |>
    dplyr::pull(group_id)
  if (length(invalid_group_ids) > 0) {
    cw_abort(paste(
      "group_id values must contain only letters, numbers, underscores, or hyphens:",
      paste(invalid_group_ids, collapse = ", ")
    ))
  }

  # Check that burn_streams is consistent within each group (group-level flag)
  burn_consistency <- sites |>
    dplyr::group_by(group_id) |>
    dplyr::summarise(
      n_burn_values = dplyr::n_distinct(burn_streams),
      .groups = "drop"
    ) |>
    dplyr::filter(n_burn_values > 1)
  if (nrow(burn_consistency) > 0) {
    cw_abort(paste(
      "burn_streams must be consistent within each group_id. Conflicting groups:",
      paste(burn_consistency$group_id, collapse = ", ")
    ))
  }

  # Check that aoi_buffer_m is consistent within each group where provided
  buffer_consistency <- sites |>
    dplyr::filter(!is.na(aoi_buffer_m)) |>
    dplyr::group_by(group_id) |>
    dplyr::summarise(
      n_buffer_values = dplyr::n_distinct(aoi_buffer_m),
      .groups = "drop"
    ) |>
    dplyr::filter(n_buffer_values > 1)
  if (nrow(buffer_consistency) > 0) {
    cw_abort(paste(
      "aoi_buffer_m must be consistent within each group_id. Conflicting groups:",
      paste(buffer_consistency$group_id, collapse = ", ")
    ))
  }

  # Check that coordinates are plausible decimal degrees.
  # Checks for: lon outside -180:0 (not Canadian), lat outside 40:90,
  # and the common mistake of swapped lat/lon.
  bad_lon <- sites |>
    dplyr::filter(lon < -180 | lon > 0) |>
    dplyr::pull(site_id)
  if (length(bad_lon) > 0) {
    cw_abort(paste(
      "Longitude values out of range (expected negative decimal degrees for Canada):",
      paste(bad_lon, collapse = ", ")
    ))
  }

  bad_lat <- sites |>
    dplyr::filter(lat < 40 | lat > 90) |>
    dplyr::pull(site_id)
  if (length(bad_lat) > 0) {
    cw_abort(paste(
      "Latitude values out of range (expected 40-90 for Canada):",
      paste(bad_lat, collapse = ", ")
    ))
  }

  # Warn if lat and lon appear swapped (lon more positive than lat)
  swapped <- sites |>
    dplyr::filter(lon > lat) |>
    dplyr::pull(site_id)
  if (length(swapped) > 0) {
    cw_warn(paste(
      "Possible swapped lat/lon — longitude is greater than latitude for:",
      paste(swapped, collapse = ", "),
      "-- verify coordinates are (lon, lat) in WGS84 decimal degrees."
    ))
  }

  cw_inform(paste(
    "Validated",
    nrow(sites),
    "sites across",
    dplyr::n_distinct(sites$group_id),
    "groups."
  ))

  sites
}

# -- Group helpers -----------------------------------------------------------

#' Extract the unique group-level attributes from the validated sites tibble
#'
#' Returns one row per group with group_id, burn_streams, and aoi_buffer_m.
#'
#' @param sites Validated sites tibble from validate_sites() or
#'   validate_sites_tibble()
#' @return A tibble with columns: group_id, burn_streams, aoi_buffer_m
get_groups <- function(sites) {
  sites |>
    dplyr::distinct(group_id, burn_streams, aoi_buffer_m)
}

#' Convert sites tibble to an sf point object (WGS84)
#'
#' @param sites Validated sites tibble
#' @return An sf object with point geometry in EPSG:4326
sites_to_sf <- function(sites) {
  sf::st_as_sf(sites, coords = c("lon", "lat"), crs = 4326)
}

# -- Spatial helpers ---------------------------------------------------------

#' Resolve the HydroBasins region tag for a set of points
#'
#' HydroBasins splits Canada between the North America (na) and Arctic (ar)
#' datasets. This is determined spatially by checking which region's level 1
#' file contains the points, rather than using a latitude threshold which can
#' misclassify sites near the boundary.
#'
#' @param group_pts sf point object for the group, in any CRS
#' @param hydrobasins_dir Character. Path to HydroBasins root directory
#' @return Character: "na" or "ar"
resolve_hydrobasins_region <- function(group_pts, hydrobasins_dir) {
  # Read both level-1 files, repair geometries, then determine which region
  # contains the points. Level 1 files are tiny so this is fast.
  # st_make_valid() is applied immediately after reading as HydroBasins
  # files frequently contain invalid geometries.
  na_lev1 <- sf::st_read(
    fs::path(hydrobasins_dir, "north_america", "hybas_na_lev01_v1c.shp"),
    quiet = TRUE
  ) |>
    sf::st_make_valid()

  ar_lev1 <- sf::st_read(
    fs::path(hydrobasins_dir, "arctic", "hybas_ar_lev01_v1c.shp"),
    quiet = TRUE
  ) |>
    sf::st_make_valid()

  pts_wgs84 <- sf::st_transform(group_pts, 4326)

  # Count how many points fall in each region
  n_na <- sum(
    lengths(sf::st_intersects(
      pts_wgs84,
      sf::st_make_valid(sf::st_union(na_lev1))
    )) >
      0
  )
  n_ar <- sum(
    lengths(sf::st_intersects(
      pts_wgs84,
      sf::st_make_valid(sf::st_union(ar_lev1))
    )) >
      0
  )

  # Use the region containing the majority of points. Ties go to na.
  if (n_ar > n_na) "ar" else "na"
}

#' Build a HydroBasins-derived AOI for a group of sites
#'
#' Finds all HydroBasins level 6 polygons containing any site in the group,
#' unions them into a single polygon, then adds a small buffer. This produces
#' a hydrologically meaningful AOI that ensures the full upstream area is
#' captured for each site. Falls back to a buffered bounding box if no
#' HydroBasins polygons are found.
#'
#' The na/ar region is resolved automatically by spatial intersection with
#' the level 1 boundaries.
#'
#' @param sites_sf         sf point object from sites_to_sf()
#' @param group            Character. The group_id to build an AOI for
#' @param hydrobasins_dir  Character. Path to HydroBasins root directory
#' @param default_buffer_m Numeric. Buffer in metres added on top of the
#'   HydroBasins polygon(s). Group-level aoi_buffer_m overrides this.
#' @param hybas_level      Integer. HydroBasins level to use. Default 6.
#' @return An sfc polygon in EPSG:3979
build_group_aoi <- function(
  sites_sf,
  group,
  hydrobasins_dir,
  default_buffer_m = 1000,
  hybas_level = 6
) {
  group_pts <- sites_sf |>
    dplyr::filter(group_id == group) |>
    sf::st_transform(3979)

  # Resolve buffer: group override if present, else global default
  group_buffer <- group_pts |>
    sf::st_drop_geometry() |>
    dplyr::pull(aoi_buffer_m) |>
    unique() |>
    stats::na.omit()

  buffer_m <- if (length(group_buffer) == 1) group_buffer else default_buffer_m

  # Determine which HydroBasins region file to use (na or ar)
  region <- resolve_hydrobasins_region(group_pts, hydrobasins_dir)

  # Build path to the appropriate level file
  level_str <- formatC(hybas_level, width = 2, flag = "0")
  hybas_path <- fs::path(
    hydrobasins_dir,
    if (region == "ar") "arctic" else "north_america",
    glue::glue("hybas_{region}_lev{level_str}_v1c.shp")
  )

  if (!fs::file_exists(hybas_path)) {
    cw_warn(glue::glue(
      "Group '{group}': HydroBasins file not found at {hybas_path}. ",
      "Falling back to buffered bounding box."
    ))
    return(
      group_pts |>
        sf::st_bbox() |>
        sf::st_as_sfc() |>
        sf::st_buffer(buffer_m)
    )
  }

  # Read HydroBasins, repair any invalid geometries, and find polygons
  # containing any group site. Geometry issues (duplicate vertices,
  # self-intersections) are common in HydroBasins data and must be fixed
  # before st_union() is called.
  hybas <- sf::st_read(hybas_path, quiet = TRUE) |>
    sf::st_transform(3979) |>
    sf::st_make_valid()

  pts_3979 <- sf::st_transform(group_pts, 3979)
  intersecting <- hybas[
    lengths(sf::st_intersects(hybas, pts_3979)) > 0,
  ]

  if (nrow(intersecting) == 0) {
    cw_warn(glue::glue(
      "Group '{group}': no HydroBasins level {hybas_level} polygons found ",
      "containing group sites. Falling back to buffered bounding box."
    ))
    return(
      group_pts |>
        sf::st_bbox() |>
        sf::st_as_sfc() |>
        sf::st_buffer(buffer_m)
    )
  }

  cw_inform(glue::glue(
    "Group '{group}': using {nrow(intersecting)} HydroBasins level ",
    "{hybas_level} polygon(s) [{region}] + {buffer_m} m buffer."
  ))

  # Union all intersecting basins, re-validate after union (union can
  # introduce new geometry issues), then add buffer. A zero-width buffer
  # via st_buffer(0) is a common alternative validity fix but
  # st_make_valid() is more robust for HydroBasins data.
  aoi <- intersecting |>
    sf::st_union() |>
    sf::st_make_valid() |>
    sf::st_buffer(buffer_m) |>
    sf::st_make_valid()

  aoi
}

# -- File system helpers -----------------------------------------------------

#' Ensure a directory exists, creating it if necessary
#' @param path Character. Directory path to create
ensure_dir <- function(path) {
  fs::dir_create(path, recurse = TRUE)
}

#' Return the cache directory path for a given group
#' @param cache_dir Character. Root cache directory
#' @param group_id  Character. Group identifier
group_cache_dir <- function(cache_dir, group_id) {
  fs::path(cache_dir, group_id)
}

#' Return the output directory path for a given site
#' @param output_dir Character. Root output directory
#' @param site_id    Character. Site identifier
site_output_dir <- function(output_dir, site_id) {
  fs::path(output_dir, site_id)
}

#' Check whether a cached file exists and is non-empty
#'
#' Used throughout the workflow to skip steps that have already been completed.
#'
#' @param path Character. File path to check
#' @return Logical
cache_exists <- function(path) {
  fs::file_exists(path) && fs::file_size(path) > 0
}

# -- Package checks ----------------------------------------------------------

#' Check that all required packages are installed
#'
#' Called at the top of delineate_catchments.R. Lists all packages needed
#' across all modules so the user gets a single informative error up front
#' rather than a failure mid-run.
check_packages <- function() {
  required <- c(
    "sf",
    "terra",
    "whitebox",
    "dplyr",
    "tidyr",
    "purrr",
    "readr",
    "tibble",
    "fs",
    "cli",
    "httr2",
    "glue"
  )

  missing <- required[
    !vapply(required, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing) > 0) {
    cw_abort(paste(
      "The following required packages are not installed:",
      paste(missing, collapse = ", "),
      "\nInstall with: install.packages(c(",
      paste0('"', missing, '"', collapse = ", "),
      "))"
    ))
  }

  # Check WhiteboxTools binary is available
  if (!whitebox::wbt_init()) {
    cw_abort(paste(
      "WhiteboxTools binary not found.",
      "Install with: whitebox::install_whitebox()"
    ))
  }

  invisible(TRUE)
}
