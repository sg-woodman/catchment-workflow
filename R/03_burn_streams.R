# 03_burn_streams.R
# ---------------------------------------------------------------------------
# Identifies which NHN GDB files intersect each group AOI, extracts the
# relevant flowline and waterbody layers, and optionally burns the stream
# network into the group DEM using wbt_fill_burn().
#
# Stream burning is controlled by the burn_streams flag in the group manifest
# (set per group in sites CSV). Groups with burn_streams = FALSE skip the
# burn step but NHN layers are still extracted and cached for use in
# 06_extract_lakes.R.
#
# NHN file structure assumed on disk:
#   <nhn_dir>/
#     nhn_rhn_<sheet>_gdb_en/
#       <anything>.gdb          # GDB filename varies — matched by pattern
#
# NHN layers used:
#   NHN_HN_PrimaryDirectedNLFlow_1  : Flowlines (used for stream burning)
#   NHN_HD_WATERBODY_2              : Lakes/ponds/watercourses (for later use)
#
# Inputs:
#   group_manifest : sf tibble from build_group_manifest()
#   nhn_dir        : Path to root NHN directory containing sheet subfolders
#   nhn_index      : Path to NHN index shapefile (.shp) mapping sheet codes
#                    to geometries
#
# Outputs (per group, written to cache/<group_id>/):
#   flowlines.gpkg     : Merged NHN flowlines clipped to group AOI (EPSG:3979)
#   waterbodies.gpkg   : Merged NHN waterbodies clipped to group AOI (EPSG:3979)
#   streams_raster.tif : Rasterised flowlines aligned to dem.tif (burn groups only)
#   dem_burned.tif     : Stream-burned DEM (burn groups only)
#
# Dependencies: sf, terra, whitebox, dplyr, purrr, fs, cli (via utils.R)
# ---------------------------------------------------------------------------

# -- Layer name constants ----------------------------------------------------

# NHN GDB layer names — consistent across all GDB versions
NHN_FLOWLINE_LAYER  <- "NHN_HN_PrimaryDirectedNLFlow_1"
NHN_WATERBODY_LAYER <- "NHN_HD_WATERBODY_2"

# -- Main entry point --------------------------------------------------------

#' Prepare NHN layers and optionally burn streams for all groups
#'
#' Iterates over all groups. For each group, identifies intersecting NHN
#' sheets, reads and merges flowlines and waterbodies, clips to AOI, and
#' writes to cache. If burn_streams = TRUE for the group, also rasterises
#' flowlines and runs wbt_fill_burn() on the group DEM.
#'
#' Groups where all expected outputs already exist in cache are skipped.
#'
#' @param group_manifest sf tibble from build_group_manifest()
#' @param nhn_dir        Character. Path to root NHN directory
#' @param nhn_index      Character. Path to NHN index .shp file
#'
#' @return The input group_manifest invisibly. Called for side effects.
prepare_nhn_layers <- function(
    group_manifest,
    nhn_dir,
    nhn_index
) {

  # Validate inputs exist before starting
  if (!fs::dir_exists(nhn_dir)) {
    cw_abort(glue::glue("NHN directory not found: {nhn_dir}"))
  }
  if (!fs::file_exists(nhn_index)) {
    cw_abort(glue::glue("NHN index shapefile not found: {nhn_index}"))
  }

  # Read NHN index once — used for all groups
  cw_inform("Reading NHN index...")
  nhn_idx <- sf::st_read(nhn_index, quiet = TRUE) |>
    sf::st_transform(3979)

  purrr::walk(seq_len(nrow(group_manifest)), function(i) {

    grp          <- group_manifest$group_id[i]
    grp_cache    <- group_manifest$cache_dir[i]
    aoi          <- group_manifest$aoi[i] |> sf::st_as_sf()
    burn_streams <- group_manifest$burn_streams[i]

    # Define expected output paths for this group
    flowlines_path   <- fs::path(grp_cache, "flowlines.gpkg")
    waterbodies_path <- fs::path(grp_cache, "waterbodies.gpkg")
    dem_burned_path  <- fs::path(grp_cache, "dem_burned.tif")
    dem_path         <- fs::path(grp_cache, "dem.tif")

    # Determine which outputs are expected for this group
    # (burned DEM only needed if burn_streams = TRUE)
    vector_outputs_exist <- cache_exists(flowlines_path) &&
      cache_exists(waterbodies_path)
    burn_outputs_exist   <- cache_exists(dem_burned_path)

    # Skip entirely if all required outputs exist
    if (vector_outputs_exist && (!burn_streams || burn_outputs_exist)) {
      cw_inform(glue::glue("Group '{grp}': NHN outputs found in cache, skipping."))
      return(invisible(NULL))
    }

    cw_inform(glue::glue(
      "Group '{grp}': preparing NHN layers",
      "{if (burn_streams) ' and burning streams' else ''}..."
    ))

    # -- Step 1: Find intersecting NHN sheets --------------------------------

    sheets <- find_nhn_sheets(aoi, nhn_idx, nhn_dir, grp)

    if (length(sheets) == 0) {
      cw_warn(glue::glue(
        "Group '{grp}': no NHN sheets found intersecting AOI. ",
        "Flowlines and waterbodies will be empty. ",
        "burn_streams will be skipped for this group."
      ))
      # Write empty outputs so cache checks pass and downstream modules
      # can handle missing NHN gracefully
      write_empty_nhn_outputs(aoi, flowlines_path, waterbodies_path)
      return(invisible(NULL))
    }

    # -- Step 2: Read and merge NHN layers -----------------------------------

    if (!vector_outputs_exist) {

      flowlines   <- read_merge_nhn_layer(sheets, NHN_FLOWLINE_LAYER,  aoi, grp)
      waterbodies <- read_merge_nhn_layer(sheets, NHN_WATERBODY_LAYER, aoi, grp)

      # Write vector outputs to cache
      if (!is.null(flowlines)) {
        sf::st_write(flowlines, flowlines_path, delete_dsn = TRUE, quiet = TRUE)
        cw_inform(glue::glue(
          "Group '{grp}': flowlines written ({nrow(flowlines)} features)."
        ))
      } else {
        write_empty_nhn_outputs(aoi, flowlines_path, waterbodies_path)
        cw_warn(glue::glue(
          "Group '{grp}': no flowlines found in intersecting NHN sheets."
        ))
      }

      if (!is.null(waterbodies)) {
        sf::st_write(waterbodies, waterbodies_path, delete_dsn = TRUE, quiet = TRUE)
        cw_inform(glue::glue(
          "Group '{grp}': waterbodies written ({nrow(waterbodies)} features)."
        ))
      } else {
        write_empty_nhn_output_layer(aoi, waterbodies_path, type = "polygon")
        cw_warn(glue::glue(
          "Group '{grp}': no waterbodies found in intersecting NHN sheets."
        ))
      }
    }

    # -- Step 3: Burn streams into DEM (if requested) ------------------------

    if (burn_streams && !burn_outputs_exist) {

      # Verify DEM exists before attempting burn
      if (!cache_exists(dem_path)) {
        cw_abort(glue::glue(
          "Group '{grp}': dem.tif not found at {dem_path}. ",
          "Run prepare_dem() before prepare_nhn_layers()."
        ))
      }

      # Reload flowlines from cache (may have just been written, or existed)
      flowlines_cached <- sf::st_read(flowlines_path, quiet = TRUE)

      if (nrow(flowlines_cached) == 0) {
        cw_warn(glue::glue(
          "Group '{grp}': burn_streams = TRUE but no flowlines available. ",
          "Skipping burn — dem.tif will be used as-is for Whitebox steps."
        ))
        return(invisible(NULL))
      }

      burn_streams_into_dem(
        flowlines       = flowlines_cached,
        dem_path        = dem_path,
        dem_burned_path = dem_burned_path,
        group_id        = grp
      )
    }

    invisible(NULL)
  })

  invisible(group_manifest)
}

# -- NHN sheet discovery -----------------------------------------------------

#' Find NHN GDB paths whose sheets intersect the group AOI
#'
#' Uses the NHN index shapefile to identify which NTS sheet codes intersect
#' the group AOI, then locates the corresponding GDB files on disk. Handles
#' the inconsistent GDB naming convention by searching each subfolder for
#' any .gdb file.
#'
#' @param aoi     sf polygon. Group AOI in EPSG:3979
#' @param nhn_idx sf object. NHN index read from the index shapefile
#' @param nhn_dir Character. Path to root NHN directory
#' @param group_id Character. Group identifier (for log messages)
#'
#' @return Character vector of .gdb paths intersecting the AOI.
#'   Returns character(0) if none found.
find_nhn_sheets <- function(aoi, nhn_idx, nhn_dir, group_id) {

  # Find index features that intersect the AOI
  # Suppress warnings from non-matching geometries
  intersecting <- suppressWarnings(
    nhn_idx[sf::st_intersects(nhn_idx, aoi, sparse = FALSE)[, 1], ]
  )

  if (nrow(intersecting) == 0) {
    return(character(0))
  }

  # Extract sheet codes from the index.
  # The NHN index uses a "WSCSSDA" or "DATASETNAM" field for the sheet code —
  # detect which is present
  code_field <- intersect(
    c("WSCSSDA", "DATASETNAM", "NID", "WSCMDA"),
    names(intersecting)
  )[1]

  if (is.na(code_field)) {
    cw_warn(glue::glue(
      "Group '{group_id}': could not identify sheet code field in NHN index. ",
      "Available fields: {paste(names(intersecting), collapse = ', ')}"
    ))
    return(character(0))
  }

  sheet_codes <- intersecting[[code_field]] |>
    tolower() |>
    unique()

  cw_inform(glue::glue(
    "Group '{group_id}': found {length(sheet_codes)} intersecting NHN sheet(s): ",
    "{paste(sheet_codes, collapse = ', ')}"
  ))

  # Locate GDB files on disk for each sheet code.
  # The NHN index uses 4-character WSCSSDA codes (e.g. "10LC") but folder
  # names on disk use 7-character codes (e.g. "nhn_rhn_10lc000_gdb_en").
  # We match by prefix against all folders in nhn_dir, which handles all
  # suffix variants (000, 001, 002, etc.) robustly.

  all_subfolders      <- fs::dir_ls(nhn_dir, type = "directory")
  all_subfolder_names <- fs::path_file(all_subfolders)

  gdb_paths <- purrr::map(sheet_codes, function(code) {

    pattern <- paste0("^nhn_rhn_", tolower(code))
    matched <- all_subfolders[grepl(pattern, all_subfolder_names,
                                    ignore.case = TRUE)]

    if (length(matched) == 0) {
      cw_warn(glue::glue(
        "Group '{group_id}': NHN subfolder not found for sheet '{code}': ",
        "no folder matching 'nhn_rhn_{tolower(code)}*' in {nhn_dir}"
      ))
      return(NULL)
    }

    purrr::map_chr(matched, function(subfolder) {
      gdbs <- fs::dir_ls(subfolder, glob = "*.gdb", type = "directory")
      if (length(gdbs) == 0) {
        cw_warn(glue::glue(
          "Group '{group_id}': no .gdb found in {subfolder}"
        ))
        return(NA_character_)
      }
      gdbs[1]
    })

  }) |>
    unlist() |>
    stats::na.omit() |>
    as.character()

  gdb_paths <- unique(gdb_paths)

  cw_inform(glue::glue(
    "Group '{group_id}': located {length(gdb_paths)} GDB file(s) on disk."
  ))

  gdb_paths
}

# -- NHN layer reading -------------------------------------------------------

#' Read, merge, and clip an NHN layer from multiple GDB files
#'
#' Reads the specified layer from each GDB, reprojects to EPSG:3979, clips
#' to the group AOI, and merges into a single sf object.
#'
#' @param gdb_paths  Character vector of .gdb paths from find_nhn_sheets()
#' @param layer_name Character. NHN layer name to read
#' @param aoi        sf polygon. Group AOI in EPSG:3979
#' @param group_id   Character. Group identifier (for log messages)
#'
#' @return Merged sf object clipped to AOI, or NULL if no features found
read_merge_nhn_layer <- function(gdb_paths, layer_name, aoi, group_id) {

  layers_list <- purrr::map(gdb_paths, function(gdb) {

    lyr <- read_nhn_from_gdb(gdb, layer_name, group_id)

    if (is.null(lyr) || nrow(lyr) == 0) return(NULL)

    clipped <- tryCatch(
      sf::st_intersection(lyr, aoi),
      error = function(e) {
        cw_warn(glue::glue(
          "Group '{group_id}': error clipping '{layer_name}' from ",
          "{fs::path_file(gdb)}: {e$message}"
        ))
        return(NULL)
      }
    )

    if (is.null(clipped) || nrow(clipped) == 0) return(NULL)

    clipped
  })

  # Remove NULLs and empty results
  layers_list <- purrr::compact(layers_list)

  if (length(layers_list) == 0) return(NULL)

  # Merge all sheets, keeping only geometry column to avoid attribute
  # conflicts between GDB/SHP versions with differing schemas
  merged <- purrr::map(layers_list, function(lyr) {
    dplyr::select(lyr, geometry = dplyr::last_col())
  }) |>
    dplyr::bind_rows()

  merged
}

# -- NHN format-specific readers --------------------------------------------

#' Read an NHN layer from a GDB folder
#' @param gdb_path   Character. Path to .gdb folder
#' @param layer_name Character. Layer name to read
#' @param group_id   Character. Group identifier (for log messages)
#' @return sf object in EPSG:3979, or NULL on failure
read_nhn_from_gdb <- function(gdb_path, layer_name, group_id) {

  available_layers <- tryCatch(
    sf::st_layers(gdb_path)$name,
    error = function(e) {
      cw_warn(glue::glue(
        "Group '{group_id}': could not read layers from ",
        "{fs::path_file(gdb_path)}: {e$message}"
      ))
      return(character(0))
    }
  )

  if (!layer_name %in% available_layers) {
    cw_warn(glue::glue(
      "Group '{group_id}': layer '{layer_name}' not found in ",
      "{fs::path_file(gdb_path)} — skipping."
    ))
    return(NULL)
  }

  tryCatch(
    sf::st_read(gdb_path, layer = layer_name, quiet = TRUE) |>
      sf::st_zm(drop = TRUE, what = "ZM") |>
      sf::st_transform(3979),
    error = function(e) {
      cw_warn(glue::glue(
        "Group '{group_id}': error reading '{layer_name}' from ",
        "{fs::path_file(gdb_path)}: {e$message}"
      ))
      NULL
    }
  )
}

# -- Stream burning ----------------------------------------------------------

#' Write flowlines as a shapefile and burn into DEM using wbt_fill_burn
#'
#' WhiteboxTools' FillBurn requires a vector stream network input and only
#' accepts shapefiles (.shp), not GeoPackages. This function writes the
#' flowlines to a temporary .shp in the group cache, runs wbt_fill_burn(),
#' then removes the temporary shapefile.
#'
#' M/Z coordinate dimensions are dropped before writing as shapefiles do
#' not support them and WhiteboxTools will panic on PolyLineM geometry.
#'
#' @param flowlines       sf object. Flowlines in EPSG:3979
#' @param dem_path        Character. Path to group dem.tif
#' @param dem_burned_path Character. Output path for burned DEM
#' @param group_id        Character. Group identifier (for log messages)
#'
#' @return Invisibly returns dem_burned_path
burn_streams_into_dem <- function(
    flowlines,
    dem_path,
    dem_burned_path,
    group_id
) {

  # Write flowlines to a temporary shapefile in the group cache directory.
  # WhiteboxTools only accepts .shp for vector inputs, not .gpkg.
  streams_shp <- fs::path(fs::path_dir(dem_path), "streams_tmp.shp")

  # -- Filter boundary-crossing stream features ------------------------------
  # NHN features that cross or exit the DEM extent cause wbt_fill_burn to
  # carve anomalously large negative values (e.g. -10000) because WhiteboxTools
  # has no valid downstream context for those segments and enforces drainage
  # continuity by carving to extreme depths.
  # Fix: remove features not fully contained within a 2-cell (60 m) inset of
  # the DEM extent. This is standard preprocessing applied to every group.
  dem_rast           <- terra::rast(dem_path)
  dem_boundary_inset <- terra::buffer(
    terra::as.polygons(terra::ext(dem_rast), crs = terra::crs(dem_rast)),
    width = -60
  )

  # Cast to LINESTRING and drop M/Z before converting to SpatVector —
  # WhiteboxTools requires single-part LINESTRING geometry
  flowlines_lines <- flowlines |>
    sf::st_zm(drop = TRUE, what = "ZM") |>
    sf::st_cast("LINESTRING", warn = FALSE)

  flowlines_vect <- terra::vect(flowlines_lines)
  n_before       <- length(flowlines_vect)

  # terra::relate() returns a matrix — convert to logical vector for subsetting
  within_mat     <- terra::relate(flowlines_vect, dem_boundary_inset,
                                  relation = "within")
  within_idx     <- as.logical(within_mat[, 1])
  flowlines_clipped <- flowlines_vect[within_idx, ]
  n_after        <- length(flowlines_clipped)

  if (n_after < n_before) {
    cw_inform(glue::glue(
      "Group '{group_id}': removed {n_before - n_after} boundary-crossing ",
      "stream feature(s) before burn ({n_after} remaining)."
    ))
  }

  if (n_after == 0) {
    cw_warn(glue::glue(
      "Group '{group_id}': no stream features remain after boundary filter. ",
      "Skipping burn — dem.tif will be used as-is for Whitebox steps."
    ))
    return(invisible(dem_burned_path))
  }

  cw_inform(glue::glue("Group '{group_id}': writing temporary streams shapefile..."))

  terra::writeVector(flowlines_clipped, streams_shp, overwrite = TRUE)

  cw_inform(glue::glue("Group '{group_id}': burning streams into DEM..."))

  # wbt_fill_burn() carves the stream network into the DEM so that flow
  # direction follows the known stream network. No burn depth parameter
  # is required — the function determines the carving depth automatically.
  # Absolute paths are required — WhiteboxTools does not resolve relative
  # paths correctly on all platforms.
  whitebox::wbt_fill_burn(
    dem     = normalizePath(dem_path,        mustWork = TRUE),
    streams = normalizePath(streams_shp,     mustWork = TRUE),
    output  = normalizePath(dem_burned_path, mustWork = FALSE)
  )

  # Remove temporary shapefile components (.shp, .shx, .dbf, .prj, etc.)
  shp_components <- fs::dir_ls(
    fs::path_dir(streams_shp),
    glob = "streams_tmp.*"
  )
  fs::file_delete(shp_components)

  cw_inform(glue::glue("Group '{group_id}': dem_burned.tif written."))

  invisible(dem_burned_path)
}

# -- Empty output helpers ----------------------------------------------------

#' Write empty flowline and waterbody GeoPackages for groups with no NHN data
#'
#' Ensures downstream modules can always read from these paths without
#' needing to check for NHN availability themselves.
#'
#' @param aoi            sf polygon. Group AOI (used to infer CRS)
#' @param flowlines_path Character. Output path for empty flowlines
#' @param waterbodies_path Character. Output path for empty waterbodies
write_empty_nhn_outputs <- function(aoi, flowlines_path, waterbodies_path) {
  write_empty_nhn_output_layer(aoi, flowlines_path,   type = "linestring")
  write_empty_nhn_output_layer(aoi, waterbodies_path, type = "polygon")
}

#' Write a single empty NHN GeoPackage layer
#'
#' @param aoi      sf polygon. Used to determine CRS
#' @param out_path Character. Output path
#' @param type     Character. Geometry type: "linestring" or "polygon"
write_empty_nhn_output_layer <- function(aoi, out_path, type = "linestring") {

  geom_type <- switch(
    type,
    "linestring" = sf::st_sfc(sf::st_linestring(), crs = sf::st_crs(aoi)),
    "polygon"    = sf::st_sfc(sf::st_polygon(),    crs = sf::st_crs(aoi))
  )

  empty_sf <- sf::st_sf(geometry = geom_type[0])

  sf::st_write(empty_sf, out_path, delete_dsn = TRUE, quiet = TRUE)
}
