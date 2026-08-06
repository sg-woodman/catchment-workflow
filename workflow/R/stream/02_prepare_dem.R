# 02_prepare_dem.R
# ---------------------------------------------------------------------------
# Crops the MRDEM virtual raster (VRT) to each group's AOI and reprojects
# to EPSG:3979 (Canada Atlas Lambert NAD83(CSRS)).
#
# The MRDEM is natively in EPSG:3979 so no reprojection is required in most
# cases — the crop step simply extracts the relevant extent. The target_crs
# parameter is retained for flexibility but defaults to 3979 to match the
# MRDEM native CRS, avoiding any resampling.
#
# Inputs:
#   group_manifest : sf tibble from build_group_manifest() in 01_group_sites.R
#   mrdem_vrt      : Path to the MRDEM .vrt file on disk
#   target_crs     : EPSG code for output CRS (default 3979 = MRDEM native)
#
# Outputs (per group, written to cache/<group_id>/):
#   dem.tif        : DEM cropped to group AOI in EPSG:3979
#
# Dependencies: terra, sf, dplyr, purrr, fs, cli (via utils.R)
# ---------------------------------------------------------------------------

#' Prepare the DEM for all groups that need processing
#'
#' Iterates over groups in the manifest, skipping any where dem.tif already
#' exists in the cache. For each group, crops the MRDEM VRT to the group AOI.
#' Since MRDEM is natively in EPSG:3979, no reprojection is needed when using
#' the default target_crs = 3979.
#'
#' @param group_manifest sf tibble from build_group_manifest()
#' @param mrdem_vrt      Character. Path to the MRDEM .vrt file
#' @param target_crs     Integer. EPSG code for output CRS. Default 3979.
#'
#' @return The input group_manifest tibble, invisibly. Called for side effects
#'   (writing dem.tif to each group cache directory).
prepare_dem <- function(
    group_manifest,
    mrdem_vrt,
    target_crs = 3979
) {

  # Validate that the VRT exists before starting any processing
  if (!fs::file_exists(mrdem_vrt)) {
    cw_abort(glue::glue("MRDEM VRT not found at: {mrdem_vrt}"))
  }

  cw_inform(glue::glue("Using MRDEM VRT: {mrdem_vrt}"))
  cw_inform(glue::glue("Target CRS: EPSG:{target_crs}"))

  # Process each group
  purrr::walk(seq_len(nrow(group_manifest)), function(i) {

    grp       <- group_manifest$group_id[i]
    grp_cache <- group_manifest$cache_dir[i]
    aoi       <- group_manifest$aoi[i]
    dem_path  <- fs::path(grp_cache, "dem.tif")

    # Skip if DEM already exists in cache
    if (cache_exists(dem_path)) {
      cw_inform(glue::glue("Group '{grp}': dem.tif found in cache, skipping."))
      return(invisible(NULL))
    }

    cw_inform(glue::glue("Group '{grp}': preparing DEM..."))

    prepare_group_dem(
      aoi        = aoi,
      mrdem_vrt  = mrdem_vrt,
      dem_path   = dem_path,
      target_crs = target_crs,
      group_id   = grp
    )

    cw_inform(glue::glue("Group '{grp}': DEM written to {dem_path}"))
  })

  invisible(group_manifest)
}

#' Prepare the DEM for a single group
#'
#' Internal function. Reprojects the AOI to match the VRT CRS, crops the VRT
#' to the AOI extent, then reprojects to target_crs if needed. Since MRDEM is
#' natively in EPSG:3979 and the AOI is also in EPSG:3979, the reprojection
#' step is a no-op in normal use — no resampling occurs.
#'
#' Uses terra::vect() and terra::project() for the AOI reprojection rather
#' than sf::st_transform() to ensure exact CRS matching with the VRT.
#'
#' @param aoi        sfc polygon. Group AOI in EPSG:3979
#' @param mrdem_vrt  Character. Path to MRDEM .vrt file
#' @param dem_path   Character. Output path for the cropped DEM
#' @param target_crs Integer. EPSG code for output CRS
#' @param group_id   Character. Group identifier (used for log messages only)
#'
#' @return Invisibly returns the output path
prepare_group_dem <- function(
    aoi,
    mrdem_vrt,
    dem_path,
    target_crs,
    group_id
) {

  # Open VRT as a terra SpatRaster (lazy — does not load into memory)
  cw_inform(glue::glue("Group '{group_id}': opening VRT..."))
  vrt <- terra::rast(mrdem_vrt)

  vrt_crs_desc <- terra::crs(vrt, describe = TRUE)
  cw_inform(glue::glue(
    "Group '{group_id}': VRT CRS = ",
    "{vrt_crs_desc$authority}:{vrt_crs_desc$code} ({vrt_crs_desc$name})"
  ))

  # Reproject AOI to match the VRT CRS for cropping. Uses terra::project()
  # on a SpatVector for reliable CRS matching against the VRT definition.
  aoi_vect   <- terra::vect(aoi)
  aoi_native <- terra::project(aoi_vect, terra::crs(vrt))

  # Verify AOI overlaps VRT before attempting crop
  aoi_ext <- terra::ext(aoi_native)
  vrt_ext <- terra::ext(vrt)
  x_overlap <- aoi_ext$xmin <= vrt_ext$xmax && aoi_ext$xmax >= vrt_ext$xmin
  y_overlap <- aoi_ext$ymin <= vrt_ext$ymax && aoi_ext$ymax >= vrt_ext$ymin

  if (!x_overlap || !y_overlap) {
    cw_abort(glue::glue(
      "Group '{group_id}': AOI does not overlap with MRDEM extent. ",
      "Check that pour point coordinates are correct and the MRDEM VRT ",
      "covers the study area."
    ))
  }

  # Crop VRT to the AOI extent — GDAL reads only the requested pixels.
  # snap = 'out' ensures no edge pixels are lost.
  cw_inform(glue::glue("Group '{group_id}': cropping to AOI extent..."))
  dem_cropped <- terra::crop(vrt, aoi_ext, snap = "out")

  # Reproject to target CRS if different from VRT CRS.
  # For MRDEM (EPSG:3979) with target_crs = 3979 this is a no-op.
  # method = "bilinear" is used if reprojection is needed (continuous data).
  vrt_epsg <- as.integer(vrt_crs_desc$code)
  if (!is.na(vrt_epsg) && vrt_epsg != target_crs) {
    cw_inform(glue::glue(
      "Group '{group_id}': reprojecting from EPSG:{vrt_epsg} ",
      "to EPSG:{target_crs}..."
    ))
    dem_cropped <- terra::project(
      dem_cropped,
      paste0("EPSG:", target_crs),
      method = "bilinear"
    )
  } else {
    cw_inform(glue::glue(
      "Group '{group_id}': VRT already in EPSG:{target_crs} — ",
      "no reprojection needed."
    ))
  }

  cw_inform(glue::glue(
    "Group '{group_id}': cropped DEM — ",
    "{terra::nrow(dem_cropped)} rows x {terra::ncol(dem_cropped)} cols, ",
    "resolution: {paste(round(terra::res(dem_cropped), 2), collapse = ' x ')} m"
  ))

  # Write the full rectangular crop to cache — do NOT mask to the AOI polygon.
  # Masking introduces NoData cells along the irregular HydroBasin boundary,
  # which prevents WhiteboxTools from computing flow direction for edge cells
  # and causes incorrect or artificially small catchments near the AOI edge.
  # The rectangular extent is sufficient — watershed delineation naturally
  # constrains outputs to the correct catchment area.
  cw_inform(glue::glue("Group '{group_id}': writing dem.tif..."))
  terra::writeRaster(
    dem_cropped,
    filename  = dem_path,
    filetype  = "GTiff",
    datatype  = "FLT4S",
    overwrite = TRUE,
    gdal      = c("COMPRESS=LZW", "BIGTIFF=IF_SAFER")
  )

  invisible(dem_path)
}

#' Verify DEM outputs for all groups
#'
#' Checks that dem.tif exists and is readable for every group in the manifest.
#' Reports basic raster properties (extent, resolution, CRS) for each group
#' as a sanity check before proceeding to stream burning or Whitebox steps.
#'
#' @param group_manifest sf tibble from build_group_manifest()
#' @return A tibble with one row per group summarising DEM properties:
#'   group_id, dem_path, nrow, ncol, res_x, res_y, crs_epsg
verify_dem_outputs <- function(group_manifest) {

  purrr::map(seq_len(nrow(group_manifest)), function(i) {

    grp      <- group_manifest$group_id[i]
    dem_path <- fs::path(group_manifest$cache_dir[i], "dem.tif")

    if (!cache_exists(dem_path)) {
      cw_warn(glue::glue("Group '{grp}': dem.tif not found at {dem_path}"))
      return(tibble::tibble(
        group_id = grp,
        dem_path = dem_path,
        nrow     = NA_integer_,
        ncol     = NA_integer_,
        res_x    = NA_real_,
        res_y    = NA_real_,
        crs_epsg = NA_character_
      ))
    }

    dem <- terra::rast(dem_path)
    crs_desc <- terra::crs(dem, describe = TRUE)

    tibble::tibble(
      group_id = grp,
      dem_path = dem_path,
      nrow     = terra::nrow(dem),
      ncol     = terra::ncol(dem),
      res_x    = round(terra::res(dem)[1], 2),
      res_y    = round(terra::res(dem)[2], 2),
      crs_epsg = paste0(crs_desc$authority, ":", crs_desc$code)
    )

  }) |>
    dplyr::bind_rows()
}
