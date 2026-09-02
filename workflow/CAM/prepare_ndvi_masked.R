# prepare_ndvi_masked.R
# =============================================================================
# Lake-masked variants of the "ndvi"/"ndvi_trend" LOIs (workflow/CAM/
# prepare_ndvi.R, prepare_ndvi_trend.R) — a SEPARATE, ADDITIONAL pair of
# LOIs ("ndvi_masked"/"ndvi_trend_masked"), not a toggle applied in place.
# The plain "ndvi"/"ndvi_trend" outputs (and their shared
# hydroweight_loi/ndvi_clean/ source) are never touched by this file —
# they stay a stable, always-available, lake-inclusive baseline to compare
# the masked variant against. See workflow/CELESTE/prepare_ndvi_masked.R's
# header for the fuller rationale (why a separate LOI rather than an
# in-place toggle, and why masking commutes with per-site cropping and
# with the per-pixel Theil-Sen trend, so masking the shared per-file
# raster up front loses nothing versus masking after cropping).
#
# Masked using the OIH/OHN waterbody layer (ohn_waterbodies_valid.gpkg,
# CAM_OIH_LAKES_PATH below — same file match_lake_polygons.R/
# remove_upstream_lakes.R already use), read via an extent-filtered
# terra::vect() scoped to each raw file's own footprint (same GDAL
# bbox-pushdown pattern those two files already use against this exact
# file — 1.4M features province-wide, never loaded whole).
#
# Applied to a SEPARATE masked-clean cache
# (hydroweight_loi/ndvi_masked_clean/<basename>), built from the plain
# (unmasked) hydroweight_loi/ndvi_clean/<basename> file — not from the raw
# export directly — so the 0->NA/10000 rescale fix stays defined in
# exactly one place (prepare_ndvi.R's clean_cam_ndvi_tile()).
#
# Usage (from run_cam_streams.R, after sourcing workflow/gee_utils.R,
# workflow/raster_attributes.R, and workflow/CAM/prepare_ndvi.R/
# prepare_ndvi_trend.R):
#   source(here("workflow/CAM/prepare_ndvi_masked.R"))
#   prepare_cam_ndvi_masked_site_rasters(sites, output_dir = output_dir, cache_dir = cache_dir)
#   prepare_cam_ndvi_trend_masked_site_rasters(sites, output_dir = output_dir, cache_dir = cache_dir)
#
# Dependencies: terra, sf, purrr, fs, glue (via utils.R); build_cam_ndvi_
#   site_map()/clean_cam_ndvi_tile() from workflow/CAM/prepare_ndvi.R,
#   sens_slope_trend()/mask_out_waterbodies() from workflow/raster_
#   attributes.R must be sourced first.
# =============================================================================

# Default location of the OIH/OHN waterbody layer used to mask lakes out
# of the NDVI LOIs — same file workflow/R/lake/match_lake_polygons.R,
# remove_upstream_lakes.R, and workflow/CAM/fix_lake_bisection.R already
# use.
CAM_OIH_LAKES_PATH <- "/Users/sam/Documents/cfs/shared_data/raw/hydro/waterbodies/ON_waterbodies/ohn_waterbodies_valid.gpkg"

#' Fetch OIH/OHN waterbody polygons covering a raster's own extent, then
#' mask them out — keeps the "ndvi_masked"/"ndvi_trend_masked" LOIs
#' terrestrial-only (see this file's header).
#'
#' Reads via an extent-filtered terra::vect() (same GDAL bbox-pushdown
#' pattern workflow/R/lake/match_lake_polygons.R and
#' remove_upstream_lakes.R already use against this exact file) rather
#' than loading all ~1.4M province-wide features — safe here specifically
#' because ohn_waterbodies_valid.gpkg is pre-validated (no curved/
#' MULTISURFACE geometries), unlike the raw NHN/harvest gdbs
#' workflow/raster_attributes.R's docstring warns about.
#'
#' @param raster    SpatRaster to mask (a cleaned NDVI export tile).
#' @param lakes_path Character. Path to the OIH/OHN waterbody layer.
#' @param min_lake_area_ha Numeric. Minimum lake area (ha) to mask.
#'   Default 1, matching this project's other waterbody consumers.
#' @param exclude_waterbody_types Character vector of WATERBODY_ values to
#'   NOT mask. Default character(0) — mask every polygon type this layer
#'   returns for this footprint (Lake, Pond, Reservoir, River, etc.):
#'   unlike this project's lake-BISECTION checks (which exclude "River"/
#'   "Pond" because those aren't a "lake" for that purpose), a river or
#'   pond surface still isn't terrestrial and still depresses NDVI, so
#'   nothing is excluded by default here.
#' @param label     Character. Used in log messages only.
#' @return `raster`, unchanged if no waterbodies are found for this
#'   footprint, or with every cell touching a qualifying polygon set to NA
#'   (see mask_out_waterbodies()).
mask_cam_ndvi_lakes <- function(
  raster, lakes_path, min_lake_area_ha = 1,
  exclude_waterbody_types = character(0), label = ""
) {
  if (!exists("mask_out_waterbodies", mode = "function")) {
    cw_abort(paste(
      "mask_cam_ndvi_lakes() requires workflow/raster_attributes.R to be",
      "sourced first (defines mask_out_waterbodies())."
    ))
  }

  lakes_here <- tryCatch(
    {
      # Read one row to get the native CRS of the lake layer (same pattern
      # match_lake_polygons.R uses), then reproject the query extent to it
      # so the GDAL bbox filter below is applied in the source's own CRS.
      poly_crs <- terra::crs(
        terra::vect(lakes_path, what = "geoms", extent = terra::ext(0, 0, 0, 0))
      )
      ext_native <- terra::project(
        terra::as.polygons(terra::ext(raster), crs = terra::crs(raster)), poly_crs
      )
      cand <- terra::vect(lakes_path, extent = terra::ext(ext_native))
      if (length(exclude_waterbody_types) > 0 && "WATERBODY_" %in% names(cand)) {
        cand <- cand[!cand$WATERBODY_ %in% exclude_waterbody_types, ]
      }
      if (nrow(cand) > 0) {
        cand <- cand[terra::expanse(cand, unit = "ha") >= min_lake_area_ha, ]
      }
      terra::project(cand, terra::crs(raster))
    },
    error = function(e) {
      cw_warn(glue::glue(
        "{label}: lake fetch for NDVI masking failed — {e$message}. ",
        "Leaving unmasked."
      ))
      NULL
    }
  )

  if (is.null(lakes_here) || nrow(lakes_here) == 0) {
    cw_inform(glue::glue("{label}: no waterbodies found to mask from NDVI."))
    return(raster)
  }

  cw_inform(glue::glue(
    "{label}: masking {nrow(lakes_here)} waterbody polygon(s) out of NDVI."
  ))
  mask_out_waterbodies(raster, lakes_here)
}

#' Build (or reuse) a lake-masked version of one distinct cleaned NDVI
#' file — shared by prepare_cam_ndvi_masked_site_rasters() and
#' prepare_cam_ndvi_trend_masked_site_rasters(), so masking only needs
#' implementing once, mirroring how clean_cam_ndvi_tile() itself is shared
#' between the two plain (unmasked) LOI functions.
#'
#' @param raw_path  Character. Path to one raw NDVI export .tif (the same
#'   key clean_cam_ndvi_tile()/build_cam_ndvi_site_map() use).
#' @param cache_dir Character. Project cache root.
#' @param lakes_path Character. Passed to mask_cam_ndvi_lakes(). Default
#'   CAM_OIH_LAKES_PATH.
#' @param min_lake_area_ha Numeric. Passed to mask_cam_ndvi_lakes().
#'   Default 1.
#' @param exclude_waterbody_types Character vector. Passed to
#'   mask_cam_ndvi_lakes(). Default character(0).
#' @return Character path to the masked-clean file
#'   (cache_dir/hydroweight_loi/ndvi_masked_clean/<basename>).
clean_and_mask_cam_ndvi_tile <- function(
  raw_path, cache_dir, lakes_path = CAM_OIH_LAKES_PATH,
  min_lake_area_ha = 1, exclude_waterbody_types = character(0)
) {
  out_dir <- fs::path(cache_dir, "hydroweight_loi", "ndvi_masked_clean")
  fs::dir_create(out_dir, recurse = TRUE)
  out_path <- fs::path(out_dir, fs::path_file(raw_path))

  if (!cache_exists(out_path)) {
    cleaned_path <- clean_cam_ndvi_tile(raw_path, cache_dir) # plain, unmasked — never overwritten
    r <- terra::rast(cleaned_path)
    r <- mask_cam_ndvi_lakes(
      r, lakes_path, min_lake_area_ha, exclude_waterbody_types,
      label = fs::path_file(raw_path)
    )
    terra::writeRaster(r, out_path, overwrite = TRUE, datatype = "FLT4S")
  }

  out_path
}

#' Materialize a per-site lake-masked NDVI raster for every site — the
#' "ndvi_masked" LOI's prep function
#'
#' Mirrors prepare_ndvi.R's prepare_cam_ndvi_site_rasters() exactly
#' (clean each distinct file once, symlink per member site), except the
#' shared file it symlinks to is the MASKED variant
#' (clean_and_mask_cam_ndvi_tile()), cached and linked under
#' hydroweight_loi/ndvi_masked_clean/ and hydroweight_loi/ndvi_masked/
#' respectively — the plain hydroweight_loi/ndvi_clean/ and
#' hydroweight_loi/ndvi/ directories are untouched.
#'
#' @param sites      Sites tibble (site_id column).
#' @param output_dir Character. Root output directory.
#' @param cache_dir  Character. Project cache root.
#' @param ndvi_dir   Character. Directory of raw NDVI export files.
#'   Default CAM_NDVI_DIR.
#' @param lakes_path Character. Passed to clean_and_mask_cam_ndvi_tile().
#'   Default CAM_OIH_LAKES_PATH.
#' @param min_lake_area_ha Numeric. Passed to
#'   clean_and_mask_cam_ndvi_tile(). Default 1.
#' @param exclude_waterbody_types Character vector. Passed to
#'   clean_and_mask_cam_ndvi_tile(). Default character(0).
#'
#' @return The site map (see build_cam_ndvi_site_map()), invisibly.
prepare_cam_ndvi_masked_site_rasters <- function(
  sites, output_dir, cache_dir, ndvi_dir = CAM_NDVI_DIR,
  lakes_path = CAM_OIH_LAKES_PATH,
  min_lake_area_ha = 1, exclude_waterbody_types = character(0)
) {
  site_map <- build_cam_ndvi_site_map(sites, output_dir, ndvi_dir)

  out_dir <- fs::path(cache_dir, "hydroweight_loi", "ndvi_masked")
  fs::dir_create(out_dir, recurse = TRUE)

  matched <- dplyr::filter(site_map, !is.na(raw_path))
  distinct_raw <- unique(matched$raw_path)
  masked_lookup <- setNames(
    purrr::map_chr(
      distinct_raw, clean_and_mask_cam_ndvi_tile, cache_dir = cache_dir,
      lakes_path = lakes_path, min_lake_area_ha = min_lake_area_ha,
      exclude_waterbody_types = exclude_waterbody_types
    ),
    distinct_raw
  )

  purrr::pwalk(matched, function(site_id, leader, raw_path) {
    target <- fs::path(out_dir, paste0(site_id, ".tif"))
    if (cache_exists(target)) {
      return(invisible(NULL))
    }
    if (fs::link_exists(target) || fs::file_exists(target)) {
      fs::file_delete(target)
    }
    fs::link_create(masked_lookup[[raw_path]], target, symbolic = TRUE)
  })

  cw_inform(glue::glue(
    "NDVI masked (continuous): {nrow(matched)}/{nrow(site_map)} site(s) ",
    "linked from {length(distinct_raw)} distinct export file(s)."
  ))

  invisible(site_map)
}

#' Rasterize the lake-masked NDVI trend LOI for every site, computed once
#' per distinct source file — the "ndvi_trend_masked" LOI's prep function
#'
#' Mirrors prepare_ndvi_trend.R's prepare_cam_ndvi_trend_site_rasters()
#' exactly, except it computes sens_slope_trend() on the MASKED-clean file
#' (clean_and_mask_cam_ndvi_tile()) and caches/links under
#' hydroweight_loi/ndvi_trend_src_masked/ and hydroweight_loi/
#' ndvi_trend_masked/ — the plain hydroweight_loi/ndvi_trend_src/ and
#' hydroweight_loi/ndvi_trend/ directories are untouched.
#'
#' @inheritParams prepare_cam_ndvi_masked_site_rasters
#' @param years   Integer vector. Default 1984:2025.
#' @param min_obs Passed to sens_slope_trend(). Default 4.
#'
#' @return The site map (see build_cam_ndvi_site_map()), invisibly.
prepare_cam_ndvi_trend_masked_site_rasters <- function(
  sites, output_dir, cache_dir, ndvi_dir = CAM_NDVI_DIR,
  lakes_path = CAM_OIH_LAKES_PATH,
  min_lake_area_ha = 1, exclude_waterbody_types = character(0),
  years = 1984:2025, min_obs = 4
) {
  site_map <- build_cam_ndvi_site_map(sites, output_dir, ndvi_dir)

  shared_dir <- fs::path(cache_dir, "hydroweight_loi", "ndvi_trend_src_masked")
  out_dir <- fs::path(cache_dir, "hydroweight_loi", "ndvi_trend_masked")
  fs::dir_create(shared_dir, recurse = TRUE)
  fs::dir_create(out_dir, recurse = TRUE)

  matched <- dplyr::filter(site_map, !is.na(raw_path))
  distinct_raw <- unique(matched$raw_path)

  shared_lookup <- setNames(
    purrr::map_chr(distinct_raw, function(raw_path) {
      shared_path <- fs::path(shared_dir, fs::path_file(raw_path))
      if (!cache_exists(shared_path)) {
        masked_clean <- clean_and_mask_cam_ndvi_tile(
          raw_path, cache_dir, lakes_path = lakes_path,
          min_lake_area_ha = min_lake_area_ha,
          exclude_waterbody_types = exclude_waterbody_types
        )
        r <- terra::rast(masked_clean)
        cw_inform(glue::glue(
          "Computing lake-masked NDVI trend for {fs::path_file(raw_path)} ",
          "({terra::ncell(r)} cells) — this can take a while for a large group..."
        ))
        trend <- sens_slope_trend(r, x = years, min_obs = min_obs)
        terra::writeRaster(trend, shared_path, overwrite = TRUE)
      }
      shared_path
    }),
    distinct_raw
  )

  purrr::pwalk(matched, function(site_id, leader, raw_path) {
    target <- fs::path(out_dir, paste0(site_id, ".tif"))
    if (cache_exists(target)) {
      return(invisible(NULL))
    }
    if (fs::link_exists(target) || fs::file_exists(target)) {
      fs::file_delete(target)
    }
    fs::link_create(shared_lookup[[raw_path]], target, symbolic = TRUE)
  })

  cw_inform(glue::glue(
    "NDVI trend masked: {nrow(matched)}/{nrow(site_map)} site(s) linked ",
    "from {length(distinct_raw)} distinct export file(s)."
  ))

  invisible(site_map)
}
