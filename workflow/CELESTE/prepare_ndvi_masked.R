# prepare_ndvi_masked.R
# =============================================================================
# Lake-masked variants of the "ndvi"/"ndvi_trend" LOIs (workflow/CELESTE/
# prepare_ndvi.R, prepare_ndvi_trend.R) — a SEPARATE, ADDITIONAL pair of
# LOIs ("ndvi_masked"/"ndvi_trend_masked"), not a toggle applied in place.
# The plain "ndvi"/"ndvi_trend" outputs are never touched by this file —
# they stay a stable, always-available, lake-inclusive baseline to compare
# the masked variant against.
#
# WHY MASKING IS A SEPARATE LOI, NOT A DEFAULT: a lake pixel's NDVI is not
# a terrestrial vegetation reading, and leaving it in systematically drags
# a catchment's summarized NDVI down in proportion to how much open water
# it contains — worth excluding when the goal is a terrestrial signal. An
# earlier version of this masked "ndvi"/"ndvi_trend" was applied IN PLACE
# (overwriting the plain LOI's own cache), which turned out to make a real
# problem (a small set of sites returning no data for an unrelated reason
# — see below) look like a masking bug, since there was no unmasked
# baseline left to compare against. Keeping both as independently-named
# LOIs means nothing is ever destructively replaced, and the two can be
# compared side by side per site.
#
# MASKED AT THE GROUP LEVEL, NOT AFTER PER-SITE CROPPING — deliberate:
# terra::mask() only ever sets values to NA within the same grid; it
# doesn't touch extent or resolution. Whether a pixel is masked before or
# after cropping down to one site's catchment, the final distance-weighted
# summary over that catchment is identical — masking commutes with
# cropping for a simple region mask. Same for the per-pixel Theil-Sen/
# Mann-Kendall trend: sens_slope_trend() fits each pixel's slope from that
# pixel's own time series alone, so masking a neighboring pixel can't
# change any other pixel's fitted slope — mask-then-fit and fit-then-mask
# give byte-identical results. Masking the group mosaic up front therefore
# loses nothing versus masking each site's already-cropped raster, while
# reusing the SAME path_lazy + per-site-crop machinery
# workflow/R/stream/hydroweight_attributes.R already provides for
# "ndvi"/"ndvi_trend" completely unmodified — "ndvi_masked"/
# "ndvi_trend_masked" are just two more named LOIs pointing at their own
# per-group files, not a new mechanism. Per-site inspection
# (cache_dir/hydroweight_loi/ndvi_masked/<site_id>_reprojected.tif) works
# exactly like today's .../ndvi/<site_id>_reprojected.tif.
#
# Masked using NHN waterbody polygons (workflow/R/engine/prepare_lake_
# conditioning.R's fetch_nhn_lakes_for_aoi(), reused rather than
# reimplemented), fetched per group against that group's own NDVI mosaic
# extent — NOT the group's full HydroBasins AOI, which is typically far
# larger (same "tight AOI, not the group's full territory" reasoning
# prepare_lake_conditioning.R's own site-buffer scoping already uses).
#
# KNOWN, SEPARATE ISSUE — not related to masking: MOR-CRANE and MOR-KEN
# have no `canlcc` (CAN_LLC_2020.tif's real extent stops short of their
# location — see README.md's canlcc note) or `harvest_regen` (MOR has no
# coverage from either harvest source at all — see prepare_harvest_
# regen.R) — both real, pre-existing, per-LOI gaps, confirmed directly
# (a genuine "no overlap"/"file not found" per LOI, not a crash). `ndvi`/
# `ndvi_masked` are NOT affected — both compute fine for these 2 sites.
# A first attempt at computing all 4 ndvi-family LOIs together for the
# full site set intermittently returned no data at all for a handful of
# sites scattered across both projects (not consistently the same ones) —
# a heavier-batch instance of the terra temp-file/GC-related flakiness
# already documented elsewhere in this project for large multi-LOI runs.
# Re-running calculate_hydroweight_attributes_stream() for just the
# affected site_ids resolved every case. If a fresh run ever shows a
# similarly scattered set of unexplained NA rows, retry that subset before
# assuming a real bug.
#
# Usage (from run_celeste.R, after sourcing workflow/raster_attributes.R,
# workflow/R/engine/prepare_lake_conditioning.R, and workflow/CELESTE/
# prepare_ndvi.R/prepare_ndvi_trend.R):
#   source(here("workflow/CELESTE/prepare_ndvi_masked.R"))
#   prepare_ndvi_masked_group_rasters(
#     group_manifest, cache_dir = cache_dir,
#     nhn_index_path = nhn_index, nhn_raw_dir = nhn_dir
#   )
#   prepare_ndvi_trend_masked_rasters(
#     group_manifest, cache_dir = cache_dir,
#     nhn_index_path = nhn_index, nhn_raw_dir = nhn_dir
#   )
#
# Dependencies: terra, sf, trend (via sens_slope_trend()), fs, glue;
#   match_group_tiles()/sens_slope_trend()/mask_out_waterbodies() from
#   workflow/raster_attributes.R, clean_ndvi_tiles() from
#   workflow/CELESTE/prepare_ndvi.R, and fetch_nhn_lakes_for_aoi() from
#   workflow/R/engine/prepare_lake_conditioning.R must be sourced first.
# =============================================================================

#' Fetch a group's NHN lake/reservoir polygons covering an NDVI raster's
#' own extent, then mask them out — shared by
#' prepare_ndvi_masked_group_rasters() and prepare_ndvi_trend_masked_rasters()
#' so the fix lives in one place.
#'
#' Thin wrapper around workflow/R/engine/prepare_lake_conditioning.R's
#' fetch_nhn_lakes_for_aoi() — reused rather than reimplemented, since it
#' already does exactly "read NHN waterbodies for an AOI, minimum area,
#' excluded types" correctly (unclipped reads, per-sheet nid-type
#' coercion, etc.). The AOI passed to it is the RASTER's own extent (this
#' group's actual NDVI tile footprint), not the group's full HydroBasins
#' AOI — see this file's header for why that distinction matters for
#' runtime.
#'
#' @param raster    SpatRaster to mask (a group's mosaicked NDVI raster or
#'   NDVI-trend output).
#' @param group_id  Character. Group identifier (log messages, NHN sheet
#'   lookup).
#' @param nhn_index_path Character. Path to the NHN index shapefile.
#' @param nhn_raw_dir    Character. Path to the directory of NHN GDB sheet
#'   folders.
#' @param min_lake_area_ha Numeric. Minimum lake area (ha) to mask —
#'   smaller features are typically shoreline-generalization slivers, not
#'   worth excising. Default 1, matching prepare_lake_conditioning.R's own
#'   default.
#' @param exclude_waterbody_types Character vector of NHN
#'   waterDefinitionText values to NOT mask (i.e. not treated as open
#'   water for this purpose). Default character(0) — mask every polygon
#'   NHN's waterbody layer returns for this footprint, deliberately
#'   broader than the c("Watercourse") exclusion this project's lake-
#'   BISECTION checks use elsewhere: a watercourse is still open water and
#'   still depresses terrestrial NDVI, even though it isn't a "lake" for
#'   bisection purposes.
#' @return `raster`, unchanged if no lakes are found for this footprint,
#'   or with every cell touching a qualifying lake polygon set to NA (see
#'   mask_out_waterbodies()).
mask_ndvi_lakes_celeste <- function(
  raster, group_id, nhn_index_path, nhn_raw_dir,
  min_lake_area_ha = 1, exclude_waterbody_types = character(0)
) {
  if (!exists("fetch_nhn_lakes_for_aoi", mode = "function")) {
    cw_abort(paste(
      "mask_ndvi_lakes_celeste() requires workflow/R/engine/prepare_lake_",
      "conditioning.R to be sourced first (defines fetch_nhn_lakes_for_aoi())."
    ))
  }
  if (!exists("mask_out_waterbodies", mode = "function")) {
    cw_abort(paste(
      "mask_ndvi_lakes_celeste() requires workflow/raster_attributes.R to",
      "be sourced first (defines mask_out_waterbodies())."
    ))
  }

  aoi_sf <- sf::st_as_sf(terra::as.polygons(terra::ext(raster), crs = terra::crs(raster)))

  lakes <- tryCatch(
    fetch_nhn_lakes_for_aoi(
      aoi = aoi_sf, nhn_index_path = nhn_index_path, nhn_raw_dir = nhn_raw_dir,
      group_id = group_id, min_area_ha = min_lake_area_ha,
      exclude_types = exclude_waterbody_types
    ),
    error = function(e) {
      cw_warn(glue::glue(
        "Group '{group_id}': lake fetch for NDVI masking failed — ",
        "{e$message}. Leaving NDVI unmasked for this group."
      ))
      NULL
    }
  )

  if (is.null(lakes) || nrow(lakes) == 0) {
    cw_inform(glue::glue("Group '{group_id}': no lakes found to mask from NDVI."))
    return(raster)
  }

  cw_inform(glue::glue("Group '{group_id}': masking {nrow(lakes)} lake(s) out of NDVI."))
  mask_out_waterbodies(raster, lakes)
}

#' Materialize one lake-masked multi-band NDVI GeoTIFF per group — the
#' "ndvi_masked" LOI's prep function
#'
#' Identical to prepare_ndvi.R's prepare_ndvi_per_group_rasters() (same
#' mosaicking of a group's own source tiles) plus a mask_ndvi_lakes_celeste()
#' call before writing, and a different output directory
#' (hydroweight_loi/ndvi_masked/ instead of hydroweight_loi/ndvi/) so the
#' plain, unmasked "ndvi" cache is never touched.
#'
#' @param group_manifest sf tibble from build_group_manifest() — only
#'   group_id is used.
#' @param ndvi_dir  Character. Directory of raw source NDVI tiles. Default
#'   here("data/ndvi").
#' @param cache_dir Character. Project cache root — rasters are written to
#'   cache_dir/hydroweight_loi/ndvi_masked/<group_id>.tif.
#' @param nhn_index_path  Character. Path to the NHN index shapefile.
#' @param nhn_raw_dir     Character. Path to the directory of NHN GDB
#'   sheet folders.
#' @param min_lake_area_ha Numeric. Passed to mask_ndvi_lakes_celeste().
#'   Default 1.
#' @param exclude_waterbody_types Character vector. Passed to
#'   mask_ndvi_lakes_celeste(). Default character(0).
#'
#' @return Character vector of the (existing-or-newly-written) per-group
#'   raster paths, invisibly. A group with no matching tiles is skipped
#'   with a warning, same as prepare_ndvi_per_group_rasters().
prepare_ndvi_masked_group_rasters <- function(
  group_manifest,
  ndvi_dir = here::here("data/ndvi"),
  cache_dir,
  nhn_index_path,
  nhn_raw_dir,
  min_lake_area_ha = 1,
  exclude_waterbody_types = character(0)
) {
  out_dir <- fs::path(cache_dir, "hydroweight_loi", "ndvi_masked")
  fs::dir_create(out_dir)
  files <- clean_ndvi_tiles(ndvi_dir = ndvi_dir, cache_dir = cache_dir)
  written <- character(0)

  for (grp in unique(group_manifest$group_id)) {
    group_raster_path <- fs::path(out_dir, paste0(grp, ".tif"))
    written <- c(written, group_raster_path)
    if (cache_exists(group_raster_path)) {
      next
    }

    matched <- match_group_tiles(files, grp)
    if (length(matched) == 0) {
      cw_warn(glue::glue(
        "No NDVI tiles matched for group '{grp}'; skipping per-group ",
        "ndvi_masked raster (will read as NA downstream)."
      ))
      next
    }

    mini_vrt <- terra::vrt(
      matched,
      filename = tempfile(fileext = ".vrt"),
      options = c("-vrtnodata", "nan"),
      overwrite = TRUE
    )
    # Same band-naming pin as prepare_ndvi_per_group_rasters() — see that
    # function's comment for why terra::vrt()'s own default (derived from
    # a random tempfile() name) is non-deterministic and must be overridden.
    names(mini_vrt) <- paste0("ndvi_mosaic_", seq_len(terra::nlyr(mini_vrt)))

    mini_vrt <- mask_ndvi_lakes_celeste(
      mini_vrt, grp, nhn_index_path, nhn_raw_dir,
      min_lake_area_ha, exclude_waterbody_types
    )

    cw_inform(glue::glue(
      "Materializing lake-masked NDVI raster for group '{grp}' ",
      "({length(matched)} tile(s), {terra::ncell(mini_vrt)} cells)..."
    ))

    terra::writeRaster(
      mini_vrt,
      group_raster_path,
      overwrite = TRUE,
      datatype = "FLT4S"
    )
  }

  invisible(written)
}

#' Rasterize the lake-masked NDVI trend (Theil-Sen slope + Mann-Kendall
#' p-value) LOI, once per group — the "ndvi_trend_masked" LOI's prep
#' function
#'
#' Identical to prepare_ndvi_trend.R's prepare_ndvi_trend_rasters() plus a
#' mask_ndvi_lakes_celeste() call before sens_slope_trend() runs, and a
#' different output directory (hydroweight_loi/ndvi_trend_masked/ instead
#' of hydroweight_loi/ndvi_trend/) so the plain, unmasked "ndvi_trend"
#' cache is never touched. Masking BEFORE fitting the trend (rather than
#' masking the trend output afterward) costs nothing extra in wall time —
#' see this file's header for why the two orders give identical results —
#' and, if anything, is slightly cheaper (fewer valid pixels for
#' sens_slope_trend() to fit, which is where this function's entire
#' runtime cost lives).
#'
#' @inheritParams prepare_ndvi_masked_group_rasters
#' @param years     Integer vector, same length/order as each tile's band
#'   count. Default 1984:2025.
#' @param min_obs   Passed to sens_slope_trend(). Default 4.
#'
#' @return Character vector of the (existing-or-newly-written) per-group
#'   raster paths, invisibly.
prepare_ndvi_trend_masked_rasters <- function(
  group_manifest,
  ndvi_dir = here::here("data/ndvi"),
  cache_dir,
  nhn_index_path,
  nhn_raw_dir,
  min_lake_area_ha = 1,
  exclude_waterbody_types = character(0),
  years = 1984:2025,
  min_obs = 4
) {
  out_dir <- fs::path(cache_dir, "hydroweight_loi", "ndvi_trend_masked")
  files <- clean_ndvi_tiles(ndvi_dir = ndvi_dir, cache_dir = cache_dir)
  written <- character(0)

  for (grp in unique(group_manifest$group_id)) {
    group_raster_path <- fs::path(out_dir, paste0(grp, ".tif"))
    written <- c(written, group_raster_path)
    if (cache_exists(group_raster_path)) next

    matched <- match_group_tiles(files, grp)
    if (length(matched) == 0) {
      cw_warn(glue::glue(
        "No NDVI tiles matched for group '{grp}'; skipping ndvi_trend_masked ",
        "raster (will read as NA downstream)."
      ))
      next
    }

    fs::dir_create(out_dir)

    mini_vrt <- terra::vrt(
      matched,
      filename = tempfile(fileext = ".vrt"),
      options = c("-vrtnodata", "nan"),
      overwrite = TRUE
    )

    mini_vrt <- mask_ndvi_lakes_celeste(
      mini_vrt, grp, nhn_index_path, nhn_raw_dir,
      min_lake_area_ha, exclude_waterbody_types
    )

    cw_inform(glue::glue(
      "Computing lake-masked NDVI trend for group '{grp}' ",
      "({length(matched)} tile(s), {terra::ncell(mini_vrt)} cells) — ",
      "this can take several minutes."
    ))

    trend <- sens_slope_trend(mini_vrt, x = years, min_obs = min_obs)

    terra::writeRaster(trend, group_raster_path, overwrite = TRUE)
  }

  invisible(written)
}
