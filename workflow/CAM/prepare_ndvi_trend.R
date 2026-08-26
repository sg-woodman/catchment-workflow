# prepare_ndvi_trend.R
# =============================================================================
# Prepares the CAM stream-site NDVI trend LOI raster (Theil-Sen slope +
# Mann-Kendall p-value, 1984-2025) for the stream hydroweight stage
# (workflow/CAM/run_cam_streams.R, Stage 7). Same approach as
# workflow/CELESTE/prepare_ndvi_trend.R (calculate a per-pixel Sen's slope
# trend on the cleaned NDVI time series), adapted for CAM's per-group
# (not per-region) export layout — see workflow/CAM/prepare_ndvi.R's
# header for the full explanation of that difference and the site->file
# mapping it builds.
#
# One raster per DISTINCT raw export file (computed once, even though it
# may cover several sites — see below), then exposed per-site via the same
# symlink mechanism prepare_ndvi.R uses for the continuous LOI. 2 bands
# ("slope", "p_value") — matches sens_slope_trend()'s return shape
# exactly, so no reshaping is needed before handing it to the hydroweight
# stage as a path_template LOI.
#
# Operates on the CLEANED files from prepare_ndvi.R's clean_cam_ndvi_tile()
# (0 -> NA masked-pixel/outside-AOI fix + /10000 rescale), not the raw
# data/ndvi/CAM/*.tif exports directly — same reasoning as CELESTE's: a
# raw 0 is Earth Engine's export fill value, not a genuine annual NDVI
# observation, and treating it as real would corrupt the Theil-Sen slope
# with false "crashes to zero" instead of correctly excluding that
# pixel-year via min_obs.
#
# COMPUTED ONCE PER FILE, NOT PER SITE: sens_slope_trend()'s cost is
# ~0.5 ms per VALID pixel (CELESTE's own measured rate — see that
# project's prepare_ndvi_trend.R and CLAUDE.md's "known quirks" section),
# essentially independent of which/how-many sites a file happens to
# cover. Since a CAM group file (e.g. the 5-member NCMN/SUD101/SUD102/
# SUD103/SUD200 group) already covers every one of its member sites'
# catchments in one raster, computing the trend once and symlinking the
# SAME output to every member is both correct (no site loses coverage —
# member sites' own catchments are all inside the shared file's extent)
# and the obviously cheaper choice (recomputing an identical trend 5x for
# 5 sites sharing one file would be pure waste).
#
# ESTIMATED COST (confirmed directly, summing each file's non-zero/non-NA
# band-1 cell count as a valid-pixel proxy across all 27 distinct CAM NDVI
# files, at CELESTE's measured ~0.5 ms/valid-pixel rate): ~43 minutes
# total, one-time, cached to disk per distinct file — not repeated on
# subsequent run_cam_streams.R sources. Dominated by the largest merged-
# group files (SUD12+VER01, SUD17+SUD22, the NCMN cluster, and the
# ILD02-exported-as-"idl02" file — each in the ~125K-2.2M valid-cell
# range); every standalone site's own file is comparatively tiny (tens to
# low-thousands of cells) and finishes in well under a second. Stays
# single-threaded (cores = 1, sens_slope_trend()'s default) — CELESTE's
# own testing found cores > 1 made this dramatically SLOWER (a PSOCK-
# cluster overhead regression, not a parallelization win); no reason to
# expect that to differ here, so not re-tested.
#
# Usage (from run_cam_streams.R, after sourcing workflow/gee_utils.R,
# workflow/raster_attributes.R, and workflow/CAM/prepare_ndvi.R):
#   source(here("workflow/CAM/prepare_ndvi_trend.R"))
#   prepare_cam_ndvi_trend_site_rasters(sites, output_dir = output_dir, cache_dir = cache_dir)
#
# Dependencies: terra, sf, purrr, fs, glue (via utils.R); sens_slope_trend()
#   from workflow/raster_attributes.R and build_cam_ndvi_site_map()/
#   clean_cam_ndvi_tile() from workflow/CAM/prepare_ndvi.R must be sourced
#   first.
# =============================================================================

#' Rasterize the NDVI trend (Theil-Sen slope + Mann-Kendall p-value) LOI
#' for every site, computed once per distinct source file
#'
#' Skips any distinct file whose trend raster already exists
#' (cache_exists()) — safe to call on every run_cam_streams.R source.
#'
#' @param sites      Sites tibble (site_id column).
#' @param output_dir Character. Root output directory (catchments read
#'   from here — see build_cam_ndvi_site_map()).
#' @param cache_dir  Character. Project cache root — the shared per-file
#'   trend raster is written to cache_dir/hydroweight_loi/ndvi_trend_src/
#'   <raw file's basename>, and each site gets a symlink to it at
#'   cache_dir/hydroweight_loi/ndvi_trend/<site_id>.tif (the path
#'   06_hydroweight_attributes.R's resolve_site_loi_raster() expects for a
#'   path_template LOI).
#' @param ndvi_dir   Character. Directory of raw NDVI export files.
#'   Default CAM_NDVI_DIR (from prepare_ndvi.R).
#' @param years      Integer vector, same length/order as each file's band
#'   count — the calendar year each band represents. Default 1984:2025
#'   (matches these exports' native range — confirmed directly, 42 bands
#'   named "{index}_NDVI_{year}").
#' @param min_obs    Passed to sens_slope_trend() — minimum non-NA years
#'   required to fit a trend for a given pixel. Default 4.
#'
#' @return The site map (see build_cam_ndvi_site_map()), invisibly.
prepare_cam_ndvi_trend_site_rasters <- function(
  sites, output_dir, cache_dir, ndvi_dir = CAM_NDVI_DIR,
  years = 1984:2025, min_obs = 4
) {
  if (!exists("build_cam_ndvi_site_map", mode = "function")) {
    cw_abort(paste(
      "prepare_cam_ndvi_trend_site_rasters() requires",
      "workflow/CAM/prepare_ndvi.R to be sourced first."
    ))
  }
  if (!exists("sens_slope_trend", mode = "function")) {
    cw_abort(paste(
      "prepare_cam_ndvi_trend_site_rasters() requires",
      "workflow/raster_attributes.R to be sourced first."
    ))
  }

  site_map <- build_cam_ndvi_site_map(sites, output_dir, ndvi_dir)

  shared_dir <- fs::path(cache_dir, "hydroweight_loi", "ndvi_trend_src")
  out_dir <- fs::path(cache_dir, "hydroweight_loi", "ndvi_trend")
  fs::dir_create(shared_dir, recurse = TRUE)
  fs::dir_create(out_dir, recurse = TRUE)

  matched <- dplyr::filter(site_map, !is.na(raw_path))
  distinct_raw <- unique(matched$raw_path)

  shared_lookup <- setNames(
    purrr::map_chr(distinct_raw, function(raw_path) {
      shared_path <- fs::path(shared_dir, fs::path_file(raw_path))
      if (!cache_exists(shared_path)) {
        cleaned <- clean_cam_ndvi_tile(raw_path, cache_dir)
        r <- terra::rast(cleaned)
        cw_inform(glue::glue(
          "Computing NDVI trend for {fs::path_file(raw_path)} ",
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
    "NDVI trend: {nrow(matched)}/{nrow(site_map)} site(s) linked ",
    "from {length(distinct_raw)} distinct export file(s)."
  ))

  invisible(site_map)
}
