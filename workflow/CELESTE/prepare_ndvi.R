# prepare_ndvi.R
# =============================================================================
# Prepares the CELESTE NDVI (continuous) LOI raster for the stream
# hydroweight stage (workflow/run_celeste.R, Stage 9). Kept as its own
# script — same reasoning as workflow/prepare_harvest_regen.R: this is a
# workflow-specific data-prep concern (source data, how it's tiled, what
# needs mosaicking/reprojecting), and how it's handled will genuinely
# differ between projects (a continuous time-series mosaic here; a
# rasterized-per-group categorical layer with temporal precedence for
# harvest/regen; something else entirely for a future project) — none of
# that belongs baked into the standardized Stage 1-9 runner or the generic
# hydroweight module.
#
# Source: data/ndvi/*.tif — 12 regional tiles (COC, KEN1-4, MOR, NIP1-5,
# TUR), 1984-2025 annual composites, all already EPSG:3979/30m/42-band-
# aligned (verified with check_tile_consistency() — no reprojection or band
# alignment needed here; a future batch that isn't this consistent would
# need that handled in this script before mosaicking, not downstream).
#
# RAW TILE ENCODING (Google Earth Engine export, code/landsat_annual_comp_
# slc.js — see get_indices() and exportVariableTimeSeries() there): NDVI is
# scaled by 10000 and stored as int16 (".multiply(10000).round().short()")
# — real file-size savings for a 42-band time series. A raw value of 7394
# is a true NDVI of 0.7394; a raw value of exactly 0 is Earth Engine's
# masked-pixel fill value on export (no formatOptions.noData was set in the
# export call, so a fully-masked composite pixel is written as literal 0,
# not tagged as NoData in the GeoTIFF), NOT a genuine NDVI of 0.
#
# Reviewed the GEE script for whether these 0/missing pixels are more
# likely a processing error or an expected consequence of cloud masking —
# it's the latter, and by design, not a bug:
#   - maskL457sr()/maskL89sr() mask fill, dilated cloud, cloud, cloud
#     shadow (and cirrus for L8/9) pixels via QA_PIXEL bits PER SCENE,
#     before compositing — standard practice, not overly aggressive.
#   - Each year's composite is ee.Reducer.median() over every unmasked
#     Landsat 4/5/7/8/9 scene from May-Sept (mrange) of that year. Masked
#     inputs are excluded from the median at each pixel individually — if
#     literally every scene available for a pixel in that year's window is
#     masked (persistent cloud, or simply few/no clear passes), the
#     composite pixel is masked, which is a legitimate "no valid
#     observation this year" outcome, not an error.
#   - Landsat 7's SLC-off scan-line gaps get an explicit, deliberate
#     mitigation (fillGaps(): a 10-iteration focal_mean blend) before
#     compositing — evidence this was already anticipated and handled as
#     far as is practical, not overlooked.
#   - The scene-level ee.Filter.lt('CLOUD_COVER', 20) filter drops whole
#     scenes with >20% cloud cover from the collection outright (even if
#     the specific catchment happened to be clear in an otherwise-cloudy
#     scene), which thins the pool of candidate observations per pixel-
#     year further — a reasonable simplifying tradeoff (much cheaper than
#     per-pixel-only filtering) but one that will produce more empty
#     pixel-years, especially early in the record (1984-1998, Landsat 4/5
#     only — half the sensor coverage of later years) or over
#     persistently-cloudy areas/years.
# Net: the fraction of 0/missing pixels is expected to vary by group and
# year for real physical/data-availability reasons (regional cloud
# climatology, sensor era), not indicate a bug to chase down. Worth
# confirming empirically per group if a particular group's missing
# fraction looks implausibly high, but nothing in the script itself looks
# wrong.
#
# clean_ndvi_tiles() (below) applies both fixes — 0 -> NA, then /10000 —
# once per raw tile, cached, and is shared by this script's
# prepare_ndvi_vrt() (whole-catchment continuous LOI) AND workflow/CELESTE/
# prepare_ndvi_trend.R's prepare_ndvi_trend_rasters() (per-group Sen's-
# slope trend), so the fix lives in exactly one place rather than needing
# to be kept in sync across both.
#
# PRIMARY entry point (used by run_celeste.R Stage 9's "ndvi" LOI):
# prepare_ndvi_per_group_rasters() — materializes ONE physical GeoTIFF per
# group, mosaicked from ONLY that group's own source tiles (via
# match_group_tiles() in workflow/raster_attributes.R), not a crop of a
# shared multi-province VRT. This is the same fix already applied to the
# NDVI trend LOI (prepare_ndvi_trend.R) and to harvest_regen
# (prepare_harvest_regen.R): per-site hydroweight processing crops+
# reprojects whichever raster a LOI's path_lazy points at
# (resolve_site_loi_raster() in workflow/R/stream/06_hydroweight_
# attributes.R), and a VRT built from many source tiles across a huge
# combined bounding box makes GDAL consider far more candidate source data
# per crop() than a small per-group file does — confirmed empirically for
# the trend LOI (~3x slower cropping the full VRT vs. a group's own tiles)
# and structurally consistent with why harvest_regen (one file per group)
# was already much faster per-site than this LOI's old VRT-based approach.
#
# prepare_ndvi_vrt() (below) is kept as a SECONDARY, ad hoc tool — a
# lightweight whole-project mosaic, handy for manual QGIS inspection across
# every group at once — but Stage 9 no longer uses it by default.
#
# Mosaicking (either function) sets vrtnodata so areas outside every tile
# read as NA, not 0 — without that fix a gap reads as valid-looking zeros
# indistinguishable from real data. The tiles' combined bounding box spans
# three widely-separated provinces (Ontario: KEN/NIP/TUR; New Brunswick:
# NBE/COC; PEI: MOR — confirmed by reprojecting each tile's extent to
# lon/lat) even though actual coverage per group is a small clustered
# patch — exactly why prepare_ndvi_per_group_rasters() materializes per
# group rather than caching that whole bounding box as one reprojected
# file (which would be both enormous and, per the above, slower per-site
# to boot).
#
# NBE now has its own 2 tiles (added after the initial NDVI batch) — no
# longer a coverage gap.
#
# Usage (from run_celeste.R, after sourcing workflow/raster_attributes.R):
#   source(here("workflow/CELESTE/prepare_ndvi.R"))
#   prepare_ndvi_per_group_rasters(group_manifest, cache_dir = cache_dir)
#
# Dependencies: terra, fs (via utils.R); match_group_tiles()/
#   build_mosaic_vrt() from workflow/raster_attributes.R must be sourced
#   first.
# =============================================================================

#' Clean one raw NDVI export tile: replace the 0 masked/no-data fill value
#' with true NA, and rescale from Earth Engine's int16 x10000 encoding back
#' to true NDVI (-1 to 1). See the header above for why both are needed and
#' why 0 is safe to treat as missing rather than a genuine NDVI reading.
#'
#' Cached per tile (by original filename) under cache_dir/hydroweight_loi/
#' ndvi_clean/ — cheap to skip on rerun, unlike the VRT (see
#' prepare_ndvi_vrt()'s docstring for why THAT stays uncached).
#'
#' @param ndvi_dir  Character. Directory of raw source NDVI tiles. Default
#'   here("data/ndvi").
#' @param cache_dir Character. Project cache root.
#'
#' @return Character vector of cleaned tile paths, one per raw tile in
#'   ndvi_dir, invisibly.
clean_ndvi_tiles <- function(ndvi_dir = here::here("data/ndvi"), cache_dir) {
  out_dir <- fs::path(cache_dir, "hydroweight_loi", "ndvi_clean")
  fs::dir_create(out_dir)

  raw_files <- list.files(ndvi_dir, pattern = "[.]tif$", full.names = TRUE)

  vapply(
    raw_files,
    function(f) {
      out_path <- fs::path(out_dir, basename(f))
      if (!cache_exists(out_path)) {
        r <- terra::rast(f)
        band_names <- names(r)
        r <- terra::subst(r, from = 0, to = NA) # EE masked-pixel export fill, not a real NDVI of 0
        r <- r / 10000 # EE export encoding: int16, scaled by 10000
        names(r) <- band_names # subst()/arithmetic can rename layers; restore the originals
        terra::writeRaster(r, out_path, overwrite = TRUE, datatype = "FLT4S")
      }
      out_path
    },
    character(1),
    USE.NAMES = FALSE
  )
}

#' Materialize one multi-band NDVI GeoTIFF per group, from ONLY that
#' group's own source tiles — the PRIMARY NDVI LOI prep function
#'
#' Reuses clean_ndvi_tiles() (0 -> NA masked-pixel fix + /10000 rescale)
#' and match_group_tiles() (workflow/raster_attributes.R) to mosaic just
#' one group's own tiles, then — unlike prepare_ndvi_vrt() — WRITES the
#' mosaic to a real per-group GeoTIFF rather than leaving it a VRT, cached
#' like prepare_harvest_regen_rasters()/prepare_ndvi_trend_rasters(). A
#' materialized single file avoids per-site crop() having to resolve which
#' of several candidate source tiles cover a given window every time —
#' genuinely worth the one-time write cost here (unlike a VRT), since it's
#' read from repeatedly (once per site).
#'
#' @param group_manifest sf tibble from build_group_manifest() — only
#'   group_id is used.
#' @param ndvi_dir  Character. Directory of raw source NDVI tiles. Default
#'   here("data/ndvi").
#' @param cache_dir Character. Project cache root — rasters are written to
#'   cache_dir/hydroweight_loi/ndvi/<group_id>.tif (the same directory
#'   resolve_site_loi_raster() already uses for this LOI's per-site
#'   "<site_id>_reprojected.tif" cache — no collision, per-site files
#'   always carry that suffix).
#'
#' @return Character vector of the (existing-or-newly-written) per-group
#'   raster paths, invisibly. A group with no matching tiles is skipped
#'   with a warning (same pattern as harvest_regen's/ndvi_trend's coverage
#'   gaps — the hydroweight module's "all NA after crop/mask" check
#'   excludes it downstream rather than erroring).
prepare_ndvi_per_group_rasters <- function(
  group_manifest,
  ndvi_dir = here::here("data/ndvi"),
  cache_dir
) {
  out_dir <- fs::path(cache_dir, "hydroweight_loi", "ndvi")
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
        "No NDVI tiles matched for group '{grp}'; skipping per-group NDVI ",
        "raster (will read as NA downstream, same as harvest_regen's MOR gap)."
      ))
      next
    }

    mini_vrt <- terra::vrt(
      matched,
      filename = tempfile(fileext = ".vrt"),
      options = c("-vrtnodata", "nan"),
      overwrite = TRUE
    )
    # terra::vrt() defaults each band's name from the VRT's OWN file
    # basename (here, a random tempfile() name) — confirmed directly: band
    # names came out as "file68e4e8de083_1" etc, different every run.
    # Since hydroweight_attributes()'s output column names embed the raw
    # raster's own band names (see run_loi_attributes_stream_multilayer()'s
    # docstring), that non-determinism would both break reproducibility and
    # produce column names inconsistent with every already-computed site
    # (which used prepare_ndvi_vrt()'s "ndvi_mosaic_N" convention). Pin it
    # explicitly to the same convention.
    names(mini_vrt) <- paste0("ndvi_mosaic_", seq_len(terra::nlyr(mini_vrt)))

    cw_inform(glue::glue(
      "Materializing NDVI raster for group '{grp}' ({length(matched)} tile(s), ",
      "{terra::ncell(mini_vrt)} cells)..."
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

#' Mosaic the (cleaned) NDVI tiles into one whole-project VRT — SECONDARY,
#' ad hoc tool, not used by Stage 9 by default
#'
#' Kept for manual QA (e.g. opening the whole-catchment NDVI mosaic in
#' QGIS across every group at once, something a per-group file can't do as
#' conveniently) — see prepare_ndvi_per_group_rasters() for the function
#' Stage 9 actually calls, and why.
#'
#' Always rebuilt (not cache_exists()-gated like
#' prepare_ndvi_per_group_rasters()'s per-group rasters) — a VRT is a
#' lightweight XML pointer, not materialized data, so rebuilding is cheap,
#' and skipping would mean a newly-added tile in ndvi_dir (this has already
#' happened once: two NBE tiles added after the initial NDVI batch) doesn't
#' get picked up until something manually deletes the stale VRT. The
#' underlying cleaned tiles ARE individually cached (see
#' clean_ndvi_tiles()), so this stays cheap even though it's unconditional.
#'
#' @param ndvi_dir  Character. Directory of raw source NDVI tiles. Default
#'   here("data/ndvi").
#' @param cache_dir Character. Project cache root — VRT is written to
#'   cache_dir/hydroweight_loi/ndvi_mosaic.vrt.
#'
#' @return Character path to the mosaic VRT, invisibly.
prepare_ndvi_vrt <- function(ndvi_dir = here::here("data/ndvi"), cache_dir) {
  vrt_path <- fs::path(cache_dir, "hydroweight_loi", "ndvi_mosaic.vrt")
  fs::dir_create(fs::path_dir(vrt_path))

  clean_files <- clean_ndvi_tiles(ndvi_dir = ndvi_dir, cache_dir = cache_dir)

  build_mosaic_vrt(
    files = clean_files,
    vrt_path = vrt_path
  )

  invisible(vrt_path)
}
