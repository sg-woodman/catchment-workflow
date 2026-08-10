# prepare_harvest_regen.R
# =============================================================================
# Prepares the harvest/regen disturbance LOI raster for the CELESTE stream
# hydroweight stage (workflow/run_celeste.R, Stage 9). Kept as its own
# script rather than inline in the runner — it's a distinct, sometimes slow
# (minutes per group, dominated by reading full source layers + rasterizing
# over each group's full DEM extent) one-time data-prep step, not part of
# the interactive delineation stage sequence.
#
# One shared LOI ("harvest_regen"), two independent source datasets — Note
# some sites fall in the Ontario tenure, others in the Irving/NB tenure, and
# the two forestry agencies track harvest/regen completely separately (own
# gdb, own field names, own CRS, own data-quality quirks). Both get
# rasterized into the SAME class scheme (harvest_regen_levels below) so
# they combine into one LOI across every group, per
# cache/CELESTE/hydroweight_loi/harvest_regen/<group_id>.tif — run_celeste.R
# doesn't need to know or care which source a given group came from.
#
#   Ontario (shared_data/raw/harvest/ontario_harvest.gdb): NIP, TUR, KEN.
#     AR_YEAR field directly. CC-only harvest scope (Harvest_CC02 +
#     Harvest_CC17 = the full 2002-2024 clear-cut record, two non-
#     overlapping vintage exports of the same product).
#   New Brunswick / Irving (data/irving_harvest/
#     LB_HarvCuHi_RefCuHi_SICuHi.gdb): NBE. HARVYR (harvest) / RFYR
#     (regen) fields, staged (cleaned + renamed to AR_YEAR) before
#     rasterizing — see stage_nb_harvest_layer() for why. CC-only harvest
#     scope too (HARVTRT == "CC"), for consistency with Ontario.
#
# To extend either source's harvest scope (e.g. add partial-cut / selection
# methods), edit that source's `buckets`/`harv_trt` below — no changes
# needed elsewhere. To add a third source/region, add another block
# following the same pattern.
#
# Both sources: two competing classes (harvest, regen), most-recent-year-
# wins where they'd otherwise overlap. Ties — which is EVERY within-year
# overlap, since a per-year band only ever includes same-year features
# from both classes by construction — go to whichever bucket is named LAST
# in `buckets`; regen is last, so regen wins ties, matching "harvested
# then regenerated -> regen" as the more ecologically current state.
#
# Usage (from run_celeste.R, after sourcing workflow/raster_attributes.R):
#   source(here("workflow/prepare_harvest_regen.R"))
#   prepare_harvest_regen_rasters(group_manifest, cache_dir = cache_dir)
#
# Dependencies: terra, sf, dplyr, fs (via utils.R); rasterize_competing_
#   classes() from workflow/raster_attributes.R must be sourced first.
# =============================================================================

#' Class levels for the harvest/regen categorical layer (shared by every
#' source — this is what makes them combinable into one LOI)
harvest_regen_levels <- data.frame(ID = 0:2, Class = c("other", "harvest", "regen"))

#' Clean, filter, and rename-to-AR_YEAR one Irving/NB harvest or regen
#' layer, caching the result as a small GPKG for rasterize_competing_
#' classes() to read directly (a plain file path, same as any other
#' bucket entry).
#'
#' Two real data-quality issues found by inspecting the source, neither of
#' which apply to the Ontario source:
#'
#'   1. `UNTREATI` (untreated indicator, documented in the accompanying
#'      data dictionary): non-zero means this polygon is a mapped sub-area
#'      WITHIN a treated block that was NOT actually harvested/
#'      regenerated (unmapped water buffers, wildlife buffers, wet areas,
#'      steep/rocky ground, immature/low-volume stands not actually cut,
#'      etc.) — these are not disturbance events and must be excluded, not
#'      just filtered by a plausible year (see #2 — many of them carry a
#'      literal sentinel year specifically because they weren't treated).
#'      Confirmed directly: EVERY record with year == 2099 has UNTREATI
#'      != 0; the reverse doesn't quite hold (a small number of untreated-
#'      indicator-0 records still carry a garbage year — #2 catches those).
#'
#'   2. A residual few records have an implausible/placeholder year (0,
#'      100, 1368, 2099, ...) even after the UNTREATI filter above — can't
#'      be assigned to any real year band, so excluded via a plausible-
#'      range check (min_year to max_year).
#'
#' @param gdb_path   Character. Path to the Irving/NB gdb.
#' @param layer      Character. Layer name (e.g. "HarvCurrent").
#' @param year_field Character. The layer's own year field name (e.g.
#'   "HARVYR" or "RFYR") — renamed to "AR_YEAR" in the output so it can
#'   share rasterize_competing_classes()'s year_field with the Ontario
#'   buckets.
#' @param harv_trt   Optional character vector. If supplied, keeps only
#'   features whose HARVTRT is in this set (e.g. "CC" for clear-cut-only,
#'   matching the Ontario harvest scope). NULL (default) keeps everything
#'   — used for the regen layers, which have no equivalent scope
#'   restriction on the Ontario side either.
#' @param cache_path Character. Output GPKG path. Skipped (returned as-is)
#'   if it already exists.
#' @param min_year,max_year Integer. Plausible year bounds for the
#'   residual-bad-year backstop. Defaults 1950 to next calendar year.
#'
#' @return cache_path, invisibly.
stage_nb_harvest_layer <- function(
  gdb_path,
  layer,
  year_field,
  harv_trt = NULL,
  cache_path,
  min_year = 1950,
  max_year = as.integer(format(Sys.Date(), "%Y")) + 1
) {
  if (cache_exists(cache_path)) {
    return(invisible(cache_path))
  }
  fs::dir_create(fs::path_dir(cache_path))

  x <- sf::st_read(gdb_path, layer = layer, quiet = TRUE)

  keep <- x[["UNTREATI"]] == 0 &
    !is.na(x[[year_field]]) &
    x[[year_field]] >= min_year &
    x[[year_field]] <= max_year
  if (!is.null(harv_trt)) {
    keep <- keep & x[["HARVTRT"]] %in% harv_trt
  }
  x <- x[keep, year_field]
  names(x)[names(x) == year_field] <- "AR_YEAR"

  sf::st_write(x, cache_path, delete_dsn = TRUE, quiet = TRUE)
  invisible(cache_path)
}

#' Stage all four Irving/NB layers (2 harvest, 2 regen) into cleaned GPKGs
#'
#' @param gdb_path  Character. Path to the Irving/NB gdb.
#' @param cache_dir Character. Where staged GPKGs are written
#'   (cache_dir/hydroweight_loi/harvest_regen/nb_staged/).
#' @return Named list with $harvest and $regen character vectors of staged
#'   file paths, ready to pass as rasterize_competing_classes() buckets.
stage_nb_harvest_regen <- function(gdb_path, cache_dir) {
  staged_dir <- fs::path(cache_dir, "hydroweight_loi", "harvest_regen", "nb_staged")

  harvest <- vapply(
    c("HarvCurrent", "HarvHistory"),
    function(lyr) {
      stage_nb_harvest_layer(
        gdb_path, lyr,
        year_field = "HARVYR",
        harv_trt = "CC", # CC-only, matching the Ontario harvest scope
        cache_path = fs::path(staged_dir, paste0(lyr, "_cc.gpkg"))
      )
    },
    character(1)
  )

  regen <- vapply(
    c("ReforCurrent", "ReforHistory"),
    function(lyr) {
      stage_nb_harvest_layer(
        gdb_path, lyr,
        year_field = "RFYR",
        cache_path = fs::path(staged_dir, paste0(lyr, ".gpkg"))
      )
    },
    character(1)
  )

  list(harvest = unname(harvest), regen = unname(regen))
}

#' Rasterize the harvest/regen disturbance layer, once per group, from
#' whichever source (Ontario or NB) that group belongs to
#'
#' Skips any group whose output raster already exists (cache_exists()) —
#' safe to call on every run_celeste.R source.
#'
#' @param group_manifest sf tibble from build_group_manifest() — used to
#'   resolve each group's AOI (for spatially filtering source layers) and
#'   cache directory (for the group's DEM template).
#' @param cache_dir      Character. Project cache root — rasters are written
#'   to cache_dir/hydroweight_loi/harvest_regen/<group_id>.tif, and each
#'   group's DEM template is read from cache_dir/<group_id>/dem_breached.tif.
#' @param sources        Named list of source configs, each with $groups
#'   (character vector of group_ids this source covers), $gdb_path, and
#'   $buckets (passed to rasterize_competing_classes(); NB's are resolved
#'   from stage_nb_harvest_regen() by default, not hardcoded, since they're
#'   cache paths that depend on cache_dir). Default: Ontario (NIP/TUR/KEN)
#'   + NB (NBE) — the two sources confirmed by real spatial-filter feature
#'   counts (not just bounding-box overlap — that check gave false
#'   positives for both sources during development) to cover every CELESTE
#'   group with harvest/regen data available. COC/MOR sites get no
#'   coverage from either and are silently excluded downstream by the
#'   hydroweight module's existing "all NA after crop/mask" check, same
#'   pattern as the NDVI/NBE gap this fills.
#'
#' @return Character vector of the (existing-or-newly-written) per-group
#'   raster paths, invisibly.
prepare_harvest_regen_rasters <- function(
  group_manifest,
  cache_dir,
  sources = list(
    ontario = list(
      groups = c("NIP", "TUR", "KEN"),
      gdb_path = "~/Documents/cfs/shared_data/raw/harvest/ontario_harvest.gdb",
      buckets = list(
        harvest = c("Harvest_CC02", "Harvest_CC17"),
        regen   = c("Regen_Seed", "Regen_Natural", "Regen_Plant")
      )
    ),
    nb = list(
      groups = "NBE",
      gdb_path = here::here("data/irving_harvest/LB_HarvCuHi_RefCuHi_SICuHi.gdb"),
      buckets = NULL # resolved below via stage_nb_harvest_regen()
    )
  )
) {
  written <- character(0)

  for (src in sources) {
    gdb_path <- path.expand(src$gdb_path)
    buckets <- src$buckets %||% stage_nb_harvest_regen(gdb_path, cache_dir)

    for (grp in src$groups) {
      group_raster_path <- fs::path(
        cache_dir, "hydroweight_loi", "harvest_regen", paste0(grp, ".tif")
      )
      written <- c(written, group_raster_path)
      if (cache_exists(group_raster_path)) next

      fs::dir_create(fs::path_dir(group_raster_path))
      grp_aoi <- group_manifest[group_manifest$group_id == grp, ]
      grp_template <- terra::rast(fs::path(cache_dir, grp, "dem_breached.tif"))

      hr <- rasterize_competing_classes(
        buckets = buckets,
        gdb_path = if (is.null(src$buckets)) NULL else gdb_path,
        template = grp_template,
        crop_to = grp_aoi
      )
      terra::writeRaster(
        hr,
        group_raster_path,
        overwrite = TRUE,
        datatype = "INT1U",
        gdal = "PHOTOMETRIC=MINISBLACK" # avoid a real GDAL quirk: an
        # exactly-3-band Byte GeoTIFF gets auto-tagged as RGB, silently
        # renaming bands to red/green/blue on every subsequent read
      )
    }
  }

  invisible(written)
}
