# raster_attributes.R
# =============================================================================
# General-purpose raster preparation utilities, independent of both the
# stream and lake catchment-delineation pipelines (no output_dir/cache_dir/
# sites/group_manifest conventions — pure raster-in, raster-out functions).
# Source this ad hoc when preparing a LOI (layer-of-interest) raster before
# adding it to a loi_layers list for the hydroweight stage, e.g. to build a
# cleaned land-cover raster or an NDVI trend raster.
#
#   reclassify_categorical() — collapse/relabel raw categorical raster codes
#     into cleaner class groups via a from/to lookup table
#   sens_slope_trend()       — per-pixel Theil-Sen slope + Mann-Kendall
#     p-value across a raster time-series stack (e.g. annual NDVI composites)
#
# Dependencies: terra, trend
# =============================================================================

library(terra)

# -- Categorical raster reclassification --------------------------------------

#' Reclassify a categorical raster using a from -> to lookup table
#'
#' Wraps terra::classify() with a two-column ("is", "becomes") matrix, which
#' uses an optimized hash lookup for integer values — the right tool for
#' collapsing many raw class codes into fewer, cleaner groups (e.g. merging
#' several raw land-cover subclasses into one "wetland" class) or simply
#' relabeling codes. Unlike terra::subst(), this also lets you optionally
#' attach a level table to the *output* classes in one step.
#'
#' @param r             SpatRaster. Categorical raster (integer codes).
#' @param reclass_table data.frame with columns matching `from_col`/`to_col`.
#'   One row per raw code being reclassified; multiple rows may map to the
#'   same `to` value to collapse classes. Codes not listed are handled per
#'   `others`.
#' @param from_col      Character. Column name for raw values. Default "from".
#' @param to_col        Character. Column name for new values. Default "to".
#' @param others        "na" sets any code not listed in `reclass_table` to
#'   NA (default — safest, makes silently-unhandled codes visible as
#'   missing data rather than passing them through unchanged). "keep"
#'   leaves unmatched codes at their original value.
#' @param class_levels  Optional data.frame with ID/Class columns (same
#'   convention as the `loi_layers[[i]]$class_levels` field used by the
#'   hydroweight stage) to attach as factor levels on the *output* raster.
#'   Default NULL — output stays a plain integer raster.
#'
#' @return SpatRaster, reclassified (and, if class_levels supplied, releveled
#'   as a factor)
reclassify_categorical <- function(
  r,
  reclass_table,
  from_col = "from",
  to_col = "to",
  others = c("na", "keep"),
  class_levels = NULL
) {
  others <- match.arg(others)

  if (!all(c(from_col, to_col) %in% names(reclass_table))) {
    stop(sprintf(
      "reclass_table must have columns '%s' and '%s'.",
      from_col, to_col
    ))
  }

  rcl <- as.matrix(reclass_table[, c(from_col, to_col)])
  storage.mode(rcl) <- "numeric"

  out <- terra::classify(
    r,
    rcl,
    others = if (others == "na") NA else NULL
  )

  if (!is.null(class_levels)) {
    out <- terra::as.factor(out)
    levels(out) <- class_levels
  }

  out
}

# -- Sen's slope trend on a raster time series --------------------------------

#' Per-pixel Theil-Sen slope + Mann-Kendall significance across a raster stack
#'
#' Applies the Theil-Sen slope estimator and a Mann-Kendall trend test to
#' each pixel's time series (one value per raster layer), returning a
#' 2-layer raster: `slope` (change per unit of `x`, e.g. per year) and
#' `p_value` (Mann-Kendall two-sided test).
#'
#' No autocorrelation prewhitening (e.g. Yue-Pilon/Zhang) is applied. That
#' was evaluated and deliberately rejected: on short/noisy series (typical
#' of annual satellite composites) prewhitening's own detrending step can
#' conflate the true trend with the autocorrelation estimate and distort the
#' slope magnitude — demonstrated empirically during development (see PR/
#' commit history) with a synthetic AR(1)-noise + known-trend series where
#' Yue-Pilon underestimated the true slope 10-fold and Zhang recovered the
#' wrong sign entirely. Plain Theil-Sen + Mann-Kendall is the standard
#' default in environmental trend analysis generally and doesn't carry that
#' risk.
#'
#' The slope is computed directly from pairwise (y_j - y_i) / (x_j - x_i)
#' using the *actual* `x` values for every valid pair, not layer position —
#' this matters because trend::sens.slope() assumes evenly-spaced,
#' gap-free input and would silently mis-estimate the slope for any pixel
#' where per-pixel NAs (e.g. cloud-masked years) leave irregular gaps after
#' filtering. The Mann-Kendall p-value (via trend::mk.test()) is unaffected
#' by this — it's a purely rank/order-based test, insensitive to actual
#' spacing — so it's used as-is.
#'
#' @param r       SpatRaster. Multi-layer time series (e.g. one layer per
#'   year of NDVI composites). Layer order must match `x` order.
#' @param x       Numeric vector, same length as terra::nlyr(r). Time index
#'   (e.g. calendar years). Default: 1:nlyr(r) (equivalent to assuming
#'   regular, unit spacing).
#' @param min_obs Integer. Minimum non-NA layers required to fit a trend for
#'   a given pixel; below this, both outputs are NA. Default 4.
#' @param cores   Integer. Passed to terra::app() for parallel evaluation
#'   across pixels. Default 1.
#'
#' @return 2-layer SpatRaster named "slope" and "p_value"
sens_slope_trend <- function(r, x = NULL, min_obs = 4, cores = 1) {
  if (!requireNamespace("trend", quietly = TRUE)) {
    stop("Package 'trend' is required. Install with install.packages('trend').")
  }
  if (terra::nlyr(r) < min_obs) {
    stop(sprintf(
      "Raster has %d layer(s); min_obs = %d requires at least that many time steps.",
      terra::nlyr(r), min_obs
    ))
  }

  x <- x %||% seq_len(terra::nlyr(r))
  if (length(x) != terra::nlyr(r)) {
    stop("length(x) must equal terra::nlyr(r).")
  }

  fit_pixel <- function(y) {
    y <- as.numeric(y)
    ok <- !is.na(y)
    if (sum(ok) < min_obs) {
      return(c(slope = NA_real_, p_value = NA_real_))
    }

    xo <- x[ok]
    yo <- y[ok]
    n <- length(yo)

    pairs <- utils::combn(n, 2)
    dx <- xo[pairs[2, ]] - xo[pairs[1, ]]
    dy <- yo[pairs[2, ]] - yo[pairs[1, ]]
    slopes <- dy / dx
    slope <- stats::median(slopes[is.finite(slopes)], na.rm = TRUE)

    p_value <- tryCatch(
      trend::mk.test(yo)$p.value,
      error = function(e) NA_real_
    )
    if (is.nan(p_value)) p_value <- NA_real_

    c(slope = slope, p_value = p_value)
  }

  out <- terra::app(r, fit_pixel, cores = cores)
  names(out) <- c("slope", "p_value")
  out
}

# -- Null coalescing operator (mirrors standard R 4.4+ behaviour) ----
`%||%` <- function(x, y) if (!is.null(x)) x else y
