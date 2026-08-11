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
#   build_mosaic_vrt()       — mosaic scattered raster tiles (e.g. per-region
#     NDVI exports) into one lightweight VRT, safe to pass as a hydroweight
#     loi_layers `path_lazy` (see workflow/R/stream/06_hydroweight_attributes.R)
#   rasterize_competing_classes() — rasterize several "class bucket" vector
#     layers (e.g. harvest vs. regen disturbance polygons) into one
#     categorical raster with temporal precedence between classes, one band
#     per year plus one all-years-combined band
#
# Dependencies: terra, trend, sf
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

# -- Mosaic VRT construction ---------------------------------------------------

#' Check that a set of raster tiles are safe to mosaic together
#'
#' A VRT does not reproject or reconcile mismatched inputs — it just points
#' at pixel windows across files, assuming they share a common grid/CRS and
#' (for multi-layer tiles) band structure. This reports per-file CRS,
#' resolution, layer count, and first/last layer names so mismatches (e.g.
#' one tile in a different CRS, or a shorter time series) are caught before
#' build_mosaic_vrt(), not after.
#'
#' @param files Character vector of raster tile file paths.
#' @return A tibble, one row per file, with columns: file, crs, res,
#'   nlyr, first_layer, last_layer, and a logical `consistent` column
#'   flagging rows that differ from the first file's crs/res/nlyr.
check_tile_consistency <- function(files) {
  info <- lapply(files, function(f) {
    r <- terra::rast(f)
    list(
      file = basename(f),
      crs = terra::crs(r, describe = TRUE)$name,
      res = terra::res(r)[1],
      nlyr = terra::nlyr(r),
      first_layer = names(r)[1],
      last_layer = names(r)[terra::nlyr(r)]
    )
  })
  tbl <- do.call(rbind.data.frame, c(info, stringsAsFactors = FALSE))
  tbl$consistent <- tbl$crs == tbl$crs[1] & tbl$res == tbl$res[1] & tbl$nlyr == tbl$nlyr[1]
  tibble::as_tibble(tbl)
}

#' Mosaic scattered raster tiles into one VRT, with a nodata safety fix
#'
#' Thin wrapper around terra::vrt() for the common case of several
#' non-overlapping (or partially overlapping) raster tiles — e.g. NDVI time-
#' series exports done per region/group rather than one file covering the
#' full site extent — that need to be readable as a single raster.
#'
#' By default, terra::vrt()/gdalbuildvrt leaves any area NOT covered by any
#' input tile at the format's fill value (typically 0), NOT NoData. This is
#' a real correctness hazard: cropping/masking that area for a site whose
#' catchment falls outside every tile's coverage returns valid-looking
#' zeros indistinguishable from real data, rather than failing or returning
#' NA (confirmed directly: a VRT built from CELESTE's data/ndvi/*.tif tiles,
#' which do not cover the NBE group at all, returned exactly 0 for every
#' NBE catchment until this fix was applied). This function always passes
#' `-vrtnodata` so uncovered areas correctly read as NA, letting
#' the existing "all NA after crop/mask" checks in the hydroweight modules
#' catch and skip them instead of silently computing wrong results.
#'
#' Does NOT reproject or align band/year structure across input tiles — all
#' inputs must already share a common CRS, resolution, and (for multi-layer
#' tiles) identical band count/order. Verify this first (e.g. loop
#' terra::rast() over the files and compare crs()/res()/nlyr()/names()); a
#' VRT cannot correct for mismatches between its source files.
#'
#' @param files    Character vector of raster tile file paths.
#' @param vrt_path Character. Output .vrt file path.
#' @param nodata   Value to use for both the vrtnodata safety fix and (if not
#'   already set on the sources) srcnodata. Default "nan", matching typical
#'   float raster NoData. Use a numeric sentinel instead if your tiles are
#'   an integer type without a NaN representation.
#' @param overwrite Logical. Overwrite vrt_path if it exists. Default TRUE.
#'
#' @return SpatRaster (the VRT, opened) — the mosaic will not be
#'   materialized to disk except at vrt_path (a lightweight XML pointer);
#'   later crop()/project() calls only read the pixels they need from the
#'   underlying tiles.
build_mosaic_vrt <- function(files, vrt_path, nodata = "nan", overwrite = TRUE) {
  missing_files <- files[!file.exists(files)]
  if (length(missing_files) > 0) {
    stop(sprintf(
      "File(s) not found: %s",
      paste(missing_files, collapse = ", ")
    ))
  }

  terra::vrt(
    files,
    filename = vrt_path,
    options = c("-vrtnodata", as.character(nodata)),
    overwrite = overwrite
  )
}

# -- Competing-class disturbance rasterization ---------------------------------

#' Rasterize competing "class bucket" vector layers with temporal precedence
#'
#' Built for cases like harvest vs. regeneration disturbance tracking: one or
#' more vector layers represent each of several classes (e.g. every harvest-
#' method layer -> "harvest", every regen-method layer -> "regen"), each
#' feature tagged with a year, and the same location can appear in more than
#' one class's source layers across time — most recent year wins. Produces
#' one band per requested year plus one all-years-combined band, all using
#' the identical precedence rule: a location's class is whichever bucket's
#' MOST RECENT covering feature is latest; a per-year band only ever sees
#' same-year features from every bucket, so ties there are the norm, not the
#' exception (`buckets` order settles them — see below).
#'
#' Geometry hygiene: features are cast to MULTIPOLYGON (fixes mixed
#' MULTIPOLYGON/MULTISURFACE layers — a real issue found in
#' ontario_harvest.gdb's newer-vintage layers, which terra::vect() cannot
#' read directly) and passed through st_make_valid() before rasterizing.
#'
#' @param buckets    Named list, IN PRIORITY ORDER — when two buckets' most
#'   recent year is exactly equal at a cell (always true for a per-year
#'   band, since every included feature already shares that year), the
#'   LATER-named bucket in this list wins. Each name is a class label; each
#'   value is a character vector of gdb layer names (resolved via
#'   `gdb_path`) or file paths readable by sf::st_read(), unioned together
#'   for that class.
#' @param gdb_path   Character. Geodatabase path for bucket entries that are
#'   bare layer names. NULL if every bucket entry is already a full file
#'   path.
#' @param year_field Character. Field holding each feature's year. Default
#'   "AR_YEAR".
#' @param template   SpatRaster defining the output grid (crs/extent/res) —
#'   e.g. a group's DEM.
#' @param years      Integer vector of years to produce per-year bands for.
#'   Default: every year present across all buckets' `year_field` values
#'   within `crop_to` and `year_range` (if supplied) — i.e. computed from
#'   the data itself.
#' @param year_range Optional length-2 integer vector `c(min, max)`.
#'   Features outside this range are dropped BEFORE any rasterization —
#'   affects both the per-year bands and "combined" (unlike `years`, which
#'   only controls which per-year bands get produced but leaves "combined"
#'   using every year present). Use this to make two sources with
#'   different native date ranges comparable — e.g. one region's record
#'   going back to 1961 and another's only to 2002; without this, the
#'   longer-record region's "combined" band reflects a much longer
#'   disturbance history and isn't a fair comparison. Default NULL (no
#'   filtering).
#' @param crop_to    sf/SpatVector (any CRS) to spatially filter source
#'   layers to — e.g. a group AOI. Recommended for province/national-scale
#'   sources so the rest of the pipeline only deals with the relevant
#'   subset. Buffered by crop_buffer_m. NOT pushed down to OGR at read
#'   time (deliberately — see below); every layer is read in full, then
#'   filtered in R after geometry cleanup. Also used to mask the OUTPUT:
#'   cells outside crop_to are set to NA (not the "other"/0 background),
#'   since we never evaluated whether they're disturbed or not — 0 would
#'   misrepresent "not checked" as "confirmed absent". Confirmed this
#'   distinction matters in practice: `template`'s extent is a rectangular
#'   DEM crop, usually noticeably larger than the (irregular) AOI actually
#'   passed as crop_to, so without this mask a real disturbance polygon
#'   just outside the AOI — never touched by rasterization at all — reads
#'   as indistinguishable from a location confirmed to have no harvest or
#'   regen.
#' @param crop_buffer_m Numeric. Buffer applied to crop_to (after
#'   transforming it to `template`'s CRS), in that CRS's units (metres for
#'   every CRS used in this project). Default 1000.
#'
#' @return SpatRaster: one layer per year in `years` (named "y<year>") plus
#'   one "combined" layer, integer-coded — 0 = none of the buckets ("other"),
#'   1 = first bucket, 2 = second bucket, etc. (matches `seq_along(buckets)`
#'   in the order `buckets` was given, NOT priority order); NA outside
#'   crop_to (see above), if supplied. Pair with `class_levels =
#'   data.frame(ID = 0:length(buckets), Class = c("other", names(buckets)))`
#'   when passing to the hydroweight stage.
rasterize_competing_classes <- function(
  buckets,
  gdb_path = NULL,
  year_field = "AR_YEAR",
  template,
  years = NULL,
  year_range = NULL,
  crop_to = NULL,
  crop_buffer_m = 1000
) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required.")
  }
  if (is.null(names(buckets)) || any(!nzchar(names(buckets)))) {
    stop("buckets must be a named list (one name per class).")
  }

  # Computed once, in template's CRS — reused both to filter input features
  # (in read_one()) and to mask the final output (see combine_bucket_years()
  # / crop_to's docstring for why masking, not just filtering, matters).
  crop_geom <- if (!is.null(crop_to)) {
    crop_to |>
      sf::st_transform(terra::crs(template)) |>
      sf::st_buffer(crop_buffer_m)
  } else {
    NULL
  }

  read_one <- function(layer_ref) {
    is_bare_layer <- is.null(gdb_path) == FALSE && !grepl("[.]", basename(layer_ref))
    src <- if (is_bare_layer) gdb_path else layer_ref
    lyr <- if (is_bare_layer) layer_ref else NULL

    # NOT filtered at read time (wkt_filter / -spat), even though that
    # would be faster for a province-scale source — deliberately, after a
    # real, confirmed failure: wkt_filter is pushed down to OGR BEFORE any
    # R-side geometry cleanup, and curved geometries (MULTISURFACE — the
    # same issue rasterize_competing_classes() already casts away below)
    # make that push-down unreliable. Confirmed on two real datasets, two
    # different severities: Ontario's harvest gdb silently dropped some
    # (not all) curved features near a filter boundary (a handful out of
    # ~4100 province-wide); New Brunswick's gdb returned ZERO features for
    # every crop_to tried, including a 100km buffer around a real site
    # coordinate confirmed well inside the layer's own reported extent.
    # Reading in full and filtering in R afterward, once geometries are
    # already cleaned up, is slower but correct regardless of the source's
    # curve usage.
    # sf::st_read(src, layer = NULL, ...) is NOT the same as omitting
    # `layer` — sf calls enc2utf8(layer) unconditionally internally, which
    # errors on NULL ("argument is not a character vector"). Confirmed
    # directly. Must branch, not just pass lyr through when it's NULL.
    x <- if (is.null(lyr)) {
      sf::st_read(src, quiet = TRUE)
    } else {
      sf::st_read(src, layer = lyr, quiet = TRUE)
    }

    if (nrow(x) == 0) {
      return(x)
    }
    if (!year_field %in% names(x)) {
      stop(sprintf("Layer '%s' has no field '%s'.", layer_ref, year_field))
    }

    # Drop empty geometries before casting — a mix of MULTISURFACE (curve)
    # and even one empty GEOMETRYCOLLECTION feature in the same column
    # makes st_cast("MULTIPOLYGON") fail outright for the WHOLE column
    # ("MULTISURFACE can only be converted into MULTIPOLYGON"), not just
    # skip the bad feature. Confirmed directly against real data: New
    # Brunswick's HarvHistory layer has 9 such empty-geometry records.
    x <- x[!sf::st_is_empty(x), ]

    x <- x |>
      dplyr::select(dplyr::all_of(year_field)) |>
      sf::st_cast("MULTIPOLYGON") |>
      sf::st_make_valid() |>
      sf::st_transform(terra::crs(template)) |>
      dplyr::filter(!is.na(.data[[year_field]]))

    if (!is.null(crop_geom) && nrow(x) > 0) {
      x <- sf::st_filter(x, crop_geom)
    }

    if (!is.null(year_range) && nrow(x) > 0) {
      x <- x[x[[year_field]] >= year_range[1] & x[[year_field]] <= year_range[2], ]
    }

    x
  }

  read_bucket <- function(layer_refs) {
    parts <- lapply(layer_refs, read_one)
    parts <- parts[vapply(parts, nrow, integer(1)) > 0]
    if (length(parts) == 0) {
      return(NULL)
    }
    dplyr::bind_rows(parts)
  }

  cw_or_message <- function(msg) {
    if (exists("cw_inform", mode = "function")) cw_inform(msg) else message(msg)
  }

  cw_or_message("Reading and cleaning bucket layers...")
  bucket_sf <- lapply(buckets, read_bucket)

  if (is.null(years)) {
    all_years <- unlist(lapply(bucket_sf, function(x) {
      if (is.null(x)) NULL else x[[year_field]]
    }))
    years <- sort(unique(all_years))
    if (length(years) == 0) {
      stop("No features found in any bucket (after crop_to filtering) — nothing to rasterize.")
    }
  }

  rasterize_year_field <- function(x) {
    if (is.null(x) || nrow(x) == 0) {
      return(NULL)
    }
    terra::rasterize(
      terra::vect(x),
      template,
      field = year_field,
      fun = "max",
      background = NA
    )
  }

  # Boolean vector, TRUE for cells inside crop_geom — computed once (a
  # single rasterize() call), applied to every band in combine_bucket_
  # years() below. Kept as a plain vector, not a SpatRaster, for the same
  # in-memory-vector-arithmetic reason as combine_bucket_years() itself.
  inside_mask <- if (is.null(crop_geom)) {
    NULL
  } else {
    m <- terra::rasterize(terra::vect(crop_geom), template, field = 1, background = NA)
    !is.na(terra::values(m)[, 1])
  }

  # Operates on plain in-memory numeric vectors (terra::values()), not
  # SpatRaster ifel() — each ifel() call is disk-backed and, measured
  # directly against this, ~4x slower per bucket than the vector
  # equivalent (confirmed: cut a 12-group real-data run from ~200s to well
  # under a minute for TUR). Writes the result straight to out_path and
  # returns only the path — critical for large groups, not just tidy:
  # KEN's grid is 139.5M cells (3.5x TUR's); the original version of this
  # function returned a SpatRaster and the caller accumulated all 24
  # (23 years + combined) of them in a list before combining at the end,
  # which meant ~24 in-memory full-size rasters alive simultaneously and
  # OOM-killed a real run on KEN (confirmed: exit 137). Writing
  # immediately and keeping only the (tiny) file path means at most one
  # band's worth of full-size vectors is ever alive at a time.
  combine_bucket_years <- function(year_rasters, out_path) {
    n_cell <- terra::ncell(template)
    best_year <- rep(NA_real_, n_cell)
    best_class <- rep(0L, n_cell)
    for (i in seq_along(year_rasters)) {
      yr <- year_rasters[[i]]
      if (is.null(yr)) next
      v <- terra::values(yr)[, 1]
      take <- !is.na(v) & (is.na(best_year) | v >= best_year)
      best_year[take] <- v[take]
      best_class[take] <- i
    }
    # NA (not 0/"other") outside crop_to — see crop_to's docstring for why
    # this distinction is real, not cosmetic: a cell here was never
    # evaluated, so 0 would misrepresent "not checked" as "confirmed
    # absent" (confirmed in practice: a real disturbance polygon just
    # outside a group's AOI, inside its rectangular DEM extent, otherwise
    # reads identically to a location with genuinely no harvest/regen).
    if (!is.null(inside_mask)) {
      best_class[!inside_mask] <- NA_integer_
    }
    out <- template[[1]]
    terra::values(out) <- best_class
    terra::writeRaster(out, out_path, overwrite = TRUE, datatype = "INT1U")
    invisible(out_path)
  }

  cw_or_message(glue::glue("Rasterizing {length(years)} year(s) + combined..."))

  # NOT cleaned up on exit — the returned raster is a LAZY reference to
  # these files (terra::rast() on paths doesn't materialize data), so
  # they must survive until the caller actually reads from it (e.g. their
  # own writeRaster() call). Left for the R session's normal temp-
  # directory cleanup, same as terra's own internal temp files.
  tmp_dir <- tempfile("rasterize_competing_classes_")
  dir.create(tmp_dir)

  year_paths <- purrr::map_chr(years, function(yr) {
    year_bucket_rasters <- lapply(bucket_sf, function(x) {
      if (is.null(x)) NULL else rasterize_year_field(x[x[[year_field]] == yr, ])
    })
    out_path <- combine_bucket_years(
      year_bucket_rasters,
      fs::path(tmp_dir, glue::glue("y{yr}.tif"))
    )
    rm(year_bucket_rasters)
    gc(FALSE)
    out_path
  })

  all_bucket_rasters <- lapply(bucket_sf, rasterize_year_field)
  combined_path <- combine_bucket_years(
    all_bucket_rasters,
    fs::path(tmp_dir, "combined.tif")
  )
  rm(all_bucket_rasters)
  gc(FALSE)

  # Lazy read — terra::rast() on file paths doesn't materialize the data,
  # only the caller's eventual writeRaster()/values() call does, and that
  # happens in terra's own chunked fashion.
  out <- terra::rast(c(year_paths, combined_path))
  names(out) <- c(paste0("y", years), "combined")
  out
}

# -- Null coalescing operator (mirrors standard R 4.4+ behaviour) ----
`%||%` <- function(x, y) if (!is.null(x)) x else y
