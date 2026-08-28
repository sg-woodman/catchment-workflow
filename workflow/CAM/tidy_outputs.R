# tidy_outputs.R
# =============================================================================
# Reshapes run_cam_streams.R's two wide deliverables — catchment_metrics.csv
# and CAM_streams_hydroweight.csv — into a handful of small, purpose-shaped
# long tables for plotting/analysis, rather than one wide table per LOI x
# scheme x stat combination.
#
# DESIGN, v2 (2026-08-26): an earlier version of this script produced ONE
# combined hydroweight long table (site, version, loi, year, class, variable,
# scheme, stat, value) — technically long, but not actually plot-ready: a
# single `value` column mixed proportions (0-1), NDVI means (~-1 to 1), trend
# slopes, and p-values, so every plot still needed a filter-then-pivot first.
# Replaced with ONE FILE PER LOI, each shaped for that LOI's natural plot
# type — `mean`/`sd` etc. are real columns, not melted into a generic
# stat/value pair. See each tidy_hydroweight_<loi>() docstring for its exact
# shape and row count.
#
# The ORIGINAL wide CSVs (catchment_metrics.csv, CAM_streams_hydroweight.csv)
# are untouched by this script and remain the right choice for statistical
# modelling (regression, PCA, ...) — one row per site x version, one column
# per covariate. This script only adds plotting-shaped long tables alongside
# them; it doesn't replace them.
#
# MERGING UNWEIGHTED + DISTANCE-WEIGHTED STATS INTO ONE COLUMN: every
# continuous LOI (ndvi, ndvi_trend) computes a flat/unweighted mean+sd
# (scheme "lumped") AND a distance-weighted mean+sd per weighted scheme
# (`<scheme>_distwtd_mean/sd`). These are the same statistical concept
# (the LOI's central estimate under a given weighting scheme) computed two
# different ways — `scheme` already tells you which, so both live in one
# `mean`/`sd` column pair rather than two column-pairs where one is always
# NA depending on the row's scheme.
#
# YEAR HANDLING PER TABLE (unchanged reasoning from v1):
#   - catchment_metrics_long: NO year — static morphometric properties, not
#     tied to a data-collection year at all.
#   - canlcc_long: ONE fixed year (canlcc_year, default 2020L, matching
#     CAN_LLC_2020.tif) — the raster is a single-year product; there's no
#     year in the raw column name to parse.
#   - ndvi_timeseries: real per-year values, 1984-2025, parsed from each
#     column name.
#   - harvest_regen_timeseries: real per-year values, 2002-2024, parsed from
#     each column name.
#   - ndvi_trend_summary, harvest_regen_combined_summary: NO year column at
#     all — both are fit/aggregated across their entire record (1984-2025
#     for the trend, 2002-2024 for the "combined" harvest/regen band), not
#     tied to any single year. Earlier draft of this script gave these a
#     year = NA row instead of dropping the column — inconsistent (NA
#     invites "why is this here" for a table that's supposed to be a clean
#     per-year time series) and confirmed as a real problem in review: a
#     value spanning a whole record belongs in its own summary table with
#     no year column, not a fake-year row in the per-year table.
#
# Usage (from an R session, after CAM_streams_hydroweight.csv and
# catchment_metrics.csv already exist — i.e. after run_cam_streams.R has
# been run):
#   source(here("workflow/R/utils.R"))
#   source(here("workflow/CAM/tidy_outputs.R"))
#   tidy_cam_outputs(output_dir = here("output/CAM/stream_delineation"))
# Writes 6 CSVs into output_dir/tidy/: catchment_metrics_long.csv,
# canlcc_long.csv, ndvi_timeseries.csv, ndvi_trend_summary.csv,
# harvest_regen_timeseries.csv, harvest_regen_combined_summary.csv. Each
# tidy_*() function can also be called directly on an already-loaded data
# frame if you don't want the files written.
#
# Dependencies: dplyr, tidyr, readr, fs, glue, tibble (via utils.R);
#   workflow/R/utils.R (cw_inform/cw_warn/cw_abort) must be sourced first.
# =============================================================================

#' Canonical hydroweight weighting-scheme names, in the case the pipeline
#' itself uses for `lumped`/`iEucO`/... — every raw column uses either this
#' exact case (continuous LOIs' distance-weighted columns) or an all-
#' lowercase version of it (categorical LOIs' "_prop" columns; see
#' CLAUDE.md's "known quirks" section on hydroweight_attributes()'s raw
#' categorical naming). normalize_scheme() maps either back to this list.
HYDROWEIGHT_SCHEMES <- c("lumped", "iEucO", "iFLO", "HAiFLO", "iEucS", "iFLS", "HAiFLS")

#' Map a raw scheme token (either case) to its canonical HYDROWEIGHT_SCHEMES
#' spelling. Aborts on any token that doesn't match one of the seven known
#' schemes — a silent NA here would be a mis-parsed column, not a genuine
#' missing value, so this fails loudly rather than producing a quietly
#' wrong tidy row.
#'
#' @param x Character vector of raw scheme tokens (any case).
#' @return Character vector, same length, canonical case.
normalize_scheme <- function(x) {
  idx <- match(tolower(x), tolower(HYDROWEIGHT_SCHEMES))
  if (anyNA(idx)) {
    cw_abort(glue::glue(
      "normalize_scheme(): unrecognized scheme token(s): ",
      "{paste(unique(x[is.na(idx)]), collapse = ', ')}"
    ))
  }
  HYDROWEIGHT_SCHEMES[idx]
}

# -- catchment_metrics.csv ---------------------------------------------------

#' Catchment metric columns that are structurally always-NA for CAM's stream
#' sites (lake morphometrics + the lake-dependent catchment-area:lake-area
#' ratio — carried over from the shared stream+lake catchment_metrics.R
#' module, but stream sites have no lake polygon to compute them from).
#' Dropped from the numeric long output — see tidy_catchment_metrics().
CATCHMENT_METRICS_LAKE_ONLY_COLS <- c(
  "lake_area_ha", "lake_fetch_m", "shoreline_length_m", "shoreline_development",
  "lake_sphericity", "lake_compactness", "lake_area_pct", "ca_la_ratio"
)

#' Tidy catchment_metrics.csv into a numeric long format for faceted plots
#' (e.g. facet_wrap(~metric, scales = "free"))
#'
#' One row per (site, version, metric). Only genuinely numeric metrics are
#' melted into `value` (a real `double` column) — CATCHMENT_METRICS_LAKE_ONLY_COLS
#' (confirmed 100% NA in this project's actual output — not assumed) are
#' dropped rather than melted, with a warning if a future run ever has real
#' data there (so a genuine future lake-adjacent metric isn't silently lost).
#' `aspect_class` (the one categorical metric) is kept as a repeated
#' companion column instead of its own `metric` row, so a faceted plot can
#' still color/filter by it without a join.
#'
#' @param metrics Data frame read from catchment_metrics.csv (site_id,
#'   version, + one column per metric).
#' @return Long tibble: site, version, aspect_class, metric, value (double).
tidy_catchment_metrics <- function(metrics) {
  dropped <- intersect(CATCHMENT_METRICS_LAKE_ONLY_COLS, names(metrics))
  non_na_dropped <- dropped[vapply(dropped, function(col) any(!is.na(metrics[[col]])), logical(1))]
  if (length(non_na_dropped) > 0) {
    cw_warn(glue::glue(
      "tidy_catchment_metrics(): {paste(non_na_dropped, collapse = ', ')} ",
      "expected to be always-NA for stream sites but has real data — ",
      "dropped from catchment_metrics_long.csv anyway; update ",
      "CATCHMENT_METRICS_LAKE_ONLY_COLS in tidy_outputs.R if this LOI ",
      "should now be kept."
    ))
  }

  metrics |>
    dplyr::rename(site = site_id) |>
    dplyr::select(-dplyr::any_of(dropped)) |>
    tidyr::pivot_longer(
      cols = -c(site, version, aspect_class),
      names_to = "metric",
      values_to = "value",
      values_transform = as.double
    ) |>
    dplyr::select(site, version, aspect_class, metric, value) |>
    dplyr::arrange(site, version, metric)
}

# -- CAM_streams_hydroweight.csv ---------------------------------------------

#' Tidy the "canlcc" (land cover) block into a composition-ready long table
#'
#' Raw columns: canlcc_<class>_<scheme>_prop (105 = 15 classes x 7 schemes).
#' `year` is fixed at `year` (the function argument) for every row — see
#' file header. Suited for stacked-bar/composition plots
#' (site on x, prop on y, fill = class, facet/filter by scheme).
#'
#' @param hw   Full hydroweight data frame (site, version, + LOI columns).
#' @param year Integer. Fixed year for every row (default 2020L, matching
#'   CAN_LLC_2020.tif).
#' @return Long tibble: site, version, year, class, scheme, prop.
tidy_hydroweight_canlcc <- function(hw, year = 2020L) {
  cols <- grep("^canlcc_", names(hw), value = TRUE)
  if (length(cols) == 0) {
    return(tibble::tibble(
      site = character(), version = character(), year = integer(),
      class = character(), scheme = character(), prop = double()
    ))
  }

  hw |>
    dplyr::select(site, version, dplyr::all_of(cols)) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(cols),
      names_to = c("class", "scheme"),
      names_pattern = "^canlcc_(.+)_(lumped|ieuco|iflo|haiflo|ieucs|ifls|haifls)_prop$",
      values_to = "prop"
    ) |>
    dplyr::mutate(year = as.integer(year), scheme = normalize_scheme(scheme)) |>
    dplyr::select(site, version, year, class, scheme, prop) |>
    dplyr::arrange(site, version, class, scheme)
}

#' Tidy the "ndvi" (continuous, per-year) block into a year x scheme time
#' series, one row per (site, version, year, scheme)
#'
#' Raw columns: ndvi_<band idx>_NDVI_<year>_<stat>. `mean`/`sd` merge each
#' scheme's own stat: for "lumped", the flat/unweighted mean/sd; for the six
#' distance-weighted schemes, `<scheme>_distwtd_mean/sd` (see file header for
#' why these share one column pair rather than two). `median`/`min`/`max`/
#' `sum`/`cell_count`/`na_cell_count` are only ever computed unweighted, so
#' they're NA for every non-"lumped" row — not a bug, just not computed per
#' weighted scheme. Suited for time-series line/ribbon plots (x = year,
#' y = mean, ymin/ymax from sd, facet/filter by scheme).
#'
#' @param hw Full hydroweight data frame.
#' @return Long tibble: site, version, year, scheme, mean, sd, median, min,
#'   max, sum, cell_count, na_cell_count.
tidy_hydroweight_ndvi_timeseries <- function(hw) {
  cols <- grep("^ndvi_[0-9]+_NDVI_[0-9]{4}_", names(hw), value = TRUE)
  if (length(cols) == 0) {
    return(tibble::tibble(
      site = character(), version = character(), year = integer(), scheme = character(),
      mean = double(), sd = double(), median = double(), min = double(),
      max = double(), sum = double(), cell_count = double(), na_cell_count = double()
    ))
  }

  long <- hw |>
    dplyr::select(site, version, dplyr::all_of(cols)) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(cols),
      names_to = c("year", "scheme_stat"),
      names_pattern = "^ndvi_[0-9]+_NDVI_([0-9]{4})_(.+)$",
      values_to = "raw_value"
    ) |>
    dplyr::mutate(
      year = as.integer(year),
      is_distwtd = grepl("distwtd", scheme_stat),
      scheme = normalize_scheme(dplyr::if_else(
        is_distwtd, sub("_distwtd_(mean|sd)$", "", scheme_stat), "lumped"
      )),
      # flat (lumped) stat name, or "mean"/"sd" for a distwtd column
      flat_stat = dplyr::if_else(
        is_distwtd, sub(".*_distwtd_(mean|sd)$", "\\1", scheme_stat), tolower(scheme_stat)
      )
    )

  long |>
    dplyr::select(site, version, year, scheme, flat_stat, raw_value) |>
    tidyr::pivot_wider(names_from = flat_stat, values_from = raw_value) |>
    dplyr::select(
      site, version, year, scheme,
      mean, sd, median, min, max, sum, cell_count, na_cell_count
    ) |>
    dplyr::arrange(site, version, year, scheme)
}

#' Tidy the "ndvi_trend" (Sen's-slope trend) block into one row per
#' (site, version, scheme)
#'
#' Raw columns: ndvi_trend_<slope|p_value>_<stat>. `<variable>_mean`/
#' `<variable>_sd` merge each scheme's own stat the same way
#' tidy_hydroweight_ndvi_timeseries() does. `<variable>_median`/`_min`/
#' `_max` are only computed unweighted (NA for non-"lumped" rows). No `year`
#' column — the trend spans the entire 1984-2025 record, not any single
#' year (see file header). Suited for scatter/significance-filter plots
#' (e.g. filter(scheme == "lumped", p_value_mean < 0.05)).
#'
#' @param hw Full hydroweight data frame.
#' @return Long tibble: site, version, scheme, slope_mean, slope_sd,
#'   slope_median, slope_min, slope_max, p_value_mean, p_value_sd,
#'   p_value_median, p_value_min, p_value_max.
tidy_hydroweight_ndvi_trend_summary <- function(hw) {
  cols <- grep("^ndvi_trend_", names(hw), value = TRUE)
  empty_cols <- c(
    "slope_mean", "slope_sd", "slope_median", "slope_min", "slope_max",
    "p_value_mean", "p_value_sd", "p_value_median", "p_value_min", "p_value_max"
  )
  if (length(cols) == 0) {
    empty <- tibble::tibble(site = character(), version = character(), scheme = character())
    for (nm in empty_cols) empty[[nm]] <- double()
    return(empty)
  }

  long <- hw |>
    dplyr::select(site, version, dplyr::all_of(cols)) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(cols),
      names_to = "rest",
      names_pattern = "^ndvi_trend_(.+)$",
      values_to = "raw_value"
    ) |>
    dplyr::mutate(
      variable = dplyr::if_else(startsWith(rest, "slope_"), "slope", "p_value"),
      rest2 = sub("^(slope|p_value)_", "", rest),
      is_distwtd = grepl("distwtd", rest2),
      scheme = normalize_scheme(dplyr::if_else(
        is_distwtd, sub("_distwtd_(mean|sd)$", "", rest2), "lumped"
      )),
      flat_stat = dplyr::if_else(
        is_distwtd, sub(".*_distwtd_(mean|sd)$", "\\1", rest2), rest2
      ),
      col_name = paste0(variable, "_", flat_stat)
    )

  long |>
    dplyr::select(site, version, scheme, col_name, raw_value) |>
    tidyr::pivot_wider(names_from = col_name, values_from = raw_value) |>
    dplyr::select(site, version, scheme, dplyr::all_of(empty_cols)) |>
    dplyr::arrange(site, version, scheme)
}

#' Shared first step for both harvest_regen tidy functions below: pivot the
#' raw wide columns into one row per (site, version, year_tok, class,
#' scheme), year_tok still raw ("y2002".."y2024" or the literal "combined").
#' Not exported/used directly — see tidy_hydroweight_harvest_regen_timeseries()
#' and tidy_hydroweight_harvest_regen_combined_summary(), which each filter
#' this to their own year_tok subset.
harvest_regen_pivot_raw <- function(hw) {
  cols <- grep("^harvest_regen_", names(hw), value = TRUE)
  if (length(cols) == 0) {
    return(tibble::tibble(
      site = character(), version = character(), year_tok = character(),
      class = character(), scheme = character(), prop = double()
    ))
  }

  hw |>
    dplyr::select(site, version, dplyr::all_of(cols)) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(cols),
      names_to = c("year_tok", "class", "scheme"),
      names_pattern = "^harvest_regen_(y[0-9]{4}|combined)_(harvest|regen|other)_(lumped|ieuco|iflo|haiflo|ieucs|ifls|haifls)_prop$",
      values_to = "prop"
    ) |>
    dplyr::mutate(scheme = normalize_scheme(scheme))
}

#' Tidy the "harvest_regen" per-year columns into a disturbance-history time
#' series, one row per (site, version, year, class, scheme)
#'
#' ONLY real per-year rows (2002-2024) — the "combined" (all-years) rows
#' live in tidy_hydroweight_harvest_regen_combined_summary() instead, same
#' split as ndvi_timeseries vs. ndvi_trend_summary: a value that isn't tied
#' to one year doesn't belong in a per-year time series table with a fake
#' NA year, it belongs in its own summary table with no year column at all.
#' Suited for stacked-area disturbance-history plots (x = year, y = prop,
#' fill = class).
#'
#' @param hw Full hydroweight data frame.
#' @return Long tibble: site, version, year, class, scheme, prop.
tidy_hydroweight_harvest_regen_timeseries <- function(hw) {
  harvest_regen_pivot_raw(hw) |>
    dplyr::filter(year_tok != "combined") |>
    dplyr::mutate(year = as.integer(sub("^y", "", year_tok))) |>
    dplyr::select(site, version, year, class, scheme, prop) |>
    dplyr::arrange(site, version, class, scheme, year)
}

#' Tidy the "harvest_regen" "combined" (all-years 2002-2024) columns into a
#' summary table, one row per (site, version, class, scheme)
#'
#' No year column — the combined band spans the entire 2002-2024 record,
#' not any single year (see tidy_hydroweight_harvest_regen_timeseries()).
#' Suited for a single-snapshot composition plot, analogous to canlcc_long
#' but for cumulative disturbance instead of current land cover.
#'
#' @param hw Full hydroweight data frame.
#' @return Long tibble: site, version, class, scheme, prop.
tidy_hydroweight_harvest_regen_combined_summary <- function(hw) {
  harvest_regen_pivot_raw(hw) |>
    dplyr::filter(year_tok == "combined") |>
    dplyr::select(site, version, class, scheme, prop) |>
    dplyr::arrange(site, version, class, scheme)
}

# -- Orchestrator -------------------------------------------------------------

#' Read catchment_metrics.csv and CAM_streams_hydroweight.csv from
#' output_dir, tidy both into 5 purpose-shaped long tables, and write the
#' results to output_dir/tidy/
#'
#' Warns (does not abort) if any hydroweight column besides site/version
#' wasn't claimed by one of the 4 known LOI prefixes — a safety net for when
#' loi_layers in run_cam_streams.R gains a new LOI this script doesn't know
#' how to parse yet; those columns are silently excluded from every tidy
#' output until this script is extended, rather than crashing the reshape.
#'
#' @param output_dir  Character. Directory containing catchment_metrics.csv
#'   and CAM_streams_hydroweight.csv (e.g.
#'   here("output/CAM/stream_delineation")).
#' @param canlcc_year Integer. Passed to tidy_hydroweight_canlcc(). Default
#'   2020L.
#' @param metrics_file,hydroweight_file Character. Filenames within
#'   output_dir. Defaults match run_cam_streams.R's actual output names.
#'
#' @return Invisibly, a named list of the 5 tidy tibbles, already written to
#'   disk.
tidy_cam_outputs <- function(
  output_dir,
  canlcc_year = 2020L,
  metrics_file = "catchment_metrics.csv",
  hydroweight_file = "CAM_streams_hydroweight.csv"
) {
  metrics_path <- fs::path(output_dir, metrics_file)
  hydroweight_path <- fs::path(output_dir, hydroweight_file)
  if (!fs::file_exists(metrics_path)) {
    cw_abort(glue::glue("tidy_cam_outputs(): not found: {metrics_path}"))
  }
  if (!fs::file_exists(hydroweight_path)) {
    cw_abort(glue::glue("tidy_cam_outputs(): not found: {hydroweight_path}"))
  }

  metrics <- readr::read_csv(metrics_path, show_col_types = FALSE)
  hw <- readr::read_csv(hydroweight_path, show_col_types = FALSE)

  claimed_cols <- c(
    grep("^canlcc_", names(hw), value = TRUE),
    grep("^ndvi_[0-9]+_NDVI_[0-9]{4}_", names(hw), value = TRUE),
    grep("^ndvi_trend_", names(hw), value = TRUE),
    grep("^harvest_regen_", names(hw), value = TRUE)
  )
  unclaimed <- setdiff(names(hw), c("site", "version", claimed_cols))
  if (length(unclaimed) > 0) {
    cw_warn(glue::glue(
      "tidy_cam_outputs(): {length(unclaimed)} column(s) not recognized by ",
      "any known LOI parser (excluded from every tidy output) — a new LOI ",
      "was likely added to run_cam_streams.R's loi_layers without a ",
      "matching tidy_hydroweight_<loi>() here. First few: ",
      "{paste(utils::head(unclaimed, 5), collapse = ', ')}"
    ))
  }

  result <- list(
    catchment_metrics_long = tidy_catchment_metrics(metrics),
    canlcc_long = tidy_hydroweight_canlcc(hw, year = canlcc_year),
    ndvi_timeseries = tidy_hydroweight_ndvi_timeseries(hw),
    ndvi_trend_summary = tidy_hydroweight_ndvi_trend_summary(hw),
    harvest_regen_timeseries = tidy_hydroweight_harvest_regen_timeseries(hw),
    harvest_regen_combined_summary = tidy_hydroweight_harvest_regen_combined_summary(hw)
  )

  tidy_dir <- fs::path(output_dir, "tidy")
  fs::dir_create(tidy_dir)

  # Remove files superseded by earlier design iterations, if present, so a
  # stale and current tidy output never silently coexist.
  stale_files <- c("CAM_streams_hydroweight_long.csv", "harvest_regen_long.csv")
  for (f in stale_files) {
    stale_path <- fs::path(tidy_dir, f)
    if (fs::file_exists(stale_path)) {
      fs::file_delete(stale_path)
    }
  }

  purrr::iwalk(result, function(df, nm) {
    readr::write_csv(df, fs::path(tidy_dir, paste0(nm, ".csv")))
  })

  cw_inform(glue::glue(
    "tidy_cam_outputs(): wrote {length(result)} file(s) to {tidy_dir}/ ",
    "({paste(names(result), vapply(result, nrow, integer(1)), sep = ': ', collapse = ' rows, ')} rows)."
  ))

  invisible(result)
}
