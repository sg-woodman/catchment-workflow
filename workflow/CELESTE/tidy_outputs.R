# tidy_outputs.R
# =============================================================================
# Reshapes CELESTE_engine_hydroweight.csv (wide, one row per site x version)
# into purpose-shaped long tables for plotting/analysis — one file per LOI,
# real stat columns (mean/sd/prop), not a generic melted stat/value pair.
# Same design as workflow/CAM/tidy_outputs.R (see that file's header for the
# full rationale — unweighted + distance-weighted stats merge into one
# column pair per scheme, scheme stays a filterable row, record-spanning
# values get their own no-year summary table rather than a fake NA year).
# CELESTE has no "canlcc" LOI (that's CAM-only), so this script covers
# exactly 4 output tables: ndvi_timeseries, ndvi_trend_summary,
# harvest_regen_timeseries, harvest_regen_combined_summary.
#
# ONE REAL DIFFERENCE FROM CAM'S ndvi COLUMNS — confirmed by inspection, not
# assumed: CELESTE's raw NDVI columns are "ndvi_ndvi_mosaic_<N>_<stat>"
# (N = 1-based band index into the per-group mosaic VRT), NOT CAM's
# "ndvi_<idx>_NDVI_<year>_<stat>" — the calendar year isn't in the column
# name at all. workflow/CELESTE/prepare_ndvi.R's prepare_ndvi_per_group_
# rasters() renames each mosaicked band to "ndvi_mosaic_<N>" (1-based,
# discarding the year the RAW source tile's own band name carried — e.g.
# "0_NDVI_1984" -> mosaic band 1). Confirmed directly against a raw source
# tile (data/ndvi/Landsat_NDVI_1984_2025_Celeste_Coc.tif): band 1 =
# "0_NDVI_1984", band 42 = "41_NDVI_2025" — sequential, no gaps, no
# reordering. So year = 1983 + N for CELESTE's ndvi_mosaic_N.
#
# Usage (from an R session, after CELESTE_engine_hydroweight.csv already
# exists — i.e. after run_celeste.R has been run):
#   source(here("workflow/R/utils.R"))
#   source(here("workflow/CELESTE/tidy_outputs.R"))
#   tidy_celeste_outputs(output_dir = here("output/CELESTE_engine"))
# Writes 4 CSVs into output_dir/tidy/: ndvi_timeseries.csv,
# ndvi_trend_summary.csv, harvest_regen_timeseries.csv,
# harvest_regen_combined_summary.csv.
#
# Dependencies: dplyr, tidyr, readr, fs, glue, tibble, purrr (via utils.R);
#   workflow/R/utils.R (cw_inform/cw_warn/cw_abort) must be sourced first.
# =============================================================================

#' Canonical hydroweight weighting-scheme names — see workflow/CAM/
#' tidy_outputs.R's HYDROWEIGHT_SCHEMES for the full rationale (raw columns
#' use either this exact case or an all-lowercase version).
HYDROWEIGHT_SCHEMES <- c("lumped", "iEucO", "iFLO", "HAiFLO", "iEucS", "iFLS", "HAiFLS")

#' Map a raw scheme token (either case) to its canonical spelling. Aborts
#' on an unrecognized token — see workflow/CAM/tidy_outputs.R for why a
#' silent NA here would be worse than failing loudly.
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

#' Tidy the "ndvi" (continuous, per-year) block into a year x scheme time
#' series, one row per (site, version, year, scheme)
#'
#' Raw columns: ndvi_ndvi_mosaic_<band>_<stat> (band 1-based, year = 1983 +
#' band — see file header). `mean`/`sd` merge each scheme's own stat, same
#' as workflow/CAM/tidy_outputs.R's ndvi function. `median`/`min`/`max`/
#' `sum`/`cell_count`/`na_cell_count` are only ever computed unweighted (NA
#' for non-"lumped" rows) — NOTE: CELESTE's own loi_layers config restricts
#' `stats` to exclude sum/cell_count/NA_cell_count entirely (see
#' run_celeste.R), so those 3 won't be present in the source data at all
#' here — included in the output schema
#' as NA-only columns for consistency with CAM's shape, not because
#' CELESTE actually computed them.
#'
#' @param hw Full hydroweight data frame.
#' @return Long tibble: site, version, year, scheme, mean, sd, median, min,
#'   max, sum, cell_count, na_cell_count.
tidy_hydroweight_ndvi_timeseries <- function(hw) {
  cols <- grep("^ndvi_ndvi_mosaic_[0-9]+_", names(hw), value = TRUE)
  empty_cols <- c("mean", "sd", "median", "min", "max", "sum", "cell_count", "na_cell_count")
  if (length(cols) == 0) {
    empty <- tibble::tibble(site = character(), version = character(), year = integer(), scheme = character())
    for (nm in empty_cols) empty[[nm]] <- double()
    return(empty)
  }

  long <- hw |>
    dplyr::select(site, version, dplyr::all_of(cols)) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(cols),
      names_to = c("band", "scheme_stat"),
      names_pattern = "^ndvi_ndvi_mosaic_([0-9]+)_(.+)$",
      values_to = "raw_value"
    ) |>
    dplyr::mutate(
      year = 1983L + as.integer(band),
      is_distwtd = grepl("distwtd", scheme_stat),
      scheme = normalize_scheme(dplyr::if_else(
        is_distwtd, sub("_distwtd_(mean|sd)$", "", scheme_stat), "lumped"
      )),
      flat_stat = dplyr::if_else(
        is_distwtd, sub(".*_distwtd_(mean|sd)$", "\\1", scheme_stat), tolower(scheme_stat)
      )
    )

  long |>
    dplyr::select(site, version, year, scheme, flat_stat, raw_value) |>
    tidyr::pivot_wider(names_from = flat_stat, values_from = raw_value) |>
    dplyr::select(site, version, year, scheme, dplyr::any_of(empty_cols)) |>
    dplyr::arrange(site, version, year, scheme)
}

#' Tidy the "ndvi_trend" (Sen's-slope trend) block into one row per
#' (site, version, scheme) — column shape identical to CAM's, see
#' workflow/CAM/tidy_outputs.R's version for the full docstring.
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

#' Shared first step for both harvest_regen tidy functions — see
#' workflow/CAM/tidy_outputs.R's version for the full docstring.
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

#' Tidy the "harvest_regen" per-year columns into a disturbance-history
#' time series — real per-year rows only (2002-2024). See workflow/CAM/
#' tidy_outputs.R's version for the full docstring.
tidy_hydroweight_harvest_regen_timeseries <- function(hw) {
  harvest_regen_pivot_raw(hw) |>
    dplyr::filter(year_tok != "combined") |>
    dplyr::mutate(year = as.integer(sub("^y", "", year_tok))) |>
    dplyr::select(site, version, year, class, scheme, prop) |>
    dplyr::arrange(site, version, class, scheme, year)
}

#' Tidy the "harvest_regen" "combined" (all-years) columns into a summary
#' table — no year column, spans the whole 2002-2024 record. See
#' workflow/CAM/tidy_outputs.R's version for the full docstring.
tidy_hydroweight_harvest_regen_combined_summary <- function(hw) {
  harvest_regen_pivot_raw(hw) |>
    dplyr::filter(year_tok == "combined") |>
    dplyr::select(site, version, class, scheme, prop) |>
    dplyr::arrange(site, version, class, scheme)
}

# -- Orchestrator -------------------------------------------------------------

#' Read CELESTE_engine_hydroweight.csv from output_dir, tidy it into 4
#' purpose-shaped long tables, and write the results to output_dir/tidy/
#'
#' Warns (does not abort) if any column besides site/version wasn't claimed
#' by one of the 3 known LOI prefixes (ndvi, ndvi_trend, harvest_regen) — a
#' safety net for a future LOI added to loi_layers without a matching
#' tidy_hydroweight_<loi>() here.
#'
#' @param output_dir      Character. Directory containing
#'   CELESTE_engine_hydroweight.csv (e.g. here("output/CELESTE_engine")).
#' @param hydroweight_file Character. Filename within output_dir. Default
#'   matches run_celeste.R's actual output name.
#'
#' @return Invisibly, a named list of the 4 tidy tibbles, already written
#'   to disk.
tidy_celeste_outputs <- function(output_dir, hydroweight_file = "CELESTE_engine_hydroweight.csv") {
  hydroweight_path <- fs::path(output_dir, hydroweight_file)
  if (!fs::file_exists(hydroweight_path)) {
    cw_abort(glue::glue("tidy_celeste_outputs(): not found: {hydroweight_path}"))
  }

  hw <- readr::read_csv(hydroweight_path, show_col_types = FALSE)

  claimed_cols <- c(
    grep("^ndvi_ndvi_mosaic_[0-9]+_", names(hw), value = TRUE),
    grep("^ndvi_trend_", names(hw), value = TRUE),
    grep("^harvest_regen_", names(hw), value = TRUE)
  )
  unclaimed <- setdiff(names(hw), c("site", "version", claimed_cols))
  if (length(unclaimed) > 0) {
    cw_warn(glue::glue(
      "tidy_celeste_outputs(): {length(unclaimed)} column(s) not recognized ",
      "by any known LOI parser (excluded from every tidy output) — a new ",
      "LOI was likely added to loi_layers without a matching ",
      "tidy_hydroweight_<loi>() here. First few: ",
      "{paste(utils::head(unclaimed, 5), collapse = ', ')}"
    ))
  }

  result <- list(
    ndvi_timeseries = tidy_hydroweight_ndvi_timeseries(hw),
    ndvi_trend_summary = tidy_hydroweight_ndvi_trend_summary(hw),
    harvest_regen_timeseries = tidy_hydroweight_harvest_regen_timeseries(hw),
    harvest_regen_combined_summary = tidy_hydroweight_harvest_regen_combined_summary(hw)
  )

  tidy_dir <- fs::path(output_dir, "tidy")
  fs::dir_create(tidy_dir)

  purrr::iwalk(result, function(df, nm) {
    readr::write_csv(df, fs::path(tidy_dir, paste0(nm, ".csv")))
  })

  cw_inform(glue::glue(
    "tidy_celeste_outputs(): wrote {length(result)} file(s) to {tidy_dir}/ ",
    "({paste(names(result), vapply(result, nrow, integer(1)), sep = ': ', collapse = ' rows, ')} rows)."
  ))

  invisible(result)
}
