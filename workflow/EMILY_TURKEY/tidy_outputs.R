# tidy_outputs.R
# =============================================================================
# Reshapes EMILY_TURKEY_hydroweight.csv (wide, one row per site x version)
# into a purpose-shaped long table for plotting/analysis. Same design as
# workflow/CAM/tidy_outputs.R and workflow/CELESTE/tidy_outputs.R (see
# CAM's file header for the full rationale) — real stat columns (`prop`
# here), not a generic melted stat/value pair.
#
# EMILY_TURKEY currently has exactly ONE hydroweight LOI — CANLCC land
# cover (see workflow/EMILY_TURKEY/README.md's "Hydroweight LOIs" section:
# NDVI and harvest/regen prep don't exist yet for this project's AOI), so
# this script covers a single output table: canlcc_long. tidy_hydroweight_
# canlcc() itself is unchanged from workflow/CAM/tidy_outputs.R / workflow/
# CELESTE/tidy_outputs.R's version — run_emily_turkey.R's loi_layers reuses
# the identical CAN_LLC_2020.tif source and canlcc_levels class scheme as
# both of those projects (see run_emily_turkey.R's Stage 7), so the raw
# column format is identical: canlcc_<class>_<scheme>_prop.
#
# YEAR: fixed at `canlcc_year` (default 2020L, matching CAN_LLC_2020.tif) —
# the raster is a single-year product; there's no year in the raw column
# name to parse.
#
# NOT covered (unlike CAM's tidy_outputs.R): catchment_metrics_long. CAM's
# version tidies catchment_metrics.csv too; CELESTE's doesn't. This file
# follows CELESTE's precedent — hydroweight only — since nothing about
# EMILY_TURKEY's own catchment_metrics.csv differs from either project's
# in a way that would need its own tidy pass; add tidy_catchment_metrics()
# here (copied from workflow/CAM/tidy_outputs.R) if/when that's wanted.
#
# When NDVI/harvest-regen LOI prep is added for this project's AOI, extend
# this file the same way workflow/CELESTE/tidy_outputs.R extends CAM's
# pattern — one tidy_hydroweight_<loi>() per LOI, added to the orchestrator
# below and to its claimed_cols check.
#
# Usage (from an R session, after EMILY_TURKEY_hydroweight.csv already
# exists — i.e. after run_emily_turkey.R has been run):
#   source(here("workflow/R/utils.R"))
#   source(here("workflow/EMILY_TURKEY/tidy_outputs.R"))
#   tidy_emily_turkey_outputs(output_dir = here("output/EMILY_TURKEY"))
# Writes 1 CSV into output_dir/tidy/: canlcc_long.csv.
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
#' silent NA here would be a mis-parsed column, not a genuine missing value.
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

#' Tidy the "canlcc" (land cover) block into a composition-ready long table
#'
#' Raw columns: canlcc_<class>_<scheme>_prop (15 classes x 7 schemes, per
#' run_emily_turkey.R's canlcc_levels — fewer classes than CAM/CELESTE's
#' output columns show only because not every class is actually present
#' across this project's AOI, same "omitted means truly absent, not this
#' function's doing" behavior CLAUDE.md's hydroweight quirks section
#' describes). `year` is fixed at `year` (the function argument) for every
#' row — see file header. Suited for stacked-bar/composition plots (site on
#' x, prop on y, fill = class, facet/filter by scheme).
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

# -- Orchestrator -------------------------------------------------------------

#' Read EMILY_TURKEY_hydroweight.csv from output_dir, tidy its (sole)
#' canlcc LOI into a composition-ready long table, and write the result to
#' output_dir/tidy/
#'
#' Warns (does not abort) if any hydroweight column besides site/version
#' wasn't claimed by the canlcc parser — a safety net for when loi_layers
#' in run_emily_turkey.R gains a new LOI this script doesn't know how to
#' parse yet; those columns are silently excluded from the tidy output
#' until this script is extended, rather than crashing the reshape.
#'
#' @param output_dir  Character. Directory containing
#'   EMILY_TURKEY_hydroweight.csv (e.g. here("output/EMILY_TURKEY")).
#' @param canlcc_year Integer. Passed to tidy_hydroweight_canlcc(). Default
#'   2020L.
#' @param hydroweight_file Character. Filename within output_dir. Default
#'   matches run_emily_turkey.R's actual output name.
#'
#' @return Invisibly, a named list holding the one tidy tibble (canlcc_long),
#'   already written to disk.
tidy_emily_turkey_outputs <- function(
  output_dir,
  canlcc_year = 2020L,
  hydroweight_file = "EMILY_TURKEY_hydroweight.csv"
) {
  hydroweight_path <- fs::path(output_dir, hydroweight_file)
  if (!fs::file_exists(hydroweight_path)) {
    cw_abort(glue::glue("tidy_emily_turkey_outputs(): not found: {hydroweight_path}"))
  }

  hw <- readr::read_csv(hydroweight_path, show_col_types = FALSE)

  claimed_cols <- grep("^canlcc_", names(hw), value = TRUE)
  unclaimed <- setdiff(names(hw), c("site", "version", claimed_cols))
  if (length(unclaimed) > 0) {
    cw_warn(glue::glue(
      "tidy_emily_turkey_outputs(): {length(unclaimed)} column(s) not ",
      "recognized by the canlcc parser (excluded from the tidy output) — ",
      "a new LOI was likely added to run_emily_turkey.R's loi_layers ",
      "without a matching tidy_hydroweight_<loi>() here. First few: ",
      "{paste(utils::head(unclaimed, 5), collapse = ', ')}"
    ))
  }

  result <- list(
    canlcc_long = tidy_hydroweight_canlcc(hw, year = canlcc_year)
  )

  tidy_dir <- fs::path(output_dir, "tidy")
  fs::dir_create(tidy_dir)

  purrr::iwalk(result, function(df, nm) {
    readr::write_csv(df, fs::path(tidy_dir, paste0(nm, ".csv")))
  })

  cw_inform(glue::glue(
    "tidy_emily_turkey_outputs(): wrote {length(result)} file(s) to ",
    "{tidy_dir}/ ({paste(names(result), vapply(result, nrow, integer(1)), sep = ': ', collapse = ' rows, ')} rows)."
  ))

  invisible(result)
}
