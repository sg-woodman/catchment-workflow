# fix_lake_bisection.R
# =============================================================================
# Corrects CAM stream-site catchments that bisect a lake — the boundary cuts
# through a lake polygon instead of fully containing it or fully excluding
# it, which is hydrologically impossible.
#
# Root cause: point-based delineation has zero lake awareness —
# wbt_watershed() traces strictly by per-cell D8 flow direction, and
# nothing guarantees a lake polygon it happens to cross is flow-consistent.
# See workflow/R/lake_containment.R's header for the general writeup, and
# workflow/CAM/README.md's "Lake-bisection fix" section for the scan
# results, the affected sites, and why CAM streams has no upfront-
# prevention equivalent (unlike CELESTE — see workflow/CELESTE/README.md).
#
# WORKFLOW: validation is automatic and side-effect-free; correction is a
# deliberate, manually-gated second step — sourcing this file never rewrites
# a catchment on its own.
#
#   1. source(here::here("workflow/CAM/run_cam_streams.R"))   # sites, group_manifest, config, output_dir, cache_dir, loi_layers
#   2. source(here::here("workflow/CAM/fix_lake_bisection.R")) # runs validate_catchment_lake_intersections(), prints/writes the report, defines correct_lake_bisected_sites() -- corrects NOTHING yet
#   3. Review output/CAM/stream_delineation/lake_intersection_report.csv
#      (site_id, lake_name, lake_area_ha, pct_lake_outside) AND the
#      catchments in QGIS against ohn_waterbodies_valid.gpkg -- not every
#      flagged pair is worth correcting (e.g. a lake sitting right at a
#      site's edge may be a real, if unusual, drainage pattern rather than
#      a delineation defect).
#   4. Once you've decided which site_ids actually need fixing:
#        correct_lake_bisected_sites(
#          site_ids = c("Tilton", "SUD17"),
#          sites = sites, group_manifest = group_manifest, config = config,
#          output_dir = output_dir, cache_dir = cache_dir, loi_layers = loi_layers
#        )
#      Re-validates fresh at call time (not against the possibly-stale
#      `flagged` object from step 2), so this is safe to call in a later
#      session too, as long as run_cam_streams.R has been re-sourced first.
#
# Safe to call correct_lake_bisected_sites() repeatedly, including as a
# second/third pass after reviewing residual flags — accumulates across
# passes, never delete cache_dir/lake_corrected/ between them. Some
# densely-lake-packed clusters can show a whack-a-mole pattern across
# passes; see workflow/CAM/README.md's "Lake-bisection fix" section for
# the confirmed cases and why a small residual may be worth leaving as
# manual-review-only rather than chasing with another pass.
# =============================================================================

source(here::here("workflow/R/lake_containment.R"))

OIH_LAKES_PATH <- "/Users/sam/Documents/cfs/shared_data/raw/hydro/waterbodies/ON_waterbodies/ohn_waterbodies_valid.gpkg"

# -- Validate: always runs, corrects nothing ---------------------------------

flagged <- validate_catchment_lake_intersections(
  sites = sites, output_dir = output_dir, lakes_path = OIH_LAKES_PATH,
  catchment_file = "catchment.gpkg"
)

if (nrow(flagged) == 0) {
  cw_inform("No bisected lakes found.")
} else {
  cw_inform(glue::glue(
    "{nrow(flagged)} bisected (site, lake) pair(s) across ",
    "{dplyr::n_distinct(flagged$site_id)} site(s) found: ",
    "{paste(unique(flagged$site_id), collapse = ', ')}. Review ",
    "{fs::path(output_dir, 'lake_intersection_report.csv')} and the ",
    "catchments in QGIS, then call correct_lake_bisected_sites(site_ids = c(...)) ",
    "for the ones confirmed to need fixing."
  ))
}

#' Correct catchments for a manually-confirmed list of site IDs
#'
#' Re-validates fresh (rather than trusting a possibly-stale `flagged`
#' object from when this file was sourced), filters to the requested
#' site_ids, flattens their bisected lakes (accumulated with any lake
#' flattened by a previous call — see this file's header), redelineates
#' just those sites via rerun_engine_sites() (which cascades remove-
#' upstream/reclip/metrics/hydroweight for them and any cascaded
#' neighbor), and re-validates to report a before/after summary.
#'
#' @param site_ids     Character vector. The site_id(s) to correct — a
#'   subset of what validate_catchment_lake_intersections() currently
#'   flags, chosen after manual review (report + QGIS).
#' @param sites          Validated sites tibble (as built by
#'   build_engine_group_manifest())
#' @param group_manifest sf tibble from build_engine_group_manifest()
#' @param config         Resolved config from resolve_engine_config()
#' @param output_dir     Character. Root output directory
#' @param cache_dir      Character. Project cache directory
#' @param loi_layers     LOI layers list for hydroweight, or NULL to skip
#'   re-running it
#' @param lakes_path     Character. Path to the OHN/OIH waterbody layer.
#'   Default OIH_LAKES_PATH.
#'
#' @return Invisibly, the before/after summary tibble
correct_lake_bisected_sites <- function(
  site_ids, sites, group_manifest, config, output_dir, cache_dir,
  loi_layers = NULL, lakes_path = OIH_LAKES_PATH
) {
  unknown_ids <- setdiff(site_ids, sites$site_id)
  if (length(unknown_ids) > 0) {
    cw_abort(glue::glue("Unknown site_id(s): {paste(unknown_ids, collapse = ', ')}"))
  }

  current_flagged <- validate_catchment_lake_intersections(
    sites = sites, output_dir = output_dir, lakes_path = lakes_path,
    catchment_file = "catchment.gpkg"
  )

  to_fix <- dplyr::filter(current_flagged, site_id %in% site_ids)
  not_flagged <- setdiff(site_ids, to_fix$site_id)
  if (length(not_flagged) > 0) {
    cw_warn(glue::glue(
      "Site(s) requested but not currently flagged (already resolved, or ",
      "below threshold) — skipping: {paste(not_flagged, collapse = ', ')}"
    ))
  }
  if (nrow(to_fix) == 0) {
    cw_inform("Nothing to correct for the requested site_id(s).")
    return(invisible(NULL))
  }

  target_ids <- unique(to_fix$site_id)
  cw_inform(glue::glue(
    "Correcting {length(target_ids)} site(s): {paste(target_ids, collapse = ', ')} ",
    "({nrow(to_fix)} bisected (site, lake) pair(s))."
  ))

  before_areas <- purrr::map_dfr(target_ids, function(sid) {
    p <- fs::path(site_output_dir(output_dir, sid), "catchment.gpkg")
    tibble::tibble(
      site_id = sid,
      area_km2_before = as.numeric(sf::st_area(sf::st_read(p, quiet = TRUE))) / 1e6
    )
  })

  lakes_full <- sf::st_read(lakes_path, quiet = TRUE)
  flagged_lake_polys <- lakes_full[lakes_full$OGF_ID %in% unique(to_fix$OGF_ID), ]
  cw_inform(glue::glue("Flattening {nrow(flagged_lake_polys)} distinct lake(s)."))

  corrected_pointer <- prepare_lake_corrected_flow_pointer(
    cache_dir = cache_dir, lake_polys = flagged_lake_polys
  )

  rerun_engine_sites(
    edited_snap_site_ids = target_ids,
    sites = sites, group_manifest = group_manifest, config = config,
    output_dir = output_dir, cache_dir = cache_dir, loi_layers = loi_layers,
    flow_pointer_override_path = corrected_pointer
  )

  after_unclipped <- validate_catchment_lake_intersections(
    sites = sites, output_dir = output_dir, lakes_path = lakes_path,
    catchment_file = "catchment.gpkg"
  )
  after_clipped <- validate_catchment_lake_intersections(
    sites = sites, output_dir = output_dir, lakes_path = lakes_path,
    catchment_file = "catchment_clipped.gpkg"
  )

  still_flagged <- intersect(
    union(after_unclipped$site_id, after_clipped$site_id), target_ids
  )
  fixed <- setdiff(target_ids, still_flagged)

  after_areas <- purrr::map_dfr(target_ids, function(sid) {
    p <- fs::path(site_output_dir(output_dir, sid), "catchment.gpkg")
    tibble::tibble(
      site_id = sid,
      area_km2_after = as.numeric(sf::st_area(sf::st_read(p, quiet = TRUE))) / 1e6
    )
  })

  summary_tbl <- dplyr::left_join(before_areas, after_areas, by = "site_id") |>
    dplyr::mutate(status = dplyr::if_else(site_id %in% fixed, "fixed", "still flagged"))

  cw_inform("\n-- before/after summary --")
  print(summary_tbl)

  if (length(still_flagged) > 0) {
    cw_warn(glue::glue(
      "{length(still_flagged)} site(s) still bisect a lake after correction ",
      "— review the residual in {fs::path(output_dir, 'lake_intersection_report.csv')} ",
      "before deciding whether another pass is worth it: ",
      "{paste(still_flagged, collapse = ', ')}"
    ))
  } else {
    cw_inform("All requested site(s) corrected.")
  }

  invisible(summary_tbl)
}
