# fix_lake_bisection.R
# =============================================================================
# Corrects CAM stream-site catchments that bisect a lake — the boundary cuts
# through a lake polygon instead of fully containing it or fully excluding
# it, which is hydrologically impossible.
#
# Root cause: point-based delineation (workflow/R/engine/04_delineate_site.R's
# delineate_engine_point_site()) has zero lake awareness — wbt_watershed()
# traces strictly by per-cell D8 flow direction from a single pour point,
# and nothing guarantees a lake polygon it happens to cross is
# flow-consistent (all interior cells routing to one outlet). OIH's own
# hydro-enforcement only "improves flow direction" for waterbodies
# intersecting the mapped Enhanced Watercourse network — not universally.
# See workflow/R/lake_containment.R's header and the plan this was built
# from for the full writeup.
#
# Confirmed on the delivered CAM streams output (2026-08-27): 14 (site,
# lake) pairs across 7 sites — Tilton, Daisy, NCMN, SUD17, SUD102 (x3),
# SUD103 (x3), SUD200 (x3) — genuinely bisect a lake >= 1 ha, excluding
# River/Pond waterbody types and boundary-touch noise (<2%/>98% overlap).
#
# SCOPE: CAM streams only. CELESTE uses a different terrain tier (raw MRDEM
# + per-group HydroBasins breach, not OIH's flow_direction tier) and a
# different lake source (NHN, not scanned) — a separate follow-on, not
# covered here.
#
# USAGE: source workflow/CAM/run_cam_streams.R once first (cache-checked,
# so this is cheap/idempotent even on a completed project) to get `sites`,
# `group_manifest`, `config`, `output_dir`, `cache_dir`, `loi_layers` in
# memory, then source this file:
#
#   source(here::here("workflow/CAM/run_cam_streams.R"))
#   source(here::here("workflow/CAM/fix_lake_bisection.R"))
#
# Safe to re-run: validate_catchment_lake_intersections() and
# rerun_engine_sites() are both idempotent — a clean second run finds
# nothing left to flag and does nothing further, except that
# prepare_lake_corrected_flow_pointer()'s cache under
# cache_dir/lake_corrected/ is NOT auto-invalidated if the flagged-lake set
# changes on a later run (e.g. new sites added) — delete that directory
# manually first if so.
# =============================================================================

source(here::here("workflow/R/lake_containment.R"))

OIH_LAKES_PATH <- "/Users/sam/Documents/cfs/shared_data/raw/hydro/waterbodies/ON_waterbodies/ohn_waterbodies_valid.gpkg"

# -- Step 1: find bisected (site, lake) pairs --------------------------------

flagged <- validate_catchment_lake_intersections(
  sites = sites, output_dir = output_dir, lakes_path = OIH_LAKES_PATH,
  catchment_file = "catchment.gpkg"
)

if (nrow(flagged) == 0) {
  cw_inform("No bisected lakes found — nothing to correct.")
} else {
  flagged_site_ids <- unique(flagged$site_id)
  cw_inform(glue::glue(
    "{nrow(flagged)} bisected (site, lake) pair(s) across ",
    "{length(flagged_site_ids)} site(s): {paste(flagged_site_ids, collapse = ', ')}"
  ))

  # Snapshot "before" catchment areas for the before/after summary below.
  before_areas <- purrr::map_dfr(flagged_site_ids, function(sid) {
    p <- fs::path(site_output_dir(output_dir, sid), "catchment.gpkg")
    tibble::tibble(
      site_id = sid,
      area_km2_before = as.numeric(sf::st_area(sf::st_read(p, quiet = TRUE))) / 1e6
    )
  })

  # -- Step 2: look up the flagged lakes' geometries --------------------------

  lakes_full <- sf::st_read(OIH_LAKES_PATH, quiet = TRUE)
  flagged_lake_polys <- lakes_full[lakes_full$OGF_ID %in% unique(flagged$OGF_ID), ]
  cw_inform(glue::glue("Flattening {nrow(flagged_lake_polys)} distinct lake(s)."))

  # -- Step 3: build the lake-corrected flow pointer ---------------------------

  corrected_pointer <- prepare_lake_corrected_flow_pointer(
    cache_dir = cache_dir, lake_polys = flagged_lake_polys
  )

  # -- Step 4: redelineate the affected sites, cascade downstream stages ------

  rerun_engine_sites(
    edited_snap_site_ids = flagged_site_ids,
    sites = sites, group_manifest = group_manifest, config = config,
    output_dir = output_dir, cache_dir = cache_dir, loi_layers = loi_layers,
    flow_pointer_override_path = corrected_pointer
  )

  # -- Step 5: re-validate, report before/after --------------------------------

  after_unclipped <- validate_catchment_lake_intersections(
    sites = sites, output_dir = output_dir, lakes_path = OIH_LAKES_PATH,
    catchment_file = "catchment.gpkg"
  )
  after_clipped <- validate_catchment_lake_intersections(
    sites = sites, output_dir = output_dir, lakes_path = OIH_LAKES_PATH,
    catchment_file = "catchment_clipped.gpkg"
  )

  still_flagged <- union(after_unclipped$site_id, after_clipped$site_id)
  fixed <- setdiff(flagged_site_ids, still_flagged)

  after_areas <- purrr::map_dfr(flagged_site_ids, function(sid) {
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
      "— needs manual review (edit pour_point_snapped.shp in QGIS, then ",
      "rerun_engine_site_watershed()): {paste(still_flagged, collapse = ', ')}"
    ))
  } else {
    cw_inform("All flagged sites corrected.")
  }
}
