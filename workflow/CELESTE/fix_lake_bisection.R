# fix_lake_bisection.R (CELESTE)
# =============================================================================
# Corrects CELESTE catchments that bisect a lake — the boundary cuts through
# a lake polygon instead of fully containing it or fully excluding it, which
# is hydrologically impossible. Same root cause and same fix as the CAM
# streams version (workflow/CAM/fix_lake_bisection.R, workflow/R/
# lake_containment.R's header) — point-based delineation has zero lake
# awareness, so nothing guarantees a lake a site's flow path happens to
# cross is flow-consistent under raw D8.
#
# Confirmed via a full systematic scan (this session): 31 (site, lake)
# pairs across 27 of 132 sites, in 5 of 6 HydroBasins groups (COC, KEN,
# NBE, NIP, TUR — MOR's 2 sites are clean). Matches and exceeds a set of 12
# sites spotted manually in QGIS beforehand (11 confirmed genuine;
# BLOXHAM1_1K_1 checks out clean — likely conflated with the genuinely
# bisected BLOXHAM1_10K_1).
#
# TWO STRUCTURAL DIFFERENCES FROM CAM STREAMS, both handled here rather
# than by changing the shared workflow/R/lake_containment.R machinery:
#
#   1. Lake source is NHN (NHN_HD_WATERBODY_2), read per-HydroBasins-group
#      from local GDBs, not one province-wide file like OHN. See
#      fetch_nhn_lakes_for_group() below.
#
#      IMPORTANT: this deliberately does NOT reuse stream/burn_streams.R's
#      read_merge_nhn_layer() — that function clips each waterbody's
#      geometry to its query AOI, correct for its original job (trimming
#      burn-in streams to a group's conditioning extent) but wrong here:
#      it distorts the containment percentage this check depends on.
#      Confirmed directly during scoping: BLOXHAM1_1K_1 read as 0.17%
#      outside with a tight per-site AOI and 7.3% with a larger group-wide
#      AOI, off the SAME lake, purely because of how much of it each AOI
#      happened to clip off. This file uses read_nhn_from_gdb() directly
#      instead (already reads+reprojects with NO clipping) and only uses
#      the AOI for a non-destructive sf::st_filter() row subset, exactly
#      matching how the CAM version filters OHN.
#
#   2. Terrain tier is raw MRDEM + per-group breach, with NHN streams
#      burned in BEFORE breaching for 5 of 6 groups (dem_burned.tif exists
#      for KEN/MOR/NBE/NIP/TUR — exactly the 5 groups with confirmed
#      bisections; COC has burn_streams = FALSE per run_celeste.R's own
#      override and only has dem.tif). prepare_lake_corrected_flow_pointer()
#      (workflow/R/lake_containment.R) already auto-detects dem_burned.tif
#      over dem.tif for exactly this reason — no CELESTE-specific change
#      needed there.
#
#   A third, purely mechanical difference: a batch of site_ids to correct
#   can span MULTIPLE HydroBasins groups at once, each with its own D8
#   pointer grid, whereas rerun_engine_sites(flow_pointer_override_path = ...)
#   takes exactly one pointer per call. Handled by looping per group in
#   correct_lake_bisected_sites_celeste() below — workflow/R/engine/
#   99_rerun_sites is NOT modified further; it's already correct per
#   individual call.
#
# WORKFLOW: same manually-gated pattern as CAM streams — sourcing this file
# only validates and reports, never corrects on its own.
#
#   1. source(here::here("workflow/CELESTE/run_celeste.R"))     # sites, group_manifest, config, output_dir, cache_dir, loi_layers
#   2. source(here::here("workflow/CELESTE/fix_lake_bisection.R")) # runs validate_catchment_lake_intersections_by_group(), prints/writes the report, defines correct_lake_bisected_sites_celeste() -- corrects NOTHING yet
#   3. Review output/CELESTE_engine/lake_intersection_report.csv AND the
#      catchments in QGIS against NHN waterbodies.
#   4. correct_lake_bisected_sites_celeste(
#        site_ids = c("CF1", "SN1UP"),   # can span multiple groups
#        sites = sites, group_manifest = group_manifest, config = config,
#        output_dir = output_dir, cache_dir = cache_dir, loi_layers = loi_layers
#      )
#
# Safe to call correct_lake_bisected_sites_celeste() repeatedly — same
# accumulation guarantee as CAM: prepare_lake_corrected_flow_pointer()
# persists its lake set per (group) cache_dir, so never delete
# <group_cache_dir>/lake_corrected/ between passes.
# =============================================================================

source(here::here("workflow/R/lake_containment.R"))
source(here::here("workflow/R/stream/burn_streams.R")) # find_nhn_sheets(), read_nhn_from_gdb()

NHN_DIR   <- "/Users/sam/Documents/cfs/shared_data/raw/hydro/networks/NHN/gdb"
NHN_INDEX <- "/Users/sam/Documents/cfs/shared_data/raw/hydro/networks/NHN/NHN_INDEX_WORKUNIT_LIMIT_2/NHN_INDEX_22_INDEX_WORKUNIT_LIMIT_2.shp"
nhn_idx   <- sf::st_read(NHN_INDEX, quiet = TRUE)

#' Fetch a HydroBasins group's NHN lake/reservoir polygons, unclipped
#'
#' @param group_catchments_sf sf object, that group's site catchments
#'   (used only to build a query AOI — union of their bbox, buffered)
#' @param group_id             Character. Group identifier (log messages)
#' @return sf object with OGF_ID, OFFICIAL_N, WATERBODY_ columns (renamed
#'   from NHN's nid/lakeName1/waterDefinitionText), in group_catchments_sf's
#'   CRS. 0-row if nothing found.
fetch_nhn_lakes_for_group <- function(group_catchments_sf, group_id) {
  aoi_idx_crs <- sf::st_transform(
    sf::st_as_sfc(sf::st_bbox(group_catchments_sf)), sf::st_crs(nhn_idx)
  ) |> sf::st_buffer(0.05) # ~5 km in degrees, generous margin for sheet lookup

  sheets <- find_nhn_sheets(aoi_idx_crs, nhn_idx, NHN_DIR, group_id)
  if (length(sheets) == 0) {
    return(NULL)
  }

  # read_nhn_from_gdb() reads + reprojects to EPSG:3979 with NO clipping —
  # deliberately not read_merge_nhn_layer(), see this file's header.
  #
  # nid is coerced to character immediately, per sheet, before any
  # bind_rows() — confirmed empirically that NHN GDBs are NOT consistent
  # about nid's column type across sheets (some read as character, at
  # least one as numeric), which otherwise breaks bind_rows() either here
  # (combining sheets within a group) or in
  # validate_catchment_lake_intersections_by_group()'s own combine (across
  # groups' already-renamed OGF_ID column). nid is only ever used for
  # identity/membership (%in%, ==), never arithmetic, so character is the
  # right type regardless.
  lakes <- purrr::map(sheets, function(gdb) {
    read_nhn_from_gdb(gdb, "NHN_HD_WATERBODY_2", group_id) |>
      dplyr::mutate(nid = as.character(nid))
  }) |>
    purrr::compact() |>
    dplyr::bind_rows()

  if (nrow(lakes) == 0) {
    return(NULL)
  }

  lakes <- dplyr::distinct(lakes, nid, .keep_all = TRUE) |>
    sf::st_transform(sf::st_crs(group_catchments_sf)) |>
    sf::st_make_valid()
  lakes <- lakes[!sf::st_is_empty(lakes), ]

  # Non-destructive row filter to the group's own bbox (performance only —
  # geometry itself is untouched, unlike read_merge_nhn_layer()'s clip).
  bbox_poly <- sf::st_as_sfc(sf::st_bbox(group_catchments_sf) + c(-2000, -2000, 2000, 2000))
  sf::st_crs(bbox_poly) <- sf::st_crs(group_catchments_sf)
  lakes <- sf::st_filter(lakes, bbox_poly)

  lakes |>
    dplyr::rename(OGF_ID = nid, OFFICIAL_N = lakeName1, WATERBODY_ = waterDefinitionText)
}

#' Load a group's site catchments (shared helper for validation +
#' correction, both need "this group's catchments as one sf object")
load_group_catchments <- function(sites, output_dir, group_id, catchment_file = "catchment.gpkg") {
  grp_sites <- dplyr::filter(sites, group_id == !!group_id)
  purrr::map(grp_sites$site_id, function(sid) {
    p <- fs::path(site_output_dir(output_dir, sid), catchment_file)
    if (!cache_exists(p)) {
      return(NULL)
    }
    x <- sf::st_read(p, quiet = TRUE)
    x$site_id <- sid
    x[, "site_id"]
  }) |>
    purrr::compact() |>
    dplyr::bind_rows() |>
    sf::st_make_valid()
}

# -- Validate: always runs, corrects nothing ---------------------------------

flagged <- validate_catchment_lake_intersections_by_group(
  sites = sites, output_dir = output_dir, fetch_lakes_fn = fetch_nhn_lakes_for_group,
  catchment_file = "catchment.gpkg", exclude_waterbody_types = c("Watercourse")
)

if (nrow(flagged) == 0) {
  cw_inform("No bisected lakes found.")
} else {
  cw_inform(glue::glue(
    "{nrow(flagged)} bisected (site, lake) pair(s) across ",
    "{dplyr::n_distinct(flagged$site_id)} site(s) found. Review ",
    "{fs::path(output_dir, 'lake_intersection_report.csv')} and the ",
    "catchments in QGIS, then call correct_lake_bisected_sites_celeste(site_ids = c(...)) ",
    "for the ones confirmed to need fixing."
  ))
}

#' Correct catchments for a manually-confirmed list of site IDs, across
#' however many HydroBasins groups they span
#'
#' Re-validates fresh, filters to the requested site_ids, splits by
#' group_id (since each group has its own D8 pointer grid), and for each
#' group represented: fetches that group's flagged lakes' full geometries,
#' builds/accumulates that group's corrected pointer, and redelineates just
#' that group's target sites via rerun_engine_sites() — one call per group,
#' not one call for the whole cross-group batch.
#'
#' @param site_ids       Character vector. The site_id(s) to correct.
#' @param sites          Validated sites tibble (site_id, group_id, ...)
#' @param group_manifest sf tibble from build_engine_group_manifest()
#' @param config         Resolved config from resolve_engine_config()
#' @param output_dir     Character. Root output directory
#' @param cache_dir      Character. Project-level cache directory (NOT a
#'   group cache dir — passed through to rerun_engine_sites() unchanged for
#'   its hydroweight call, which resolves per-group internally)
#' @param loi_layers     LOI layers list for hydroweight, or NULL to skip
#'
#' @return Invisibly, the before/after summary tibble
correct_lake_bisected_sites_celeste <- function(
  site_ids, sites, group_manifest, config, output_dir, cache_dir, loi_layers = NULL
) {
  unknown_ids <- setdiff(site_ids, sites$site_id)
  if (length(unknown_ids) > 0) {
    cw_abort(glue::glue("Unknown site_id(s): {paste(unknown_ids, collapse = ', ')}"))
  }

  current_flagged <- validate_catchment_lake_intersections_by_group(
    sites = sites, output_dir = output_dir, fetch_lakes_fn = fetch_nhn_lakes_for_group,
    catchment_file = "catchment.gpkg", exclude_waterbody_types = c("Watercourse")
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
  target_groups <- sites |>
    dplyr::filter(site_id %in% target_ids) |>
    dplyr::distinct(site_id, group_id)
  to_fix <- dplyr::left_join(to_fix, target_groups, by = "site_id")

  cw_inform(glue::glue(
    "Correcting {length(target_ids)} site(s) across ",
    "{dplyr::n_distinct(to_fix$group_id)} group(s): {paste(target_ids, collapse = ', ')} ",
    "({nrow(to_fix)} bisected (site, lake) pair(s))."
  ))

  before_areas <- purrr::map_dfr(target_ids, function(sid) {
    p <- fs::path(site_output_dir(output_dir, sid), "catchment.gpkg")
    tibble::tibble(
      site_id = sid,
      area_km2_before = as.numeric(sf::st_area(sf::st_read(p, quiet = TRUE))) / 1e6
    )
  })

  for (grp in unique(to_fix$group_id)) {
    grp_target_ids <- unique(to_fix$site_id[to_fix$group_id == grp])
    grp_lake_ids <- unique(to_fix$OGF_ID[to_fix$group_id == grp])
    grp_cache <- group_manifest$cache_dir[group_manifest$group_id == grp]

    grp_catchments <- load_group_catchments(sites, output_dir, grp)
    grp_lakes_all <- fetch_nhn_lakes_for_group(grp_catchments, grp)
    grp_lake_polys <- grp_lakes_all[grp_lakes_all$OGF_ID %in% grp_lake_ids, ]

    cw_inform(glue::glue(
      "Group '{grp}': correcting {length(grp_target_ids)} site(s), ",
      "flattening {nrow(grp_lake_polys)} lake(s)."
    ))

    corrected_pointer <- prepare_lake_corrected_flow_pointer(
      cache_dir = grp_cache, lake_polys = grp_lake_polys
    )

    rerun_engine_sites(
      edited_snap_site_ids = grp_target_ids,
      sites = sites, group_manifest = group_manifest, config = config,
      output_dir = output_dir, cache_dir = cache_dir, loi_layers = loi_layers,
      flow_pointer_override_path = corrected_pointer
    )
  }

  after_unclipped <- validate_catchment_lake_intersections_by_group(
    sites = sites, output_dir = output_dir, fetch_lakes_fn = fetch_nhn_lakes_for_group,
    catchment_file = "catchment.gpkg", exclude_waterbody_types = c("Watercourse")
  )
  after_clipped <- validate_catchment_lake_intersections_by_group(
    sites = sites, output_dir = output_dir, fetch_lakes_fn = fetch_nhn_lakes_for_group,
    catchment_file = "catchment_clipped.gpkg", exclude_waterbody_types = c("Watercourse")
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
