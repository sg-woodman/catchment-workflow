# lake_containment.R
# ---------------------------------------------------------------------------
# Shared "does this catchment fully contain this lake" primitives, used by
# both project kinds:
#
#   - validate_lake_containment() — CAM lakes (lake-polygon pour point):
#     each site has exactly one designated lake (config$lake_polygons'
#     matched_lake join key). Currently dormant/unwired — a scan of all 45
#     CAM lakes sites found 0 bisected, so there's nothing to fix there
#     today, but the check is cheap QA infrastructure worth keeping.
#
#   - validate_catchment_lake_intersections() — point-based projects (CAM
#     streams; potentially CELESTE later) with no single designated lake per
#     site: spatially joins each catchment against a project-wide waterbody
#     layer and flags any intersecting lake that's only partially contained.
#     This is the active case — see workflow/CAM/fix_lake_bisection.R.
#
# Root cause (both cases): wbt_watershed() traces strictly by per-cell D8
# flow direction. Nothing guarantees a lake polygon is flow-consistent (all
# interior cells routing to one outlet) unless something explicitly enforces
# it — OIH's own hydro-enforcement only "improves flow direction" for
# waterbodies intersecting the mapped Enhanced Watercourse network, not
# universally (OIH User Guide, Hydro Enforcement Overview). Where that
# doesn't hold, a traced watershed boundary can slice through a lake instead
# of being a strict superset (or fully excluding it).
#
# Dependencies: sf, dplyr, purrr, tibble, fs, readr, glue (via utils.R)
# ---------------------------------------------------------------------------

#' Percent of a lake polygon's area falling outside a catchment polygon
#'
#' @param lake_sf       sf object, single lake polygon (any CRS)
#' @param catchment_sf  sf object, single catchment polygon (any CRS)
#' @return Numeric, 0-100. 0 = catchment fully contains the lake.
lake_containment_pct <- function(lake_sf, catchment_sf) {
  catchment_sf <- sf::st_transform(catchment_sf, sf::st_crs(lake_sf))
  lake_sf <- sf::st_make_valid(lake_sf)
  catchment_sf <- sf::st_make_valid(catchment_sf)

  lake_area <- as.numeric(sf::st_area(sf::st_union(lake_sf)))
  if (lake_area == 0) {
    return(NA_real_)
  }

  outside <- sf::st_difference(sf::st_union(lake_sf), sf::st_union(catchment_sf))
  outside_area <- if (length(outside) == 0 || sf::st_is_empty(outside)) {
    0
  } else {
    as.numeric(sf::st_area(outside))
  }

  100 * outside_area / lake_area
}

#' Validate lake containment for every site against its own designated lake
#'
#' CAM lakes case: each site has exactly one matched lake
#' (config$lake_polygons' matched_lake column, keyed to sites$lake_name).
#'
#' @param sites          tibble with site_id, lake_name columns
#' @param lake_polys     SpatVector or sf from match_lake_polygons(), with a
#'                       matched_lake column keyed to sites$lake_name
#' @param output_dir     Character. Root output directory (reads
#'                       <site_id>/catchment.gpkg, writes
#'                       lake_containment_report.csv here)
#' @param threshold_pct  Numeric. % of lake area outside the catchment above
#'                       which a site is flagged as bisecting its lake.
#'                       Default 1 (shoreline tolerance at 30 m resolution).
#'
#' @return A tibble with columns: site_id, lake_name, pct_lake_outside,
#'   bisects_lake, status
validate_lake_containment <- function(sites, lake_polys, output_dir, threshold_pct = 1) {
  if (inherits(lake_polys, "SpatVector")) {
    lake_polys <- sf::st_as_sf(lake_polys)
  }

  results <- purrr::map(seq_len(nrow(sites)), function(i) {
    sid <- sites$site_id[i]
    lake_name <- sites$lake_name[i]
    catchment_path <- fs::path(site_output_dir(output_dir, sid), "catchment.gpkg")

    if (!cache_exists(catchment_path)) {
      return(tibble::tibble(
        site_id = sid, lake_name = lake_name, pct_lake_outside = NA_real_,
        bisects_lake = NA, status = "skipped (no catchment.gpkg)"
      ))
    }

    lake_sf <- lake_polys[lake_polys$matched_lake == lake_name, ]
    if (nrow(lake_sf) == 0) {
      return(tibble::tibble(
        site_id = sid, lake_name = lake_name, pct_lake_outside = NA_real_,
        bisects_lake = NA, status = "skipped (no matched lake polygon)"
      ))
    }
    if (nrow(lake_sf) > 1) {
      lake_sf <- lake_sf[1, ]
    }

    catchment_sf <- sf::st_read(catchment_path, quiet = TRUE)

    pct_outside <- tryCatch(
      lake_containment_pct(lake_sf, catchment_sf),
      error = function(e) {
        cw_warn(glue::glue("Site '{sid}': containment check failed — {e$message}"))
        NA_real_
      }
    )

    bisects <- !is.na(pct_outside) && pct_outside > threshold_pct
    if (bisects) {
      cw_warn(glue::glue(
        "Site '{sid}': catchment bisects lake '{lake_name}' — ",
        "{round(pct_outside, 2)}% of lake area outside catchment"
      ))
    }

    tibble::tibble(
      site_id = sid, lake_name = lake_name, pct_lake_outside = pct_outside,
      bisects_lake = bisects, status = "checked"
    )
  }) |>
    dplyr::bind_rows()

  report_path <- fs::path(output_dir, "lake_containment_report.csv")
  readr::write_csv(results, report_path)

  n_bisecting <- sum(results$bisects_lake, na.rm = TRUE)
  cw_inform(glue::glue(
    "Lake containment check: {n_bisecting} of {nrow(results)} site(s) bisect ",
    "their lake (threshold {threshold_pct}%). Report written to {report_path}."
  ))

  results
}

#' Validate every catchment against every lake it intersects, for projects
#' with no single designated lake per site (point pour points)
#'
#' Reads the full lakes_path layer WITHOUT wkt_filter — GDAL's spatial
#' pushdown is pushed to OGR before any R-side geometry cleanup and has been
#' observed to silently drop features on this kind of source (see
#' rasterize_competing_classes()'s docstring in workflow/raster_attributes.R
#' for the confirmed New Brunswick case). Filters to the sites' combined
#' bbox in R instead via sf::st_filter(), the same fix pattern.
#'
#' @param sites                    tibble with site_id column
#' @param output_dir               Character. Root output directory (reads
#'   <site_id>/<catchment_file>, writes lake_intersection_report.csv here)
#' @param lakes_path                Character. Path to OHN/OIH waterbody layer
#' @param catchment_file           Character. Which per-site file to check —
#'   "catchment.gpkg" (root-cause check, before any correction) or
#'   "catchment_clipped.gpkg" (post-correction final QA, also catches
#'   remove_upstream_catchments()'s erase step reintroducing a bisection).
#'   Default "catchment.gpkg".
#' @param min_lake_area_ha         Numeric. Minimum lake area (ha) to
#'   consider — smaller features are typically shoreline-generalization
#'   slivers or wetland fragments, not real lakes. Default 1.
#' @param exclude_waterbody_types  Character vector of WATERBODY_ values to
#'   exclude. Default c("River", "Pond") — matching match_lake_polygons.R's/
#'   remove_upstream_lakes.R's existing conventions.
#' @param partial_band             Numeric length-2. pct_lake_outside range
#'   considered a genuine partial bisection; outside this range is treated
#'   as boundary-touch noise (near 0%) or "lake legitimately not part of
#'   this catchment" (near 100%) rather than a bug. Default c(2, 98).
#'
#' @return A long tibble (site_id, OGF_ID, lake_name, lake_area_ha,
#'   pct_lake_outside), one row per flagged (site, lake) pair
validate_catchment_lake_intersections <- function(
  sites, output_dir, lakes_path, catchment_file = "catchment.gpkg",
  min_lake_area_ha = 1, exclude_waterbody_types = c("River", "Pond"),
  partial_band = c(2, 98)
) {
  catchments <- purrr::map(sites$site_id, function(sid) {
    p <- fs::path(site_output_dir(output_dir, sid), catchment_file)
    if (!cache_exists(p)) {
      return(NULL)
    }
    x <- sf::st_read(p, quiet = TRUE)
    x$site_id <- sid
    x[, "site_id"]
  }) |>
    purrr::compact() |>
    dplyr::bind_rows()

  if (nrow(catchments) == 0) {
    cw_warn(glue::glue("No {catchment_file} files found under {output_dir} — nothing to validate."))
    empty <- tibble::tibble(
      site_id = character(0), OGF_ID = numeric(0), lake_name = character(0),
      lake_area_ha = numeric(0), pct_lake_outside = numeric(0)
    )
    readr::write_csv(empty, fs::path(output_dir, "lake_intersection_report.csv"))
    return(empty)
  }

  cw_inform(glue::glue("Validating lake intersections for {nrow(catchments)} site(s) [{catchment_file}]..."))

  catchments <- sf::st_make_valid(catchments)
  bbox_poly <- sf::st_as_sfc(sf::st_bbox(catchments) + c(-2000, -2000, 2000, 2000))
  sf::st_crs(bbox_poly) <- sf::st_crs(catchments)

  lakes_full <- sf::st_read(lakes_path, quiet = TRUE) |>
    sf::st_transform(sf::st_crs(catchments))
  lakes <- sf::st_filter(lakes_full, bbox_poly) |>
    sf::st_make_valid()
  lakes <- lakes[!sf::st_is_empty(lakes), ]
  cw_inform(glue::glue("{nrow(lakes)} candidate lake(s) within site bounding area."))

  results <- purrr::map(seq_len(nrow(catchments)), function(i) {
    sid <- catchments$site_id[i]
    cat_geom <- catchments[i, ]
    hits <- sf::st_intersects(lakes, cat_geom, sparse = FALSE)[, 1]
    if (sum(hits) == 0) {
      return(NULL)
    }
    cand <- lakes[hits, ]
    purrr::map(seq_len(nrow(cand)), function(j) {
      lake_row <- cand[j, ]
      lake_area_ha <- as.numeric(sf::st_area(lake_row)) / 1e4
      if (lake_area_ha < min_lake_area_ha) {
        return(NULL)
      }
      if (!is.na(lake_row$WATERBODY_) && lake_row$WATERBODY_ %in% exclude_waterbody_types) {
        return(NULL)
      }
      pct_outside <- tryCatch(
        lake_containment_pct(lake_row, cat_geom),
        error = function(e) NA_real_
      )
      if (is.na(pct_outside) || pct_outside <= partial_band[1] || pct_outside >= partial_band[2]) {
        return(NULL)
      }
      tibble::tibble(
        site_id = sid, OGF_ID = lake_row$OGF_ID, lake_name = lake_row$OFFICIAL_N,
        lake_area_ha = lake_area_ha, pct_lake_outside = pct_outside
      )
    }) |>
      purrr::compact() |>
      dplyr::bind_rows()
  }) |>
    purrr::compact() |>
    dplyr::bind_rows()

  if (nrow(results) == 0) {
    results <- tibble::tibble(
      site_id = character(0), OGF_ID = numeric(0), lake_name = character(0),
      lake_area_ha = numeric(0), pct_lake_outside = numeric(0)
    )
  } else {
    results <- dplyr::arrange(results, dplyr::desc(lake_area_ha))
  }

  report_path <- fs::path(output_dir, "lake_intersection_report.csv")
  readr::write_csv(results, report_path)

  cw_inform(glue::glue(
    "Lake intersection check: {nrow(results)} bisected (site, lake) pair(s) ",
    "across {dplyr::n_distinct(results$site_id)} site(s). ",
    "Report written to {report_path}."
  ))
  if (nrow(results) > 0) {
    print(results)
  }

  results
}

#' Build a lake-flattened, fill-and-D8-corrected flow pointer for a
#' point-based project's group, to redelineate sites whose catchment
#' bisects a lake
#'
#' Reuses cache_dir/dem.tif directly rather than re-cropping/reprojecting
#' from the raw DEM source — for a "flow_direction" terrain tier
#' (engine/02_prepare_terrain.R's prepare_flow_direction_tier()), that file
#' is already the project's DEM cropped + reprojected to the exact same
#' extent/CRS/alignment as flow_pointer.tif, so this guarantees the
#' corrected pointer aligns with the one it's replacing without duplicating
#' any crop/reproject logic.
#'
#' Follows FlattenLakes with wbt_fill_depressions(fix_flats = TRUE), NOT
#' engine/02_prepare_terrain.R's prepare_raw_dem_tier() breach step
#' (wbt_breach_depressions_least_cost()) — confirmed by direct testing to
#' matter: breach is built for small, localized pits (it carves a channel
#' through them) and does not resolve the large, genuinely-flat region
#' FlattenLakes creates. Run against Daisy Lake, breach produced a
#' catchment collapsed to a single ~900 sq m cell (the pour point landing in
#' the flattened lake's now-ambiguous flat with no assigned downstream
#' direction). fill_depressions(fix_flats = TRUE) imposes the tiny gradient
#' standard practice uses to guarantee flow converges to one outlet across a
#' flat — the standard FlattenLakes -> FillDepressions -> D8Pointer sequence.
#'
#' ACCUMULATES across repeated calls, keyed by lake_polys' OGF_ID column,
#' via the persisted lakes_to_flatten.shp itself: a lake flattened in an
#' earlier call stays flattened in every later rebuild, even if that later
#' call's lake_polys doesn't include it. This matters because a corrected
#' pointer is shared by every site redelineated against it — a second
#' correction pass that only knew about its own newly-flagged lakes was
#' confirmed to silently UN-flatten the first pass's lakes for any site
#' redelineated in that second pass, reverting them back toward their
#' original bisected state. Because of this, the caller should never delete
#' output_subdir's contents to "force a rebuild" — the function detects a
#' grown lake set on its own (comparing OGF_IDs) and rebuilds automatically;
#' deleting the directory would just discard the accumulated history and
#' reintroduce the bug above.
#'
#' @param cache_dir     Character. Project cache directory (must contain
#'   dem.tif, already cropped/reprojected by the engine's terrain prep)
#' @param lake_polys    sf object with an OGF_ID column. The lake polygon(s)
#'   to flatten (typically just the ones flagged by
#'   validate_catchment_lake_intersections()) — merged with, not replacing,
#'   whatever's already accumulated from previous calls.
#' @param output_subdir Character. Subdirectory of cache_dir to write
#'   corrected products into. Default "lake_corrected".
#'
#' @return Character. Path to the corrected flow_pointer.tif.
prepare_lake_corrected_flow_pointer <- function(cache_dir, lake_polys, output_subdir = "lake_corrected") {
  dem_path <- fs::path(cache_dir, "dem.tif")
  if (!cache_exists(dem_path)) {
    cw_abort(glue::glue("prepare_lake_corrected_flow_pointer(): {dem_path} not found."))
  }

  out_dir <- fs::path(cache_dir, output_subdir)
  ensure_dir(out_dir)
  pointer_path <- fs::path(out_dir, "flow_pointer.tif")

  # The accumulation manifest is a GPKG, not a Shapefile — confirmed
  # empirically that writing OGF_ID (9-digit integers) to a .shp's .dbf
  # triggers a GDAL "not successfully written, possibly due to too large a
  # number" warning; it happened to round-trip correctly in one test, but
  # that's an unreliable numeric-field-width edge case to depend on for
  # the exact-match logic below, not something to build accumulation
  # correctness on. WhiteboxTools' own vector reader is the thing that
  # requires Shapefile specifically (confirmed empirically: a .gpkg input
  # crashes FlattenLakes with "Unrecognized ShapeType" — its shapefile
  # reader is being handed GPKG bytes) — so a separate, disposable .shp
  # (attributes don't matter for that call, only geometry) is written just
  # before each wbt_flatten_lakes() call, from whatever the accumulated
  # GPKG says at that point. Same Shapefile constraint
  # workflow/R/engine/04_delineate_site.R's write_pour_point_shp_dynamic()/
  # stream/delineate_sites.R's write_pour_point_shp() already work around
  # for wbt_watershed()'s pour_pts argument.
  lakes_gpkg <- fs::path(out_dir, "lakes_to_flatten.gpkg")
  lakes_shp  <- fs::path(out_dir, "lakes_to_flatten_wbt_input.shp")

  dem_rast <- terra::rast(dem_path)
  lake_polys_proj <- sf::st_transform(lake_polys, terra::crs(dem_rast)) |>
    sf::st_make_valid()

  if (cache_exists(lakes_gpkg)) {
    previous <- sf::st_read(lakes_gpkg, quiet = TRUE)
    new_ids <- setdiff(lake_polys_proj$OGF_ID, previous$OGF_ID)
    if (length(new_ids) == 0 && cache_exists(pointer_path)) {
      cw_inform(glue::glue(
        "Lake-corrected flow pointer already covers all {nrow(lake_polys_proj)} ",
        "requested lake(s) — skipping: {pointer_path}"
      ))
      return(pointer_path)
    }
    # Reduce both to just OGF_ID + geometry and force a common geometry
    # column name before combining — lake_polys_proj's geometry column name
    # depends on lakes_path's own schema, and may not match previous's;
    # dplyr::bind_rows()/rbind.sf() need matching names to align.
    previous_min <- previous["OGF_ID"]
    new_min <- lake_polys_proj["OGF_ID"]
    sf::st_geometry(previous_min) <- "geometry"
    sf::st_geometry(new_min) <- "geometry"
    lake_polys_proj <- rbind(previous_min, new_min) |>
      dplyr::distinct(OGF_ID, .keep_all = TRUE)
    cw_inform(glue::glue(
      "Accumulating {length(new_ids)} newly-flagged lake(s) with ",
      "{nrow(previous)} previously-flattened lake(s) -> {nrow(lake_polys_proj)} total."
    ))
  } else if (cache_exists(pointer_path)) {
    cw_inform(glue::glue("Lake-corrected flow pointer found in cache — skipping: {pointer_path}"))
    return(pointer_path)
  }

  sf::st_write(lake_polys_proj, lakes_gpkg, delete_dsn = TRUE, quiet = TRUE)
  # Disposable — WBT input only, attribute values (including OGF_ID) don't
  # need to survive this round-trip, only geometry does.
  sf::st_write(lake_polys_proj, lakes_shp, delete_dsn = TRUE, quiet = TRUE)

  cw_inform(glue::glue(
    "Flattening {nrow(lake_polys_proj)} lake(s) (accumulated) into {dem_path}..."
  ))

  flattened_path <- fs::path(out_dir, "dem_lakes_flattened.tif")
  whitebox::wbt_flatten_lakes(
    dem    = normalizePath(dem_path, mustWork = TRUE),
    lakes  = normalizePath(lakes_shp, mustWork = TRUE),
    output = normalizePath(flattened_path, mustWork = FALSE)
  )
  if (!cache_exists(flattened_path)) {
    cw_abort("prepare_lake_corrected_flow_pointer(): wbt_flatten_lakes() did not produce output.")
  }

  filled_path <- fs::path(out_dir, "dem_filled.tif")
  cw_inform("Filling depressions (fix_flats = TRUE) in the lake-flattened DEM...")
  whitebox::wbt_fill_depressions(
    dem        = normalizePath(flattened_path, mustWork = TRUE),
    output     = normalizePath(filled_path, mustWork = FALSE),
    fix_flats  = TRUE
  )
  if (!cache_exists(filled_path)) {
    cw_abort("prepare_lake_corrected_flow_pointer(): wbt_fill_depressions() did not produce output.")
  }

  cw_inform("Computing corrected D8 flow pointer...")
  whitebox::wbt_d8_pointer(
    dem    = normalizePath(filled_path, mustWork = TRUE),
    output = normalizePath(pointer_path, mustWork = FALSE)
  )
  if (!cache_exists(pointer_path)) {
    cw_abort("prepare_lake_corrected_flow_pointer(): wbt_d8_pointer() did not produce output.")
  }

  cw_inform(glue::glue("Lake-corrected flow pointer written: {pointer_path}"))
  pointer_path
}
