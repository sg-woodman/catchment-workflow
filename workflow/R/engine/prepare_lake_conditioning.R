# prepare_lake_conditioning.R
# ---------------------------------------------------------------------------
# Resolves config$lake_conditioning for the raw-dem terrain tier only
# (flattening has to happen before D8 derivation is finalized — see
# 00_resolve_config.R's tier-vs-lake_conditioning validation, deliberately
# mirroring streams_burn's own restriction). Called from
# 02_prepare_terrain.R::prepare_raw_dem_tier(), after breaching, before D8
# pointer derivation.
#
# WHY THIS EXISTS: reactively detecting and correcting individual bisected
# catchments after delineation (workflow/R/lake_containment.R,
# workflow/CELESTE/fix_lake_bisection.R) works, but is fixing a symptom.
# Correcting one site's D8 pointer without touching an overlapping
# neighbor's breaks remove_upstream_catchments()'s implicit assumption that
# every site in a group is traced from the SAME flow field, so nested
# catchments are always cleanly nested — confirmed directly to silently
# corrupt output (workflow/R/remove_upstream.R's integrity check, added
# after finding a 3-part disconnected catchment and a clip that excised its
# own pour point). Conditioning known lakes into the DEM up front, as a
# normal part of terrain preparation, means every site in a group shares
# one internally-consistent, lake-aware flow field from the start — the
# whole bug class this file exists to prevent can't occur.
#
#   source = "nhn_auto"  — fetch_nhn_lakes_for_aoi() below, generalizing
#     workflow/CELESTE/fix_lake_bisection.R's fetch_nhn_lakes_for_group()
#     to key off the group's AOI directly (group_manifest$aoi) instead of
#     already-delineated catchments, since at this point in the pipeline
#     (Stage 2, terrain prep) no catchment exists yet.
#   source = "supplied"  — reads config$lake_conditioning$path directly,
#     clips to the group AOI (any vector format sf::st_read() handles).
#   source = "none"      — never reaches this file (00_resolve_config.R
#     only allows lake_conditioning to matter for the dem tier, and
#     prepare_raw_dem_tier() only calls resolve_lake_conditioning() when
#     source != "none").
#
# ORDERING, deliberate — see 02_prepare_terrain.R's prepare_raw_dem_tier():
# breach runs on the real terrain FIRST (preserves the existing validated
# breach-based behavior for every non-lake depression exactly as before
# this feature existed), THEN lakes are flattened (overwriting whatever
# breach did to those cells, which is fine — that elevation was noise
# anyway), THEN wbt_fill_depressions(fix_flats = TRUE) resolves the flat
# plateau flattening just created. fill (not breach) is required for that
# last step specifically — confirmed directly (see workflow/R/
# lake_containment.R's header) that breaching a freshly-flattened lake
# collapsed a catchment to a single cell, since breach has no local lower
# cell to carve toward across a genuinely flat region.
#
# Requires workflow/R/stream/burn_streams.R already sourced for
# find_nhn_sheets()/read_nhn_from_gdb() (reused unmodified — same NHN
# reading primitives streams_burn's nhn_auto path and the reactive
# CELESTE fix both already use).
#
# Dependencies: sf, terra, whitebox, dplyr, purrr, fs, glue, cli (via
# utils.R)
# ---------------------------------------------------------------------------

#' Condition known lakes into a group's (already-breached) DEM
#'
#' Deliberately scopes the lake search to a buffer around the group's own
#' SITES, not the group's full terrain-conditioning AOI — confirmed
#' directly (2026-08-29) that the latter is wildly disproportionate: for
#' CELESTE's KEN group, the HydroBasins group AOI is 58,928 sq km (the
#' group's 21 sites' own bare bounding box is 2,663 sq km, ~22x smaller),
#' and a first attempt scoped to the full group AOI pulled in 31,525
#' candidate lakes and took ~68 minutes for one group. A lake far from
#' every site in the group can't affect any site's catchment or its
#' nesting relationship with a neighbor — the only property this feature
#' needs to preserve — so there's no correctness reason to flatten it.
#'
#' @param config          Resolved config from resolve_engine_config()
#'   (lake_conditioning$source/min_area_ha/exclude_types/site_buffer_m
#'   already defaulted; config$sites has lon/lat/group_id for every site)
#' @param input_dem_path  Character. Path to the breached DEM to condition
#'   (prepare_raw_dem_tier()'s breach output — NOT the raw dem.tif)
#' @param output_dem_path Character. Output path — prepare_raw_dem_tier()
#'   passes its final dem_breached.tif here, so every downstream consumer
#'   of that canonical filename (D8 pointer, hillshade, hydroweight) sees
#'   the fully lake-conditioned surface with no changes needed on their end
#' @param group_id        Character. Group identifier (log messages)
#' @return Invisibly, output_dem_path
resolve_lake_conditioning <- function(config, input_dem_path, output_dem_path, group_id) {
  if (!exists("find_nhn_sheets", mode = "function")) {
    cw_abort(paste(
      "resolve_lake_conditioning() requires workflow/R/stream/burn_streams.R",
      "to be sourced first (defines find_nhn_sheets()/read_nhn_from_gdb())."
    ))
  }

  min_area_ha <- config$lake_conditioning$min_area_ha
  exclude_types <- config$lake_conditioning$exclude_types

  fetch_aoi_sf <- build_site_buffer_aoi(
    sites = config$sites, group_id = group_id, working_crs = config$working_crs,
    buffer_m = config$lake_conditioning$site_buffer_m
  )
  cw_inform(glue::glue(
    "Group '{group_id}': lake search scoped to a {config$lake_conditioning$site_buffer_m / 1000} km ",
    "buffer around its own sites ({round(as.numeric(sf::st_area(fetch_aoi_sf)) / 1e6)} sq km) — ",
    "not the full group terrain AOI."
  ))

  lakes <- switch(
    config$lake_conditioning$source,
    nhn_auto = fetch_nhn_lakes_for_aoi(
      aoi = fetch_aoi_sf, nhn_index_path = config$nhn_index_path,
      nhn_raw_dir = config$nhn_raw_dir, group_id = group_id,
      min_area_ha = min_area_ha, exclude_types = exclude_types
    ),
    supplied = {
      cw_inform(glue::glue("Group '{group_id}': reading supplied lake polygons..."))
      raw <- sf::st_read(config$lake_conditioning$path, quiet = TRUE) |>
        sf::st_transform(sf::st_crs(fetch_aoi_sf))
      clipped <- tryCatch(
        sf::st_filter(raw, fetch_aoi_sf),
        error = function(e) {
          cw_warn(glue::glue(
            "Group '{group_id}': error filtering supplied lakes — {e$message}"
          ))
          raw[0, ]
        }
      )
      areas_ha <- as.numeric(sf::st_area(clipped)) / 1e4
      clipped[areas_ha >= min_area_ha, ]
    }
  )

  if (is.null(lakes) || nrow(lakes) == 0) {
    cw_inform(glue::glue(
      "Group '{group_id}': no lakes found for conditioning — ",
      "using breached DEM as-is."
    ))
    fs::file_copy(input_dem_path, output_dem_path, overwrite = TRUE)
    return(invisible(output_dem_path))
  }

  cw_inform(glue::glue(
    "Group '{group_id}': conditioning {nrow(lakes)} lake(s) into the DEM..."
  ))

  out_dir <- fs::path_dir(output_dem_path)

  # WhiteboxTools' vector reader requires Shapefile for FlattenLakes
  # (confirmed empirically — a .gpkg input crashes with "Unrecognized
  # ShapeType" — see workflow/R/lake_containment.R's
  # prepare_lake_corrected_flow_pointer() for the same finding). Disposable
  # — attribute values don't need to survive this round-trip, only geometry
  # does, so dissolved into one (possibly multi-part) feature rather than
  # written as N separate rows — FlattenLakes only needs "which cells fall
  # inside any lake", and collapsing potentially thousands of individual
  # polygons into a single feature avoids whatever per-feature overhead
  # scales with row count in WBT's own reader/algorithm.
  lakes_shp <- fs::path(out_dir, glue::glue("lakes_to_condition_{group_id}.shp"))
  lakes_dissolved <- sf::st_union(sf::st_make_valid(lakes)) |> sf::st_as_sf()
  sf::st_write(lakes_dissolved, lakes_shp, delete_dsn = TRUE, quiet = TRUE)

  flattened_path <- fs::path(out_dir, "dem_lakes_flattened.tif")
  whitebox::wbt_flatten_lakes(
    dem    = normalizePath(input_dem_path, mustWork = TRUE),
    lakes  = normalizePath(lakes_shp, mustWork = TRUE),
    output = normalizePath(flattened_path, mustWork = FALSE)
  )
  if (!cache_exists(flattened_path)) {
    cw_abort(glue::glue("Group '{group_id}': wbt_flatten_lakes() did not produce output."))
  }

  cw_inform(glue::glue("Group '{group_id}': filling depressions (fix_flats = TRUE)..."))
  whitebox::wbt_fill_depressions(
    dem       = normalizePath(flattened_path, mustWork = TRUE),
    output    = normalizePath(output_dem_path, mustWork = FALSE),
    fix_flats = TRUE
  )
  if (!cache_exists(output_dem_path)) {
    cw_abort(glue::glue("Group '{group_id}': wbt_fill_depressions() did not produce output."))
  }

  invisible(output_dem_path)
}

#' Build a buffered AOI around a group's own sites (not its full
#' terrain-conditioning AOI)
#'
#' @param sites      tibble with site_id, lon, lat, group_id columns
#'   (config$sites — lon/lat in EPSG:4326)
#' @param group_id   Character. Group identifier
#' @param working_crs Character. "EPSG:####" to build the AOI in
#' @param buffer_m   Numeric. Buffer distance in metres around the sites'
#'   own bounding box
#' @return sf polygon, in working_crs
build_site_buffer_aoi <- function(sites, group_id, working_crs, buffer_m) {
  grp_sites <- dplyr::filter(sites, group_id == !!group_id)
  pts <- sf::st_as_sf(grp_sites, coords = c("lon", "lat"), crs = 4326) |>
    sf::st_transform(working_crs)
  sf::st_as_sfc(sf::st_bbox(pts)) |>
    sf::st_buffer(buffer_m) |>
    sf::st_as_sf()
}

#' Fetch NHN lake/reservoir polygons for a group's AOI, unclipped
#'
#' Generalizes workflow/CELESTE/fix_lake_bisection.R's
#' fetch_nhn_lakes_for_group() to key off a group's AOI directly rather
#' than already-delineated site catchments — at terrain-prep time (Stage
#' 2) no catchment exists yet. Deliberately does NOT reuse stream/
#' burn_streams.R's read_merge_nhn_layer(): that function clips each
#' waterbody's geometry to the query AOI, correct for its original job
#' (trimming burn-in streams to a group's conditioning extent) but wrong
#' here, since it would distort which lakes get flattened and how much of
#' each — confirmed as a real accuracy bug in the reactive fix (this
#' file's sibling, fix_lake_bisection.R). Uses read_nhn_from_gdb()
#' directly (reads + reprojects with no clipping) instead, exactly like
#' the reactive fix does.
#'
#' @param aoi             sf polygon. Group AOI, in any CRS (transformed
#'   internally to the NHN index's CRS for the sheet lookup)
#' @param nhn_index_path  Character. Path to the NHN index shapefile
#' @param nhn_raw_dir     Character. Path to the directory of NHN GDB
#'   sheet folders
#' @param group_id        Character. Group identifier (log messages)
#' @param min_area_ha     Numeric. Minimum lake area (ha) to include
#' @param exclude_types    Character vector of waterDefinitionText values
#'   to exclude (e.g. "Watercourse")
#' @return sf object with OGF_ID, OFFICIAL_N, WATERBODY_ columns (renamed
#'   from NHN's nid/lakeName1/waterDefinitionText, same normalization the
#'   reactive fix uses), in aoi's CRS. NULL if nothing found.
fetch_nhn_lakes_for_aoi <- function(aoi, nhn_index_path, nhn_raw_dir, group_id,
                                     min_area_ha = 1, exclude_types = c("Watercourse")) {
  nhn_idx <- sf::st_read(nhn_index_path, quiet = TRUE)
  aoi_idx_crs <- sf::st_transform(aoi, sf::st_crs(nhn_idx))

  sheets <- find_nhn_sheets(aoi_idx_crs, nhn_idx, nhn_raw_dir, group_id)
  if (length(sheets) == 0) {
    return(NULL)
  }

  # nid is coerced to character immediately, per sheet, before any
  # bind_rows() — NHN GDBs are not consistent about nid's column type
  # across sheets (confirmed empirically: some read as character, at
  # least one as numeric), which otherwise breaks bind_rows() combining
  # sheets. nid is only ever used for identity/dedup, never arithmetic.
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
    sf::st_transform(sf::st_crs(aoi)) |>
    sf::st_make_valid()
  lakes <- lakes[!sf::st_is_empty(lakes), ]

  # Non-destructive row filter to the group's own bbox (performance only —
  # geometry itself is untouched).
  bbox_poly <- sf::st_as_sfc(sf::st_bbox(aoi) + c(-2000, -2000, 2000, 2000))
  sf::st_crs(bbox_poly) <- sf::st_crs(aoi)
  lakes <- sf::st_filter(lakes, bbox_poly)

  if (nrow(lakes) == 0) {
    return(NULL)
  }

  areas_ha <- as.numeric(sf::st_area(lakes)) / 1e4
  lakes <- lakes[areas_ha >= min_area_ha, ]
  lakes <- lakes[is.na(lakes$waterDefinitionText) | !lakes$waterDefinitionText %in% exclude_types, ]

  if (nrow(lakes) == 0) {
    return(NULL)
  }

  lakes |>
    dplyr::rename(OGF_ID = nid, OFFICIAL_N = lakeName1, WATERBODY_ = waterDefinitionText)
}
