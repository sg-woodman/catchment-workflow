# 04_remove_upstream_lakes.R
# ---------------------------------------------------------------------------
# Clips each target lake catchment to remove the contributing area of any
# upstream lakes found within the catchment boundary.
#
# Strategy:
#   1. For each target lake site, find all OHN lake polygons intersecting its
#      catchment (excluding the target lake itself and small ponds)
#   2. Delineate and cache catchments for those upstream lakes using the same
#      lake-raster pour-point approach as 03_delineate_lakes.R
#   3. Build a pool of all catchments (target sites + upstream lakes)
#   4. For each target site, erase all smaller intersecting catchments
#   5. Write catchment_clipped.gpkg per site (same convention as the stream
#      06_remove_upstream.R) so downstream shared modules work identically
#
# Run after 03_delineate_lakes.R.
#
# Inputs:
#   sites            : tibble with site_id and lake_name columns
#   polys_all        : SpatVector from match_lake_polygons() — target lake
#                      polygons (used to exclude target lakes from upstream
#                      candidates)
#   lakes_path       : Path to OHN/OIH waterbody layer (same as used in
#                      01_match_lake_polygons.R)
#   cache_dir        : Project cache directory
#   output_dir       : Root output directory
#
# Outputs:
#   output_dir/<site_id>/catchment_clipped.gpkg  — per-site clipped catchment
#   output_dir/all_catchments_clipped.gpkg        — combined for QA in QGIS
#   cache_dir/upstream_catchments/               — cached upstream lake
#                                                   catchments (reused across
#                                                   runs)
#
# Dependencies: terra, sf, dplyr, purrr, fs, cli (via utils.R)
# ---------------------------------------------------------------------------

#' Remove upstream lake catchments from all target site catchments
#'
#' @param sites                    tibble with site_id, lake_name columns
#' @param polys_all                SpatVector from match_lake_polygons()
#' @param lakes_path               Character. Path to OHN/OIH waterbody layer
#' @param cache_dir                Character. Project cache directory
#' @param output_dir               Character. Root output directory
#' @param lake_buffer_m            Numeric. Buffer (m) for upstream lake pour
#'                                 points. Default 30 m.
#' @param min_upstream_area_m2     Numeric. Minimum area (m²) for an upstream
#'                                 lake to be included. Default 10000 (1 ha).
#' @param exclude_waterbody_types  Character vector of waterbody type values
#'                                 to exclude. Default c("Pond").
#' @param manual_exclusions        Named list of OGF_ID vectors to exclude
#'                                 for specific target lakes. Names are
#'                                 site_ids. E.g. list(CAM01 = c(12345L)).
#'
#' @return A tibble summarising clipping results with columns:
#'   site_id, status, n_erased, area_km2_before, area_km2_after
remove_upstream_lake_catchments <- function(
  sites,
  polys_all,
  lakes_path,
  cache_dir,
  output_dir,
  lake_buffer_m            = 30,
  min_upstream_area_m2     = 10000,
  exclude_waterbody_types  = c("Pond"),
  manual_exclusions        = list()
) {
  upstream_cache_dir <- fs::path(cache_dir, "upstream_catchments")
  ensure_dir(upstream_cache_dir)

  # Load D8 pointer once as rasterization template (full extent)
  d8_template <- terra::rast(fs::path(cache_dir, "d8_pntr.tif"))
  poly_crs    <- terra::crs(
    terra::vect(lakes_path, what = "geoms", extent = terra::ext(0, 0, 0, 0))
  )

  cw_inform(glue::glue(
    "Removing upstream lake catchments for {nrow(sites)} site(s)..."
  ))

  # -- Step 1: Delineate and cache upstream lake catchments --------------------
  cw_inform("Step 1: Delineating upstream lake catchments...")

  purrr::walk(seq_len(nrow(sites)), function(i) {
    sid         <- sites$site_id[i]
    lake_name   <- sites$lake_name[i]
    site_dir    <- site_output_dir(output_dir, sid)
    catchment_path <- fs::path(site_dir, "catchment.gpkg")

    if (!cache_exists(catchment_path)) {
      cw_warn(glue::glue(
        "Site '{sid}': catchment.gpkg not found — upstream step skipped"
      ))
      return(invisible(NULL))
    }

    target_catchment <- terra::vect(catchment_path) |>
      terra::project(terra::crs(d8_template))
    target_ogf_ids <- polys_all[polys_all$matched_lake == lake_name, ]$OGF_ID

    # Manually excluded OGF_IDs for this site
    manual_excl <- manual_exclusions[[sid]]

    # Load OHN lakes intersecting this catchment extent (bbox filter)
    catchment_filter <- terra::project(target_catchment, poly_crs)
    upstream_cands <- terra::vect(lakes_path, extent = terra::ext(catchment_filter))

    # Filter to intersecting, non-excluded, non-target, area >= minimum
    # Always exclude rivers (same as original upstream_clipping_lake.R)
    exclude_all_types    <- unique(c("River", exclude_waterbody_types))
    exclude_types_filter <- !upstream_cands$WATERBODY_ %in% exclude_all_types
    area_filter          <- terra::expanse(upstream_cands) >= min_upstream_area_m2
    target_filter        <- !upstream_cands$OGF_ID %in% target_ogf_ids
    manual_filter        <- if (!is.null(manual_excl)) {
      !upstream_cands$OGF_ID %in% manual_excl
    } else {
      rep(TRUE, nrow(upstream_cands))
    }

    upstream_cands <- upstream_cands[
      exclude_types_filter & area_filter & target_filter & manual_filter,
    ]
    upstream_cands <- terra::project(upstream_cands, terra::crs(d8_template))

    # Precise intersection test with catchment polygon
    intersects <- terra::relate(upstream_cands, target_catchment, relation = "intersects")
    upstream_lakes <- upstream_cands[rowSums(intersects) > 0, ]

    if (nrow(upstream_lakes) == 0) {
      cw_inform(glue::glue("Site '{sid}': no upstream lakes found"))
      return(invisible(NULL))
    }

    cw_inform(glue::glue(
      "Site '{sid}': {nrow(upstream_lakes)} upstream lake(s) — delineating..."
    ))

    purrr::walk(seq_len(nrow(upstream_lakes)), function(j) {
      delineate_upstream_lake(
        ogf_id          = upstream_lakes$OGF_ID[j],
        lake_poly       = upstream_lakes[j, ],
        d8_template     = d8_template,
        cache_dir       = cache_dir,
        upstream_dir    = upstream_cache_dir,
        lake_buffer_m   = lake_buffer_m
      )
    })
  })

  # -- Step 2: Build catchment pool (target sites + upstream lakes) -----------
  cw_inform("Step 2: Building catchment pool...")

  target_pool <- purrr::map(seq_len(nrow(sites)), function(i) {
    sid  <- sites$site_id[i]
    path <- fs::path(site_output_dir(output_dir, sid), "catchment.gpkg")
    if (!cache_exists(path)) return(NULL)
    v <- terra::vect(path) |> terra::project(terra::crs(d8_template))
    v$catchment_id <- sid
    v$area_m2      <- sum(terra::expanse(v))
    v
  }) |> purrr::compact()

  upstream_files <- fs::dir_ls(upstream_cache_dir, glob = "*_catchment.gpkg")
  upstream_pool <- purrr::map(upstream_files, function(p) {
    terra::vect(p) |> terra::project(terra::crs(d8_template))
  }) |> purrr::compact()

  # Remove upstream catchments whose catchment_id (OGF_ID) matches a CAM
  # target site polygon. Without this, a target lake found as "upstream" of a
  # neighbouring site gets cached with its OGF_ID as catchment_id, then
  # appears in `others` during Step 3 and erases the site from itself.
  cam_ogf_ids_char <- as.character(polys_all$OGF_ID)
  upstream_pool <- purrr::keep(
    upstream_pool,
    function(v) !v$catchment_id[1] %in% cam_ogf_ids_char
  )

  if (length(target_pool) == 0) {
    cw_warn("No target catchments found — aborting upstream removal")
    return(tibble::tibble())
  }

  all_catchments <- do.call(
    rbind,
    c(target_pool, upstream_pool)
  )

  cw_inform(glue::glue(
    "Pool: {length(target_pool)} target + {length(upstream_pool)} upstream = ",
    "{nrow(all_catchments)} total"
  ))

  # -- Step 3: Clip each target catchment -------------------------------------
  cw_inform("Step 3: Clipping target catchments...")

  results <- purrr::map(seq_len(nrow(sites)), function(i) {
    sid      <- sites$site_id[i]
    site_dir <- site_output_dir(output_dir, sid)
    out_path <- fs::path(site_dir, "catchment_clipped.gpkg")

    focal <- all_catchments[all_catchments$catchment_id == sid, ]
    if (nrow(focal) == 0) {
      cw_warn(glue::glue("Site '{sid}': not in catchment pool — skipping"))
      return(tibble::tibble(
        site_id = sid, status = "skipped (no catchment)",
        n_erased = NA_integer_, area_km2_before = NA_real_, area_km2_after = NA_real_
      ))
    }

    focal_area <- focal$area_m2[1]

    # Find smaller catchments that intersect the focal catchment
    others <- all_catchments[all_catchments$catchment_id != sid, ]
    if (nrow(others) == 0) {
      terra::writeVector(focal, out_path, overwrite = TRUE)
      return(tibble::tibble(
        site_id = sid, status = "no upstream",
        n_erased = 0L, area_km2_before = round(focal_area / 1e6, 4),
        area_km2_after = round(focal_area / 1e6, 4)
      ))
    }

    intersects_focal <- terra::relate(others, focal, relation = "intersects")
    smaller          <- others[
      rowSums(intersects_focal) > 0 & others$area_m2 < focal_area,
    ]

    if (nrow(smaller) == 0) {
      cw_inform(glue::glue("Site '{sid}': no smaller intersecting catchments"))
      terra::writeVector(focal, out_path, overwrite = TRUE)
      return(tibble::tibble(
        site_id = sid, status = "headwater (no clipping)",
        n_erased = 0L, area_km2_before = round(focal_area / 1e6, 4),
        area_km2_after = round(focal_area / 1e6, 4)
      ))
    }

    cw_inform(glue::glue("Site '{sid}': erasing {nrow(smaller)} smaller catchment(s)"))

    erase_mask <- terra::aggregate(smaller) |> terra::makeValid()
    clipped    <- terra::erase(focal, erase_mask)

    if (nrow(clipped) == 0 || sum(terra::expanse(clipped)) == 0) {
      cw_warn(glue::glue(
        "Site '{sid}': clipping resulted in empty geometry — saving unclipped"
      ))
      clipped <- focal
    }
    if (nrow(clipped) > 1) {
      clipped <- terra::aggregate(clipped)
    }

    area_after <- round(sum(terra::expanse(clipped)) / 1e6, 4)
    terra::writeVector(clipped, out_path, overwrite = TRUE)

    cw_inform(glue::glue(
      "Site '{sid}': {round(focal_area/1e6,4)} km² → {area_after} km²"
    ))

    tibble::tibble(
      site_id = sid, status = "clipped",
      n_erased = nrow(smaller),
      area_km2_before = round(focal_area / 1e6, 4),
      area_km2_after  = area_after
    )
  }) |>
    dplyr::bind_rows()

  # -- Step 4: Write combined output for QA -----------------------------------
  combined_path <- fs::path(output_dir, "all_catchments_clipped.gpkg")
  combined <- purrr::map(sites$site_id, function(sid) {
    p <- fs::path(site_output_dir(output_dir, sid), "catchment_clipped.gpkg")
    if (!cache_exists(p)) return(NULL)
    v <- terra::vect(p) |> sf::st_as_sf()
    v$site_id <- sid
    v
  }) |> purrr::compact()

  if (length(combined) > 0) {
    combined_sf <- dplyr::bind_rows(combined)
    sf::st_write(combined_sf, combined_path, delete_dsn = TRUE, quiet = TRUE)
    cw_inform(glue::glue(
      "Combined {nrow(combined_sf)} clipped catchment(s) → {combined_path}"
    ))
  }

  results
}

# -- Upstream lake delineation helper ----------------------------------------

#' Delineate and cache a single upstream lake catchment
#'
#' Internal function. Caches the result in upstream_dir using OGF_ID as the
#' filename. Returns the path to the cached vector file, or NULL on failure.
#'
#' @param ogf_id       Integer. OHN OGF_ID for the upstream lake
#' @param lake_poly    SpatVector. Single upstream lake polygon feature
#' @param d8_template  SpatRaster. D8 pointer used as rasterization template
#' @param cache_dir    Character. Project cache directory
#' @param upstream_dir Character. Upstream cache subdirectory
#' @param lake_buffer_m Numeric. Buffer (m) applied to lake polygon
#'
#' @return Character path to cached catchment.gpkg, or NULL on failure
delineate_upstream_lake <- function(
  ogf_id,
  lake_poly,
  d8_template,
  cache_dir,
  upstream_dir,
  lake_buffer_m
) {
  rast_path <- fs::path(upstream_dir, paste0(ogf_id, "_catchment.tif"))
  vect_path <- fs::path(upstream_dir, paste0(ogf_id, "_catchment.gpkg"))

  if (cache_exists(rast_path) && cache_exists(vect_path)) {
    return(invisible(vect_path))
  }

  cw_inform(glue::glue("  Upstream OGF_ID {ogf_id}..."))

  tryCatch({
    d8_pntr_path    <- fs::path_abs(fs::path(cache_dir, "d8_pntr.tif"))
    flow_accum_path <- fs::path(cache_dir, "flow_accum.tif")

    lake_proj     <- terra::project(lake_poly, terra::crs(d8_template))
    lake_buffered <- terra::buffer(lake_proj, width = lake_buffer_m)
    lake_rast     <- terra::rasterize(lake_buffered, d8_template, field = 1, background = NA)

    pour_path <- fs::path(upstream_dir, paste0(ogf_id, "_pourpoint.tif"))
    terra::writeRaster(lake_rast, pour_path, overwrite = TRUE)

    whitebox::wbt_watershed(
      d8_pntr  = d8_pntr_path,
      pour_pts = fs::path_abs(pour_path),
      output   = fs::path_abs(rast_path)
    )

    r <- terra::rast(rast_path) |> terra::trim()
    if (any(is.nan(as.vector(terra::ext(r))))) {
      cw_warn(glue::glue("  Upstream OGF_ID {ogf_id}: degenerate raster — skipping"))
      return(invisible(NULL))
    }

    terra::writeRaster(r, rast_path, overwrite = TRUE)

    v <- r |> terra::as.polygons() |> terra::project(terra::crs(d8_template))
    v$catchment_id <- as.character(ogf_id)
    v$area_m2      <- sum(terra::expanse(v))

    if (cache_exists(vect_path)) fs::file_delete(vect_path)
    terra::writeVector(v, vect_path)

    invisible(vect_path)
  }, error = function(e) {
    cw_warn(glue::glue("  Upstream OGF_ID {ogf_id}: failed — {e$message}"))
    invisible(NULL)
  })
}
