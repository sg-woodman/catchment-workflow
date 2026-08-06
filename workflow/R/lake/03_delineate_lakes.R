# 03_delineate_lakes.R
# ---------------------------------------------------------------------------
# Delineates catchments for lake sites using a buffered lake polygon as the
# pour point. Unlike the stream workflow, pour point snapping is not needed —
# the entire lake boundary drives delineation, producing catchments that
# fully contain the lake polygon.
#
# Per-site outputs use the same file names as the stream workflow so that all
# shared downstream modules (07_reclip_outputs.R, 08_catchment_metrics.R)
# work identically for lake and stream projects.
#
# Key differences from the stream workflow:
#   - Pour point: buffered lake polygon rasterized to full D8 pointer extent
#     (not a snapped single point)
#   - D8 pointer source: OIH Enhanced Flow Direction (pre-conditioned, no
#     breaching needed) — cached in cache_dir/d8_pntr.tif
#   - No hillshade, dem_breached — not applicable to the OIH dataset
#
# Outputs (per site, written to output_dir/<site_id>/):
#   catchment.gpkg     : Catchment polygon (CRS matches D8 pointer)
#   watershed.tif      : Binary watershed raster (1 = in catchment)
#   lake_pourpoint.tif : Buffered lake raster used as pour point (editable
#                        in QGIS to correct catchment boundaries)
#   dem.tif            : OIH DEM clipped and masked to catchment
#   flow_pointer.tif   : D8 flow pointer clipped to catchment
#   flow_accum.tif     : Flow accumulation clipped to catchment
#   streams.tif        : Stream network clipped to catchment
#
# Run after 02_prepare_oih_dem.R and 01_match_lake_polygons.R.
#
# Dependencies: terra, whitebox, sf, dplyr, purrr, fs, cli (via utils.R)
# ---------------------------------------------------------------------------

#' Delineate catchments for all lake sites
#'
#' Iterates over all sites in the sites tibble, matching each to its lake
#' polygon from polys_all, then running lake-raster pour-point delineation.
#' Sites where catchment.gpkg already exists in the output directory are
#' skipped.
#'
#' @param sites        tibble with at least site_id and lake_name columns.
#'                     lake_name must match the matched_lake column in
#'                     polys_all.
#' @param polys_all    SpatVector from match_lake_polygons() with a
#'                     matched_lake column
#' @param cache_dir    Character. Project cache directory containing
#'                     d8_pntr.tif, flow_accum.tif, streams.tif, dem.tif
#' @param output_dir   Character. Root output directory (per-site folders
#'                     are created as output_dir/<site_id>/)
#' @param lake_buffer_m Numeric. Buffer (m) applied to lake polygon before
#'                     rasterizing. Ensures full lake boundary is captured.
#'                     Default 30 m.
#'
#' @return A tibble with one row per site and columns:
#'   site_id, status, catchment_km2, flagged, flag_reason
delineate_lake_catchments <- function(
  sites,
  polys_all,
  cache_dir,
  output_dir,
  lake_buffer_m = 30
) {
  # Validate required cache files before starting any sites
  required_cache <- c("d8_pntr.tif", "flow_accum.tif", "streams.tif", "dem.tif")
  missing_cache <- required_cache[
    !vapply(fs::path(cache_dir, required_cache), cache_exists, logical(1))
  ]
  if (length(missing_cache) > 0) {
    cw_abort(glue::glue(
      "Missing cache files in {cache_dir}: {paste(missing_cache, collapse = ', ')}. ",
      "Run prepare_oih_products() first."
    ))
  }

  # Load D8 pointer once as the rasterization template (full extent, no trim)
  d8_template <- terra::rast(fs::path(cache_dir, "d8_pntr.tif"))

  cw_inform(glue::glue(
    "Delineating lake catchments for {nrow(sites)} site(s)..."
  ))

  results <- purrr::map(seq_len(nrow(sites)), function(i) {
    sid        <- sites$site_id[i]
    lake_name  <- sites$lake_name[i]
    site_dir   <- site_output_dir(output_dir, sid)
    ensure_dir(site_dir)

    catchment_path <- fs::path(site_dir, "catchment.gpkg")
    if (cache_exists(catchment_path)) {
      cw_inform(glue::glue("Site '{sid}': catchment.gpkg found — skipping"))
      return(tibble::tibble(
        site_id = sid, status = "skipped (cached)",
        catchment_km2 = NA_real_, flagged = FALSE, flag_reason = NA_character_
      ))
    }

    # Match lake polygon
    lake_poly <- polys_all[polys_all$matched_lake == lake_name, ]
    if (nrow(lake_poly) == 0) {
      cw_warn(glue::glue(
        "Site '{sid}': no polygon matched for lake_name = '{lake_name}' — skipping"
      ))
      return(tibble::tibble(
        site_id = sid, status = "skipped (no polygon)",
        catchment_km2 = NA_real_, flagged = TRUE,
        flag_reason = glue::glue("no polygon matched for lake '{lake_name}'")
      ))
    }
    if (nrow(lake_poly) > 1) {
      lake_poly <- lake_poly[1, ]
      cw_warn(glue::glue(
        "Site '{sid}': multiple polygons matched — using first only"
      ))
    }

    tryCatch(
      delineate_single_lake(
        site_id      = sid,
        lake_poly    = lake_poly,
        d8_template  = d8_template,
        cache_dir    = cache_dir,
        site_dir     = site_dir,
        lake_buffer_m = lake_buffer_m
      ),
      error = function(e) {
        cw_warn(glue::glue("Site '{sid}': delineation failed — {e$message}"))
        tibble::tibble(
          site_id = sid, status = paste("failed:", e$message),
          catchment_km2 = NA_real_, flagged = TRUE, flag_reason = e$message
        )
      }
    )
  }) |>
    dplyr::bind_rows()

  n_ok      <- sum(results$status %in% c("success", "skipped (cached)"))
  n_flagged <- sum(results$flagged, na.rm = TRUE)
  cw_inform(glue::glue(
    "Lake delineation complete: {n_ok}/{nrow(sites)} sites OK, ",
    "{n_flagged} flagged."
  ))

  results
}

# -- Single-site delineation ------------------------------------------------

#' Delineate catchment for a single lake site
#'
#' Internal function called by delineate_lake_catchments(). Writes all
#' standardized per-site outputs.
#'
#' @param site_id      Character. Site identifier
#' @param lake_poly    SpatVector. Single lake polygon feature
#' @param d8_template  SpatRaster. D8 pointer used as rasterization template
#' @param cache_dir    Character. Project cache directory
#' @param site_dir     Character. Per-site output directory
#' @param lake_buffer_m Numeric. Lake polygon buffer (m)
#'
#' @return Single-row tibble with delineation result
delineate_single_lake <- function(
  site_id,
  lake_poly,
  d8_template,
  cache_dir,
  site_dir,
  lake_buffer_m
) {
  cw_inform(glue::glue("Site '{site_id}': delineating..."))

  d8_pntr_path    <- fs::path_abs(fs::path(cache_dir, "d8_pntr.tif"))
  flow_accum_path <- fs::path_abs(fs::path(cache_dir, "flow_accum.tif"))
  streams_path    <- fs::path_abs(fs::path(cache_dir, "streams.tif"))
  dem_path        <- fs::path_abs(fs::path(cache_dir, "dem.tif"))

  watershed_path    <- fs::path(site_dir, "watershed.tif")
  pourpoint_path    <- fs::path(site_dir, "lake_pourpoint.tif")
  catchment_vect    <- fs::path(site_dir, "catchment.gpkg")

  # -- Step 1: Buffer lake polygon and rasterize to full D8 extent ------------
  # CRITICAL: use d8_template (full extent, no trim) as the raster template.
  # wbt_watershed() requires the pour point raster to share extent and
  # resolution with the D8 pointer — any trimming causes a mismatch.
  lake_proj     <- terra::project(lake_poly, terra::crs(d8_template))
  lake_buffered <- terra::buffer(lake_proj, width = lake_buffer_m)
  lake_rast     <- terra::rasterize(lake_buffered, d8_template, field = 1, background = NA)

  terra::writeRaster(lake_rast, pourpoint_path, overwrite = TRUE)

  # -- Step 2: Delineate watershed using lake raster as pour point ------------
  whitebox::wbt_watershed(
    d8_pntr  = d8_pntr_path,
    pour_pts = fs::path_abs(pourpoint_path),
    output   = fs::path_abs(watershed_path)
  )

  # Trim AFTER delineation — the untrimmed raster was needed by wbt_watershed()
  watershed_rast <- terra::rast(watershed_path) |> terra::trim()
  terra::writeRaster(watershed_rast, watershed_path, overwrite = TRUE)

  # -- Step 3: Convert watershed raster to polygon ----------------------------
  catchment_sf <- watershed_rast |>
    terra::as.polygons() |>
    terra::project(terra::crs(d8_template)) |>
    sf::st_as_sf() |>
    dplyr::mutate(site_id = site_id)

  if (nrow(catchment_sf) == 0 || sf::st_is_empty(sf::st_union(catchment_sf))) {
    cw_abort(glue::glue("Site '{site_id}': watershed polygon is empty"))
  }

  sf::st_write(catchment_sf, catchment_vect, delete_dsn = TRUE, quiet = TRUE)

  # -- Step 4: Clip group rasters to catchment --------------------------------
  catchment_vect_terra <- terra::vect(catchment_sf)

  rasters_to_clip <- list(
    dem          = terra::rast(dem_path),
    flow_pointer = terra::rast(d8_pntr_path),
    flow_accum   = terra::rast(flow_accum_path),
    streams      = terra::rast(streams_path)
  )

  purrr::iwalk(rasters_to_clip, function(r, name) {
    out <- r |>
      terra::crop(catchment_vect_terra, snap = "out") |>
      terra::mask(catchment_vect_terra)

    terra::writeRaster(
      out,
      fs::path(site_dir, paste0(name, ".tif")),
      overwrite = TRUE,
      gdal = c("COMPRESS=LZW", "BIGTIFF=IF_SAFER")
    )
  })

  # -- Step 5: Report catchment size ------------------------------------------
  area_km2 <- round(
    sum(as.numeric(sf::st_area(catchment_sf))) / 1e6,
    4
  )
  n_cells  <- sum(terra::values(watershed_rast) == 1, na.rm = TRUE)
  flagged  <- n_cells < 10
  flag_reason <- if (flagged) {
    glue::glue("catchment only {n_cells} cells ({area_km2} km²)")
  } else {
    NA_character_
  }

  if (flagged) {
    cw_warn(glue::glue("Site '{site_id}': FLAGGED — {flag_reason}"))
  } else {
    cw_inform(glue::glue(
      "Site '{site_id}': {n_cells} cells, {area_km2} km²"
    ))
  }

  tibble::tibble(
    site_id       = site_id,
    status        = "success",
    catchment_km2 = area_km2,
    flagged       = flagged,
    flag_reason   = flag_reason
  )
}

# -- Pour point correction -----------------------------------------------------

#' Rerun watershed delineation for a single lake site using an edited pour point
#'
#' Use after editing lake_pourpoint.tif in QGIS to correct a catchment that
#' does not fully contain the lake or is otherwise wrong. Re-runs watershed
#' delineation, polygon conversion, and raster clipping using the edited
#' lake_pourpoint.tif. Does NOT regenerate the pour point raster from the lake
#' polygon — it reads the existing (edited) lake_pourpoint.tif directly.
#'
#' Workflow:
#'   1. Edit output_dir/<site_id>/lake_pourpoint.tif in QGIS and save
#'   2. Run rerun_watershed_lake(site_id, cache_dir, output_dir)
#'   3. Re-inspect catchment.gpkg in QGIS
#'
#' @param site_id    Character. Site identifier
#' @param cache_dir  Character. Project cache directory
#' @param output_dir Character. Root output directory
#'
#' @return Single-row tibble with delineation result
rerun_watershed_lake <- function(site_id, cache_dir, output_dir) {
  site_dir      <- site_output_dir(output_dir, site_id)
  pourpoint_path <- fs::path(site_dir, "lake_pourpoint.tif")

  if (!fs::file_exists(pourpoint_path)) {
    cw_abort(glue::glue(
      "Site '{site_id}': lake_pourpoint.tif not found at {site_dir}. ",
      "Has this site been delineated yet?"
    ))
  }

  cw_inform(glue::glue(
    "Site '{site_id}': rerunning watershed delineation with edited pour point..."
  ))

  d8_template <- terra::rast(fs::path(cache_dir, "d8_pntr.tif"))

  tryCatch(
    {
      watershed_path <- fs::path(site_dir, "watershed.tif")
      catchment_vect <- fs::path(site_dir, "catchment.gpkg")

      whitebox::wbt_watershed(
        d8_pntr  = fs::path_abs(fs::path(cache_dir, "d8_pntr.tif")),
        pour_pts = fs::path_abs(pourpoint_path),
        output   = fs::path_abs(watershed_path)
      )

      watershed_rast <- terra::rast(watershed_path) |> terra::trim()
      terra::writeRaster(watershed_rast, watershed_path, overwrite = TRUE)

      catchment_sf <- watershed_rast |>
        terra::as.polygons() |>
        terra::project(terra::crs(d8_template)) |>
        sf::st_as_sf() |>
        dplyr::mutate(site_id = site_id)

      sf::st_write(catchment_sf, catchment_vect, delete_dsn = TRUE, quiet = TRUE)

      # Force-overwrite clipped rasters
      catchment_vect_terra <- terra::vect(catchment_sf)
      rasters_to_clip <- list(
        dem          = terra::rast(fs::path(cache_dir, "dem.tif")),
        flow_pointer = terra::rast(fs::path(cache_dir, "d8_pntr.tif")),
        flow_accum   = terra::rast(fs::path(cache_dir, "flow_accum.tif")),
        streams      = terra::rast(fs::path(cache_dir, "streams.tif"))
      )
      purrr::iwalk(rasters_to_clip, function(r, name) {
        out_path <- fs::path(site_dir, paste0(name, ".tif"))
        r |>
          terra::crop(catchment_vect_terra, snap = "out") |>
          terra::mask(catchment_vect_terra) |>
          terra::writeRaster(out_path, overwrite = TRUE,
            gdal = c("COMPRESS=LZW", "BIGTIFF=IF_SAFER"))
      })

      area_km2 <- round(sum(as.numeric(sf::st_area(catchment_sf))) / 1e6, 4)
      n_cells  <- sum(terra::values(watershed_rast) == 1, na.rm = TRUE)
      flagged  <- n_cells < 10

      cw_inform(glue::glue(
        "Site '{site_id}': rerun complete — {n_cells} cells, {area_km2} km²"
      ))

      tibble::tibble(
        site_id       = site_id,
        status        = "rerun success",
        catchment_km2 = area_km2,
        flagged       = flagged,
        flag_reason   = if (flagged) glue::glue("{n_cells} cells") else NA_character_
      )
    },
    error = function(e) {
      cw_warn(glue::glue("Site '{site_id}': rerun failed — {e$message}"))
      tibble::tibble(
        site_id = site_id, status = paste("rerun failed:", e$message),
        catchment_km2 = NA_real_, flagged = TRUE, flag_reason = e$message
      )
    }
  )
}
