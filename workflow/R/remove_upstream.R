# 06_remove_upstream.R
# ---------------------------------------------------------------------------
# Clips each site's catchment to remove the contributing area of any nested
# upstream sites. This is the river-site equivalent of the upstream lake
# clipping step in the CAM workflow — instead of erasing upstream lake
# catchments, smaller site catchments that fall within a larger downstream
# site's catchment are erased.
#
# Steps:
#   1. Build a pool of all site catchments (catchment.gpkg from each site
#      output directory), regardless of group_id — clipping is performed
#      across the full study area
#   2. For each site, find all OTHER site catchments that are smaller and
#      intersect the focal catchment — these are nested upstream sites
#   3. Erase the union of those smaller catchments from the focal catchment
#   4. Write the clipped catchment to catchment_clipped.gpkg in the site's
#      output directory
#   5. Combine all clipped catchments into a single GeoPackage for inspection
#
# Run after 05_delineate_sites.R.
#
# Inputs:
#   sites      : Validated sites tibble from validate_sites()
#   output_dir : Root output directory (contains one subfolder per site_id,
#                each with catchment.gpkg from 05_delineate_sites.R)
#
# Outputs (per site, written to output/<site_id>/):
#   catchment_clipped.gpkg : Catchment polygon with nested upstream site
#                            catchments erased
#
# Output (combined, written to output_dir):
#   all_catchments_clipped.gpkg : All clipped catchments in a single file,
#                                  with a site_id column
#
# Dependencies: sf, dplyr, purrr, fs (via utils.R)
# ---------------------------------------------------------------------------

#' Remove upstream nested catchments from all site catchments
#'
#' For each site, identifies other site catchments that are smaller and
#' intersect the focal catchment (i.e. nested upstream sites), erases their
#' combined extent from the focal catchment, and writes the result.
#'
#' @param sites      Validated sites tibble from validate_sites()
#' @param output_dir Character. Root output directory
#'
#' @return A tibble summarising clipping results for all sites, with columns:
#'   site_id, status, n_erased, erased_site_ids, area_km2_before,
#'   area_km2_after
remove_upstream_catchments <- function(sites, output_dir) {
  sf::sf_use_s2(FALSE)

  cw_inform(glue::glue(
    "Removing upstream nested catchments for {nrow(sites)} site(s)..."
  ))

  # -- Step 1: Build pool of all site catchments -----------------------------
  catchment_pool <- build_catchment_pool(sites, output_dir)

  if (nrow(catchment_pool) == 0) {
    cw_abort(
      "No catchment.gpkg files found. Run 05_delineate_sites.R first."
    )
  }

  cw_inform(glue::glue(
    "Catchment pool: {nrow(catchment_pool)} site(s) loaded."
  ))

  # Erase/intersection math below runs in a fixed EPSG:3979 working CRS
  # (build_catchment_pool() transforms every catchment to it) — but the
  # WRITTEN catchment_clipped.gpkg must come back out in whatever CRS the
  # rest of THIS project actually uses, not silently stay in 3979. Fine
  # (a no-op) for a 3979-native project (CELESTE); a real, previously
  # unhandled mismatch for any other working CRS — confirmed directly on
  # CAM streams (EPSG:3161 catchments): catchment_clipped.gpkg came out in
  # EPSG:3979, which then made every downstream terra::crop() against a
  # 3161-native raster fail with "[crop] extents do not overlap" — the
  # geometry was valid, just numerically in the wrong CRS relative to
  # everything else.
  native_crs <- attr(catchment_pool, "native_crs")

  # -- Steps 2-4: Clip each site's catchment ----------------------------------
  results <- purrr::map(sites$site_id, function(sid) {
    clip_site_catchment(
      site_id = sid,
      catchment_pool = catchment_pool,
      output_dir = output_dir,
      native_crs = native_crs
    )
  }) |>
    dplyr::bind_rows()

  # Report sites with nested upstream catchments removed
  with_erasures <- dplyr::filter(results, n_erased > 0)
  if (nrow(with_erasures) > 0) {
    cw_inform(glue::glue(
      "{nrow(with_erasures)} site(s) had nested upstream catchments removed:\n",
      "{paste(with_erasures$site_id, '-', with_erasures$n_erased, ",
      "'nested site(s):', with_erasures$erased_site_ids, collapse = '\n')}"
    ))
  }

  failed <- dplyr::filter(results, status != "success")
  if (nrow(failed) > 0) {
    cw_warn(glue::glue(
      "{nrow(failed)} site(s) failed:\n",
      "{paste(failed$site_id, ':', failed$status, collapse = '\n')}"
    ))
  }

  # -- Step 5: Combine all clipped catchments ---------------------------------
  combine_clipped_catchments(sites, output_dir)

  results
}

# -- Catchment pool -----------------------------------------------------------

#' Build a pool of all site catchments with area
#'
#' Loads catchment.gpkg for every site that has one, tags each with its
#' site_id and area, and combines into a single sf object. Sites without a
#' catchment.gpkg are skipped with a warning.
#'
#' @param sites      Validated sites tibble from validate_sites()
#' @param output_dir Character. Root output directory
#'
#' @return sf polygon object with columns: site_id, area_m2, geometry.
#'   Returns an empty sf object (0 rows) if no catchments are found.
build_catchment_pool <- function(sites, output_dir) {
  catchments <- purrr::map(sites$site_id, function(sid) {
    catchment_path <- fs::path(
      site_output_dir(output_dir, sid),
      "catchment.gpkg"
    )

    if (!cache_exists(catchment_path)) {
      cw_warn(glue::glue(
        "Site '{sid}': catchment.gpkg not found — excluded from pool."
      ))
      return(NULL)
    }

    sf::st_read(catchment_path, quiet = TRUE) |>
      sf::st_union() |>
      sf::st_as_sf() |>
      dplyr::rename(geometry = x) |>
      dplyr::mutate(
        site_id = sid,
        area_m2 = as.numeric(sf::st_area(geometry))
      ) |>
      dplyr::select(site_id, area_m2, geometry)
  })

  catchments <- catchments[!sapply(catchments, is.null)]

  if (length(catchments) == 0) {
    return(sf::st_sf(
      site_id = character(),
      area_m2 = numeric(),
      geometry = sf::st_sfc(crs = 3979)
    ))
  }

  # Capture the catchments' own native CRS (all sites in one project share
  # one working CRS — EPSG:3979 for a CELESTE-style project, but not
  # necessarily for others, e.g. EPSG:3161 for CAM streams) BEFORE
  # transforming to a fixed EPSG:3979 working CRS for the erase/
  # intersection math below — clip_site_catchment() transforms the final
  # output back to this before writing, so catchment_clipped.gpkg always
  # matches whatever CRS the rest of the project actually uses.
  native_crs <- sf::st_crs(catchments[[1]])

  catchments <- purrr::map(catchments, sf::st_transform, crs = 3979)

  pool <- dplyr::bind_rows(catchments)
  attr(pool, "native_crs") <- native_crs
  pool
}

# -- Single-site clipping ------------------------------------------------------

#' Clip a single site's catchment by erasing smaller intersecting catchments
#'
#' Identifies all OTHER site catchments in the pool that are smaller than the
#' focal catchment and intersect it — these are nested upstream sites. Their
#' combined extent is erased from the focal catchment. If no nested sites are
#' found (e.g. a headwater site), the unclipped catchment is written as-is.
#'
#' @param site_id        Character. Focal site identifier
#' @param catchment_pool sf polygon object from build_catchment_pool()
#' @param output_dir     Character. Root output directory
#' @param native_crs     crs object (from build_catchment_pool()'s
#'   "native_crs" attribute) to transform the final clipped geometry back
#'   to before writing — catchment_pool itself is in a fixed EPSG:3979
#'   working CRS for the erase math, which may not match the project's
#'   actual CRS. Defaults to EPSG:3979 (a no-op) if not supplied, for
#'   backward compatibility with any caller that doesn't pass it.
#'
#' @return Single-row tibble with columns: site_id, status, n_erased,
#'   erased_site_ids, area_km2_before, area_km2_after
clip_site_catchment <- function(site_id, catchment_pool, output_dir, native_crs = sf::st_crs(3979)) {
  site_dir <- site_output_dir(output_dir, site_id)
  out_path <- fs::path(site_dir, "catchment_clipped.gpkg")

  focal <- dplyr::filter(catchment_pool, site_id == !!site_id)

  if (nrow(focal) == 0) {
    cw_warn(glue::glue(
      "Site '{site_id}': not found in catchment pool — skipping."
    ))
    return(tibble::tibble(
      site_id = site_id,
      status = "skipped (no catchment.gpkg)",
      n_erased = NA_integer_,
      erased_site_ids = NA_character_,
      area_km2_before = NA_real_,
      area_km2_after = NA_real_
    ))
  }

  focal_area <- focal$area_m2[1]

  others <- dplyr::filter(catchment_pool, site_id != !!site_id)

  # Find other catchments that are smaller AND intersect the focal catchment
  intersects_focal <- sf::st_intersects(others, focal, sparse = FALSE)[, 1]
  smaller <- others$area_m2 < focal_area

  nested <- others[intersects_focal & smaller, ]

  area_km2_before <- round(focal_area / 1e6, 4)

  tryCatch(
    {
      if (nrow(nested) == 0) {
        cw_inform(glue::glue(
          "Site '{site_id}': no nested upstream sites — saving unclipped."
        ))

        clipped <- focal
        area_km2_after <- area_km2_before
        n_erased <- 0L
        erased_ids <- NA_character_
      } else {
        cw_inform(glue::glue(
          "Site '{site_id}': erasing {nrow(nested)} nested upstream site(s): ",
          "{paste(nested$site_id, collapse = ', ')}"
        ))

        erase_mask <- sf::st_union(sf::st_make_valid(nested))

        clipped_geom <- sf::st_difference(
          sf::st_make_valid(focal),
          erase_mask
        )

        # Guard against empty result from erasure
        if (nrow(clipped_geom) == 0 || all(sf::st_is_empty(clipped_geom))) {
          cw_warn(glue::glue(
            "Site '{site_id}': clipping produced empty geometry — ",
            "saving unclipped instead."
          ))
          clipped <- focal
        } else {
          clipped <- clipped_geom |>
            sf::st_union() |>
            sf::st_as_sf() |>
            dplyr::rename(geometry = x) |>
            dplyr::mutate(site_id = site_id)
        }

        area_km2_after <- round(
          as.numeric(sf::st_area(clipped$geometry[1])) / 1e6,
          4
        )
        n_erased <- nrow(nested)
        erased_ids <- paste(nested$site_id, collapse = ", ")
      }

      # Attach area metadata before writing
      clipped <- clipped |>
        dplyr::mutate(
          area_km2_before = area_km2_before,
          area_km2_after = area_km2_after,
          n_erased = n_erased
        ) |>
        dplyr::select(
          site_id,
          area_km2_before,
          area_km2_after,
          n_erased,
          geometry
        )

      # Transform back to the project's actual working CRS before writing —
      # the erase math above ran in a fixed EPSG:3979 CRS regardless of
      # what CRS this project actually uses (see build_catchment_pool()).
      clipped <- sf::st_transform(clipped, native_crs)

      sf::st_write(clipped, out_path, delete_dsn = TRUE, quiet = TRUE)

      cw_inform(glue::glue(
        "Site '{site_id}': catchment_clipped.gpkg written ",
        "({area_km2_before} -> {area_km2_after} km2)."
      ))

      tibble::tibble(
        site_id = site_id,
        status = "success",
        n_erased = n_erased,
        erased_site_ids = erased_ids,
        area_km2_before = area_km2_before,
        area_km2_after = area_km2_after
      )
    },
    error = function(e) {
      cw_warn(glue::glue("Site '{site_id}': clipping failed — {e$message}"))

      tibble::tibble(
        site_id = site_id,
        status = paste("failed:", e$message),
        n_erased = NA_integer_,
        erased_site_ids = NA_character_,
        area_km2_before = area_km2_before,
        area_km2_after = NA_real_
      )
    }
  )
}

# -- Combine clipped catchments -------------------------------------------------

#' Combine all clipped catchments into a single GeoPackage
#'
#' Reads catchment_clipped.gpkg from each site's output directory and writes
#' a combined file to the root output directory. Sites missing
#' catchment_clipped.gpkg are skipped with a warning.
#'
#' @param sites      Validated sites tibble from validate_sites()
#' @param output_dir Character. Root output directory
#'
#' @return Invisibly NULL. Called for side effects.
combine_clipped_catchments <- function(sites, output_dir) {
  out_path <- fs::path(output_dir, "all_catchments_clipped.gpkg")

  combined <- purrr::map(sites$site_id, function(sid) {
    p <- fs::path(site_output_dir(output_dir, sid), "catchment_clipped.gpkg")

    if (!cache_exists(p)) {
      cw_warn(glue::glue(
        "Site '{sid}': catchment_clipped.gpkg not found — excluded from combined output."
      ))
      return(NULL)
    }

    sf::st_read(p, quiet = TRUE)
  })

  combined <- combined[!sapply(combined, is.null)]

  if (length(combined) == 0) {
    cw_warn("No clipped catchments found — combined output not written.")
    return(invisible(NULL))
  }

  combined_sf <- dplyr::bind_rows(combined)

  sf::st_write(combined_sf, out_path, delete_dsn = TRUE, quiet = TRUE)

  cw_inform(glue::glue(
    "Combined {nrow(combined_sf)} clipped catchment(s) -> {out_path}"
  ))

  invisible(NULL)
}

