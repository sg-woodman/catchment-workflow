# match_lake_polygons.R
# ---------------------------------------------------------------------------
# Matches a table of lake sites to lake polygons from an OHN/OIH waterbody
# dataset using a three-pass strategy:
#
#   Pass 1 — Spatial:     buffer each site point, intersect with polygons
#   Pass 2 — Fuzzy name:  Jaro-Winkler match against OFFICIAL_N within a
#                         search radius, for any site not matched spatially
#   Pass 3 — Manual ID:   direct OGF_ID lookup for confirmed hard cases
#
# The three passes are combined and deduplicated on OGF_ID. Sites still
# unmatched after all passes are reported but do not cause an error.
#
# Inputs:
#   sites            : tibble with columns site_id, lake_name, lon, lat
#                      (lake_name is matched against OFFICIAL_N in the
#                      waterbody dataset)
#   lakes_path       : Path to OHN/OIH waterbody GeoPackage or shapefile
#   buffer_m         : Buffer radius (metres) for spatial point matching.
#                      Default 10 m (shoreline snap tolerance).
#   name_dist_max    : Maximum Jaro-Winkler distance (0–1) to accept a name
#                      match. Default 0.15 (≈ minor spelling differences).
#   name_search_km   : Search radius (km) around unmatched points for name
#                      matching candidates. Default 40 km.
#   manual_id_lookup : Optional tibble with columns OGF_ID (integer) and
#                      lake_name (character) for forced manual matches.
#   excluded_ogf_ids : Optional integer vector of OGF_IDs to exclude from
#                      spatial and name matches (e.g. known wrong polygons).
#   target_crs       : CRS string for output SpatVector. Default EPSG:3161.
#
# Returns:
#   A SpatVector in target_crs containing matched lake polygons, with added
#   columns: matched_lake (site_id), match_method, jw_distance.
#
# Dependencies: terra, sf, dplyr, purrr, stringdist, cli (via utils.R)
# ---------------------------------------------------------------------------

#' Match lake sites to OHN waterbody polygons
#'
#' @param sites            tibble with site_id, lake_name, lon, lat
#' @param lakes_path       Character. Path to OHN/OIH waterbody layer
#' @param buffer_m         Numeric. Shoreline buffer for spatial match (m)
#' @param name_dist_max    Numeric. Max Jaro-Winkler distance for name match
#' @param name_search_km   Numeric. Search radius for name candidates (km)
#' @param manual_id_lookup Optional tibble with OGF_ID and lake_name columns
#' @param excluded_ogf_ids Optional integer vector of OGF_IDs to exclude
#' @param target_crs       Character. Output CRS. Default "EPSG:3161"
#'
#' @return SpatVector of matched lake polygons with matched_lake, match_method,
#'   jw_distance columns added
match_lake_polygons <- function(
  sites,
  lakes_path,
  buffer_m         = 10,
  name_dist_max    = 0.15,
  name_search_km   = 40,
  manual_id_lookup = NULL,
  excluded_ogf_ids = NULL,
  target_crs       = "EPSG:3161"
) {
  if (!requireNamespace("stringdist", quietly = TRUE)) {
    cw_abort("Package 'stringdist' is required for fuzzy name matching.")
  }

  pts_sf <- sf::st_as_sf(sites, coords = c("lon", "lat"), crs = 4326)
  pts    <- terra::vect(pts_sf)
  pts$lake_name <- sites$lake_name

  # Project to target_crs for metric buffering
  pts_proj     <- terra::project(pts, target_crs)
  pts_buffered <- terra::buffer(pts_proj, width = buffer_m)

  # Read one row to get the native CRS of the lake layer
  poly_crs <- terra::crs(
    terra::vect(lakes_path, what = "geoms", extent = terra::ext(0, 0, 0, 0))
  )

  # Reproject buffer to polygon CRS for GDAL bbox filter
  pts_for_filter <- if (!terra::same.crs(pts_buffered, poly_crs)) {
    terra::project(pts_buffered, poly_crs)
  } else {
    pts_buffered
  }

  # -- Pass 1: Spatial match --------------------------------------------------
  cw_inform("Lake matching — Pass 1: spatial intersection")

  polys_candidate <- terra::vect(lakes_path, extent = terra::ext(pts_for_filter))
  polys_candidate <- polys_candidate[polys_candidate$WATERBODY_ != "River", ]
  polys_candidate <- terra::project(polys_candidate, target_crs)

  intersect_matrix <- terra::relate(
    polys_candidate,
    pts_buffered,
    relation = "intersects"
  )

  poly_match_idx <- which(rowSums(intersect_matrix) > 0)
  polys_spatial  <- polys_candidate[poly_match_idx, ]

  polys_spatial$matched_lake  <- apply(
    intersect_matrix[poly_match_idx, , drop = FALSE],
    1,
    function(row) pts_buffered$lake_name[which(row)[1]]
  )
  polys_spatial$match_method <- "spatial"
  polys_spatial$jw_distance  <- NA_real_

  # Remove manually excluded OGF_IDs
  if (!is.null(excluded_ogf_ids) && length(excluded_ogf_ids) > 0) {
    polys_spatial <- polys_spatial[
      !polys_spatial$OGF_ID %in% excluded_ogf_ids,
    ]
  }

  spatially_matched <- unique(polys_spatial$matched_lake)
  cw_inform(glue::glue(
    "  {nrow(polys_spatial)} polygon(s) matched spatially; ",
    "{sum(!sites$lake_name %in% spatially_matched)} site(s) unmatched"
  ))

  # -- Pass 2: Fuzzy name match -----------------------------------------------
  cw_inform("Lake matching — Pass 2: fuzzy name match")

  unmatched_sites <- dplyr::filter(sites, !lake_name %in% spatially_matched)
  polys_name_dedup <- NULL

  if (nrow(unmatched_sites) > 0) {
    unmatched_proj <- terra::project(
      terra::vect(sf::st_as_sf(unmatched_sites, coords = c("lon", "lat"), crs = 4326)),
      target_crs
    )
    search_ext <- terra::ext(terra::buffer(unmatched_proj, name_search_km * 1000))
    search_poly <- terra::project(
      terra::as.polygons(search_ext, crs = target_crs),
      poly_crs
    )

    polys_name_cand <- terra::vect(lakes_path, extent = terra::ext(search_poly))
    polys_name_cand <- polys_name_cand[polys_name_cand$WATERBODY_ != "River", ]
    polys_name_cand <- terra::project(polys_name_cand, target_crs)

    name_matches <- purrr::map(seq_len(nrow(unmatched_sites)), function(i) {
      sname    <- unmatched_sites$lake_name[i]
      site_pt  <- terra::project(
        terra::vect(sf::st_as_sf(unmatched_sites[i, ], coords = c("lon", "lat"), crs = 4326)),
        target_crs
      )
      nearby_idx <- which(
        terra::relate(
          polys_name_cand,
          terra::buffer(site_pt, name_search_km * 1000),
          relation = "intersects"
        )
      )

      if (length(nearby_idx) == 0) {
        cw_warn(glue::glue("  '{sname}' — no candidates within search radius"))
        return(NULL)
      }

      nearby        <- polys_name_cand[nearby_idx, ]
      official      <- nearby$OFFICIAL_N
      jw_dist       <- stringdist::stringdist(tolower(sname), tolower(official), method = "jw")
      best_idx      <- which.min(jw_dist)
      best_dist     <- jw_dist[best_idx]

      if (best_dist <= name_dist_max) {
        cw_inform(glue::glue(
          "  '{sname}' → '{official[best_idx]}' (JW: {round(best_dist, 3)})"
        ))
        m                <- nearby[best_idx, ]
        m$matched_lake   <- sname
        m$match_method   <- "name"
        m$jw_distance    <- best_dist
        return(m)
      } else {
        cw_warn(glue::glue(
          "  '{sname}' — best match '{official[best_idx]}' too dissimilar ",
          "(JW: {round(best_dist, 3)})"
        ))
        return(NULL)
      }
    })

    polys_name_list <- name_matches[!vapply(name_matches, is.null, logical(1))]
    if (length(polys_name_list) > 0) {
      polys_name_all  <- do.call(rbind, polys_name_list)
      polys_name_dedup <- polys_name_all[
        !polys_name_all$matched_lake %in% spatially_matched,
      ]
    }
    cw_inform(glue::glue(
      "  {if (is.null(polys_name_dedup)) 0 else nrow(polys_name_dedup)} ",
      "polygon(s) matched by name"
    ))
  } else {
    cw_inform("  All sites matched spatially — name pass skipped")
  }

  # -- Pass 3: Manual OGF_ID lookup -------------------------------------------
  polys_manual <- NULL
  if (!is.null(manual_id_lookup) && nrow(manual_id_lookup) > 0) {
    cw_inform("Lake matching — Pass 3: manual OGF_ID lookup")

    ogf_sql <- paste0(
      "SELECT * FROM ohn_waterbodies_valid WHERE OGF_ID IN (",
      paste(manual_id_lookup$OGF_ID, collapse = ", "),
      ")"
    )
    polys_manual <- sf::st_read(lakes_path, query = ogf_sql, quiet = TRUE) |>
      terra::vect() |>
      terra::project(target_crs)

    polys_manual$matched_lake <- manual_id_lookup$lake_name[
      match(polys_manual$OGF_ID, manual_id_lookup$OGF_ID)
    ]
    polys_manual$match_method <- "manual_id"
    polys_manual$jw_distance  <- NA_real_

    cw_inform(glue::glue("  {nrow(polys_manual)} polygon(s) matched by OGF_ID"))
  }

  # -- Combine and deduplicate ------------------------------------------------
  polys_list <- list(polys_spatial, polys_name_dedup, polys_manual)
  polys_list <- polys_list[!vapply(polys_list, is.null, logical(1))]
  polys_all  <- do.call(rbind, polys_list)
  polys_all  <- polys_all[!duplicated(polys_all$OGF_ID), ]

  cw_inform(glue::glue("{nrow(polys_all)} total polygon(s) matched"))

  still_unmatched <- sites$lake_name[!sites$lake_name %in% unique(polys_all$matched_lake)]
  if (length(still_unmatched) > 0) {
    cw_warn(glue::glue(
      "The following sites were not matched — manual review needed:\n",
      "{paste(still_unmatched, collapse = ', ')}"
    ))
  } else {
    cw_inform(glue::glue("All {nrow(sites)} site(s) matched successfully"))
  }

  polys_all
}
