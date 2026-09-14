# assign_hydrobasins_groups.R
# ---------------------------------------------------------------------------
# EMILY_NHN-specific: derives the group_id column the engine's
# grouping$strategy = "hydrobasins" requires, by spatially assigning each
# site to the HydroBasins level-6 polygon it physically falls in —
# group_id = "HYBAS_<HYBAS_ID>". Unlike CELESTE (group_id pre-assigned as
# semantic project labels like "COC"/"NIP" when its sites gpkg was first
# built), EMILY_NHN's sites have no pre-existing group structure suited to
# this engine, so this derives one directly from the same HydroBasins data
# build_group_aoi() (workflow/R/utils.R) will use to build each group's AOI
# — "one group per polygon actually touched by a site", not one group per
# TOUCHED-POLYGON-CLUSTER: checked directly, all 5 level-6 polygons the
# 177 sites fall into are mutually adjacent (st_touches() forms a single
# connected chain across all 5), so merging adjacent polygons into one
# group would collapse everything into a single 177-site group — the
# opposite of the point of grouping (amortizing DEM crop/breach/NHN
# burn-in cost across a manageable AOI, not merging distinct clusters that
# happen to share a border). Kept simple: no merging, one group per
# distinct polygon.
#
# Verified against the real data before use: all 177 sites match exactly
# one polygon (zero unmatched, zero site matching more than one), producing
# 5 groups sized 2 / 4 / 7 / 23 / 141 — a sensible split matching the real
# geographic/project clusters (visible in the source data's own
# location_project column): a Northshore/Ranger Lake cluster, a Carpenter
# Lake cluster, a Tabor Creek cluster, and one large basin covering the
# Goulais River / Harmony River / Christina Mine / Garden River / Wabos /
# Achigan / Batchewana / TLW / Turkey Lake / Will / Murphy sites.
# ---------------------------------------------------------------------------

#' Assign each site a group_id from the HydroBasins level-6 polygon it
#' falls in
#'
#' Resolves the na/ar HydroBasins region once for the whole site set (via
#' the existing resolve_hydrobasins_region(), workflow/R/utils.R — same
#' spatial na-vs-ar check build_group_aoi() uses, not a latitude
#' threshold), then spatially joins every site to that region's level-6
#' polygons. Aborts loudly — rather than silently dropping or
#' mis-grouping — if any site matches zero polygons (outside HydroBasins
#' coverage) or more than one (sitting exactly on a shared boundary
#' vertex; st_intersects() can return multiple matches there).
#'
#' @param sites tibble with site_id, lon, lat (WGS84 decimal degrees).
#'   Must NOT already have a group_id column.
#' @param hydrobasins_dir Path to HydroBasins root directory (north_america/
#'   arctic subfolders — see root CLAUDE.md's "External data" table)
#' @param hybas_level HydroBasins level to assign groups from. Default 6,
#'   matching build_group_aoi()'s own default so each group's later AOI
#'   union operates on the exact same polygons used here to assign it.
#' @return `sites` with a new group_id column ("HYBAS_<HYBAS_ID>")
assign_hydrobasins_group_ids <- function(sites, hydrobasins_dir, hybas_level = 6) {
  if ("group_id" %in% names(sites)) {
    cw_abort("assign_hydrobasins_group_ids(): sites already has a group_id column — remove it before calling this.")
  }

  sites_sf <- sf::st_as_sf(sites, coords = c("lon", "lat"), crs = 4326)

  region <- resolve_hydrobasins_region(sites_sf, hydrobasins_dir)
  level_str <- formatC(hybas_level, width = 2, flag = "0")
  hybas_path <- fs::path(
    hydrobasins_dir,
    if (region == "ar") "arctic" else "north_america",
    glue::glue("hybas_{region}_lev{level_str}_v1c.shp")
  )
  if (!fs::file_exists(hybas_path)) {
    cw_abort(glue::glue("assign_hydrobasins_group_ids(): HydroBasins file not found at {hybas_path}."))
  }

  hybas <- sf::st_read(hybas_path, quiet = TRUE) |>
    sf::st_transform(3979) |>
    sf::st_make_valid()

  pts_3979 <- sf::st_transform(sites_sf, 3979)
  joined <- sf::st_join(pts_3979, hybas["HYBAS_ID"], join = sf::st_intersects) |>
    sf::st_drop_geometry()

  # A site matching >1 polygon (sitting exactly on a shared boundary
  # vertex) produces >1 row for that site_id after st_join — caught here
  # via a row-count mismatch, rather than silently duplicating that site
  # downstream with two different group_ids.
  if (nrow(joined) != nrow(sites)) {
    dup_ids <- joined$site_id[duplicated(joined$site_id) | duplicated(joined$site_id, fromLast = TRUE)]
    cw_abort(glue::glue(
      "assign_hydrobasins_group_ids(): {length(unique(dup_ids))} site(s) matched more than one ",
      "level-{hybas_level} polygon (likely sitting exactly on a shared boundary vertex) — ",
      "group assignment is ambiguous: {paste(unique(dup_ids), collapse = ', ')}"
    ))
  }

  unmatched <- joined$site_id[is.na(joined$HYBAS_ID)]
  if (length(unmatched) > 0) {
    cw_abort(glue::glue(
      "assign_hydrobasins_group_ids(): {length(unmatched)} site(s) matched no level-{hybas_level} ",
      "polygon (outside HydroBasins coverage): {paste(unmatched, collapse = ', ')}"
    ))
  }

  group_lookup <- joined |>
    dplyr::transmute(site_id, group_id = paste0("HYBAS_", HYBAS_ID))

  sites |>
    dplyr::left_join(group_lookup, by = "site_id")
}
