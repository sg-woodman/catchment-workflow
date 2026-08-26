# prepare_ndvi.R
# =============================================================================
# Prepares the CAM stream-site NDVI (continuous) LOI raster for the stream
# hydroweight stage (workflow/CAM/run_cam_streams.R, Stage 7). Same
# reasoning as workflow/CELESTE/prepare_ndvi.R for keeping this as its own
# script: source data / how it's organized is a project-specific data-prep
# concern, not part of the standardized engine stage sequence or the
# generic hydroweight module.
#
# SOURCE DATA AND HOW IT DIFFERS FROM CELESTE'S:
#   data/ndvi/CAM/*.tif — one GeoTIFF per catchment GROUP (NOT per
#   HydroBasins region like CELESTE's tiles), already masked to that
#   group's own catchment polygon(s), exported from Google Earth Engine
#   after uploading the adjacency-grouped catchments via
#   workflow/CAM/upload_catchments_to_gee.R (Example 2b —
#   group_by_adjacency = TRUE: sites whose catchments touch/share an edge
#   were merged into one polygon before export, everything else exported
#   standalone). 42 bands, "{index}_NDVI_{year}" naming, 1984-2025,
#   EPSG:3161, 30 m — confirmed directly (terra::rast() on 4 representative
#   files, sizes from 64 to ~2.2M cells).
#
#   Consequence: a CAM NDVI file can cover MULTIPLE sites (every member of
#   an adjacency group), unlike CELESTE's per-HydroBasins-region tiles
#   which needed match_group_tiles() + mosaicking to assemble one group's
#   coverage from several small tiles. Here it's the opposite direction —
#   one file already covers several SITES, and per-site resolution means
#   figuring out which file a given site_id's data lives in, then exposing
#   that shared file under a per-site name so the existing path_template
#   contract (06_hydroweight_attributes.R's resolve_site_loi_raster())
#   works unmodified. build_cam_ndvi_site_map() + prepare_cam_ndvi_site_
#   rasters() below do exactly that (materializing a per-site SYMLINK, not
#   a full copy, to the shared per-group cleaned raster — cheap, and
#   avoids duplicating the two ~180-190 MB SUD12/SUD17 group files once
#   per member site).
#
# FILE NAMING CONVENTION (confirmed directly against every file in
# data/ndvi/CAM/, cross-checked against workflow/CAM/upload_catchments_to_
# gee.R's actual adjacency grouping re-run on output/CAM/stream_delineation/
# all_catchments_clipped.gpkg — see build_cam_ndvi_site_map()):
#   Landsat_NDVI_1984_2025_[Cc]am_<name>.tif, where <name> is either:
#     - a standalone site's own site_id, OR
#     - the ALPHABETICALLY-FIRST site_id among an adjacency group's
#       members (e.g. group {SUD12, VER01} -> "sud12"; group {NCMN,
#       SUD101, SUD102, SUD103, SUD200} -> "ncmn")
#   normalized: lowercased, non-alphanumeric characters stripped (matches
#   gee_utils.R's own normalization used when the groups were built).
#
#   Two confirmed exceptions — one-off manual naming choices/typos on the
#   export side, not systematic, and not re-derivable from the site_id
#   alone — handled via CAM_NDVI_FILENAME_OVERRIDES below:
#     - Group led by "Aurora_Whitepine" (+ Marina, Sunnywater) was
#       exported as "..._cam_aurora.tif", not the full normalized leader
#       name "aurorawhitepine".
#     - Standalone site "ILD02" was exported as "..._cam_idl02.tif" — a
#       letter-swap typo (I-L-D vs I-D-L).
#
# RAW TILE ENCODING (Google Earth Engine export): confirmed identical to
# CELESTE's — int16, NDVI x10000 (".multiply(10000).round().short()"), and
# a masked/no-valid-observation pixel written as literal 0 on export, not
# tagged NoData. Confirmed two DISTINCT sources of 0 in these files,
# directly (Wolf lake's file):
#   1. Cells outside the exported polygon's own footprint (the raster is a
#      rectangular crop, like every other export in this repo) — read 0 in
#      EVERY band, all 42 years (Wolf: 15/64 cells, every year).
#   2. A genuine "no valid composite this year" cell WITHIN the polygon
#      footprint — reads 0 (indistinguishable from case 1) on top of a
#      real 0->NA masked/no-observation year (Wolf's 1985 band: an
#      ADDITIONAL 49 cells beyond the 15 permanently-outside-AOI ones are
#      NA specifically in 1985 — a real "not enough cloud-free imagery
#      that year" gap, matching the example the user described).
#   clean_cam_ndvi_tile() applies the same fix CELESTE's clean_ndvi_tiles()
#   does — subst(0, NA) then /10000 — since both cases above are
#   legitimately "no real NDVI reading," not a genuine value of 0.
#
# PRIMARY entry point (for the "ndvi" LOI, if used — see run_cam_streams.R
# Stage 7's comment on why this isn't wired in by default alongside the
# new "ndvi_trend" LOI):
#   prepare_cam_ndvi_site_rasters(sites, output_dir, cache_dir)
#
# Usage (from run_cam_streams.R, after sourcing workflow/gee_utils.R for
# group_polygons_by_adjacency() and workflow/raster_attributes.R is NOT
# required by this file specifically, but IS required by prepare_ndvi_
# trend.R which sources this file):
#   source(here("workflow/gee_utils.R"))
#   source(here("workflow/CAM/prepare_ndvi.R"))
#   prepare_cam_ndvi_site_rasters(sites, output_dir = output_dir, cache_dir = cache_dir)
#
# Dependencies: terra, sf, purrr, fs, glue (via utils.R); workflow/
#   gee_utils.R's prepare_polygons_for_gee()/group_polygons_by_adjacency()
#   must be sourced first (used to reconstruct which sites were merged
#   into which exported group).
# =============================================================================

# Default location of the CAM stream-site NDVI exports.
CAM_NDVI_DIR <- here::here("data/ndvi/CAM")

# Confirmed one-off naming exceptions on the GEE export side (see header) —
# keyed by the NORMALIZED (lowercased, non-alphanumeric stripped) leader
# name a filename would be EXPECTED to carry; value is the normalized
# suffix the file was ACTUALLY saved under. Verified directly against
# every file in CAM_NDVI_DIR as of 2026-08-24 — if a future export batch
# renames these consistently, this table becomes a no-op (the direct match
# succeeds before this is ever consulted) rather than a wrong override.
CAM_NDVI_FILENAME_OVERRIDES <- c(
  aurorawhitepine = "aurora", # group {Aurora_Whitepine, Marina, Sunnywater} -> "..._cam_aurora.tif"
  ild02           = "idl02"   # standalone site "ILD02" -> "..._cam_idl02.tif" (letter-swap typo)
)

#' Normalize an id/name for filename matching: lowercase, strip everything
#' but letters and digits. Same normalization gee_utils.R's adjacency
#' grouping and this repo's site_id-from-raw-name derivation both use.
#'
#' @param x Character vector
#' @return Character vector, same length
cam_ndvi_normalize_id <- function(x) {
  tolower(gsub("[^A-Za-z0-9]", "", x))
}

#' Build the site_id -> raw NDVI export file mapping
#'
#' Re-derives the SAME adjacency grouping used when the catchments were
#' uploaded to GEE (workflow/CAM/upload_catchments_to_gee.R's Example 2b —
#' prepare_polygons_for_gee(group_by_adjacency = TRUE) on
#' output/CAM/stream_delineation/all_catchments_clipped.gpkg), to recover
#' which sites were merged into which exported group, then matches each
#' resulting leader name against the actual files present in ndvi_dir.
#'
#' Re-derivation (rather than a hardcoded mapping) is deliberate — it self-
#' updates if sites are added/removed/re-delineated later — but is only
#' correct as long as this is run against the SAME catchment geometry that
#' was actually exported. If you re-run pour-point corrections
#' (rerun_engine_sites()) AFTER a new NDVI export batch, re-verify this
#' mapping before trusting new results — a shifted catchment could change
#' which sites are adjacent.
#'
#' @param sites       Sites tibble (site_id column) — the mapping is
#'   restricted to just these sites, catchments for any other historical
#'   site_id present in the gpkg are ignored.
#' @param output_dir  Character. Root output directory — catchments_file
#'   is read from here.
#' @param ndvi_dir    Character. Directory of raw NDVI export files.
#'   Default CAM_NDVI_DIR.
#' @param catchments_file Character. Filename (within output_dir) of the
#'   clipped combined catchments used for the GEE upload. Default
#'   "all_catchments_clipped.gpkg" (matches upload_catchments_to_gee.R's
#'   Example 2b).
#'
#' @return A tibble with columns: site_id, leader (the adjacency group's
#'   alphabetically-first member — itself, for a standalone site),
#'   raw_path (character path to the matching file in ndvi_dir, or NA if
#'   no file matched — a warning is emitted for every such site, same
#'   "real, expected gap" pattern as harvest_regen's MOR gap / CELESTE's
#'   NDVI tile gaps).
build_cam_ndvi_site_map <- function(
  sites, output_dir, ndvi_dir = CAM_NDVI_DIR,
  catchments_file = "all_catchments_clipped.gpkg"
) {
  if (!exists("prepare_polygons_for_gee", mode = "function")) {
    cw_abort(paste(
      "build_cam_ndvi_site_map() requires workflow/gee_utils.R to be",
      "sourced first (defines prepare_polygons_for_gee())."
    ))
  }

  catchments_path <- fs::path(output_dir, catchments_file)
  if (!fs::file_exists(catchments_path)) {
    cw_abort(glue::glue("Catchments file not found: {catchments_path}"))
  }

  catchments <- sf::st_read(catchments_path, quiet = TRUE)
  catchments <- catchments[catchments$site_id %in% sites$site_id, ]

  groups <- prepare_polygons_for_gee(catchments, group_by_adjacency = TRUE, id_col = "site_id")

  site_leader <- purrr::map_dfr(groups, function(g) {
    if ("adjacency_members" %in% names(g)) {
      members <- trimws(strsplit(g$adjacency_members, ",")[[1]])
      leader <- sort(members)[1]
    } else {
      members <- g$site_id
      leader <- members
    }
    tibble::tibble(site_id = members, leader = leader)
  })

  # Resolve actual files, applying the confirmed naming exceptions first.
  files <- fs::dir_ls(ndvi_dir, glob = "*.tif")
  file_suffix <- sub("^Landsat_NDVI_[0-9]+_[0-9]+_[Cc]am_", "", fs::path_file(files))
  file_suffix <- cam_ndvi_normalize_id(sub("[.]tif$", "", file_suffix, ignore.case = TRUE))
  file_lookup <- setNames(as.character(files), file_suffix)

  site_leader$leader_norm <- cam_ndvi_normalize_id(site_leader$leader)
  overridden <- site_leader$leader_norm %in% names(CAM_NDVI_FILENAME_OVERRIDES)
  site_leader$leader_norm[overridden] <- CAM_NDVI_FILENAME_OVERRIDES[site_leader$leader_norm[overridden]]

  site_leader$raw_path <- unname(file_lookup[site_leader$leader_norm])

  unmatched <- dplyr::filter(site_leader, is.na(raw_path))
  if (nrow(unmatched) > 0) {
    cw_warn(glue::glue(
      "No NDVI export file found for {nrow(unmatched)} site(s): ",
      "{paste(unmatched$site_id, collapse = ', ')} ",
      "(expected filename suffix(es): {paste(unique(unmatched$leader_norm), collapse = ', ')}). ",
      "These sites will have no NDVI/NDVI-trend data — will read as NA downstream."
    ))
  }

  dplyr::select(site_leader, site_id, leader, raw_path)
}

#' Clean one raw CAM NDVI export file: replace the 0 masked/no-data fill
#' value with true NA, and rescale from Earth Engine's int16 x10000
#' encoding back to true NDVI. Same fix as CELESTE's clean_ndvi_tiles(),
#' applied to a single file rather than a directory (CAM's exports are
#' already per-group, not per-region tiles needing mosaicking).
#'
#' Cached under cache_dir/hydroweight_loi/ndvi_clean/ by basename — shared
#' across every site whose leader points at the same raw file, so a
#' 5-member group's raster is only ever cleaned once, not once per member.
#'
#' @param raw_path  Character. Path to one raw NDVI export .tif.
#' @param cache_dir Character. Project cache root.
#' @return Character path to the cleaned file.
clean_cam_ndvi_tile <- function(raw_path, cache_dir) {
  out_dir <- fs::path(cache_dir, "hydroweight_loi", "ndvi_clean")
  fs::dir_create(out_dir, recurse = TRUE)
  out_path <- fs::path(out_dir, fs::path_file(raw_path))

  if (!cache_exists(out_path)) {
    r <- terra::rast(raw_path)
    band_names <- names(r)
    r <- terra::subst(r, from = 0, to = NA) # EE masked-pixel/outside-AOI export fill, not a real NDVI of 0
    r <- r / 10000 # EE export encoding: int16, scaled by 10000
    names(r) <- band_names
    terra::writeRaster(r, out_path, overwrite = TRUE, datatype = "FLT4S")
  }

  out_path
}

#' Materialize a per-site cleaned, rescaled NDVI raster for every site —
#' the PRIMARY entry point for the continuous "ndvi" LOI, if used
#'
#' Cleans each DISTINCT raw export file once (clean_cam_ndvi_tile()), then
#' exposes it under every member site's own name via a symlink (not a
#' copy) at cache_dir/hydroweight_loi/ndvi/<site_id>.tif — satisfies the
#' path_template contract documented in workflow/R/stream/06_hydroweight_
#' attributes.R (resolve_site_loi_raster()) exactly: "one pre-clipped
#' raster per site." A member site's file is a superset of its own
#' catchment (it also covers its group-mates' catchments), which is fine —
#' run_loi_attributes_stream() crops/masks it down to that site's own
#' catchment regardless.
#'
#' Symlinked, not copied, deliberately — SUD12/SUD17's shared group files
#' are ~180-190 MB each; copying them once per member site (up to 5, for
#' the NCMN/SUD101/SUD102/SUD103/SUD200 group) would waste real disk space
#' for zero benefit (confirmed terra::rast() reads through a symlink
#' identically to a real file).
#'
#' @param sites      Sites tibble (site_id column).
#' @param output_dir Character. Root output directory (catchments read
#'   from here — see build_cam_ndvi_site_map()).
#' @param cache_dir  Character. Project cache root.
#' @param ndvi_dir   Character. Directory of raw NDVI export files.
#'   Default CAM_NDVI_DIR.
#'
#' @return The site map (see build_cam_ndvi_site_map()), invisibly.
prepare_cam_ndvi_site_rasters <- function(
  sites, output_dir, cache_dir, ndvi_dir = CAM_NDVI_DIR
) {
  site_map <- build_cam_ndvi_site_map(sites, output_dir, ndvi_dir)

  out_dir <- fs::path(cache_dir, "hydroweight_loi", "ndvi")
  fs::dir_create(out_dir, recurse = TRUE)

  matched <- dplyr::filter(site_map, !is.na(raw_path))
  distinct_raw <- unique(matched$raw_path)
  cleaned_lookup <- setNames(
    purrr::map_chr(distinct_raw, clean_cam_ndvi_tile, cache_dir = cache_dir),
    distinct_raw
  )

  purrr::pwalk(matched, function(site_id, leader, raw_path) {
    target <- fs::path(out_dir, paste0(site_id, ".tif"))
    if (cache_exists(target)) {
      return(invisible(NULL))
    }
    if (fs::link_exists(target) || fs::file_exists(target)) {
      fs::file_delete(target) # stale/broken link from a source that's since changed
    }
    fs::link_create(cleaned_lookup[[raw_path]], target, symbolic = TRUE)
  })

  cw_inform(glue::glue(
    "NDVI (continuous): {nrow(matched)}/{nrow(site_map)} site(s) linked ",
    "from {length(distinct_raw)} distinct export file(s)."
  ))

  invisible(site_map)
}
