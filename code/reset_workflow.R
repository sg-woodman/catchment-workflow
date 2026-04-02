# reset_workflow.R
# ---------------------------------------------------------------------------
# Clears cache and/or output directories to allow a clean rerun of the
# catchment delineation workflow.
#
# Three reset levels are available:
#
#   reset_all()         — deletes everything in cache/ and output/
#                         Use when starting fresh or after major code changes
#
#   reset_cache()       — deletes only cache/ (group-level rasters)
#                         Preserves site output folders
#                         Use when group processing needs to rerun but
#                         site outputs are still valid
#
#   reset_sites()       — deletes only output/<site_id>/ folders
#                         Preserves cache/ (group rasters)
#                         Use when only delineation needs to rerun
#
#   reset_group(group)  — deletes cache/<group_id>/ and all site output
#                         folders belonging to that group
#                         Use when one group needs reprocessing
#
#   reset_site(site)    — deletes output/<site_id>/ for a single site
#                         Use when one site's delineation needs rerunning
#                         (e.g. after manually editing pour_point_snapped.shp)
#
# USAGE:
#   source("reset_workflow.R")
#   reset_all()           # full reset
#   reset_sites(sites)    # rerun delineation only
#   reset_site("L03")     # rerun one site
# ---------------------------------------------------------------------------

#' Full reset — delete all cache and output content
#'
#' @param cache_dir  Character. Root cache directory. Default "cache".
#' @param output_dir Character. Root output directory. Default "output".
#' @param confirm    Logical. If TRUE (default), asks for confirmation first.
reset_all <- function(
    cache_dir  = "cache",
    output_dir = "output",
    confirm    = TRUE
) {
  if (confirm) {
    answer <- readline(glue::glue(
      "This will delete ALL contents of '{cache_dir}/' and '{output_dir}/'. ",
      "Type 'yes' to confirm: "
    ))
    if (tolower(trimws(answer)) != "yes") {
      message("Reset cancelled.")
      return(invisible(NULL))
    }
  }

  .reset_dir(cache_dir,  "cache")
  .reset_dir(output_dir, "output")
  message("Full reset complete.")
}

#' Reset cache only — delete group-level raster cache
#'
#' @param cache_dir Character. Root cache directory. Default "cache".
#' @param confirm   Logical. If TRUE (default), asks for confirmation first.
reset_cache <- function(cache_dir = "cache", confirm = TRUE) {
  if (confirm) {
    answer <- readline(glue::glue(
      "This will delete ALL contents of '{cache_dir}/'. ",
      "Type 'yes' to confirm: "
    ))
    if (tolower(trimws(answer)) != "yes") {
      message("Reset cancelled.")
      return(invisible(NULL))
    }
  }
  .reset_dir(cache_dir, "cache")
  message("Cache reset complete.")
}

#' Reset site outputs only — delete all per-site output folders
#'
#' @param sites      Validated sites tibble. Used to find site_ids.
#' @param output_dir Character. Root output directory. Default "output".
#' @param confirm    Logical. If TRUE (default), asks for confirmation first.
reset_sites <- function(sites, output_dir = "output", confirm = TRUE) {
  if (confirm) {
    answer <- readline(glue::glue(
      "This will delete output folders for {nrow(sites)} site(s). ",
      "Type 'yes' to confirm: "
    ))
    if (tolower(trimws(answer)) != "yes") {
      message("Reset cancelled.")
      return(invisible(NULL))
    }
  }

  purrr::walk(sites$site_id, function(sid) {
    site_dir <- fs::path(output_dir, sid)
    if (fs::dir_exists(site_dir)) {
      fs::dir_delete(site_dir)
      fs::dir_create(site_dir)
      message(glue::glue("  Reset: {site_dir}"))
    }
  })

  message(glue::glue("Reset complete for {nrow(sites)} site(s)."))
}

#' Reset a single group — delete group cache and its sites' output folders
#'
#' @param group_id     Character. The group_id to reset.
#' @param sites        Validated sites tibble. Used to find sites in the group.
#' @param cache_dir    Character. Root cache directory. Default "cache".
#' @param output_dir   Character. Root output directory. Default "output".
#' @param confirm      Logical. If TRUE (default), asks for confirmation first.
reset_group <- function(
    group_id,
    sites,
    cache_dir  = "cache",
    output_dir = "output",
    confirm    = TRUE
) {
  grp_sites <- dplyr::filter(sites, group_id == !!group_id)

  if (confirm) {
    answer <- readline(glue::glue(
      "This will delete cache for group '{group_id}' and output folders for ",
      "{nrow(grp_sites)} site(s): {paste(grp_sites$site_id, collapse = ', ')}. ",
      "Type 'yes' to confirm: "
    ))
    if (tolower(trimws(answer)) != "yes") {
      message("Reset cancelled.")
      return(invisible(NULL))
    }
  }

  # Delete group cache
  grp_cache <- fs::path(cache_dir, group_id)
  if (fs::dir_exists(grp_cache)) {
    fs::dir_delete(grp_cache)
    fs::dir_create(grp_cache)
    message(glue::glue("  Reset cache: {grp_cache}"))
  }

  # Delete site output folders for this group
  purrr::walk(grp_sites$site_id, function(sid) {
    site_dir <- fs::path(output_dir, sid)
    if (fs::dir_exists(site_dir)) {
      fs::dir_delete(site_dir)
      fs::dir_create(site_dir)
      message(glue::glue("  Reset site: {site_dir}"))
    }
  })

  message(glue::glue("Reset complete for group '{group_id}'."))
}

#' Reset a single site — delete its output folder and recreate it empty
#'
#' Use after manually editing pour_point_snapped.shp to rerun delineation
#' for a specific site without affecting others.
#'
#' @param site_id    Character. The site_id to reset.
#' @param output_dir Character. Root output directory. Default "output".
reset_site <- function(site_id, output_dir = "output") {
  site_dir <- fs::path(output_dir, site_id)

  if (!fs::dir_exists(site_dir)) {
    message(glue::glue("Site '{site_id}': output folder not found, nothing to reset."))
    return(invisible(NULL))
  }

  fs::dir_delete(site_dir)
  fs::dir_create(site_dir)
  message(glue::glue("Reset complete for site '{site_id}'."))
}

# -- Internal helpers --------------------------------------------------------

#' Delete and recreate a directory
#' @param path Character. Directory path
#' @param label Character. Label for log messages
.reset_dir <- function(path, label) {
  if (fs::dir_exists(path)) {
    fs::dir_delete(path)
    message(glue::glue("  Deleted: {path}/"))
  }
  fs::dir_create(path)
  message(glue::glue("  Recreated: {path}/"))
}
