# 06_hydroweight_attributes.R
# ---------------------------------------------------------------------------
# Computes inverse distance-weighted catchment attributes for stream sites
# using the hydroweight package, with the site's pour point as the O-scheme
# target (in place of the lake polygon used by workflow/R/lake/
# 05_hydroweight_attributes.R — see that file for the shared conceptual
# design; this module differs in three structural ways):
#
#   1. target_O is a POINT (output/<site_id>/pour_point.gpkg) rather than a
#      polygon, so there is no remove_region — a point has no area to
#      exclude from the LOI summary.
#   2. dem/flow_accum are resolved PER SITE from that site's group cache
#      (cache/<PROJECT_ID>/<group_id>/), not a single project-wide raster —
#      CELESTE groups span separate HydroBasins regions, each with its own
#      DEM, unlike CAM's single contiguous lake project.
#   3. Attributes are computed for BOTH catchment versions — the unclipped
#      catchment (catchment.gpkg / streams.tif from 05_delineate_sites.R)
#      and the upstream-removed catchment (catchment_clipped.gpkg /
#      streams_clipped.tif from 06_remove_upstream.R + 07_reclip_outputs.R)
#      — tagged by a "version" column, matching 08_catchment_metrics.R's
#      convention.
#
# Supports the same categorical/continuous loi_layers descriptors as the
# lake module, PLUS a `path_template` field (containing a literal "{site_id}"
# placeholder) for LOI layers that are already provided as one pre-clipped
# raster per site rather than a single project-wide raster — e.g.
# "data/celeste_loi/ndvi_trend/{site_id}.tif". A per-site raster only needs
# to cover the UNCLIPPED catchment extent; the clipped-version pass crops/
# masks it down further, so one file per site serves both versions. Exactly
# one of `path` / `path_template` must be set per layer.
#
# Run after Stage 7 (reclip_outputs) so *_clipped rasters exist. Requires:
#   - Unclipped catchment: output_dir/<site_id>/catchment.gpkg
#   - Clipped catchment:   output_dir/<site_id>/catchment_clipped.gpkg
#   - Pour point:          output_dir/<site_id>/pour_point.gpkg
#   - Group rasters:       cache_dir/<group_id>/dem_breached.tif,
#                           cache_dir/<group_id>/flow_accum.tif
#
# Dependencies: terra, sf, hydroweight, dplyr, purrr, fs, cli (via utils.R)
# ---------------------------------------------------------------------------

#' Compute distance-weighted catchment attributes for all stream sites
#'
#' @param sites             tibble with site_id and group_id columns
#'   (validate_sites_tibble() output)
#' @param group_manifest    sf tibble from build_group_manifest() — used to
#'   resolve each site's group-level dem/flow_accum
#' @param output_dir        Character. Root output directory
#' @param cache_dir         Character. Project-level cache directory (root;
#'   per-group rasters live in cache_dir/<group_id>/)
#' @param loi_layers        Named list of LOI layer descriptors. Each element
#'   is a list with:
#'   \describe{
#'     \item{path}{Character. Path to a single project-wide raster (single
#'       or multi-layer .tif). Mutually exclusive with path_template.}
#'     \item{path_template}{Character containing "{site_id}". Path to a
#'       per-site pre-clipped raster. Mutually exclusive with path.}
#'     \item{name}{Character. Short label used as column prefix in output}
#'     \item{type}{"categorical" or "continuous"}
#'     \item{class_levels}{data.frame with ID/Class columns. Categorical
#'       only. Optional.}
#'     \item{stats}{Character vector of stat columns to keep. Continuous
#'       only. Optional.}
#'   }
#' @param catchment_versions Character vector, subset of c("unclipped",
#'   "clipped"). Default both.
#' @param weighting_schemes Character vector. Default all seven schemes;
#'   stream-based schemes (iEucS/iFLS/HAiFLS) are dropped per-site-version
#'   when no stream raster is available for that version.
#' @param inv_function      Function mapping distance (m) to weight. Default
#'   (d_km + 1)^-1.
#' @param hw_dir            Character. Hydroweight intermediate raster root.
#'   Default output_dir/hydroweight.
#' @param raster_crs        Character. Target CRS for LOI rasters. Default
#'   "EPSG:3979" (stream project native CRS).
#'
#' @return A tibble with one row per site x catchment_version successfully
#'   processed, a "version" column, and one column per LOI attribute x
#'   weighting scheme combination.
calculate_hydroweight_attributes_stream <- function(
  sites,
  group_manifest,
  output_dir,
  cache_dir,
  loi_layers,
  catchment_versions = c("unclipped", "clipped"),
  weighting_schemes = c("lumped", "iEucO", "iFLO", "HAiFLO", "iEucS", "iFLS", "HAiFLS"),
  inv_function      = function(x) (x * 0.001 + 1)^-1,
  hw_dir            = fs::path(output_dir, "hydroweight"),
  raster_crs        = "EPSG:3979"
) {
  if (!requireNamespace("hydroweight", quietly = TRUE)) {
    cw_abort(paste(
      "Package 'hydroweight' is required.",
      "Install from GitHub: remotes::install_github('bencuddihee/hydroweight')"
    ))
  }
  if (!is.list(loi_layers) || length(loi_layers) == 0) {
    cw_abort("loi_layers must be a non-empty list of LOI layer descriptors.")
  }
  catchment_versions <- match.arg(
    catchment_versions,
    c("unclipped", "clipped"),
    several.ok = TRUE
  )

  validate_loi_layers_stream(loi_layers)
  loi_prepared <- prepare_loi_layers_stream(loi_layers, cache_dir, raster_crs)

  fs::dir_create(hw_dir, recurse = TRUE)

  cw_inform(glue::glue(
    "hydroweight (stream): {nrow(sites)} site(s), {length(loi_layers)} LOI ",
    "layer(s), {length(catchment_versions)} catchment version(s)"
  ))

  jobs <- tidyr::expand_grid(site_id = sites$site_id, version = catchment_versions)

  results <- purrr::pmap(jobs, function(site_id, version) {
    grp <- sites$group_id[sites$site_id == site_id][1]
    grp_manifest <- dplyr::filter(group_manifest, group_id == grp)
    if (nrow(grp_manifest) == 0) {
      cw_warn(glue::glue(
        "Site '{site_id}' [{version}]: group '{grp}' not in group_manifest — skipping"
      ))
      return(NULL)
    }
    grp_cache <- grp_manifest$cache_dir[1]

    process_hw_site_stream(
      site_id           = site_id,
      version           = version,
      output_dir        = output_dir,
      cache_dir         = cache_dir,
      hw_dir            = hw_dir,
      dem_path          = fs::path(grp_cache, "dem_breached.tif"),
      flow_accum_path   = fs::path(grp_cache, "flow_accum.tif"),
      loi_prepared      = loi_prepared,
      loi_layers        = loi_layers,
      weighting_schemes = weighting_schemes,
      inv_function      = inv_function,
      raster_crs        = raster_crs
    )
  })
  names(results) <- paste(jobs$site_id, jobs$version, sep = "__")

  n_failed <- sum(vapply(results, is.null, logical(1)))
  if (n_failed > 0) {
    failed <- names(results)[vapply(results, is.null, logical(1))]
    cw_warn(glue::glue(
      "{n_failed} site x version combination(s) failed: {paste(failed, collapse = ', ')}"
    ))
  }

  results_clean <- results[!vapply(results, is.null, logical(1))]
  if (length(results_clean) == 0) {
    cw_warn("No sites produced hydroweight results.")
    return(tibble::tibble())
  }

  cw_inform(glue::glue(
    "hydroweight (stream) complete: {length(results_clean)}/{nrow(jobs)} ",
    "site x version combination(s) processed."
  ))

  dplyr::bind_rows(results_clean)
}

# -- LOI validation and preparation -------------------------------------------

#' Validate loi_layers list structure (stream variant: path / path_template /
#' path_lazy — exactly one)
validate_loi_layers_stream <- function(loi_layers) {
  for (i in seq_along(loi_layers)) {
    lyr <- loi_layers[[i]]
    lbl <- if (!is.null(lyr$name)) lyr$name else paste0("layer ", i)

    # NOTE: use [[ (exact) not $ (partial-matches "path" against
    # "path_template"/"path_lazy" on a list — a real bug caught during
    # testing)
    has_path <- !is.null(lyr[["path"]])
    has_template <- !is.null(lyr[["path_template"]])
    has_lazy <- !is.null(lyr[["path_lazy"]])
    if (sum(has_path, has_template, has_lazy) != 1) {
      cw_abort(glue::glue(
        "loi_layers[[{i}]] ({lbl}): exactly one of 'path', 'path_template', ",
        "or 'path_lazy' must be set."
      ))
    }
    if (has_template && !grepl("\\{site_id\\}", lyr[["path_template"]])) {
      cw_abort(glue::glue(
        "loi_layers[[{i}]] ({lbl}): path_template must contain '{{site_id}}'."
      ))
    }
    if (has_path && !cache_exists(lyr[["path"]])) {
      cw_abort(glue::glue("loi_layers[[{i}]] ({lbl}): file not found — {lyr[['path']]}"))
    }
    if (has_lazy && !cache_exists(lyr[["path_lazy"]])) {
      cw_abort(glue::glue("loi_layers[[{i}]] ({lbl}): file not found — {lyr[['path_lazy']]}"))
    }
    if (is.null(lyr$type) || !lyr$type %in% c("categorical", "continuous")) {
      cw_abort(glue::glue(
        "loi_layers[[{i}]] ({lbl}): 'type' must be 'categorical' or 'continuous'."
      ))
    }
    if (is.null(lyr$name) || !nzchar(lyr$name)) {
      cw_abort(glue::glue("loi_layers[[{i}]]: 'name' is required (used as column prefix)."))
    }
    if (lyr$type == "categorical" && !is.null(lyr$class_levels)) {
      if (!all(c("ID", "Class") %in% names(lyr$class_levels))) {
        cw_abort(glue::glue(
          "loi_layers[[{i}]] ({lbl}): class_levels must have 'ID' and 'Class' columns."
        ))
      }
    }
  }
}

#' Prepare project-wide LOI rasters (path-based layers only)
#'
#' Mirrors workflow/R/lake/05_hydroweight_attributes.R's prepare_loi_layers().
#' Layers declared with path_template or path_lazy (both resolved per site)
#' are left NULL here — see resolve_site_loi_raster().
#'
#' @param loi_layers Validated list of LOI layer descriptors
#' @param cache_dir  Character. Project cache directory
#' @param raster_crs Character. Target CRS
#' @return Named list, one element per LOI layer: a prepared SpatRaster for
#'   path-based layers, or NULL for path_template/path_lazy-based layers
prepare_loi_layers_stream <- function(loi_layers, cache_dir, raster_crs) {
  hw_cache_dir <- fs::path(cache_dir, "hydroweight_loi")
  fs::dir_create(hw_cache_dir, recurse = TRUE)

  purrr::map(loi_layers, function(lyr) {
    if (!is.null(lyr[["path_template"]]) || !is.null(lyr[["path_lazy"]])) {
      return(NULL)
    }
    prepare_one_loi_raster(
      raw_path   = lyr[["path"]],
      cache_path = fs::path(hw_cache_dir, paste0(lyr$name, "_reprojected.tif")),
      type       = lyr$type,
      class_levels = lyr$class_levels,
      raster_crs = raster_crs,
      label      = lyr$name
    )
  }) |> setNames(purrr::map_chr(loi_layers, "name"))
}

#' Resolve (crop, reproject/cache, apply levels) one per-site LOI raster
#'
#' Handles both per-site source shapes:
#'   path_template — one real file per site (already small; no crop needed
#'     before caching).
#'   path_lazy     — one shared source too large to materialize whole (e.g.
#'     a VRT mosaicing scattered regional tiles) — cropped to the site's
#'     catchment (in the source's own CRS, buffered) BEFORE reprojecting/
#'     caching, so only a small per-catchment piece is ever written to disk
#'     or passed through terra::project(). See build_mosaic_vrt() in
#'     workflow/raster_attributes.R for why the source must have its nodata
#'     set correctly first — otherwise areas outside every source tile
#'     silently read as 0, not NA.
#'
#' @param loi_desc     LOI layer descriptor with path_template or path_lazy
#' @param site_id      Character. Site identifier
#' @param catchment_sf sf polygon. The site's catchment (current version),
#'   used to crop path_lazy sources. Ignored for path_template.
#' @param cache_dir    Character. Project cache directory
#' @param raster_crs   Character. Target CRS
#' @return SpatRaster, or NULL if the site's raster file doesn't exist
resolve_site_loi_raster <- function(loi_desc, site_id, catchment_sf, cache_dir, raster_crs) {
  is_lazy <- !is.null(loi_desc[["path_lazy"]])
  raw_path <- if (is_lazy) {
    loi_desc[["path_lazy"]]
  } else {
    gsub("\\{site_id\\}", site_id, loi_desc[["path_template"]])
  }

  if (!cache_exists(raw_path)) {
    cw_warn(glue::glue(
      "Site '{site_id}', LOI '{loi_desc$name}': ",
      "{if (is_lazy) 'path_lazy source' else 'per-site raster'} not found — {raw_path}"
    ))
    return(NULL)
  }

  hw_cache_dir <- fs::path(cache_dir, "hydroweight_loi", loi_desc$name)
  fs::dir_create(hw_cache_dir, recurse = TRUE)
  cache_path <- fs::path(hw_cache_dir, paste0(site_id, "_reprojected.tif"))

  prepare_one_loi_raster(
    raw_path     = raw_path,
    cache_path   = cache_path,
    type         = loi_desc$type,
    class_levels = loi_desc$class_levels,
    raster_crs   = raster_crs,
    label        = glue::glue("{loi_desc$name} [{site_id}]"),
    crop_to      = if (is_lazy) catchment_sf else NULL
  )
}

#' Shared reprojection/caching/factor-level logic for one LOI raster
#'
#' @param crop_to Optional sf polygon (any CRS). If supplied, the raw raster
#'   is cropped to this geometry (transformed to the raster's own CRS,
#'   buffered by crop_buffer_m) BEFORE reprojecting — critical for large/
#'   lazy sources (e.g. a mosaic VRT spanning a huge, mostly-empty extent)
#'   so terra::project()/writeRaster() only ever touch a small window.
#' @param crop_buffer_m Numeric. Buffer applied to crop_to before cropping,
#'   in crop_to's units (metres for a projected CRS). Guards against edge
#'   pixels being lost to grid-rotation effects if raw and raster_crs
#'   differ. Default 500.
prepare_one_loi_raster <- function(
  raw_path,
  cache_path,
  type,
  class_levels,
  raster_crs,
  label,
  crop_to = NULL,
  crop_buffer_m = 500
) {
  if (!cache_exists(cache_path)) {
    cw_inform(glue::glue("Preparing LOI '{label}' — reprojecting to {raster_crs}..."))
    raw <- terra::rast(raw_path)

    if (!is.null(crop_to)) {
      crop_geom <- crop_to |>
        sf::st_transform(terra::crs(raw)) |>
        sf::st_buffer(crop_buffer_m)
      raw <- terra::crop(raw, terra::vect(crop_geom), snap = "out")
    }

    if (!terra::same.crs(raw, raster_crs)) {
      method <- if (type == "categorical") "near" else "bilinear"
      raw <- terra::project(raw, raster_crs, method = method)
    }
    dtype <- if (type == "categorical") "INT1U" else "FLT4S"
    terra::writeRaster(raw, cache_path, overwrite = TRUE, datatype = dtype)
    cw_inform(glue::glue("  -> Cached: {cache_path}"))
  }

  r <- terra::rast(cache_path)

  if (type == "categorical") {
    r <- terra::as.factor(r)
    if (!is.null(class_levels)) {
      levels(r) <- class_levels
    }
  }

  r
}

# -- Per-site processing -------------------------------------------------------

#' Resolve the effective streams raster path for a site/version
#'
#' @param output_dir Character. Root output directory
#' @param site_id    Character. Site identifier
#' @param version    "unclipped" or "clipped"
#' @return Character path, or NULL if not found
resolve_hw_streams_path_stream <- function(output_dir, site_id, version) {
  site_dir <- site_output_dir(output_dir, site_id)
  candidate <- if (version == "clipped") {
    fs::path(site_dir, "streams_clipped.tif")
  } else {
    fs::path(site_dir, "streams.tif")
  }
  if (cache_exists(candidate)) candidate else NULL
}

#' Process hydroweight for a single site x catchment-version combination
#'
#' @return Single-row tibble with site/version columns and all LOI
#'   attributes, or NULL
process_hw_site_stream <- function(
  site_id,
  version,
  output_dir,
  cache_dir,
  hw_dir,
  dem_path,
  flow_accum_path,
  loi_prepared,
  loi_layers,
  weighting_schemes,
  inv_function,
  raster_crs
) {
  cw_inform(glue::glue("\nhw: {site_id} [{version}]"))

  for (p in c(dem_path, flow_accum_path)) {
    if (!cache_exists(p)) {
      cw_warn(glue::glue("Site '{site_id}' [{version}]: required raster not found — {p}"))
      return(NULL)
    }
  }

  site_dir <- site_output_dir(output_dir, site_id)
  catchment_file <- if (version == "clipped") "catchment_clipped.gpkg" else "catchment.gpkg"
  catchment_path <- fs::path(site_dir, catchment_file)
  pour_point_path <- fs::path(site_dir, "pour_point.gpkg")

  if (!cache_exists(catchment_path)) {
    cw_warn(glue::glue("Site '{site_id}' [{version}]: {catchment_file} not found — skipping"))
    return(NULL)
  }
  if (!cache_exists(pour_point_path)) {
    cw_warn(glue::glue("Site '{site_id}' [{version}]: pour_point.gpkg not found — skipping"))
    return(NULL)
  }

  site_catch_sf <- sf::st_read(catchment_path, quiet = TRUE)
  pour_point_sf <- sf::st_read(pour_point_path, quiet = TRUE)

  hw_site_dir <- fs::path(hw_dir, site_id, version)
  fs::dir_create(hw_site_dir, recurse = TRUE)

  # Load and extend stream raster to full DEM extent (same fix as the lake
  # module — hydroweight() needs target_S to span the same domain as
  # dem/flow_accum to trace flow paths to streams from all catchment cells)
  streams_path <- resolve_hw_streams_path_stream(output_dir, site_id, version)
  target_S <- NULL
  schemes_final <- weighting_schemes
  stream_schemes <- c("iEucS", "iFLS", "HAiFLS")

  if (!is.null(streams_path)) {
    target_S <- tryCatch(
      {
        s <- terra::rast(streams_path)
        dem_template <- terra::rast(dem_path)
        terra::extend(s, dem_template)
      },
      error = function(e) {
        cw_warn(glue::glue(
          "Site '{site_id}' [{version}]: failed to load stream raster — ",
          "{e$message}. Dropping stream-based schemes."
        ))
        NULL
      }
    )
    if (is.null(target_S)) {
      schemes_final <- setdiff(weighting_schemes, stream_schemes)
    }
  } else {
    schemes_final <- setdiff(weighting_schemes, stream_schemes)
  }

  if (length(schemes_final) == 0) {
    cw_warn(glue::glue("Site '{site_id}' [{version}]: no weighting schemes available — skipping"))
    return(NULL)
  }

  cw_inform(glue::glue(
    "  Running hydroweight() [{paste(schemes_final, collapse=', ')}]..."
  ))

  hw <- tryCatch(
    hydroweight::hydroweight(
      hydroweight_dir  = hw_site_dir,
      target_O         = pour_point_sf,
      target_S         = target_S,
      target_uid       = paste(site_id, version, sep = "_"),
      clip_region      = site_catch_sf,
      OS_combine       = !is.null(target_S),
      dem              = dem_path,
      flow_accum       = flow_accum_path,
      weighting_scheme = schemes_final,
      inv_function     = inv_function
    ),
    error = function(e) {
      cw_warn(glue::glue("Site '{site_id}' [{version}]: hydroweight() failed — {e$message}"))
      NULL
    }
  )

  if (is.null(hw)) return(NULL)

  hw <- lapply(hw, terra::rast)
  cw_inform(glue::glue("  Distance weights: {paste(names(hw), collapse=', ')}"))

  attr_tables <- purrr::map(seq_along(loi_layers), function(i) {
    lyr_desc <- loi_layers[[i]]
    lyr_rast <- loi_prepared[[i]] %||%
      resolve_site_loi_raster(lyr_desc, site_id, site_catch_sf, cache_dir, raster_crs)

    if (is.null(lyr_rast)) return(NULL)

    cw_inform(glue::glue("  LOI '{lyr_desc$name}'..."))

    run_loi_attributes_stream(
      loi_rast   = lyr_rast,
      loi_desc   = lyr_desc,
      site_catch = site_catch_sf,
      hw         = hw,
      site_id    = site_id,
      version    = version
    )
  })

  attr_tables <- attr_tables[!vapply(attr_tables, is.null, logical(1))]
  if (length(attr_tables) == 0) return(NULL)

  result <- dplyr::bind_cols(attr_tables)
  dplyr::mutate(
    result,
    site = site_id,
    version = version,
    .before = dplyr::everything()
  )
}

# -- LOI attribute computation ------------------------------------------------

#' Compute hydroweight attributes for one LOI layer at one site/version
#'
#' Same core logic as the lake module's run_loi_attributes(), minus
#' remove_region (no polygon to exclude for a point target).
run_loi_attributes_stream <- function(
  loi_rast,
  loi_desc,
  site_catch,
  hw,
  site_id,
  version
) {
  loi_name    <- loi_desc$name
  loi_numeric <- (loi_desc$type == "continuous")

  catch_vect <- terra::vect(site_catch)
  loi_site <- tryCatch(
    {
      loi_rast |>
        terra::crop(catch_vect) |>
        terra::mask(catch_vect)
    },
    error = function(e) {
      cw_warn(glue::glue(
        "Site '{site_id}' [{version}], LOI '{loi_name}': crop/mask failed — {e$message}"
      ))
      NULL
    }
  )

  if (is.null(loi_site)) return(NULL)
  if (all(is.na(terra::values(loi_site)))) {
    cw_warn(glue::glue(
      "Site '{site_id}' [{version}], LOI '{loi_name}': all NA after crop/mask — skipping"
    ))
    return(NULL)
  }

  if (loi_numeric && terra::nlyr(loi_site) > 1) {
    return(run_loi_attributes_stream_multilayer(
      loi_site   = loi_site,
      loi_desc   = loi_desc,
      site_catch = site_catch,
      hw         = hw,
      site_id    = site_id,
      version    = version
    ))
  }

  hwa <- tryCatch(
    hydroweight::hydroweight_attributes(
      loi              = loi_site,
      loi_numeric      = loi_numeric,
      roi              = site_catch,
      roi_uid          = "1",
      roi_uid_col      = site_id,
      distance_weights = hw,
      remove_region    = NULL,
      return_products  = FALSE
    ),
    error = function(e) {
      cw_warn(glue::glue(
        "Site '{site_id}' [{version}], LOI '{loi_name}': ",
        "hydroweight_attributes() failed — {e$message}"
      ))
      NULL
    }
  )

  if (is.null(hwa)) return(NULL)

  attr_tbl <- hwa$attribute_table |>
    dplyr::select(-dplyr::any_of(c(site_id, "1")))

  if (!loi_numeric) {
    attr_tbl <- clean_categorical_columns_stream(attr_tbl, loi_name, loi_desc$class_levels)
  } else {
    attr_tbl <- clean_continuous_columns_stream(attr_tbl, loi_name, loi_desc$stats)
  }

  attr_tbl
}

#' Process a multi-layer continuous LOI one layer at a time (see lake module
#' for why: hydroweight_attributes()'s internal pivot_wider() cannot handle
#' multi-layer rasters directly).
run_loi_attributes_stream_multilayer <- function(
  loi_site,
  loi_desc,
  site_catch,
  hw,
  site_id,
  version
) {
  loi_name    <- loi_desc$name
  layer_names <- names(loi_site)

  if (anyDuplicated(layer_names) > 0 || any(!nzchar(layer_names))) {
    layer_names <- paste0(loi_name, "_", seq_len(terra::nlyr(loi_site)))
  }

  layer_tbls <- purrr::map(seq_len(terra::nlyr(loi_site)), function(k) {
    single <- loi_site[[k]]
    names(single) <- layer_names[k]

    hwa_k <- tryCatch(
      hydroweight::hydroweight_attributes(
        loi              = single,
        loi_numeric      = TRUE,
        roi              = site_catch,
        roi_uid          = "1",
        roi_uid_col      = site_id,
        distance_weights = hw,
        remove_region    = NULL,
        return_products  = FALSE
      ),
      error = function(e) {
        cw_warn(glue::glue(
          "Site '{site_id}' [{version}], LOI '{loi_name}' layer {k}: ",
          "hydroweight_attributes() failed — {e$message}"
        ))
        NULL
      }
    )

    if (is.null(hwa_k)) return(NULL)
    hwa_k$attribute_table |>
      dplyr::select(-dplyr::any_of(c(site_id, "1")))
  })

  layer_tbls <- layer_tbls[!vapply(layer_tbls, is.null, logical(1))]
  if (length(layer_tbls) == 0) return(NULL)

  combined <- dplyr::bind_cols(layer_tbls)
  clean_continuous_columns_stream(combined, loi_name, loi_desc$stats)
}

# -- Column cleaning ----------------------------------------------------------

#' Rename categorical LOI columns from Class_{id}_{scheme} to {class}_{scheme}
clean_categorical_columns_stream <- function(tbl, loi_name, class_levels) {
  if (!is.null(class_levels)) {
    levels_map <- setNames(
      as.character(class_levels$Class),
      as.character(class_levels$ID)
    )
    tbl <- dplyr::rename_with(tbl, function(nms) {
      id   <- regmatches(nms, regexpr("(?<=^Class_)\\d+", nms, perl = TRUE))
      rest <- sub("^Class_\\d+_", "", nms)
      cls  <- levels_map[id]
      changed <- !is.na(cls) & nzchar(cls)
      nms[changed] <- paste0(loi_name, "_", cls[changed], "_", rest[changed])
      nms
    })
  } else {
    tbl <- dplyr::rename_with(tbl, ~ paste0(loi_name, "_", .x))
  }

  dplyr::rename_with(tbl, ~ gsub("[^A-Za-z0-9_]", "_", .x) |>
    gsub("_+", "_", x = _) |>
    tolower() |>
    sub("_$", "", x = _))
}

#' Prefix continuous LOI columns with loi_name and optionally filter stats
clean_continuous_columns_stream <- function(tbl, loi_name, stats) {
  tbl <- dplyr::rename_with(tbl, ~ paste0(loi_name, "_", .x))

  if (!is.null(stats)) {
    keep <- grep(paste(stats, collapse = "|"), names(tbl), value = TRUE)
    tbl <- dplyr::select(tbl, dplyr::any_of(keep))
  }

  tbl
}

# -- Null coalescing operator (mirrors standard R 4.4+ behaviour) ----
`%||%` <- function(x, y) if (!is.null(x)) x else y
