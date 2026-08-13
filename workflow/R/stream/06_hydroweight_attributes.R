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
      group_id          = grp,
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
    # path_lazy may itself be a template (containing "{group_id}" and/or
    # "{site_id}") for a source that varies per group/site (e.g. one
    # rasterized harvest/regen layer per HydroBasins group) — only check
    # existence up front for a literal (non-templated) path; templated ones
    # are checked per site/group in resolve_site_loi_raster().
    if (has_lazy && !grepl("\\{(group_id|site_id)\\}", lyr[["path_lazy"]]) &&
      !cache_exists(lyr[["path_lazy"]])) {
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
#' @param group_id     Character. The site's group identifier — substituted
#'   into a "{group_id}" placeholder in path_lazy, e.g. for one shared
#'   source per HydroBasins group (rasterized once per group, cropped once
#'   more per site here since a group-scale source is itself often too
#'   large to cache whole — see build_mosaic_vrt()'s docstring for the same
#'   reasoning at province scale).
#' @param catchment_sf sf polygon. The site's catchment (current version),
#'   used to crop path_lazy sources. Ignored for path_template.
#' @param cache_dir    Character. Project cache directory
#' @param raster_crs   Character. Target CRS
#' @return SpatRaster, or NULL if the site's raster file doesn't exist
resolve_site_loi_raster <- function(loi_desc, site_id, group_id, catchment_sf, cache_dir, raster_crs) {
  is_lazy <- !is.null(loi_desc[["path_lazy"]])
  raw_path <- if (is_lazy) {
    gsub("\\{site_id\\}", site_id, gsub("\\{group_id\\}", group_id, loi_desc[["path_lazy"]]))
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
    # PHOTOMETRIC=MINISBLACK guards against a real, confirmed GDAL quirk:
    # an exactly-3-band Byte GeoTIFF gets auto-tagged as RGB imagery on
    # write, which silently renames the bands to "red"/"green"/"blue" on
    # every subsequent read — corrupting layer names for e.g. a group with
    # exactly 2 years + 1 combined band. Harmless for other band counts/
    # types; always applied for consistency.
    gdal_opts <- if (type == "categorical") "PHOTOMETRIC=MINISBLACK" else NULL
    terra::writeRaster(raw, cache_path, overwrite = TRUE, datatype = dtype, gdal = gdal_opts)
    cw_inform(glue::glue("  -> Cached: {cache_path}"))
  }

  r <- terra::rast(cache_path)

  if (type == "categorical" && !is.null(class_levels)) {
    # Multi-layer SpatRasters have unreliable factor-level semantics in
    # terra (confirmed by direct testing, not assumed):
    #   - `levels(r) <- class_levels` on the whole object only applies to
    #     layer 1; layers 2+ silently keep as.factor()'s own auto-generated
    #     labels (observed: "green"/"blue"-style defaults, not an error).
    #   - Setting it per layer (`levels(r[[k]]) <- ...`) works, but each
    #     call renames THAT layer to the class_levels attribute column's
    #     name (e.g. "Class"), discarding whatever it was named before.
    #   - Restoring names afterwards with `names(r) <- orig_names` on the
    #     already-leveled multi-layer object then *wipes the levels back
    #     out* — silently reverting to as.factor()'s defaults again.
    # The only combination that reliably keeps both name and levels intact
    # is doing it per layer on independent single-layer rasters (each
    # individually verified safe), then recombining. Cheap enough — this
    # runs once per LOI per cache miss, not per site.
    band_names <- names(r)
    bands <- lapply(seq_len(terra::nlyr(r)), function(k) {
      b <- terra::as.factor(r[[k]])
      levels(b) <- class_levels
      names(b) <- band_names[k]
      b
    })
    r <- terra::rast(bands)
  } else if (type == "categorical") {
    r <- terra::as.factor(r)
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
  group_id,
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

  # st_make_valid() defensively — hydroweight::hydroweight() does its own
  # internal GEOS geometry processing (union/buffer/etc. against the clip
  # region) and a self-touching artifact in a raster-to-polygon conversion
  # upstream (watershed_to_polygon() in 05_delineate_sites.R) is enough to
  # trip a "TopologyException: side location conflict" there even though
  # the polygon reads back fine on its own. Confirmed directly: 3 sites
  # (NBI1/NBI4/NBI5, unclipped version only — the largest unclipped
  # catchments in their group, so most likely to reach whatever pixel
  # produces the artifact) failed with this exact error at the identical
  # coordinate until this fix was added.
  site_catch_sf <- sf::st_read(catchment_path, quiet = TRUE) |> sf::st_make_valid()
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
      resolve_site_loi_raster(lyr_desc, site_id, group_id, site_catch_sf, cache_dir, raster_crs)

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

  if (terra::nlyr(loi_site) > 1) {
    return(run_loi_attributes_stream_multilayer(
      loi_site    = loi_site,
      loi_desc    = loi_desc,
      loi_numeric = loi_numeric,
      site_catch  = site_catch,
      hw          = hw,
      site_id     = site_id,
      version     = version
    ))
  }

  present_id <- if (!loi_numeric) single_class_value(loi_site) else NULL

  if (!is.null(present_id)) {
    attr_tbl <- degenerate_categorical_table(present_id, loi_desc$class_levels, names(hw))
  } else {
    hwa <- tryCatch(
      hydroweight::hydroweight_attributes(
        loi              = if (loi_numeric) {
          loi_site
        } else {
          prep_categorical_for_hydroweight(loi_site, loi_desc$class_levels)
        },
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
  }

  if (!loi_numeric) {
    attr_tbl <- clean_categorical_columns_stream(attr_tbl, loi_name, loi_desc$class_levels)
    attr_tbl <- ensure_full_categorical_schema(attr_tbl, loi_name, loi_desc$class_levels, names(hw))
  } else {
    attr_tbl <- clean_continuous_columns_stream(attr_tbl, loi_name, loi_desc$stats)
  }

  attr_tbl
}

#' Process a multi-layer LOI (continuous or categorical)
#'
#' Continuous and categorical LOIs are handled completely differently here,
#' and the difference is not cosmetic — verified directly against real
#' hydroweight() output AND the installed hydroweight package's own source
#' (v2.0.0, deparse(body(hydroweight::hydroweight_attributes))):
#'
#'   continuous: hydroweight_attributes() has genuine native multi-layer
#'     support — it vectorizes stat computation across every layer in ONE
#'     call (terra::global(), raster arithmetic against the distance-weight
#'     rasters), producing already-distinct per-layer column names straight
#'     from the raster's own band names (e.g. "y2015_mean"). Confirmed
#'     ~10x faster than the old per-band R loop (42.6s -> ~4.2s on a real
#'     42-band NDVI raster) with byte-identical output. See
#'     run_loi_attributes_stream_multilayer_continuous() for the one
#'     documented exception (median) and why it's handled separately.
#'   categorical: raw column names are ALWAYS "Class_{id}_{scheme}_prop" —
#'     the LOI's names() has no effect on this whatsoever (confirmed: an
#'     explicitly-renamed layer still produced "Class_..." columns, not
#'     "<name>_..."). Every year's raw output is therefore byte-identical
#'     to every other year's, and a single multi-layer call would silently
#'     collide/mangle columns — genuinely needs the per-layer loop, with
#'     clean_categorical_columns_stream() run PER LAYER (layer-specific
#'     name prefix) BEFORE binding.
run_loi_attributes_stream_multilayer <- function(
  loi_site,
  loi_desc,
  loi_numeric,
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

  if (loi_numeric) {
    combined <- run_loi_attributes_stream_multilayer_continuous(
      loi_site      = loi_site,
      layer_names   = layer_names,
      loi_name      = loi_name,
      site_catch    = site_catch,
      hw            = hw,
      site_id       = site_id,
      version       = version,
      numeric_stats = loi_desc$stats
    )
    if (is.null(combined)) return(NULL)
    return(clean_continuous_columns_stream(combined, loi_name, loi_desc$stats))
  }

  layer_tbls <- purrr::map(seq_len(terra::nlyr(loi_site)), function(k) {
    single <- loi_site[[k]]
    if (loi_numeric) {
      names(single) <- layer_names[k]
    }

    if (all(is.na(terra::values(single)))) {
      return(NULL)
    }

    present_id <- if (!loi_numeric) single_class_value(single) else NULL

    if (!is.null(present_id)) {
      tbl_k <- degenerate_categorical_table(present_id, loi_desc$class_levels, names(hw))
    } else {
      hwa_k <- tryCatch(
        hydroweight::hydroweight_attributes(
          loi              = if (loi_numeric) {
            single
          } else {
            prep_categorical_for_hydroweight(single, loi_desc$class_levels)
          },
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
            "Site '{site_id}' [{version}], LOI '{loi_name}' layer {k}: ",
            "hydroweight_attributes() failed — {e$message}"
          ))
          NULL
        }
      )

      if (is.null(hwa_k)) return(NULL)
      tbl_k <- hwa_k$attribute_table |>
        dplyr::select(-dplyr::any_of(c(site_id, "1")))
    }

    if (!loi_numeric) {
      layer_prefix <- paste0(loi_name, "_", layer_names[k])
      tbl_k <- clean_categorical_columns_stream(tbl_k, layer_prefix, loi_desc$class_levels)
      tbl_k <- ensure_full_categorical_schema(tbl_k, layer_prefix, loi_desc$class_levels, names(hw))
    }
    tbl_k
  })

  layer_tbls <- layer_tbls[!vapply(layer_tbls, is.null, logical(1))]
  if (length(layer_tbls) == 0) return(NULL)

  combined <- dplyr::bind_cols(layer_tbls)
  if (loi_numeric) {
    clean_continuous_columns_stream(combined, loi_name, loi_desc$stats)
  } else {
    combined
  }
}

#' Compute a multi-layer CONTINUOUS LOI's attributes in one native
#' hydroweight_attributes() call, working around a confirmed package bug
#' in its "median" stat
#'
#' hydroweight_attributes() (v2.0.0) vectorizes every "lumped" stat
#' (mean/sd/median/min/max/sum/cell_count/NA_cell_count) and every
#' distance-weighted stat (distwtd_mean/distwtd_sd, per weighting scheme)
#' across all layers of `loi` in one call — EXCEPT median, which has a
#' real, reproducible bug in this package version: unlike every other
#' lumped stat (computed via a shared `fun_global()` helper that leaves
#' terra::global()'s own per-layer row names — e.g. "y2015" — intact),
#' median's computation explicitly does
#' `stats::setNames(unlist(terra::global(loi, function(x)
#' stats::median(x, na.rm = TRUE))), "median")` — renaming the result to
#' the literal string "median" for EVERY layer. The function's downstream
#' column-naming step extracts a per-layer index from the trailing digits
#' of each stat's names() (to tell "y2015_mean" apart from "y2016_mean");
#' for median, every layer's name is now the same non-numeric string
#' "median", so that extraction fails (produces NA) and pivot_wider()
#' throws "'list' object cannot be coerced to type 'double'" the moment
#' more than one layer is present. Confirmed directly (real 42-band NDVI
#' raster): fails every time with median included; a call with
#' loi_numeric_stats excluding "median" alone succeeds.
#'
#' Fix: exclude "median" from the one-shot call, then compute it
#' separately using the EXACT same method hydroweight_attributes() itself
#' uses internally (terra::global(loi, function(x) stats::median(x, na.rm
#' = TRUE))) — confirmed byte-for-byte identical to the old per-band-loop
#' path's median output (diff at floating-point noise, ~1e-16) and fast
#' (~0.3s for 42 layers — one vectorized terra call, not a loop), then
#' merged back in as "{layer_name}_median" — bare, no loi_name prefix,
#' matching every other stat's raw naming here; the caller
#' (run_loi_attributes_stream_multilayer()) adds the loi_name prefix once,
#' uniformly, via clean_continuous_columns_stream().
#'
#' All-NA layers are dropped up front (e.g. a year with no coverage for
#' this site's tiles) — matches the old per-band loop's behavior of
#' silently excluding a degenerate layer rather than erroring or passing
#' an empty layer into hydroweight_attributes().
#'
#' @param numeric_stats Character vector of grep patterns (same semantics
#'   as loi_desc$stats elsewhere in this file — see
#'   clean_continuous_columns_stream()) restricting which of the 10
#'   possible stats (mean/sd/median/min/max/sum/cell_count/NA_cell_count/
#'   distwtd_mean/distwtd_sd) are actually computed — not just filtered
#'   from the output afterward, since the caller may not want the cost of
#'   computing them either (e.g. sum/cell_count are meaningless for NDVI
#'   analysis). NULL (default) computes every stat, matching prior
#'   behavior. Unlike clean_continuous_columns_stream()'s post-hoc filter,
#'   this must match against the bare stat name only (not yet loi_name-
#'   prefixed), so it's resolved once here via the same "grep any pattern
#'   against a fixed candidate list" approach.
#'
#' @return A tibble of raw (not yet cleaned) attribute columns, or NULL if
#'   every layer was all-NA or the one-shot call failed.
run_loi_attributes_stream_multilayer_continuous <- function(
  loi_site,
  layer_names,
  loi_name,
  site_catch,
  hw,
  site_id,
  version,
  numeric_stats = NULL
) {
  names(loi_site) <- layer_names

  full_stats <- c(
    "mean", "sd", "median", "min", "max", "sum", "cell_count",
    "NA_cell_count", "distwtd_mean", "distwtd_sd"
  )
  wanted_stats <- if (is.null(numeric_stats)) {
    full_stats
  } else {
    full_stats[grepl(paste(numeric_stats, collapse = "|"), full_stats)]
  }
  compute_median <- "median" %in% wanted_stats
  hwa_stats <- setdiff(wanted_stats, "median") # median always computed separately — see below

  valid_layers <- which(!vapply(seq_len(terra::nlyr(loi_site)), function(k) {
    all(is.na(terra::values(loi_site[[k]])))
  }, logical(1)))
  if (length(valid_layers) == 0) return(NULL)
  loi_valid <- loi_site[[valid_layers]]

  attr_tbl <- NULL
  if (length(hwa_stats) > 0) {
    hwa <- tryCatch(
      hydroweight::hydroweight_attributes(
        loi               = loi_valid,
        loi_numeric       = TRUE,
        loi_numeric_stats = hwa_stats,
        roi               = site_catch,
        roi_uid           = "1",
        roi_uid_col       = site_id,
        distance_weights  = hw,
        remove_region     = NULL,
        return_products   = FALSE
      ),
      error = function(e) {
        cw_warn(glue::glue(
          "Site '{site_id}' [{version}], LOI '{loi_name}': one-shot ",
          "hydroweight_attributes() failed — {e$message}"
        ))
        NULL
      }
    )
    if (is.null(hwa)) return(NULL)
    attr_tbl <- hwa$attribute_table |> dplyr::select(-dplyr::any_of(c(site_id, "1")))
  }

  if (!compute_median) {
    return(attr_tbl %||% tibble::tibble())
  }

  # No loi_name prefix here (unlike the docstring's naming description) —
  # hydroweight_attributes()'s own per-layer output for every other stat is
  # bare "{layer_name}_{stat}" (e.g. "ndvi_mosaic_1_mean"), and the OUTER
  # run_loi_attributes_stream_multilayer() wrapper adds the loi_name prefix
  # once, uniformly, via clean_continuous_columns_stream() after this
  # function returns. Prefixing it here too double-prefixes median columns
  # only (confirmed directly: "ndvi_ndvi_ndvi_mosaic_1_median" instead of
  # "ndvi_ndvi_mosaic_1_median" — caught by the end-to-end column-name diff
  # against the old per-band-loop output).
  med_raw <- terra::global(loi_valid, function(x) stats::median(x, na.rm = TRUE))
  med_named <- stats::setNames(
    as.list(med_raw[[1]]),
    paste0(names(loi_valid), "_median")
  )
  med_tbl <- tibble::as_tibble(med_named)

  if (is.null(attr_tbl)) med_tbl else dplyr::bind_cols(attr_tbl, med_tbl)
}

# -- Degenerate (single-class) categorical layers -----------------------------

#' Check whether a categorical LOI has only one distinct non-NA value
#'
#' Common for a sparse categorical time series clipped to a small
#' catchment — e.g. most years of a harvest/regen record having zero
#' recorded activity within a given site's catchment, leaving every cell
#' "other".
#'
#' @param loi SpatRaster, single layer
#' @return The single present value (numeric), or NULL if 2+ distinct
#'   values are present (the normal, non-degenerate case) or the layer is
#'   entirely NA (handled separately by the caller).
single_class_value <- function(loi) {
  u <- unique(terra::values(loi))
  u <- u[!is.na(u)]
  if (length(u) == 1) u[1] else NULL
}

#' Build a "Class_{id}_{scheme}_prop" attribute table directly for a
#' degenerate (single-class) categorical layer, bypassing
#' hydroweight_attributes() entirely
#'
#' hydroweight_attributes() produces an ambiguous, id-less column name
#' ("Class_{scheme}_prop", no class number) when only one class is present
#' in the LOI — confirmed empirically against real output, not documented —
#' which clean_categorical_columns_stream() cannot parse (there's no id to
#' look up in class_levels). Sidestepped here with a shortcut that is
#' EXACT, not approximate: when one class covers 100% of the (already
#' catchment-masked) LOI extent, its distance-weighted proportion is
#' exactly 1.0 under ANY weighting scheme — the class indicator is 1
#' everywhere in the ROI, so sum(indicator * weight) / sum(weight) =
#' sum(weight) / sum(weight) = 1 regardless of the actual weight values —
#' and every other class is exactly 0.0. No need to touch the distance-
#' weight rasters, or even call hydroweight_attributes(), at all.
#'
#' @param present_id   Numeric. The one class ID present (from
#'   single_class_value()).
#' @param class_levels data.frame with ID/Class columns.
#' @param scheme_names Character vector of weighting scheme names (e.g.
#'   names(hw)) — matches what hydroweight_attributes() would have used.
#' @return tibble with one row, columns "Class_{id}_{scheme}_prop" for
#'   every id x scheme combination — the same raw format
#'   clean_categorical_columns_stream() expects from a real
#'   hydroweight_attributes() call.
degenerate_categorical_table <- function(present_id, class_levels, scheme_names) {
  grid <- tidyr::expand_grid(id = class_levels$ID, scheme = scheme_names)
  vals <- as.list(as.numeric(grid$id == present_id))
  names(vals) <- paste0("Class_", grid$id, "_", grid$scheme, "_prop")
  tibble::as_tibble(vals)
}

# -- Column cleaning ----------------------------------------------------------

#' Force a categorical LOI's layer name to its class_levels attribute
#' column name (e.g. "Class"), right before hydroweight_attributes()
#'
#' hydroweight_attributes()'s raw categorical output column format depends
#' on the LOI raster's names() in ways that are not documented and were
#' confirmed empirically, not assumed:
#'   - names() == the class_levels attribute column name (e.g. "Class") ->
#'     "Class_{id}_{scheme}_prop" (scheme included even with 7 real
#'     schemes — verified against real hydroweight() output).
#'   - names() == anything else (e.g. "y2015") -> "{that_name}_{id}_
#'     {scheme}_prop" instead — a DIFFERENT format, breaking
#'     clean_categorical_columns_stream()'s parsing.
#' Forcing the name here — on a throwaway copy, right before the call —
#' guarantees the first (parseable) format regardless of what the LOI was
#' named upstream (e.g. a per-year band explicitly named "y2015" so it
#' stays distinguishable through crop/mask). The year/layer-specific prefix
#' in the FINAL column names comes from the `loi_name` argument callers
#' pass to clean_categorical_columns_stream() afterwards, not from this.
#'
#' @param loi          SpatRaster, single categorical layer
#' @param class_levels data.frame with ID + one other column (e.g. Class),
#'   or NULL (in which case `loi` is returned unchanged — no known attribute
#'   column name to force it to).
#' @return SpatRaster with names() set to the class_levels attribute column
prep_categorical_for_hydroweight <- function(loi, class_levels) {
  if (is.null(class_levels)) {
    return(loi)
  }
  attr_col <- setdiff(names(class_levels), "ID")[1]
  names(loi) <- attr_col
  loi
}

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

  dplyr::rename_with(tbl, snake_case_cols)
}

#' Shared column-name normalization: lowercase snake_case, no repeated/
#' trailing underscores. Used by clean_categorical_columns_stream() and
#' ensure_full_categorical_schema() — kept as one function so the two never
#' silently drift apart and fail to recognize each other's column names.
snake_case_cols <- function(x) {
  x |>
    gsub("[^A-Za-z0-9_]", "_", x = _) |>
    gsub("_+", "_", x = _) |>
    tolower() |>
    sub("_$", "", x = _)
}

#' Ensure every {class} x {scheme} column is present after cleaning
#'
#' hydroweight_attributes() only reports columns for classes it actually
#' sees occurring in the ROI — a class with zero pixels present is simply
#' omitted, not reported as 0. Harmless for a single site processed alone,
#' but once results from multiple sites/layers are bound together (which
#' calculate_hydroweight_attributes_stream() always does), a column absent
#' for one row and present for another becomes NA after dplyr::bind_rows()
#' — indistinguishable from genuinely missing data. Fill any such column
#' explicitly with 0 (the class is *known* to be absent, not unmeasured).
#'
#' @param tbl          Cleaned tibble from clean_categorical_columns_stream()
#' @param loi_name     Character. Prefix passed to
#'   clean_categorical_columns_stream() for this call (already includes any
#'   layer-specific suffix, e.g. "harvest_regen_y2015").
#' @param class_levels data.frame with ID/Class columns, or NULL (no-op).
#' @param scheme_names Character vector of weighting scheme names (e.g.
#'   names(hw)).
#' @return tbl with any missing expected column added, filled with 0.
ensure_full_categorical_schema <- function(tbl, loi_name, class_levels, scheme_names) {
  if (is.null(class_levels)) {
    return(tbl)
  }

  expected <- tidyr::expand_grid(
    cls = as.character(class_levels$Class),
    scheme = scheme_names
  ) |>
    dplyr::mutate(col = paste0(loi_name, "_", cls, "_", scheme, "_prop")) |>
    dplyr::pull(col) |>
    snake_case_cols()

  missing <- setdiff(expected, names(tbl))
  if (length(missing) > 0) {
    tbl[missing] <- 0
  }
  tbl
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
