# 05_hydroweight_attributes.R
# ---------------------------------------------------------------------------
# Computes inverse distance-weighted catchment attributes for lake sites using
# the hydroweight package.
#
# Supports two types of LOI (layer of interest) rasters:
#   Categorical (e.g., land cover): returns per-class proportions under each
#     weighting scheme via hydroweight_attributes(..., loi_numeric = FALSE)
#   Continuous (e.g., NDVI, elevation): returns summary statistics under each
#     weighting scheme via hydroweight_attributes(..., loi_numeric = TRUE).
#     Multi-layer rasters (e.g., annual NDVI stacks) are also supported.
#
# All seven weighting schemes are available:
#   lumped   — equal weighting (standard catchment mean)
#   iEucO    — inverse Euclidean distance to lake
#   iFLO     — inverse flow-path distance to lake
#   HAiFLO   — hydrologically-active inverse flow-path distance to lake
#   iEucS    — inverse Euclidean distance to nearest stream
#   iFLS     — inverse flow-path distance to nearest stream
#   HAiFLS   — hydrologically-active inverse flow-path distance to stream
#
# Stream-based schemes (iEucS, iFLS, HAiFLS) require a stream raster. Sites
# missing stream data fall back to lake-only schemes automatically.
#
# Run after Stage 6 (catchment metrics). Requires:
#   - Clipped catchment polygons: output_dir/<site_id>/catchment_clipped.gpkg
#   - Lake polygons: polys_all from match_lake_polygons() with matched_lake col
#   - Breached DEM:             cache_dir/dem_breached.tif
#   - Flow accumulation raster: cache_dir/flow_accum.tif
#
# Dependencies: terra, sf, hydroweight, dplyr, purrr, fs, cli (via utils.R)
# ---------------------------------------------------------------------------

#' Compute distance-weighted catchment attributes for all lake sites
#'
#' @param sites            tibble with site_id and lake_name columns
#' @param output_dir       Character. Root output directory
#' @param cache_dir        Character. Project-level cache directory
#' @param loi_layers       Named list of LOI layer descriptors. Each element is
#'   a list with:
#'   \describe{
#'     \item{path}{Character. Path to raster file (single or multi-layer .tif)}
#'     \item{name}{Character. Short label used as column prefix in output}
#'     \item{type}{"categorical" or "continuous"}
#'     \item{class_levels}{data.frame with ID (integer) and Class (character)
#'       columns. For categorical rasters only: replaces numeric class IDs in
#'       output column names with readable class names. Optional.}
#'     \item{stats}{Character vector. For continuous rasters: which statistic
#'       columns to keep in the final table. NULL keeps all. Optional.}
#'   }
#' @param weighting_schemes Character vector. Weighting schemes to compute.
#'   Stream-based schemes (iEucS, iFLS, HAiFLS) are silently dropped for sites
#'   without a stream raster. Default: all seven schemes.
#' @param inv_function     Function. Maps distance in metres to a weight.
#'   Default: (d_km + 1)^-1, giving weight 1.0 at 0 m, 0.5 at 1 km.
#' @param lake_polys       SpatVector from match_lake_polygons(). Used as
#'   target_O (lake outflow polygon for flow-path distance calculations).
#'   Required.
#' @param streams_path     One of:
#'   \describe{
#'     \item{NULL}{Auto-detect per site: prefers streams_clipped.tif in the
#'       site output folder, then falls back to streams.tif. No stream target
#'       if neither exists.}
#'     \item{Character path}{A single project-level stream raster used for all
#'       sites (extended to DEM extent automatically).}
#'     \item{NA}{Disable target_S entirely; only lake-based schemes are run.}
#'   }
#' @param hw_dir           Character. Directory for hydroweight intermediate
#'   rasters. Default: output_dir/hydroweight
#' @param dem_path         Character. Path to breached DEM. Default:
#'   cache_dir/dem_breached.tif
#' @param flow_accum_path  Character. Path to flow accumulation raster.
#'   Default: cache_dir/flow_accum.tif
#'
#' @return A tibble with one row per successfully processed site and one column
#'   per LOI attribute × weighting scheme combination, plus a leading site
#'   column.
calculate_hydroweight_attributes <- function(
  sites,
  output_dir,
  cache_dir,
  loi_layers,
  weighting_schemes = c("lumped", "iEucO", "iFLO", "HAiFLO", "iEucS", "iFLS", "HAiFLS"),
  inv_function      = function(x) (x * 0.001 + 1)^-1,
  lake_polys        = NULL,
  streams_path      = NULL,
  hw_dir            = fs::path(output_dir, "hydroweight"),
  dem_path          = NULL,
  flow_accum_path   = NULL
) {
  if (!requireNamespace("hydroweight", quietly = TRUE)) {
    cw_abort(paste(
      "Package 'hydroweight' is required.",
      "Install from GitHub: remotes::install_github('bencuddihee/hydroweight')"
    ))
  }
  if (is.null(lake_polys)) {
    cw_abort("lake_polys is required — pass the SpatVector from match_lake_polygons().")
  }
  if (!is.list(loi_layers) || length(loi_layers) == 0) {
    cw_abort("loi_layers must be a non-empty list of LOI layer descriptors.")
  }

  dem_path        <- dem_path        %||% fs::path(cache_dir, "dem_breached.tif")
  flow_accum_path <- flow_accum_path %||% fs::path(cache_dir, "flow_accum.tif")

  for (p in c(dem_path, flow_accum_path)) {
    if (!cache_exists(p)) {
      cw_abort(glue::glue("Required raster not found: {p}"))
    }
  }

  # Validate LOI layer descriptors up front
  validate_loi_layers(loi_layers)

  # Prepare (reproject/cache) LOI rasters once — expensive step done globally
  loi_prepared <- prepare_loi_layers(loi_layers, cache_dir)

  fs::dir_create(hw_dir, recurse = TRUE)

  cw_inform(glue::glue(
    "hydroweight: {nrow(sites)} site(s), {length(loi_layers)} LOI layer(s)"
  ))

  results <- purrr::map(sites$site_id, function(sid) {
    lake_name <- sites$lake_name[sites$site_id == sid]

    lake_poly_sf <- extract_lake_poly(lake_polys, lake_name, sid)
    if (is.null(lake_poly_sf)) return(NULL)

    effective_streams <- resolve_hw_streams_path(streams_path, output_dir, sid)

    process_hw_site(
      site_id           = sid,
      lake_name         = lake_name,
      output_dir        = output_dir,
      hw_dir            = hw_dir,
      lake_poly_sf      = lake_poly_sf,
      loi_prepared      = loi_prepared,
      loi_layers        = loi_layers,
      weighting_schemes = weighting_schemes,
      inv_function      = inv_function,
      streams_path      = effective_streams,
      dem_path          = dem_path,
      flow_accum_path   = flow_accum_path
    )
  })

  names(results) <- sites$site_id
  n_failed <- sum(vapply(results, is.null, logical(1)))
  if (n_failed > 0) {
    failed <- names(results)[vapply(results, is.null, logical(1))]
    cw_warn(glue::glue(
      "{n_failed} site(s) failed: {paste(failed, collapse = ', ')}"
    ))
  }

  results_clean <- results[!vapply(results, is.null, logical(1))]
  if (length(results_clean) == 0) {
    cw_warn("No sites produced hydroweight results.")
    return(tibble::tibble())
  }

  cw_inform(glue::glue(
    "hydroweight complete: {length(results_clean)}/{nrow(sites)} sites processed."
  ))

  dplyr::bind_rows(results_clean)
}

# -- LOI validation and preparation -------------------------------------------

#' Validate loi_layers list structure
validate_loi_layers <- function(loi_layers) {
  for (i in seq_along(loi_layers)) {
    lyr <- loi_layers[[i]]
    lbl <- if (!is.null(lyr$name)) lyr$name else paste0("layer ", i)
    if (is.null(lyr$path)) {
      cw_abort(glue::glue("loi_layers[[{i}]] ({lbl}): 'path' is required."))
    }
    if (!cache_exists(lyr$path)) {
      cw_abort(glue::glue("loi_layers[[{i}]] ({lbl}): file not found — {lyr$path}"))
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

#' Prepare LOI rasters: reproject to RASTER_CRS and cache to disk
#'
#' Each LOI raster is reprojected once to EPSG:3161 and written to a permanent
#' cache file to avoid terra temp file cleanup issues during per-site processing.
#' Multi-layer rasters are supported; all layers are written to a single cached
#' file.
#'
#' @param loi_layers Validated list of LOI layer descriptors
#' @param cache_dir  Character. Project cache directory
#' @param raster_crs Character. Target CRS. Default: EPSG:3161
#' @return Named list of SpatRasters (one per LOI layer, in raster_crs)
prepare_loi_layers <- function(loi_layers, cache_dir, raster_crs = "EPSG:3161") {
  hw_cache_dir <- fs::path(cache_dir, "hydroweight_loi")
  fs::dir_create(hw_cache_dir, recurse = TRUE)

  purrr::map(loi_layers, function(lyr) {
    cached_path <- fs::path(hw_cache_dir, paste0(lyr$name, "_reprojected.tif"))

    if (!cache_exists(cached_path)) {
      cw_inform(glue::glue("Preparing LOI '{lyr$name}' — reprojecting to {raster_crs}..."))
      raw <- terra::rast(lyr$path)
      if (!terra::same.crs(raw, raster_crs)) {
        method <- if (lyr$type == "categorical") "near" else "bilinear"
        raw <- terra::project(raw, raster_crs, method = method)
      }
      dtype <- if (lyr$type == "categorical") "INT1U" else "FLT4S"
      terra::writeRaster(raw, cached_path, overwrite = TRUE, datatype = dtype)
      cw_inform(glue::glue("  -> Cached: {cached_path}"))
    }

    r <- terra::rast(cached_path)

    # Apply factor levels for categorical rasters so hydroweight_attributes()
    # uses class labels rather than numeric IDs in output column names
    if (lyr$type == "categorical") {
      r <- terra::as.factor(r)
      if (!is.null(lyr$class_levels)) {
        levels(r) <- lyr$class_levels
      }
    }

    r
  }) |> setNames(purrr::map_chr(loi_layers, "name"))
}

# -- Per-site processing -------------------------------------------------------

#' Extract and convert lake polygon for one site to sf in EPSG:3161
#'
#' @param lake_polys SpatVector with matched_lake column
#' @param lake_name  Character. Site lake name to look up
#' @param site_id    Character. Site identifier (for log messages)
#' @return sf polygon in EPSG:3161, or NULL if not found
extract_lake_poly <- function(lake_polys, lake_name, site_id) {
  poly <- lake_polys[lake_polys$matched_lake == lake_name, ]
  if (nrow(poly) == 0) {
    cw_warn(glue::glue("Site '{site_id}': no lake polygon for '{lake_name}' — skipping"))
    return(NULL)
  }
  poly |> sf::st_as_sf() |> sf::st_transform("EPSG:3161")
}

#' Resolve the effective streams raster path for a site
#'
#' @param streams_path User-supplied streams_path argument
#' @param output_dir   Character. Root output directory
#' @param site_id      Character. Site identifier
#' @return Character path, or NULL if disabled/not found
resolve_hw_streams_path <- function(streams_path, output_dir, site_id) {
  if (isTRUE(is.na(streams_path))) return(NULL)
  if (!is.null(streams_path)) return(streams_path)

  site_dir <- site_output_dir(output_dir, site_id)
  clipped  <- fs::path(site_dir, "streams_clipped.tif")
  raw      <- fs::path(site_dir, "streams.tif")
  oih      <- fs::path(site_dir, paste0(site_id, "_streams_oih.tif"))

  if (cache_exists(clipped)) return(clipped)
  if (cache_exists(oih))     return(oih)
  if (cache_exists(raw))     return(raw)
  NULL
}

#' Process hydroweight for a single site
#'
#' @param site_id           Character. Site identifier
#' @param lake_name         Character. Lake name for this site
#' @param output_dir        Character. Root output directory
#' @param hw_dir            Character. Hydroweight intermediate raster directory
#' @param lake_poly_sf      sf polygon. Target lake (unbuffered), EPSG:3161
#' @param loi_prepared      Named list of SpatRasters from prepare_loi_layers()
#' @param loi_layers        Original LOI layer descriptor list (for metadata)
#' @param weighting_schemes Character vector of schemes to run
#' @param inv_function      Inverse distance function
#' @param streams_path      Character path or NULL
#' @param dem_path          Character. Path to breached DEM
#' @param flow_accum_path   Character. Path to flow accumulation raster
#' @return Single-row tibble with site column and all LOI attributes, or NULL
process_hw_site <- function(
  site_id,
  lake_name,
  output_dir,
  hw_dir,
  lake_poly_sf,
  loi_prepared,
  loi_layers,
  weighting_schemes,
  inv_function,
  streams_path,
  dem_path,
  flow_accum_path
) {
  cw_inform(glue::glue("\nhw: {site_id}"))
  hw_site_dir <- fs::path(hw_dir, site_id)
  fs::dir_create(hw_site_dir, recurse = TRUE)

  site_dir <- site_output_dir(output_dir, site_id)
  catchment_path <- fs::path(site_dir, "catchment_clipped.gpkg")

  if (!cache_exists(catchment_path)) {
    cw_warn(glue::glue("Site '{site_id}': catchment_clipped.gpkg not found — skipping"))
    return(NULL)
  }

  site_catch_sf <- sf::st_read(catchment_path, quiet = TRUE) |>
    sf::st_transform("EPSG:3161")

  # Load and extend stream raster to full DEM extent (Bug 4 fix)
  # hydroweight() needs target_S to span the same domain as dem/flow_accum so
  # it can compute flow-path distances to streams from all catchment cells.
  target_S <- NULL
  schemes_final <- weighting_schemes
  stream_schemes <- c("iEucS", "iFLS", "HAiFLS")

  if (!is.null(streams_path)) {
    if (!cache_exists(streams_path)) {
      cw_warn(glue::glue(
        "Site '{site_id}': stream raster not found ({streams_path}) — ",
        "dropping stream-based schemes"
      ))
      schemes_final <- setdiff(weighting_schemes, stream_schemes)
    } else {
      target_S <- tryCatch({
        s <- terra::rast(streams_path)
        # Extend to full DEM extent so hydroweight() can trace flow paths
        # from all catchment cells to the nearest stream
        dem_template <- terra::rast(dem_path)
        terra::extend(s, dem_template)
      }, error = function(e) {
        cw_warn(glue::glue(
          "Site '{site_id}': failed to load stream raster — {e$message}. ",
          "Dropping stream-based schemes."
        ))
        NULL
      })
      if (is.null(target_S)) {
        schemes_final <- setdiff(weighting_schemes, stream_schemes)
      } else {
        cw_inform(glue::glue(
          "  Stream raster loaded: ",
          "{sum(!is.na(terra::values(terra::crop(target_S, terra::vect(site_catch_sf)))))} ",
          "stream cells in catchment"
        ))
      }
    }
  } else {
    schemes_final <- setdiff(weighting_schemes, stream_schemes)
  }

  if (length(schemes_final) == 0) {
    cw_warn(glue::glue("Site '{site_id}': no weighting schemes available — skipping"))
    return(NULL)
  }

  # Run hydroweight() — generates one distance-weight raster per scheme
  cw_inform(glue::glue(
    "  Running hydroweight() [{paste(schemes_final, collapse=', ')}]..."
  ))

  hw <- tryCatch(
    hydroweight::hydroweight(
      hydroweight_dir  = hw_site_dir,
      target_O         = lake_poly_sf,
      target_S         = target_S,
      target_uid       = site_id,
      clip_region      = site_catch_sf,
      OS_combine       = !is.null(target_S),
      dem              = dem_path,
      flow_accum       = flow_accum_path,
      weighting_scheme = schemes_final,
      inv_function     = inv_function
    ),
    error = function(e) {
      cw_warn(glue::glue("Site '{site_id}': hydroweight() failed — {e$message}"))
      NULL
    }
  )

  if (is.null(hw)) return(NULL)

  # Unpack PackedSpatRaster → SpatRaster (required in hydroweight >= 2.0)
  hw <- lapply(hw, terra::rast)
  cw_inform(glue::glue(
    "  Distance weights: {paste(names(hw), collapse=', ')}"
  ))

  # Run hydroweight_attributes() for each LOI layer and column-bind
  attr_tables <- purrr::map(seq_along(loi_prepared), function(i) {
    lyr_name <- names(loi_prepared)[i]
    lyr_rast <- loi_prepared[[i]]
    lyr_desc <- loi_layers[[i]]

    cw_inform(glue::glue("  LOI '{lyr_name}'..."))

    run_loi_attributes(
      loi_rast     = lyr_rast,
      loi_desc     = lyr_desc,
      site_catch   = site_catch_sf,
      lake_poly    = lake_poly_sf,
      hw           = hw,
      site_id      = site_id
    )
  })

  # Drop failed LOI results, column-bind survivors
  attr_tables <- attr_tables[!vapply(attr_tables, is.null, logical(1))]
  if (length(attr_tables) == 0) return(NULL)

  result <- dplyr::bind_cols(attr_tables)
  dplyr::mutate(result, site = site_id, .before = dplyr::everything())
}

# -- LOI attribute computation ------------------------------------------------

#' Compute hydroweight attributes for one LOI layer at one site
#'
#' @param loi_rast   SpatRaster. Prepared LOI raster (in EPSG:3161)
#' @param loi_desc   List. LOI descriptor from loi_layers (path, name, type, ...)
#' @param site_catch sf polygon. Site catchment (EPSG:3161)
#' @param lake_poly  sf polygon. Target lake (EPSG:3161) — excluded from LOI summary
#' @param hw         Named list of SpatRasters. Distance weight rasters from hydroweight()
#' @param site_id    Character. Site identifier (for log messages)
#' @return Single-row tibble of LOI attributes, or NULL on failure
run_loi_attributes <- function(
  loi_rast,
  loi_desc,
  site_catch,
  lake_poly,
  hw,
  site_id
) {
  loi_name     <- loi_desc$name
  loi_numeric  <- (loi_desc$type == "continuous")

  # Crop and mask to catchment extent
  catch_vect <- terra::vect(site_catch)
  loi_site <- tryCatch({
    loi_rast |>
      terra::crop(catch_vect) |>
      terra::mask(catch_vect)
  }, error = function(e) {
    cw_warn(glue::glue(
      "Site '{site_id}', LOI '{loi_name}': crop/mask failed — {e$message}"
    ))
    NULL
  })

  if (is.null(loi_site)) return(NULL)
  if (all(is.na(terra::values(loi_site)))) {
    cw_warn(glue::glue(
      "Site '{site_id}', LOI '{loi_name}': all NA after crop/mask — skipping"
    ))
    return(NULL)
  }

  # hydroweight_attributes() uses pivot_wider() internally and cannot handle
  # multi-layer rasters reliably — it produces duplicate name errors regardless
  # of layer naming. For continuous multi-layer rasters (e.g. annual NDVI
  # stacks), call hydroweight_attributes() once per layer and bind columns.
  if (loi_numeric && terra::nlyr(loi_site) > 1) {
    return(run_loi_attributes_multilayer(
      loi_site   = loi_site,
      loi_desc   = loi_desc,
      site_catch = site_catch,
      lake_poly  = lake_poly,
      hw         = hw,
      site_id    = site_id
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
      remove_region    = lake_poly,
      return_products  = FALSE
    ),
    error = function(e) {
      cw_warn(glue::glue(
        "Site '{site_id}', LOI '{loi_name}': hydroweight_attributes() failed — {e$message}"
      ))
      NULL
    }
  )

  if (is.null(hwa)) return(NULL)

  attr_tbl <- hwa$attribute_table |>
    dplyr::select(-dplyr::any_of(c(site_id, "1")))

  # Rename and clean columns by type
  if (!loi_numeric) {
    attr_tbl <- clean_categorical_columns(attr_tbl, loi_name, loi_desc$class_levels)
  } else {
    attr_tbl <- clean_continuous_columns(attr_tbl, loi_name, loi_desc$stats)
  }

  attr_tbl
}

#' Process a multi-layer continuous LOI one layer at a time
#'
#' hydroweight_attributes() uses pivot_wider() internally and cannot correctly
#' pivot multi-layer rasters — it produces duplicate column errors. This helper
#' calls hydroweight_attributes() once per layer and column-binds the results.
#'
#' @param loi_site   SpatRaster. Cropped/masked LOI (multiple layers)
#' @param loi_desc   List. LOI descriptor (name, type, stats, ...)
#' @param site_catch sf polygon. Site catchment
#' @param lake_poly  sf polygon. Target lake
#' @param hw         Named list of SpatRasters. Distance weight rasters
#' @param site_id    Character. Site identifier
#' @return Single-row tibble of all layer attributes, or NULL on failure
run_loi_attributes_multilayer <- function(
  loi_site,
  loi_desc,
  site_catch,
  lake_poly,
  hw,
  site_id
) {
  loi_name    <- loi_desc$name
  layer_names <- names(loi_site)

  # Enforce unique layer names — use loi_name_N if any are duplicated or blank
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
        remove_region    = lake_poly,
        return_products  = FALSE
      ),
      error = function(e) {
        cw_warn(glue::glue(
          "Site '{site_id}', LOI '{loi_name}' layer {k}: ",
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
  clean_continuous_columns(combined, loi_name, loi_desc$stats)
}

# -- Column cleaning ----------------------------------------------------------

#' Rename categorical LOI columns from Class_{id}_{scheme} to {class}_{scheme}
#'
#' @param tbl          tibble. Raw output from hydroweight_attributes()
#' @param loi_name     Character. LOI name prefix
#' @param class_levels data.frame with ID and Class columns, or NULL
#' @return tibble with cleaned column names
clean_categorical_columns <- function(tbl, loi_name, class_levels) {
  if (!is.null(class_levels)) {
    levels_map <- setNames(
      as.character(class_levels$Class),
      as.character(class_levels$ID)
    )
    tbl <- dplyr::rename_with(tbl, function(nms) {
      id   <- regmatches(nms, regexpr("(?<=^Class_)\\d+", nms, perl = TRUE))
      rest <- sub("^Class_\\d+_", "", nms)
      cls  <- levels_map[id]
      # Replace Class_{id}_ prefix; non-matching columns pass through unchanged
      changed <- !is.na(cls) & nzchar(cls)
      nms[changed] <- paste0(loi_name, "_", cls[changed], "_", rest[changed])
      nms
    })
  } else {
    # No class_levels supplied — just prefix with loi_name
    tbl <- dplyr::rename_with(tbl, ~ paste0(loi_name, "_", .x))
  }

  # Sanitize to snake_case (removes spaces/hyphens, lowercases)
  dplyr::rename_with(tbl, ~ gsub("[^A-Za-z0-9_]", "_", .x) |>
    gsub("_+", "_", x = _) |>
    tolower() |>
    sub("_$", "", x = _))
}

#' Prefix continuous LOI columns with loi_name and optionally filter stats
#'
#' @param tbl      tibble. Raw output from hydroweight_attributes()
#' @param loi_name Character. LOI name prefix
#' @param stats    Character vector of stat substrings to keep, or NULL for all
#' @return tibble with prefixed column names
clean_continuous_columns <- function(tbl, loi_name, stats) {
  tbl <- dplyr::rename_with(tbl, ~ paste0(loi_name, "_", .x))

  if (!is.null(stats)) {
    keep <- grep(paste(stats, collapse = "|"), names(tbl), value = TRUE)
    tbl <- dplyr::select(tbl, dplyr::any_of(keep))
  }

  tbl
}

# -- Null coalescing operator (mirrors standard R 4.4+ behaviour) ----
`%||%` <- function(x, y) if (!is.null(x)) x else y
