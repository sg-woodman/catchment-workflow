# 02_prepare_terrain.R
# ---------------------------------------------------------------------------
# Conditional terrain preparation, per group. Materializes the SAME
# canonical set of per-group cache files every existing reused function
# already expects (workflow/R/stream/05_delineate_sites.R's
# load_group_rasters(), workflow/R/stream/06_hydroweight_attributes.R's
# dem_breached.tif/flow_accum.tif lookup) regardless of which terrain tier
# was supplied — so those functions can be reused completely unmodified,
# for both grouping strategies:
#
#   cache/<group>/dem.tif            best-available elevation surface
#   cache/<group>/dem_burned.tif     (only if streams were burned in)
#   cache/<group>/dem_breached.tif   depression-breached DEM, OR an alias
#                                     of the best available DEM when the
#                                     terrain tier is already conditioned
#                                     (no real breaching happened, but the
#                                     filename invariant is kept so every
#                                     reused function finds what it expects)
#   cache/<group>/flow_pointer.tif   WhiteboxTools-encoded D8 pointer
#   cache/<group>/flow_accum.tif     D8 flow accumulation (cells)
#   cache/<group>/streams.tif        binary stream raster (pour-point
#                                     snapping + hydroweight target_S)
#   cache/<group>/hillshade.tif      cheap, kept for load_group_rasters()
#                                     compatibility
#
# Terrain tier (resolved in 00_resolve_config.R, highest-conditioning wins):
#   flow_pointer   supplied → crop, use directly. No breach.
#   flow_direction supplied → crop, recode to WBT encoding (config$flow_direction$recode
#                    matrix, e.g. OIH's rotation table — generic here, not
#                    OIH-specific). No breach.
#   dem            (raw) only → crop, run 02_prepare_streams_burn.R's
#                    resolve_streams_burn() [config$streams_burn], then
#                    wbt_breach_depressions_least_cost() + wbt_d8_pointer().
#
# For every tier, `dem` — if separately supplied (e.g. CAM streams' OIH
# case: flow_direction for conditioning, dem for the elevation surface used
# in per-site clipping/output/hydroweight) — is cropped and cached as
# dem.tif regardless of which tier drove pointer generation. Cropping to
# each group's AOI is applied uniformly for both grouping strategies —
# for "whole_domain" the AOI already equals the source's own full extent
# (see 01_build_group_manifest.R::whole_domain_aoi()), so the crop is a
# harmless identity operation; no strategy branching needed here.
#
# Dependencies: terra, whitebox, fs, glue, cli (via utils.R)
# ---------------------------------------------------------------------------

#' Prepare terrain products for every group in the manifest
#'
#' @param config         Resolved config from resolve_engine_config()
#' @param group_manifest sf tibble from build_engine_group_manifest()
#' @return group_manifest, invisibly. Called for side effects.
prepare_engine_terrain <- function(config, group_manifest) {
  purrr::walk(seq_len(nrow(group_manifest)), function(i) {
    grp       <- group_manifest$group_id[i]
    grp_cache <- group_manifest$cache_dir[i]
    aoi       <- group_manifest$aoi[i]
    burn      <- group_manifest$burn_streams[i]

    ensure_dir(grp_cache)

    pointer_path  <- fs::path(grp_cache, "flow_pointer.tif")
    breached_path <- fs::path(grp_cache, "dem_breached.tif")
    accum_path    <- fs::path(grp_cache, "flow_accum.tif")
    streams_path  <- fs::path(grp_cache, "streams.tif")
    hillshade_path <- fs::path(grp_cache, "hillshade.tif")
    dem_path      <- fs::path(grp_cache, "dem.tif")

    if (cache_exists(pointer_path) && cache_exists(breached_path) &&
      cache_exists(accum_path) && cache_exists(streams_path)) {
      cw_inform(glue::glue("Group '{grp}': terrain products found in cache, skipping."))
      return(invisible(NULL))
    }

    cw_inform(glue::glue("Group '{grp}': preparing terrain (tier = '{config$terrain_tier}')..."))

    # -- dem.tif: cache the best-available elevation surface if supplied,
    # independent of which tier drives conditioning (decoupled on purpose —
    # e.g. CAM streams' flow_direction supersedes dem for conditioning, but
    # dem is still what per-site output/hydroweight clip against).
    if (!is.null(config$dem[["path"]]) && !cache_exists(dem_path)) {
      crop_and_reproject_source(
        source_path = config$dem[["path"]], aoi = aoi,
        target_crs = config$working_crs, out_path = dem_path,
        method = "bilinear", datatype = "FLT4S", group_id = grp, label = "dem.tif"
      )
    }

    if (config$terrain_tier == "flow_pointer") {
      if (!cache_exists(pointer_path)) {
        crop_and_reproject_source(
          source_path = config$flow_pointer[["path"]], aoi = aoi,
          target_crs = config$working_crs, out_path = pointer_path,
          method = "near", datatype = "INT1U", group_id = grp, label = "flow_pointer.tif"
        )
      }
      alias_dem_breached(dem_path, breached_path, grp)
    } else if (config$terrain_tier == "flow_direction") {
      prepare_flow_direction_tier(config, aoi, grp_cache, grp)
      alias_dem_breached(dem_path, breached_path, grp)
    } else {
      prepare_raw_dem_tier(config, aoi, grp_cache, grp, burn)
    }

    # -- flow_accum.tif: supplied → crop; else derive from flow_pointer.tif
    if (!cache_exists(accum_path)) {
      if (!is.null(config$flow_accum[["path"]])) {
        crop_and_reproject_source(
          source_path = config$flow_accum[["path"]], aoi = aoi,
          target_crs = config$working_crs, out_path = accum_path,
          method = "near", datatype = "FLT4S", group_id = grp, label = "flow_accum.tif"
        )
      } else {
        cw_inform(glue::glue("Group '{grp}': computing D8 flow accumulation..."))
        whitebox::wbt_d8_flow_accumulation(
          input    = normalizePath(pointer_path, mustWork = TRUE),
          output   = normalizePath(accum_path, mustWork = FALSE),
          out_type = "cells",
          pntr     = TRUE
        )
        if (!cache_exists(accum_path)) {
          cw_abort(glue::glue(
            "Group '{grp}': wbt_d8_flow_accumulation() did not produce output."
          ))
        }
      }
    }

    # -- streams.tif: always extracted from flow_accum (needed for
    # pour-point snapping in point mode, and as hydroweight's target_S in
    # both modes)
    if (!cache_exists(streams_path)) {
      cw_inform(glue::glue(
        "Group '{grp}': extracting streams (threshold = {config$stream_threshold} cells)..."
      ))
      whitebox::wbt_extract_streams(
        flow_accum = normalizePath(accum_path, mustWork = TRUE),
        output     = normalizePath(streams_path, mustWork = FALSE),
        threshold  = config$stream_threshold
      )
      if (!cache_exists(streams_path)) {
        cw_abort(glue::glue("Group '{grp}': wbt_extract_streams() did not produce output."))
      }
    }

    # -- hillshade.tif: cheap, kept only so load_group_rasters() (reused
    # unmodified from workflow/R/stream/05_delineate_sites.R) finds every
    # raster it expects.
    if (!cache_exists(hillshade_path) && cache_exists(breached_path)) {
      cw_inform(glue::glue("Group '{grp}': computing hillshade..."))
      whitebox::wbt_hillshade(
        dem    = normalizePath(breached_path, mustWork = TRUE),
        output = normalizePath(hillshade_path, mustWork = FALSE)
      )
    }

    cw_inform(glue::glue("Group '{grp}': terrain preparation complete."))
  })

  invisible(group_manifest)
}

# -- Per-tier helpers ---------------------------------------------------------

#' flow_direction tier: crop, then recode to WhiteboxTools D8 encoding
prepare_flow_direction_tier <- function(config, aoi, grp_cache, grp) {
  pointer_path <- fs::path(grp_cache, "flow_pointer.tif")
  if (cache_exists(pointer_path)) {
    return(invisible(pointer_path))
  }

  cropped <- crop_source_to_aoi(
    source_path = config$flow_direction[["path"]], aoi = aoi,
    target_crs = config$working_crs, method = "near", group_id = grp,
    label = "flow_direction"
  )

  recode <- config$flow_direction[["recode"]]
  if (!is.null(recode)) {
    cw_inform(glue::glue(
      "Group '{grp}': recoding external flow direction to WhiteboxTools encoding..."
    ))
    cropped <- terra::classify(cropped, recode, others = NA)
  }

  terra::writeRaster(
    cropped, pointer_path,
    overwrite = TRUE, datatype = "INT1U", gdal = c("COMPRESS=LZW")
  )
  cw_inform(glue::glue("Group '{grp}': flow_pointer.tif written (from flow_direction)."))
  invisible(pointer_path)
}

#' Raw dem tier: crop, optionally burn streams in, breach, D8 pointer
prepare_raw_dem_tier <- function(config, aoi, grp_cache, grp, burn) {
  dem_path      <- fs::path(grp_cache, "dem.tif")
  breached_path <- fs::path(grp_cache, "dem_breached.tif")
  pointer_path  <- fs::path(grp_cache, "flow_pointer.tif")

  if (!cache_exists(dem_path)) {
    crop_and_reproject_source(
      source_path = config$dem[["path"]], aoi = aoi,
      target_crs = config$working_crs, out_path = dem_path,
      method = "bilinear", datatype = "FLT4S", group_id = grp, label = "dem.tif"
    )
  }

  input_dem <- dem_path
  if (burn && config$streams_burn$source != "none") {
    dem_burned_path <- fs::path(grp_cache, "dem_burned.tif")
    if (!cache_exists(dem_burned_path)) {
      resolve_streams_burn(config, aoi, dem_path, dem_burned_path, grp)
    }
    if (cache_exists(dem_burned_path)) {
      input_dem <- dem_burned_path
    }
  }

  if (!cache_exists(breached_path)) {
    cw_inform(glue::glue("Group '{grp}': breaching depressions..."))
    whitebox::wbt_breach_depressions_least_cost(
      dem    = normalizePath(input_dem, mustWork = TRUE),
      output = normalizePath(breached_path, mustWork = FALSE),
      dist   = 10,
      fill   = TRUE
    )
    if (!cache_exists(breached_path)) {
      cw_abort(glue::glue("Group '{grp}': wbt_breach_depressions_least_cost() did not produce output."))
    }
  }

  if (!cache_exists(pointer_path)) {
    cw_inform(glue::glue("Group '{grp}': computing D8 flow pointer..."))
    whitebox::wbt_d8_pointer(
      dem    = normalizePath(breached_path, mustWork = TRUE),
      output = normalizePath(pointer_path, mustWork = FALSE)
    )
    if (!cache_exists(pointer_path)) {
      cw_abort(glue::glue("Group '{grp}': wbt_d8_pointer() did not produce output."))
    }
  }

  invisible(pointer_path)
}

#' Alias dem_breached.tif to the best-available dem.tif for tiers that skip
#' real breaching (flow_pointer/flow_direction already conditioned) — kept
#' so workflow/R/stream/06_hydroweight_attributes.R's hardcoded
#' `fs::path(grp_cache, "dem_breached.tif")` lookup always finds a file,
#' without implying any actual depression-breaching happened.
alias_dem_breached <- function(dem_path, breached_path, grp) {
  if (cache_exists(breached_path)) {
    return(invisible(breached_path))
  }
  if (!cache_exists(dem_path)) {
    cw_warn(glue::glue(
      "Group '{grp}': no dem.tif available to alias as dem_breached.tif — ",
      "dem-dependent downstream steps (hydroweight, dem.tif output) will ",
      "be unavailable for this group."
    ))
    return(invisible(NULL))
  }
  fs::file_copy(dem_path, breached_path, overwrite = TRUE)
  invisible(breached_path)
}

# -- Shared crop/reproject helper ---------------------------------------------

#' Crop a source raster to a group AOI, in the source's own native CRS
#' (cheap — GDAL reads only the requested pixels), without reprojecting.
#' Matched the retired workflow/R/stream/02_prepare_dem.R's crop step
#' exactly (same terra::crop(..., snap = "out") logic) — confirmed via
#' direct parity testing before that file was removed; see git history if
#' you need to compare against it again.
crop_source_to_aoi <- function(source_path, aoi, target_crs, method, group_id, label) {
  src <- terra::rast(source_path)
  aoi_vect   <- terra::vect(aoi)
  aoi_native <- terra::project(aoi_vect, terra::crs(src))
  aoi_ext    <- terra::ext(aoi_native)
  src_ext    <- terra::ext(src)

  x_overlap <- aoi_ext$xmin <= src_ext$xmax && aoi_ext$xmax >= src_ext$xmin
  y_overlap <- aoi_ext$ymin <= src_ext$ymax && aoi_ext$ymax >= src_ext$ymin
  if (!x_overlap || !y_overlap) {
    cw_abort(glue::glue(
      "Group '{group_id}': AOI does not overlap {label} extent ({source_path})."
    ))
  }

  cropped <- terra::crop(src, aoi_ext, snap = "out")

  src_crs_desc <- terra::crs(src, describe = TRUE)
  src_epsg <- paste0(src_crs_desc$authority, ":", src_crs_desc$code)
  if (!identical(src_epsg, target_crs)) {
    cw_inform(glue::glue(
      "Group '{group_id}': reprojecting {label} from {src_epsg} to {target_crs}..."
    ))
    cropped <- terra::project(cropped, target_crs, method = method)
  }

  cropped
}

#' Crop + reproject + write a source raster to a group cache file
crop_and_reproject_source <- function(
  source_path, aoi, target_crs, out_path, method, datatype, group_id, label
) {
  cropped <- crop_source_to_aoi(source_path, aoi, target_crs, method, group_id, label)
  terra::writeRaster(
    cropped, out_path,
    overwrite = TRUE, datatype = datatype,
    gdal = c("COMPRESS=LZW", "BIGTIFF=IF_SAFER")
  )
  cw_inform(glue::glue("Group '{group_id}': {label} written."))
  invisible(out_path)
}
