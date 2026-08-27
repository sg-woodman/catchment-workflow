# 08_catchment_metrics.R
# ---------------------------------------------------------------------------
# Calculates catchment morphometric metrics for all delineated sites and
# compiles results into a summary tibble and a gt reference table.
#
# Computes metrics for BOTH the unclipped catchment (catchment.gpkg, dem.tif
# from 05_delineate_sites.R) and the clipped catchment (catchment_clipped.gpkg,
# dem_clipped.tif from 07_reclip_outputs.R), producing two rows per site
# tagged by a "version" column ("unclipped" / "clipped"). If a site has no
# clipped outputs, a warning is issued and that row is omitted (the unclipped
# row is still produced).
#
# Metric categories:
#   Geometry  : area, perimeter, basin length, catchment length, mean width
#   Areal     : form factor, circularity ratio, elongation ratio,
#               compactness coefficient, shape index, lemniscate ratio
#   Relief    : basin relief, relief ratio, relative relief, ruggedness number
#               (uses actual drainage density when streams provided, else proxy),
#               dissection index, hypsometric integral, slope, aspect
#   Lake      : lake area, fetch, shoreline length, shoreline development,
#               sphericity, compactness, lake %, drainage-area ratio
#               (only computed when lake_polys supplied)
#   Stream    : drainage density, stream frequency
#               (only computed when streams are detected or supplied)
#
# Metrics follow:
#   Shekar, P.R. & Mathew, A. (2024). Morphometric analysis of watersheds.
#   Watershed Ecology and the Environment, 6, 13-25.
#   https://doi.org/10.1016/j.wsee.2023.12.001
#
#   CFS Catchment and Lake Metric Methods (internal methods document)
#
# Inputs:
#   sites       : Validated sites tibble from validate_sites()
#   output_dir  : Root output directory containing per-site folders
#   lake_polys  : Optional SpatVector from match_lake_polygons() with
#                 matched_lake column — enables lake metric category
#   streams_path: Optional path to a project-level stream layer (.gpkg or
#                 .tif). When NULL, per-site streams_clipped.* files are
#                 auto-detected in each site output directory.
#
# Run after 07_reclip_outputs.R (clipped rows will be skipped with a warning
# if 07 has not been run for a site).
#
# Dependencies: sf, terra, dplyr, purrr, tibble, readr, gt, fs, cli (via utils.R)
# ---------------------------------------------------------------------------

#' Calculate morphometric metrics for all delineated sites
#'
#' @param sites        Validated sites tibble from validate_sites()
#' @param output_dir   Character. Root output directory
#' @param lake_polys   Optional SpatVector from match_lake_polygons() with a
#'                     matched_lake column. Enables the Lake metric category.
#'                     Requires sites to have a lake_name column.
#' @param streams_path Optional character path to a project-level stream layer
#'                     (.gpkg vector or .tif raster). NULL (default) triggers
#'                     auto-detection of per-site streams_clipped.* files.
#'                     Set to NA to disable stream metrics entirely.
#' @return A tibble with up to two rows per site (version = "unclipped" /
#'   "clipped") containing all computed metrics
calculate_catchment_metrics <- function(
  sites,
  output_dir,
  lake_polys   = NULL,
  streams_path = NULL
) {
  cw_inform(glue::glue(
    "Calculating morphometric metrics for {nrow(sites)} site(s)..."
  ))
  if (!is.null(lake_polys))
    cw_inform("Lake metrics enabled (lake_polys supplied).")
  if (!is.null(streams_path) && !is.na(streams_path))
    cw_inform(glue::glue("Stream metrics: using project-level layer — {streams_path}"))
  else if (is.null(streams_path))
    cw_inform("Stream metrics: auto-detecting per-site streams_clipped.* files.")
  else
    cw_inform("Stream metrics: disabled (streams_path = NA).")

  has_lake_name <- "lake_name" %in% names(sites)

  results <- purrr::map(seq_len(nrow(sites)), function(i) {
    site    <- sites[i, ]
    sid     <- site$site_id
    site_dir <- site_output_dir(output_dir, sid)

    # Resolve lake polygon for this site
    lake_poly <- if (!is.null(lake_polys) && has_lake_name) {
      lp <- lake_polys[lake_polys$matched_lake == site$lake_name, ]
      if (nrow(lp) > 0) lp else NULL
    } else NULL

    site_rows <- list()

    # -- Unclipped version ----------------------------------------------------
    unclipped_row <- compute_site_metrics(
      site_dir      = site_dir,
      site_id       = sid,
      version       = "unclipped",
      catchment_file = "catchment.gpkg",
      dem_file       = "dem.tif",
      lake_poly      = lake_poly,
      streams_path   = streams_path
    )

    if (is.null(unclipped_row)) {
      cw_warn(glue::glue(
        "Site '{sid}': catchment.gpkg not found — site skipped entirely."
      ))
      return(NULL)
    }

    site_rows[["unclipped"]] <- unclipped_row

    # -- Clipped version -------------------------------------------------------
    clipped_catchment_path <- fs::path(site_dir, "catchment_clipped.gpkg")

    if (!cache_exists(clipped_catchment_path)) {
      cw_warn(glue::glue(
        "Site '{sid}': catchment_clipped.gpkg not found — ",
        "clipped metrics skipped (run 06/07 first if expected)."
      ))
    } else {
      clipped_row <- compute_site_metrics(
        site_dir       = site_dir,
        site_id        = sid,
        version        = "clipped",
        catchment_file = "catchment_clipped.gpkg",
        dem_file       = "dem_clipped.tif",
        lake_poly      = lake_poly,
        streams_path   = streams_path
      )
      if (!is.null(clipped_row)) {
        site_rows[["clipped"]] <- clipped_row
      }
    }

    dplyr::bind_rows(site_rows)
  }) |>
    purrr::compact() |>
    dplyr::bind_rows()

  n_sites <- dplyr::n_distinct(results$site_id)
  cw_inform(glue::glue(
    "Metrics computed for {n_sites} site(s) ({nrow(results)} row(s) total)."
  ))

  results
}

# -- Per-version site metrics --------------------------------------------------

#' Compute all metric categories for a single site/version combination
#'
#' @param site_dir       Character. Site output directory path
#' @param site_id        Character. Site identifier (for log messages)
#' @param version        Character. "unclipped" or "clipped"
#' @param catchment_file Character. Filename of the catchment polygon
#' @param dem_file       Character. Filename of the DEM
#' @param lake_poly      Optional SpatVector. Single lake polygon for this site
#' @param streams_path   Optional character path to project-level stream layer,
#'                       NULL to auto-detect per-site files, or NA to disable
#' @return A single-row tibble with site_id, version, and all metric columns,
#'   or NULL if catchment_file does not exist
compute_site_metrics <- function(
  site_dir,
  site_id,
  version,
  catchment_file,
  dem_file,
  lake_poly    = NULL,
  streams_path = NULL
) {
  catchment_path <- fs::path(site_dir, catchment_file)
  dem_path       <- fs::path(site_dir, dem_file)

  if (!cache_exists(catchment_path)) return(NULL)

  cw_inform(glue::glue(
    "Site '{site_id}' [{version}]: computing metrics from {catchment_file}..."
  ))

  catchment <- sf::st_read(catchment_path, quiet = TRUE)

  dem <- if (cache_exists(dem_path)) {
    terra::rast(dem_path)
  } else {
    cw_warn(glue::glue(
      "Site '{site_id}' [{version}]: {dem_file} not found — ",
      "relief metrics will be NA."
    ))
    NULL
  }

  # Auto-detect stream layer when streams_path = NULL
  effective_streams_path <- if (!is.null(streams_path)) {
    streams_path  # supplied (project-level path or NA)
  } else {
    # Look for per-site clipped streams (prefer vector, fall back to raster)
    if (version == "clipped") {
      gpkg <- fs::path(site_dir, "streams_clipped.gpkg")
      tif  <- fs::path(site_dir, "streams_clipped.tif")
      if (cache_exists(gpkg)) gpkg else if (cache_exists(tif)) tif else NULL
    } else {
      tif <- fs::path(site_dir, "streams.tif")
      if (cache_exists(tif)) tif else NULL
    }
  }

  # Compute metric categories
  geom   <- compute_geometry_metrics(catchment, site_id)
  areal  <- compute_areal_metrics(geom)
  stream <- compute_stream_metrics(geom$area_km2, effective_streams_path, catchment, site_id)
  relief <- compute_relief_metrics(
    dem, geom, site_id,
    drainage_density = stream$drainage_density_km_km2
  )
  lake   <- compute_lake_metrics(geom$area_km2, lake_poly, site_id)

  tibble::tibble(site_id = site_id, version = version) |>
    dplyr::bind_cols(geom, areal, relief, lake, stream)
}

# -- Geometry metrics ----------------------------------------------------------

#' Compute basic geometry metrics from catchment polygon
#'
#' All measurements in metres and km². Catchment is reprojected to EPSG:3979
#' if not already in that CRS to ensure metric units.
#'
#' @param catchment sf polygon. Catchment boundary
#' @param site_id   Character. Site identifier (for log messages)
#' @return Named tibble row with geometry metrics
compute_geometry_metrics <- function(catchment, site_id) {
  catchment <- sf::st_transform(catchment, 3979)

  # Area (km²) and perimeter (km)
  area_m2  <- as.numeric(sf::st_area(catchment))
  area_km2 <- area_m2 / 1e6
  perim_m  <- as.numeric(sf::st_length(sf::st_cast(
    sf::st_union(catchment), "MULTILINESTRING"
  )))
  perim_km <- perim_m / 1000

  # Basin length (Lb): diagonal of axis-aligned minimum bounding rectangle
  bbox     <- sf::st_bbox(catchment)
  width_m  <- bbox["xmax"] - bbox["xmin"]
  height_m <- bbox["ymax"] - bbox["ymin"]
  lb_m     <- sqrt(width_m^2 + height_m^2)
  lb_km    <- lb_m / 1000

  # Mean basin width (Wb = A / Lb)
  wb_km <- area_km2 / lb_km

  # Catchment length: longer side of the minimum rotated rectangle (MRR)
  # of the convex hull — distinct from lb_km which uses the axis-aligned bbox.
  catchment_length_m <- tryCatch({
    ch     <- sf::st_convex_hull(sf::st_union(catchment))
    mrr    <- sf::st_minimum_rotated_rectangle(ch)
    coords <- sf::st_coordinates(mrr)[1:4, 1:2]
    side1  <- sqrt(sum((coords[1, ] - coords[2, ])^2))
    side2  <- sqrt(sum((coords[2, ] - coords[3, ])^2))
    max(side1, side2)
  }, error = function(e) {
    cw_warn(glue::glue(
      "Site '{site_id}': catchment_length_m failed — {e$message}"
    ))
    NA_real_
  })

  tibble::tibble(
    area_km2           = round(area_km2, 4),
    perim_km           = round(perim_km, 4),
    lb_km              = round(lb_km, 4),
    wb_km              = round(wb_km, 4),
    catchment_length_m = round(catchment_length_m, 1)
  )
}

# -- Areal metrics -------------------------------------------------------------

#' Compute areal (shape) morphometric metrics
#'
#' @param geom Named tibble row from compute_geometry_metrics()
#' @return Named tibble row with areal metrics
compute_areal_metrics <- function(geom) {
  A  <- geom$area_km2
  P  <- geom$perim_km
  Lb <- geom$lb_km

  Rf <- A / Lb^2
  Rc <- (4 * pi * A) / P^2
  Re <- (2 / Lb) * sqrt(A / pi)
  Cc <- P / (2 * sqrt(pi * A))
  Sw <- Lb^2 / A
  k  <- Lb^2 / (4 * A)

  tibble::tibble(
    form_factor        = round(Rf, 4),
    circularity_ratio  = round(Rc, 4),
    elongation_ratio   = round(Re, 4),
    compactness_coeff  = round(Cc, 4),
    shape_index        = round(Sw, 4),
    lemniscate_ratio   = round(k, 4)
  )
}

# -- Relief metrics ------------------------------------------------------------

#' Compute relief morphometric metrics from the catchment DEM
#'
#' @param dem              SpatRaster or NULL. Site DEM clipped to catchment
#' @param geom             Named tibble row from compute_geometry_metrics()
#' @param site_id          Character. Site identifier (for log messages)
#' @param drainage_density Numeric or NULL/NA. Actual drainage density (km/km²)
#'                         from compute_stream_metrics(). When supplied and
#'                         non-NA, used for ruggedness_number; otherwise the
#'                         geometric proxy Lb/A is used.
#' @return Named tibble row with relief metrics
compute_relief_metrics <- function(dem, geom, site_id, drainage_density = NULL) {
  na_row <- tibble::tibble(
    elev_min_m           = NA_real_,
    elev_max_m           = NA_real_,
    elev_sd_m            = NA_real_,
    basin_relief_m       = NA_real_,
    relief_ratio         = NA_real_,
    relative_relief      = NA_real_,
    ruggedness_number    = NA_real_,
    dissection_index     = NA_real_,
    hypsometric_integral = NA_real_,
    mean_slope_deg       = NA_real_,
    mean_aspect_deg      = NA_real_,
    aspect_class         = NA_character_
  )

  if (is.null(dem)) return(na_row)

  vals    <- terra::values(dem, na.rm = TRUE)
  elev_min <- min(vals)
  elev_max <- max(vals)
  H        <- elev_max - elev_min

  Rh  <- (H / 1000) / geom$lb_km
  Rhp <- (H * 100) / (geom$perim_km * 1000)

  # Ruggedness number: use actual Dd when available, else geometric proxy
  Dd  <- if (!is.null(drainage_density) && !is.na(drainage_density)) {
    drainage_density
  } else {
    geom$lb_km / geom$area_km2
  }
  Rn <- (H / 1000) * Dd

  Dis <- H / elev_max

  mean_elev <- mean(vals)
  HI <- (mean_elev - elev_min) / (elev_max - elev_min)

  slope_rast   <- terra::terrain(dem, v = "slope",  unit = "degrees")
  mean_slope   <- mean(terra::values(slope_rast, na.rm = TRUE))

  aspect_rast  <- terra::terrain(dem, v = "aspect", unit = "degrees")
  aspect_vals  <- terra::values(aspect_rast, na.rm = TRUE)
  mean_sin     <- mean(sin(aspect_vals * pi / 180))
  mean_cos     <- mean(cos(aspect_vals * pi / 180))
  mean_aspect  <- (atan2(mean_sin, mean_cos) * 180 / pi) %% 360

  aspect_class <- dplyr::case_when(
    mean_aspect >= 337.5 | mean_aspect <  22.5 ~ "N",
    mean_aspect >=  22.5 & mean_aspect <  67.5 ~ "NE",
    mean_aspect >=  67.5 & mean_aspect < 112.5 ~ "E",
    mean_aspect >= 112.5 & mean_aspect < 157.5 ~ "SE",
    mean_aspect >= 157.5 & mean_aspect < 202.5 ~ "S",
    mean_aspect >= 202.5 & mean_aspect < 247.5 ~ "SW",
    mean_aspect >= 247.5 & mean_aspect < 292.5 ~ "W",
    mean_aspect >= 292.5 & mean_aspect < 337.5 ~ "NW",
    TRUE ~ NA_character_
  )

  elev_sd <- stats::sd(vals)

  tibble::tibble(
    elev_min_m           = round(elev_min,   2),
    elev_max_m           = round(elev_max,   2),
    elev_sd_m            = round(elev_sd,    2),
    basin_relief_m       = round(H,          2),
    relief_ratio         = round(Rh,         4),
    relative_relief      = round(Rhp,        4),
    ruggedness_number    = round(Rn,         4),
    dissection_index     = round(Dis,        4),
    hypsometric_integral = round(HI,         4),
    mean_slope_deg       = round(mean_slope, 2),
    mean_aspect_deg      = round(mean_aspect, 1),
    aspect_class         = aspect_class
  )
}

# -- Lake metrics --------------------------------------------------------------

#' Compute lake morphometric metrics from the matched OHN lake polygon
#'
#' Returns all-NA when lake_poly is NULL. The lake polygon is projected to
#' EPSG:3979 before computation so all results are in metric units.
#'
#' @param catchment_area_km2 Numeric. Catchment area (km²) from geometry metrics
#' @param lake_poly          SpatVector or NULL. Matched lake polygon for this site
#' @param site_id            Character. Site identifier (for log messages)
#' @return Named tibble row with lake metrics (all NA when lake_poly = NULL)
compute_lake_metrics <- function(catchment_area_km2, lake_poly, site_id) {
  na_row <- tibble::tibble(
    lake_area_ha          = NA_real_,
    lake_fetch_m          = NA_real_,
    shoreline_length_m    = NA_real_,
    shoreline_development = NA_real_,
    lake_sphericity       = NA_real_,
    lake_compactness      = NA_real_,
    lake_area_pct         = NA_real_,
    ca_la_ratio           = NA_real_
  )

  if (is.null(lake_poly) || nrow(lake_poly) == 0) return(na_row)

  tryCatch({
    lake_sf <- sf::st_as_sf(lake_poly) |> sf::st_transform(3979)

    lake_area_m2 <- sum(as.numeric(sf::st_area(lake_sf)))
    lake_area_ha <- lake_area_m2 / 10000

    perim_m <- as.numeric(sf::st_length(
      sf::st_cast(sf::st_union(lake_sf), "MULTILINESTRING")
    ))

    # Fetch: longer side of MRR of convex hull
    ch     <- sf::st_convex_hull(sf::st_union(lake_sf))
    mrr    <- sf::st_minimum_rotated_rectangle(ch)
    coords <- sf::st_coordinates(mrr)[1:4, 1:2]
    side1  <- sqrt(sum((coords[1, ] - coords[2, ])^2))
    side2  <- sqrt(sum((coords[2, ] - coords[3, ])^2))
    lake_fetch_m <- max(side1, side2)
    D_c    <- sqrt(side1^2 + side2^2)   # MRR diagonal ≈ circumscribed circle diameter

    # Shoreline development index (SDI) — 1 = perfect circle, increases with irregularity
    SDI <- perim_m / (2 * sqrt(pi * lake_area_m2))

    # Sphericity — ratio of equivalent-circle diameter to circumscribed circle diameter
    d_s <- 2 * sqrt(lake_area_m2 / pi)
    sphericity <- d_s / D_c

    # Lake compactness — isoperimetric ratio (complement to SDI²)
    compactness <- (4 * pi * lake_area_m2) / perim_m^2

    # Catchment–lake ratios
    catchment_area_ha <- catchment_area_km2 * 100
    lake_area_pct     <- (lake_area_ha / catchment_area_ha) * 100
    ca_la_ratio       <- catchment_area_ha / lake_area_ha

    tibble::tibble(
      lake_area_ha          = round(lake_area_ha,  2),
      lake_fetch_m          = round(lake_fetch_m,  1),
      shoreline_length_m    = round(perim_m,        1),
      shoreline_development = round(SDI,            4),
      lake_sphericity       = round(sphericity,     4),
      lake_compactness      = round(compactness,    4),
      lake_area_pct         = round(lake_area_pct,  2),
      ca_la_ratio           = round(ca_la_ratio,    2)
    )
  }, error = function(e) {
    cw_warn(glue::glue(
      "Site '{site_id}': lake metrics failed — {e$message}"
    ))
    na_row
  })
}

# -- Stream metrics ------------------------------------------------------------

#' Compute stream network metrics from a stream vector or raster layer
#'
#' Accepts either a vector path (.gpkg/.shp) or a raster path (.tif).
#' For vector input: clips to catchment and computes drainage density and
#' stream frequency. For raster input: counts stream cells (value == 1),
#' computes drainage density only (stream_frequency returns NA for rasters).
#'
#' @param catchment_area_km2 Numeric. Catchment area (km²) from geometry metrics
#' @param streams_path       Character path to stream layer, or NULL/NA to skip
#' @param catchment          sf polygon. Catchment boundary (any CRS)
#' @param site_id            Character. Site identifier (for log messages)
#' @return Named tibble row with drainage_density_km_km2 and
#'   stream_frequency_per_km2 (all NA when streams_path is NULL or NA)
compute_stream_metrics <- function(
  catchment_area_km2,
  streams_path,
  catchment,
  site_id
) {
  na_row <- tibble::tibble(
    drainage_density_km_km2  = NA_real_,
    stream_frequency_per_km2 = NA_real_
  )

  if (is.null(streams_path) || is.na(streams_path) || catchment_area_km2 <= 0) {
    return(na_row)
  }

  tryCatch({
    ext <- tolower(fs::path_ext(streams_path))

    if (ext %in% c("gpkg", "shp", "geojson")) {
      # -- Vector path ----------------------------------------------------------
      streams_raw <- sf::st_read(streams_path, quiet = TRUE)
      catchment_3979 <- sf::st_transform(sf::st_union(catchment), 3979)
      streams_sf <- streams_raw |>
        sf::st_transform(3979) |>
        sf::st_intersection(catchment_3979)
      # Keep only line geometries (intersection can produce points at boundaries)
      streams_sf <- streams_sf[
        sf::st_geometry_type(streams_sf) %in% c("LINESTRING", "MULTILINESTRING"),
        ,
        drop = FALSE
      ]

      if (nrow(streams_sf) == 0) {
        cw_warn(glue::glue(
          "Site '{site_id}': no stream features within catchment — stream metrics NA"
        ))
        return(na_row)
      }

      total_length_km  <- sum(as.numeric(sf::st_length(streams_sf))) / 1000
      drainage_density <- total_length_km / catchment_area_km2
      stream_frequency <- nrow(streams_sf) / catchment_area_km2

    } else if (ext == "tif") {
      # -- Raster path ----------------------------------------------------------
      r <- terra::rast(streams_path)
      catchment_vect <- terra::vect(
        sf::st_transform(sf::st_union(catchment), terra::crs(r))
      )
      r_clip <- terra::mask(terra::crop(r, catchment_vect), catchment_vect)

      # Stream cells have value 1; cell size gives approximate reach length
      stream_cells     <- sum(terra::values(r_clip) == 1, na.rm = TRUE)
      cell_size_m      <- sqrt(prod(terra::res(r_clip)))
      total_length_km  <- (stream_cells * cell_size_m) / 1000
      drainage_density <- total_length_km / catchment_area_km2
      stream_frequency <- NA_real_  # segment count not meaningful from raster

    } else {
      cw_warn(glue::glue(
        "Site '{site_id}': unrecognised stream file extension '.{ext}' — skipping"
      ))
      return(na_row)
    }

    tibble::tibble(
      drainage_density_km_km2  = round(drainage_density, 4),
      stream_frequency_per_km2 = round(stream_frequency, 4)
    )
  }, error = function(e) {
    cw_warn(glue::glue(
      "Site '{site_id}': stream metrics failed — {e$message}"
    ))
    na_row
  })
}

# -- Reference table -----------------------------------------------------------

#' Build a gt reference table of all morphometric metrics
#'
#' @return A gt table object
build_metrics_reference_table <- function() {
  ref <- tibble::tribble(
    ~category, ~parameter, ~symbol, ~formula, ~interpretation, ~source,

    # -- Geometry ---------------------------------------------------------------
    "Geometry", "Catchment Area", "A",
    "Planimetric area of catchment polygon (km²)",
    "Total area draining to the pour point. Larger catchments generally produce higher total discharge but lower peak runoff per unit area.",
    "—",

    "Geometry", "Perimeter", "P",
    "Length of catchment boundary (km)",
    "Longer perimeters relative to area indicate more irregular, elongated shapes.",
    "—",

    "Geometry", "Basin Length", "Lb",
    "Diagonal of axis-aligned minimum bounding rectangle (km)",
    "Approximates the longest axis of the basin using the bounding-box diagonal. Used in shape and relief ratio calculations.",
    "Schumm (1956)",

    "Geometry", "Mean Basin Width", "Wb",
    "Wb = A / Lb",
    "Average width perpendicular to the basin length. Narrow basins produce more peaked hydrographs.",
    "Horton (1932)",

    "Geometry", "Catchment Length", "—",
    "Longer side of minimum rotated rectangle of catchment convex hull (m)",
    "Longest physical dimension of the catchment. Computed from the tightest-fitting rectangle around the convex hull, giving the dominant elongation axis independent of orientation.",
    "—",

    # -- Areal ------------------------------------------------------------------
    "Areal", "Form Factor", "Rf",
    "Rf = A / Lb²",
    "Values near 1 indicate a circular basin with high, peaked flood discharge. Values < 0.45 indicate elongated basins with lower, flatter hydrographs and lower erosion potential.",
    "Horton (1932)",

    "Areal", "Circularity Ratio", "Rc",
    "Rc = 4πA / P²",
    "Range 0–1. Values approaching 1 indicate a circular, compact basin. Low values (<0.4) indicate elongated basins with low discharge and high permeability.",
    "Miller (1953)",

    "Areal", "Elongation Ratio", "Re",
    "Re = (2/Lb) × √(A/π)",
    "Range 0–1. Values near 1 indicate circular basins. <0.5 = highly elongated; 0.5–0.7 = elongated; 0.7–0.9 = oval; 0.9–1.0 = circular.",
    "Schumm (1956)",

    "Areal", "Compactness Coefficient", "Cc",
    "Cc = P / (2√(πA))",
    "Equals 1 for a perfect circle; increases with elongation. Higher values indicate more elongated basins with lower flood susceptibility.",
    "Gravelius (1914)",

    "Areal", "Shape Index", "Sw",
    "Sw = Lb² / A",
    "Inverse of form factor. Higher values indicate more elongated basins. Values < 2 = circular; 2–4 = slightly elongated; > 4 = strongly elongated.",
    "Horton (1932)",

    "Areal", "Lemniscate Ratio", "k",
    "k = Lb² / (4A)",
    "Values near 1 indicate circular basins. Values > 1 indicate elongated basins with lower peak discharge and reduced erosion susceptibility.",
    "Chorley et al. (1957)",

    # -- Relief -----------------------------------------------------------------
    "Relief", "Basin Relief", "H",
    "H = Elevₘₐₓ − Elevₘᴵₙ (m)",
    "Vertical distance between highest and lowest points. High relief indicates steep terrain, greater erosive energy, and higher peak flows.",
    "Hadley & Schumm (1961)",

    "Relief", "Relief Ratio", "Rh",
    "Rh = H (km) / Lb (km)",
    "Measures overall steepness of the watershed. Higher values indicate steeper terrain and greater intensity of erosion processes.",
    "Schumm (1956)",

    "Relief", "Relative Relief", "Rhp",
    "Rhp = H × 100 / P (m/km)",
    "Normalises relief by perimeter length. Higher values indicate more rugged, steeply incised terrain.",
    "Melton (1957)",

    "Relief", "Ruggedness Number", "Rn",
    "Rn = H (km) × Dd (km/km²)¹",
    "Combines relief and drainage density to indicate structural complexity and erosion potential. High values indicate steep, well-dissected terrain prone to high runoff and erosion. ¹Uses actual Dd from stream data when available; otherwise approximated as Lb/A.",
    "Patton & Baker (1976)",

    "Relief", "Dissection Index", "Dis",
    "Dis = H / Elevₘₐₓ",
    "Range 0–1. Measures the degree of terrain dissection relative to maximum elevation. Higher values indicate more deeply incised terrain.",
    "Singh & Dubey (1994)",

    "Relief", "Hypsometric Integral", "HI",
    "HI = (Elevₘₑₐₙ − Elevₘᴵₙ) / (Elevₘₐₓ − Elevₘᴵₙ)",
    "Range 0–1. HI > 0.6 = youthful stage (high erosion potential); 0.35–0.6 = mature stage; < 0.35 = old stage (highly eroded, low relief).",
    "Pike & Wilson (1971)",

    "Relief", "Mean Slope", "—",
    "Mean of per-cell slope (degrees) across catchment DEM",
    "Average gradient across the catchment. Larger slopes provide greater opportunity for surface runoff and nutrient transport.",
    "—",

    "Relief", "Mean Aspect", "—",
    "Circular mean of per-cell aspect (degrees, 0–360° from N)",
    "Dominant solar exposure direction. Affects snowmelt timing, evapotranspiration, and vegetation composition.",
    "—",

    "Relief", "Aspect Class", "—",
    "8-point compass class derived from mean aspect (N, NE, E, SE, S, SW, W, NW)",
    "Categorical summary of dominant catchment orientation.",
    "—",

    "Relief", "Elevation Min / Max / SD", "—",
    "Min, max, and standard deviation of DEM values within catchment (m)",
    "Elevation range and variability. Higher SD indicates more topographically complex, rugged terrain with greater potential for nutrient runoff.",
    "—",

    # -- Lake -------------------------------------------------------------------
    "Lake", "Lake Area", "—",
    "Lake polygon area (ha)",
    "Lake surface area. Larger lakes dilute inputs and may stratify more readily.",
    "Hutchinson (1957)",

    "Lake", "Lake Fetch", "—",
    "Longer side of minimum rotated rectangle of lake convex hull (m)",
    "Longest open-water distance. Larger fetch drives stronger wind mixing and wave action, affecting thermal stratification and nutrient cycling.",
    "Hutchinson (1957)",

    "Lake", "Shoreline Length", "—",
    "Lake polygon perimeter (m)",
    "Total shoreline length including islands. Longer shorelines increase the interface between terrestrial and aquatic zones, enabling greater nutrient exchange.",
    "Hutchinson (1957)",

    "Lake", "Shoreline Development", "SDI",
    "SDI = Pₗₐₖₑ / (2√(π × Aₗₐₖₑ))",
    "Ratio of actual shoreline to the circumference of a circle with equal area. SDI = 1 for a perfect circle; higher values indicate more irregular shorelines with larger potential littoral zones.",
    "Hakanson (1981)",

    "Lake", "Lake Sphericity", "dₛ/Dₙ",
    "dₛ / Dₙ: diameter of equal-area circle / MRR diagonal",
    "Range 0–1. Values near 1 indicate nearly circular lakes. Lower values indicate more elongated or irregular lake shapes. Dₙ approximated as the diagonal of the minimum rotated rectangle.",
    "Hakanson (1981)",

    "Lake", "Lake Compactness", "—",
    "4π × Aₗₐₖₑ / Pₗₐₖₑ²",
    "Isoperimetric ratio. Range 0–1; equals 1 for a perfect circle, decreases with elongation or shoreline irregularity. Complement to SDI².",
    "—",

    "Lake", "Lake Area as % of Catchment", "—",
    "Lake area (ha) / Catchment area (ha) × 100",
    "Smaller values indicate a larger terrestrial catchment relative to the lake, implying greater nutrient runoff potential and shorter water residence times.",
    "—",

    "Lake", "Catchment-to-Lake Area Ratio", "DA ratio",
    "Catchment area (ha) / Lake area (ha)",
    "Drainage-area ratio. Larger ratios (≥20) increase surface water inflows carrying nutrients and reduce water residence time. Smaller ratios indicate more groundwater-dominated systems.",
    "—",

    # -- Stream -----------------------------------------------------------------
    "Stream", "Drainage Density", "Dd",
    "Dd = total stream length (km) / catchment area (km²)",
    "Length of streams per unit area. Higher values indicate better-drained, less permeable terrain with faster runoff and greater sediment/nutrient transport. Also used to compute the true ruggedness number.",
    "Horton (1932)",

    "Stream", "Stream Frequency", "Fu",
    "Fu = number of stream reaches / catchment area (km²)¹",
    "Number of distinct stream segments per unit area. Higher values indicate finer drainage texture and more active erosion. ¹Requires vector stream input; NA for raster-derived streams.",
    "Horton (1932)"
  )

  ref |>
    gt::gt(groupname_col = "category") |>
    gt::tab_header(
      title    = "Catchment Morphometric Metrics",
      subtitle = "Geometry, areal, relief, lake, and stream parameters"
    ) |>
    gt::cols_label(
      parameter      = "Parameter",
      symbol         = "Symbol",
      formula        = "Formula / Method",
      interpretation = "Interpretation",
      source         = "Source"
    ) |>
    gt::cols_width(
      parameter      ~ gt::px(170),
      symbol         ~ gt::px(55),
      formula        ~ gt::px(200),
      interpretation ~ gt::px(340),
      source         ~ gt::px(130)
    ) |>
    gt::tab_style(
      style     = gt::cell_text(weight = "bold"),
      locations = gt::cells_row_groups()
    ) |>
    gt::tab_style(
      style     = gt::cell_text(
        font  = gt::google_font("Source Code Pro"),
        size  = "small"
      ),
      locations = gt::cells_body(columns = formula)
    ) |>
    gt::tab_style(
      style     = gt::cell_text(weight = "bold"),
      locations = gt::cells_body(columns = symbol)
    ) |>
    gt::tab_footnote(
      footnote = gt::md(
        "**References:** Shekar & Mathew (2024). *Watershed Ecology and the Environment*, 6, 13–25. https://doi.org/10.1016/j.wsee.2023.12.001; Hutchinson, G.E. (1957). *A Treatise on Limnology*, Vol. 1. Wiley; Hakanson, L. (1981). *A Manual of Lake Morphometry*. Springer; Horton, R.E. (1932). *Trans. Am. Geophys. Union*, 13, 350–361."
      )
    ) |>
    gt::opt_row_striping() |>
    gt::opt_table_font(font = gt::google_font("Source Sans Pro")) |>
    gt::tab_options(
      heading.background.color        = "#2C3E50",
      heading.title.font.size         = gt::px(16),
      heading.subtitle.font.size      = gt::px(12),
      row_group.background.color      = "#ECF0F1",
      row_group.font.weight           = "bold",
      column_labels.background.color  = "#34495E",
      column_labels.font.weight       = "bold",
      table.border.top.color          = "#2C3E50",
      table.border.bottom.color       = "#2C3E50"
    ) |>
    gt::tab_style(
      style     = gt::cell_text(color = "white"),
      locations = gt::cells_column_labels()
    )
}

# -- Output writers ------------------------------------------------------------

#' Write metrics results to a CSV and the reference table to an HTML file
#'
#' @param metrics    Tibble from calculate_catchment_metrics()
#' @param ref_table  gt table from build_metrics_reference_table()
#' @param output_dir Character. Root output directory
write_metrics_outputs <- function(metrics, ref_table, output_dir) {
  metrics_path <- fs::path(output_dir, "catchment_metrics.csv")
  readr::write_csv(metrics, metrics_path)
  cw_inform(glue::glue("Metrics written to: {metrics_path}"))

  ref_path <- fs::path(output_dir, "metrics_reference_table.html")
  gt::gtsave(ref_table, ref_path)
  cw_inform(glue::glue("Reference table written to: {ref_path}"))

  invisible(list(metrics = metrics_path, ref_table = ref_path))
}
