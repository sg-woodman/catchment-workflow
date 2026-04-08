# 06_catchment_metrics.R
# ---------------------------------------------------------------------------
# Calculates catchment morphometric metrics for all delineated sites and
# compiles results into a summary tibble and a gt reference table.
#
# Metrics are drawn from:
#   Shekar, P.R. & Mathew, A. (2024). Morphometric analysis of watersheds:
#   A comprehensive review of data sources, quality, and geospatial techniques.
#   Watershed Ecology and the Environment, 6, 13-25.
#   https://doi.org/10.1016/j.wsee.2023.12.001
#
# Three parameter categories are computed (stream/linear metrics excluded):
#   Geometry  : area, perimeter, basin length, mean basin width
#   Areal     : form factor, circularity ratio, elongation ratio,
#               compactness coefficient, shape index, lemniscate ratio
#   Relief    : basin relief, relief ratio, relative relief, ruggedness number,
#               dissection index, gradient ratio, hypsometric integral
#               (hypsometric integral requires the DEM)
#
# Inputs:
#   sites       : Validated sites tibble from validate_sites()
#   output_dir  : Root output directory containing per-site folders
#
# Outputs:
#   metrics tibble : One row per site, one column per metric
#   gt table       : Reference table of metric definitions, formulas,
#                    interpretations, and source
#
# Dependencies: sf, terra, dplyr, purrr, tibble, gt, fs, cli (via utils.R)
# ---------------------------------------------------------------------------

#' Calculate morphometric metrics for all delineated sites
#'
#' Iterates over all sites, loading catchment.gpkg and the site DEM, and
#' computing all available metrics. Sites with missing catchment or DEM
#' outputs are skipped with a warning.
#'
#' @param sites      Validated sites tibble from validate_sites()
#' @param output_dir Character. Root output directory
#' @return A tibble with one row per site containing all computed metrics
calculate_catchment_metrics <- function(sites, output_dir) {

  cw_inform(glue::glue(
    "Calculating morphometric metrics for {nrow(sites)} site(s)..."
  ))

  results <- purrr::map(seq_len(nrow(sites)), function(i) {

    site     <- sites[i, ]
    sid      <- site$site_id
    site_dir <- site_output_dir(output_dir, sid)

    catchment_path <- fs::path(site_dir, "catchment.gpkg")
    dem_path       <- fs::path(site_dir, "dem.tif")

    # Check required inputs exist
    if (!cache_exists(catchment_path)) {
      cw_warn(glue::glue("Site '{sid}': catchment.gpkg not found, skipping."))
      return(NULL)
    }

    cw_inform(glue::glue("Site '{sid}': computing metrics..."))

    catchment <- sf::st_read(catchment_path, quiet = TRUE)

    # Load DEM if available (needed for relief metrics)
    dem <- if (cache_exists(dem_path)) {
      terra::rast(dem_path)
    } else {
      cw_warn(glue::glue(
        "Site '{sid}': dem.tif not found — relief metrics will be NA."
      ))
      NULL
    }

    # Compute all metric categories
    geom    <- compute_geometry_metrics(catchment, sid)
    areal   <- compute_areal_metrics(geom)
    relief  <- compute_relief_metrics(dem, geom, sid)

    # Combine into one row
    tibble::tibble(site_id = sid) |>
      dplyr::bind_cols(geom) |>
      dplyr::bind_cols(areal) |>
      dplyr::bind_cols(relief)

  }) |>
    purrr::compact() |>
    dplyr::bind_rows()

  cw_inform(glue::glue(
    "Metrics computed for {nrow(results)} site(s)."
  ))

  results
}

# -- Geometry metrics --------------------------------------------------------

#' Compute basic geometry metrics from catchment polygon
#'
#' All measurements in metres and km². Catchment is reprojected to EPSG:3979
#' if not already in that CRS to ensure metric units.
#'
#' @param catchment sf polygon. Catchment boundary
#' @param site_id   Character. Site identifier (for log messages)
#' @return Named tibble row with geometry metrics
compute_geometry_metrics <- function(catchment, site_id) {

  # Ensure metric CRS
  catchment <- sf::st_transform(catchment, 3979)

  # Area (km²) and perimeter (km)
  area_m2   <- as.numeric(sf::st_area(catchment))
  area_km2  <- area_m2 / 1e6
  perim_m   <- as.numeric(sf::st_length(sf::st_cast(
    sf::st_union(catchment), "MULTILINESTRING"
  )))
  perim_km  <- perim_m / 1000

  # Basin length (Lb): longest dimension of the catchment bounding box
  # (approximated as the diagonal of the minimum bounding rectangle)
  bbox      <- sf::st_bbox(catchment)
  width_m   <- bbox["xmax"] - bbox["xmin"]
  height_m  <- bbox["ymax"] - bbox["ymin"]
  lb_m      <- sqrt(width_m^2 + height_m^2)
  lb_km     <- lb_m / 1000

  # Mean basin width (Wb = A / Lb)
  wb_km <- area_km2 / lb_km

  tibble::tibble(
    area_km2  = round(area_km2,  4),
    perim_km  = round(perim_km,  4),
    lb_km     = round(lb_km,     4),
    wb_km     = round(wb_km,     4)
  )
}

# -- Areal metrics -----------------------------------------------------------

#' Compute areal (shape) morphometric metrics
#'
#' All metrics derived from geometry metrics only — no DEM required.
#' Formulas follow Shekar & Mathew (2024).
#'
#' @param geom Named tibble row from compute_geometry_metrics()
#' @return Named tibble row with areal metrics
compute_areal_metrics <- function(geom) {

  A  <- geom$area_km2
  P  <- geom$perim_km
  Lb <- geom$lb_km

  # Form factor (Rf) — Horton (1932)
  # Rf = A / Lb²
  # Rf close to 1 → circular basin (high peak flow)
  # Rf << 1 → elongated basin (lower, more sustained flow)
  Rf <- A / Lb^2

  # Circularity ratio (Rc) — Miller (1953)
  # Rc = 4π × A / P²
  # Range 0–1; values approaching 1 indicate circular basin
  Rc <- (4 * pi * A) / P^2

  # Elongation ratio (Re) — Schumm (1956)
  # Re = (2 / Lb) × √(A / π)
  # Range 0–1; <0.5 = highly elongated, 0.5–0.7 = elongated,
  # 0.7–0.8 = less elongated, 0.8–0.9 = oval, 0.9–1.0 = circular
  Re <- (2 / Lb) * sqrt(A / pi)

  # Compactness coefficient (Cc) — Gravelius (1914)
  # Cc = P / (2√(πA))
  # Cc = 1 for a perfect circle; increases with elongation
  Cc <- P / (2 * sqrt(pi * A))

  # Shape index (Sw) — inverse of form factor (Horton 1932)
  # Sw = Lb² / A
  # Higher value → more elongated
  Sw <- Lb^2 / A

  # Lemniscate ratio (k) — Chorley et al. (1957)
  # k = Lb² / (4A)
  # k approaching 1 → circular basin; k > 1 → elongated
  k <- Lb^2 / (4 * A)

  tibble::tibble(
    form_factor          = round(Rf, 4),
    circularity_ratio    = round(Rc, 4),
    elongation_ratio     = round(Re, 4),
    compactness_coeff    = round(Cc, 4),
    shape_index          = round(Sw, 4),
    lemniscate_ratio     = round(k,  4)
  )
}

# -- Relief metrics ----------------------------------------------------------

#' Compute relief morphometric metrics from the catchment DEM
#'
#' Requires the site-level DEM (clipped to catchment). Returns NA for all
#' metrics if DEM is not available.
#'
#' @param dem     SpatRaster or NULL. Site DEM clipped to catchment
#' @param geom    Named tibble row from compute_geometry_metrics()
#' @param site_id Character. Site identifier (for log messages)
#' @return Named tibble row with relief metrics
compute_relief_metrics <- function(dem, geom, site_id) {

  # Return all-NA if DEM unavailable
  if (is.null(dem)) {
    return(tibble::tibble(
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
    ))
  }

  vals     <- terra::values(dem, na.rm = TRUE)
  elev_min <- min(vals)
  elev_max <- max(vals)

  # Basin relief (H) — difference between highest and lowest elevation (m)
  H <- elev_max - elev_min

  # Relief ratio (Rh) — Schumm (1956)
  # Rh = H / Lb (H in km, Lb in km)
  # Indicates overall steepness and intensity of erosion
  Lb_m  <- geom$lb_km * 1000
  Rh    <- (H / 1000) / geom$lb_km

  # Relative relief (Rhp) — Melton (1957)
  # Rhp = H × 100 / P (H in m, P in km → multiply P by 1000)
  Rhp <- (H * 100) / (geom$perim_km * 1000)

  # Ruggedness number (Rn) — Patton & Baker (1976)
  # Rn = H × Dd  where Dd = drainage density (km/km²)
  # Without stream data we approximate Dd from basin geometry:
  # Dd ≈ Lb / A (a simplified proxy — flag in output)
  # Full Dd requires stream network data
  Dd_proxy <- geom$lb_km / geom$area_km2
  Rn       <- (H / 1000) * Dd_proxy

  # Dissection index (Dis) — Singh & Dubey (1994)
  # Dis = H / elev_max (both in same units)
  # Range 0–1; higher values indicate more dissected terrain
  Dis <- H / elev_max

  # Hypsometric integral (HI) — Pike & Wilson (1971) elevation-relief ratio method
  # HI = (mean_elev - elev_min) / (elev_max - elev_min)
  # HI > 0.6 → monadnock phase (youthful, high erosion potential)
  # HI 0.35–0.6 → equilibrium (mature)
  # HI < 0.35 → peneplain phase (old, low erosion)
  mean_elev <- mean(vals)
  HI        <- (mean_elev - elev_min) / (elev_max - elev_min)

  # -- Slope (degrees) -------------------------------------------------------
  # Mean slope across the catchment in degrees, derived from the DEM using
  # a 3x3 neighbourhood gradient. Higher values indicate steeper terrain.
  slope_rast      <- terra::terrain(dem, v = "slope", unit = "degrees")
  slope_vals      <- terra::values(slope_rast, na.rm = TRUE)
  mean_slope_deg  <- mean(slope_vals)

  # -- Aspect (degrees, circular mean) ----------------------------------------
  # Catchment mean aspect using circular statistics to handle the 0/360°
  # discontinuity. Aspect is decomposed into sin/cos unit vectors, averaged,
  # then converted back to degrees. Also returns a dominant aspect class.
  aspect_rast  <- terra::terrain(dem, v = "aspect", unit = "degrees")
  aspect_vals  <- terra::values(aspect_rast, na.rm = TRUE)

  # Circular mean: average unit vectors then convert back
  mean_sin     <- mean(sin(aspect_vals * pi / 180))
  mean_cos     <- mean(cos(aspect_vals * pi / 180))
  mean_aspect  <- (atan2(mean_sin, mean_cos) * 180 / pi) %% 360

  # Dominant aspect class (8-point compass)
  aspect_class <- dplyr::case_when(
    mean_aspect >= 337.5 | mean_aspect < 22.5   ~ "N",
    mean_aspect >= 22.5  & mean_aspect < 67.5   ~ "NE",
    mean_aspect >= 67.5  & mean_aspect < 112.5  ~ "E",
    mean_aspect >= 112.5 & mean_aspect < 157.5  ~ "SE",
    mean_aspect >= 157.5 & mean_aspect < 202.5  ~ "S",
    mean_aspect >= 202.5 & mean_aspect < 247.5  ~ "SW",
    mean_aspect >= 247.5 & mean_aspect < 292.5  ~ "W",
    mean_aspect >= 292.5 & mean_aspect < 337.5  ~ "NW",
    TRUE ~ NA_character_
  )

  # -- Elevation standard deviation -------------------------------------------
  # Variability of elevation across the catchment. Higher values indicate
  # more heterogeneous, rugged terrain.
  elev_sd <- stats::sd(vals)

  tibble::tibble(
    elev_min_m           = round(elev_min,       2),
    elev_max_m           = round(elev_max,       2),
    elev_sd_m            = round(elev_sd,        2),
    basin_relief_m       = round(H,              2),
    relief_ratio         = round(Rh,             4),
    relative_relief      = round(Rhp,            4),
    ruggedness_number    = round(Rn,             4),
    dissection_index     = round(Dis,            4),
    hypsometric_integral = round(HI,             4),
    mean_slope_deg       = round(mean_slope_deg, 2),
    mean_aspect_deg      = round(mean_aspect,    1),
    aspect_class         = aspect_class
  )
}

# -- Reference table ---------------------------------------------------------

#' Build a gt reference table of all morphometric metrics
#'
#' Returns a formatted gt table with one row per metric showing the parameter
#' name, symbol, formula, interpretation, and source citation.
#'
#' @return A gt table object
build_metrics_reference_table <- function() {

  ref <- tibble::tribble(

    ~category, ~parameter, ~symbol, ~formula, ~interpretation, ~source,

    # -- Geometry -------------------------------------------------------------
    "Geometry", "Catchment Area",    "A",   "Planimetric area of catchment polygon",
    "Total area draining to the pour point. Larger catchments generally produce higher total discharge but lower peak runoff per unit area.",
    "—",

    "Geometry", "Perimeter",         "P",   "Length of catchment boundary",
    "Longer perimeters relative to area indicate more irregular, elongated shapes.",
    "—",

    "Geometry", "Basin Length",      "Lb",  "Diagonal of minimum bounding rectangle (km)",
    "Approximates the longest axis of the basin. Used in shape and relief ratio calculations.",
    "Schumm (1956)",

    "Geometry", "Mean Basin Width",  "Wb",  "Wb = A / Lb",
    "Average width perpendicular to the basin length. Narrow basins produce more peaked hydrographs.",
    "Horton (1932)",

    # -- Areal ----------------------------------------------------------------
    "Areal", "Form Factor",          "Rf",  "Rf = A / Lb\u00b2",
    "Values near 1 indicate a circular basin with high, peaked flood discharge. Values < 0.45 indicate elongated basins with lower, flatter hydrographs and lower erosion potential.",
    "Horton (1932)",

    "Areal", "Circularity Ratio",    "Rc",  "Rc = 4\u03c0A / P\u00b2",
    "Range 0\u20131. Values approaching 1 indicate a circular, compact basin. Low values (<0.4) indicate elongated basins with low discharge and high permeability. Values 0.4\u20130.5 = strongly elongated; 0.5\u20130.6 = elongated; 0.6\u20130.7 = slightly elongated; 0.7\u20130.8 = oval; 0.8\u20131.0 = circular.",
    "Miller (1953)",

    "Areal", "Elongation Ratio",     "Re",  "Re = (2/Lb) \u00d7 \u221a(A/\u03c0)",
    "Range 0\u20131. Values near 1 indicate circular basins. <0.5 = highly elongated; 0.5\u20130.7 = elongated; 0.7\u20130.8 = less elongated; 0.8\u20130.9 = oval; 0.9\u20131.0 = circular. Circular basins have higher runoff and erosion susceptibility.",
    "Schumm (1956)",

    "Areal", "Compactness Coefficient", "Cc", "Cc = P / (2\u221a(\u03c0A))",
    "Equals 1 for a perfect circle; increases with elongation. Higher values indicate more elongated basins with lower flood susceptibility.",
    "Gravelius (1914)",

    "Areal", "Shape Index",          "Sw",  "Sw = Lb\u00b2 / A",
    "Inverse of form factor. Higher values indicate more elongated basins. Values < 2 = circular; 2\u20134 = slightly elongated; > 4 = strongly elongated.",
    "Horton (1932)",

    "Areal", "Lemniscate Ratio",     "k",   "k = Lb\u00b2 / (4A)",
    "Values near 1 indicate circular basins. Values > 1 indicate elongated basins with lower peak discharge and reduced erosion susceptibility.",
    "Chorley et al. (1957)",

    # -- Relief ---------------------------------------------------------------
    "Relief", "Basin Relief",        "H",   "H = Elev\u2098\u2090\u2093 - Elev\u2098\u1d35\u2099 (m)",
    "Vertical distance between highest and lowest points. High relief indicates steep terrain, greater erosive energy, and higher peak flows.",
    "Hadley & Schumm (1961)",

    "Relief", "Relief Ratio",        "Rh",  "Rh = H (km) / Lb (km)",
    "Measures overall steepness of the watershed. Higher values indicate steeper terrain and greater intensity of erosion processes. Low: <0.1; Medium: 0.1\u20130.3; High: >0.3.",
    "Schumm (1956)",

    "Relief", "Relative Relief",     "Rhp", "Rhp = H \u00d7 100 / P (m/km)",
    "Normalises relief by perimeter length. Higher values indicate more rugged, steeply incised terrain.",
    "Melton (1957)",

    "Relief", "Ruggedness Number",   "Rn",  "Rn = H (km) \u00d7 Dd (km/km\u00b2)\u00b9",
    "Combines relief and drainage density to indicate the structural complexity and erosion potential of the basin. High values indicate steep, well-dissected terrain prone to high runoff and erosion. \u00b9Dd approximated as Lb/A without stream network data.",
    "Patton & Baker (1976)",

    "Relief", "Dissection Index",    "Dis", "Dis = H / Elev\u2098\u2090\u2093",
    "Range 0\u20131. Measures the degree of terrain dissection relative to maximum elevation. Higher values indicate more deeply incised, dissected terrain.",
    "Singh & Dubey (1994)",

    "Relief", "Hypsometric Integral", "HI", "HI = (Elev\u2098\u2091\u2090\u2099 - Elev\u2098\u1d35\u2099) / (Elev\u2098\u2090\u2093 - Elev\u2098\u1d35\u2099)",
    "Range 0\u20131. HI > 0.6 = monadnock/youthful stage (high erosion potential, most mass still uneroded); 0.35\u20130.6 = equilibrium/mature stage; < 0.35 = peneplain/old stage (highly eroded, low relief).",
    "Pike & Wilson (1971)"
  )

  ref |>
    gt::gt(groupname_col = "category") |>
    gt::tab_header(
      title    = "Catchment Morphometric Metrics",
      subtitle = "Areal and relief parameters following Shekar & Mathew (2024)"
    ) |>
    gt::cols_label(
      parameter      = "Parameter",
      symbol         = "Symbol",
      formula        = "Formula",
      interpretation = "Interpretation",
      source         = "Source"
    ) |>
    gt::cols_width(
      parameter      ~ gt::px(160),
      symbol         ~ gt::px(50),
      formula        ~ gt::px(180),
      interpretation ~ gt::px(340),
      source         ~ gt::px(130)
    ) |>
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_row_groups()
    ) |>
    gt::tab_style(
      style = gt::cell_text(font = gt::google_font("Source Code Pro"), size = "small"),
      locations = gt::cells_body(columns = formula)
    ) |>
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_body(columns = symbol)
    ) |>
    gt::tab_footnote(
      footnote = gt::md(
        "**Reference:** Shekar, P.R. & Mathew, A. (2024). Morphometric analysis of watersheds: A comprehensive review of data sources, quality, and geospatial techniques. *Watershed Ecology and the Environment*, 6, 13\u201325. https://doi.org/10.1016/j.wsee.2023.12.001"
      )
    ) |>
    gt::opt_row_striping() |>
    gt::opt_table_font(font = gt::google_font("Source Sans Pro")) |>
    gt::tab_options(
      heading.background.color       = "#2C3E50",
      heading.title.font.size        = gt::px(16),
      heading.subtitle.font.size     = gt::px(12),
      row_group.background.color     = "#ECF0F1",
      row_group.font.weight          = "bold",
      column_labels.background.color = "#34495E",
      column_labels.font.weight      = "bold",
      table.border.top.color         = "#2C3E50",
      table.border.bottom.color      = "#2C3E50"
    ) |>
    gt::tab_style(
      style     = gt::cell_text(color = "white"),
      locations = gt::cells_column_labels()
    )
}

#' Write metrics results to a CSV and the reference table to an HTML file
#'
#' @param metrics    Tibble from calculate_catchment_metrics()
#' @param ref_table  gt table from build_metrics_reference_table()
#' @param output_dir Character. Root output directory
write_metrics_outputs <- function(metrics, ref_table, output_dir) {

  # Write metrics CSV
  metrics_path <- fs::path(output_dir, "catchment_metrics.csv")
  readr::write_csv(metrics, metrics_path)
  cw_inform(glue::glue("Metrics written to: {metrics_path}"))

  # Write reference table as HTML
  ref_path <- fs::path(output_dir, "metrics_reference_table.html")
  gt::gtsave(ref_table, ref_path)
  cw_inform(glue::glue("Reference table written to: {ref_path}"))

  invisible(list(metrics = metrics_path, ref_table = ref_path))
}
