# =============================================================================
# config.R — Shared paths and constants
# =============================================================================
# Edit this file when moving the project to a new machine or data directory.
# No paths or magic numbers should appear anywhere else in the codebase.

# -----------------------------------------------------------------------------
# Shared data paths (downloaded once, used across projects)
# -----------------------------------------------------------------------------

# MRDEM Virtual Raster — seamless 30 m DEM for all of Canada
MRDEM_PATH <- "/Users/sam/Documents/cfs/shared_data/raw/dem/mrdem-30-dtm.vrt"

# HydroBasins regional directories (hybas_na_lev0N_v1c.shp, hybas_ar_lev0N_v1c.shp)
# Download from: https://www.hydrosheds.org/products/hydrobasins
HYDROBASINS_DIRS <- list(
  na = "/Users/sam/Documents/cfs/shared_data/raw/hydro/watersheds/HydroBasins/north_america",
  ar = "/Users/sam/Documents/cfs/shared_data/raw/hydro/watersheds/HydroBasins/arctic"
)

# NHN work unit index shapefile — download once from:
#   https://ftp.maps.canada.ca/pub/nrcan_rncan/vector/geobase_nhn_rhn/index/
#   NHN_INDEX_WORKUNIT_LIMIT_2.zip
NHN_INDEX_PATH <- "/Users/sam/Documents/cfs/shared_data/raw/hydro/networks/NHN/NHN_INDEX_WORKUNIT_LIMIT_2/NHN_INDEX_22_INDEX_WORKUNIT_LIMIT_2.shp"

# Directory where NHN GDB zips are extracted permanently (shared across projects)
NHN_RAW_DIR <- "/Users/sam/Documents/cfs/shared_data/raw/hydro/networks/NHN/gdb"

# -----------------------------------------------------------------------------
# Project output paths
# -----------------------------------------------------------------------------

# All derived products — targets manages this directory
CACHE_DIR <- here::here("data/processed/cache")

# Final catchment polygons
OUTPUT_DIR <- here::here("data/processed/catchments")

# -----------------------------------------------------------------------------
# HydroBasins constants
# -----------------------------------------------------------------------------

# Latitude boundary between North America and Arctic regional files
HYBAS_BOUNDARY_LAT  <- 60.0

# Sites within this many degrees of the boundary are checked against both files
HYBAS_BOUNDARY_ZONE <- 1.0

# Pfafstetter levels to try when iteratively expanding the AOI, coarsest last.
# Level 6 covers most catchments; 4 and 3 handle very large drainages.
HYBAS_LEVELS <- c(6L, 5L, 4L, 3L)

# Buffer (metres) added around the HydroBasins polygon to ensure the full
# upstream area is captured even when the pour point sits near a basin edge
HYBAS_BUFFER_M <- 10000L

# -----------------------------------------------------------------------------
# WhiteboxTools constants
# -----------------------------------------------------------------------------

# Stream burn depth (metres) — deep enough to override DEM noise, shallow
# enough to avoid unrealistic incisions
WBT_BURN_DEPTH_M <- 5L

# Maximum breach distance (raster cells) for least-cost depression breaching
WBT_BREACH_DIST <- 10L

# Pour point snap distance (metres) — snaps to nearest stream cell
WBT_SNAP_DIST_M <- 200L

# Margin (metres) inset from AOI edge used to detect boundary-touching catchments.
# ~17 MRDEM pixels at 30 m resolution — large enough to be unambiguous
BOUNDARY_MARGIN_M <- 500L

# -----------------------------------------------------------------------------
# Coordinate reference systems
# -----------------------------------------------------------------------------

# Working CRS for all vector and raster operations
CRS_PROJECTED <- "EPSG:3979"

# Output CRS for final catchment polygons (equal-area for correct area calc)
CRS_OUTPUT <- "EPSG:3347"
