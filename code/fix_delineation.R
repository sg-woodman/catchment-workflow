# fix_delineation.R
# ---------------------------------------------------------------------------
# Script designed to manually correct catchment delineation for sites where the
# pour point was not accurately snapped to the stream layer. Incorrect
# delineation is often the result of a bad stream threshold in
# wbt_extract_streams.
#
# This approach is designed to use after the full delineation workflow. The
# expected inputs the D8 pntr file and a manually adjusted snapped pout point.
# When erroneous catchments (i.e. too small) are inspected in QGIS, the flow
# accumulation and streams raster should be loaded. An appropriate flow
# accumlation pixel (i.e. larger number) should be selected and the snapped
# point is moved to that location. Once edited, this script can be rerun to
# generate a proper catchment boundary.
#
# SETUP: NA
# ---------------------------------------------------------------------------

# Load required packages
library(here)
library(tidyverse)
library(sf)
library(terra)
library(whitebox)

# Set catchment ID and cache ID
site_id <- "PHPP06"
cache_id <- "Sask"

# Redelineate the watershed using the updated pour point
wbt_watershed(
  d8_pntr = paste0(here("cache/"), cache_id, "/flow_pointer.tif"),
  pour_pts = paste0(here("output/"), site_id, "/pour_point_snapped.shp"),
  output = paste0(here("output/"), site_id, "/watershed_corrected.tif")
)

# Convert the corrected watershed to a polygon
watershed <- rast(paste0(here("output/"), site_id, "/watershed_corrected.tif"))
catchment <- as.polygons(watershed)
plot(catchment)

# Save output
writeVector(
  catchment,
  paste0(here("output/"), site_id, "/catchment.gpkg"),
  overwrite = T
)
