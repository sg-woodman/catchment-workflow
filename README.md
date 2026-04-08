# Catchment Delineation Workflow
 
A modular, reproducible R workflow for delineating catchments and computing morphometric metrics across multiple sites in Canada. Designed for project-agnostic use — bring a sites table, point it at your data, and run.
 
**Author:** Sam W.  
Natural Resources Canada
 
> Caching logic, data download helpers, and directory structure components were developed with assistance from Claude (Anthropic).
 
---
 
## Overview
 
The workflow delineates catchments from pour points using the [MRDEM](https://open.canada.ca/data/en/dataset/18752265-bda3-498c-a4ba-9dfe68cb98da) (30 m resolution) and [WhiteboxTools](https://www.whiteboxgeo.com/), with optional stream burning from the [National Hydro Network (NHN)](https://www.nrcan.gc.ca/science-and-data/science-and-research/earth-sciences/geography/topographic-information/geobase-surface-water-program-geeau/national-hydrographic-network/21361). Sites are grouped using [HydroBasins](https://www.hydrosheds.org/products/hydrobasins) so that expensive DEM conditioning steps run once per hydrologically coherent region rather than once per site.
 
For each site the workflow produces a self-contained output folder containing the catchment polygon, snapped pour point, all hydrological raster layers clipped to the catchment, NHN flowlines, and hillshade. An optional final step computes catchment morphometric metrics following Shekar & Mathew (2024).
 
---
 
## Requirements
 
### R packages
 
```r
install.packages(c(
  "sf", "terra", "whitebox",
  "dplyr", "tidyr", "purrr", "readr", "tibble",
  "fs", "cli", "glue", "gt"
))
```
 
### WhiteboxTools binary
 
```r
whitebox::install_whitebox()
```
 
### External data (user-supplied)
 
| Dataset | Source | Notes |
|---|---|---|
| MRDEM DTM `.vrt` | [NRCan Open Government](https://open.canada.ca/data/en/dataset/18752265-bda3-498c-a4ba-9dfe68cb98da) | Stream, not download |
| NHN GDB files | [GeoBase NHN](https://open.canada.ca/data/en/dataset/a4b190fe-e090-4e6d-881e-b87956c07977) | Per NTS sheet, pre-downloaded |
| NHN index shapefile | Same source | `NHN_INDEX_WORKUNIT_LIMIT_2.shp` |
| HydroBasins | [HydroSHEDS](https://www.hydrosheds.org/products/hydrobasins) | `na` and `ar` level 1–12 `.shp` files |
 
---
 
## Directory Structure
 
```
catchment-delineation/
  R/
    utils.R                  # shared helpers, validation, logging
    01_group_sites.R         # group AOI construction from HydroBasins
    02_prepare_dem.R         # MRDEM crop to group AOI
    03_burn_streams.R        # NHN stream burn-in
    04_run_whitebox.R        # breach, pointer, accumulation, streams, hillshade
    05_delineate_sites.R     # per-site watershed delineation
    06_catchment_metrics.R   # morphometric metrics
  sites_template.R           # define sites as a tibble (recommended)
  sites_template.csv         # alternative CSV input
  delineate_catchments.R     # top-level runner
  reset_workflow.R           # cache/output reset helpers
  cache/                     # group-level rasters (auto-created)
  output/                    # per-site outputs (auto-created)
```
 
---
 
## Quick Start
 
**1. Define your sites** in `sites_template.R`:
 
```r
sites <- tribble(
  ~site_id,  ~site_name,   ~lon,      ~lat,    ~group_id,   ~burn_streams, ~aoi_buffer_m, ~stream_threshold,
  "L03",     "Lake 03",    -81.234,   46.123,  "sudbury",   TRUE,          NA,            NA,
  "L04",     "Lake 04",    -81.198,   46.145,  "sudbury",   TRUE,          NA,            NA,
  "YK01",    "Yukon 01",   -114.567,  62.456,  "yukon",     FALSE,         25000,         2000
)
```
 
**2. Set paths** in `delineate_catchments.R`:
 
```r
mrdem_vrt       <- "/path/to/mrdem-30-dtm.vrt"
nhn_dir         <- "/path/to/NHN/gdb"
nhn_index       <- "/path/to/NHN_INDEX_WORKUNIT_LIMIT_2.shp"
hydrobasins_dir <- "/path/to/HydroBasins"
```
 
**3. Run**:
 
```r
source("delineate_catchments.R")
```
 
---
 
## Site Table Columns
 
| Column | Required | Description |
|---|---|---|
| `site_id` | ✓ | Unique identifier, used as output folder name |
| `site_name` | ✓ | Human-readable name |
| `lon` | ✓ | Pour point longitude (WGS84 decimal degrees, negative for Canada) |
| `lat` | ✓ | Pour point latitude (WGS84 decimal degrees) |
| `group_id` | ✓ | Processing group — sites sharing a group share a DEM and flow products |
| `burn_streams` | ✓ | `TRUE`/`FALSE` — whether to burn NHN streams into DEM (group-level) |
| `aoi_buffer_m` | | Buffer in metres added to HydroBasins polygon. Default 1000 m |
| `stream_threshold` | | Flow accumulation threshold for stream extraction. Default 1000 cells |
 
---
 
## Per-Site Outputs
 
Each site produces a folder `output/<site_id>/` containing:
 
| File | Description |
|---|---|
| `catchment.gpkg` | Catchment polygon (EPSG:3979) |
| `pour_point.gpkg` | Snapped pour point (EPSG:3979) |
| `pour_point.shp` | Original pour point from coordinates |
| `pour_point_snapped.shp` | Snapped pour point — **edit this to correct misaligned pour points** |
| `streams_tmp.tif` | Streams raster used for snapping |
| `flow_pointer_tmp.tif` | Flow pointer used for watershed delineation |
| `watershed.tif` | Binary watershed raster |
| `dem.tif` | DEM clipped to catchment |
| `dem_breached.tif` | Breached DEM clipped to catchment |
| `flow_pointer.tif` | D8 flow pointer clipped to catchment |
| `flow_accum.tif` | Flow accumulation clipped to catchment |
| `streams.tif` | Stream network raster clipped to catchment |
| `hillshade.tif` | Hillshade clipped to catchment |
| `streams.gpkg` | NHN flowlines clipped to catchment |
 
Group-level cache files live in `cache/<group_id>/` and are shared across all sites in the group.
 
---
 
## Correcting a Pour Point
 
If a catchment looks wrong, the pour point likely snapped to the wrong stream cell. To fix:
 
1. Open `pour_point_snapped.shp` and `streams_tmp.tif` in QGIS for the affected site
2. Edit `pour_point_snapped.shp` — move the point to the correct stream cell
3. In R:
```r
source("reset_workflow.R")
reset_site("site_id")   # clears only that site's output folder
source("delineate_catchments.R")  # re-runs from Stage 5 only
```
 
---
 
## Resetting the Workflow
 
```r
source("reset_workflow.R")
 
reset_all()              # delete everything in cache/ and output/
reset_cache()            # delete group rasters only
reset_sites(sites)       # delete all site output folders
reset_group("sudbury", sites)  # delete one group cache + its sites
reset_site("L03")        # delete one site output folder
```
 
---
 
## Morphometric Metrics
 
After delineation, compute catchment metrics:
 
```r
source("R/06_catchment_metrics.R")
 
metrics   <- calculate_catchment_metrics(sites, output_dir = "output")
ref_table <- build_metrics_reference_table()
write_metrics_outputs(metrics, ref_table, output_dir = "output")
```
 
Outputs: `output/catchment_metrics.csv` and `output/metrics_reference_table.html`.
 
Metrics follow: Shekar, P.R. & Mathew, A. (2024). Morphometric analysis of watersheds: A comprehensive review of data sources, quality, and geospatial techniques. *Watershed Ecology and the Environment*, 6, 13–25.
 
---
 
## Key Dependencies and Data Sources
 
- **MRDEM**: Canada's 30 m national DEM in EPSG:3979, replacing the CDEM. Streamed via VRT — no tile downloads required.
- **NHN**: National Hydro Network flowlines and waterbodies used for stream burning and lake extraction. Distributed per NTS sheet as GDB files.
- **HydroBasins**: HydroSHEDS watershed polygons used to define processing group AOIs. Both North America (`na`) and Arctic (`ar`) datasets required for Canadian sites.
- **WhiteboxTools**: Open-source geospatial analysis library used for DEM conditioning, flow routing, and watershed delineation.
 
---
 
## License
 
For internal use. Contact Sam W. (Natural Resources Canada) for inquiries.
 