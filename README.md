# Catchment Delineation Workflow

A modular, reproducible R workflow for delineating stream and lake catchments and computing morphometric and distance-weighted landscape metrics across multiple sites in Canada. Project-agnostic: bring a sites table and paths to your terrain/hydrology data, and it produces self-contained per-site output folders.

**Author:** Sam Woodman.
Natural Resources Canada

> Caching logic, the modular delineation engine, and directory structure components were developed with assistance from Claude (Anthropic).

---

## Overview

At the core is a single, input-driven delineation **engine** (`workflow/R/engine/`) rather than one fixed pipeline. It runs off whatever terrain a project actually supplies — a raw DEM needing depression-breaching, or an already-conditioned flow direction/flow pointer raster — and off whatever pour points a project uses — plain point coordinates (Jenson-snapped to a stream network) or lake polygons. Optional NHN stream burn-in and HydroBasins-based regional grouping are both configurable, not hardcoded. See `workflow/R/engine/00_resolve_config.R`'s header for the full config schema.

Three projects currently run on this engine:

| Project | Runner | Pour points | Terrain source | Grouping |
|---|---|---|---|---|
| CELESTE (stream/river sites) | `workflow/CELESTE/run_celeste.R` | Point (Jenson snap) | MRDEM (raw, needs breach) + NHN burn-in | HydroBasins regions |
| CAM streams | `workflow/CAM/run_cam_streams.R` | Point (Jenson snap) | OIH Enhanced Flow Direction (pre-conditioned) | Whole domain |
| CAM lakes | `workflow/run_cam_lakes.R` | Lake polygon | OIH Enhanced Flow Direction (pre-conditioned) | Whole domain |

For each site the workflow produces a self-contained output folder containing the catchment polygon (unclipped and, where nested sites exist, clipped), snapped pour point (or lake pour-point raster), every hydrological raster clipped to the catchment, NHN flowlines where applicable, and hillshade. Standard downstream stages compute catchment morphometric metrics (Shekar & Mathew, 2024) and distance-weighted landscape attributes (the `hydroweight` package — % land cover, mean NDVI, etc. under 7 weighting schemes) for whatever land-cover-type layers of interest (LOIs) the project defines.

`CLAUDE.md` is the full architecture reference — module-by-module purpose tables, the `group_manifest` object, caching/CRS conventions, and a running list of real `terra`/`hydroweight` gotchas found while building this. This README is the shorter orientation; read `CLAUDE.md` for anything not covered here.

---

## Requirements

### R packages

```r
install.packages(c(
  "sf", "terra", "whitebox", "hydroweight",
  "dplyr", "tidyr", "purrr", "readr", "tibble",
  "fs", "cli", "glue", "here", "gt"
))
```

### WhiteboxTools binary

```r
whitebox::install_whitebox()
```

### External data (user-supplied)

| Dataset | Used by | Notes |
|---|---|---|
| MRDEM DTM `.vrt` | CELESTE | [NRCan Open Government](https://open.canada.ca/data/en/dataset/18752265-bda3-498c-a4ba-9dfe68cb98da) — streamed via VRT, no tile downloads |
| NHN GDB files | CELESTE | [GeoBase NHN](https://open.canada.ca/data/en/dataset/a4b190fe-e090-4e6d-881e-b87956c07977), per NTS sheet (`nhn_rhn_<sheet>_gdb_en/`) |
| NHN index shapefile | CELESTE | Same source — `NHN_INDEX_WORKUNIT_LIMIT_2.shp`, maps sheet codes to geometries |
| HydroBasins | CELESTE | [HydroSHEDS](https://www.hydrosheds.org/products/hydrobasins) — `north_america/` and `arctic/` subdirectories, levels 1–12 |
| OIH Enforced DEM + Enhanced Flow Direction | CAM streams, CAM lakes | Ontario Integrated Hydrology — pre-conditioned, no breach needed |
| OIH waterbody polygons | CAM lakes | Used to match sites to lake polygons and as the hydroweight `target_O` |
| CANLCC 2020 land cover | CAM lakes (hydroweight LOI) | Reprojected once and cached on first use |
| Harvest/regen source data (Ontario MNR gdb, NB/Irving gdb) | CELESTE, CAM streams (hydroweight LOI) | See `workflow/CELESTE/prepare_harvest_regen.R` / `workflow/CAM/prepare_harvest_regen.R` |
| NDVI tiles / GEE exports | CELESTE, CAM streams (hydroweight LOI) | See the respective `prepare_ndvi.R` |

On Sam's machine these live under `~/Documents/cfs/shared_data/raw/`; exact paths are set at the top of each runner.

---

## Directory Structure

```
workflow/
├── R/
│   ├── utils.R                   shared: logging, validation, AOI helpers
│   ├── stream/                   reusable stream/river building blocks —
│   │   │                         used directly by engine/, not a pipeline
│   │   │                         of their own
│   │   ├── group_sites.R
│   │   ├── burn_streams.R
│   │   ├── delineate_sites.R
│   │   ├── 99_rerun_sites        pre-engine only; kept for an unrelated script
│   │   └── hydroweight_attributes.R
│   ├── lake/                     reusable lake building blocks, same role
│   │   ├── match_lake_polygons.R
│   │   ├── remove_upstream_lakes.R
│   │   └── hydroweight_attributes.R
│   ├── engine/                   the modular, input-driven delineation
│   │   │                         engine — standard approach for every
│   │   │                         project now
│   │   ├── 00_resolve_config.R
│   │   ├── 01_build_group_manifest.R
│   │   ├── 02_prepare_terrain.R
│   │   ├── 03_prepare_streams_burn.R
│   │   ├── 04_delineate_site.R
│   │   └── 99_rerun_sites        rerun_engine_site_watershed(), etc.
│   ├── remove_upstream.R         shared (stream nesting)
│   ├── reclip_outputs.R          shared (stream + lake)
│   └── catchment_metrics.R       shared (stream + lake)
├── raster_attributes.R           generic LOI-prep toolbox (project-agnostic)
├── CELESTE/                      CELESTE runner (engine-based) + LOI prep
│   ├── run_celeste.R
│   ├── prepare_ndvi.R
│   ├── prepare_ndvi_trend.R
│   ├── prepare_harvest_regen.R
│   └── tidy_outputs.R            long-format plotting-ready tables
├── CAM/                          CAM streams runner (engine-based) + LOI prep
│   ├── run_cam_streams.R
│   ├── prepare_ndvi.R
│   ├── prepare_ndvi_trend.R
│   ├── prepare_harvest_regen.R
│   ├── upload_catchments_to_gee.R
│   └── tidy_outputs.R
├── templates/
│   └── run_engine_template.R     blank, TODO-annotated skeleton for a new project
└── run_cam_lakes.R                CAM lakes runner (engine-based)

code/
└── reset_workflow.R              cache/output reset helpers

cache/<group_id>/  or  cache/<PROJECT>/       group-level rasters (auto-created)
output/<site_id>/  or  output/<PROJECT>/…     per-site outputs (auto-created)
```

---

## Running an existing project

Each runner is sourced interactively in R (or run section-by-section), not executed as a CLI script:

```r
# From an R session in the project root:
source("workflow/CELESTE/run_celeste.R")     # or workflow/CAM/run_cam_streams.R, workflow/run_cam_lakes.R
```

Before running, set the path variables near the top of that file (DEM/flow direction source, NHN or OIH data, HydroBasins directory as applicable) and point it at your sites source.

---

## Starting a new project

Copy `workflow/templates/run_engine_template.R` to `workflow/<PROJECT>/run_<project>.R` and fill in the TODOs. The config you build is a plain list, resolved by `resolve_engine_config()`:

```r
run_config <- list(
  project_id = "MY_PROJECT",
  output_dir = here("output/MY_PROJECT"),
  cache_dir  = here("cache/MY_PROJECT"),
  sites      = sites,                    # tibble: site_id, lon, lat, ...

  dem            = list(path = "/path/to/dem.vrt"),  # OR flow_direction / flow_pointer
  flow_direction = NULL,
  flow_pointer   = NULL,
  crs              = NULL,               # NULL = match the terrain's own native CRS
  stream_threshold = 1000,

  streams_burn = list(source = "none"),  # "nhn_auto" | "supplied" | "none"

  lake_polygons = NULL,                  # NULL = point pour points; supply an sf/SpatVector for lake mode
  lake_buffer_m = 30,

  grouping = list(strategy = "whole_domain"),  # or "hydrobasins" (+ hydrobasins_dir)

  loi_layers = NULL                      # set later, for hydroweight
)
```

`workflow/CAM/run_cam_streams.R` is the leanest filled-in example (OIH DEM, point pour points, whole-domain grouping); `workflow/CELESTE/run_celeste.R` shows the HydroBasins + NHN-burn-in case, and `workflow/run_cam_lakes.R` the lake-polygon case.

---

## Site table columns

Exact requirements depend on `grouping$strategy` and pour-point mode, but in general:

| Column | Required | Description |
|---|---|---|
| `site_id` | ✓ | Unique identifier, used as output folder name |
| `site_name` | ✓ | Human-readable name |
| `lon`, `lat` | ✓ (point mode) | Pour point coordinates (WGS84 decimal degrees) |
| `lake_name` | ✓ (lake mode) | Matched against lake polygon source data |
| `group_id` | only for `"hydrobasins"` grouping | Sites sharing a group share a DEM crop and flow products |
| `burn_streams` | only when `streams_burn` applies | `TRUE`/`FALSE`, group-level |
| `aoi_buffer_m` | optional | Buffer (m) added to the HydroBasins AOI polygon. Default 1000 m |
| `stream_threshold` | optional | Flow accumulation threshold for stream extraction. Default from config |

---

## Per-site outputs

Each site produces a folder `output/<PROJECT>/<site_id>/` (point mode) or the lake-mode equivalent, containing:

| File | Description |
|---|---|
| `catchment.gpkg` | Unclipped catchment polygon |
| `catchment_clipped.gpkg` | Catchment with nested upstream sites' area erased (identical to unclipped if none) |
| `pour_point.gpkg` / `pour_point_snapped.shp` | Snapped pour point (point mode) — **edit the `.shp` to correct a misaligned snap** |
| `lake_pourpoint.tif` | Buffered lake polygon rasterized as the pour point (lake mode) |
| `dem.tif`, `dem_breached.tif` | DEM / breached DEM, clipped to catchment |
| `flow_pointer.tif`, `flow_accum.tif` | D8 flow pointer and flow accumulation, clipped |
| `streams.tif` / `streams_clipped.tif` | Stream network raster, unclipped/clipped versions |
| `hillshade.tif` | Hillshade, clipped |
| `streams.gpkg` | NHN flowlines clipped to catchment (empty where no NHN coverage) |

Group-level cache files live in `cache/<group_id>/` (or the flat `cache/<PROJECT>/` for whole-domain projects) and are shared across every site in the group.

---

## Correcting a pour point

If a catchment looks wrong, the pour point likely snapped to the wrong stream cell (point mode) or the lake buffer needs adjusting (lake mode):

**Point mode:**
1. Open `output/<PROJECT>/<site_id>/pour_point_snapped.shp` and `streams_tmp.tif` in QGIS
2. Edit `pour_point_snapped.shp` to move the point to the correct stream cell
3. In R:
```r
source("workflow/R/engine/99_rerun_sites")
rerun_engine_site_watershed("site_id", sites, group_manifest, config, output_dir)
```
`rerun_engine_sites()` is the higher-level orchestrator — handles a raw-coordinate fix and/or an edited-snap fix in one call, then cascades remove-upstream/reclip/metrics/hydroweight and merges into the existing combined CSVs.

**Lake mode:**
1. Open `output/<PROJECT>/<site_id>/lake_pourpoint.tif` in QGIS and paint/correct the pour-point cells
2. In R:
```r
source("workflow/R/engine/99_rerun_sites")
rerun_engine_lake_site_watershed("site_id", sites, group_manifest, output_dir)
```

---

## Resetting the workflow

```r
source("code/reset_workflow.R")

reset_all()                     # delete everything in cache/ and output/
reset_sites(sites)               # delete all site output folders, keep cache
reset_group("group_id", sites)   # delete one group cache + its sites
reset_site("site_id")            # delete one site output folder
```

---

## Morphometric metrics

```r
source("workflow/R/catchment_metrics.R")

metrics   <- calculate_catchment_metrics(sites, output_dir = output_dir)  # + lake_polys= for lake mode
ref_table <- build_metrics_reference_table()
write_metrics_outputs(metrics, ref_table, output_dir = output_dir)
```

Computes geometry, areal, and relief metrics for every site (both catchment versions), plus lake-specific metrics (fetch, shoreline development, lake area %, catchment:lake area ratio) when a lake polygon is present. Outputs: `catchment_metrics.csv` and `metrics_reference_table.html` in `output_dir`.

Metrics follow: Shekar, P.R. & Mathew, A. (2024). Morphometric analysis of watersheds: A comprehensive review of data sources, quality, and geospatial techniques. *Watershed Ecology and the Environment*, 6, 13–25.

---

## Hydroweight (distance-weighted landscape attributes)

Each project computes inverse-distance-weighted landscape attributes (e.g. % land cover, mean NDVI) via the `hydroweight` package, under 7 weighting schemes (`lumped`, `iEucO`, `iFLO`, `HAiFLO`, `iEucS`, `iFLS`, `HAiFLS`). The "O" (outflow) target is the site's pour point for stream projects or the lake polygon for lake projects.

Each layer of interest (LOI) — land cover, NDVI, harvest/regen disturbance history, etc. — has its own `workflow/<PROJECT>/prepare_*.R` script that produces a path (or per-group/per-site template) to feed into `loi_layers`, rather than being hardcoded into a runner or the shared hydroweight module. `workflow/raster_attributes.R` is the generic, project-agnostic toolbox those scripts build on (reclassification, Sen's-slope trend rasters, mosaic VRTs, competing-class rasterization). See `CLAUDE.md`'s Hydroweight section for the exact `loi_layers` shape and a list of real `terra`/`hydroweight` quirks worth knowing before adding a new categorical LOI.

---

## Key dependencies and data sources

- **MRDEM**: Canada's 30 m national DEM in EPSG:3979, replacing the CDEM. Streamed via VRT — no tile downloads required. (CELESTE)
- **OIH**: Ontario Integrated Hydrology's pre-conditioned Enforced DEM and Enhanced Flow Direction, EPSG:3161. No breach step needed. (CAM streams, CAM lakes)
- **NHN**: National Hydro Network flowlines and waterbodies used for stream burning. Distributed per NTS sheet as GDB files. (CELESTE)
- **HydroBasins**: HydroSHEDS watershed polygons used to define processing group AOIs. Both North America (`na`) and Arctic (`ar`) datasets required for Canadian sites. (CELESTE)
- **WhiteboxTools**: Open-source geospatial analysis library used for DEM conditioning, flow routing, and watershed delineation.
- **hydroweight**: R package for inverse-distance-weighted landscape attribute summaries.

---

## License

For internal use. Contact Sam W. (Natural Resources Canada) for inquiries.
