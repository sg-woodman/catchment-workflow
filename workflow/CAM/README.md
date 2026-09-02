# CAM

Two related but distinct delineation runs for the CAM project, both on the
shared, project-agnostic engine described in the repo root `CLAUDE.md`:
**CAM lakes** (45 lake-polygon pour points) and **CAM streams** (39 point
pour points). Both use the same underlying OIH terrain data. This file
records what makes CAM's runs of that engine what they are — data sources,
config values, one-off corrections, and the reasoning behind each — so the
runs can be explained (e.g. for a methods section) or reproduced from
scratch. It is the CAM-specific counterpart to the generic architecture
documented in `CLAUDE.md`; nothing here should be needed to understand or
modify the shared modules under `workflow/R/`.

## At a glance

| | CAM lakes | CAM streams |
|---|---|---|
| Sites | 45 lakes | 39 (25 "Summer"/CRADLES-lakes-adjacent + 14 "Fall"/MOE long-term monitoring gauges) |
| Runner | `workflow/run_cam_lakes.R` | `workflow/CAM/run_cam_streams.R` |
| Delineation | Lake-polygon pour point (buffered 30 m) | Point pour point (Jenson snap) |
| Terrain | OIH Enforced DEM + Enhanced Flow Direction (pre-conditioned) | Same OIH data |
| Grouping | `"whole_domain"` (one flat group) | `"whole_domain"` (one flat group) |
| Working CRS | EPSG:3161 (NAD83 / Ontario MNR Lambert), OIH's native CRS | EPSG:3161 |
| Stream threshold | 100 flow-accumulation cells | 100 cells |
| Output | `output/CAM/lake_delineation_engine/` | `output/CAM/stream_delineation/` |
| Cache | `cache/CAM_lakes_engine/` | `cache/CAM/stream_delineation/` |

## Data inputs

| Input | Path (Sam's machine) | Notes |
|---|---|---|
| OIH Enforced DEM | `/Users/sam/Downloads/IntegratedHydrologyNE/EnforcedDEM.tif` | Shared by both CAM lakes and CAM streams, untouched by either |
| OIH Enhanced Flow Direction | `/Users/sam/Downloads/IntegratedHydrologyNE/EnhancedFlowDirection.tif` | Pre-conditioned — no breach step needed |
| OIH waterbodies | `/Users/sam/Documents/cfs/shared_data/raw/hydro/waterbodies/ON_waterbodies/ohn_waterbodies_valid.gpkg` | OHN lake polygons, lake matching + lake-bisection checks |
| CAM streams sites | `data/cam_stream_sites_raw.csv` | Converted from `Stream Coordinates Summer and Fall - Cameron Lefebvre.xlsx`, Sheet1 |
| CANLCC 2020 | `/Users/sam/Documents/cfs/shared_data/raw/landcover/CAN_LLC_2020.tif` | Land-cover categorical LOI, both CAM lakes and CAM streams |
| NDVI (lakes) | `/Users/sam/Downloads/NDVI_CAM_sites_BAP.tif` | Built for lake-site locations only — see "NDVI" below |
| NDVI (streams) | `data/ndvi/CAM/*.tif` | Per-catchment/group Google Earth Engine exports |
| Harvest (streams) | `shared_data/raw/harvest/ontario_harvest.gdb` | Ontario MNR only — CAM's AOI has no NB/Irving relevance |

### OIH → WhiteboxTools flow-direction recoding

Both CAM lakes and CAM streams recode OIH's Enhanced Flow Direction to
WhiteboxTools' D8 encoding via the same one-step clockwise rotation
matrix:

```
OIH 128 (NE) -> WBT 1     OIH 1  (E)  -> WBT 2
OIH 2   (SE) -> WBT 4     OIH 4  (S)  -> WBT 8
OIH 8   (SW) -> WBT 16    OIH 16 (W)  -> WBT 32
OIH 32  (NW) -> WBT 64    OIH 64 (N)  -> WBT 128
```

## CAM lakes (`workflow/run_cam_lakes.R`)

### Site list

45 lakes, hardcoded as a `tibble::tribble()` in the runner (`cam_sites_raw`)
— this is project data, not pipeline-specific, so it's the same list
across every pipeline version this project has run on. Two supporting
tables handle imperfect name/spatial matching against the OHN polygon
layer:

- `manual_id_lookup` — 8 lakes matched to a specific OGF_ID by hand
  (Marina, Tilton, Middle, Hannah, Lohi, Jerry, Daisy, Chiniguchi).
- `excluded_ogf_ids` — 2 OGF_IDs excluded from automatic matching.

Matching itself (`match_lake_polygons()`, `workflow/R/lake/`) uses a 10 m
shoreline snap buffer, Jaro-Winkler name distance ≤ 0.15, and a 40 km name
search radius.

### Excluded waterbody types

`c("Pond")` for upstream-lake removal (`remove_upstream_lake_catchments()`)
— narrower than CAM streams' lake-bisection check (`c("River", "Pond")`,
see below), since this is filtering *upstream lake* candidates for the
target lakes themselves, not filtering candidate bisected lakes.

### Migration provenance

Migrated from a non-modular `workflow/R/lake/`-only pipeline (retired —
see git history) once the engine was validated as an exact geometric
match: 45/45 sites, IoU = 1.0000, 0% area difference on every site — OIH
needs no breach step, so this migration avoided the MRDEM-breach
nondeterminism CELESTE's own migration had to work through (see
`CLAUDE.md`'s note on WhiteboxTools breach nondeterminism). Hydroweight
canlcc values matched production to within ~0.1% (normal
raster-resampling noise).

Reused **completely unmodified** from the old pipeline (confirmed by
reading each function's actual file dependencies, not assumed):
`match_lake_polygons()` (pure spatial/name matching, no cache_dir
dependency), and `calculate_hydroweight_attributes()` (defaults to
`cache_dir/dem_breached.tif` + `cache_dir/flow_accum.tif`, both already
materialized by the engine's terrain prep). The **one** real gap:
`remove_upstream_lake_catchments()` hardcodes `cache_dir/d8_pntr.tif` (the
old pipeline's own naming), while the engine writes `flow_pointer.tif` —
fixed with a cheap symlink alias created right after Stage 2, rather than
modifying the shared function itself.

### Reconciling with `run_cam_streams.R`

`run_cam_streams.R`'s stream-monitoring sites are near, but do not match,
these 45 lake sites — they were delineated separately (point pour points,
not lake polygons) against the same physical OIH data, under
`workflow/CAM/` rather than at the repo's top level, purely because CAM
lakes doesn't have its own `workflow/CAM_lakes/` subdirectory yet (see
`CLAUDE.md`).

### For correcting a single bad lake catchment

Edit `output/<site_id>/lake_pourpoint.tif` in QGIS, then use
`rerun_engine_lake_site_watershed()` (`workflow/R/engine/99_rerun_sites`)
— the same tool the old pipeline's `rerun_watershed_lake()` provided,
ported to the engine specifically so that capability wasn't lost when the
old pipeline was retired.

## CAM streams (`workflow/CAM/run_cam_streams.R`)

### Why this isn't on either existing pipeline

These 39 sites needed point pour points (Jenson snap), like CELESTE, but
against the *same* OIH DEM/flow-direction CAM lakes uses (EPSG:3161), not
MRDEM (EPSG:3979) — a combination neither the old CELESTE pipeline nor the
old CAM-lakes pipeline supported directly. This combination is exactly
what the shared engine (`workflow/R/engine/`) was generalized to handle.

### Site list and the SUD11 exclusion

Converted from `Stream Coordinates Summer and Fall - Cameron
Lefebvre.xlsx`, Sheet1, via `data/cam_stream_sites_raw.csv` (39 sites: 25
"Summer"/CRADLES-lakes-adjacent + 14 "Fall"/MOE long-term monitoring
gauges).

**`SUD11` is excluded** (`EXCLUDED_SITE_IDS`) — source row has
`lon = -51.2`, 2500+ km off from its neighbors `SUD12`/`VER01` (both
~-81.0), almost certainly a typo in the source workbook. Verify against
Cameron Lefebvre's original data before re-including.

### Terrain

No stream burn-in (`streams_burn = list(source = "none")`) — OIH's flow
direction is already hydrologically conditioned, and burning would have
to happen before conditioning, not after. `dem` (OIH Enforced DEM) is
supplied purely as the elevation surface for per-site clipping/output +
hydroweight, decoupled from `flow_direction` (which drives conditioning).

### Lake-bisection fix (`fix_lake_bisection.R`) — CAM streams only

Root cause: point-based delineation has zero lake awareness —
`wbt_watershed()` traces strictly by per-cell D8 flow direction, and
OIH's own hydro-enforcement only improves flow direction for waterbodies
intersecting its mapped Enhanced Watercourse network, not universally.

Confirmed on the delivered CAM streams output (2026-08-27): 11 (site,
lake) pairs across 7 sites — Tilton, Daisy, NCMN, SUD17, SUD102 (×3),
SUD103 (×1), SUD200 (×3) — genuinely bisect a lake ≥ 1 ha, excluding
River/Pond waterbody types and boundary-touch noise (<2%/>98% overlap).

**Scope note**: this fix predates CELESTE's own upfront lake-conditioning
approach (see `workflow/CELESTE/README.md`) — CAM streams has **no
upfront-prevention equivalent by design**, only this reactive,
manually-gated fix. CAM streams uses a different terrain tier (OIH's
pre-conditioned `flow_direction`, not a raw DEM needing breach) — upfront
lake conditioning only applies to the raw-`dem` terrain tier (flattening
must happen before D8 derivation), so it isn't available for this
project's terrain source.

Usage:

```r
source(here::here("workflow/CAM/run_cam_streams.R"))
source(here::here("workflow/CAM/fix_lake_bisection.R")) # validates + reports only
# review output/CAM/stream_delineation/lake_intersection_report.csv and
# the catchments in QGIS against ohn_waterbodies_valid.gpkg
correct_lake_bisected_sites(
  site_ids = c("Tilton", "SUD17"),
  sites = sites, group_manifest = group_manifest, config = config,
  output_dir = output_dir, cache_dir = cache_dir, loi_layers = loi_layers
)
```

Safe to call repeatedly, including as a second/third pass — accumulates
flattened lakes across passes (never delete `cache_dir/lake_corrected/`
between passes; confirmed doing so once reverted 4 sites' fixes).
Densely-lake-packed clusters (confirmed: SUD17/SUD102/SUD103/SUD200) can
show a whack-a-mole pattern across passes — review the report after each
pass rather than assuming convergence to zero; a small residual (e.g. a
sub-2-ha lake right at the edge of the 98% "acceptable" band) may be worth
leaving as manual-review-only.

### 2026-08-31/09-01 clean rerun and manual corrections

The originally-delivered CAM streams output (2026-08-26) was superseded
by a clean-cache rerun, after which 7 of the lake-bisection sites reverted
to their pre-correction state (expected — a clean rerun starts from the
uncorrected flow field again). Sam's manual review of that revert
(`workflow/CAM/one_off/rerun_corrected_sites_20260901.R`) resolved them as
follows:

- **NCMN, SUD17, SUD200** — left as-is (Sam's call after manual QGIS
  review; not every flagged bisection is a genuine defect).
- **Daisy, Wolf** — `pour_point_snapped.shp` manually corrected by Sam in
  QGIS (confirmed via file mtime); rerun via `edited_snap_site_ids`, no
  lake correction involved.
- **SUD102, SUD103, Tilton** — genuine lake bisection; rerun via
  `fix_lake_bisection.R`'s `correct_lake_bisected_sites()`.

### NDVI: superseded source (2026-08-26)

The `"ndvi"` LOI was originally wired to
`/Users/sam/Downloads/NDVI_CAM_sites_BAP.tif`, a raster built for the
*lake*-site locations — that source genuinely doesn't cover `SUD22` (~8 km
outside its extent), so those columns read NA under it. Rather than leave
`SUD22` uncovered, `workflow/CAM/prepare_ndvi.R`'s
`prepare_cam_ndvi_site_rasters()` now builds a properly
stream-catchment-matched *continuous* NDVI layer from the same
per-catchment/group Google Earth Engine exports (`data/ndvi/CAM/*.tif`)
that `"ndvi_trend"` uses — confirmed directly that `SUD22` maps to
`Landsat_NDVI_1984_2025_cam_sud17.tif` (shared with `SUD17`) with full
coverage. `"ndvi_trend"` is a separate Sen's-slope trend computed from the
same export data (mirrors `workflow/CELESTE/prepare_ndvi_trend.R`), kept
alongside `"ndvi"` since the two answer different questions (a static
mean/distribution vs. a directional trend).

### NDVI lake masking — a separate, additional LOI, not a default

`ndvi_masked`/`ndvi_trend_masked` (`workflow/CAM/prepare_ndvi_masked.R`)
are lake-masked counterparts of the continuous `"ndvi"` and `"ndvi_trend"`
LOIs — separate LOIs, not a toggle applied in place. `"ndvi"`/
`"ndvi_trend"` themselves are never masked; masking writes to its own
`hydroweight_loi/ndvi_masked_clean/` cache, built from the plain
(unmasked) `hydroweight_loi/ndvi_clean/` file — not from the raw export
directly — so the 0→NA/10000 rescale fix stays defined in exactly one
place (`clean_cam_ndvi_tile()`).

**Deliberately a separate LOI, not an in-place toggle** — an earlier
version masked `"ndvi"`/`"ndvi_trend"` directly, overwriting their cache
in place. See `workflow/CELESTE/README.md`'s equivalent note for why that
turned out to be the wrong call: it makes a real, unrelated problem look
like a masking bug, since there's no unmasked baseline left to compare
against.

Masked using the OIH/OHN waterbody layer (`ohn_waterbodies_valid.gpkg`,
`CAM_OIH_LAKES_PATH` in `prepare_ndvi_masked.R`), read via an
extent-filtered `terra::vect()` scoped to each raw file's own footprint —
the same GDAL bbox-pushdown pattern `workflow/R/lake/
match_lake_polygons.R` and `remove_upstream_lakes.R` already use against
this exact 1.4M-feature province-wide file — rather than loading it
whole. Unlike this project's lake-*bisection* check (which excludes
`"River"`/`"Pond"` because those aren't a "lake" for that purpose),
masking excludes nothing by type by default (`exclude_waterbody_types =
character(0)`) — a pond or river surface is still open water and still
depresses NDVI just like a lake does.

**A real operational finding, not a masking bug**: an initial full-batch
run of all 4 ndvi-family LOIs together across all 38 sites returned no
`ndvi`/`ndvi_trend`/`ndvi_masked`/`ndvi_trend_masked` data at all for 17
scattered site × version rows (canlcc/harvest_regen were fine for every
one of them, so `hydroweight::hydroweight()` itself was succeeding). This
looked at first like a genuine "lake-outlet catchment masked down to
nothing" case — several of the affected sites (George, Manitou, Wolf,
etc.) are pour points right at or below a named lake's outlet (the
"CRADLES-lakes-adjacent" sites — see "Site list" above), a plausible
mechanism. It wasn't: re-running `calculate_hydroweight_attributes_stream()`
for just those 17 site × version combinations resolved every single one,
including ones with no lake-outlet story at all (`Bell`, `Marina`,
`Smoothwater`, `SUD12`). Retrying the exact same computation on a smaller
batch producing correct results rules out a real per-site data gap — this
is the same terra temp-file/GC-race class of issue already documented
elsewhere in this project for large multi-LOI runs, here manifesting as
silently-dropped rows in a `path_template`-based call rather than a hard
crash. If a run ever produces a similarly scattered, small set of
unexplained NA rows, re-running that subset is the fix — not a sign of a
real coverage gap.

### LOI layers use `path_lazy`/`path_template`, not `path`

CAM streams' project-wide LOI rasters (CANLCC, harvest_regen) are declared
as `path_lazy`, and the per-site NDVI/NDVI-trend rasters as
`path_template`, deliberately — not just for consistency with CELESTE.
With `path`, `calculate_hydroweight_attributes_stream()` prepares ONE
shared SpatRaster per LOI up front and reuses that SAME object across
every site × version job (76 here: 38 sites × 2 versions); each job's
crop()/mask() then operates on the full project-wide extent every time.
Confirmed in practice to eventually crash with a terra `"[readStart] file
does not exist"` error partway through a run (a temp-file GC race from
heavy repeated crop/mask churn on a long-lived shared object).
`path_lazy`/`path_template` avoid this structurally: each site's LOI is
pre-cropped down to just its own catchment once, cached to a real per-site
file, before further processing touches it.

### Rerunning after a correction

```r
# Raw coordinate fixed in data/cam_stream_sites_raw.csv:
rerun_engine_sites(
  resnap_site_ids = "Bell", sites = sites, group_manifest = group_manifest,
  config = config, output_dir = output_dir, cache_dir = cache_dir,
  loi_layers = loi_layers # or NULL to skip re-running hydroweight
)

# Snapped pour point manually edited in QGIS
# (edit output/<site_id>/pour_point_snapped.shp directly first):
rerun_engine_sites(
  edited_snap_site_ids = "Wolf", sites = sites, group_manifest = group_manifest,
  config = config, output_dir = output_dir, cache_dir = cache_dir,
  loi_layers = loi_layers
)
```

## The CAM lakes Stage 5 reclip bug (found 2026-08-31)

`reclip_outputs()`'s lake branch called `load_lake_cache_rasters(cache_dir)`
— a function that, at the time, was never defined anywhere in the repo.
Every call raised `"could not find function"`, caught by `reclip_site()`'s
own `tryCatch` and downgraded to a per-site warning rather than an abort,
so `run_cam_lakes.R` "completed successfully" while silently producing
**zero** `*_clipped.tif` rasters for all 45 sites, every run. Found while
checking why CAM lakes' output had no `dem_clipped.tif` at all. Fixed in
`workflow/R/reclip_outputs.R` (`load_lake_cache_rasters()` now defined —
see `CLAUDE.md`'s "Reclip systemic-failure detection" note for the general
hardening this prompted across every project). Backfilled via
`workflow/one_off/backfill_cam_lakes_reclip.R`, which reruns Stage 5
(fixed), Stage 6 (metrics — the "clipped" version rows were reading a
`dem_clipped.tif` that never existed), and Stage 7 (hydroweight — the
"clipped" version was silently missing its stream-based weighting
schemes for the same reason), skipping the unaffected Stage 1–4 outputs.

## Output shapes

`workflow/CAM/tidy_outputs.R`'s `tidy_cam_outputs()` reshapes
`catchment_metrics.csv` and `CAM_streams_hydroweight.csv` (both wide, one
row per site × version) into 5 small, purpose-shaped long tables — one per
LOI (canlcc/ndvi/ndvi_trend/harvest_regen) plus one for the morphometric
metrics — rather than one wide row per site or one generic melt. The
original wide CSVs are untouched and remain the design matrix for
statistical modelling. See that file's header for the exact year-handling
rules per table (canlcc is a single fixed year 2020; ndvi/harvest_regen
carry real per-year columns; ndvi_trend/harvest_regen_combined have no
year at all, being fit across their whole record).

## Reproducing from scratch

```r
# Lakes:
source("workflow/run_cam_lakes.R")

# Streams:
source("workflow/CAM/run_cam_streams.R")
source("workflow/CAM/fix_lake_bisection.R") # check for residual lake bisection
```

Before running either: confirm the OIH DEM/flow-direction paths and (for
streams) `data/cam_stream_sites_raw.csv` match your machine. To rerun a
single site after a manual pour-point correction in QGIS, see
`workflow/R/engine/99_rerun_sites` — same pattern described in the repo
root `CLAUDE.md`.
