# EMILY_TURKEY

Stream-site catchment delineation for the EMILY_TURKEY project, on the
shared, project-agnostic engine described in the repo root `CLAUDE.md`.
Same point-pour-point (Jenson snap), OIH pre-conditioned terrain, EPSG:3161
approach as CAM streams (`workflow/CAM/run_cam_streams.R`) — **except
grouping**: this project's 177 sites split across **two disjoint OIH
terrain tiles with zero overlap**, so it uses the engine's
`"manual_groups"` strategy (one group per tile, each with its own terrain
source) instead of CAM's `"whole_domain"` (a single shared raster set).
This file records what makes this project's run of that engine what it is
— data source, the two-tile discovery, one-off corrections, and the
reasoning behind each — so the run can be explained or reproduced from
scratch. Nothing here should be needed to understand or modify the shared
modules under `workflow/R/`.

## At a glance

| | EMILY_TURKEY |
|---|---|
| Sites | 177 |
| Runner | `workflow/EMILY_TURKEY/run_emily_turkey.R` |
| Delineation | Point pour point (Jenson snap) |
| Terrain | OIH Enforced DEM + Enhanced Flow Direction (pre-conditioned) — **two tiles**, see below |
| Grouping | `"manual_groups"` — group `"NE"` (24 sites) + group `"NC"` (153 sites), auto-detected by terrain coverage |
| Working CRS | EPSG:3161 (NAD83 / Ontario MNR Lambert), both tiles' native CRS |
| Stream threshold | 100 flow-accumulation cells |
| Output | `output/EMILY_TURKEY/` |
| Cache | `cache/EMILY_TURKEY/NE/`, `cache/EMILY_TURKEY/NC/` |
| Hydroweight LOIs | CANLCC land cover only (no NDVI/harvest prep exists yet for this AOI) |

## The two-tile discovery

The first full run failed identically on every site processed
(`wbt_watershed()` producing an all-NA `watershed.tif`, then WhiteboxTools'
`RasterToVectorPolygons` erroring with "The file does not currently
contain any record data"). Traced to: `IntegratedHydrologyNE` (the OIH
tile CAM streams/CAM lakes use) does **not** actually cover most of this
project's AOI — its real data has a large hole there, even though the
sites fall inside its raster's bounding box. A second OIH tile on Sam's
machine, `IntegratedHydrologyNC`
(`/Users/sam/Downloads/IntegratedHydrologyNC/`), covers the rest. Checked
systematically (extract DEM value at every site's coordinate, both tiles):

| | NE tile | NC tile |
|---|---|---|
| Sites covered | 24 | 153 |
| Overlap between tiles | 0 | 0 |

Every site is covered by exactly one tile, none by both, none by neither
— a clean, total split:

- **NE (24 sites):** Garden River (13), Carpenter Lake (4), Northshore (4), Ranger Lake (3)
- **NC (153 sites):** Goulais River (44), Harmony River (34), Christina Mine (14), Wabos (12), Achigan (8), Tabor Creek (7), Batchewana (5), TLW (4), Batchwana (2), Turkey Lake (2), Will (2), Murphy (1), unlabeled (18)

This is why the engine gained a new `grouping$strategy = "manual_groups"`
(`workflow/R/engine/00_resolve_config.R` / `01_build_group_manifest.R` /
`02_prepare_terrain.R`) — the existing `"whole_domain"` assumes one raster
covers every site, and `"hydrobasins"` assumes one shared national mosaic
cropped per region; neither fit a project needing two genuinely
independent, non-overlapping terrain products. Group assignment is
auto-detected (which tile's raster has real data at each site's
coordinate) rather than a manual `sites$group_id` column — deliberately,
since hand-assigning 177 rows from a one-time coverage check would go
stale silently if a site were added later. Validated on a 6-site subset
(3 NE, 3 NC) before the full 177-site rerun: all six delineated with real,
sane catchments (0.16–87.3 km²), none flagged.

## Data inputs

| Input | Path (Sam's machine) | Notes |
|---|---|---|
| OIH Enforced DEM (NE) | `/Users/sam/Downloads/IntegratedHydrologyNE/EnforcedDEM.tif` | Shared with both CAM projects, untouched by this one — covers 24 of 177 sites |
| OIH Enhanced Flow Direction (NE) | `/Users/sam/Downloads/IntegratedHydrologyNE/EnhancedFlowDirection.tif` | Pre-conditioned — no breach step needed |
| OIH Enforced DEM (NC) | `/Users/sam/Downloads/IntegratedHydrologyNC/EnforcedDEM.tif` | Covers the remaining 153 of 177 sites — not used by any other project in this repo yet |
| OIH Enhanced Flow Direction (NC) | `/Users/sam/Downloads/IntegratedHydrologyNC/EnhancedFlowDirection.tif` | Pre-conditioned — no breach step needed |
| Site coordinates | `/Users/sam/Downloads/DataCodes.csv` | Source spreadsheet — N/W decimal-degree WGS84, `Sample Site` as raw identifier |
| CANLCC 2020 | `/Users/sam/Documents/cfs/shared_data/raw/landcover/CAN_LLC_2020.tif` | Land-cover categorical LOI, reused verbatim from CAM |

### OIH → WhiteboxTools flow-direction recoding

Same one-step clockwise rotation matrix as both CAM projects:

```
OIH 128 (NE) -> WBT 1     OIH 1  (E)  -> WBT 2
OIH 2   (SE) -> WBT 4     OIH 4  (S)  -> WBT 8
OIH 8   (SW) -> WBT 16    OIH 16 (W)  -> WBT 32
OIH 32  (NW) -> WBT 64    OIH 64 (N)  -> WBT 128
```

## Site list

177 sites, converted from `DataCodes.csv` into `data/emily_turkey_sites_raw.csv`
(`sample_site, location_project, lat, lon, date_collected, sample_code`) —
the raw `Sample Site` text is kept unmodified in that file; `site_id`
derivation happens in the runner itself so the sanitizing/disambiguation
logic stays visible rather than being baked silently into a pre-cleaned CSV.

### Coverage check

Site coordinates (lat 46.60–47.85, lon -84.69 to -83.26 — Sault Ste
Marie/Algoma, Ontario) confirmed comfortably inside the OIH
`EnforcedDEM`/`EnhancedFlowDirection` extent (lon -84.78 to -78.19, lat
45.21–48.55, EPSG:3161) before committing to reusing OIH for this project.

### site_id derivation and the HARM030 disambiguation

`site_id` is derived from `sample_site` by stripping apostrophes, replacing
whitespace with `_`, then stripping any remaining character outside
`[A-Za-z0-9_-]` (the same sanitize CAM streams applies to its `stream_id`).

One genuine collision survives that sanitize: `"HARM030'"` (Project Sample
Code 149, 46.889617/-84.29985, 29-07-2019) collapses onto `"HARM030"`
(Project Sample Code 148, 46.81364/-84.36151, 26-07-2019). These are two
distinct sites — different coordinates, different collection dates — not a
duplicate entry; the apostrophe in the source spreadsheet was evidently
being used to distinguish a second, nearby site sharing the same base name.
Resolved by appending the Project Sample Code to the second row only,
producing `HARM030_149`. A `stopifnot(!anyDuplicated(site_id))` check in
the runner guards against this kind of collision recurring silently if the
source data changes.

No other `sample_site` value collides after the same sanitize (176 unique
values from 177 rows, fully explained by this one case). One other value,
`"WAB011 MCDCRK"`, needed the whitespace step but produced no collision.

## Hydroweight LOIs

CANLCC (`CAN_LLC_2020.tif`) only, reusing the exact path and
`canlcc_levels` class lookup table already defined in
`workflow/CAM/run_cam_streams.R`. NDVI and harvest/regen LOIs are
deliberately deferred — no per-project `prepare_*.R` data exists yet for
this AOI (Ontario harvest coverage would need to be confirmed against this
project's actual site footprint, same as CAM streams' own harvest-regen
prep did for its AOI, before being added here).

`workflow/EMILY_TURKEY/tidy_outputs.R`'s `tidy_emily_turkey_outputs()`
reshapes `EMILY_TURKEY_hydroweight.csv`'s CANLCC block into a
composition-ready long table (`output/EMILY_TURKEY/tidy/canlcc_long.csv`
— site, version, year, class, scheme, prop), same design as `workflow/
CAM/tidy_outputs.R` / `workflow/CELESTE/tidy_outputs.R`. Only one output
table for now, since CANLCC is this project's only LOI — extend the same
way CELESTE's version extends CAM's when NDVI/harvest-regen prep is
added. See the "newly-found issue" note above: 4 sites' rows in this
table currently report 0% for every class instead of their real
(single-class, 100%) composition.

### hydroweight() crash on an all-NA pour-point column (engine bug, fixed)

The first hydroweight pass silently produced rows for only 155 of 177
sites — no error in the run log, just missing rows, discovered only by
comparing site counts against `catchment_metrics.csv`. Traced to a real,
undocumented behavior of `hydroweight::hydroweight()`: it rasterizes
**every non-geometry column** of `target_O` (the pour point) onto the DEM
grid via `process_input()` -> `terra::rasterize(field = varname)`, one
column at a time — not just its geometry. `pour_point.gpkg` carries
whatever extra columns a project's own sites tibble had; here that's
`location_project`, added purely for traceability (see "Site list" above).
For the 18 sites whose source data had no Location/Project value at all,
that column is entirely `NA`, which `terra` treats as a 0-level factor —
and `hydroweight()` crashes deep inside `levels<-`/`set.cats()` with an
opaque `arguments imply differing number of rows: 2, 0`, caught by this
project's own `tryCatch` only at the outer call, with no per-site
diagnostic beyond a deferred R warning invisible in a non-interactive
`Rscript` log.

Fixed in `workflow/R/stream/hydroweight_attributes.R` (shared by CELESTE
and CAM streams too — they simply never carried an NA-eligible extra
column through to `pour_point.gpkg`, not because they were structurally
immune) by stripping `pour_point_sf` down to just its `site_id` column
(always present, always non-NA, unique per file) before passing it as
`target_O`. Verified two ways before shipping: a regression check
(`ST78`'s already-correct output is bit-identical before/after the fix)
and a recovery check (previously-crashing sites now succeed). Re-running
Stage 7 in full recovered 17 of the 22 originally-missing sites — all of
the `MIC` cluster except `MIC013`.

**5 sites still had no hydroweight row, for an unrelated reason (as of
2026-09-04):** `WABO7161`, `WABO7164`, `CHRIST003`, `CHRIST005`, `MIC013`.
These had genuinely tiny or degenerate catchments (1–3 cells, `CHRIST005`
even computed a *negative* catchment area) and hit a GDAL TIFF
encode/decode failure (`TIFFReadEncodedStrip() failed`) inside
`hydroweight()`'s internal cost-distance raster, or a `[writeRaster]
there are no cell values` failure — confirmed to reproduce even after the
fix above, on catchments too small for a meaningful distance-weighted
landscape summary regardless. All 5 (plus `CHRIST004`/`CHRIST025`/
`ACHI001`/others) overlapped the 11 sites Stage 3 flagged for suspiciously
small catchments.

**2026-09-06: all 11 flagged sites' pour points manually corrected in
QGIS** (`MIC013, CHRIST005, CHRIST003, HWCAR08, CHRIST025, WABO7161,
CUGOUL1, GOUL011, HWGOU3, StreamAspen, HWSTREAM2, WABO7164, ACHI004,
GAR010, HARM015, GOUL005` — 16 total, the 11 flagged plus 5 more caught by
the same review), rerun via `rerun_engine_sites(edited_snap_site_ids =
...)` per root `CLAUDE.md`'s standard workflow. All 16 re-delineated to
healthy 14–58 cell catchments; the correction cascaded to 11 downstream
neighbors via `remove_upstream_catchments()` (`BRGOUL2, GOUL003, GOUL007,
GOUL009, GAR009, GAR011, CHRIST023, HARM004, HARM005, HARM016, HARM017`),
all reclipped/remetriced/rehydroweighted and merged into the combined
CSVs. `ACHI001` (3 cells) and `CHRIST004` (7 cells) remain flagged —
outside this round's correction list.

**WABO7161/WABO7164 still failed hydroweight after the pour-point fix**
(`[writeRaster] there are no cell values`), even at a healthy 27–28
cells — contradicting the "genuinely tiny catchment" explanation above.
Traced to a real, separate bug in `workflow/R/engine/04_delineate_site.R`'s
`watershed_to_polygon()` (used by every project on the engine — CELESTE,
CAM streams, EMILY_TURKEY alike): both sites' D8 trace connects into the
pour-point cell only **diagonally** (a valid 8-connected flow path) —
`whitebox::wbt_raster_to_vector_polygons()` groups such a cell into the
SAME polygon feature as the rest of the catchment, producing a single
self-touching ("bowtie") ring at the corner pinch, invalid per OGC rules.
The existing `sf::st_make_valid()` call resolves that invalidity, but
(confirmed directly — this GEOS version only supports `geos_method =
"valid_structure"`, no alternative) by silently **dropping the smaller
lobe entirely** rather than preserving it as a separate polygon part.
Confirmed on real data: `watershed.tif` had 28 valid cells for
`WABO7161`, the resulting `catchment.gpkg` covered only 27 (24300 m² vs.
the correct 25200 m²) — missing exactly the pour point's own cell,
leaving it ~25m outside its own catchment polygon. Silent at delineation
time; only surfaced later as this opaque `hydroweight()` crash, since the
DEM crop `hydroweight()` aligns to the (undersized) polygon didn't cover
the pour point at all.

Fixed by replacing the whitebox-RTV-+-`st_union()`-+-`st_make_valid()`
approach with `terra::trim()` + `terra::as.polygons(dissolve = TRUE)`
directly on the watershed raster — groups cells by the same
8-connectivity rule but returns a proper, already-valid MULTIPOLYGON
preserving every connected lobe (confirmed: 28/28 cells recovered for
`WABO7161`, pour point correctly contained). Regression-tested against 6
already-correct EMILY_TURKEY sites spanning 1 cell to 790 km² — areas
bit-identical before/after. `WABO7161`/`WABO7164` rerun via the same
`rerun_engine_sites(edited_snap_site_ids = ...)` path (no raw-coordinate
or snap change needed — only the polygonization step was buggy); both now
have complete 4/4 hydroweight rows.

**Separate issue, found and fixed the same day:** 4 sites (`BATCH004`,
`HARM002`, `MIC006`, `MIC016`) had a `class_<scheme>_prop` column (no
class number, un-prefixed) leaking into `EMILY_TURKEY_hydroweight.csv`
alongside a full set of `canlcc_*_prop` columns that were incorrectly
all-zero — i.e. these sites' real land-cover composition (100% of one
class, per the leaked value) was being reported as 0% for every class
instead. Root cause: `single_class_value()`
(`workflow/R/stream/hydroweight_attributes.R`) — meant to intercept
exactly this single-class-ROI case and produce a clean, correctly-named
1.0/0.0 table via `degenerate_categorical_table()` — evaluated "which
classes are present" over a different (looser) mask than
`hydroweight::hydroweight_attributes()` uses internally (the LOI's own
native-resolution crop/mask vs. the package's internal nearest-neighbor
resample onto the distance-weight grid), so the two could disagree right
at a catchment boundary; when they did, the raw ambiguous
`"Class_{scheme}_prop"` (no id) output described in root `CLAUDE.md`'s
"known quirks" section leaked straight through un-prefixed, and
`ensure_full_categorical_schema()` then zero-filled every real class
column for that row since none of them existed under the expected name.

Fixed by adding `align_loi_to_distance_weight_grid()` — resamples/crops/
masks the LOI onto one of the site's own distance-weight rasters' exact
grid (the same grid `hydroweight_attributes()` uses internally) *before*
`single_class_value()` checks it, so our own degenerate-ROI detection now
fires whenever the package's fallthrough would otherwise trigger.
Verified: all 4 sites now report the correct 1.0 for their true sole
class (ID 5, `broadleaf_deciduous_forest`) and 0.0 elsewhere, zero `class_`
columns remain, and two unrelated multi-class sites (`ST78`, `GAR012`)
are bit-identical before/after. Full Stage 7 rerun for all 177 sites
confirmed no other site was hitting this.

This lives in `workflow/R/stream/hydroweight_attributes.R` — shared by
CELESTE and CAM streams too, both automatically covered by the fix going
forward. Checked both projects' existing hydroweight CSVs directly: no
leaked `Class_`/`class_` columns in either, so neither has actually hit
this in production. `workflow/R/lake/hydroweight_attributes.R` (CAM
lakes) is a separate file with **no equivalent protection at all** — no
`single_class_value()`/degenerate bypass, no `ensure_full_categorical_
schema()` zero-fill — structurally exposed the same way, currently no
evidence of having hit it (45/45 sites clean), deliberately left
unfixed for now since CAM lakes is still a work in progress (Sam,
2026-09-06) — porting this fix there is a known follow-up for whenever
that project is revisited.

## The "clipped" drainage_density/stream_frequency-always-NA bug (found + fixed 2026-09-06)

Every single "clipped"-version row in `catchment_metrics.csv` had
`drainage_density_km_km2`/`stream_frequency_per_km2` = `NA` (68/68 rows,
100%) — indistinguishable at a glance from "genuinely no stream here,"
but actually a real bug, not a physical result. Root cause:
`catchment_metrics.R`'s `compute_site_metrics()` auto-detects a per-site
stream layer for the `clipped` version, preferring `streams_clipped.gpkg`
(vector) over `streams_clipped.tif` (raster) whenever the `.gpkg` merely
*exists* — it never checks feature count. `workflow/R/reclip_outputs.R`'s
`reclip_site()` was unconditionally calling
`clip_flowlines_to_catchment_clipped()` regardless of whether the group
actually has burned-in NHN flowlines cached — unlike the engine's own
delineation stage and `rerun_engine_sites()`, which both already gate
this on `cache_exists(flowlines.gpkg)`. EMILY_TURKEY has no NHN burn-in
(`streams_burn = list(source = "none")` — OIH terrain only), so every
site got an **empty** `streams_clipped.gpkg` (0 features) written
regardless — which then won the auto-detect over the real, non-empty
`streams_clipped.tif` sitting right next to it (confirmed directly:
`HWGOU04`/`MGOU1`/`GOULST001` each had 0 gpkg features but 166/147/550
real stream cells in their `.tif`), discarding perfectly good data for
every clipped catchment.

Fixed by gating `reclip_site()`'s flowline-clip call on
`cache_exists(flowlines.gpkg)` (matching the guard used elsewhere), then
deleting the 177 stale empty `streams_clipped.gpkg` files (confirmed
every single one was 0-feature before deleting — not assumed) and
rerunning Stage 6 (metrics only; Stage 5's raster outputs were never
wrong, only the flowlines vector was). Verified: 0/68 clipped rows NA
afterward (67 real nonzero values, 1 legitimate 0 for `CHRIST003`'s tiny
clipped remainder — its stream cells genuinely fall below the flow-
accumulation threshold, a real physical result, not this bug). The
`unclipped` version's own 22 legitimate zeros (tiny headwater catchments
below the 100-cell stream threshold) are unrelated and unaffected either
way. Checked CAM streams (identical no-burn-in setup) and found the
identical bug there too (8/8 clipped rows NA) — fixed and reran in the
same pass, see `workflow/CAM/README.md`. CELESTE (real NHN burn-in) was
never affected — its `streams_clipped.gpkg` genuinely has features
wherever a real flowline exists, so the vector-preferred design is
correct there, not a bug.

## Reproducing this run from scratch

1. `source("workflow/EMILY_TURKEY/run_emily_turkey.R")` through Stage 1 —
   confirm 177 sites resolve into 2 groups (`"NE"`: 24, `"NC"`: 153), no
   CRS-suitability warning, no aborted/ambiguous-site error.
2. Stage 2 (terrain prep) — spot-check `cache/EMILY_TURKEY/NE/streams.tif` /
   `dem.tif` and `cache/EMILY_TURKEY/NC/streams.tif` / `dem.tif` in QGIS.
   Adjust `STREAM_THRESHOLD` and re-run if the extracted stream network
   looks wrong (delete the relevant `streams.tif` first to force
   regeneration).
3. Stage 3 (delineation) — review the `flagged` sites table (small
   catchments / edge cases) before trusting Stages 4–7.
4. Stages 4–6 (remove-upstream, reclip, metrics) — no project-specific
   decisions, same as CAM streams.
5. Stage 7 (hydroweight) — CANLCC only, as above.
