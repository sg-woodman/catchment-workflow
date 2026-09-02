# CELESTE

River/stream-site catchment delineation for the CELESTE project: 132 sites
across 6 HydroBasins groups (COC, KEN, MOR, NBE, NIP, TUR), on the shared,
project-agnostic engine described in the repo root `CLAUDE.md`. This file
records what makes CELESTE's *specific* run of that engine what it is —
data sources, config values, one-off corrections, and the reasoning behind
each — so the run can be explained (e.g. for a methods section) or
reproduced from scratch. It is the CELESTE-specific counterpart to the
generic architecture documented in `CLAUDE.md`; nothing here should be
needed to understand or modify the shared modules under `workflow/R/`.

## At a glance

| | |
|---|---|
| Sites | 132, across 6 HydroBasins level-6 groups |
| Delineation | Point pour point (Jenson snap), `workflow/CELESTE/run_celeste.R` |
| Terrain | MRDEM 30 m DEM (raw, needs breach), `~/Documents/cfs/shared_data/raw/dem/mrdem-30-dtm.vrt` |
| Grouping | `"hydrobasins"` — HydroBasins level-6 union per group, 1000 m default buffer |
| Stream burn-in | NHN (`streams_burn = list(source = "nhn_auto")`), per-group AOI, auto-downloaded/cached |
| Lake conditioning | Upfront, opt-in (`lake_conditioning = list(source = "nhn_auto")`) — see below |
| Working CRS | EPSG:3979 (MRDEM's native CRS, Canada Atlas Lambert NAD83(CSRS)) — no override |
| Stream threshold | 1000 flow-accumulation cells |
| Site source | `data/celeste_milli_sites_clean_corrected.gpkg` |
| Output | `output/CELESTE/`, cache `cache/CELESTE/` |
| Runtime | ~90 min on a cold cache (6 groups; MOR's DEM crop is ~5x the smallest group's) |

## Data inputs

| Input | Path (Sam's machine) | Notes |
|---|---|---|
| MRDEM DEM | `~/Documents/cfs/shared_data/raw/dem/mrdem-30-dtm.vrt` | 30 m national DEM, streamed via VRT, native EPSG:3979 |
| HydroBasins | `/Users/sam/Documents/cfs/shared_data/raw/hydro/watersheds/HydroBasins` | `north_america/` + `arctic/` subfolders, level 1–12 |
| NHN GDBs | `/Users/sam/Documents/cfs/shared_data/raw/hydro/networks/NHN/gdb` | Per-work-unit File Geodatabases, auto-downloaded/cached here |
| NHN index | `.../NHN/NHN_INDEX_WORKUNIT_LIMIT_2/NHN_INDEX_22_INDEX_WORKUNIT_LIMIT_2.shp` | Maps NHN work-unit sheet codes to geometries |
| Sites | `data/celeste_milli_sites_clean_corrected.gpkg` | See "Site coordinates" below |
| NDVI tiles | `data/ndvi/*.tif` | Regional Landsat NDVI composites, 1984–2025, per-group filename match |
| Harvest (Ontario) | `shared_data/raw/harvest/ontario_harvest.gdb` | NIP/TUR/KEN coverage, `AR_YEAR` field |
| Harvest (NB/Irving) | `data/irving_harvest/LB_HarvCuHi_RefCuHi_SICuHi.gdb` | NBE + COC coverage, `HARVYR`/`RFYR` fields |
| CANLCC 2020 | `/Users/sam/Documents/cfs/shared_data/raw/landcover/CAN_LLC_2020.tif` | National land-cover product, same source + class scheme as both CAM projects; already EPSG:3979, no reprojection needed |

## Site coordinates

Uses `data/celeste_milli_sites_clean_corrected.gpkg` — confirmed by the
user to match what they were actually sent. Includes 6 sites whose pour
points were manually corrected in QGIS after the initial engine run
flagged them as near-zero catchments; see
`data/celeste_milli_sites_corrections_2026-08-26.csv` for the audit trail
of exactly what changed and why.

## How it was run (`run_celeste.R`, Stages 1–8)

1. **Config + group manifest** — `grouping$strategy = "hydrobasins"`,
   delegating to `workflow/R/stream/group_sites.R` unmodified.
2. **Terrain prep** — crop MRDEM per group, burn NHN streams in (see
   below), breach, upfront lake conditioning (see below), D8 pointer, flow
   accumulation, extract streams.
3. **Delineate** — point pour point, Jenson snap, `snap_dist = 200`,
   `min_cells = 10`.
4. **Remove upstream** nested catchments.
5. **Re-clip** rasters/flowlines to clipped catchments — **required**
   before Stage 7 hydroweighting, not optional: the clipped version's
   stream-based weighting schemes (iEucS/iFLS/HAiFLS) silently come out
   missing (not erroring) if `streams_clipped.tif` doesn't exist yet.
6. **Morphometric metrics**, unclipped + clipped.
7. **LOI prep + hydroweighting** — see "LOI layers" below.
8. **Tidy** — long-format plotting tables via `tidy_celeste_outputs()`.

### COC runs without stream burn-in

`COC-MAIN`/`COC-NWB` are overridden to `burn_streams = FALSE` in the sites
tibble inside `run_celeste.R` (not in the source `.gpkg` — that file is
the user's authoritative reference data). This is a real, confirmed
correctness fix, not cosmetic: `COC-NWB`'s catchment was 248.07 km²
burned-in (wrong — connected to a drainage path that isn't actually
there) vs. 66.33 km² without burn-in (correct).

### Stream burn-in (NHN)

`streams_burn = list(source = "nhn_auto")`, resolved **per group** from
each group's own HydroBasins-derived AOI rather than a fixed per-group
vector — necessary because a corrected site coordinate can, in principle,
shift which NHN work units apply. `download_nhn()` checks the local NHN
directory for an already-downloaded GDB matching the required WSCSSDA
code before ever hitting the FTP server, so a rerun with everything
already cached needs no network access.

### Upfront lake conditioning (opt-in, `lake_conditioning = list(source = "nhn_auto")`)

Flattens known NHN lakes into the DEM before D8 derivation, so every site
in a group is delineated from one internally-consistent, lake-aware flow
field — instead of reactively correcting individual bisected catchments
after the fact (`fix_lake_bisection.R`, still kept as a fallback for any
residual — see below). Verified clean on 2026-08-29:

- KEN: 3/3 originally-bisected sites fully resolved.
- NIP + TUR: went from 31 bisected (site, lake) pairs across 28 sites down
  to 2, both at the documented noise-band edge (S3/S4, Little Turkey
  Lake).
- COC/MOR/NBE: 0 bisections either before or after (don't need
  conditioning, but still pay a small, bounded lake-search cost since
  `lake_conditioning` applies uniformly per group — no per-group opt-out
  exists; not worth building for groups this small).

The lake search is deliberately scoped to a buffer around a group's own
sites (default 20 km), not its full HydroBasins AOI — confirmed directly
(2026-08-29) that the latter is wildly disproportionate: for KEN, the
HydroBasins group AOI is 58,928 km² vs. 2,663 km² for the group's own 21
sites' bare bounding box (~22x smaller), and a first attempt scoped to the
full group AOI pulled in 31,525 candidate lakes and took ~68 minutes for
one group.

### Reactive lake-bisection fix (`fix_lake_bisection.R`) — fallback only

Superseded as the *primary* fix by upfront lake conditioning above, but
kept as a fallback for any residual bisection. A full systematic scan
found 33 (site, lake) pairs across 28 of 132 sites, in 3 of 6 HydroBasins
groups (KEN, NIP, TUR — COC/MOR/NBE confirmed clean). This exceeded a set
of 12 sites spotted manually in QGIS beforehand (11 confirmed genuine;
`BLOXHAM1_1K_1` checked out clean — likely conflated with the genuinely
bisected `BLOXHAM1_10K_1`).

Two structural differences from the CAM streams version of this fix (see
`workflow/CAM/README.md`), both handled inside `fix_lake_bisection.R`
rather than by changing the shared `workflow/R/lake_containment.R`:

1. **Lake source is NHN** (`NHN_HD_WATERBODY_2`), read per-HydroBasins-group
   from local GDBs, not one province-wide file. Deliberately does **not**
   reuse `stream/burn_streams.R`'s `read_merge_nhn_layer()` — that
   function clips each waterbody's geometry to its query AOI, which
   distorts the containment percentage this check depends on. Confirmed
   directly: `BLOXHAM1_1K_1` read as 0.17% outside with a tight per-site
   AOI and 7.3% outside with a larger group-wide AOI, off the *same* lake
   — purely an artifact of how much of it each AOI happened to clip off.
   `fetch_nhn_lakes_for_group()` reads unclipped and only uses the AOI for
   a non-destructive row filter.
2. **Terrain tier is raw MRDEM + per-group breach**, with NHN streams
   burned in before breaching for 5 of the 6 groups (`dem_burned.tif`
   exists for KEN/MOR/NBE/NIP/TUR — exactly the 5 groups with confirmed
   bisections; COC has `burn_streams = FALSE`, see above, and only has
   `dem.tif`).

Usage:

```r
source(here::here("workflow/CELESTE/run_celeste.R"))
source(here::here("workflow/CELESTE/fix_lake_bisection.R")) # validates + reports only
# review output/CELESTE/lake_intersection_report.csv and the catchments in QGIS
correct_lake_bisected_sites_celeste(
  site_ids = c("CF1", "SN1UP"), # can span multiple HydroBasins groups
  sites = sites, group_manifest = group_manifest, config = config,
  output_dir = output_dir, cache_dir = cache_dir, loi_layers = loi_layers
)
```

`exclude_waterbody_types = c("Watercourse")` — NHN's own waterbody-type
value for a feature that isn't a real lake, the CELESTE equivalent of
CAM's `c("River", "Pond")` (OHN's own vocabulary).

### `remove_upstream.R` integrity-check case studies

The clip-integrity check in `workflow/R/remove_upstream.R` (fragmentation
/ pour-point-exclusion guard) was added after two real CELESTE cases where
a lake-bisection correction on one site left an unrelated neighbor's
nesting relationship non-clean: a correction on `BAT2NEW` produced a
3-part disconnected `catchment_clipped.gpkg`, and one on `SN1UP` produced
a clipped catchment that no longer contained its own pour point — both
silently written before this check existed.

## LOI layers (Stage 7)

| LOI | Type | Source | Cache scope |
|---|---|---|---|
| `canlcc` | categorical | `CAN_LLC_2020.tif` (same national product + class scheme as both CAM projects), declared `path_lazy` directly in `run_celeste.R` — no `prepare_*.R` script needed | Cropped straight from the source per site; nothing cached at the group/project level |
| `ndvi` | continuous | `workflow/CELESTE/prepare_ndvi.R`, mosaics `data/ndvi/*.tif` per group | Depends only on `group_id` — reused directly across reruns, no DEM dependency |
| `ndvi_trend` | continuous, 2-band (slope, p-value) | `workflow/CELESTE/prepare_ndvi_trend.R`, Theil-Sen/Mann-Kendall over the full 1984–2025 series | Same as `ndvi` — no DEM dependency |
| `ndvi_masked` | continuous | `workflow/CELESTE/prepare_ndvi_masked.R` — lake-masked counterpart of `ndvi`, same source tiles | Separate cache directory; `ndvi` itself is never touched |
| `ndvi_trend_masked` | continuous, 2-band | `workflow/CELESTE/prepare_ndvi_masked.R` — lake-masked counterpart of `ndvi_trend` | Separate cache directory; `ndvi_trend` itself is never touched |
| `harvest_regen` | categorical | `workflow/CELESTE/prepare_harvest_regen.R`, Ontario MNR (NIP/TUR/KEN) + NB/Irving (NBE/COC) | Genuinely DEM-dependent (rasterized onto a template from `dem_breached.tif`) — recomputed fresh whenever the terrain changes |

### NDVI lake masking — a separate, additional LOI, not a default

`ndvi_masked`/`ndvi_trend_masked` mask NHN lake/reservoir polygons out of
each group's mosaic before it's written (`mask_ndvi_lakes_celeste()` in
`prepare_ndvi_masked.R`) — a lake pixel's NDVI isn't a terrestrial
vegetation reading, and left in, it systematically drags a catchment's
summarized NDVI down in proportion to how much open water it contains.
Lakes are fetched per group against that group's own NDVI mosaic extent
(not the group's full HydroBasins AOI, which is typically far larger —
same tight-AOI reasoning as the upfront lake-conditioning feature above),
reusing `engine/prepare_lake_conditioning.R`'s `fetch_nhn_lakes_for_aoi()`
rather than a separate implementation. Unlike this project's lake-
*bisection* checks (which exclude `"Watercourse"` because a river isn't a
"lake" for that purpose), masking excludes nothing by type by default — a
watercourse is still open water and still depresses NDVI just like a lake
does.

**Deliberately a separate LOI, not an in-place toggle on `ndvi`/
`ndvi_trend`** — an earlier version masked those two LOIs directly,
overwriting their cache in place. That made a real, unrelated problem (a
small number of sites returning no data at all — see below) look like a
masking bug, since there was no unmasked baseline left to compare
against. `ndvi`/`ndvi_trend` are never masked and never touched by the
masked variant's prep functions.

**Masked at the group level, not after per-site cropping**, deliberately:
`terra::mask()` only sets values to NA within the same grid, so masking
before or after cropping to one site's catchment gives an identical final
weighted summary — masking commutes with cropping. Same for the per-pixel
Theil-Sen/Mann-Kendall trend: each pixel's slope is fit from that pixel's
own time series alone, so mask-then-fit and fit-then-mask are
byte-identical. Masking the group mosaic up front therefore costs nothing
extra and reuses the exact same `path_lazy` + per-site-crop machinery
`ndvi`/`ndvi_trend` already use — `ndvi_masked`/`ndvi_trend_masked` are
just two more named LOIs pointing at their own per-group files.

This is a semantic addition, not a change to already-cached output —
`ndvi`/`ndvi_trend`'s own cache directories are untouched; `ndvi_masked`/
`ndvi_trend_masked` write to their own `cache_dir/hydroweight_loi/
ndvi_masked/` / `.../ndvi_trend_masked/` directories.

**`MOR-CRANE`/`MOR-KEN` — resolved, not a masking issue**: an initial
investigation found these 2 sites returning no data under *any* LOI and
suspected a `hydroweight::hydroweight()`-level crash. Directly testing
just those 2 sites in isolation disproved that: `hydroweight()` succeeds
fine for both, `ndvi`/`ndvi_masked` compute correctly, and the *only*
genuinely missing LOIs are `canlcc` (a confirmed real "no overlap with
the source raster's extent" — see the canlcc note below) and
`harvest_regen` (MOR has no coverage from either harvest source at all —
already documented in `prepare_harvest_regen.R`). Both are real,
pre-existing, per-LOI gaps, not new. What actually produced the
"everything is NA" appearance in a full multi-site run was a separate,
genuinely transient issue (see next note) that happened to also hit
`ndvi` for these 2 rows on that particular run.

**A separate, real operational finding**: computing all 4 ndvi-family
LOIs together across the full 132-site batch intermittently returned no
data at all for a handful of scattered site × version combinations (5 for
CELESTE, e.g. `NBI2`/`NBR1`/`SN1UP` — 17 for CAM) — not consistently the
same sites, and not reproducible in isolation (every one of them computed
correctly the moment it was re-run on its own or in a small batch). This
matches the terra temp-file/GC-race class of issue already documented
elsewhere in this project for large multi-LOI runs, just manifesting
here as a handful of silently-dropped rows rather than a hard crash. If a
run ever produces a similarly scattered, small set of unexplained NA
rows across otherwise-unrelated sites, re-running `calculate_hydroweight_
attributes_stream()` for just that subset is the fix — not a sign of a
per-site data problem.

### canlcc (added on the same basis as CAM)

Uses the identical `CAN_LLC_2020.tif` source and 15-class scheme as both
CAM projects (see `workflow/CAM/README.md`) — same `ID`/`Class` lookup
table, duplicated inline in `run_celeste.R` rather than centralized,
matching the existing convention of each runner defining its own
`canlcc_levels`. Declared `path_lazy` (matching CAM streams' choice, not
CAM lakes' `path`) since it's a single ~2.9 GB national raster — each
site's LOI is cropped straight from the source once per site rather than
reprojecting/caching the whole file up front. Unlike CAM (working CRS
EPSG:3161), `CAN_LLC_2020.tif` is already in EPSG:3979 — CELESTE's own
working CRS (MRDEM's native CRS) — so this needs no reprojection at all.
No project-specific `prepare_*.R` script was needed: it's a single-year,
ready-to-use raster with nothing to mosaic, rasterize, or reproject ahead
of time. `workflow/CELESTE/tidy_outputs.R`'s `tidy_hydroweight_canlcc()`
mirrors CAM's version exactly (same raw column format
`canlcc_<class>_<scheme>_prop`, fixed `year = 2020L` for every row) and
writes `canlcc_long.csv` alongside the other 4 tidy tables.

Backfilled via `workflow/CELESTE/one_off/add_canlcc_loi.R`, which computes
only the `canlcc` LOI (Stage 1 config/group-manifest resolution only — no
terrain/delineation stage is touched) and joins the result onto the
existing `CELESTE_hydroweight.csv` by `(site, version)`, rather than
re-running the full multi-LOI Stage 7. Verified byte-identical
`catchment.gpkg`/`catchment_clipped.gpkg` checksums before and after.

**Confirmed real coverage gap**: `CAN_LLC_2020.tif`'s actual pixel extent
(EPSG:3979, x up to 1,619,940 m) stops short of New Brunswick/PEI, despite
its own CRS-area metadata nominally listing "New Brunswick;
... Prince Edward Island" as covered — every site in the **COC** (2/2),
**MOR** (2/2), and **NBE** (20/20) HydroBasins groups sits east of that
boundary and gets no `canlcc` value at all (24 of 132 sites, 18%). This
was caught correctly by the existing pre-crop overlap check in
`prepare_one_loi_raster()` (`workflow/R/stream/hydroweight_attributes.R`)
— those rows come out NA, not a wrong value like 0 — the same "genuine
no-coverage surfaces as NA, not silently wrong data" behavior already
established for the NDVI/NBE gap (see `CLAUDE.md`'s
`build_mosaic_vrt()` note). Every KEN/NIP/TUR site (Ontario, within the
raster's real extent) got real `canlcc` values. This is a property of the
shared `CAN_LLC_2020.tif` source file, not a bug in this run — CAM never
surfaced it because CAM's own AOI never reaches that far east.

`ndvi`/`ndvi_trend` group-tile matching uses a case-insensitive filename
match on `group_id`, with an optional literal `"Celeste_"` prefix and
optional trailing tile number — e.g. `NIP1`..`NIP5`, `KEN1`..`KEN4` — the
actual naming convention across `data/ndvi/*.tif` (verified via
`check_tile_consistency()` in `workflow/raster_attributes.R`).

`ndvi_trend` mosaics **only** a group's own source tiles, not a crop of
the full multi-province mosaic VRT (source tiles span three separate
provinces — Ontario: KEN/NIP/TUR; New Brunswick: NBE/COC; PEI: MOR — not
transcontinental). Measured total across all 6 groups: ~29 min one-time
(COC 248.5s, MOR 162.0s, TUR 226.3s, NIP 306.5s, NBE 673.7s, KEN 123.9s).
Stays single-threaded (`cores = 1`) — `cores = 8` was tested on NIP and
made it ~9.5x **worse** (2896.8s vs. 306.5s), most likely `terra::app()`'s
parallel path using a PSOCK cluster rather than fork.

`harvest_regen`'s `rasterize_competing_classes()` call OOM-killed a real
run on KEN (139.5M-cell grid, 3.5x TUR's) before the per-band-to-disk fix
described in `CLAUDE.md`; post-fix, KEN completes in ~8 min. MOR has zero
harvest/regen coverage from either source (confirmed against real gdb
features, not bbox overlap) and stays excluded — this is a pre-existing,
documented gap, not caching-related. (Separately, `cache/CELESTE/MOR/`'s
own NHN coverage was found to be independently stale — roughly half of
current on-disk NHN data — at the time of this investigation; unrelated
to `harvest_regen`, but worth knowing if MOR's numbers ever look
surprising.)

`ndvi`/`ndvi_trend` raw column names are `"ndvi_ndvi_mosaic_<N>_<stat>"`
(N = 1-based mosaic band index), **not** a calendar year — CELESTE's
`prepare_ndvi_per_group_rasters()` renames each band to `ndvi_mosaic_<N>`,
discarding the source tile's own year-bearing band name (e.g.
`"0_NDVI_1984"` → mosaic band 1). Confirmed against a raw source tile
(`data/ndvi/Landsat_NDVI_1984_2025_Celeste_Coc.tif`): band 1 =
`"0_NDVI_1984"`, band 42 = `"41_NDVI_2025"` — sequential, no gaps, no
reordering. So year = 1983 + N. `workflow/CELESTE/tidy_outputs.R`'s
`tidy_celeste_outputs()` handles this conversion; see its header for the
full column-name contract (this differs from CAM's `ndvi` column format,
which does carry the year).

## Known one-off corrections (see git history / `workflow/CELESTE/one_off/` for the scripts)

- **`backfill_streams_gpkg_nbe.R`** — NBE was delineated before an
  `engine/03_prepare_streams_burn.R` fix that makes `resolve_streams_burn()`
  persist `flowlines.gpkg` to the group cache; since that step only runs
  when `dem_burned.tif` doesn't exist yet, the fix didn't retroactively
  backfill NBE (needed for `catchment_delineation_guide.qmd`'s worked
  example). Backfills `cache/CELESTE/NBE/flowlines.gpkg`,
  `waterbodies.gpkg`, and each NBE site's `streams.gpkg`.
- **`rerun_stage7_ken_nip_tur.R`** — after the 2026-08-29 upfront
  lake-conditioning terrain fix (commit `ac5f908`) reconditioned
  KEN/NIP/TUR's DEMs and redelineated their sites, Stage 7 was never
  rerun for those 3 groups (confirmed via file timestamps —
  `CELESTE_hydroweight.csv` predated the reconditioned DEMs). Deletes and
  recomputes the per-site and group-level LOI caches that depend on the
  changed catchment geometry/DEM for KEN/NIP/TUR only, and merges the
  result into the existing combined CSV rather than overwriting it.

## Reproducing from scratch

```r
source("workflow/CELESTE/run_celeste.R")
```

Before running: confirm the MRDEM VRT, HydroBasins directory, NHN
directory, and NHN index shapefile paths near the top of the file match
your machine, and that `data/celeste_milli_sites_clean_corrected.gpkg`
exists. Then, optionally:

```r
source("workflow/CELESTE/fix_lake_bisection.R") # check for any residual lake bisection
```

To rerun a single site after a manual pour-point correction in QGIS, see
`workflow/R/engine/99_rerun_sites`'s `rerun_engine_sites()` — same pattern
described in the repo root `CLAUDE.md`.

## Migration provenance

CELESTE was migrated from a non-modular, `workflow/R/stream/`-only
pipeline (retired — see git history) onto the shared engine once the
engine was validated as an exact behavioral match: breach/D8/threshold
parameters confirmed identical by direct code comparison, ~88% of sites
reproduced pre-migration catchments near-exactly, and the remaining
differences were traced to real, explained causes (stale site
coordinates, stale cached NHN data, a since-fixed no-burn-in config
error) rather than engine bugs.
