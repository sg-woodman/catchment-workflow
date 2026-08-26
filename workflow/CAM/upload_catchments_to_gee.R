# upload_catchments_to_gee.R
# =============================================================================
# Worked examples: take CAM's catchment polygons (lake or stream, sf or
# terra), filter/group them, and upload as a Google Earth Engine asset.
#
# Uses the generic toolbox in workflow/gee_utils.R — see that file's header
# for the one-time rgee/GCS/earthengine-CLI setup, and for full parameter
# docs on every function called here. Nothing in this file runs
# automatically; every upload call below is a worked example to copy/adapt,
# same convention as gee_utils.R's own lineage
# (/Users/sam/Documents/collaborations/sommer_hydroweights/code/
# upload_gee_assets.R's "EXAMPLE USAGE" section).
#
# Source catchment files (see CLAUDE.md for what each contains):
#   output/CAM/lake_delineation/all_catchments.gpkg           (45 lake sites, clean: watershed, site_id, geom)
#   output/CAM/lake_delineation/all_catchments_clipped.gpkg   (has a messy wide X<OGF_ID>_catchment column block from upstream-lake removal — select() down to site_id/catchment_id/area_m2/geom before uploading)
#   output/CAM/stream_delineation/all_catchments.gpkg         (38 stream sites: site_id, n_cells, area_km2, flagged, geom)
#   output/CAM/stream_delineation/all_catchments_clipped.gpkg (site_id, area_km2_before, area_km2_after, n_erased, geom)
#
# Neither catchment file carries campaign/grouping metadata on its own —
# that lives in data/cam_stream_sites_raw.csv (stream sites: campaign,
# source_note) — Example 2 below joins it on before grouping, which is the
# general pattern for grouping by anything not already a column on the
# polygon layer itself.
# =============================================================================

library(sf)
library(dplyr)
library(readr)
library(here)

source(here("workflow/gee_utils.R"))

# =============================================================================
# CONFIGURATION — one-time machine setup, then session init
# =============================================================================
# One-time per machine (installs a dedicated Python + earthengine-api conda
# environment rgee calls into; NOT needed every session — uncomment only if
# rgee/GCS/the earthengine CLI aren't already set up on this machine):
# rgee::ee_install()

# One-time (or whenever the stored token expires/is revoked) — Earth
# Engine's OWN OAuth, separate from both ee_install() above and the GCS
# credentials gee_init()'s gcs argument controls. Without this, gee_init()
# fails with: "Please authorize access to your Earth Engine account by
# running earthengine authenticate ... or ee.Authenticate() in Python."
# Opens a browser for you to sign in and grant EE access — uncomment,
# run once, then leave commented again:
# rgee::ee_Authenticate(user = "samuel.woodman@gmail.com")

# GEE_USER/GEE_PROJECT below are the values already used successfully
# against this repo's own CELESTE catchments (see workflow/gee_utils.R's
# adaptation notes) — adjust if your account or Cloud project has changed.
# GEE_PROJECT must be the Cloud project's STRING ID (e.g. "gee-woodman"),
# not its numeric project number — see workflow/gee_utils.R's header for
# why that distinction matters (a numeric ID makes uploads fail with a
# misleading "asset does not exist" error).
GEE_USER <- "samuel.woodman@gmail.com"
GEE_PROJECT <- "gee-woodman"

# Run once per R session, before any upload call below. Safe/idempotent —
# unlike the actual upload calls further down, this has no side effects
# beyond authenticating (may open a browser for OAuth the first time).
#
# gcs = FALSE: Examples 1/2/4 (vector uploads via = "getInfo") need only
# Earth Engine auth, not Google Cloud Storage — a separate credential
# chain (rgee/googleCloudStorageR normally expect GCS Application Default
# Credentials, set up via `gcloud auth application-default login`, which
# most rgee setups haven't done and don't need for plain EE use). Forcing
# gcs = TRUE here is exactly what was causing a "missing GCS credentials"
# failure with nothing actually using GCS yet. Only flip this to TRUE
# once you're ready to run Example 3 (rasterized Image upload) or a
# via = "gcs_to_asset" vector upload — both genuinely need GCS, at which
# point run `gcloud auth application-default login` first.
gee_init(user = GEE_USER, project = GEE_PROJECT, gcs = FALSE)

# =============================================================================
# Example 1 — CAM lake catchments, filter by an explicit site list, upload
# as a vector FeatureCollection
# =============================================================================

lake_catchments <- st_read(
  here("output/CAM/lake_delineation/all_catchments.gpkg"),
  quiet = TRUE
)

some_lakes <- c("Chief", "Wolf", "Manitou") # TODO: your actual site_ids

lake_subset <- prepare_polygons_for_gee(
  lake_catchments,
  filter_ids = some_lakes,
  id_col = "site_id"
)
plot(st_geometry(lake_subset))

# upload_polygons_to_gee(
#   polygon = lake_subset,
#   asset_id = "projects/gee-woodman/assets/cam_lake_catchments_subset",
#   mode = "vector"
# )

# =============================================================================
# Example 2 — CAM stream catchments, grouped by campaign (metadata joined
# in from the raw sites CSV, since it isn't a column on the catchment
# layer itself), one asset uploaded per group
# =============================================================================

stream_catchments <- st_read(
  here("output/CAM/stream_delineation/all_catchments_clipped.gpkg"),
  quiet = TRUE
)

stream_meta <- read_csv(
  here("data/cam_stream_sites_raw.csv"),
  show_col_types = FALSE
) |>
  mutate(
    site_id = gsub("\\s+", "_", stream_id) |> gsub("[^A-Za-z0-9_-]", "", x = _)
  ) |>
  select(site_id, campaign, source_note)

stream_catchments_with_meta <- left_join(
  stream_catchments,
  stream_meta,
  by = "site_id"
)

campaign_groups <- prepare_polygons_for_gee(
  stream_catchments_with_meta,
  group_by = "campaign"
)
# campaign_groups is a named list: campaign_groups[["Summer"]], campaign_groups[["Fall"]]

# Optional: write each group locally first, to inspect in QGIS before upload
write_polygon_groups(
  campaign_groups,
  out_dir = here("tmp/gee_upload"),
  prefix = "cam_streams_"
)

# purrr::iwalk(campaign_groups, function(g, campaign_name) {
#   upload_polygons_to_gee(
#     polygon = g,
#     asset_id = paste0("projects/gee-woodman/assets/cam_streams_", tolower(campaign_name)),
#     mode = "vector"
#   )
# })

# =============================================================================
# Example 2b — CAM stream catchments, grouped by SPATIAL ADJACENCY instead
# of an attribute — polygons that touch/share an edge (e.g. a downstream
# site with a nested upstream site erased from it, and that upstream
# site's own catchment) get dissolved into one feature; everything else
# stays standalone. Use this when a single combined-campaign upload (as
# in Example 2) is too large/complex a polygon for your GEE script —
# grouping only genuinely-touching neighbors keeps each uploaded feature
# smaller than "everything in one campaign" while still cutting down the
# asset count vs. uploading all 38 sites individually.
# =============================================================================

adjacency_groups <- prepare_polygons_for_gee(
  stream_catchments,
  group_by_adjacency = TRUE,
  id_col = "site_id"
)
# Standalone sites are named by their own site_id (e.g. "Bell"); merged
# groups get a camelCase combined name (e.g. "acidGeorgeLumsden") — see
# group_polygons_by_adjacency()'s docstring in workflow/gee_utils.R for
# the exact naming rule, including the acronym-style-id caveat (SUD22,
# NCMN, etc. only get their very first letter lowercased, e.g.
# "sUD17SUD22" — flag if you'd rather those stay fully uppercase).

write_polygon_groups(
  adjacency_groups,
  out_dir = here("tmp/gee_upload"),
  prefix = "cam_streams_adj_"
)

purrr::iwalk(adjacency_groups, function(g, group_name) {
  upload_polygons_to_gee(
    polygon = g,
    asset_id = paste0("projects/gee-woodman/assets/cam_streams_", group_name),
    mode = "vector"
  )
})

# =============================================================================
# Example 3 — rasterize a catchment set to a binary AOI mask and upload as
# an EE Image asset (e.g. for a BAP-style compositing app that needs a
# rasterized mask rather than a vector)
# =============================================================================

# upload_polygons_to_gee(
#   polygon = stream_catchments,
#   asset_id = "projects/gee-woodman/assets/cam_streams_mask",
#   mode = "raster",
#   gcs_bucket = "your-gcs-bucket-name",
#   res = 30, # match your target sensor's pixel size, e.g. Landsat=30, Sentinel-2=10
#   field = NULL # NULL -> binary 0/1 mask; or e.g. "site_id" is character here, so rasterize a numeric column instead if you want per-site IDs burned in
# )

# =============================================================================
# Example 4 — same thing, starting from a terra SpatVector instead of sf
# (prepare_polygons_for_gee()/upload_polygons_to_gee() accept either)
# =============================================================================

# lake_catchments_terra <- terra::vect(here("output/CAM/lake_delineation/all_catchments.gpkg"))
# lake_subset_terra <- prepare_polygons_for_gee(
#   lake_catchments_terra, filter_ids = some_lakes, id_col = "site_id"
# )
# upload_polygons_to_gee(
#   polygon = lake_subset_terra,
#   asset_id = "projects/gee-woodman/assets/cam_lake_catchments_subset",
#   mode = "vector"
# )
