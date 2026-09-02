# gee_utils.R
# ---------------------------------------------------------------------------
# Standalone, project-agnostic toolbox for getting polygon (or rasterized
# polygon) data into Google Earth Engine as an asset — same positioning as
# workflow/raster_attributes.R: not sourced into either project's Stage
# sequence by default, sourced ad hoc when you need it.
#
# Adapted from two sources:
#   - /Users/sam/Documents/collaborations/sommer_hydroweights/code/
#     upload_gee_assets.R (internally titled gee_polygon_upload.R) — the
#     upload mechanics (gee_init(), upload_polygon_vector(),
#     rasterize_polygon(), upload_raster_to_gee(), upload_polygons_to_gee())
#     are carried over close to verbatim, generalized to accept a terra
#     SpatVector as well as an sf object (the original only accepted sf).
#   - /Users/sam/Documents/cfs/catchment_delineation/code/
#     export_shp_for_gee.R — that script's actual logic (filter a
#     catchment layer down to a handful of hardcoded per-region site-code
#     lists, write each as its own shapefile) is generalized here into
#     prepare_polygons_for_gee() (filter/group on ANY id/group column, not
#     hardcoded site lists) + write_polygon_groups() (write out however
#     many groups result, any vector driver) — the write-to-local-file step
#     is now optional, not required: rgee::sf_as_ee()'s default "getInfo"
#     upload path uploads directly from an in-memory sf object, no
#     shapefile roundtrip needed for small-to-moderate feature counts.
#
# Purpose
# -------
# Three composable steps:
#   1. prepare_polygons_for_gee() — filter and/or split a polygon layer
#      (sf or terra SpatVector) into the group(s) you actually want to
#      upload, by any id/group column already on it (join on whatever
#      project-specific metadata you need BEFORE calling this — see a
#      calling project's own upload_catchments_to_gee.R for a worked example).
#   2. write_polygon_groups() — OPTIONAL. Write the result(s) of step 1 to
#      local vector files, if you want to inspect/QA a group before
#      uploading, or need shapefiles for something else.
#   3. upload_polygons_to_gee() — upload one polygon layer (sf, SpatVector,
#      or a file path) to GEE, as a vector FeatureCollection or (rasterized
#      first) an Image asset.
#
# Dependencies
# ------------
#   R packages : rgee, sf, terra, jsonlite, purrr
#   External   : Google Cloud SDK (`gsutil`) and the Earth Engine CLI
#                (`earthengine`) must be on PATH — both are installed
#                automatically the first time you run rgee::ee_install(),
#                which sets up a dedicated Python + earthengine-api conda
#                environment that rgee calls into.
#
# One-time setup (run once per machine, NOT every session)
# ---------------------------------------------------------------------------
#   install.packages(c("rgee", "sf", "terra", "jsonlite", "purrr"))
#   rgee::ee_install()
#
# Per-session init — gee_init(user, project) below (gcs = FALSE by
# default; only needed for a raster upload or via = "gcs_to_asset", see
# gee_init()'s own docstring):
#   gee_init(user = "your_email@gmail.com", project = "your-cloud-project-id")
#
# `project` must be the Cloud project's string ID (e.g. "gee-woodman"), not
# its numeric project number. Read the string ID off any existing asset
# path: projects/<this part>/assets/.... Earth Engine's asset-path
# validation silently mangles a numeric project ID into a broken
# "earthengine-legacy" path on export/create calls — reads (getAsset,
# getAssetRoots, listing) will look fine, but uploads fail with a confusing
# "asset does not exist" error.
#
# Raster uploads specifically require write access to a Google Cloud
# Storage (GCS) bucket that your GEE-linked Cloud project can also read
# from — GEE ingests images from GCS, not directly from your local disk.
# Vector uploads do not need a bucket for small-to-moderate polygon counts
# (rgee stages these via the EE REST "getInfo" upload path automatically).
# ---------------------------------------------------------------------------

# =============================================================================
# 1. INITIALIZATION
# =============================================================================

#' Initialize the Earth Engine session
#'
#' Thin wrapper around rgee::ee_Initialize() so the rest of the script has a
#' single, obvious entry point.
#'
#' `project` must be the Cloud project's string ID (e.g. "gee-woodman"), not
#' its numeric project number — see the note at the top of this file. If
#' omitted, ee_Initialize() falls back to project = NULL, which makes
#' vector and raster asset uploads fail server-side with a misleading
#' "asset does not exist" error, even though read-only calls
#' (getAssetRoots, listing) succeed.
#'
#' `gcs` defaults to FALSE, not TRUE — GCS auth is a separate credential
#' chain from Earth Engine's own OAuth (rgee/googleCloudStorageR normally
#' expect GCS Application Default Credentials, set up via `gcloud auth
#' application-default login`), and most of what this toolbox does
#' (upload_polygon_vector() via = "getInfo", the default) never touches
#' GCS at all. Forcing it on unconditionally produces a "missing GCS
#' credentials" failure for a plain EE-only session that was never going
#' to need it — confirmed directly. Only pass gcs = TRUE once you're
#' actually about to rasterize+upload an Image asset (upload_raster_to_gee())
#' or use via = "gcs_to_asset" for a vector upload.
#'
#' @param user Character. Your GEE-linked Google account email.
#' @param project Character. Cloud project string ID that assets are
#'   uploaded to, e.g. "gee-woodman".
#' @param gcs Logical. Whether to also initialize GCS credentials. Default
#'   FALSE — see above.
#' @return Invisibly returns TRUE on success.
gee_init <- function(user = NULL, project, gcs = FALSE) {
  rgee::ee_Initialize(user = user, project = project, gcs = gcs)
  invisible(TRUE)
}

# =============================================================================
# 2. FILTER / GROUP (generalizes catchment_delineation's export_shp_for_gee.R)
# =============================================================================

#' Filter and/or split a polygon layer into named groups for GEE upload
#'
#' Accepts either an sf object or a terra SpatVector — returns the same
#' class it was given. Generalizes export_shp_for_gee.R's pattern (filter a
#' catchment layer down to several hardcoded per-region site-code lists)
#' into filtering/grouping by any column already on the layer, instead of
#' a fixed set of site lists baked into the script. To group by metadata
#' that ISN'T already a column on your polygon layer (e.g. a campaign or
#' watershed label that lives in a separate sites table), join it on
#' first — see a calling project's own upload_catchments_to_gee.R for a
#' worked example — this function only ever operates on columns already
#' present.
#'
#' @param polygons  sf object or terra SpatVector, with an identifying
#'   column (e.g. "site_id") and, if grouping by attribute, a grouping
#'   column.
#' @param filter_ids Character vector of ids to keep (matched against
#'   id_col), or NULL (default) to keep everything. Applied before either
#'   form of grouping below.
#' @param id_col    Character. Column filter_ids is matched against, and
#'   (for group_by_adjacency) the column standalone/combined names are
#'   built from. Default "site_id".
#' @param group_by  Character or NULL (default). Column name to split by
#'   VALUE. If supplied, returns a NAMED LIST of sf/SpatVector objects, one
#'   per distinct value of that column (list names are that value, coerced
#'   to character). Mutually exclusive with group_by_adjacency.
#' @param group_by_adjacency Logical. Default FALSE. If TRUE, splits by
#'   SPATIAL adjacency instead of an attribute value — groups polygons
#'   that touch/share an edge into one dissolved feature per connected
#'   component, leaving polygons with no touching neighbor standalone. See
#'   group_polygons_by_adjacency() for the full contract (naming,
#'   snap_tolerance_m). Mutually exclusive with group_by.
#' @param snap_tolerance_m Numeric. Only used when group_by_adjacency =
#'   TRUE — passed through to group_polygons_by_adjacency(). Default 1.
#'
#' @return An sf/SpatVector object (group_by = NULL and
#'   group_by_adjacency = FALSE) or a named list of them (either grouping
#'   mode) — same class as the input either way.
prepare_polygons_for_gee <- function(
  polygons, filter_ids = NULL, id_col = "site_id", group_by = NULL,
  group_by_adjacency = FALSE, snap_tolerance_m = 1
) {
  is_terra <- inherits(polygons, "SpatVector")
  if (!is_terra && !inherits(polygons, "sf")) {
    stop("polygons must be an sf object or a terra SpatVector.")
  }
  if (!is.null(group_by) && group_by_adjacency) {
    stop("group_by and group_by_adjacency are mutually exclusive — pick one grouping mode.")
  }

  get_col <- function(x, col) {
    if (!col %in% names(x)) {
      stop(glue::glue("Column '{col}' not found on polygons. Available: {paste(names(x), collapse = ', ')}"))
    }
    if (is_terra) terra::values(x)[[col]] else sf::st_drop_geometry(x)[[col]]
  }

  if (!is.null(filter_ids)) {
    keep <- get_col(polygons, id_col) %in% filter_ids
    n_before <- nrow(polygons)
    polygons <- polygons[keep, ]
    message(glue::glue("Filtered {n_before} -> {nrow(polygons)} feature(s) by {id_col}."))
  }

  if (group_by_adjacency) {
    return(group_polygons_by_adjacency(polygons, id_col = id_col, snap_tolerance_m = snap_tolerance_m))
  }

  if (is.null(group_by)) {
    return(polygons)
  }

  group_vals <- unique(get_col(polygons, group_by))
  groups <- purrr::map(group_vals, function(v) {
    polygons[get_col(polygons, group_by) == v, ]
  })
  names(groups) <- as.character(group_vals)

  message(glue::glue(
    "Split into {length(groups)} group(s) by '{group_by}': {paste(names(groups), collapse = ', ')}"
  ))

  groups
}

#' Group polygons by spatial adjacency (touching / sharing an edge)
#'
#' Alternative to attribute-based grouping — instead of splitting by a
#' column value, finds connected components of polygons that touch (share
#' a boundary/edge; sf::st_touches()'s definition — interiors don't
#' overlap, so this is genuine adjacency, not "nearby" or "overlapping")
#' and dissolves each component into one feature. A standalone polygon
#' (no touching neighbor) passes through unchanged, named by its own
#' id_col value. A component of 2+ touching polygons is unioned into a
#' single feature, named by combining their id_col values — sorted
#' alphabetically first (so the name is stable regardless of input row
#' order), camelCase, no separators: "Chief" + "Wolf" -> "chiefWolf".
#' Non-alphanumeric characters (spaces, underscores) in an id are stripped
#' before combining, e.g. "Sam_Martin" -> "SamMartin" mid-name.
#'
#' Called by prepare_polygons_for_gee(group_by_adjacency = TRUE); also
#' fine to call directly.
#'
#' @param polygons sf object or terra SpatVector.
#' @param id_col   Character. Column identifying each feature. Default
#'   "site_id".
#' @param snap_tolerance_m Numeric. Buffer distance (metres, in the
#'   polygons' own CRS — must be projected/metric) applied to a COPY of
#'   the geometry before testing adjacency, to treat near-touching
#'   polygons (separated only by rasterization/floating-point noise, not
#'   a real gap) as touching. Does not affect the geometry actually
#'   written to the output — only the adjacency test. Set to 0 for exact-
#'   touching only. Default 1.
#'
#' @return A named list of sf/SpatVector objects (same class as input),
#'   one per connected component. A merged feature carries only id_col
#'   (the combined name) and `adjacency_members` (comma-separated list of
#'   the original ids it replaces) — its other attribute columns are
#'   dropped rather than carried over from one arbitrary member, since
#'   they'd no longer describe the unioned geometry (e.g. an area column
#'   computed for one member's footprint, not the merged one's). A
#'   standalone feature keeps all of its original columns unchanged.
group_polygons_by_adjacency <- function(polygons, id_col = "site_id", snap_tolerance_m = 1) {
  is_terra <- inherits(polygons, "SpatVector")
  sf_polygons <- if (is_terra) sf::st_as_sf(polygons) else polygons
  if (!inherits(sf_polygons, "sf")) {
    stop("polygons must be an sf object or a terra SpatVector.")
  }
  if (!id_col %in% names(sf_polygons)) {
    stop(glue::glue("Column '{id_col}' not found on polygons. Available: {paste(names(sf_polygons), collapse = ', ')}"))
  }

  n <- nrow(sf_polygons)
  # st_make_valid() defensively — confirmed directly against real
  # catchment output that a clipped (upstream-erased) polygon can come out
  # geometrically invalid (self-intersection from the erase operation),
  # which risks unpredictable GEOS topology results downstream.
  test_geom <- sf::st_geometry(sf::st_make_valid(sf_polygons))

  if (snap_tolerance_m > 0) {
    # Buffering turns "touching" into "overlapping" — two polygons that
    # exactly share an edge no longer merely TOUCH once either side is
    # expanded into the other's territory, so st_touches() would (and,
    # confirmed directly, does) go FALSE the moment any buffer is applied,
    # the opposite of what a tolerance is supposed to do. st_intersects()
    # is the correct predicate once buffered — but its single-argument
    # form includes each geometry as its own neighbor (a geometry always
    # intersects itself), which st_touches() never does, so self-matches
    # are stripped explicitly to keep the same adjacency-list shape either
    # way.
    test_geom <- sf::st_buffer(test_geom, snap_tolerance_m)
    touch_list <- sf::st_intersects(test_geom)
    touch_list <- purrr::imap(touch_list, function(nbrs, i) setdiff(nbrs, i))
  } else {
    touch_list <- sf::st_touches(test_geom)
  }

  comp_id <- gee_connected_components(touch_list)
  ids <- as.character(sf_polygons[[id_col]])
  comp_ids <- sort(unique(comp_id))
  n_merged <- sum(tabulate(comp_id) > 1)

  groups <- purrr::map(comp_ids, function(cid) {
    idx <- which(comp_id == cid)
    member_ids <- sort(ids[idx])

    if (length(idx) == 1) {
      return(sf_polygons[idx, ])
    }

    dissolved_geom <- sf::st_union(sf::st_geometry(sf::st_make_valid(sf_polygons[idx, ])))
    combined_name <- gee_to_camel_case(member_ids)
    new_row <- tibble::tibble(!!id_col := combined_name, adjacency_members = paste(member_ids, collapse = ", "))
    sf::st_sf(new_row, geometry = dissolved_geom)
  })

  group_names <- purrr::map_chr(comp_ids, function(cid) {
    member_ids <- sort(ids[comp_id == cid])
    if (length(member_ids) == 1) member_ids else gee_to_camel_case(member_ids)
  })
  names(groups) <- group_names

  if (is_terra) {
    groups <- purrr::map(groups, terra::vect)
  }

  message(glue::glue(
    "Adjacency grouping: {n} feature(s) -> {length(groups)} group(s) ",
    "({n_merged} merged from touching neighbors, {length(groups) - n_merged} standalone)."
  ))

  groups
}

#' Connected components of an adjacency list (BFS) — no igraph dependency
#'
#' @param adj_list List, one element per node, each an integer vector of
#'   neighbor indices (e.g. sf::st_touches()'s return shape).
#' @return Integer vector, one component id per node (1-indexed, in
#'   discovery order).
gee_connected_components <- function(adj_list) {
  n <- length(adj_list)
  visited <- rep(FALSE, n)
  comp_id <- integer(n)
  current_comp <- 0L

  for (i in seq_len(n)) {
    if (visited[i]) next
    current_comp <- current_comp + 1L
    queue <- i
    visited[i] <- TRUE
    while (length(queue) > 0) {
      node <- queue[1]
      queue <- queue[-1]
      comp_id[node] <- current_comp
      nbrs <- adj_list[[node]]
      new_nodes <- nbrs[!visited[nbrs]]
      if (length(new_nodes) > 0) {
        visited[new_nodes] <- TRUE
        queue <- c(queue, new_nodes)
      }
    }
  }

  comp_id
}

#' Combine a set of ids into one camelCase name
#'
#' First id lowercased at its first letter; every subsequent id
#' uppercased at its first letter; non-alphanumeric characters stripped
#' from every id first (so "Sam_Martin" contributes "SamMartin", not
#' "Sam_Martin"); no separators between parts. Caller is responsible for
#' sorting `ids` into a stable order first if determinism across runs
#' matters (group_polygons_by_adjacency() sorts alphabetically).
#'
#' @param ids Character vector, in the order to combine them.
#' @return Character, length 1.
gee_to_camel_case <- function(ids) {
  clean <- gsub("[^A-Za-z0-9]", "", ids)
  parts <- purrr::imap_chr(clean, function(id, i) {
    if (nchar(id) == 0) {
      return(id)
    }
    first <- if (i == 1) tolower(substr(id, 1, 1)) else toupper(substr(id, 1, 1))
    paste0(first, substr(id, 2, nchar(id)))
  })
  paste(parts, collapse = "")
}

#' Write one or more polygon groups to local vector files
#'
#' Optional — you don't need this before uploading (upload_polygons_to_gee()
#' accepts an sf object or SpatVector directly), but useful for inspecting/
#' QA-ing a group in QGIS first, or if you need shapefiles for something
#' else. Generalizes export_shp_for_gee.R's per-region
#' writeVector(..., here("tmp/<region>_catch.shp")) calls.
#'
#' @param groups  A single sf/SpatVector object, or a named list of them
#'   (e.g. from prepare_polygons_for_gee()'s group_by path).
#' @param out_dir Character. Destination directory (created if needed).
#' @param prefix  Character. Filename prefix; each group is written as
#'   "<prefix><group_name>.<ext>". If `groups` isn't a named list, `prefix`
#'   alone (no suffix) is used as the filename stem.
#' @param driver  Character. "shp" or "gpkg". Default "gpkg" (fewer file-
#'   count/field-name-length gotchas than shapefile; export_shp_for_gee.R
#'   used "shp" only because that's what its era's workflow expected).
#'
#' @return Character vector of paths written, invisibly.
write_polygon_groups <- function(groups, out_dir, prefix = "", driver = c("gpkg", "shp")) {
  driver <- match.arg(driver)
  ext <- if (driver == "gpkg") "gpkg" else "shp"
  fs::dir_create(out_dir, recurse = TRUE)

  is_list_of_groups <- is.list(groups) && !inherits(groups, "sf") && !inherits(groups, "SpatVector")

  if (!is_list_of_groups) {
    groups <- setNames(list(groups), prefix)
    prefix <- ""
  }

  paths <- purrr::imap_chr(groups, function(g, name) {
    out_path <- fs::path(out_dir, paste0(prefix, name, ".", ext))
    if (inherits(g, "SpatVector")) {
      terra::writeVector(g, out_path, overwrite = TRUE)
    } else {
      sf::st_write(g, out_path, delete_dsn = TRUE, quiet = TRUE)
    }
    message(glue::glue("Wrote {nrow(g)} feature(s) -> {out_path}"))
    out_path
  })

  invisible(unname(paths))
}

# =============================================================================
# 3. VECTOR UPLOAD (polygon -> EE FeatureCollection asset)
# =============================================================================

#' Upload a polygon object to Earth Engine as a FeatureCollection asset
#'
#' @param polygon sf object or terra SpatVector (POLYGON / MULTIPOLYGON) to
#'   upload.
#' @param asset_id Character. Full destination asset path, e.g.
#'   "projects/gee-woodman/assets/my_polygon" (string project ID, not the
#'   numeric project number). The parent folder must already exist in your
#'   EE asset tree.
#' @param via Character. Upload path used internally by rgee::sf_as_ee().
#'   "getInfo" works for small-to-moderate feature collections (no bucket
#'   needed); "gcs_to_asset" routes through GCS and scales better for large
#'   or complex geometries.
#' @param gcs_bucket Character. Required only when via = "gcs_to_asset".
#' @param overwrite Logical. Overwrite an existing asset at asset_id.
#' @param monitor Logical. Block and print progress until the export task
#'   finishes.
#'
#' @return The rgee Task object (invisibly).
upload_polygon_vector <- function(
  polygon, asset_id, via = c("getInfo", "gcs_to_asset"),
  gcs_bucket = NULL, overwrite = FALSE, monitor = TRUE
) {
  via <- match.arg(via)

  if (inherits(polygon, "SpatVector")) {
    polygon <- sf::st_as_sf(polygon)
  } else if (!inherits(polygon, "sf")) {
    stop("polygon must be an sf object or a terra SpatVector. Use sf::st_read()/terra::vect() first if you have a file path.")
  }

  if (via == "gcs_to_asset" && is.null(gcs_bucket)) {
    stop("gcs_bucket is required when via = 'gcs_to_asset'.")
  }

  # sf_as_ee() converts the sf object into an ee$FeatureCollection, staging
  # it via the chosen method under the hood.
  ee_fc <- rgee::sf_as_ee(
    x = polygon, via = via,
    bucket = if (via == "gcs_to_asset") gcs_bucket else NULL
  )

  # ee_table_to_asset() persists that FeatureCollection as a permanent asset
  # (sf_as_ee alone only creates a transient/temp EE object).
  task <- rgee::ee_table_to_asset(
    collection = ee_fc, description = paste0("upload_", basename(asset_id)),
    assetId = asset_id, overwrite = overwrite
  )

  task$start()

  if (monitor) {
    rgee::ee_monitoring(task) # blocks and prints task status until done
  } else {
    message("Task started. Check progress with rgee::ee_monitoring(task) or the EE Task Manager.")
  }

  invisible(task)
}

# =============================================================================
# 4. RASTERIZATION (polygon -> local GeoTIFF)
# =============================================================================

#' Rasterize a polygon object to a local GeoTIFF
#'
#' @param polygon sf object or terra SpatVector to rasterize.
#' @param template Optional SpatRaster defining output extent/resolution/
#'   CRS. If NULL, a template is built automatically from the polygon
#'   extent at `res`, in the polygon's own CRS.
#' @param res Numeric. Output pixel size (map units of the CRS), used only
#'   when template is NULL.
#' @param field Character or NULL. Attribute column to burn into the raster
#'   (e.g. a numeric class/ID field). If NULL, all polygon cells are burned
#'   with value 1 — the typical case for a binary AOI mask.
#' @param background Numeric. Value assigned to cells outside the
#'   polygon(s).
#' @param touches Logical. If TRUE, any pixel touched by the polygon
#'   boundary is included (useful for AOI masks); if FALSE, only pixel
#'   centers inside the polygon are included.
#' @param datatype Character. terra/GDAL datatype for the output, e.g.
#'   "INT1U" for a 0/1 mask, "INT2S" for small signed integers.
#' @param out_path Character. Destination path for the GeoTIFF.
#'
#' @return Character path to the written GeoTIFF, invisibly.
rasterize_polygon <- function(
  polygon, template = NULL, res = 30, field = NULL, background = 0,
  touches = TRUE, datatype = "INT1U", out_path = tempfile(fileext = ".tif")
) {
  if (!inherits(polygon, "sf") && !inherits(polygon, "SpatVector")) {
    stop("polygon must be an sf object or a terra SpatVector. Use sf::st_read()/terra::vect() first if you have a file path.")
  }

  polygon_vect <- terra::vect(polygon)

  if (is.null(template)) {
    template <- terra::rast(polygon_vect, resolution = res)
  }

  # terra::rasterize()'s field = NULL does NOT mean "burn 1" in the
  # installed terra version (1.9.27) — confirmed directly: it burns
  # `background` for every cell instead (silently producing an all-zero
  # raster, no error). The documented-working way to get a binary 1/0
  # mask is field = 1 (a literal constant value, not a column reference).
  rasterize_field <- if (is.null(field)) 1 else field

  r <- terra::rasterize(
    x = polygon_vect, y = template,
    field = rasterize_field, # NULL `field` arg -> literal 1, burns all polygon cells; a column name burns that column's values instead
    background = background, touches = touches
  )

  terra::writeRaster(
    r, filename = out_path, overwrite = TRUE,
    datatype = datatype, gdal = c("COMPRESS=DEFLATE")
  )

  message("Rasterized polygon written to: ", out_path)
  invisible(out_path)
}

# =============================================================================
# 5. RASTER UPLOAD (local GeoTIFF -> GCS -> EE Image asset)
# =============================================================================

#' Upload a local GeoTIFF to Earth Engine as an Image asset
#'
#' GEE cannot ingest images directly from local disk — the file must first
#' be staged in a GCS bucket, then registered with EE via an ingestion
#' manifest submitted through the `earthengine` CLI.
#'
#' @param tif_path Character. Path to the local GeoTIFF (e.g. output of
#'   rasterize_polygon()).
#' @param asset_id Character. Full destination asset path, e.g.
#'   "projects/gee-woodman/assets/my_raster" (string project ID, not the
#'   numeric project number). The parent folder must already exist in your
#'   EE asset tree.
#' @param gcs_bucket Character. Name of a GCS bucket you can write to and
#'   that your EE Cloud project can read from.
#' @param gcs_path Character or NULL. Object path within the bucket. If
#'   NULL, defaults to "gee_uploads/<filename>".
#' @param pyramiding_policy Character. One of "MEAN", "MODE", "MIN", "MAX",
#'   "SAMPLE". Use "MODE" for categorical/mask rasters (typical for a
#'   rasterized AOI), "MEAN" for continuous data.
#' @param nodata_value Numeric. Value in the GeoTIFF to treat as missing
#'   data.
#' @param wait Logical. If TRUE, poll the earthengine CLI until the
#'   ingestion task finishes and report success/failure.
#'
#' @return Character EE task ID, invisibly.
upload_raster_to_gee <- function(
  tif_path, asset_id, gcs_bucket, gcs_path = NULL,
  pyramiding_policy = c("MODE", "MEAN", "MIN", "MAX", "SAMPLE"),
  nodata_value = 0, wait = TRUE
) {
  pyramiding_policy <- match.arg(pyramiding_policy)

  if (!file.exists(tif_path)) {
    stop("tif_path does not exist: ", tif_path)
  }

  if (is.null(gcs_path)) {
    gcs_path <- paste0("gee_uploads/", basename(tif_path))
  }
  gcs_uri <- paste0("gs://", gcs_bucket, "/", gcs_path)

  # --- Step 1: stage the GeoTIFF in GCS ---------------------------------
  message("Uploading to GCS: ", gcs_uri)
  gsutil_status <- system2(
    "gsutil", args = c("cp", shQuote(tif_path), shQuote(gcs_uri)),
    stdout = TRUE, stderr = TRUE
  )
  message(paste(gsutil_status, collapse = "\n"))

  # --- Step 2: build the EE ingestion manifest ---------------------------
  # See: https://developers.google.com/earth-engine/guides/image_upload
  manifest <- list(
    name = asset_id,
    tilesets = list(list(sources = list(list(uris = list(gcs_uri))))),
    pyramiding_policy = pyramiding_policy,
    missing_data = list(values = list(nodata_value))
  )

  manifest_path <- tempfile(fileext = ".json")
  jsonlite::write_json(manifest, manifest_path, auto_unbox = TRUE)

  # --- Step 3: kick off ingestion via the earthengine CLI -----------------
  message("Submitting ingestion task to Earth Engine...")
  upload_result <- system2(
    "earthengine", args = c("upload", "image", "--manifest", shQuote(manifest_path)),
    stdout = TRUE, stderr = TRUE
  )
  message(paste(upload_result, collapse = "\n"))

  # CLI output typically looks like: "Started upload task with ID: XXXXXXXX"
  task_id_line <- grep("task with ID", upload_result, value = TRUE)
  task_id <- if (length(task_id_line) > 0) {
    sub(".*ID:\\s*", "", task_id_line[1])
  } else {
    NA_character_
  }

  if (wait && !is.na(task_id)) {
    message("Waiting on ingestion task ", task_id, " ...")
    system2("earthengine", args = c("task", "wait", task_id))
    message("Ingestion complete: ", asset_id)
  } else if (wait) {
    message("Could not parse task ID from CLI output; check status manually with `earthengine task list`.")
  }

  invisible(task_id)
}

# =============================================================================
# 6. WRAPPER: single entry point for either mode
# =============================================================================

#' Upload polygon(s) to Earth Engine, either as a vector or a rasterized
#' image
#'
#' @param polygon sf object, terra SpatVector, OR a file path (shapefile/
#'   GeoPackage/etc, read via sf::st_read()).
#' @param asset_id Character. Destination EE asset path.
#' @param mode Character. "vector" for a FeatureCollection asset, "raster"
#'   for a rasterized Image asset.
#' @param gcs_bucket Character. Required for mode = "raster" (and for
#'   mode = "vector" only if vector_via = "gcs_to_asset").
#' @param res Numeric. Pixel size for rasterization (mode = "raster" only).
#' @param field Character or NULL. Attribute to burn when rasterizing
#'   (mode = "raster" only); NULL burns a binary 1/0 mask.
#' @param vector_via Character. Passed through to upload_polygon_vector()
#'   as `via` (mode = "vector" only).
#' @param overwrite Logical. Overwrite an existing asset (mode = "vector"
#'   only; raster ingestion will fail if the asset already exists).
#' @param monitor Logical. Block and print progress until the vector
#'   export task finishes (mode = "vector" only).
#'
#' @return The result of the underlying upload function (Task object for
#'   vector, task ID string for raster), invisibly.
upload_polygons_to_gee <- function(
  polygon, asset_id, mode = c("vector", "raster"), gcs_bucket = NULL,
  res = 30, field = NULL, vector_via = "getInfo", overwrite = FALSE,
  monitor = TRUE
) {
  mode <- match.arg(mode)

  # Accept an sf object, a SpatVector, or a path to read
  if (is.character(polygon)) {
    polygon <- sf::st_read(polygon, quiet = TRUE)
  }

  if (mode == "vector") {
    result <- upload_polygon_vector(
      polygon = polygon, asset_id = asset_id, via = vector_via,
      gcs_bucket = gcs_bucket, overwrite = overwrite, monitor = monitor
    )
  } else {
    tif_path <- rasterize_polygon(polygon = polygon, res = res, field = field)
    result <- upload_raster_to_gee(
      tif_path = tif_path, asset_id = asset_id, gcs_bucket = gcs_bucket
    )
  }

  invisible(result)
}
