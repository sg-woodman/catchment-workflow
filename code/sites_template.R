# sites_template.R
# ---------------------------------------------------------------------------
# Template for defining sites as an R tibble, as an alternative to
# sites_template.csv. Use this approach when building site lists
# programmatically or when working entirely within R.
#
# Once populated, the tibble can be passed directly to validate_sites_tibble()
# (see utils.R) instead of validate_sites() which reads from a CSV file.
#
# Column descriptions:
#   site_id      : Unique site identifier. Used as the output folder name.
#                  No spaces or special characters (e.g., "L03", "YK_01").
#   site_name    : Human-readable site name. Can contain spaces.
#                  Carried through to outputs but not used in processing.
#   lon          : Longitude of pour point in decimal degrees (WGS84 / EPSG:4326).
#   lat          : Latitude of pour point in decimal degrees (WGS84 / EPSG:4326).
#   group_id     : Processing group identifier. Sites sharing a group_id are
#                  processed over a single continuous DEM extent. Sites that
#                  are hydrologically connected or geographically close should
#                  share a group_id. No spaces or special characters.
#   burn_streams : Whether to burn NHN stream network into the DEM before
#                  hydrological conditioning. Applied at the group level —
#                  all sites in a group must share this value. TRUE or FALSE.
#   aoi_buffer_m      : Optional. Buffer distance in metres around the group
#                        bounding box when defining the processing AOI.
#                        Overrides the global default (typically 1000 m).
#                        Applied at the group level. Use NA to fall back to
#                        the global default.
#   stream_threshold  : Optional. Flow accumulation threshold (in cells)
#                        used to define streams via wbt_extract_streams().
#                        Overrides the global default (typically 1000 cells).
#                        Applied at the group level. Use NA to fall back to
#                        the global default. Increase for flatter landscapes
#                        or sparse drainage networks.
# ---------------------------------------------------------------------------

library(tibble)

sites <- tribble(
  # -------------------------------------------------------------------------
  # Column names — do not modify
  ~site_id  , ~site_name , ~lat        , ~lon          , ~group_id , ~burn_streams , ~aoi_buffer_m , ~stream_threshold ,
  # -------------------------------------------------------------------------
  # Add your sites below. One row per site.
  # Examples:
  "BNWT1"   , "BNWT1"    , 67.85839    , -133.66695    , "Inuvik"  , TRUE          , NA            ,               500 ,
  "BNWT2"   , "BNWT2"    , 67.87995    , -133.62755    , "Inuvik"  , TRUE          , NA            ,               500 ,
  "CW-S-06" , "CW-S-06"  , 68.2608     , -133.26308    , "Inuvik"  , TRUE          , NA            ,              1000 ,
  "CW-S-07" , "CW-S-07"  , 68.08668    , -133.49118    , "Inuvik"  , TRUE          , NA            ,              1000 ,
  "CW-S-08" , "CW-S-08"  , 67.84236    , -133.69424    , "Inuvik"  , TRUE          , NA            ,              1000 ,
  "Km43"    , "Km43"     , 64.301      , -138.477      , "Yukon"   , TRUE          , NA            ,              1000 ,
  "Km44"    , "Km44"     , 64.31       , -138.473      , "Yukon"   , TRUE          , NA            ,              1000 ,
  "Km47"    , "Km47"     , 64.334      , -138.447      , "Yukon"   , TRUE          , NA            ,              1000 ,
  "Km69-5"  , "Km69.5"   , 64.489      , -138.206      , "Yukon"   , TRUE          , NA            ,              1000 ,
  "Km185"   , "Km185"    , 65.288      , -138.225      , "Yukon"   , TRUE          , NA            ,              1000 ,
  "PHPP05"  , "PHPP05"   , 52.39034    , -102.17277    , "Sask"    , TRUE          , NA            ,              1000 ,
  "PHPP06"  , "PHPP06"   , 52.37478    , -102.18974    , "Sask"    , TRUE          , NA            ,              1000 ,
  "U01"     , "U01"      , 48.77769603 ,  -66.13483225 , "Gaspe"   , TRUE          , NA            ,              1000 ,
  "U02"     , "U02"      , 48.84714641 ,  -65.95554781 , "Gaspe"   , TRUE          , NA            ,              1000 ,
  "U03"     , "U03"      , 48.87196147 ,  -65.82999172 , "Gaspe"   , TRUE          , NA            ,              1000 ,
  "C04"     , "C04"      , 48.52289617 ,  -66.10337881 , "Gaspe"   , TRUE          , NA            ,              1000 ,
  "C05"     , "C05"      , 48.53412982 ,  -65.8718441  , "Gaspe"   , TRUE          , NA            ,              1000 ,
  "C06"     , "C06"      , 48.55861841 ,  -65.78111336 , "Gaspe"   , TRUE          , NA            ,              1000 ,
  "C07"     , "C07"      , 48.58631372 ,  -65.76315142 , "Gaspe"   , TRUE          , NA            ,              1000 ,
  "L08"     , "L08"      , 48.40597698 ,  -65.61372834 , "Gaspe"   , TRUE          , NA            ,              1000 ,
  "L09"     , "L09"      , 48.30161709 ,  -65.7212379  , "Gaspe"   , TRUE          , NA            ,              1000 ,
  "L10"     , "L10"      , 48.28311055 ,  -65.67445472 , "Gaspe"   , TRUE          , NA            ,              1000 ,
  "L11"     , "L11"      , 48.31098102 ,  -65.54432336 , "Gaspe"   , TRUE          , NA            ,              1000 ,
  "L12"     , "L12"      , 48.31115015 ,  -65.54230006 , "Gaspe"   , TRUE          , NA            ,              1000 ,
  "HB1"     , "HB1"      , 49.088413   ,  -57.901861   , "NFLD"    , TRUE          , NA            ,              1000 ,
  "PA1"     , "PA1"      , 49.007639   ,  -57.566861   , "NFLD"    , TRUE          , NA            ,              1000 ,
  "LS1"     , "LS1"      , 48.859666   ,  -57.946887   , "NFLD"    , TRUE          , NA            ,              1000
  # -------------------------------------------------------------------------
)


# Preview the sites tibble
print(sites)
