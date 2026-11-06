
# ==============================================================
# FIGURE 4 — INDEPENDENT, REVISED, COMPLETE SCRIPT
#
# Complementarity between Nutrient Loading Patterns
# and Information-Network Structure
#
# a/b: CAST loading context + F-Net PCS/backbone
# c:   Node importance vs NTN/CAST loading
# d:   Discordant PCS/NTN monitoring priorities
# e:   Directional TE across CAST loading contexts
# f:   Candidate complementary monitoring priorities
#
# Additional:
# - CAST unique-region correlation
# - PCS vs in/out-strength
# - Common-sample top-20% overlap
# - Wilson descriptive intervals
# - TE threshold/coverage QA
# - SPARROW 2012 supplemental comparison
#
# NO old Fig4_* intermediate files required,
# except the explicitly recomputed PCS.
# ==============================================================


# ==============================================================
# 0. PACKAGES
# ==============================================================

PKGS <- c(
  "readr", "dplyr", "tidyr", "tibble",
  "ggplot2", "sf", "patchwork",
  "scales", "stringr"
)

for (p in PKGS) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

invisible(
  lapply(PKGS, library, character.only = TRUE)
)


# ==============================================================
# 1. PATHS — EDIT ROOT ONLY IF NEEDED
# ==============================================================

source("config.R")
ROOT <- get_data_root()

OUT <- file.path(
  ROOT, "Fig4_final_original_style_yield_CASTDM"
)

dir.create(
  OUT, recursive = TRUE, showWarnings = FALSE
)

PCS_FILE <- file.path(
  ROOT,
  "Fig4_PCS_recomputed",
  "Global_PCS_TN_TP_recomputed.csv"
)

FUSION_DIR <- file.path(
  ROOT, "Fusion_DualLag_2011_2023"
)

NTN_FILE <- file.path(
  ROOT, "NTN_2023_MonthlyLoadTable.csv"
)

# Verified metric export converted from official "CAST_DM, Mlbs".
CASTDM_FILE <- file.path(ROOT, "CAST_DM_NTN_annual.csv")
# Optional verified full upstream area at EACH monitoring station (km2).
# If absent, area is estimated from ntnwatersheds.shp and flagged for QA.
AREA_FILE <- file.path(ROOT, "NTN_watershed_areas.csv")

STATION_FILE <- file.path(
  ROOT, "NTN-station.txt"
)

RIVER_FILE <- file.path(
  ROOT, "水系", "水系.shp"
)

BASIN_FILE <- file.path(
  ROOT, "ntnwatersheds",
  "ntnwatersheds.shp"
)

CAST_FILE <- file.path(
  ROOT,
  "22_STEP4C_CAST",
  "00_CAST_RAW",
  "CAST_92segments_2006_2012_2019_2024_TNTP_long.csv"
)

CAST_SHP <- file.path(
  ROOT,
  "CAST_Targeting",
  "P6Beta_v3_LRSegs_081516.shp"
)

TE_FILES <- c(
  TN = file.path(
    ROOT,
    "TE_network_results",
    "TN",
    "main_2011_2023_all_edges.csv"
  ),
  TP = file.path(
    ROOT,
    "TE_network_results",
    "TP",
    "main_2011_2023_all_edges.csv"
  )
)

SP_DIR <- file.path(
  ROOT, "17_STEP4B2_SPARROW64",
  "00_USGS_RAW"
)

SP_INPUT <- file.path(
  SP_DIR, "input_clean.csv"
)

SP_FILES <- c(
  TN = file.path(SP_DIR, "predict_rnd_N_f.csv"),
  TP = file.path(SP_DIR, "predict_rnd_P_f.csv")
)

# If a river-topology-verified crosswalk exists,
# set its path here:
#
# Columns: STAID, comid_stack
#
SP_CROSSWALK <- NA_character_


# ==============================================================
# 2. SETTINGS
# ==============================================================

MAP_CRS <- 5070

CAST_YEAR <- 2019
CAST_SCENARIO <- "2019 Progress"

SP_MODELSPACE <- 2012

NTN_START_YEAR <- 2012
NTN_END_YEAR <- 2022

# IMPORTANT:
# This code assumes KLoad is a MONTHLY INTEGRATED LOAD.
# Confirm units and time aggregation in NTN metadata.

NTN_LOAD_FIELD <- "KLoad"
# Provisional numeric-consistency result, not primary metadata confirmation.
# 1000 lb per unit of NTN KLoad (2012 monthly Q*C check mean ratio 1.0038).
NTN_KLOAD_UNIT <- "1000lb"
NTN_KG_FACTOR <- switch(NTN_KLOAD_UNIT,
                        "kg"=1, "lb"=0.45359237,
                        "1000lb"=453.59237,
                        stop("Verify NTN KLoad units"))

MIN_MONTHS <- 12L
MIN_YEARS <- 5L

HOT_FRAC <- 0.20

TE_P <- 0.10
TE_P_SENS <- 0.05

CAST_HIGH <- 0.80

MIN_COR_N <- 8L
MIN_TE_GROUP_N <- 30L

# Restore original displayed F-Net lines
MAX_DRAW_EDGES <- 180L

# NULL = show complete backbone
# MAX_DRAW_EDGES <- NULL

COL_TN <- "#287F9C"
COL_TP <- "#CE7953"

COL_DUAL <- "#287C6C"
COL_INFO <- "#CE7953"
COL_LOAD <- "#377F9E"
COL_NEITHER <- "#AFBBC1"
COL_MISSING <- "#555B60"

FIG_WIDTH <- 19
FIG_HEIGHT <- 18


# ==============================================================
# 3. HELPERS
# ==============================================================

need <- function(path) {
  
  if (!file.exists(path)) {
    stop(
      "\nRequired input does not exist:\n",
      path
    )
  }
}

save_csv <- function(x, name) {
  
  readr::write_csv(
    x,
    file.path(OUT, name),
    na = ""
  )
}

clean_id <- function(x) {
  
  x <- trimws(as.character(x))
  
  x <- sub("\\.0+$", "", x)
  
  x[x %in% c("", "NA", "NaN")] <- NA_character_
  
  x
}

bool <- function(x) {
  
  as.character(x) %in%
    c("TRUE", "True", "true", "T", "1")
}

optional_col <- function(d, name, default = NA) {
  
  if (name %in% names(d)) {
    d[[name]]
  } else {
    rep(default, nrow(d))
  }
}

pick_col <- function(d, choices, what) {
  
  found <- choices[choices %in% names(d)]
  
  if (!length(found)) {
    
    stop(
      "\nMissing ", what, "\nColumns:\n",
      paste(names(d), collapse = ", ")
    )
  }
  
  found[1]
}

rho_safe <- function(x, y, min_n = MIN_COR_N) {
  
  good <- is.finite(x) & is.finite(y)
  
  x <- x[good]
  y <- y[good]
  
  if (
    length(x) < min_n ||
    length(unique(x)) < 2 ||
    length(unique(y)) < 2
  ) {
    return(NA_real_)
  }
  
  suppressWarnings(
    cor(x, y, method = "spearman")
  )
}

top_fraction <- function(x, frac = HOT_FRAC) {
  
  ans <- rep(NA, length(x))
  
  valid <- which(is.finite(x))
  
  if (!length(valid)) return(ans)
  
  ans[valid] <- FALSE
  
  k <- ceiling(length(valid) * frac)
  
  ranked <- valid[
    order(-x[valid], valid)
  ]
  
  ans[head(ranked, k)] <- TRUE
  
  ans
}

theme_fig <- function(size = 9) {
  
  theme_bw(base_size = size) +
    
    theme(
      panel.grid.minor = element_blank(),
      
      strip.text = element_text(
        face = "bold",
        size = 9
      ),
      
      plot.title = element_text(
        face = "bold",
        size = 11
      ),
      
      plot.subtitle = element_text(
        size = 8,
        color = "grey40"
      ),
      
      legend.position = "bottom",
      
      legend.text = element_text(size = 8),
      
      legend.title = element_text(size = 8)
    )
}

read_sf_safe <- function(path) {
  
  need(path)
  
  x <- sf::st_read(path, quiet = TRUE)
  
  if (is.na(sf::st_crs(x))) {
    
    bb <- sf::st_bbox(x)
    
    plausible_lonlat <- (
      all(is.finite(bb)) &&
        bb["xmin"] >= -180 &&
        bb["xmax"] <= 180 &&
        bb["ymin"] >= -90 &&
        bb["ymax"] <= 90
    )
    
    if (!plausible_lonlat) {
      stop(
        "Missing CRS cannot be inferred:\n",
        path
      )
    }
    
    warning(
      "EPSG:4326 assumed for ",
      basename(path),
      ". Verify original CRS."
    )
    
    sf::st_crs(x) <- 4326
  }
  
  sf::st_transform(x, MAP_CRS)
}

wilson <- function(k, n, level = 0.95) {
  
  if (n <= 0 || !is.finite(n)) {
    return(c(NA_real_, NA_real_))
  }
  
  z <- qnorm(1 - (1 - level) / 2)
  
  p <- k / n
  
  den <- 1 + z^2 / n
  
  center <- (
    p + z^2 / (2 * n)
  ) / den
  
  rad <- z / den * sqrt(
    p * (1 - p) / n +
      z^2 / (4 * n^2)
  )
  
  c(
    max(0, center - rad),
    min(1, center + rad)
  )
}


# ==============================================================
# 4. CHECK CORE SOURCE FILES
# ==============================================================

CORE_FILES <- c(
  PCS_FILE,
  NTN_FILE,
  CASTDM_FILE,
  STATION_FILE,
  BASIN_FILE,
  RIVER_FILE,
  CAST_FILE,
  CAST_SHP,
  TE_FILES,
  file.path(
    FUSION_DIR, "TN",
    "Fusion_DualLag_backbone.csv"
  ),
  file.path(
    FUSION_DIR, "TP",
    "Fusion_DualLag_backbone.csv"
  ),
  file.path(
    FUSION_DIR, "TN",
    "Fusion_all_edges.csv"
  ),
  file.path(
    FUSION_DIR, "TP",
    "Fusion_all_edges.csv"
  )
)

missing_core <- CORE_FILES[
  !file.exists(CORE_FILES)
]

if (length(missing_core)) {
  
  stop(
    "\nMissing required files:\n",
    paste(missing_core, collapse = "\n")
  )
}

cat("\nAll core source files found.\n")


# ==============================================================
# 5. READ PCS
# ==============================================================

pcs_raw <- readr::read_csv(
  PCS_FILE,
  
  col_types = readr::cols(
    STAID = readr::col_character(),
    .default = readr::col_guess()
  ),
  
  show_col_types = FALSE
)

pcs <- pcs_raw %>%
  
  transmute(
    STAID = clean_id(STAID),
    
    nutrient = toupper(
      as.character(nutrient)
    ),
    
    PCS = as.numeric(
      PCS_recomputed
    ),
    
    out_strength = as.numeric(out_strength),
    in_strength = as.numeric(in_strength),
    
    out_degree = as.numeric(out_degree),
    in_degree = as.numeric(in_degree)
  ) %>%
  
  filter(
    nutrient %in% c("TN", "TP"),
    is.finite(PCS),
    PCS >= 0
  ) %>%
  
  group_by(nutrient) %>%
  
  mutate(
    PCS_top20_full = top_fraction(PCS),
    
    PCS_pct_full = percent_rank(PCS)
  ) %>%
  
  ungroup()

if (
  anyDuplicated(
    pcs[c("STAID", "nutrient")]
  )
) {
  stop("Duplicate PCS station/nutrient.")
}

save_csv(
  pcs,
  "01_PCS_nodes.csv"
)

cat("\nPCS nodes:\n")

print(
  pcs %>% count(nutrient)
)


# ==============================================================
# 6. READ ORIGINAL F-NET EDGES
# ==============================================================

read_fusion_edges <- function(
    nu,
    filename
) {
  
  path <- file.path(
    FUSION_DIR, nu, filename
  )
  
  x <- readr::read_csv(
    path,
    
    col_types = readr::cols(
      from = readr::col_character(),
      to = readr::col_character(),
      .default = readr::col_guess()
    ),
    
    show_col_types = FALSE
  )
  
  tibble(
    nutrient = nu,
    
    from = clean_id(x$from),
    to = clean_id(x$to),
    
    W_fusion = as.numeric(x$W_fusion),
    
    R = suppressWarnings(
      as.numeric(optional_col(x, "R"))
    ),
    
    physical_status = as.character(
      optional_col(x, "physical_status")
    ),
    
    lag_assigned = bool(
      optional_col(
        x, "lag_assigned", FALSE
      )
    )
  ) %>%
    
    filter(
      is.finite(W_fusion),
      W_fusion > 0,
      !is.na(from),
      !is.na(to),
      from != to
    )
}

backbone <- bind_rows(
  read_fusion_edges(
    "TN", "Fusion_DualLag_backbone.csv"
  ),
  
  read_fusion_edges(
    "TP", "Fusion_DualLag_backbone.csv"
  )
)

full_edges <- bind_rows(
  read_fusion_edges(
    "TN", "Fusion_all_edges.csv"
  ),
  
  read_fusion_edges(
    "TP", "Fusion_all_edges.csv"
  )
)

if (
  anyDuplicated(
    full_edges[c("nutrient", "from", "to")]
  )
) {
  stop("Duplicate full F-Net edges.")
}

save_csv(
  backbone,
  "02_FNet_backbone.csv"
)

save_csv(
  full_edges,
  "02_FNet_full_edges.csv"
)


# ==============================================================
# 6B. FULL-UPSTREAM CATCHMENT AREAS FOR NTN AND CAST-DM YIELD
# ==============================================================
# IMPORTANT: An ntnwatersheds polygon is only suitable when it represents
# the whole upstream catchment of the station, NOT a local incremental polygon.
# For publication, prefer NTN_watershed_areas.csv with STAID,area_km2,
# derived from VERIFIED full upstream delineations.

if (file.exists(AREA_FILE)) {
  area_raw <- readr::read_csv(AREA_FILE,
                              col_types=readr::cols(STAID=readr::col_character(), .default=readr::col_guess()),
                              show_col_types=FALSE)
  if (!all(c("STAID","area_km2") %in% names(area_raw))) {
    stop("AREA_FILE must contain STAID and area_km2")
  }
  areas <- area_raw %>%
    transmute(STAID=clean_id(STAID), area_km2=as.numeric(area_km2),
              area_source="verified_upstream_csv")
} else {
  area_poly <- read_sf_safe(BASIN_FILE)
  area_id <- intersect(c("SITE_NO","STAID","site_no","SiteNo","station_id"),
                       names(area_poly))
  if (!length(area_id)) {
    stop("No site identifier in ntnwatersheds.shp. Provide NTN_watershed_areas.csv")
  }
  area_poly <- area_poly %>%
    mutate(STAID=clean_id(.data[[area_id[1]]])) %>%
    filter(!is.na(STAID)) %>%
    sf::st_transform(5070) %>%
    group_by(STAID) %>%
    summarise(.groups="drop")
  areas <- area_poly %>% sf::st_drop_geometry() %>%
    mutate(area_km2=as.numeric(sf::st_area(area_poly))/1e6,
           area_source="shapefile_geometry_UNVERIFIED_upstream")
  warning("Full upstream watershed area NOT verified. Inspect the area QA; ",
          "a partial drainage polygon is NOT suitable for yield.")
}
if (anyDuplicated(areas$STAID)) stop("Duplicate catchment area station IDs")
if (any(!is.finite(areas$area_km2) | areas$area_km2 <= 0)) {
  stop("Invalid catchment areas; inspect NTN_watershed_areas.csv")
}
save_csv(areas, "01_station_upstream_areas_QA.csv")

# ==============================================================
# 7. READ NTN MONTHLY LOADS
# ==============================================================

ntn_raw <- readr::read_csv(
  NTN_FILE,
  
  col_types = readr::cols(
    STAID = readr::col_character(),
    PCODE = readr::col_character(),
    .default = readr::col_guess()
  ),
  
  show_col_types = FALSE
)

if (!NTN_LOAD_FIELD %in% names(ntn_raw)) {
  
  stop(
    "NTN load field missing: ",
    NTN_LOAD_FIELD,
    "\nAvailable columns: ",
    paste(names(ntn_raw), collapse = ", ")
  )
}

ntn <- ntn_raw %>%
  
  transmute(
    STAID = clean_id(STAID),
    
    PCODE = trimws(as.character(PCODE)),
    
    year = as.integer(Year),
    month = as.integer(Month),
    
    monthly_load = as.numeric(
      .data[[NTN_LOAD_FIELD]]
    )
  ) %>%
  
  mutate(
    nutrient = case_when(
      PCODE == "P00600" ~ "TN",
      PCODE == "P00665" ~ "TP",
      TRUE ~ NA_character_
    )
  ) %>%
  
  filter(
    nutrient %in% c("TN", "TP"),
    year >= NTN_START_YEAR,
    year <= NTN_END_YEAR,
    month %in% 1:12
  )

if (
  anyDuplicated(
    ntn[c(
      "STAID", "nutrient",
      "year", "month"
    )]
  )
) {
  stop("Duplicate NTN site-month nutrient loads.")
}

ntn_annual <- ntn %>%
  group_by(STAID, nutrient, year) %>%
  summarise(
    valid_months = sum(is.finite(monthly_load) & monthly_load >= 0),
    annual_load = if (sum(is.finite(monthly_load) & monthly_load >= 0) == MIN_MONTHS) {
      sum(monthly_load, na.rm = TRUE) * NTN_KG_FACTOR
    } else NA_real_,
    .groups="drop"
  ) %>%
  left_join(areas, by="STAID") %>%
  mutate(yield_kg_km2_yr=annual_load/area_km2)

ntn_site <- ntn_annual %>%
  group_by(STAID, nutrient) %>%
  summarise(
    valid_years=sum(is.finite(yield_kg_km2_yr)),
    # Keep field used throughout original Fig. 4, but now contents are NTN yield.
    mean_annual_load=if (sum(is.finite(yield_kg_km2_yr)) >= MIN_YEARS) {
      mean(yield_kg_km2_yr,na.rm=TRUE)
    } else NA_real_,
    .groups="drop"
  )

save_csv(ntn_annual,"03_NTN_annual_yield.csv")
save_csv(ntn_site,"03_NTN_mean_annual_yield.csv")
save_csv(tibble::tibble(
  field="NTN KLoad", assumed_unit=NTN_KLOAD_UNIT,
  multiplier_to_kg=NTN_KG_FACTOR,
  check="01486000 2012 Q*C vs KLoad; ratio mean=1.0038",
  metadata_confirmed=FALSE
),"03_NTN_unit_assumption_QA.csv")

# The official METRIC CAST-DM station-year load is stored in pounds/year.
# Convert to kg/km2/year with the SAME upstream area denominator as NTN.
dm_raw <- readr::read_csv(CASTDM_FILE,
                          col_types=readr::cols(STAID=readr::col_character(),.default=readr::col_guess()),
                          show_col_types=FALSE)
if (!all(c("STAID","nutrient","year","annual_load_lb") %in% names(dm_raw))) {
  stop("CAST-DM table must have STAID,nutrient,year,annual_load_lb")
}
dm_annual <- dm_raw %>%
  transmute(STAID=clean_id(STAID),nutrient=toupper(trimws(nutrient)),
            year=as.integer(year),annual_load_lb=as.numeric(annual_load_lb)) %>%
  filter(nutrient %in% c("TN","TP"),
         year >= NTN_START_YEAR, year <= NTN_END_YEAR,
         is.finite(annual_load_lb),annual_load_lb >= 0) %>%
  left_join(areas, by="STAID") %>%
  mutate(CASTDM_yield=annual_load_lb*0.45359237/area_km2)
if (anyDuplicated(dm_annual[c("STAID","nutrient","year")])) {
  stop("Duplicate CAST-DM site/nutrient/year")
}
dm_site <- dm_annual %>%
  group_by(STAID,nutrient) %>%
  summarise(CASTDM_valid_years=sum(is.finite(CASTDM_yield)),
            CASTDM_yield=if (sum(is.finite(CASTDM_yield)) >= MIN_YEARS) {
              mean(CASTDM_yield,na.rm=TRUE)
            } else NA_real_,
            .groups="drop") %>%
  group_by(nutrient) %>%
  mutate(CASTDM_pct=percent_rank(CASTDM_yield)) %>%
  ungroup()
save_csv(dm_annual,"03_CASTDM_annual_yield.csv")
save_csv(dm_site,"03_CASTDM_mean_yield.csv")

# Cross-check data availability and the common site-year coverage.
annual_compare <- inner_join(
  ntn_annual %>% select(STAID,nutrient,year,NTN_yield=yield_kg_km2_yr),
  dm_annual %>% select(STAID,nutrient,year,CASTDM_yield),
  by=c("STAID","nutrient","year"))
save_csv(annual_compare,"03_NTN_CASTDM_common_year_yield_QA.csv")


# ==============================================================
# 8. READ CAST — 92 REGIONS
# ==============================================================

cast_raw <- readr::read_csv(
  CAST_FILE,
  show_col_types = FALSE
)

cast <- cast_raw %>%
  
  mutate(
    Geography = trimws(
      as.character(Geography)
    ),
    
    # Example:
    # ANATF_DC-ANACOSTIA... -> ANATF_DC
    CBSEG_92 = trimws(
      sub("-.*$", "", Geography)
    ),
    
    nutrient = toupper(
      as.character(nutrient)
    )
  ) %>%
  
  filter(
    scenario_year == CAST_YEAR,
    scenario == CAST_SCENARIO,
    nutrient %in% c("TN", "TP")
  )

if (
  anyDuplicated(
    cast[c("CBSEG_92", "nutrient")]
  )
) {
  stop("Duplicate CAST region/nutrient rows.")
}

cast_shapes <- read_sf_safe(
  CAST_SHP
) %>%
  
  mutate(
    CBSEG_92 = trimws(
      as.character(CBSEG_92)
    )
  ) %>%
  
  filter(
    !is.na(CBSEG_92),
    CBSEG_92 != ""
  )

# Dissolve land-river segments to 92 CAST regions

cast_regions <- cast_shapes %>%
  
  group_by(CBSEG_92) %>%
  
  summarise(
    .groups = "drop"
  )

unmatched_cast <- anti_join(
  cast %>% distinct(CBSEG_92),
  
  cast_regions %>%
    st_drop_geometry() %>%
    distinct(CBSEG_92),
  
  by = "CBSEG_92"
)

save_csv(
  unmatched_cast,
  "QA_CAST_unmatched_regions.csv"
)

if (nrow(unmatched_cast)) {
  
  stop(
    "CAST Geography cannot be matched ",
    "to CBSEG_92 shapefile. ",
    "Inspect QA_CAST_unmatched_regions.csv."
  )
}

cast_map <- cast_regions %>%
  
  left_join(
    cast %>%
      select(
        CBSEG_92,
        nutrient,
        EOS_lb_yr,
        EOT_lb_yr,
        EOS_lb_acre_yr,
        area_acres
      ),
    
    by = "CBSEG_92"
  ) %>%
  
  filter(
    nutrient %in% c("TN", "TP")
  ) %>%
  
  group_by(nutrient) %>%
  
  mutate(
    CAST_pct = percent_rank(
      EOS_lb_acre_yr
    )
  ) %>%
  
  ungroup()

if (
  any(!is.finite(cast_map$EOS_lb_acre_yr))
) {
  stop(
    "CAST EOS loading rate contains missing ",
    "or non-finite values."
  )
}

save_csv(
  st_drop_geometry(cast_map),
  "04_CAST_regions.csv"
)

st_write(
  cast_map,
  
  file.path(
    OUT, "04_CAST_2019_regions.gpkg"
  ),
  
  delete_dsn = TRUE,
  quiet = TRUE
)


# ==============================================================
# 9. STATION COORDINATES
# ==============================================================

stations_raw <- readr::read_csv(
  STATION_FILE,
  
  col_types = readr::cols(
    SITE_NO = readr::col_character(),
    .default = readr::col_guess()
  ),
  
  show_col_types = FALSE
)

stations <- stations_raw %>%
  
  transmute(
    STAID = clean_id(SITE_NO),
    
    lon = as.numeric(lon),
    lat = as.numeric(lat)
  ) %>%
  
  filter(
    !is.na(STAID),
    is.finite(lon),
    is.finite(lat),
    
    lon >= -180,
    lon <= 180,
    
    lat >= -90,
    lat <= 90
  ) %>%
  
  distinct(
    STAID,
    .keep_all = TRUE
  )

# Only pcs has no coordinates.
# This avoids the old lon.x / lon.y issue.

pcs_geo <- pcs %>%
  
  left_join(
    stations,
    by = "STAID"
  )

save_csv(
  pcs_geo %>%
    
    filter(
      !is.finite(lon) |
        !is.finite(lat)
    ),
  
  "QA_PCS_missing_coordinates.csv"
)

pcs_geo <- pcs_geo %>%
  
  filter(
    is.finite(lon),
    is.finite(lat)
  )

pcs_sf <- sf::st_as_sf(
  pcs_geo,
  
  coords = c("lon", "lat"),
  
  crs = 4326,
  remove = FALSE
) %>%
  
  sf::st_transform(MAP_CRS)

station_sf <- sf::st_as_sf(
  stations,
  
  coords = c("lon", "lat"),
  
  crs = 4326,
  remove = FALSE
) %>%
  
  sf::st_transform(MAP_CRS)

basin <- read_sf_safe(BASIN_FILE)

rivers <- read_sf_safe(RIVER_FILE)


# ==============================================================
# 10. MATCH PCS STATIONS TO CAST REGIONS
# ==============================================================

# Only one unique CAST polygon assignment
# is accepted for station-level context.

hits <- sf::st_intersects(
  pcs_sf,
  cast_regions
)

pcs_cast <- pcs_sf %>%
  
  mutate(
    n_region_matches = lengths(hits),
    
    CBSEG_92 = vapply(
      hits,
      
      function(idx) {
        
        if (length(idx) == 1L) {
          
          as.character(
            cast_regions$CBSEG_92[idx]
          )
          
        } else {
          
          NA_character_
        }
      },
      
      character(1)
    )
  ) %>%
  
  st_drop_geometry() %>%
  
  left_join(
    st_drop_geometry(cast_map) %>%
      
      select(
        nutrient,
        CBSEG_92,
        EOS_lb_acre_yr,
        CAST_pct
      ),
    
    by = c("nutrient", "CBSEG_92")
  )

save_csv(
  pcs_cast,
  "05_PCS_CAST_station_context.csv"
)

save_csv(
  pcs_cast %>%
    
    count(
      nutrient,
      n_region_matches
    ),
  
  "05_CAST_station_matching_QA.csv"
)


# ==============================================================
# 11. UNIFIED NODE-LEVEL TABLE
# ==============================================================

station_data <- pcs_cast %>%
  left_join(ntn_site,by=c("STAID","nutrient")) %>%
  left_join(dm_site,by=c("STAID","nutrient"))

if (
  anyDuplicated(
    station_data[c("STAID", "nutrient")]
  )
) {
  stop("Duplicate rows in unified station table.")
}


# ==============================================================
# 12. COMMON-SAMPLE PRIORITIES
# ==============================================================

station_data <- station_data %>%
  
  group_by(nutrient) %>%
  
  mutate(
    common_valid = (
      is.finite(PCS) &
        is.finite(mean_annual_load) &
        mean_annual_load >= 0
    ),
    
    PCS_common = if_else(
      common_valid, PCS, NA_real_
    ),
    
    load_common = if_else(
      common_valid,
      mean_annual_load,
      NA_real_
    ),
    
    high_PCS = top_fraction(PCS_common),
    
    high_load = top_fraction(load_common),
    
    PCS_pct_common = percent_rank(
      PCS_common
    ),
    
    NTN_pct_common = percent_rank(
      load_common
    ),
    
    priority_class = case_when(
      !common_valid ~ "Not evaluated",
      
      high_PCS & high_load ~
        "Dual priority",
      
      high_PCS & !high_load ~
        "Information priority",
      
      !high_PCS & high_load ~
        "Loading priority",
      
      TRUE ~ "Neither"
    )
  ) %>%
  
  ungroup()

priority_levels <- c(
  "Dual priority",
  "Information priority",
  "Loading priority",
  "Neither",
  "Not evaluated"
)

station_data <- station_data %>%
  
  mutate(
    priority_class = factor(
      priority_class,
      levels = priority_levels
    )
  )

save_csv(
  station_data,
  "06_unified_station_metrics.csv"
)

priority_counts <- station_data %>%
  
  count(
    nutrient,
    priority_class,
    name = "n_stations",
    .drop = FALSE
  )

save_csv(
  priority_counts,
  "06_priority_counts.csv"
)

priority_overlap <- station_data %>%
  
  filter(common_valid) %>%
  
  group_by(nutrient) %>%
  
  summarise(
    n = n(),
    
    n_high_PCS = sum(high_PCS),
    
    n_high_load = sum(high_load),
    
    n_both = sum(high_PCS & high_load),
    
    overlap_fraction =
      n_both / n_high_PCS,
    
    expected_fraction =
      n_high_load / n,
    
    enrichment =
      overlap_fraction / expected_fraction,
    
    .groups = "drop"
  )

save_csv(
  priority_overlap,
  "06_priority_overlap.csv"
)


# ==============================================================
# 13. CAST REGION-LEVEL STATISTICS
# ==============================================================

# One observation per unique CAST region.
# Stations inside the same region do not
# become independent CAST model outputs.

cast_region_nodes <- station_data %>%
  
  filter(
    n_region_matches == 1,
    !is.na(CBSEG_92),
    is.finite(EOS_lb_acre_yr)
  ) %>%
  
  group_by(
    nutrient,
    CBSEG_92
  ) %>%
  
  summarise(
    n_stations = n(),
    
    mean_PCS = mean(PCS),
    
    mean_out_strength = mean(
      out_strength
    ),
    
    mean_in_strength = mean(
      in_strength
    ),
    
    CAST_rate = first(
      EOS_lb_acre_yr
    ),
    
    CAST_pct = first(CAST_pct),
    
    .groups = "drop"
  )

cast_region_cor <- cast_region_nodes %>%
  
  group_by(nutrient) %>%
  
  summarise(
    n_regions = n(),
    
    rho_PCS = rho_safe(
      mean_PCS,
      CAST_rate
    ),
    
    rho_out = rho_safe(
      mean_out_strength,
      CAST_rate
    ),
    
    rho_in = rho_safe(
      mean_in_strength,
      CAST_rate
    ),
    
    .groups = "drop"
  )

save_csv(
  cast_region_nodes,
  "07_CAST_unique_region_nodes.csv"
)

save_csv(
  cast_region_cor,
  "07_CAST_unique_region_correlations.csv"
)


# ==============================================================
# 14. PCS VS STRENGTH
# ==============================================================

pcs_strength_QA <- pcs %>%
  
  group_by(nutrient) %>%
  
  summarise(
    n = n(),
    
    rho_PCS_out = rho_safe(
      PCS,
      out_strength
    ),
    
    rho_PCS_in = rho_safe(
      PCS,
      in_strength
    ),
    
    .groups = "drop"
  )

save_csv(
  pcs_strength_QA,
  "07_PCS_strength_QA.csv"
)


# ==============================================================
# 15. MAP EXTENT & ORIGINAL NETWORK LINES
# ==============================================================

bb <- sf::st_bbox(basin)

dx <- as.numeric(
  bb["xmax"] - bb["xmin"]
)

dy <- as.numeric(
  bb["ymax"] - bb["ymin"]
)

XLIM <- c(
  as.numeric(bb["xmin"]) - 0.035 * dx,
  as.numeric(bb["xmax"]) + 0.035 * dx
)

YLIM <- c(
  as.numeric(bb["ymin"]) - 0.035 * dy,
  as.numeric(bb["ymax"]) + 0.035 * dy
)

xy <- sf::st_coordinates(station_sf)

station_xy <- station_sf %>%
  
  st_drop_geometry() %>%
  
  mutate(
    x = xy[, 1],
    y = xy[, 2]
  ) %>%
  
  select(STAID, x, y)

edge_xy <- backbone %>%
  
  left_join(
    station_xy %>%
      rename(
        from = STAID,
        x1 = x,
        y1 = y
      ),
    
    by = "from"
  ) %>%
  
  left_join(
    station_xy %>%
      rename(
        to = STAID,
        x2 = x,
        y2 = y
      ),
    
    by = "to"
  ) %>%
  
  filter(
    is.finite(x1),
    is.finite(y1),
    is.finite(x2),
    is.finite(y2)
  )

if (is.null(MAX_DRAW_EDGES)) {
  
  draw_edges <- edge_xy
  
} else {
  
  draw_edges <- edge_xy %>%
    
    group_by(nutrient) %>%
    
    arrange(
      desc(W_fusion),
      .by_group = TRUE
    ) %>%
    
    slice_head(
      n = MAX_DRAW_EDGES
    ) %>%
    
    ungroup()
}

draw_edges <- draw_edges %>%
  
  mutate(
    lag_class = if_else(
      lag_assigned,
      "Lag assigned",
      "Lag unassigned"
    )
  )

save_csv(
  draw_edges,
  "08_displayed_FNet_backbone.csv"
)


# ============================================================
# 16. PANELS A-B — BASIN + RIVER NETWORK ONLY
# ============================================================


# ============================================================
# FIX FIG4(a-b) COLUMN NAMES
# ============================================================

# 1. PCS节点字段统一
pcs_sf <- pcs_sf %>%
  mutate(
    PCS_top20 = PCS_top20_full
  )

# 2. 网络边字段统一
# 兼容独立版中的 draw_edges / lag_class
edge_draw <- draw_edges %>%
  mutate(
    lag_status = if_else(
      lag_assigned,
      "Lag assigned",
      "Lag unassigned"
    )
  )

# 3. 检查字段
cat("\nPCS columns:\n")
print(names(pcs_sf))

cat("\nEdge columns:\n")
print(names(edge_draw))

cat("\nPCS top-20% counts:\n")
print(
  pcs_sf %>%
    sf::st_drop_geometry() %>%
    count(nutrient, PCS_top20)
)

cat("\nDisplayed backbone edges:\n")
print(
  edge_draw %>%
    count(nutrient)
)


make_network_map <- function(nu, letter) {
  
  nodes <- pcs_sf %>%
    filter(nutrient == nu)
  
  links <- edge_draw %>%
    filter(nutrient == nu)
  
  ggplot() +
    
    # Basin background
    geom_sf(
      data = basin,
      fill = "#F5F5F5",
      color = "#CFCFCF",
      linewidth = 0.25
    ) +
    
    # River network
    geom_sf(
      data = rivers,
      color = "#9FC7D9",
      linewidth = 0.18,
      alpha = 0.75
    ) +
    
    # Original F-Net backbone lines
    geom_segment(
      data = links,
      aes(
        x = x1,
        y = y1,
        xend = x2,
        yend = y2,
        linetype = lag_status
      ),
      color = "#407B8E",
      linewidth = 0.30,
      alpha = 0.42
    ) +
    
    # PCS nodes
    geom_sf(
      data = nodes,
      aes(
        size = PCS,
        color = PCS_top20
      ),
      shape = 21,
      fill = "white",
      stroke = 0.55
    ) +
    
    scale_color_manual(
      values = c(
        "FALSE" = "#73848B",
        "TRUE" = "#C46D48"
      ),
      labels = c(
        "Other nodes",
        "Top-20% PCS"
      ),
      name = "Information hub"
    ) +
    
    scale_size_continuous(
      range = c(1, 5),
      name = "PCS"
    ) +
    
    scale_linetype_manual(
      values = c(
        "Lag assigned" = "solid",
        "Lag unassigned" = "dashed"
      ),
      name = "Dual-lag"
    ) +
    
    coord_sf(
      xlim = XLIM,
      ylim = YLIM,
      expand = FALSE,
      datum = NA
    ) +
    
    labs(
      title = paste0(
        "(", letter, ") ",
        nu,
        " F-Net structure on basin river network"
      ),
      x = NULL,
      y = NULL
    ) +
    
    theme_fig(9) +
    
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
}

p_a <- make_network_map("TN", "a")
p_b <- make_network_map("TP", "b")
# 先单独测试
print(p_a)
print(p_b)
# ==============================================================
# 17. PANEL (c): NETWORK-LOADING CORRELATIONS
# ==============================================================

c_data <- station_data %>%
  
  select(
    STAID,
    nutrient,
    PCS,
    out_strength,
    in_strength,
    mean_annual_load,
    CASTDM_yield
  ) %>%
  
  pivot_longer(
    cols = c(
      PCS,
      out_strength,
      in_strength
    ),
    
    names_to = "network_metric",
    values_to = "network_value"
  ) %>%
  
  pivot_longer(
    cols = c(
      mean_annual_load,
      CASTDM_yield
    ),
    
    names_to = "load_metric",
    values_to = "load_value"
  ) %>%
  
  mutate(
    load_metric = recode(
      load_metric,
      
      mean_annual_load = "NTN yield",
      
      CASTDM_yield = "CAST-DM yield"
    )
  )

c_stats <- c_data %>%
  
  group_by(
    nutrient,
    network_metric,
    load_metric
  ) %>%
  
  summarise(
    n = sum(
      is.finite(network_value) &
        is.finite(load_value)
    ),
    
    rho = rho_safe(
      network_value,
      load_value
    ),
    
    .groups = "drop"
  )

save_csv(
  c_stats,
  "09_node_loading_correlations.csv"
)

c_plot <- c_stats %>%
  
  mutate(
    row_label = paste(
      nutrient,
      network_metric
    ),
    
    row_label = factor(
      row_label,
      
      levels = rev(c(
        "TN PCS",
        "TN out_strength",
        "TN in_strength",
        "TP PCS",
        "TP out_strength",
        "TP in_strength"
      ))
    ),
    
    label = if_else(
      is.finite(rho),
      
      sprintf(
        "%.2f\nn=%d",
        rho, n
      ),
      
      paste0(
        "NA\nn=", n
      )
    )
  )

p_c <- ggplot(
  c_plot,
  
  aes(
    x = load_metric,
    y = row_label,
    fill = rho
  )
) +
  
  geom_tile(
    color = "white",
    linewidth = 0.75,
    width = 0.96,
    height = 0.94
  ) +
  
  geom_text(
    aes(label = label),
    size = 3.1
  ) +
  
  scale_fill_gradient2(
    low = "#3C82A3",
    mid = "#F8F8F8",
    high = "#CC7953",
    midpoint = 0,
    limits = c(-1, 1),
    na.value = "#D8D8D8",
    name = "Spearman rho"
  ) +
  
  labs(
    title = "(c) Network importance and loading context",
    
    subtitle = paste0(
      "NTN and CAST-DM station yields"
    ),
    
    x = NULL,
    y = NULL
  ) +
  
  theme_fig() +
  
  theme(
    panel.grid = element_blank(),
    
    axis.text.x = element_text(
      angle = 20,
      hjust = 1
    )
  )


# ==============================================================
# 18. PANEL (d): PRIORITY DISCORDANCE
# ==============================================================

d_data <- station_data %>%
  
  filter(common_valid) %>%
  
  mutate(
    log_load = log10(
      mean_annual_load + 1
    )
  )

d_cuts <- d_data %>%
  
  group_by(nutrient) %>%
  
  summarise(
    xcut = min(
      log_load[high_load]
    ),
    
    ycut = min(
      PCS[high_PCS]
    ),
    
    .groups = "drop"
  )

d_labels <- d_data %>%
  
  filter(
    priority_class %in% c(
      "Information priority",
      "Loading priority"
    )
  ) %>%
  
  mutate(
    label_score = if_else(
      priority_class == "Information priority",
      
      PCS_pct_common,
      
      NTN_pct_common
    )
  ) %>%
  
  group_by(
    nutrient,
    priority_class
  ) %>%
  
  arrange(
    desc(label_score),
    .by_group = TRUE
  ) %>%
  
  slice_head(n = 3) %>%
  
  ungroup()

p_d <- ggplot(
  d_data,
  
  aes(
    x = log_load,
    y = PCS
  )
) +
  
  geom_vline(
    data = d_cuts,
    
    aes(xintercept = xcut),
    
    inherit.aes = FALSE,
    
    linetype = "dashed",
    color = "grey60"
  ) +
  
  geom_hline(
    data = d_cuts,
    
    aes(yintercept = ycut),
    
    inherit.aes = FALSE,
    
    linetype = "dashed",
    color = "grey60"
  ) +
  
  geom_point(
    aes(fill = priority_class),
    
    shape = 21,
    color = "grey35",
    stroke = 0.30,
    
    size = 2.2,
    alpha = 0.92
  ) +
  
  geom_text(
    data = d_labels,
    
    aes(label = STAID),
    
    size = 2.2,
    vjust = -0.65,
    check_overlap = TRUE
  ) +
  
  facet_wrap(
    ~nutrient,
    
    ncol = 2,
    scales = "free"
  ) +
  
  scale_fill_manual(
    values = c(
      "Dual priority" = COL_DUAL,
      "Information priority" = COL_INFO,
      "Loading priority" = COL_LOAD,
      "Neither" = COL_NEITHER
    ),
    
    drop = TRUE,
    name = "Candidate priority"
  ) +
  
  scale_x_continuous(
    breaks = scales::breaks_pretty(n = 4)
  ) +
  
  labs(
    title = "(d) Information–loading priority discordance",
    
    subtitle = paste0(
      "Common matched sample; ",
      "top 20% within TN and TP"
    ),
    
    x = "log10(mean annual NTN yield + 1)",
    y = "Structural PCS"
  ) +
  
  theme_fig()


# ==============================================================
# 19. READ TE RESULTS
# ==============================================================

read_te <- function(file, nu) {
  
  need(file)
  
  x <- readr::read_csv(
    file,
    show_col_types = FALSE
  )
  
  from_name <- pick_col(
    x,
    c("from", "source", "Source"),
    "TE source"
  )
  
  to_name <- pick_col(
    x,
    c("to", "target", "Target"),
    "TE target"
  )
  
  p_name <- pick_col(
    x,
    c("p_value", "p", "p_val"),
    "TE p-value"
  )
  
  te_name <- pick_col(
    x,
    c(
      "TE_bias_corrected",
      "TE_corrected",
      "TE",
      "te"
    ),
    "TE value"
  )
  
  lag_name <- pick_col(
    x,
    c(
      "lag_months",
      "lag",
      "best_lag"
    ),
    "TE lag"
  )
  
  x %>%
    
    transmute(
      nutrient = nu,
      
      from = clean_id(
        .data[[from_name]]
      ),
      
      to = clean_id(
        .data[[to_name]]
      ),
      
      p_value = as.numeric(
        .data[[p_name]]
      ),
      
      TE_value = as.numeric(
        .data[[te_name]]
      ),
      
      TE_lag = as.numeric(
        .data[[lag_name]]
      )
    ) %>%
    
    filter(
      !is.na(from),
      !is.na(to),
      from != to
    )
}

te_all <- bind_rows(
  read_te(TE_FILES[["TN"]], "TN"),
  read_te(TE_FILES[["TP"]], "TP")
)

if (
  anyDuplicated(
    te_all[c("nutrient", "from", "to")]
  )
) {
  stop(
    "Duplicate TE directed pairs. ",
    "Check whether all lags were stored separately."
  )
}


# ==============================================================
# 20. TE CANDIDATE COVERAGE QA
# ==============================================================

te_coverage <- te_all %>%
  
  group_by(nutrient) %>%
  
  summarise(
    n_rows = n(),
    
    n_nodes = n_distinct(
      c(from, to)
    ),
    
    n_possible_directed =
      n_nodes * (n_nodes - 1),
    
    n_evaluated = sum(
      is.finite(p_value) &
        p_value >= 0 &
        p_value <= 1 &
        is.finite(TE_value)
    ),
    
    coverage = n_rows /
      n_possible_directed,
    
    .groups = "drop"
  )

save_csv(
  te_coverage,
  "10_TE_candidate_coverage.csv"
)

if (
  any(te_coverage$coverage < 0.90)
) {
  
  warning(
    "TE edge coverage below 90%. ",
    "Confirm that the source contains all ",
    "evaluated candidate directed pairs."
  )
}


# ==============================================================
# 21. MATCH TE EDGES TO CAST BACKGROUND
# ==============================================================

contexts <- station_data %>%
  
  select(
    nutrient,
    STAID,
    CASTDM_pct
  ) %>%
  
  distinct()

from_context <- contexts %>%
  
  rename(
    from = STAID,
    
    from_pct = CASTDM_pct
  )

to_context <- contexts %>%
  
  rename(
    to = STAID,
    
    to_pct = CASTDM_pct
  )

te_context <- te_all %>%
  
  left_join(
    from_context,
    
    by = c("nutrient", "from")
  ) %>%
  
  left_join(
    to_context,
    
    by = c("nutrient", "to")
  ) %>%
  
  mutate(
    loading_context = case_when(
      is.finite(from_pct) &
        is.finite(to_pct) &
        from_pct >= CAST_HIGH &
        to_pct >= CAST_HIGH ~
        "Both high",
      
      is.finite(from_pct) &
        is.finite(to_pct) &
        from_pct < CAST_HIGH &
        to_pct < CAST_HIGH ~
        "Both lower",
      
      is.finite(from_pct) &
        is.finite(to_pct) ~
        "Mixed",
      
      TRUE ~ NA_character_
    ),
    
    evaluated = (
      is.finite(p_value) &
        p_value >= 0 &
        p_value <= 1 &
        is.finite(TE_value)
    ),
    
    sig_010 = (
      evaluated &
        p_value < TE_P &
        TE_value > 0
    ),
    
    sig_005 = (
      evaluated &
        p_value < TE_P_SENS &
        TE_value > 0
    )
  )

save_csv(
  te_context,
  "10_TE_all_pairs_CASTDM_yield_context.csv"
)


# ==============================================================
# 22. TE SIGNIFICANCE RATES
# ==============================================================

te_rates <- te_context %>%
  
  filter(
    !is.na(loading_context),
    evaluated
  ) %>%
  
  group_by(
    nutrient,
    loading_context
  ) %>%
  
  summarise(
    n_evaluated = n(),
    
    n_sig_010 = sum(sig_010),
    
    n_sig_005 = sum(sig_005),
    
    n_sources = n_distinct(from),
    
    n_targets = n_distinct(to),
    
    rate_010 =
      n_sig_010 / n_evaluated,
    
    rate_005 =
      n_sig_005 / n_evaluated,
    
    .groups = "drop"
  )

ci <- lapply(
  seq_len(nrow(te_rates)),
  
  function(i) {
    
    wilson(
      te_rates$n_sig_010[i],
      te_rates$n_evaluated[i]
    )
  }
)

te_rates$ci_low <- vapply(
  ci, function(x) x[1], numeric(1)
)

te_rates$ci_high <- vapply(
  ci, function(x) x[2], numeric(1)
)

te_rates <- te_rates %>%
  
  mutate(
    loading_context = factor(
      loading_context,
      
      levels = c(
        "Both lower",
        "Mixed",
        "Both high"
      )
    ),
    
    small_sample = (
      n_evaluated < MIN_TE_GROUP_N
    ),
    
    label = sprintf(
      "%d/%d",
      n_sig_010,
      n_evaluated
    ),
    
    label = if_else(
      small_sample,
      paste0(label, " *"),
      label
    )
  )

save_csv(
  te_rates,
  "11_TE_significance_rates_CI.csv"
)


# ==============================================================
# 23. PANEL (e): TE FRACTION + DESCRIPTIVE INTERVALS
# ==============================================================

# Wilson intervals below assume independent
# Bernoulli trials mathematically.
#
# Here directed pairs share monitoring stations.
# Therefore, intervals are descriptive, NOT
# valid independent-edge inferential confidence
# intervals for differences between groups.

e_top <- max(
  te_rates$ci_high,
  na.rm = TRUE
)

e_top <- min(
  1,
  max(0.30, e_top + 0.12)
)

p_e <- ggplot(
  te_rates,
  
  aes(
    x = loading_context,
    y = rate_010,
    fill = nutrient
  )
) +
  
  geom_col(
    width = 0.65,
    alpha = 0.92
  ) +
  
  geom_errorbar(
    aes(
      ymin = ci_low,
      ymax = ci_high
    ),
    
    width = 0.18,
    linewidth = 0.6,
    color = "#30383E"
  ) +
  
  geom_text(
    aes(
      y = pmin(
        1,
        ci_high + 0.025
      ),
      
      label = label
    ),
    
    size = 2.9
  ) +
  
  facet_wrap(
    ~nutrient,
    ncol = 2
  ) +
  
  scale_fill_manual(
    values = c(
      TN = COL_TN,
      TP = COL_TP
    ),
    
    guide = "none"
  ) +
  
  scale_y_continuous(
    labels = scales::label_percent(),
    
    breaks = scales::breaks_pretty(n = 5)
  ) +
  
  coord_cartesian(
    ylim = c(0, e_top)
  ) +
  
  labs(
    title = "(e) Directional TE across loading contexts",
    
    subtitle = paste0(
      "Nominal p < 0.10; ",
      "descriptive Wilson intervals; ",
      "* fewer than 30 evaluated pairs"
    ),
    
    x = "CAST-DM yield context",
    
    y = "Significant TE fraction"
  ) +
  
  theme_fig() +
  
  theme(
    panel.grid.major.x = element_blank()
  )


# ==============================================================
# 24. TE LAG AND THRESHOLD SENSITIVITY
# ==============================================================

te_lag_summary <- te_context %>%
  
  filter(
    !is.na(loading_context),
    sig_010,
    is.finite(TE_lag)
  ) %>%
  
  group_by(
    nutrient,
    loading_context
  ) %>%
  
  summarise(
    n_significant = n(),
    
    median_lag = median(TE_lag),
    
    mean_lag = mean(TE_lag),
    
    median_TE = median(TE_value),
    
    .groups = "drop"
  )

save_csv(
  te_lag_summary,
  "11_TE_lag_summary.csv"
)

save_csv(
  te_rates %>%
    
    select(
      nutrient,
      loading_context,
      n_evaluated,
      n_sig_010,
      rate_010,
      n_sig_005,
      rate_005
    ),
  
  "11_TE_p_threshold_sensitivity.csv"
)


# ==============================================================
# 25. PANEL (f): CANDIDATE MONITORING PRIORITIES
# ==============================================================

# The same matched station set is used for
# PCS and NTN top-20% rankings.
#
# CAST remains background context.

priority_summary <- priority_counts %>%
  
  filter(
    priority_class %in% c(
      "Dual priority",
      "Information priority",
      "Loading priority"
    )
  ) %>%
  
  select(
    nutrient,
    priority_class,
    n_stations
  ) %>%
  
  pivot_wider(
    names_from = priority_class,
    values_from = n_stations,
    values_fill = 0
  ) %>%
  
  mutate(
    count_label = paste0(
      "Dual: ", `Dual priority`,
      " | Info: ", `Information priority`,
      " | Load: ", `Loading priority`
    )
  ) %>%
  
  select(
    nutrient,
    count_label
  )


f_map <- pcs_sf %>%
  
  left_join(
    station_data %>%
      select(
        STAID,
        nutrient,
        priority_class
      ),
    
    by = c("STAID", "nutrient")
  )

# Since pcs_sf has no priority_class already,
# this join is safe and retains every node.

f_map <- f_map %>%
  
  mutate(
    priority_class = factor(
      priority_class,
      levels = priority_levels
    )
  )

p_f <- ggplot() +
  
  geom_sf(
    data = cast_map,
    aes(fill = CAST_pct),
    color = NA
  ) +
  
  geom_sf(
    data = rivers,
    color = "#B8D5DF",
    linewidth = 0.12,
    alpha = 0.75
  ) +
  
  # white halo
  geom_sf(
    data = f_map,
    aes(shape = priority_class),
    color = "white",
    fill = NA,
    size = 3.1,
    stroke = 1.3,
    show.legend = FALSE
  ) +
  
  # 第二层：在相同位置叠加较小的彩色形状
  geom_sf(
    data = f_map,
    aes(
      shape = priority_class,
      color = priority_class
    ),
    size = 2.2,
    stroke = 0.65,
    alpha = 0.98
  ) +
  
  
  facet_wrap(
    ~nutrient,
    ncol = 2
  ) +
  
  scale_fill_viridis_c(
    option = "viridis",
    limits = c(0, 1),
    labels = scales::label_percent(),
    name = "CAST EOS rate\npercentile"
  ) +
  
  scale_shape_manual(
    values = c(
      "Dual priority" = 17,
      "Information priority" = 16,
      "Loading priority" = 15,
      "Neither" = 1,
      "Not evaluated" = 4
    ),
    drop = TRUE,
    name = "Monitoring class"
  ) +
  
  scale_color_manual(
    values = c(
      "Dual priority" = COL_DUAL,
      "Information priority" = COL_INFO,
      "Loading priority" = COL_LOAD,
      "Neither" = COL_NEITHER,
      "Not evaluated" = COL_MISSING
    ),
    drop = TRUE,
    name = "Monitoring class"
  ) +
  
  coord_sf(
    xlim = XLIM,
    ylim = YLIM,
    datum = NA,
    expand = FALSE
  ) +
  
  labs(
    title = "(f) Candidate complementary monitoring priorities",
    subtitle = paste0(
      "Common-sample PCS and NTN top-20%; ",
      "counts reported by nutrient"
    ),
    x = NULL,
    y = NULL
  ) +
  
  theme_fig() +
  
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

# ==============================================================
# 26. SPARROW SUPPLEMENT AND CROSSWALK QA
# ==============================================================

# SPARROW is retained as a supplemental model
# reference until a verified station–reach
# crosswalk is available.
#
# Do not treat ntn_site_stack as a confirmed
# STAID mapping without topology QA.

if (
  all(file.exists(c(SP_INPUT, SP_FILES)))
) {
  
  sp_input <- readr::read_csv(
    SP_INPUT,
    show_col_types = FALSE
  ) %>%
    
    filter(
      modelspace == SP_MODELSPACE
    )
  
  if (
    "ntn_site_stack" %in% names(sp_input)
  ) {
    
    sp_candidate <- sp_input %>%
      
      transmute(
        comid_stack = clean_id(
          comid_stack
        ),
        
        ntn_site_stack = clean_id(
          ntn_site_stack
        )
      ) %>%
      
      filter(
        ntn_site_stack %in% pcs$STAID
      )
    
    save_csv(
      sp_candidate,
      "QA_SPARROW_exact_ID_candidates.csv"
    )
  }
  
  read_sp <- function(nu) {
    
    x <- readr::read_csv(
      SP_FILES[[nu]],
      show_col_types = FALSE
    )
    
    if (anyDuplicated(x$comid_stack)) {
      stop("Duplicate SPARROW prediction key.")
    }
    
    input_small <- sp_input %>%
      select(
        comid_stack,
        comid
      )
    
    y <- x %>%
      
      select(
        comid_stack,
        TOTAL_YIELD,
        INC_TOTAL_YIELD,
        PLOAD_TOTAL,
        PLOAD_INC_TOTAL
      ) %>%
      
      inner_join(
        input_small,
        by = "comid_stack"
      ) %>%
      
      mutate(nutrient = nu)
    
    if (nrow(y) != nrow(sp_input)) {
      stop("Incomplete SPARROW join: ", nu)
    }
    
    y
  }
  
  sparrow <- bind_rows(
    read_sp("TN"),
    read_sp("TP")
  )
  
  save_csv(
    sparrow,
    "12_SPARROW_2012_yields.csv"
  )
  
  p_sp <- sparrow %>%
    
    filter(
      is.finite(TOTAL_YIELD),
      is.finite(INC_TOTAL_YIELD),
      TOTAL_YIELD >= 0,
      INC_TOTAL_YIELD >= 0
    ) %>%
    
    ggplot(
      aes(
        x = log10(INC_TOTAL_YIELD + 1),
        y = log10(TOTAL_YIELD + 1)
      )
    ) +
    
    geom_bin2d(bins = 40) +
    
    facet_wrap(
      ~nutrient,
      ncol = 2
    ) +
    
    scale_fill_viridis_c(
      trans = "log1p",
      name = "Reaches"
    ) +
    
    labs(
      title = "SPARROW 2012 nutrient yield structure",
      
      subtitle = paste0(
        "Model-reach distributions; ",
        "not yet verified against NTN stations"
      ),
      
      x = "log10(incremental yield + 1)",
      y = "log10(total yield + 1)"
    ) +
    
    theme_fig()
  
  ggsave(
    file.path(
      OUT, "FigS4_SPARROW.pdf"
    ),
    
    p_sp,
    width = 10,
    height = 6,
    bg = "white"
  )
  
} else {
  
  warning(
    "SPARROW files incomplete. ",
    "Skipping SPARROW supplemental plot."
  )
}


# ==============================================================
# 27. CAST UNIQUE-REGION SCATTER SUPPLEMENT
# ==============================================================

cast_cor_labels <- cast_region_cor %>%
  
  mutate(
    label = paste0(
      "rho = ",
      
      if_else(
        is.finite(rho_PCS),
        
        sprintf("%.2f", rho_PCS),
        "NA"
      ),
      
      "; regions = ", n_regions
    )
  )

p_cast_region <- ggplot(
  cast_region_nodes,
  
  aes(
    x = CAST_pct,
    y = mean_PCS
  )
) +
  
  geom_point(
    aes(
      size = n_stations,
      color = nutrient
    ),
    
    alpha = 0.75
  ) +
  
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    
    se = FALSE,
    linewidth = 0.8,
    
    color = "#34434A"
  ) +
  
  geom_label(
    data = cast_cor_labels,
    
    aes(
      x = -Inf,
      y = Inf,
      label = label
    ),
    
    inherit.aes = FALSE,
    
    hjust = -0.05,
    vjust = 1.15,
    
    size = 3
  ) +
  
  facet_wrap(
    ~nutrient,
    ncol = 2,
    scales = "free_y"
  ) +
  
  scale_x_continuous(
    limits = c(0, 1),
    
    labels = scales::label_percent()
  ) +
  
  scale_color_manual(
    values = c(
      TN = COL_TN,
      TP = COL_TP
    ),
    
    guide = "none"
  ) +
  
  labs(
    title = "CAST regional loading versus F-Net PCS",
    
    subtitle = paste0(
      "One observation per unique CAST ",
      "region with monitoring stations"
    ),
    
    x = "CAST EOS loading-rate percentile",
    
    y = "Mean station PCS per CAST region",
    
    size = "Stations"
  ) +
  
  theme_fig()


# ==============================================================
# 27B. NEW PANEL (g) — SPARROW vs F-Net PCS
# ==============================================================
# SPARROW 2012 static reach-level predicted total load.
# The TN and TP model output BOTH use PLOAD_TOTAL;
# nutrient identity comes from the corresponding model file.
# The selected 2012 NTN station–reach correspondence is
# provided by the previously produced matching result.
# This panel does not modify the original panels (a–f).

SP_MATCH_FILE <- file.path(
  ROOT, "Fig5_FNet_SPARROW_verified",
  "02_SPARROW_NTN_station_loads.csv"
)

if (!file.exists(SP_MATCH_FILE)) {
  stop(
    "Panel (g) requires the previously validated NTN–SPARROW ",
    "crosswalk result:\n", SP_MATCH_FILE,
    "\nRun the previously completed SPARROW matching step first."
  )
}

sp_matched_raw <- readr::read_csv(
  SP_MATCH_FILE,
  col_types = readr::cols(
    STAID = readr::col_character(),
    comid_stack = readr::col_character(),
    .default = readr::col_guess()
  ),
  show_col_types = FALSE
)

sp_require <- c("STAID", "nutrient", "comid_stack", "sparrow_load")
sp_missing <- setdiff(sp_require, names(sp_matched_raw))
if (length(sp_missing)) {
  stop("SPARROW match table missing: ", paste(sp_missing, collapse=", "))
}

sp_matched <- sp_matched_raw %>%
  transmute(
    STAID = clean_id(STAID),
    nutrient = toupper(as.character(nutrient)),
    comid_stack = clean_id(comid_stack),
    model_year = if ("model_year" %in% names(sp_matched_raw))
      as.integer(model_year) else SP_MODELSPACE,
    sparrow_load = as.numeric(sparrow_load)
  ) %>%
  filter(
    nutrient %in% c("TN", "TP"),
    model_year == SP_MODELSPACE,
    !is.na(STAID),
    is.finite(sparrow_load),
    sparrow_load >= 0
  ) %>%
  distinct()

if (anyDuplicated(sp_matched[c("nutrient", "STAID")])) {
  stop("Panel (g): duplicated NTN station / nutrient SPARROW matches")
}

# Use original 2011–2023 F-Net PCS in THIS fig4.R,
# rather than importing the 2011–2020 Fig6 PCS.
sp_pcs <- pcs %>%
  select(STAID, nutrient, PCS) %>%
  filter(is.finite(PCS), PCS >= 0)

# Preserve original station-reach crosswalk but use published SPARROW
# TOTAL_YIELD rather than the unnormalized sparrow_load.
sp_yield_all <- bind_rows(lapply(c("TN","TP"), function(nu) {
  predfile <- SP_FILES[[nu]]
  if (!file.exists(predfile)) stop("Missing SPARROW prediction: ",predfile)
  z <- readr::read_csv(predfile,
                       col_types=readr::cols(comid_stack=readr::col_character(), .default=readr::col_guess()),
                       show_col_types=FALSE)
  if (!"TOTAL_YIELD" %in% names(z)) stop("SPARROW lacks TOTAL_YIELD: ",predfile)
  z %>% transmute(comid_stack=clean_id(comid_stack),nutrient=nu,
                  SPARROW_yield=as.numeric(TOTAL_YIELD))
}))
if (anyDuplicated(sp_yield_all[c("nutrient","comid_stack")])) {
  stop("Repeated SPARROW nutrient/COMID yield keys")
}
sp_g <- sp_matched %>%
  inner_join(sp_pcs, by = c("STAID", "nutrient")) %>%
  inner_join(sp_yield_all,by=c("nutrient","comid_stack")) %>%
  mutate(
    log_sparrow = log10(SPARROW_yield + 1),
    log_PCS = log10(PCS + 1)
  )

sp_all_qc <- sp_matched %>%
  count(nutrient, name = "n_sparrow_matched")

sp_join_qc <- sp_g %>%
  group_by(nutrient) %>%
  summarise(
    n_both = n(),
    n_distinct_comids = n_distinct(comid_stack),
    rho = rho_safe(PCS, SPARROW_yield),
    .groups = "drop"
  )

sp_qa <- sp_all_qc %>%
  full_join(sp_join_qc, by = "nutrient") %>%
  arrange(nutrient)

save_csv(sp_g, "13_SPARROW_PCS_matched_station_panel_g.csv")
save_csv(sp_qa, "13_SPARROW_PCS_matching_QA.csv")

if (!all(c("TN", "TP") %in% sp_g$nutrient)) {
  stop("Panel (g) has no matched observations for TN or TP")
}

sp_lab <- sp_qa %>%
  mutate(
    label = ifelse(
      is.finite(rho),
      sprintf("Spearman rho = %.2f; n = %d", rho, n_both),
      sprintf("Spearman rho = NA; n = %d", n_both)
    ),
    x = -Inf,
    y = Inf
  )

p_g <- ggplot(sp_g, aes(x = log_sparrow, y = log_PCS)) +
  geom_point(aes(color = nutrient), size = 2.0, alpha = 0.78) +
  geom_smooth(
    aes(color = nutrient),
    method = "lm", formula = y ~ x,
    se = TRUE, linewidth = 0.70, alpha = 0.12
  ) +
  geom_text(
    data = sp_lab,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = -0.05, vjust = 1.4,
    size = 3.4, color = "#364650"
  ) +
  facet_wrap(~nutrient, nrow = 1, scales = "free") +
  scale_color_manual(values = c(TN = COL_TN, TP = COL_TP), guide = "none") +
  labs(
    title = "(g) SPARROW nutrient loading versus F-Net structural importance",
    subtitle = paste0(
      "SPARROW 2012 long-term model reference; ",
      "station–reach matched subset; Spearman rho on original values"
    ),
    x = "log10(SPARROW predicted total yield + 1)",
    y = "log10(F-Net PCS + 1)"
  ) +
  theme_fig()

# Export the added panel by itself.
ggsave(
  file.path(OUT, "Fig4_panel_g_SPARROW.pdf"),
  plot = p_g, width = 15, height = 5.2, bg = "white"
)
ggsave(
  file.path(OUT, "Fig4_panel_g_SPARROW.png"),
  plot = p_g, width = 15, height = 5.2, dpi = 400, bg = "white"
)

# ==============================================================
# 28. COMBINE FINAL FIGURE
# ==============================================================

fig4 <- (
  (p_a | p_b) /
    (p_c | p_d) /
    (p_e | p_f) /
    p_g
) +
  patchwork::plot_layout(
    heights = c(1.15, 0.98, 1.00, 0.83)
  ) +
  patchwork::plot_annotation(
    title = paste0(
      "Complementarity between Nutrient Loading ",
      "Patterns and Information-Network Structure"
    ),
    subtitle = paste0(
      "F-Net (2011–2023), ",
      "NTN/CAST-DM yield (2012–2022), ",
      "CAST 2019 background; SPARROW 2012 reference"
    ),
    caption = paste0(
      "PCS is a weighted multi-step structural score ",
      "(H=5, alpha=0.85), not nutrient mass. ",
      "Static CAST remains the original regional map background; ",
      "Panel (g) adds a station-matched SPARROW static loading reference. ",
      "Directional TE significance is nominal ",
      "(p<0.10); Wilson intervals are descriptive."
    )
  )

# 先测试整图能否渲染
print(fig4)

# 导出
ggsave(
  file.path(OUT, "Fig4_final.pdf"),
  plot = fig4,
  width = FIG_WIDTH,
  height = FIG_HEIGHT + 4,
  bg = "white"
)

ggsave(
  file.path(OUT, "Fig4_final.png"),
  plot = fig4,
  width = FIG_WIDTH,
  height = FIG_HEIGHT + 4,
  dpi = 400,
  bg = "white"
)

# ==============================================================
# 30. SUMMARY AND QA
# ==============================================================

cat("\n====================================\n")
cat("FIGURE 4 FINISHED\n")
cat("====================================\n")

cat("\nPCS stations:\n")
print(pcs %>% count(nutrient))

cat("\nMatched PCS–NTN stations:\n")

print(
  station_data %>%
    
    group_by(nutrient) %>%
    
    summarise(
      n_total = n(),
      n_matched = sum(common_valid),
      .groups = "drop"
    )
)

cat("\nNode vs loading correlation:\n")
print(c_stats, n = Inf)

cat("\nCAST unique-region correlations:\n")
print(cast_region_cor, n = Inf)

cat("\nPCS vs strength:\n")
print(pcs_strength_QA, n = Inf)

cat("\nPriority overlap:\n")
print(priority_overlap, n = Inf)

cat("\nPriority counts:\n")
print(priority_counts, n = Inf)

cat("\nTE candidate coverage:\n")
print(te_coverage, n = Inf)

cat("\nTE significance rates:\n")
print(te_rates, n = Inf)

cat("\nTE lag summary:\n")
print(te_lag_summary, n = Inf)

cat("\nOUTPUT:\n")
cat(OUT, "\n")

print(fig4)

# Area provenance: confirm station polygons truly represent full upstream basin.
cat("\nCAST-DM station overlap with PCS:\n")
print(station_data %>% group_by(nutrient) %>% summarise(
  n_pcs=n(), n_castdm=sum(is.finite(CASTDM_yield)),
  n_ntn=sum(is.finite(mean_annual_load)),
  n_both=sum(is.finite(mean_annual_load)&is.finite(CASTDM_yield)), .groups="drop"))
cat("\nNOTE: Fig. 4(f) polygon fill is original static CAST EOS.\n")
