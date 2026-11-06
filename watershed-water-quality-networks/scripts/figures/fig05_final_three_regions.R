
# ================================================================
# FIGURE 5 — COMPLETE FINAL SCRIPT
#
# Anomaly-triggered complementary monitoring prioritization
# using source-conditioned F-Net
#
# Six selected cases:
#   Region 1 North:   TN + TP
#   Region 2 Central: TN + TP
#   Region 3 South:   TN + TP
#
# (a) Combined TN monitoring map
# (b) Combined TP monitoring map
# (c) Top-5 anomaly verification
# (d) All dual-lag attention windows
# (e) PCS priority vs observed outcomes
# (f) Representative observed anomaly series
#
# This script uses PRECOMPUTED results.
# No F-Net, PCS or dual-lag recalculation.
# ================================================================


# ================================================================
# 0. PACKAGES
# ================================================================

PKGS <- c(
  "sf", "readr", "dplyr", "tidyr",
  "tibble", "ggplot2", "patchwork",
  "lubridate", "scales"
)

for (p in PKGS) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

invisible(
  lapply(PKGS, library, character.only = TRUE)
)


# ================================================================
# 1. PATHS
# ================================================================

source("config.R")
ROOT <- get_data_root()

BASE <- file.path(
  ROOT, "Fig5_anomaly_trigger"
)

STEP51 <- file.path(
  BASE, "Step5_1_event_catalog"
)

REGION_DIR <- file.path(
  BASE, "Fig5_three_regions"
)

OUT <- file.path(
  BASE, "Fig5_main_final"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

# Previously completed region results

CASE_FILE <- file.path(
  REGION_DIR, "05_six_selected_events.csv"
)

TOP5_FILE <- file.path(
  REGION_DIR, "06_six_cases_Top5.csv"
)

MONTH_FILE <- file.path(
  REGION_DIR, "07_six_cases_monthly_verification.csv"
)

RESULT_FILE <- file.path(
  REGION_DIR, "08_six_cases_target_verification.csv"
)

ANOMALY_FILE <- file.path(
  STEP51, "04_all_monthly_concentration_anomalies.csv"
)

STATION_FILE <- file.path(
  ROOT, "NTN-station.txt"
)

BASIN_FILE <- file.path(
  ROOT, "ntnwatersheds", "ntnwatersheds.shp"
)

RIVER_FILE <- file.path(
  ROOT, "水系", "水系.shp"
)


# ================================================================
# 2. SETTINGS
# ================================================================

TOP_K <- 5L

MAX_LAG <- 3L

Z_THRESHOLD <- 1.5

MAP_CRS <- 5070

REGION_NAMES <- c(
  "Region 1 — North",
  "Region 2 — Central",
  "Region 3 — South"
)

REGION_SHORT <- c(
  "Region 1 — North" = "North",
  "Region 2 — Central" = "Central",
  "Region 3 — South" = "South"
)

# Only STATIONS are colored by region.
REGION_COLS <- c(
  "Region 1 — North" = "#267DA8",
  "Region 2 — Central" = "#DE8842",
  "Region 3 — South" = "#4C9B77"
)

# Verification status colors
OUTCOME_COLS <- c(
  "Hit" = "#2C9077",
  "No anomaly" = "#CE7656",
  "Not evaluable" = "#9EA8AD"
)

WINDOW_COLS <- c(
  "Attention" = "#438EAA",
  "Outside" = "#E8EDF0",
  "Unassigned" = "#BAC3C8"
)

COL_SOURCE <- "#C9554B"

FIG_WIDTH <- 17
FIG_HEIGHT <- 17


# ================================================================
# 3. HELPER FUNCTIONS
# ================================================================

clean_id <- function(x) {
  
  x <- trimws(as.character(x))
  
  x <- sub("\\.0+$", "", x)
  
  x[x %in% c("", "NA", "NaN")] <- NA_character_
  
  ii <- !is.na(x) &
    grepl("^[0-9]{1,7}$", x)
  
  x[ii] <- sprintf(
    "%08.0f",
    as.numeric(x[ii])
  )
  
  x
}


read_ids <- function(
    path,
    id_columns = character()
) {
  
  if (!file.exists(path)) {
    stop("Missing input file:\n", path)
  }
  
  spec <- readr::cols(
    .default = readr::col_guess()
  )
  
  for (nm in id_columns) {
    spec$cols[[nm]] <- readr::col_character()
  }
  
  readr::read_csv(
    path,
    col_types = spec,
    show_col_types = FALSE
  )
}


as_bool <- function(x) {
  
  if (is.logical(x)) {
    return(x)
  }
  
  y <- tolower(
    trimws(as.character(x))
  )
  
  ifelse(
    is.na(y),
    NA,
    y %in% c("true", "t", "1", "yes")
  )
}


save_csv <- function(d, name) {
  
  readr::write_csv(
    d,
    file.path(OUT, name),
    na = ""
  )
}


check_columns <- function(d, required, label) {
  
  missing <- setdiff(
    required, names(d)
  )
  
  if (length(missing)) {
    stop(
      label,
      " missing columns: ",
      paste(missing, collapse = ", ")
    )
  }
}


theme_fig <- function(base_size = 10) {
  
  theme_bw(base_size = base_size) +
    
    theme(
      panel.grid.minor = element_blank(),
      
      panel.grid.major = element_line(
        color = "#E9EDEF",
        linewidth = 0.2
      ),
      
      plot.title = element_text(
        face = "bold",
        size = base_size + 0.5,
        margin = margin(b = 5)
      ),
      
      plot.subtitle = element_text(
        color = "#647078",
        size = base_size - 1
      ),
      
      strip.background = element_rect(
        fill = "#F2F4F5",
        color = "#D9E0E3"
      ),
      
      strip.text = element_text(
        face = "bold",
        size = base_size - 1
      ),
      
      legend.title = element_text(
        face = "bold",
        size = 9
      ),
      
      legend.text = element_text(
        size = 8
      ),
      
      plot.margin = margin(7, 7, 7, 7)
    )
}


# ================================================================
# 4. CHECK ALL INPUT FILES
# ================================================================

input_files <- c(
  CASE_FILE = CASE_FILE,
  TOP5_FILE = TOP5_FILE,
  MONTH_FILE = MONTH_FILE,
  RESULT_FILE = RESULT_FILE,
  ANOMALY_FILE = ANOMALY_FILE
)

print(
  data.frame(
    file = names(input_files),
    exists = file.exists(input_files)
  )
)

if (!all(file.exists(input_files))) {
  stop("Some Fig.5 input files are missing.")
}


# ================================================================
# 5. READ THE SIX SELECTED EVENTS
# ================================================================

cases <- read_ids(
  CASE_FILE,
  "source"
) %>%
  mutate(
    event_id = as.character(event_id),
    source = clean_id(source),
    region = as.character(region),
    nutrient = toupper(
      as.character(nutrient)
    ),
    trigger_date = as.Date(trigger_date),
    source_Z = as.numeric(source_Z)
  ) %>%
  mutate(
    region_short = unname(
      REGION_SHORT[region]
    ),
    
    case_label = paste0(
      region_short, " | ", nutrient
    ),
    
    case_detail = paste0(
      region_short, " | ", nutrient,
      " | ", event_id,
      " | ", format(trigger_date, "%b %Y")
    )
  )

check_columns(
  cases,
  c(
    "event_id", "region", "nutrient",
    "source", "trigger_date", "source_Z"
  ),
  "Selected cases"
)

if (anyDuplicated(cases$event_id)) {
  stop("Duplicate event_id in CASE_FILE.")
}

if (nrow(cases) != 6) {
  stop(
    "Expected 6 selected cases; found ",
    nrow(cases)
  )
}

if (anyNA(cases$region_short)) {
  stop(
    "Unknown region names in CASE_FILE. ",
    "Check REGION_NAMES."
  )
}

if (anyDuplicated(
  cases[c("region", "nutrient")]
)) {
  stop(
    "More than one selected event per region/nutrient."
  )
}

cases <- cases %>%
  mutate(
    region = factor(
      region,
      levels = REGION_NAMES
    ),
    nutrient = factor(
      nutrient,
      levels = c("TN", "TP")
    )
  ) %>%
  arrange(region, nutrient) %>%
  mutate(
    region = as.character(region),
    nutrient = as.character(nutrient)
  )

save_csv(
  cases,
  "01_six_events_used.csv"
)


# ================================================================
# 6. READ LOCAL TOP-5
# ================================================================

top5 <- read_ids(
  TOP5_FILE,
  c("source", "target")
) %>%
  mutate(
    event_id = as.character(event_id),
    source = clean_id(source),
    target = clean_id(target),
    rank_region = as.integer(rank_region),
    source_PCS = as.numeric(source_PCS),
    region = as.character(region),
    nutrient = toupper(
      as.character(nutrient)
    ),
    window_defined = as_bool(window_defined),
    lag_min = as.integer(lag_min),
    lag_max = as.integer(lag_max)
  ) %>%
  semi_join(
    cases %>% select(event_id),
    by = "event_id"
  )

check_columns(
  top5,
  c(
    "event_id", "target", "rank_region",
    "source_PCS", "lag_min", "lag_max",
    "window_defined"
  ),
  "Top-5"
)

if (anyDuplicated(
  top5[c("event_id", "target")]
)) {
  stop("Duplicate event-target in TOP5_FILE.")
}

if (nrow(top5) != 6 * TOP_K) {
  stop(
    "Expected 30 Top-5 target records; found ",
    nrow(top5)
  )
}

if (any(
  is.na(top5$rank_region) |
  top5$rank_region < 1 |
  top5$rank_region > TOP_K
)) {
  stop("Invalid rank_region.")
}

# Add ONLY display labels to avoid .x/.y collisions.

case_labels <- cases %>%
  select(
    event_id,
    region_short,
    case_label,
    case_detail
  )

top5 <- top5 %>%
  select(
    -any_of(c(
      "region_short", "case_label", "case_detail"
    ))
  ) %>%
  left_join(
    case_labels,
    by = "event_id"
  ) %>%
  mutate(
    target_label = paste0(
      "#", rank_region, " ", target
    ),
    
    window_label = case_when(
      !window_defined ~ "Unassigned",
      is.na(lag_min) | is.na(lag_max) ~
        "Unassigned",
      lag_min > MAX_LAG ~ ">3 months",
      lag_min == lag_max ~
        paste0(lag_min, " mo"),
      TRUE ~ paste0(
        lag_min, "–", lag_max, " mo"
      )
    )
  )


# ================================================================
# 7. READ ACTUAL MONTHLY VERIFICATION
# ================================================================

months <- read_ids(
  MONTH_FILE,
  c("source", "target")
) %>%
  mutate(
    event_id = as.character(event_id),
    source = clean_id(source),
    target = clean_id(target),
    lag_months = as.integer(lag_months),
    target_observed = as_bool(target_observed),
    target_anomaly_15 = as_bool(target_anomaly_15),
    target_Z = as.numeric(target_Z),
    attention = as_bool(attention)
  ) %>%
  semi_join(
    cases %>% select(event_id),
    by = "event_id"
  )

check_columns(
  months,
  c(
    "event_id", "target",
    "lag_months", "target_observed",
    "target_anomaly_15", "target_Z",
    "attention"
  ),
  "Monthly verification"
)

if (anyDuplicated(
  months[c("event_id", "target", "lag_months")]
)) {
  stop("Duplicate event-target-month records.")
}

months <- months %>%
  select(
    -any_of(c(
      "region_short", "case_label",
      "case_detail", "target_label"
    ))
  ) %>%
  left_join(
    top5 %>%
      select(
        event_id,
        target,
        region_short,
        case_label,
        case_detail,
        target_label
      ),
    by = c("event_id", "target")
  )


# ================================================================
# 8. READ TARGET VERIFICATION
# ================================================================

results <- read_ids(
  RESULT_FILE,
  "target"
) %>%
  mutate(
    event_id = as.character(event_id),
    target = clean_id(target),
    rank_region = as.integer(rank_region),
    evaluable = as_bool(evaluable),
    outcome = as.character(outcome)
  ) %>%
  semi_join(
    cases %>% select(event_id),
    by = "event_id"
  )

check_columns(
  results,
  c(
    "event_id", "target",
    "rank_region", "evaluable", "outcome"
  ),
  "Target verification"
)

if (anyDuplicated(
  results[c("event_id", "target")]
)) {
  stop("Duplicate event-target verification.")
}

results <- results %>%
  select(
    -any_of(c(
      "source_PCS", "region_short",
      "case_label", "case_detail"
    ))
  ) %>%
  left_join(
    top5 %>%
      select(
        event_id, target, source_PCS,
        region_short, case_label, case_detail
      ),
    by = c("event_id", "target")
  )

results$outcome <- factor(
  results$outcome,
  levels = names(OUTCOME_COLS)
)

if (anyNA(results$outcome)) {
  stop(
    "Unexpected outcome label in RESULT_FILE."
  )
}


# ================================================================
# 9. READ ACTUAL ANOMALY TIME SERIES
# ================================================================

anomaly <- read_ids(
  ANOMALY_FILE,
  "STAID"
) %>%
  mutate(
    STAID = clean_id(STAID),
    nutrient = toupper(
      as.character(nutrient)
    ),
    date = as.Date(date),
    Z = as.numeric(Z)
  )

check_columns(
  anomaly,
  c("STAID", "nutrient", "date", "Z"),
  "Anomaly series"
)


# ================================================================
# 10. READ SPATIAL DATA
# ================================================================

# Prefer spatial objects already loaded in R.
# Read original files only if the object is absent.

if (!exists("basin")) {
  
  basin <- sf::st_read(
    BASIN_FILE,
    quiet = TRUE
  )
}

if (!exists("rivers")) {
  
  rivers <- sf::st_read(
    RIVER_FILE,
    quiet = TRUE
  )
}

if (!exists("stations_sf")) {
  
  stations_raw <- readr::read_csv(
    STATION_FILE,
    col_types = readr::cols(
      .default = readr::col_character()
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
      is.finite(lat)
    ) %>%
    distinct(STAID, .keep_all = TRUE)
  
  stations_sf <- sf::st_as_sf(
    stations,
    coords = c("lon", "lat"),
    crs = 4326
  )
}

# Original river CRS sometimes missing.
# Assume WGS84 ONLY if bounds are plausible.

fix_crs <- function(x, label) {
  
  if (is.na(sf::st_crs(x))) {
    
    bb <- sf::st_bbox(x)
    
    plausible <- (
      bb["xmin"] >= -180 &&
        bb["xmax"] <= 180 &&
        bb["ymin"] >= -90 &&
        bb["ymax"] <= 90
    )
    
    if (!plausible) {
      stop(
        label,
        " has missing CRS and non-geographic coordinates."
      )
    }
    
    warning(
      label,
      " CRS missing; assuming EPSG:4326. ",
      "Verify the original spatial metadata."
    )
    
    sf::st_crs(x) <- 4326
  }
  
  sf::st_transform(
    x, MAP_CRS
  )
}

basin <- fix_crs(basin, "Basin")
rivers <- fix_crs(rivers, "Rivers")
stations_sf <- fix_crs(stations_sf, "Stations")

stations_sf$STAID <- clean_id(
  stations_sf$STAID
)

if (anyDuplicated(stations_sf$STAID)) {
  stop("Duplicate STAID in stations_sf.")
}

# Spatial coordinates

station_coordinates <- sf::st_coordinates(
  stations_sf
)

station_xy <- tibble(
  STAID = stations_sf$STAID,
  x = station_coordinates[, 1],
  y = station_coordinates[, 2]
)

# Main map extent

bb <- sf::st_bbox(basin)

XLIM <- as.numeric(
  bb[c("xmin", "xmax")]
)

YLIM <- as.numeric(
  bb[c("ymin", "ymax")]
)


# ================================================================
# 11. PREPARE SOURCE AND TARGET SPATIAL OBJECTS
# ================================================================

# These are explicitly created here.
# No missing seg_df, source_pts or target_pts.

source_pts <- stations_sf %>%
  select(STAID) %>%
  inner_join(
    cases %>%
      select(
        event_id,
        source,
        nutrient,
        region,
        region_short
      ),
    by = c("STAID" = "source")
  )

target_pts <- stations_sf %>%
  select(STAID) %>%
  inner_join(
    top5 %>%
      select(
        event_id,
        target,
        nutrient,
        region,
        region_short,
        rank_region,
        source_PCS
      ),
    by = c("STAID" = "target")
  )

# Schematic straight connections:
# These are NOT actual river transport paths.

seg_df <- top5 %>%
  select(
    event_id,
    nutrient,
    region,
    source,
    target,
    rank_region
  ) %>%
  left_join(
    station_xy %>%
      transmute(
        source = STAID,
        sx = x,
        sy = y
      ),
    by = "source"
  ) %>%
  left_join(
    station_xy %>%
      transmute(
        target = STAID,
        tx = x,
        ty = y
      ),
    by = "target"
  )

if (
  nrow(source_pts) != 6 ||
  nrow(target_pts) != 30 ||
  nrow(seg_df) != 30
) {
  
  stop(
    "Spatial match incomplete. Expected ",
    "6 sources, 30 targets and 30 segments."
  )
}

if (any(
  !is.finite(seg_df$sx) |
  !is.finite(seg_df$sy) |
  !is.finite(seg_df$tx) |
  !is.finite(seg_df$ty)
)) {
  stop("Missing station coordinates in seg_df.")
}

cat("\nSpatial QA:\n")

print(
  seg_df %>%
    count(nutrient, region),
  n = Inf
)


# ================================================================
# 12. SHARED MAP BACKGROUND
# ================================================================

map_base <- function() {
  
  list(
    geom_sf(
      data = basin,
      fill = "#F8FAFB",
      color = "#CBD4D8",
      linewidth = 0.16
    ),
    
    geom_sf(
      data = rivers,
      color = "#ABCEDA",
      linewidth = 0.17,
      alpha = 0.8
    ),
    
    geom_sf(
      data = stations_sf,
      color = "#A7B1B6",
      size = 0.9,
      alpha = 0.75
    )
  )
}


# ================================================================
# 13. PANEL A-B — TN AND TP COMBINED MAPS
# ================================================================

make_map_panel <- function(nu, letter) {
  
  src <- source_pts %>%
    filter(nutrient == nu)
  
  tar <- target_pts %>%
    filter(nutrient == nu)
  
  links <- seg_df %>%
    filter(nutrient == nu)
  
  ggplot() +
    
    map_base() +
    
    # Uncolored schematic links
    geom_segment(
      data = links,
      aes(
        x = sx,
        y = sy,
        xend = tx,
        yend = ty
      ),
      color = "#92A5AE",
      linewidth = 0.32,
      alpha = 0.43,
      linetype = "dashed"
    ) +
    
    # Region-colored target stations
    geom_sf(
      data = tar,
      aes(fill = region),
      shape = 21,
      color = "white",
      stroke = 0.9,
      size = 3.8
    ) +
    
    # Rank numbers
    geom_sf_text(
      data = tar,
      aes(label = rank_region),
      color = "white",
      fontface = "bold",
      size = 2.3
    ) +
    
    # Region-colored source stations
    geom_sf(
      data = src,
      aes(color = region),
      shape = 8,
      size = 4.4,
      stroke = 0.9
    ) +
    
    scale_fill_manual(
      values = REGION_COLS,
      name = "Spatial region",
      drop = FALSE
    ) +
    
    scale_color_manual(
      values = REGION_COLS,
      name = "Spatial region",
      drop = FALSE
    ) +
    
    coord_sf(
      xlim = XLIM,
      ylim = YLIM,
      datum = NA,
      expand = FALSE
    ) +
    
    labs(
      title = paste0(
        "(", letter, ") ", nu,
        ": anomaly triggers and Top-5 targets"
      ),
      
      subtitle = paste0(
        "Three source events; numbered targets ",
        "are ranked within each spatial region"
      ),
      
      x = NULL,
      y = NULL
    ) +
    
    theme_fig(9) +
    
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      
      legend.position = "bottom"
    )
}

p_a <- make_map_panel("TN", "a")
p_b <- make_map_panel("TP", "b")


# ================================================================
# 14. PANEL C — TOP-5 VERIFICATION
# ================================================================

performance <- results %>%
  group_by(
    event_id
  ) %>%
  summarise(
    n_recommended = n(),
    
    n_evaluable = sum(
      evaluable,
      na.rm = TRUE
    ),
    
    n_hits = sum(
      outcome == "Hit",
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  left_join(
    cases %>%
      select(
        event_id,
        region,
        region_short,
        nutrient
      ),
    by = "event_id"
  ) %>%
  mutate(
    hit_rate = if_else(
      n_evaluable > 0,
      n_hits / n_evaluable,
      NA_real_
    ),
    
    coverage = n_evaluable / TOP_K,
    
    label = paste0(
      n_hits, "/", n_evaluable
    ),
    
    region_short = factor(
      region_short,
      levels = c("North", "Central", "South")
    )
  )

save_csv(
  performance,
  "02_six_case_performance.csv"
)

# Display hit fraction and evaluability together.
# No evaluation: show NA rather than 0.

p_c <- ggplot(
  performance,
  aes(
    x = region_short,
    y = hit_rate,
    fill = nutrient
  )
) +
  
  geom_col(
    position = position_dodge(
      width = 0.75
    ),
    width = 0.67,
    color = "white",
    na.rm = TRUE
  ) +
  
  geom_text(
    aes(label = label),
    position = position_dodge(
      width = 0.75
    ),
    vjust = -0.4,
    size = 3,
    na.rm = TRUE
  ) +
  
  scale_fill_manual(
    values = c(
      TN = "#4387AA",
      TP = "#D88863"
    ),
    name = "Nutrient"
  ) +
  
  scale_y_continuous(
    limits = c(0, 1.12),
    breaks = seq(0, 1, 0.25),
    labels = scales::label_percent(),
    expand = expansion(mult = c(0, 0))
  ) +
  
  labs(
    title = "(c) Observed anomaly verification",
    subtitle = "Hits / evaluable Top-5 targets; coverage reported separately",
    x = NULL,
    y = "Within-window anomaly fraction"
  ) +
  
  theme_fig(10) +
  
  theme(
    legend.position = "bottom"
  )


# ================================================================
# 15. PANEL D — ALL ATTENTION WINDOWS
# ================================================================

# Build a complete 0–3-month grid for
# all 30 selected targets.

month_info <- months %>%
  select(
    event_id,
    target,
    lag_months,
    target_observed,
    target_anomaly_15,
    target_Z
  )

attention_grid <- tidyr::crossing(
  event_id = cases$event_id,
  lag_months = 0:MAX_LAG
) %>%
  inner_join(
    top5 %>%
      select(
        event_id,
        target,
        rank_region,
        window_defined,
        lag_min,
        lag_max,
        region_short,
        region,
        case_label,
        nutrient
      ),
    by = "event_id",
    relationship = "many-to-many"
  ) %>%
  left_join(
    month_info,
    by = c(
      "event_id", "target", "lag_months"
    )
  ) %>%
  mutate(
    attention = (
      window_defined %in% TRUE &
        !is.na(lag_min) &
        !is.na(lag_max) &
        lag_months >= lag_min &
        lag_months <= lag_max
    ),
    
    tile_status = case_when(
      !window_defined | is.na(window_defined) ~
        "Unassigned",
      
      attention ~ "Attention",
      
      TRUE ~ "Outside"
    ),
    
    real_anomaly = (
      target_observed %in% TRUE &
        target_anomaly_15 %in% TRUE
    ),
    
    missing_observation = (
      is.na(target_observed) |
        target_observed == FALSE
    ),
    
    target_label = paste0(
      "#", rank_region, " ", target
    )
  )

# Each facet has its own five station labels.
attention_grid <- attention_grid %>%
  mutate(
    region_short = factor(
      region_short,
      levels = c("North", "Central", "South")
    ),
    
    nutrient = factor(
      nutrient,
      levels = c("TN", "TP")
    )
  )

p_d <- ggplot(
  attention_grid,
  aes(
    x = lag_months,
    y = target_label
  )
) +
  
  geom_tile(
    aes(fill = tile_status),
    width = 0.96,
    height = 0.92,
    color = "white",
    linewidth = 0.6
  ) +
  
  geom_point(
    data = attention_grid %>%
      filter(real_anomaly),
    
    shape = 21,
    fill = COL_SOURCE,
    color = "white",
    stroke = 0.65,
    size = 2.8
  ) +
  
  geom_point(
    data = attention_grid %>%
      filter(missing_observation),
    
    shape = 4,
    color = "#58656A",
    size = 1.8,
    stroke = 0.6
  ) +
  
  facet_grid(
    nutrient ~ region_short,
    scales = "free_y",
    space = "free_y"
  ) +
  
  scale_fill_manual(
    values = WINDOW_COLS,
    name = "Window status"
  ) +
  
  scale_x_continuous(
    breaks = 0:MAX_LAG,
    labels = paste0(0:MAX_LAG, " mo")
  ) +
  
  labs(
    title = "(d) Dual-lag attention windows",
    subtitle = "Red points: observed Z > 1.5; crosses: missing observations",
    x = "Months after source anomaly",
    y = NULL
  ) +
  
  theme_fig(9) +
  
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 7.4),
    strip.text = element_text(size = 8.5)
  )


# ================================================================
# 16. PANEL E — PCS RANKING VS OBSERVED OUTCOME
# ================================================================

pcs_plot <- results %>%
  mutate(
    nutrient = factor(
      nutrient,
      levels = c("TN", "TP")
    ),
    
    region_short = factor(
      region_short,
      levels = c("North", "Central", "South")
    ),
    
    rank_region = as.integer(
      rank_region
    )
  ) %>%
  arrange(
    nutrient,
    region_short,
    rank_region
  )

p_e <- ggplot(
  pcs_plot,
  aes(
    x = rank_region,
    y = source_PCS
  )
) +
  
  geom_line(
    aes(group = event_id),
    linewidth = 0.6,
    color = "#B7C3C8"
  ) +
  
  geom_point(
    aes(
      color = outcome
    ),
    size = 2.8
  ) +
  
  scale_color_manual(
    values = OUTCOME_COLS,
    name = "Observed outcome",
    drop = FALSE
  ) +
  
  scale_x_continuous(
    breaks = 1:TOP_K
  ) +
  
  facet_grid(
    nutrient ~ region_short,
    scales = "free_y"
  ) +
  
  labs(
    title = "(e) Source-conditioned ranking and verification",
    subtitle = "Priority ranks and observed within-window outcomes",
    x = "Within-event monitoring rank",
    y = "Source-conditioned PCS"
  ) +
  
  theme_fig(9) +
  
  theme(
    strip.text = element_text(size = 8.5)
  )


# ================================================================
# 17. PANEL F — REPRESENTATIVE OBSERVED TIME SERIES
# ================================================================

# Select one TN and one TP event based on:
# 1. number of evaluable targets
# 2. source anomaly magnitude
#
# NOT based on observed target hits.

representative_cases <- performance %>%
  left_join(
    cases %>%
      select(
        event_id,
        source,
        trigger_date,
        source_Z
      ),
    by = "event_id"
  ) %>%
  group_by(nutrient) %>%
  arrange(
    desc(n_evaluable),
    desc(source_Z),
    event_id,
    .by_group = TRUE
  ) %>%
  slice_head(n = 1) %>%
  ungroup()

save_csv(
  representative_cases,
  "03_representative_cases.csv"
)

# Use source + two highest-ranked targets
# per representative event.
#
# This keeps panel (f) compact and avoids
# selection based on actual anomaly hits.

rep_targets <- top5 %>%
  semi_join(
    representative_cases %>%
      select(event_id),
    by = "event_id"
  ) %>%
  filter(
    rank_region <= 2
  )

rep_station_roles <- bind_rows(
  
  representative_cases %>%
    transmute(
      event_id,
      nutrient = as.character(nutrient),
      STAID = source,
      station_role = "Source",
      station_order = 0L
    ),
  
  rep_targets %>%
    transmute(
      event_id,
      nutrient = as.character(nutrient),
      STAID = target,
      station_role = paste0(
        "Target #", rank_region
      ),
      station_order = rank_region
    )
  
)

# Time series from one month before trigger
# to three months after trigger.

rep_series <- rep_station_roles %>%
  left_join(
    representative_cases %>%
      select(
        event_id,
        trigger_date,
        region_short
      ),
    by = "event_id"
  ) %>%
  inner_join(
    anomaly %>%
      select(
        STAID,
        nutrient,
        date,
        Z
      ),
    by = c("STAID", "nutrient"),
    relationship = "many-to-many"
  ) %>%
  filter(
    date >= trigger_date %m-% months(1),
    date <= trigger_date %m+% months(MAX_LAG)
  ) %>%
  mutate(
    relative_month =
      (lubridate::year(date) -
         lubridate::year(trigger_date)) * 12 +
      lubridate::month(date) -
      lubridate::month(trigger_date),
    
    panel_label = paste0(
      station_role,
      " | ", STAID
    )
  )

# Add original assigned windows to target curves.
rep_windows <- rep_targets %>%
  select(
    event_id,
    target,
    rank_region,
    window_defined,
    lag_min,
    lag_max
  ) %>%
  mutate(
    panel_label = paste0(
      "Target #", rank_region,
      " | ", target
    )
  ) %>%
  filter(
    window_defined %in% TRUE,
    !is.na(lag_min),
    !is.na(lag_max),
    lag_min <= MAX_LAG
  ) %>%
  transmute(
    event_id,
    panel_label,
    xmin = lag_min - 0.45,
    xmax = lag_max + 0.45
  )

# Keep each case separate.
# Each contains Source + Target #1 + Target #2.

rep_series <- rep_series %>%
  mutate(
    nutrient = factor(
      nutrient,
      levels = c("TN", "TP")
    ),
    
    station_role = factor(
      station_role,
      levels = c(
        "Source", "Target #1", "Target #2"
      )
    )
  )

p_f <- ggplot() +
  
  geom_rect(
    data = rep_windows,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = "#D9EAF0",
    alpha = 0.8,
    inherit.aes = FALSE
  ) +
  
  geom_hline(
    yintercept = Z_THRESHOLD,
    linetype = "dashed",
    color = COL_SOURCE,
    linewidth = 0.55
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    color = "#7A878E",
    linewidth = 0.5
  ) +
  
  geom_line(
    data = rep_series,
    aes(
      x = relative_month,
      y = Z,
      group = interaction(
        event_id, STAID
      ),
      color = station_role
    ),
    linewidth = 0.75,
    na.rm = TRUE
  ) +
  
  geom_point(
    data = rep_series,
    aes(
      x = relative_month,
      y = Z,
      color = station_role
    ),
    size = 1.7,
    na.rm = TRUE
  ) +
  
  geom_point(
    data = rep_series %>%
      filter(
        is.finite(Z),
        Z > Z_THRESHOLD
      ),
    aes(
      x = relative_month,
      y = Z
    ),
    shape = 21,
    fill = COL_SOURCE,
    color = "white",
    size = 2.5,
    stroke = 0.6
  ) +
  
  facet_wrap(
    ~nutrient,
    ncol = 1,
    scales = "free_y"
  ) +
  
  scale_color_manual(
    values = c(
      "Source" = "#C9554B",
      "Target #1" = "#247F9E",
      "Target #2" = "#759BAC"
    ),
    name = "Station"
  ) +
  
  scale_x_continuous(
    breaks = -1:MAX_LAG,
    labels = c(
      "-1", "0", "1", "2", "3"
    )
  ) +
  
  labs(
    title = "(f) Representative observed anomaly trajectories",
    subtitle = "One TN and one TP event; source and two highest-ranked targets",
    x = "Months relative to trigger",
    y = "Observed concentration anomaly (Z)"
  ) +
  
  theme_fig(9) +
  
  theme(
    legend.position = "bottom"
  )


# ================================================================
# 18. BUILD FINAL FIGURE
# ================================================================

# Top row: two full-basin maps
# Middle row: performance + attention matrix
# Bottom row: PCS verification + time series

fig5 <- (
  (p_a | p_b) /
    (p_c | p_d) /
    (p_e | p_f)
) +
  
  patchwork::plot_layout(
    heights = c(1.10, 1.05, 1.0),
    widths = c(1, 1.45)
  ) +
  
  patchwork::plot_annotation(
    title = paste0(
      "Anomaly-triggered complementary monitoring ",
      "prioritization using source-conditioned F-Net"
    ),
    
    subtitle = paste0(
      "Three north–south spatial regions | ",
      "One observed TN and TP anomaly event per region | ",
      "Local Top-5 monitoring targets"
    ),
    
    caption = paste0(
      "Regional colors in (a–b) identify monitoring stations, ",
      "not watershed partitions. ",
      "Source-conditioned PCS provides relative attention ranks; ",
      "dual-lag rules provide indicative monitoring windows. ",
      "Observed anomalies do not establish causal nutrient transport. ",
      "Retrospective analysis using a full-period F-Net."
    ),
    
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 15
      ),
      
      plot.subtitle = element_text(
        size = 10,
        color = "#626E75"
      ),
      
      plot.caption = element_text(
        size = 8,
        color = "#626E75",
        hjust = 0
      ),
      
      plot.margin = margin(
        10, 12, 10, 12
      )
    )
  )


# ================================================================
# 19. EXPORT COMPLETE FIGURE
# ================================================================

ggsave(
  filename = file.path(
    OUT, "Fig5_main_final.pdf"
  ),
  plot = fig5,
  width = FIG_WIDTH,
  height = FIG_HEIGHT,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = file.path(
    OUT, "Fig5_main_final.png"
  ),
  plot = fig5,
  width = FIG_WIDTH,
  height = FIG_HEIGHT,
  dpi = 400,
  bg = "white",
  limitsize = FALSE
)


# ================================================================
# 20. EXPORT INDIVIDUAL PANELS
# ================================================================

panel_list <- list(
  a = p_a,
  b = p_b,
  c = p_c,
  d = p_d,
  e = p_e,
  f = p_f
)

for (nm in names(panel_list)) {
  
  ggsave(
    filename = file.path(
      OUT,
      paste0("Fig5_panel_", nm, ".pdf")
    ),
    plot = panel_list[[nm]],
    width = if (nm %in% c("d", "e")) 11 else 8,
    height = if (nm %in% c("d", "e")) 7 else 6,
    bg = "white"
  )
}


# ================================================================
# 21. EXPORT FIGURE DATA
# ================================================================

save_csv(
  top5,
  "04_Fig5_all_Top5_targets.csv"
)

save_csv(
  attention_grid,
  "05_Fig5_attention_matrix.csv"
)

save_csv(
  results,
  "06_Fig5_target_verification.csv"
)

save_csv(
  rep_series,
  "07_Fig5_representative_timeseries.csv"
)


# ================================================================
# 22. FINAL QA
# ================================================================

cat("\n========================================\n")
cat("FIGURE 5 COMPLETED\n")
cat("========================================\n")

cat("\nSix selected cases:\n")

print(
  cases %>%
    select(
      region,
      nutrient,
      event_id,
      source,
      trigger_date,
      source_Z
    ),
  n = Inf
)

cat("\nSpatial source and target counts:\n")

print(
  seg_df %>%
    count(region, nutrient),
  n = Inf
)

cat("\nTop-5 performance:\n")

print(
  performance %>%
    select(
      region_short,
      nutrient,
      n_recommended,
      n_evaluable,
      n_hits,
      hit_rate,
      coverage
    ),
  n = Inf
)

cat("\nRecommended time window distribution:\n")

print(
  attention_grid %>%
    distinct(
      event_id,
      target,
      tile_status,
      lag_months
    ) %>%
    count(
      tile_status
    ),
  n = Inf
)

cat("\nOutput directory:\n")
cat(OUT, "\n")

print(fig5)
