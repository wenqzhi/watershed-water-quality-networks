
# ================================================================
# FIGURE 5 — OPTIMIZED FINAL VERSION
#
# Anomaly-triggered monitoring prioritization
# using source-conditioned F-Net
#
# Panels:
# (a) Observed anomaly trigger
# (b) Source-conditioned Top-5 targets
# (c) Representative network associations
# (d) Attention windows + observed anomalies
# (e) Top-5 monitoring verification table
# (f) Observed nutrient anomaly trajectories
#
# No recalculation of F-Net, PCS or dual lags.
# ================================================================


# ================================================================
# 0. PACKAGES
# ================================================================

PKGS <- c(
  "readr", "dplyr", "tidyr", "tibble",
  "ggplot2", "sf", "patchwork",
  "scales", "lubridate", "stringr"
)

for (p in PKGS) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

invisible(lapply(PKGS, library, character.only = TRUE))

sf::sf_use_s2(FALSE)


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

STEP52 <- file.path(
  BASE, "Step5_2_source_PCS"
)

STEP53 <- file.path(
  BASE, "Step5_3_dual_lag_TopK"
)

OUT <- file.path(
  BASE, "Fig5_optimized"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

EVENT_FILE <- file.path(
  STEP51, "05_all_real_anomaly_triggers.csv"
)

CATALOG_FILE <- file.path(
  STEP51, "08_Fig5_candidate_event_catalog.csv"
)

OBS_FILE <- file.path(
  STEP51, "06_event_target_observation_matrix.csv"
)

ANOMALY_FILE <- file.path(
  STEP51, "04_all_monthly_concentration_anomalies.csv"
)

RANK_FILE <- file.path(
  STEP52, "03_all_event_target_rankings.csv"
)

WINDOW_FILE <- file.path(
  STEP53, "04_TopK_dual_lag_windows.csv"
)

VERIFY_FILE <- file.path(
  STEP53, "06_TopK_target_verification.csv"
)

SUMMARY_FILE <- file.path(
  STEP53, "08_TopK_dual_lag_summary.csv"
)

STATION_FILE <- file.path(
  ROOT, "NTN-station.txt"
)

BASIN_FILE <- file.path(
  ROOT, "ntnwatersheds",
  "ntnwatersheds.shp"
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

# NULL = automatically select using lag coverage,
# observation coverage, and source anomaly strength.
#
# Specify "EV0136" to reproduce the previous example.
SELECT_EVENT <- "EV0136"

MAP_CRS <- 5070

FIG_WIDTH <- 15.5
FIG_HEIGHT <- 15.2

# Colors
COL_SOURCE <- "#C94E4A"
COL_TOP <- "#247C9B"
COL_OTHER <- "#ABB5BA"

COL_PHYSICAL <- "#287F9F"
COL_STATISTICAL <- "#8368AD"
COL_MIXED <- "#D49A43"
COL_UNVERIFIED <- "#B28AC7"
COL_UNASSIGNED <- "#969FA5"

COL_WINDOW <- "#3D8BA7"
COL_OUTSIDE <- "#EBEFF1"
COL_HIT <- "#2F927A"
COL_MISS <- "#CF765A"
COL_UNKNOWN <- "#99A3A8"

ROUTE_COLORS <- c(
  "Physical" = COL_PHYSICAL,
  "Statistical" = COL_STATISTICAL,
  "Mixed" = COL_MIXED,
  "TE: unverified" = COL_UNVERIFIED,
  "Unassigned" = COL_UNASSIGNED
)

OUTCOME_COLORS <- c(
  "Hit" = COL_HIT,
  "No anomaly" = COL_MISS,
  "Not evaluable" = COL_UNKNOWN
)


# ================================================================
# 3. HELPERS
# ================================================================

clean_id <- function(x) {
  
  x <- trimws(as.character(x))
  x <- sub("\\.0+$", "", x)
  
  x[x %in% c("", "NA", "NaN")] <- NA_character_
  
  # Restore leading zeros in short numeric USGS IDs.
  ii <- !is.na(x) & grepl("^[0-9]{1,7}$", x)
  
  x[ii] <- sprintf(
    "%08.0f",
    as.numeric(x[ii])
  )
  
  x
}

read_ids <- function(path, ids = character()) {
  
  if (!file.exists(path)) {
    stop("Missing input file:\n", path)
  }
  
  spec <- readr::cols(
    .default = readr::col_guess()
  )
  
  for (nm in ids) {
    spec$cols[[nm]] <- readr::col_character()
  }
  
  readr::read_csv(
    path,
    col_types = spec,
    show_col_types = FALSE
  )
}

as_bool <- function(x) {
  
  tolower(as.character(x)) %in%
    c("true", "t", "1", "yes")
}

save_csv <- function(d, filename) {
  
  readr::write_csv(
    d,
    file.path(OUT, filename),
    na = ""
  )
}

theme_fig <- function(base_size = 10) {
  
  theme_bw(base_size = base_size) +
    
    theme(
      panel.grid.minor = element_blank(),
      
      panel.grid.major = element_line(
        color = "#E7EBED",
        linewidth = 0.2
      ),
      
      plot.title = element_text(
        face = "bold",
        size = base_size + 0.5,
        margin = margin(b = 5)
      ),
      
      plot.subtitle = element_text(
        color = "#68737A",
        size = base_size - 1,
        margin = margin(b = 6)
      ),
      
      strip.background = element_rect(
        fill = "#F3F5F6",
        color = "#DCE1E3"
      ),
      
      strip.text = element_text(
        face = "bold"
      ),
      
      legend.position = "bottom",
      
      legend.title = element_text(
        face = "bold",
        size = 9
      ),
      
      legend.text = element_text(
        size = 8
      ),
      
      plot.margin = margin(8, 8, 8, 8)
    )
}


# ================================================================
# 4. READ STEP 5.1–5.3 RESULTS
# ================================================================

events <- read_ids(
  EVENT_FILE,
  "source"
) %>%
  mutate(
    source = clean_id(source),
    nutrient = toupper(nutrient),
    trigger_date = as.Date(trigger_date)
  )

catalog <- read_ids(
  CATALOG_FILE,
  "source"
) %>%
  mutate(
    source = clean_id(source),
    nutrient = toupper(nutrient),
    trigger_date = as.Date(trigger_date)
  )

obs <- read_ids(
  OBS_FILE,
  c("source", "target")
) %>%
  mutate(
    source = clean_id(source),
    target = clean_id(target),
    nutrient = toupper(nutrient),
    target_date = as.Date(target_date),
    lag_months = as.integer(lag_months),
    
    target_observed = as_bool(target_observed),
    target_anomaly_15 = as_bool(target_anomaly_15)
  )

anomaly <- read_ids(
  ANOMALY_FILE,
  "STAID"
) %>%
  mutate(
    STAID = clean_id(STAID),
    nutrient = toupper(nutrient),
    date = as.Date(date),
    Z = as.numeric(Z)
  )

rankings <- read_ids(
  RANK_FILE,
  c("source", "target")
) %>%
  mutate(
    source = clean_id(source),
    target = clean_id(target),
    nutrient = toupper(nutrient),
    rank_FNet = as.integer(rank_FNet),
    source_PCS = as.numeric(source_PCS)
  )

windows <- read_ids(
  WINDOW_FILE,
  c("source", "target")
) %>%
  mutate(
    source = clean_id(source),
    target = clean_id(target),
    nutrient = toupper(nutrient),
    rank_FNet = as.integer(rank_FNet),
    
    window_defined = as_bool(window_defined),
    
    lag_min = as.integer(lag_min),
    lag_max = as.integer(lag_max),
    
    path_months = as.numeric(path_months)
  )

verification <- read_ids(
  VERIFY_FILE,
  c("source", "target")
) %>%
  mutate(
    source = clean_id(source),
    target = clean_id(target),
    nutrient = toupper(nutrient),
    
    observed_anomaly = as_bool(observed_anomaly),
    
    window_fully_observed = as_bool(
      window_fully_observed
    )
  )

summary53 <- read_ids(SUMMARY_FILE)

if (
  anyDuplicated(
    windows[c("event_id", "target")]
  ) > 0
) {
  stop("Duplicate event-target records in windows.")
}

if (
  anyDuplicated(
    obs[c("event_id", "target", "lag_months")]
  ) > 0
) {
  stop("Duplicate event-target-month observations.")
}


# ================================================================
# 5. SELECT REPRESENTATIVE EVENT
# ================================================================

event_quality <- windows %>%
  filter(rank_FNet <= TOP_K) %>%
  group_by(event_id) %>%
  summarise(
    n_assigned = sum(
      window_defined &
        is.finite(lag_min) &
        lag_min <= MAX_LAG,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  left_join(
    catalog %>%
      select(
        event_id, nutrient,
        source, trigger_date,
        source_Z, n_targets_observed
      ),
    by = "event_id"
  ) %>%
  arrange(
    desc(n_assigned),
    desc(n_targets_observed),
    desc(source_Z),
    event_id
  )

save_csv(
  event_quality,
  "01_event_selection_QA.csv"
)

if (is.null(SELECT_EVENT)) {
  SELECT_EVENT <- event_quality$event_id[1]
}

selected_event <- events %>%
  filter(event_id == SELECT_EVENT)

if (nrow(selected_event) != 1) {
  stop("Selected event is missing or duplicated.")
}

NU <- selected_event$nutrient[1]
SOURCE <- selected_event$source[1]
TRIGGER_DATE <- selected_event$trigger_date[1]
SOURCE_Z <- as.numeric(selected_event$source_Z[1])

cat("\nSelected event:\n")
print(selected_event)


# ================================================================
# 6. PREPARE TOP-5
# ================================================================

top5 <- windows %>%
  filter(
    event_id == SELECT_EVENT,
    rank_FNet <= TOP_K
  ) %>%
  left_join(
    verification %>%
      select(
        event_id,
        target,
        window_fully_observed,
        observed_anomaly,
        verification_status
      ),
    by = c("event_id", "target")
  ) %>%
  arrange(rank_FNet) %>%
  mutate(
    station_label = paste0(
      "#", rank_FNet, " ", target
    ),
    
    route_class = case_when(
      route_type == "Physical" ~ "Physical",
      
      route_type == "Statistical" ~ "Statistical",
      
      route_type == "Mixed" ~ "Mixed",
      
      route_type == "TE-unverified" ~
        "TE: unverified",
      
      TRUE ~ "Unassigned"
    ),
    
    outcome = case_when(
      window_fully_observed & observed_anomaly ~
        "Hit",
      
      window_fully_observed & !observed_anomaly ~
        "No anomaly",
      
      TRUE ~ "Not evaluable"
    ),
    
    outcome = factor(
      outcome,
      levels = c(
        "Hit", "No anomaly", "Not evaluable"
      )
    ),
    
    display_window = case_when(
      !window_defined ~ "Unassigned",
      
      is.na(lag_min) | is.na(lag_max) ~
        "Unassigned",
      
      lag_min > MAX_LAG ~ ">3 mo",
      
      lag_min == lag_max ~
        paste0(lag_min, " mo"),
      
      TRUE ~ paste0(
        lag_min, "–", lag_max, " mo"
      )
    )
  )

if (nrow(top5) != TOP_K) {
  stop("Top-5 data incomplete.")
}

save_csv(
  top5,
  "02_selected_Top5_results.csv"
)


# ================================================================
# 7. STATION COORDINATES
# ================================================================

# Important:
# NTN-station.txt contains comma-separated:
# FID,SITE_NO,lon,lat

stations_raw <- readr::read_csv(
  STATION_FILE,
  col_types = readr::cols(
    .default = readr::col_character()
  ),
  show_col_types = FALSE
)

names(stations_raw) <- trimws(
  names(stations_raw)
)

required_cols <- c(
  "SITE_NO", "lon", "lat"
)

if (
  !all(required_cols %in% names(stations_raw))
) {
  stop(
    "Missing station columns: ",
    paste(
      setdiff(required_cols, names(stations_raw)),
      collapse = ", "
    )
  )
}

stations <- stations_raw %>%
  transmute(
    STAID = clean_id(SITE_NO),
    
    lon = suppressWarnings(
      as.numeric(lon)
    ),
    
    lat = suppressWarnings(
      as.numeric(lat)
    )
  ) %>%
  filter(
    !is.na(STAID),
    is.finite(lon),
    is.finite(lat),
    abs(lon) <= 180,
    abs(lat) <= 90
  ) %>%
  distinct(
    STAID,
    .keep_all = TRUE
  )

stations_sf <- stations %>%
  sf::st_as_sf(
    coords = c("lon", "lat"),
    crs = 4326,
    remove = FALSE
  ) %>%
  sf::st_transform(MAP_CRS)


# ================================================================
# 8. BASIN AND RIVER NETWORK
# ================================================================

read_map <- function(path, label) {
  
  if (!file.exists(path)) {
    stop(label, " missing: ", path)
  }
  
  x <- sf::st_read(
    path,
    quiet = TRUE
  )
  
  if (is.na(sf::st_crs(x))) {
    
    bb <- sf::st_bbox(x)
    
    plausible_lonlat <- (
      bb["xmin"] >= -180 &&
        bb["xmax"] <= 180 &&
        bb["ymin"] >= -90 &&
        bb["ymax"] <= 90
    )
    
    if (!plausible_lonlat) {
      stop(
        label,
        " CRS missing and coordinates not lon/lat."
      )
    }
    
    warning(
      label,
      ": CRS missing; assuming EPSG:4326."
    )
    
    sf::st_crs(x) <- 4326
  }
  
  sf::st_transform(
    x, MAP_CRS
  )
}

basin <- read_map(
  BASIN_FILE, "Basin"
)

rivers <- read_map(
  RIVER_FILE, "Rivers"
)

bbox <- sf::st_bbox(basin)

XLIM <- as.numeric(
  bbox[c("xmin", "xmax")]
)

YLIM <- as.numeric(
  bbox[c("ymin", "ymax")]
)


# ================================================================
# 9. SELECTED EVENT SPATIAL DATA
# ================================================================

all_targets <- rankings %>%
  filter(
    event_id == SELECT_EVENT
  ) %>%
  select(
    STAID = target,
    source_PCS,
    rank_FNet
  ) %>%
  left_join(
    stations,
    by = "STAID"
  ) %>%
  filter(
    is.finite(lon),
    is.finite(lat)
  ) %>%
  sf::st_as_sf(
    coords = c("lon", "lat"),
    crs = 4326
  ) %>%
  sf::st_transform(MAP_CRS)

source_sf <- stations_sf %>%
  filter(STAID == SOURCE)

top5_sf <- stations_sf %>%
  inner_join(
    top5 %>%
      select(
        STAID = target,
        rank_FNet,
        source_PCS,
        route_class,
        outcome
      ),
    by = "STAID"
  )

if (nrow(source_sf) != 1) {
  stop("Source station coordinates missing.")
}

if (nrow(top5_sf) != TOP_K) {
  stop("Top-5 station coordinates missing.")
}

# Schematic straight connections only.
# These are NOT actual traced river routes.

source_xy <- sf::st_coordinates(
  source_sf
)[1, ]

top_xy <- sf::st_coordinates(top5_sf)

lines_df <- top5_sf %>%
  sf::st_drop_geometry() %>%
  mutate(
    x1 = source_xy[1],
    y1 = source_xy[2],
    x2 = top_xy[, 1],
    y2 = top_xy[, 2]
  )


# ================================================================
# 10. MAP STYLE
# ================================================================

map_base <- function() {
  
  list(
    geom_sf(
      data = basin,
      fill = "#FBFCFC",
      color = "#D5DADC",
      linewidth = 0.16
    ),
    
    geom_sf(
      data = rivers,
      color = "#AFD0DC",
      linewidth = 0.20,
      alpha = 0.83
    ),
    
    coord_sf(
      xlim = XLIM,
      ylim = YLIM,
      datum = NA,
      expand = FALSE
    ),
    
    theme_fig(9),
    
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      
      legend.position = "bottom"
    )
  )
}


# ================================================================
# 11. PANEL A — OBSERVED ANOMALY SOURCE
# ================================================================

p_a <- ggplot() +
  
  map_base() +
  
  geom_sf(
    data = all_targets,
    color = "#A6B0B6",
    size = 1.4,
    alpha = 0.8
  ) +
  
  geom_sf(
    data = source_sf,
    shape = 8,
    color = COL_SOURCE,
    size = 4.3,
    stroke = 1.15
  ) +
  
  labs(
    title = "(a) Observed anomaly trigger",
    subtitle = paste0(
      "Source ", SOURCE,
      " | Z = ",
      sprintf("%.2f", SOURCE_Z)
    )
  )


# ================================================================
# 12. PANEL B — SOURCE-CONDITIONED PCS
# ================================================================

p_b <- ggplot() +
  
  map_base() +
  
  geom_sf(
    data = all_targets,
    aes(size = source_PCS),
    shape = 21,
    fill = "#B7D6E0",
    color = "#53889A",
    alpha = 0.82,
    stroke = 0.3
  ) +
  
  geom_sf(
    data = top5_sf,
    shape = 21,
    fill = COL_TOP,
    color = "white",
    size = 4.2,
    stroke = 1.0
  ) +
  
  geom_sf_text(
    data = top5_sf,
    aes(label = rank_FNet),
    color = "white",
    fontface = "bold",
    size = 2.5
  ) +
  
  geom_sf(
    data = source_sf,
    shape = 8,
    color = COL_SOURCE,
    size = 3.6
  ) +
  
  scale_size_continuous(
    range = c(0.7, 3.1),
    name = "Source-conditioned PCS"
  ) +
  
  labs(
    title = "(b) Source-conditioned Top-5",
    subtitle = "Numbered nodes indicate monitoring rank"
  )


# ================================================================
# 13. PANEL C — REPRESENTATIVE CONNECTIONS
# ================================================================

p_c <- ggplot() +
  
  map_base() +
  
  geom_segment(
    data = lines_df %>%
      filter(route_class != "Unassigned"),
    
    aes(
      x = x1,
      y = y1,
      xend = x2,
      yend = y2,
      color = route_class
    ),
    
    linetype = "dashed",
    linewidth = 0.72,
    alpha = 0.85
  ) +
  
  geom_sf(
    data = top5_sf,
    aes(fill = route_class),
    shape = 21,
    color = "white",
    size = 3.9,
    stroke = 1.0
  ) +
  
  geom_sf_text(
    data = top5_sf,
    aes(label = rank_FNet),
    color = "white",
    fontface = "bold",
    size = 2.3
  ) +
  
  geom_sf(
    data = source_sf,
    shape = 8,
    color = COL_SOURCE,
    size = 3.7
  ) +
  
  scale_color_manual(
    values = ROUTE_COLORS,
    drop = TRUE,
    name = "Connection type"
  ) +
  
  scale_fill_manual(
    values = ROUTE_COLORS,
    drop = TRUE,
    guide = "none"
  ) +
  
  labs(
    title = "(c) Representative network links",
    subtitle = "Straight dashed links are schematic, not river routes"
  )


# ================================================================
# 14. PANEL D — ATTENTION + ACTUAL ANOMALIES
# ================================================================

# The panel uses a complete month-by-target grid.
# Actual anomaly detections are shown as red dots
# regardless of whether the month was recommended.
#
# This avoids hiding observed anomalies that fall
# outside the attention window.

attention_grid <- tidyr::crossing(
  target = top5$target,
  lag_months = 0:MAX_LAG
) %>%
  left_join(
    top5 %>%
      select(
        target,
        rank_FNet,
        station_label,
        lag_min,
        lag_max,
        window_defined
      ),
    by = "target"
  ) %>%
  left_join(
    obs %>%
      filter(event_id == SELECT_EVENT) %>%
      select(
        target,
        lag_months,
        target_observed,
        target_anomaly_15,
        target_Z
      ),
    by = c("target", "lag_months")
  ) %>%
  mutate(
    attention = (
      window_defined &
        !is.na(lag_min) &
        !is.na(lag_max) &
        lag_months >= lag_min &
        lag_months <= lag_max
    ),
    
    tile_status = case_when(
      !window_defined ~ "Unassigned",
      attention ~ "Attention",
      TRUE ~ "Outside window"
    ),
    
    observed_anomaly = (
      target_observed %in% TRUE &
        target_anomaly_15 %in% TRUE
    ),
    
    station_label = factor(
      station_label,
      levels = rev(top5$station_label)
    )
  )

p_d <- ggplot(
  attention_grid,
  aes(
    x = lag_months,
    y = station_label
  )
) +
  
  geom_tile(
    aes(fill = tile_status),
    width = 0.97,
    height = 0.92,
    color = "white",
    linewidth = 0.9
  ) +
  
  # Real observed anomaly:
  # red circle over the corresponding month.
  geom_point(
    data = attention_grid %>%
      filter(observed_anomaly),
    
    shape = 21,
    fill = COL_SOURCE,
    color = "white",
    size = 3.5,
    stroke = 0.7
  ) +
  
  # Missing observations, if any.
  geom_point(
    data = attention_grid %>%
      filter(
        is.na(target_observed) |
          !target_observed
      ),
    
    shape = 4,
    color = "#555F66",
    size = 2,
    stroke = 0.65
  ) +
  
  scale_fill_manual(
    values = c(
      "Attention" = COL_WINDOW,
      "Outside window" = COL_OUTSIDE,
      "Unassigned" = "#BCC4C8"
    ),
    name = "Attention status"
  ) +
  
  scale_x_continuous(
    breaks = 0:MAX_LAG,
    labels = paste0(0:MAX_LAG, " mo"),
    expand = expansion(add = 0.5)
  ) +
  
  labs(
    title = "(d) Attention windows and actual anomalies",
    subtitle = "Red dots: observed Z > 1.5; crosses: missing observations",
    x = "Months after source anomaly",
    y = NULL
  ) +
  
  theme_fig(10) +
  
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 9),
    legend.position = "bottom"
  )


# ================================================================
# 15. PANEL E — TOP-5 RESULT TABLE
# ================================================================

# Use a readable table-style ggplot instead of
# repeating PCS as another bar chart.

table_data <- top5 %>%
  arrange(rank_FNet) %>%
  mutate(
    row_id = TOP_K + 1L - rank_FNet,
    
    PCS_text = sprintf(
      "%.3f", source_PCS
    ),
    
    window_text = display_window,
    
    outcome_text = case_when(
      outcome == "Hit" ~ "Hit",
      outcome == "No anomaly" ~ "No anomaly",
      TRUE ~ "Not evaluable"
    )
  )

# Draw each column with aligned positions.
# x=1: rank
# x=2: station
# x=3: PCS
# x=4: window
# x=5: observed

p_e <- ggplot() +
  
  geom_rect(
    data = table_data,
    
    aes(
      xmin = 0.55,
      xmax = 5.5,
      ymin = row_id - 0.47,
      ymax = row_id + 0.47
    ),
    
    fill = rep(
      c("#F5F7F8", "#FFFFFF"),
      length.out = nrow(table_data)
    ),
    
    color = "#E5E9EB",
    linewidth = 0.3
  ) +
  
  # Column headings
  annotate(
    "text",
    x = c(1, 2, 3, 4, 5),
    y = TOP_K + 0.95,
    
    label = c(
      "Rank", "Station",
      "PCS", "Window", "Observed"
    ),
    
    fontface = "bold",
    size = 3.3,
    color = "#354047"
  ) +
  
  geom_text(
    data = table_data,
    aes(
      x = 1,
      y = row_id,
      label = paste0("#", rank_FNet)
    ),
    fontface = "bold",
    size = 3.5
  ) +
  
  geom_text(
    data = table_data,
    aes(
      x = 2,
      y = row_id,
      label = target
    ),
    size = 3.25
  ) +
  
  geom_text(
    data = table_data,
    aes(
      x = 3,
      y = row_id,
      label = PCS_text
    ),
    size = 3.3
  ) +
  
  geom_text(
    data = table_data,
    aes(
      x = 4,
      y = row_id,
      label = window_text
    ),
    size = 3.3
  ) +
  
  geom_point(
    data = table_data,
    aes(
      x = 4.62,
      y = row_id,
      fill = outcome
    ),
    shape = 21,
    size = 3,
    color = "white",
    stroke = 0.5
  ) +
  
  geom_text(
    data = table_data,
    aes(
      x = 5.12,
      y = row_id,
      label = outcome_text
    ),
    hjust = 0.5,
    size = 2.65
  ) +
  
  scale_fill_manual(
    values = OUTCOME_COLORS,
    guide = "none"
  ) +
  
  scale_x_continuous(
    limits = c(0.5, 5.65),
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    limits = c(0.4, TOP_K + 1.45),
    expand = c(0, 0)
  ) +
  
  labs(
    title = "(e) Top-5 monitoring verification",
    subtitle = paste0(
      sum(
        top5$outcome == "Hit",
        na.rm = TRUE
      ),
      "/",
      sum(
        top5$outcome != "Not evaluable",
        na.rm = TRUE
      ),
      " evaluable targets with an observed anomaly"
    )
  ) +
  
  theme_void(base_size = 10) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 10.5,
      margin = margin(b = 5)
    ),
    
    plot.subtitle = element_text(
      color = "#68737A",
      size = 9
    ),
    
    plot.margin = margin(8, 8, 8, 8)
  )


# ================================================================
# 16. PANEL F — OBSERVED ANOMALY SERIES
# ================================================================

START_PLOT <- TRIGGER_DATE %m-% months(1)
END_PLOT <- TRIGGER_DATE %m+% months(MAX_LAG)

series_ids <- c(
  SOURCE,
  top5$target
)

series_order <- c(
  paste0("Source ", SOURCE),
  paste0(
    "#", top5$rank_FNet,
    " ", top5$target
  )
)

series <- anomaly %>%
  filter(
    nutrient == NU,
    STAID %in% series_ids,
    date >= START_PLOT,
    date <= END_PLOT
  ) %>%
  mutate(
    station_label = case_when(
      STAID == SOURCE ~
        paste0("Source ", STAID),
      
      TRUE ~ paste0(
        "#",
        match(STAID, top5$target),
        " ",
        STAID
      )
    ),
    
    station_label = factor(
      station_label,
      levels = series_order
    ),
    
    station_role = if_else(
      STAID == SOURCE,
      "Source",
      "Top-5 target"
    )
  )

# Actual recommended shading for each target
# during the corresponding month.

attention_bands <- attention_grid %>%
  filter(attention) %>%
  mutate(
    start_date = TRIGGER_DATE %m+%
      months(lag_months),
    
    end_date = start_date %m+% months(1),
    
    station_label = factor(
      as.character(station_label),
      levels = series_order
    )
  )

# Ensure the source is always the first facet.
# All facets use identical x-axis limits.

p_f <- ggplot() +
  
  geom_rect(
    data = attention_bands,
    aes(
      xmin = start_date,
      xmax = end_date,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = "#D8ECF1",
    alpha = 0.7,
    inherit.aes = FALSE
  ) +
  
  geom_hline(
    yintercept = Z_THRESHOLD,
    linetype = "dashed",
    linewidth = 0.55,
    color = COL_SOURCE
  ) +
  
  geom_vline(
    xintercept = as.numeric(TRIGGER_DATE),
    linetype = "dotted",
    linewidth = 0.55,
    color = "#7B878C"
  ) +
  
  geom_line(
    data = series,
    aes(
      x = date,
      y = Z,
      group = station_label,
      color = station_role
    ),
    linewidth = 0.9,
    na.rm = TRUE
  ) +
  
  geom_point(
    data = series,
    aes(
      x = date,
      y = Z,
      color = station_role
    ),
    size = 1.9,
    na.rm = TRUE
  ) +
  
  # Show observed threshold crossings.
  geom_point(
    data = series %>%
      filter(
        is.finite(Z),
        Z > Z_THRESHOLD
      ),
    
    aes(
      x = date,
      y = Z
    ),
    
    shape = 21,
    fill = COL_SOURCE,
    color = "white",
    stroke = 0.7,
    size = 3.0
  ) +
  
  facet_wrap(
    ~station_label,
    ncol = 2,
    scales = "fixed",
    drop = FALSE
  ) +
  
  scale_color_manual(
    values = c(
      "Source" = COL_SOURCE,
      "Top-5 target" = COL_TOP
    ),
    guide = "none"
  ) +
  
  scale_x_date(
    limits = c(
      START_PLOT,
      END_PLOT
    ),
    date_breaks = "1 month",
    date_labels = "%b",
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  
  labs(
    title = "(f) Observed nutrient anomalies",
    subtitle = paste0(
      "Shaded: attention windows; dashed: Z = ",
      Z_THRESHOLD
    ),
    x = "Month",
    y = "Observed concentration anomaly (Z)"
  ) +
  
  theme_fig(9) +
  
  theme(
    strip.text = element_text(size = 8.5),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    
    panel.grid.minor = element_blank()
  )


# ================================================================
# 17. COMBINE FIGURE
# ================================================================

# Preserve two equal-width panels per row.
# Collecting all legends is avoided here because
# panel-specific legends can become overcrowded.

fig5 <- (
  (p_a | p_b) /
    (p_c | p_d) /
    (p_e | p_f)
) +
  
  patchwork::plot_layout(
    heights = c(1.08, 1.08, 1.18)
  ) +
  
  patchwork::plot_annotation(
    title = paste0(
      "Anomaly-triggered monitoring prioritization ",
      "using source-conditioned F-Net"
    ),
    
    subtitle = paste0(
      NU, " anomaly | ",
      "Event ", SELECT_EVENT,
      " | Source ", SOURCE,
      " | ",
      format(TRIGGER_DATE, "%B %Y")
    ),
    
    caption = paste0(
      "Top-5 stations are selected by source-conditioned PCS. ",
      "Attention windows use available dual-lag information; ",
      "red markers denote observed concentration anomalies. ",
      "Schematic links are not river routes. ",
      "Retrospective assessment, not an independent forecast."
    ),
    
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 15
      ),
      
      plot.subtitle = element_text(
        color = "#647078",
        size = 10
      ),
      
      plot.caption = element_text(
        color = "#647078",
        size = 8,
        hjust = 0
      ),
      
      plot.margin = margin(10, 12, 10, 12)
    )
  )


# ================================================================
# 18. EXPORT
# ================================================================

ggsave(
  file.path(
    OUT, "Fig5_optimized.pdf"
  ),
  plot = fig5,
  width = FIG_WIDTH,
  height = FIG_HEIGHT,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  file.path(
    OUT, "Fig5_optimized.png"
  ),
  plot = fig5,
  width = FIG_WIDTH,
  height = FIG_HEIGHT,
  dpi = 400,
  bg = "white",
  limitsize = FALSE
)


# ================================================================
# 19. EXPORT INDIVIDUAL PANELS
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
    file.path(
      OUT, paste0("Fig5_", nm, ".pdf")
    ),
    plot = panel_list[[nm]],
    width = 7.5,
    height = 5.5,
    bg = "white"
  )
}


# ================================================================
# 20. EXPORT SELECTED EVENT RESULT TABLE
# ================================================================

result_table <- top5 %>%
  select(
    event_id,
    nutrient,
    source,
    trigger_date,
    target,
    rank_FNet,
    source_PCS,
    route_class,
    representative_route,
    path_months,
    window_label,
    display_window,
    window_defined,
    outcome,
    verification_status
  )

save_csv(
  result_table,
  "03_Fig5_optimized_event_table.csv"
)

save_csv(
  attention_grid,
  "04_Fig5_monthly_attention_and_observation.csv"
)


# ================================================================
# 21. FINAL REPORT
# ================================================================

cat("\n========================================\n")
cat("FIGURE 5 — OPTIMIZED VERSION FINISHED\n")
cat("========================================\n")

cat("\nSelected event:\n")
print(selected_event)

cat("\nTop-5 result:\n")
print(
  result_table,
  n = Inf,
  width = Inf
)

cat("\nSummary:\n")

print(
  top5 %>%
    summarise(
      n_recommended = n(),
      
      n_assigned = sum(
        window_defined,
        na.rm = TRUE
      ),
      
      n_evaluable = sum(
        outcome != "Not evaluable",
        na.rm = TRUE
      ),
      
      n_hits = sum(
        outcome == "Hit",
        na.rm = TRUE
      )
    )
)

cat("\nFull Step 5.3 summary:\n")
print(summary53, n = Inf)

cat("\nOutput directory:\n")
cat(OUT, "\n")

print(fig5)
