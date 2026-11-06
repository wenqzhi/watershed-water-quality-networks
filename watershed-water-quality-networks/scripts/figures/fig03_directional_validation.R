
# ============================================================
# FIGURE 3 — REVISED OBSERVATIONAL VALIDATION
#
# Agreement and Divergence between Physical River
# Connectivity and Directional Nutrient Dependencies
#
# (a) Selected river-network configurations
# (b) Observed TN / TP dynamics
# (c) Directional TE strength and significance
# (d) Physical connectivity × TE significance
# (e) Lagged correspondence of nutrient anomalies
#
# 3 groups | 9 stations | 9 pairs | TN / TP
# 2011-10 to 2023-09
#
# New Panel E:
# - Uses measured concentrations directly
# - Seasonal anomalies (calendar-month z-scores)
# - Upstream -> downstream only
# - Lags: 0, 1, 2, 3 months
# - Conditional event probability minus background
#
# No interpolation
# No new TE calculation
# No TE lag used to choose event lag
# ============================================================


# ============================================================
# 0. PACKAGES
# ============================================================

PKGS <- c(
  "readr",
  "dplyr",
  "tidyr",
  "purrr",
  "stringr",
  "tibble",
  "ggplot2",
  "sf",
  "patchwork",
  "scales"
)

for (p in PKGS) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

invisible(
  lapply(PKGS, library, character.only = TRUE)
)


# ============================================================
# 1. PATHS
# ============================================================

source("config.R")
ROOT <- get_data_root()

CONC_FILE <- file.path(
  ROOT,
  "NTN_2023_MonthlyLoadTable.csv"
)

STATION_FILE <- file.path(
  ROOT,
  "NTN-station.txt"
)

RIVER_FILE <- file.path(
  ROOT,
  "水系/水系.shp"
)

BASIN_FILE <- file.path(
  ROOT,
  "ntnwatersheds/ntnwatersheds.shp"
)

TE_FILE <- file.path(
  ROOT,
  "Fig3_revised",
  "04_all_36_directed_TE_results.csv"
)

OUT <- file.path(
  ROOT,
  "Fig3_observational_validation"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. SETTINGS
# ============================================================

START_DATE <- as.Date("2011-10-01")
END_DATE   <- as.Date("2023-09-01")

CONC_FIELD <- "KConc"

TE_P_THRESHOLD <- 0.10

# Panel B: display last 36 months.
# Full record is used for standardized anomalies.
PLOT_MONTHS <- 36L

# Panel E: anomaly threshold
Z_THRESHOLD <- 1.0

# Sensitivity options can be examined separately:
# 0.5, 1.0, 1.5

LAGS <- 0:3

# Minimum high-anomaly upstream months
MIN_UPSTREAM_EVENTS <- 5L

# Standardization needs at least this many
# observations for each calendar month.
MIN_CLIMATOLOGY_YEARS <- 8L

FIG_DPI <- 400

# Station role colors
COL_A <- "#D97745"
COL_B <- "#9568AA"
COL_D <- "#2785A1"

# TE classification
COL_BOTH <- "#217F70"
COL_PHYS_ONLY <- "#8AC6AC"
COL_TE_ONLY <- "#9970B0"
COL_NEITHER <- "#D4DCE0"

# Anomaly probability colors
COL_NEG <- "#3B83A5"
COL_ZERO <- "#F7F7F7"
COL_POS <- "#C86D4B"

# TE continuous palette
TE_PALETTE <- c(
  "#F7F7F7",
  "#D7E9E8",
  "#A4D3CB",
  "#5EBAA9",
  "#238A79",
  "#075C50"
)


# ============================================================
# 3. HELPERS
# ============================================================

norm_id <- function(x) {
  
  x <- trimws(as.character(x))
  x <- sub("\\.0+$", "", x)
  
  stringr::str_pad(
    x,
    width = 8,
    pad = "0"
  )
}

check_file <- function(file) {
  
  if (!file.exists(file)) {
    stop("Missing file:\n", file)
  }
}

save_csv <- function(x, filename) {
  
  readr::write_csv(
    x,
    file.path(OUT, filename),
    na = ""
  )
}

zscore <- function(x) {
  
  n_ok <- sum(is.finite(x))
  s <- sd(x, na.rm = TRUE)
  
  if (
    n_ok < 2 ||
    !is.finite(s) ||
    s <= 0
  ) {
    return(rep(NA_real_, length(x)))
  }
  
  (
    x - mean(x, na.rm = TRUE)
  ) / s
}

theme_fig <- function(size = 10) {
  
  theme_bw(base_size = size) +
    
    theme(
      panel.grid.minor = element_blank(),
      
      strip.text = element_text(
        face = "bold",
        size = 10
      ),
      
      plot.title = element_text(
        face = "bold",
        size = 11
      ),
      
      plot.subtitle = element_text(
        size = 8.5,
        color = "grey35"
      ),
      
      legend.position = "bottom",
      
      legend.title = element_text(
        size = 8.5
      ),
      
      legend.text = element_text(
        size = 8
      )
    )
}


# ============================================================
# 4. THREE STATION GROUPS
# ============================================================

groups <- tibble::tribble(
  
  ~group, ~A, ~B, ~D,
  
  "G1",
  "01529500",
  "01562000",
  "01578310",
  
  "G2",
  "01599000",
  "01631000",
  "01646580",
  
  "G3",
  "02034000",
  "02024752",
  "02035000"
  
)

groups <- groups %>%
  mutate(
    across(c(A, B, D), norm_id)
  )

station_roles <- groups %>%
  pivot_longer(
    cols = c(A, B, D),
    names_to = "role",
    values_to = "STAID"
  )

# Each group:
#
# P1 A -> D, physical downstream
# P2 B -> D, physical downstream
# P3 A <-> B, cross-tributary
#
# River topology is user-confirmed.

pairs <- purrr::map_dfr(
  seq_len(nrow(groups)),
  function(i) {
    
    g <- groups[i, ]
    
    tibble(
      group = g$group,
      
      pair = paste0(
        g$group,
        "_P",
        1:3
      ),
      
      from = c(
        g$A,
        g$B,
        g$A
      ),
      
      to = c(
        g$D,
        g$D,
        g$B
      ),
      
      pair_type = c(
        "A-D",
        "B-D",
        "A-B"
      ),
      
      connection = c(
        "Physical",
        "Physical",
        "Cross-tributary"
      )
    )
  }
)

PAIR_LEVELS <- pairs$pair

save_csv(
  pairs,
  "01_selected_pairs.csv"
)


# ============================================================
# 5. READ EXISTING TE DATA
# ============================================================

check_file(TE_FILE)

te_raw <- readr::read_csv(
  TE_FILE,
  
  col_types = readr::cols(
    source = readr::col_character(),
    target = readr::col_character(),
    .default = readr::col_guess()
  ),
  
  show_col_types = FALSE
)

te_required <- c(
  "nutrient",
  "source",
  "target",
  "TE_raw",
  "p_value",
  "lag_months"
)

missing_te <- setdiff(
  te_required,
  names(te_raw)
)

if (length(missing_te)) {
  
  stop(
    "Missing TE columns: ",
    paste(
      missing_te,
      collapse = ", "
    )
  )
}

te_data <- te_raw %>%
  
  transmute(
    nutrient = as.character(nutrient),
    
    source = norm_id(source),
    target = norm_id(target),
    
    TE_raw = as.numeric(TE_raw),
    
    p_value = as.numeric(p_value),
    
    lag_months = as.numeric(lag_months)
  )

# Must have one record per directed TE edge.

te_duplicate <- te_data %>%
  count(
    nutrient,
    source,
    target
  ) %>%
  filter(n > 1)

if (nrow(te_duplicate)) {
  
  save_csv(
    te_duplicate,
    "QA_duplicate_TE.csv"
  )
  
  stop(
    "Duplicate directed TE results."
  )
}


# ============================================================
# 6. REBUILD ALL 36 DIRECTIONS
# ============================================================

expected <- pairs %>%
  
  tidyr::crossing(
    nutrient = c("TN", "TP"),
    direction = c("Forward", "Reverse")
  ) %>%
  
  mutate(
    source = if_else(
      direction == "Forward",
      from,
      to
    ),
    
    target = if_else(
      direction == "Forward",
      to,
      from
    ),
    
    physical_reachable = (
      connection == "Physical" &
        direction == "Forward"
    ),
    
    direction_class = case_when(
      
      connection == "Physical" &
        direction == "Forward" ~ "Downstream",
      
      connection == "Physical" &
        direction == "Reverse" ~ "Reverse",
      
      connection == "Cross-tributary" ~
        "Cross-tributary"
    )
  )

directed <- expected %>%
  
  left_join(
    te_data,
    
    by = c(
      "nutrient",
      "source",
      "target"
    )
  ) %>%
  
  mutate(
    TE_positive = if_else(
      is.finite(TE_raw),
      pmax(TE_raw, 0),
      NA_real_
    ),
    
    evaluated = (
      is.finite(TE_raw) &
        is.finite(p_value)
    ),
    
    significant = (
      evaluated &
        TE_positive > 0 &
        p_value < TE_P_THRESHOLD
    ),
    
    significant = coalesce(
      significant,
      FALSE
    ),
    
    TE_status = case_when(
      !evaluated ~ "Unavailable",
      significant ~ "Significant",
      TRUE ~ "Not significant"
    ),
    
    agreement_class = case_when(
      
      !evaluated ~ "Unavailable",
      
      physical_reachable & significant ~
        "Both",
      
      physical_reachable & !significant ~
        "Physical only",
      
      !physical_reachable & significant ~
        "TE only",
      
      TRUE ~ "Neither"
    )
  )

stopifnot(
  nrow(directed) == 36
)

save_csv(
  directed,
  "02_all_36_directional_TE.csv"
)


# ============================================================
# 7. READ MONTHLY TN / TP CONCENTRATIONS
# ============================================================

check_file(CONC_FILE)

conc_raw <- readr::read_csv(
  CONC_FILE,
  
  col_types = readr::cols(
    STAID = readr::col_character(),
    PCODE = readr::col_character(),
    .default = readr::col_guess()
  ),
  
  show_col_types = FALSE
)

selected_ids <- unique(
  station_roles$STAID
)

obs <- conc_raw %>%
  
  mutate(
    STAID = norm_id(STAID),
    
    PCODE = trimws(
      as.character(PCODE)
    ),
    
    Year = as.integer(Year),
    Month = as.integer(Month),
    
    concentration = suppressWarnings(
      as.numeric(.data[[CONC_FIELD]])
    )
  ) %>%
  
  filter(
    STAID %in% selected_ids,
    
    PCODE %in% c(
      "P00600",
      "P00665"
    ),
    
    is.finite(Year),
    is.finite(Month),
    
    Month >= 1,
    Month <= 12
  ) %>%
  
  mutate(
    nutrient = case_when(
      PCODE == "P00600" ~ "TN",
      PCODE == "P00665" ~ "TP"
    ),
    
    date = as.Date(
      sprintf(
        "%04d-%02d-01",
        Year,
        Month
      )
    )
  ) %>%
  
  filter(
    date >= START_DATE,
    date <= END_DATE
  )

duplicate_obs <- obs %>%
  count(
    STAID,
    nutrient,
    date
  ) %>%
  filter(n > 1)

if (nrow(duplicate_obs)) {
  
  save_csv(
    duplicate_obs,
    "QA_duplicate_observations.csv"
  )
  
  stop(
    "Duplicate station-month records."
  )
}

obs <- obs %>%
  filter(
    is.finite(concentration),
    concentration >= 0
  ) %>%
  select(
    STAID,
    nutrient,
    date,
    concentration
  )

# Expected monthly calendar
all_dates <- seq(
  START_DATE,
  END_DATE,
  by = "month"
)

obs_coverage <- tidyr::expand_grid(
  STAID = selected_ids,
  nutrient = c("TN", "TP")
) %>%
  
  left_join(
    obs %>%
      count(
        STAID,
        nutrient,
        name = "n_valid_months"
      ),
    
    by = c("STAID", "nutrient")
  ) %>%
  
  mutate(
    n_valid_months = coalesce(
      n_valid_months,
      0L
    ),
    
    n_expected = length(all_dates),
    
    coverage = n_valid_months / n_expected
  )

save_csv(
  obs_coverage,
  "03_station_observation_coverage.csv"
)


# ============================================================
# 8. MATCH THREE STATIONS FOR PANEL B
# ============================================================

get_group_observations <- function(
    g,
    nutrient_name
) {
  
  d <- obs %>%
    
    filter(
      nutrient == nutrient_name,
      
      STAID %in% c(
        g$A,
        g$B,
        g$D
      )
    ) %>%
    
    mutate(
      role = case_when(
        STAID == g$A ~ "A",
        STAID == g$B ~ "B",
        STAID == g$D ~ "D"
      )
    ) %>%
    
    select(
      date,
      role,
      concentration
    ) %>%
    
    pivot_wider(
      names_from = role,
      values_from = concentration
    )
  
  for (nm in c("A", "B", "D")) {
    
    if (!nm %in% names(d)) {
      d[[nm]] <- NA_real_
    }
  }
  
  d %>%
    
    filter(
      is.finite(A),
      is.finite(B),
      is.finite(D)
    ) %>%
    
    mutate(
      group = g$group,
      nutrient = nutrient_name
    ) %>%
    
    select(
      group,
      nutrient,
      date,
      A,
      B,
      D
    ) %>%
    
    arrange(date)
}

group_obs <- purrr::map_dfr(
  seq_len(nrow(groups)),
  
  function(i) {
    
    bind_rows(
      lapply(
        c("TN", "TP"),
        
        function(nu) {
          
          get_group_observations(
            groups[i, ],
            nu
          )
        }
      )
    )
  }
)

save_csv(
  group_obs,
  "04_three_station_common_months.csv"
)


# ============================================================
# 9. PANEL A — STATION LOCATION AND RIVER NETWORK
# ============================================================

check_file(STATION_FILE)
check_file(RIVER_FILE)
check_file(BASIN_FILE)

station_raw <- readr::read_csv(
  STATION_FILE,
  
  col_types = readr::cols(
    SITE_NO = readr::col_character(),
    .default = readr::col_guess()
  ),
  
  show_col_types = FALSE
)

station_xy <- station_raw %>%
  transmute(
    STAID = norm_id(SITE_NO),
    lon = as.numeric(lon),
    lat = as.numeric(lat)
  ) %>%
  filter(
    is.finite(lon),
    is.finite(lat)
  ) %>%
  distinct(
    STAID,
    .keep_all = TRUE
  )

station_info <- station_roles %>%
  left_join(
    station_xy,
    by = "STAID"
  )

if (any(
  !is.finite(station_info$lon) |
  !is.finite(station_info$lat)
)) {
  
  stop(
    "Missing station coordinates."
  )
}

station_sf <- st_as_sf(
  station_info,
  coords = c("lon", "lat"),
  crs = 4326,
  remove = FALSE
)

rivers <- st_read(
  RIVER_FILE,
  quiet = TRUE
)

basin <- st_read(
  BASIN_FILE,
  quiet = TRUE
)

if (is.na(st_crs(basin))) {
  
  stop(
    "Watershed CRS is missing."
  )
}

if (is.na(st_crs(rivers))) {
  
  bb <- st_bbox(rivers)
  
  looks_geographic <- all(
    is.finite(bb)
  ) &&
    abs(bb["xmin"]) <= 180 &&
    abs(bb["xmax"]) <= 180 &&
    abs(bb["ymin"]) <= 90 &&
    abs(bb["ymax"]) <= 90
  
  if (!looks_geographic) {
    
    stop(
      "River CRS is missing and cannot ",
      "be safely inferred."
    )
  }
  
  st_crs(rivers) <- 4326
  
  warning(
    "River CRS assigned EPSG:4326. ",
    "Verify original metadata."
  )
}

rivers <- st_transform(
  rivers,
  4326
)

basin <- st_transform(
  basin,
  4326
)

map_edges <- pairs %>%
  
  left_join(
    station_xy %>%
      rename(
        from = STAID,
        x1 = lon,
        y1 = lat
      ),
    by = "from"
  ) %>%
  
  left_join(
    station_xy %>%
      rename(
        to = STAID,
        x2 = lon,
        y2 = lat
      ),
    by = "to"
  )


# ============================================================
# 9. PANEL A — THREE MAPS WITH IDENTICAL DIMENSIONS
# ============================================================

# 统一地图框的宽高比
MAP_ASPECT <- 1.5

# 地图外围留白比例
MAP_PADDING <- 0.12


# ------------------------------------------------------------
# 9.1 Calculate geographic extents
# ------------------------------------------------------------

get_map_bounds <- function(g) {
  
  p <- station_info %>%
    filter(group == g)
  
  lon_min <- min(p$lon)
  lon_max <- max(p$lon)
  
  lat_min <- min(p$lat)
  lat_max <- max(p$lat)
  
  cx <- (lon_min + lon_max) / 2
  cy <- (lat_min + lat_max) / 2
  
  # 将经纬度近似转换为实际距离比例
  cos_lat <- cos(cy * pi / 180)
  
  # 地理范围对应的宽高
  dx <- (lon_max - lon_min) * cos_lat
  dy <- lat_max - lat_min
  
  # 加入留白
  dx <- max(dx * (1 + 2 * MAP_PADDING), 0.05)
  dy <- max(dy * (1 + 2 * MAP_PADDING), 0.05)
  
  # 统一地图范围的宽高比
  if (dx / dy < MAP_ASPECT) {
    dx <- dy * MAP_ASPECT
  } else {
    dy <- dx / MAP_ASPECT
  }
  
  half_lon <- dx / (2 * cos_lat)
  half_lat <- dy / 2
  
  list(
    xmin = cx - half_lon,
    xmax = cx + half_lon,
    ymin = cy - half_lat,
    ymax = cy + half_lat
  )
}


# ------------------------------------------------------------
# 9.2 Map function
# ------------------------------------------------------------

make_group_map <- function(g) {
  
  bb <- get_map_bounds(g)
  
  pts <- station_sf %>%
    filter(group == g)
  
  seg <- map_edges %>%
    filter(group == g)
  
  ggplot() +
    
    geom_sf(
      data = basin,
      fill = "#F8F8F8",
      color = "grey83",
      linewidth = 0.15
    ) +
    
    geom_sf(
      data = rivers,
      color = "#87B7D4",
      linewidth = 0.32
    ) +
    
    geom_segment(
      data = seg,
      aes(
        x = x1,
        y = y1,
        xend = x2,
        yend = y2,
        color = connection,
        linetype = connection
      ),
      linewidth = 0.70,
      alpha = 0.90
    ) +
    
    geom_sf(
      data = pts,
      aes(fill = role),
      shape = 21,
      size = 3,
      stroke = 0.55,
      color = "#333333"
    ) +
    
    geom_text(
      data = station_info %>%
        filter(group == g),
      aes(
        x = lon,
        y = lat,
        label = STAID
      ),
      size = 2.7,
      nudge_y = 0.015,
      check_overlap = TRUE
    ) +
    
    scale_fill_manual(
      values = c(
        A = COL_A,
        B = COL_B,
        D = COL_D
      )
    ) +
    
    scale_color_manual(
      values = c(
        Physical = COL_BOTH,
        "Cross-tributary" = COL_TE_ONLY
      )
    ) +
    
    scale_linetype_manual(
      values = c(
        Physical = "solid",
        "Cross-tributary" = "dashed"
      )
    ) +
    
    coord_sf(
      xlim = c(bb$xmin, bb$xmax),
      ylim = c(bb$ymin, bb$ymax),
      datum = NA,
      expand = FALSE,
      default_crs = sf::st_crs(4326)
    ) +
    
    labs(
      title = g,
      x = NULL,
      y = NULL
    ) +
    
    theme_fig(9) +
    
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none",
      
      # 统一地图框比例
      aspect.ratio = 1 / MAP_ASPECT,
      
      plot.margin = margin(
        4, 4, 4, 4
      ),
      
      plot.title = element_text(
        face = "bold",
        size = 11,
        hjust = 0
      )
    )
}


# ------------------------------------------------------------
# 9.3 Assemble three equal-size maps
# ------------------------------------------------------------

p_map1 <- make_group_map("G1")
p_map2 <- make_group_map("G2")
p_map3 <- make_group_map("G3")

p_a <- patchwork::wrap_plots(
  p_map1,
  p_map2,
  p_map3,
  ncol = 3,
  widths = c(1, 1, 1)
) +
  patchwork::plot_annotation(
    title = "(a) Representative river-network configurations"
  )


# ------------------------------------------------------------
# 9.4 Export map panel
# ------------------------------------------------------------

ggsave(
  file.path(
    OUT,
    "Fig3a_equal_size_maps.pdf"
  ),
  plot = p_a,
  width = 15,
  height = 4.5,
  bg = "white"
)

ggsave(
  file.path(
    OUT,
    "Fig3a_equal_size_maps.png"
  ),
  plot = p_a,
  width = 15,
  height = 4.5,
  dpi = 400,
  bg = "white"
)

print(p_a)

# ============================================================
# 10. PANEL B — OBSERVED TN / TP DYNAMICS
# ============================================================

time_long <- group_obs %>%
  
  pivot_longer(
    cols = c(A, B, D),
    
    names_to = "role",
    values_to = "concentration"
  ) %>%
  
  group_by(
    group,
    nutrient,
    role
  ) %>%
  
  mutate(
    z = zscore(concentration)
  ) %>%
  
  ungroup()

# Standardization based on the entire study period.

if (PLOT_MONTHS > 0) {
  
  display_dates <- tail(
    all_dates,
    PLOT_MONTHS
  )
  
  time_long <- time_long %>%
    filter(
      date %in% display_dates
    )
}

p_b <- ggplot(
  time_long,
  
  aes(
    x = date,
    y = z,
    color = role
  )
) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.3,
    color = "grey70"
  ) +
  
  geom_line(
    linewidth = 0.55,
    alpha = 0.9
  ) +
  
  facet_grid(
    group ~ nutrient,
    scales = "free_y"
  ) +
  
  scale_color_manual(
    values = c(
      A = COL_A,
      B = COL_B,
      D = COL_D
    ),
    
    labels = c(
      A = "Tributary A",
      B = "Tributary B",
      D = "Downstream"
    )
  ) +
  
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y"
  ) +
  
  labs(
    title = "(b) Observed TN and TP concentration dynamics",
    
    subtitle = paste0(
      "Three-station common months; ",
      "standardized over the full period"
    ),
    
    x = "Year",
    y = "Standardized concentration",
    color = NULL
  ) +
  
  theme_fig(10) +
  
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    )
  )


# ============================================================
# 11. PANEL C — BIDIRECTIONAL TE
# ============================================================

c_data <- directed %>%
  
  mutate(
    pair = factor(
      pair,
      levels = rev(PAIR_LEVELS)
    ),
    
    direction = factor(
      direction,
      levels = c(
        "Forward",
        "Reverse"
      )
    ),
    
    TE_label = case_when(
      
      !evaluated ~ "NA",
      
      significant ~ paste0(
        sprintf(
          "%.3f",
          TE_positive
        ),
        "*"
      ),
      
      TRUE ~ sprintf(
        "%.3f",
        TE_positive
      )
    )
  )

max_TE <- max(
  c_data$TE_positive,
  na.rm = TRUE
)

if (!is.finite(max_TE) || max_TE <= 0) {
  max_TE <- 1
}

p_c <- ggplot(
  c_data,
  
  aes(
    x = direction,
    y = pair,
    fill = TE_positive
  )
) +
  
  geom_tile(
    width = 0.94,
    height = 0.94,
    color = "white",
    linewidth = 0.65
  ) +
  
  geom_text(
    aes(label = TE_label),
    size = 2.8,
    color = "#1F3936"
  ) +
  
  scale_fill_gradientn(
    colours = TE_PALETTE,
    limits = c(0, max_TE),
    na.value = "#D9D9D9",
    name = "TE strength"
  ) +
  
  facet_wrap(
    ~nutrient,
    nrow = 1
  ) +
  
  labs(
    title = "(c) Directional transfer entropy",
    
    subtitle = paste0(
      "* nominal p < ",
      TE_P_THRESHOLD,
      "; Forward: U→D or A→B"
    ),
    
    x = NULL,
    y = "Station pair"
  ) +
  
  theme_fig(10) +
  
  theme(
    panel.grid = element_blank()
  )


# ============================================================
# 12. PANEL D — PHYSICAL–TE AGREEMENT
# ============================================================

d_data <- directed %>%
  
  filter(evaluated) %>%
  
  count(
    nutrient,
    agreement_class,
    name = "n"
  ) %>%
  
  complete(
    nutrient = c("TN", "TP"),
    
    agreement_class = c(
      "Both",
      "Physical only",
      "TE only",
      "Neither"
    ),
    
    fill = list(n = 0L)
  ) %>%
  
  mutate(
    agreement_class = factor(
      agreement_class,
      
      levels = c(
        "Both",
        "Physical only",
        "TE only",
        "Neither"
      )
    )
  )

save_csv(
  d_data,
  "05_connectivity_TE_agreement.csv"
)

p_d <- ggplot(
  d_data,
  
  aes(
    x = agreement_class,
    y = n,
    fill = agreement_class
  )
) +
  
  geom_col(
    width = 0.67,
    color = "white"
  ) +
  
  geom_text(
    aes(label = n),
    vjust = -0.45,
    size = 3.2
  ) +
  
  facet_wrap(
    ~nutrient,
    nrow = 1
  ) +
  
  scale_fill_manual(
    values = c(
      "Both" = COL_BOTH,
      "Physical only" = COL_PHYS_ONLY,
      "TE only" = COL_TE_ONLY,
      "Neither" = COL_NEITHER
    ),
    
    guide = "none"
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(0, 0.13)
    )
  ) +
  
  labs(
    title = "(d) Connectivity–TE agreement",
    
    subtitle = paste0(
      "Physical reachability versus ",
      "significant directional TE"
    ),
    
    x = NULL,
    y = "Directed relationships"
  ) +
  
  theme_fig(10) +
  
  theme(
    axis.text.x = element_text(
      angle = 25,
      hjust = 1,
      size = 8
    ),
    
    panel.grid.major.x = element_blank()
  )


# ============================================================
# 13. PREPARE SEASONALLY STANDARDIZED ANOMALIES
# ============================================================

# Important:
#
# Use the original observed data, not TE weights.
#
# Each station and nutrient is standardized
# separately within each calendar month.
#
# E.g. all Januaries are standardized together,
# all Februaries together, etc.
#
# Missing monthly observations are not filled.

seasonal_obs <- obs %>%
  
  mutate(
    calendar_month = as.integer(
      format(date, "%m")
    )
  ) %>%
  
  group_by(
    STAID,
    nutrient,
    calendar_month
  ) %>%
  
  mutate(
    n_climatology = sum(
      is.finite(concentration)
    ),
    
    monthly_mean = mean(
      concentration,
      na.rm = TRUE
    ),
    
    monthly_sd = sd(
      concentration,
      na.rm = TRUE
    ),
    
    seasonal_z = if_else(
      
      n_climatology >= MIN_CLIMATOLOGY_YEARS &
        is.finite(monthly_sd) &
        monthly_sd > 0,
      
      (
        concentration - monthly_mean
      ) / monthly_sd,
      
      NA_real_
    ),
    
    high_event = if_else(
      is.finite(seasonal_z),
      seasonal_z > Z_THRESHOLD,
      NA
    )
  ) %>%
  
  ungroup()

save_csv(
  seasonal_obs,
  "06_seasonal_concentration_anomalies.csv"
)


# ============================================================
# 14. DEFINE SIX PHYSICAL UPSTREAM–DOWNSTREAM PAIRS
# ============================================================

physical_pairs <- pairs %>%
  
  filter(
    connection == "Physical"
  ) %>%
  
  transmute(
    group,
    pair,
    upstream = from,
    downstream = to
  )

stopifnot(
  nrow(physical_pairs) == 6
)


# ============================================================
# 15. LAGGED EVENT CORRESPONDENCE
# ============================================================

# For a fixed station pair and lag k:
#
# p_cond =
# P[Downstream high at t+k | Upstream high at t]
#
# p_base =
# P[Downstream high at t+k]
#
# delta_p = p_cond - p_base
#
# The same matched calendar months are used
# in both probabilities.
#
# A positive delta_p indicates enrichment of
# downstream high anomalies conditional on
# upstream high anomalies.
#
# It does not establish causal transport.


calculate_lag <- function(
    group_id,
    pair_id,
    upstream_id,
    downstream_id,
    nutrient_name,
    k
) {
  
  upstream_obs <- seasonal_obs %>%
    
    filter(
      STAID == upstream_id,
      nutrient == nutrient_name
    ) %>%
    
    transmute(
      date_up = date,
      upstream_event = high_event
    ) %>%
    
    mutate(
      # Shift by actual calendar-month position.
      # This avoids treating missing records as
      # consecutive months.
      
      month_index = match(
        date_up,
        all_dates
      ),
      
      target_index = month_index + k,
      
      date_down = as.Date(
        all_dates[
          ifelse(
            target_index >= 1 &
              target_index <= length(all_dates),
            target_index,
            NA_integer_
          )
        ]
      )
    ) %>%
    
    filter(
      !is.na(date_down)
    ) %>%
    
    select(
      date_down,
      upstream_event
    )
  
  downstream_obs <- seasonal_obs %>%
    
    filter(
      STAID == downstream_id,
      nutrient == nutrient_name
    ) %>%
    
    transmute(
      date_down = date,
      downstream_event = high_event
    )
  
  matched <- inner_join(
    upstream_obs,
    downstream_obs,
    by = "date_down"
  ) %>%
    
    filter(
      !is.na(upstream_event),
      !is.na(downstream_event)
    )
  
  n_valid <- nrow(matched)
  
  n_up <- sum(
    matched$upstream_event
  )
  
  n_down <- sum(
    matched$downstream_event
  )
  
  n_joint <- sum(
    matched$upstream_event &
      matched$downstream_event
  )
  
  if (n_valid > 0) {
    
    p_background <- n_down / n_valid
    
  } else {
    
    p_background <- NA_real_
  }
  
  if (n_up > 0) {
    
    p_conditional <- n_joint / n_up
    
  } else {
    
    p_conditional <- NA_real_
  }
  
  delta_raw <- (
    p_conditional - p_background
  )
  
  # Minimum-event criterion applies to
  # the plotted estimate, not raw counts.
  
  delta_plot <- if (
    n_up >= MIN_UPSTREAM_EVENTS &&
    is.finite(delta_raw)
  ) {
    
    delta_raw
    
  } else {
    
    NA_real_
  }
  
  tibble(
    group = group_id,
    pair = pair_id,
    
    upstream = upstream_id,
    downstream = downstream_id,
    
    nutrient = nutrient_name,
    
    lag_months = k,
    
    n_valid = n_valid,
    
    n_upstream_events = n_up,
    n_downstream_events = n_down,
    n_joint_events = n_joint,
    
    p_conditional = p_conditional,
    p_background = p_background,
    
    delta_probability_raw = delta_raw,
    
    delta_probability = delta_plot,
    
    sufficient_events = (
      n_up >= MIN_UPSTREAM_EVENTS
    )
  )
}


# ============================================================
# 16. CALCULATE ALL 48 RESULTS
# ============================================================

event_design <- tidyr::crossing(
  physical_pairs,
  
  nutrient = c("TN", "TP"),
  
  lag_months = LAGS
)

event_results <- purrr::pmap_dfr(
  event_design,
  
  function(
    group,
    pair,
    upstream,
    downstream,
    nutrient,
    lag_months
  ) {
    
    calculate_lag(
      group_id = group,
      pair_id = pair,
      
      upstream_id = upstream,
      downstream_id = downstream,
      
      nutrient_name = nutrient,
      k = lag_months
    )
  }
)

stopifnot(
  nrow(event_results) == 48
)

save_csv(
  event_results,
  "07_lagged_anomaly_correspondence.csv"
)


# ============================================================
# 17. EVENT QA
# ============================================================

event_QA <- event_results %>%
  
  group_by(
    nutrient,
    lag_months
  ) %>%
  
  summarise(
    n_pairs = n(),
    
    n_sufficient = sum(
      sufficient_events
    ),
    
    median_upstream_events = median(
      n_upstream_events
    ),
    
    mean_delta = if (
      any(is.finite(delta_probability))
    ) {
      
      mean(
        delta_probability,
        na.rm = TRUE
      )
      
    } else {
      
      NA_real_
    },
    
    .groups = "drop"
  )

save_csv(
  event_QA,
  "08_lagged_anomaly_QA.csv"
)


# ============================================================
# 18. PANEL E — LAGGED ANOMALY HEATMAP
# ============================================================

e_data <- event_results %>%
  
  mutate(
    pair = factor(
      pair,
      levels = rev(physical_pairs$pair)
    ),
    
    lag_months = factor(
      lag_months,
      levels = LAGS
    ),
    
    label = case_when(
      
      !sufficient_events ~ "NA",
      
      !is.finite(delta_probability) ~ "NA",
      
      TRUE ~ sprintf(
        "%+.0f",
        100 * delta_probability
      )
    )
  )

max_abs <- max(
  abs(e_data$delta_probability),
  na.rm = TRUE
)

if (
  !is.finite(max_abs) ||
  max_abs <= 0
) {
  max_abs <- 0.10
}

# Round to a slightly wider color limit.
max_abs <- max(
  0.10,
  ceiling(max_abs * 20) / 20
)

p_e <- ggplot(
  e_data,
  
  aes(
    x = lag_months,
    y = pair,
    fill = delta_probability
  )
) +
  
  geom_tile(
    width = 0.94,
    height = 0.94,
    color = "white",
    linewidth = 0.65
  ) +
  
  geom_text(
    aes(label = label),
    size = 2.9
  ) +
  
  facet_wrap(
    ~nutrient,
    nrow = 1
  ) +
  
  scale_fill_gradient2(
    low = COL_NEG,
    mid = COL_ZERO,
    high = COL_POS,
    
    midpoint = 0,
    
    limits = c(
      -max_abs,
      max_abs
    ),
    
    oob = scales::squish,
    
    na.value = "#D9D9D9",
    
    labels = function(x) {
      sprintf(
        "%+.0f pp",
        100 * x
      )
    },
    
    name = "Probability\ndifference"
  ) +
  
  labs(
    title = "(e) Lagged correspondence of nutrient anomalies",
    
    subtitle = paste0(
      "Seasonal z > ",
      Z_THRESHOLD,
      "; ΔP in percentage points; ",
      "NA: insufficient events"
    ),
    
    x = "Lag after upstream anomaly (months)",
    
    y = "Upstream–downstream pair"
  ) +
  
  theme_fig(10) +
  
  theme(
    panel.grid = element_blank(),
    
    axis.text.x = element_text(
      size = 9
    )
  )


# ============================================================
# 19. COMBINE FINAL FIGURE
# ============================================================

# Recommended:
#
# Row 1: Spatial configuration
# Row 2: Observed time series
# Row 3: TE strength + connectivity agreement
# Row 4: Observed event correspondence
#
# Panel E is wider because it contains
# 6 pairs × 4 lags × 2 nutrients.

middle_bottom <- patchwork::wrap_plots(
  p_c,
  p_d,
  
  ncol = 2,
  
  widths = c(1.20, 1.00)
)

fig3 <- patchwork::wrap_plots(
  p_a,
  p_b,
  middle_bottom,
  p_e,
  
  ncol = 1,
  
  heights = c(
    0.95,
    1.40,
    1.15,
    1.05
  )
) +
  
  patchwork::plot_annotation(
    
    title = paste0(
      "Physical Connectivity, Directional Information ",
      "and Observed Nutrient Anomaly Correspondence"
    ),
    
    subtitle = paste0(
      "Three river configurations | ",
      "Nine station pairs | ",
      "TN and TP | ",
      "October 2011–September 2023"
    ),
    
    caption = paste0(
      "Panel (a) uses schematic station-to-station links. ",
      "TE significance is nominal (p < 0.10), ",
      "without multiple-testing correction. ",
      "Panel (e) uses observed concentration anomalies ",
      "independently of TE and reports conditional ",
      "probability differences, not causal effects."
    ),
    
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 15
      ),
      
      plot.subtitle = element_text(
        size = 10
      ),
      
      plot.caption = element_text(
        size = 8,
        color = "grey40"
      )
    )
  )


# ============================================================
# 20. EXPORT
# ============================================================

ggsave(
  file.path(
    OUT,
    "Fig3_observational_validation.pdf"
  ),
  
  plot = fig3,
  
  width = 17,
  height = 18,
  
  bg = "white"
)

ggsave(
  file.path(
    OUT,
    "Fig3_observational_validation.png"
  ),
  
  plot = fig3,
  
  width = 17,
  height = 18,
  
  dpi = FIG_DPI,
  
  bg = "white"
)

ggsave(
  file.path(
    OUT,
    "Fig3_observational_validation.tif"
  ),
  
  plot = fig3,
  
  width = 17,
  height = 18,
  
  dpi = FIG_DPI,
  compression = "lzw",
  
  bg = "white"
)


# ============================================================
# 21. INDIVIDUAL PANELS
# ============================================================

panel_plots <- list(
  Fig3a_river_map = p_a,
  Fig3b_observed_timeseries = p_b,
  Fig3c_directional_TE = p_c,
  Fig3d_connectivity_TE = p_d,
  Fig3e_anomaly_correspondence = p_e
)

panel_sizes <- list(
  c(15, 5),
  c(14, 8),
  c(11, 7),
  c(10, 7),
  c(13, 7)
)

for (i in seq_along(panel_plots)) {
  
  ggsave(
    file.path(
      OUT,
      paste0(
        names(panel_plots)[i],
        ".pdf"
      )
    ),
    
    plot = panel_plots[[i]],
    
    width = panel_sizes[[i]][1],
    height = panel_sizes[[i]][2],
    
    bg = "white"
  )
}


# ============================================================
# 22. REPORT RESULTS
# ============================================================

cat("\n")
cat("=============================================\n")
cat("FIGURE 3 COMPLETED\n")
cat("=============================================\n")

cat("\nOBSERVATION COVERAGE:\n")
print(
  obs_coverage,
  n = Inf
)

cat("\nTE SIGNIFICANCE:\n")

print(
  directed %>%
    group_by(
      nutrient,
      direction_class
    ) %>%
    summarise(
      n_expected = n(),
      n_evaluated = sum(evaluated),
      n_significant = sum(significant),
      .groups = "drop"
    ),
  n = Inf
)

cat("\nPHYSICAL–TE AGREEMENT:\n")
print(
  d_data,
  n = Inf
)

cat("\nANOMALY EVENT QA:\n")
print(
  event_QA,
  n = Inf
)

cat("\nALL 48 ANOMALY RESULTS:\n")

print(
  event_results %>%
    select(
      pair,
      nutrient,
      lag_months,
      n_valid,
      n_upstream_events,
      n_joint_events,
      p_conditional,
      p_background,
      delta_probability
    ),
  n = Inf
)

cat("\nOutput folder:\n")
cat(OUT, "\n")

print(fig3)
