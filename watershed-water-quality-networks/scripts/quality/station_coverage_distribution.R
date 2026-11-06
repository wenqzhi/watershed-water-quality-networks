source("config.R")
ROOT <- get_data_root()
# ============================================================
# Fig. 1 | 122 NTN stations and temporal coverage
# Panels:
#   A. Station map over watershed and river network
#   B. Station-level temporal coverage (gantt-style timeline)
#   C. Number of stations available under different common windows
#
# Input files:
#   1) QA_122_stations_summary.csv
#   2) NTN-station.txt   (FID,SITE_NO,lon,lat)
#   3) ntnwatersheds.shp
#   4) 水系.shp
# ============================================================

# -----------------------------
# 0. Packages
# -----------------------------
need_pkgs <- c(
  "sf", "dplyr", "readr", "stringr", "lubridate",
  "ggplot2", "patchwork", "tidyr", "forcats", "scales"
)

for (p in need_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

library(sf)
library(dplyr)
library(readr)
library(stringr)
library(lubridate)
library(ggplot2)
library(patchwork)
library(tidyr)
library(forcats)
library(scales)

# -----------------------------
# 1. File paths (modify as needed)
# -----------------------------
qa_file <- file.path(ROOT, "NTN_completeness_QA/QA_122_stations_summary.csv")
station_file <- file.path(ROOT, "NTN-station.txt")
watershed_shp <- file.path(ROOT, "ntnwatersheds/ntnwatersheds.shp")
river_shp <- file.path(ROOT, "水系/水系.shp")
out_dir <- file.path(ROOT, "figures", "Fig1_NTN_coverage")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -----------------------------
# 2. Read data
# -----------------------------
qa <- read_csv(qa_file, show_col_types = FALSE) %>%
  mutate(
    STAID = str_pad(as.character(STAID), width = 8, pad = "0"),
    start_date = ym(start_month),
    end_date   = ym(end_month)
  )

stations <- read_csv(station_file, show_col_types = FALSE) %>%
  transmute(
    STAID = str_pad(as.character(SITE_NO), width = 8, pad = "0"),
    lon = as.numeric(lon),
    lat = as.numeric(lat)
  )

watersheds <- st_read(watershed_shp, quiet = TRUE)
rivers <- st_read(river_shp, quiet = TRUE)

# Reproject to lon/lat if needed
if (is.na(st_crs(watersheds))) stop("流域 SHP 缺少 CRS，请先确认真实坐标系")
if (is.na(st_crs(rivers))) {
  st_crs(rivers) <- 4326
}

watersheds <- st_transform(watersheds, 4326)
rivers <- st_transform(rivers, 4326)

# Check spatial extents before plotting (incorrect CRS often causes invisible rivers)
cat("Watershed extent:", paste(round(st_bbox(watersheds), 3), collapse = ", "), "\n")
cat("River extent:", paste(round(st_bbox(rivers), 3), collapse = ", "), "\n")
if (!any(lengths(st_intersects(st_as_sfc(st_bbox(watersheds)), st_as_sfc(st_bbox(rivers)))) > 0)) {
  warning("河网与流域边界范围未重叠：请检查河网 CRS、投影和坐标单位")
}
# Limit rivers to the watershed bounding rectangle for clearer map plotting
rivers_plot <- suppressWarnings(st_crop(rivers, st_bbox(watersheds)))

stations_sf <- stations %>%
  left_join(qa, by = "STAID") %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# -----------------------------
# 3. Coverage groups for visualization
# -----------------------------
# Start-time cohorts often used for this dataset
stations_sf <- stations_sf %>%
  mutate(
    start_group = case_when(
      start_date <= ymd("1984-10-01") ~ "1984.10 start",
      start_date <= ymd("1997-10-01") ~ "1997.10 start",
      start_date <= ymd("2004-10-01") ~ "2004.10 start",
      start_date <= ymd("2011-10-01") ~ "2011.10 start",
      TRUE                               ~ "2014.10 start"
    ),
    start_group = factor(
      start_group,
      levels = c("1984.10 start", "1997.10 start", "2004.10 start", "2011.10 start", "2014.10 start")
    )
  )

# -----------------------------
# 4. Common windows and station counts
# -----------------------------
common_windows <- tibble(
  period = c("1984.10–2023.09", "2004.10–2023.09", "2011.10–2023.09", "2014.10–2023.09"),
  w_start = ymd(c("1984-10-01", "2004-10-01", "2011-10-01", "2014-10-01")),
  w_end   = ymd(c("2023-09-01", "2023-09-01", "2023-09-01", "2023-09-01"))
)

window_counts <- common_windows %>%
  rowwise() %>%
  mutate(
    n_station = sum(qa$start_date <= w_start & qa$end_date >= w_end, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(period = factor(period, levels = period))

# -----------------------------
# 5. Panel A - map
# -----------------------------
# Use a clean extent based on watersheds + stations
bbox_obj <- st_bbox(watersheds)

p_map <- ggplot() +
  geom_sf(data = watersheds, fill = "grey97", color = "grey95", linewidth = 0.3) +
  geom_sf(data = rivers_plot, color = "#398AC6", linewidth = 0.55, alpha = 0.95, show.legend = FALSE) +
  geom_sf(data = watersheds, fill = NA, color = "grey85", linewidth = 0.38) +
  geom_sf(data = stations_sf, color = "white", size = 3.15, show.legend = FALSE) +
  geom_sf(data = stations_sf, aes(color = start_group), size = 2.15, alpha = 0.98) +
  scale_color_manual(
    values = c(
      "1984.10 start" = "#B2182B",
      "1997.10 start" = "#EF8A62",
      "2004.10 start" = "#FDCC8A",
      "2011.10 start" = "#67A583",
      "2014.10 start" = "#7B95C6"
    ),
    name = "Station start cohort"
  ) +
  coord_sf(
    xlim = c(bbox_obj["xmin"] - 0.3, bbox_obj["xmax"] + 0.3),
    ylim = c(bbox_obj["ymin"] - 0.2, bbox_obj["ymax"] + 0.2),
    expand = FALSE
  ) +
  labs(
    title = "NTN stations and river network (122 sites)",
    subtitle = "Blue lines: river network; station colors: monitoring start cohort"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_line(color = "grey95", linewidth = 0.2),
    legend.position = c(0.16, 0.18),
    legend.background = element_rect(fill = alpha("white", 0.85), color = "grey80"),
    plot.title = element_text(face = "bold"),
    axis.title = element_blank()
  )

# -----------------------------
# 6. Panel B - station-level coverage timeline
# -----------------------------
# Sort by start date (and latitude as tie-breaker) for clearer display
coverage_df <- stations_sf %>%
  st_drop_geometry() %>%
  arrange(start_date, desc(lat)) %>%
  mutate(
    station_order = row_number(),
    station_label = factor(STAID, levels = STAID)
  )

# Optional reference windows shown as light background bands
window_bands <- tibble(
  xmin = ymd(c("1984-10-01", "2004-10-01", "2011-10-01", "2014-10-01")),
  xmax = ymd(c("2023-09-01", "2023-09-01", "2023-09-01", "2023-09-01")),
  ymin = c(0.5, 0.5, 0.5, 0.5),
  ymax = c(122.5, 122.5, 122.5, 122.5),
  band = c("1984–2023", "2004–2023", "2011–2023", "2014–2023")
)

p_timeline <- ggplot(coverage_df) +
  geom_segment(
    aes(x = start_date, xend = end_date, y = station_order, yend = station_order, color = start_group),
    linewidth = 1.0, lineend = "round"
  ) +
  scale_color_manual(
    values = c(
      "1984.10 start" = "#B2182B",
      "1997.10 start" = "#EF8A62",
      "2004.10 start" = "#FDCC8A",
      "2011.10 start" = "#67A583",
      "2014.10 start" = "#7B95C6"
    ),
    guide = "none"
  ) +
  scale_x_date(
    breaks = seq(ymd("1985-01-01"), ymd("2025-01-01"), by = "5 years"),
    labels = label_date("%Y")
  ) +
  scale_y_continuous(
    breaks = c(1, 20, 40, 60, 80, 100, 122),
    expand = c(0.01, 0.01)
  ) +
  labs(
    title = "Temporal coverage of 122 stations",
    subtitle = "Each line represents one station (ordered by start month)",
    x = NULL,
    y = "Station order"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

# -----------------------------
# 7. Panel C - count stations under common windows
# -----------------------------
p_bar <- ggplot(window_counts, aes(x = period, y = n_station)) +
  geom_col(fill = "grey45", width = 0.72) +
  geom_text(aes(label = n_station), vjust = -0.4, size = 4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Number of stations available in common windows",
    x = NULL,
    y = "No. of stations"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 15, hjust = 1),
    plot.title = element_text(face = "bold")
  )

# -----------------------------
# 8. Combine panels
# -----------------------------
fig1 <- p_map + (p_timeline / p_bar) +
  plot_layout(widths = c(1.05, 1.25)) +
  plot_annotation(tag_levels = "A")

# -----------------------------
# 9. Export
# -----------------------------
out_png <- file.path(out_dir, "Fig1_NTN_stations_and_coverage.png")
out_pdf <- file.path(out_dir, "Fig1_NTN_stations_and_coverage.pdf")
out_tif <- file.path(out_dir, "Fig1_NTN_stations_and_coverage.tif")

ggsave(out_png, fig1, width = 15, height = 9, dpi = 400, bg = "white")
ggsave(out_pdf, fig1, width = 15, height = 9, bg = "white")
ggsave(out_tif, fig1, width = 15, height = 9, dpi = 400, compression = "lzw", bg = "white")

# Also export key summary tables used in the figure
write_csv(window_counts, file.path(out_dir, "Fig1_common_window_station_counts.csv"))
write_csv(coverage_df, file.path(out_dir, "Fig1_station_coverage_table.csv"))

print(fig1)
cat("\nFigure saved to:\n", out_png, "\n", out_pdf, "\n", out_tif, "\n")
