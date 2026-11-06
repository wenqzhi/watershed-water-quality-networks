source("config.R")
# ============================================================
# Figure S4 | PCS versus NTN nutrient loads
#
# 数据：06_unified_station_metrics.csv
#
# 输出两套结果：
# 1. Original：TN 84站、TP 73站
# 2. Sensitivity：删除TN年负荷最低的1个站点
#
# (a) TN: PCS vs. annual load
# (b) TP: PCS vs. annual load
# (c) High-priority station agreement
# (d) Top-ranked station Jaccard overlap
#
# 修复：
# - high_PCS / high_load 字段冲突
# - 删除站点后重新计算排名和重合度
# - 同时保存删除前后结果
# ============================================================


# ============================================================
# 1. Packages
# ============================================================

packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "ggplot2",
  "patchwork",
  "scales"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)


# ============================================================
# 2. Paths
# ============================================================

ROOT <- get_data_root()

INPUT_FILE <- file.path(
  ROOT,
  "Fig4_final_independent",
  "06_unified_station_metrics.csv"
)

OUT_DIR <- file.path(
  ROOT,
  "Figure_S4_TN_min_sensitivity"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(INPUT_FILE)) {
  stop("Input file not found: ", INPUT_FILE)
}


# ============================================================
# 3. Parameters
# ============================================================

REMOVE_TN_MIN <- TRUE

TOP_FRACTIONS <- c(0.10, 0.20, 0.30)

N_RANDOM <- 3000L

RANDOM_SEED <- 20261101L

COLORS <- c(
  TN = "#307FA6",
  TP = "#DE8750"
)

PRIORITY_COLORS <- c(
  "Both high" = "#448477",
  "PCS only" = "#307FA6",
  "Load only" = "#DE8750",
  "Neither" = "#B5BBC2"
)


# ============================================================
# 4. Read original data
# ============================================================

d0 <- read_csv(
  INPUT_FILE,
  col_types = cols(
    STAID = col_character(),
    .default = col_guess()
  ),
  show_col_types = FALSE
)

required_columns <- c(
  "STAID",
  "nutrient",
  "PCS",
  "mean_annual_load",
  "common_valid"
)

missing_columns <- setdiff(
  required_columns,
  names(d0)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# 只保留本次分析需要的字段
# 避免原始 high_PCS / high_load 影响重新计算

d_original <- d0 %>%
  transmute(
    STAID = as.character(STAID),
    nutrient = toupper(as.character(nutrient)),
    PCS = as.numeric(PCS),
    mean_annual_load = as.numeric(mean_annual_load),
    common_valid = common_valid
  ) %>%
  filter(
    nutrient %in% c("TN", "TP"),
    common_valid %in% TRUE,
    !is.na(STAID),
    is.finite(PCS),
    is.finite(mean_annual_load),
    mean_annual_load > 0
  ) %>%
  arrange(nutrient, STAID)

if (anyDuplicated(
  paste(d_original$nutrient, d_original$STAID)
)) {
  stop("Duplicate nutrient-station combinations detected.")
}

cat("\n========== Original station counts ==========\n")

print(
  d_original %>%
    count(nutrient)
)


# ============================================================
# 5. Remove lowest TN annual load
# ============================================================

tn_min <- d_original %>%
  filter(nutrient == "TN") %>%
  arrange(mean_annual_load, STAID) %>%
  slice(1)

cat("\n========== Removed TN station ==========\n")

print(tn_min, width = Inf)

write_csv(
  tn_min,
  file.path(
    OUT_DIR,
    "00_removed_TN_lowest_load.csv"
  )
)

if (REMOVE_TN_MIN) {
  
  d_sensitivity <- d_original %>%
    filter(
      !(
        nutrient == "TN" &
          STAID == tn_min$STAID[1]
      )
    )
  
} else {
  
  d_sensitivity <- d_original
  
}

cat("\n========== Sensitivity station counts ==========\n")

print(
  d_sensitivity %>%
    count(nutrient)
)

write_csv(
  d_original,
  file.path(
    OUT_DIR,
    "01_original_stations.csv"
  )
)

write_csv(
  d_sensitivity,
  file.path(
    OUT_DIR,
    "02_sensitivity_stations.csv"
  )
)


# ============================================================
# 6. Spearman correlations
# ============================================================

calculate_correlations <- function(df) {
  
  df %>%
    group_by(nutrient) %>%
    summarise(
      n = n(),
      spearman_rho = cor(
        PCS,
        mean_annual_load,
        method = "spearman"
      ),
      .groups = "drop"
    )
  
}

corr_original <- calculate_correlations(
  d_original
)

corr_sensitivity <- calculate_correlations(
  d_sensitivity
)

corr_comparison <- bind_rows(
  corr_original %>%
    mutate(scenario = "Original"),
  corr_sensitivity %>%
    mutate(scenario = "TN minimum removed")
) %>%
  select(
    scenario,
    nutrient,
    n,
    spearman_rho
  )

cat("\n========== Spearman comparison ==========\n")

print(corr_comparison, n = Inf)

write_csv(
  corr_comparison,
  file.path(
    OUT_DIR,
    "03_before_after_spearman.csv"
  )
)


# ============================================================
# 7. Priority classification
# ============================================================

classify_priority <- function(df, fraction = 0.20) {
  
  df <- df %>%
    select(
      -any_of(
        c(
          "rank_pcs",
          "rank_load",
          "high_PCS",
          "high_load",
          "priority_class"
        )
      )
    ) %>%
    arrange(STAID)
  
  n_all <- nrow(df)
  
  n_top <- max(
    1L,
    ceiling(n_all * fraction)
  )
  
  pcs_order <- order(
    -df$PCS,
    df$STAID
  )
  
  load_order <- order(
    -df$mean_annual_load,
    df$STAID
  )
  
  pcs_rank <- integer(n_all)
  load_rank <- integer(n_all)
  
  pcs_rank[pcs_order] <- seq_len(n_all)
  load_rank[load_order] <- seq_len(n_all)
  
  df %>%
    mutate(
      rank_pcs = pcs_rank,
      rank_load = load_rank,
      high_PCS = rank_pcs <= n_top,
      high_load = rank_load <= n_top,
      priority_class = case_when(
        high_PCS & high_load ~ "Both high",
        high_PCS & !high_load ~ "PCS only",
        !high_PCS & high_load ~ "Load only",
        TRUE ~ "Neither"
      )
    )
  
}


calculate_priority <- function(df) {
  
  df %>%
    group_by(nutrient) %>%
    group_modify(
      ~ classify_priority(
        .x,
        fraction = 0.20
      )
    ) %>%
    ungroup()
  
}


priority_original <- calculate_priority(
  d_original
)

priority_sensitivity <- calculate_priority(
  d_sensitivity
)

write_csv(
  priority_original,
  file.path(
    OUT_DIR,
    "04_priority_original.csv"
  )
)

write_csv(
  priority_sensitivity,
  file.path(
    OUT_DIR,
    "05_priority_sensitivity.csv"
  )
)


# ============================================================
# 8. Priority QA
# ============================================================

priority_QA <- function(df, scenario) {
  
  df %>%
    group_by(nutrient) %>%
    summarise(
      total_stations = n(),
      expected_top20 = ceiling(n() * 0.20),
      PCS_top20 = sum(high_PCS),
      Load_top20 = sum(high_load),
      Both_high = sum(
        high_PCS & high_load
      ),
      .groups = "drop"
    ) %>%
    mutate(
      scenario = scenario,
      .before = 1
    )
  
}

qa_priority <- bind_rows(
  priority_QA(
    priority_original,
    "Original"
  ),
  priority_QA(
    priority_sensitivity,
    "TN minimum removed"
  )
)

print(qa_priority, width = Inf)

write_csv(
  qa_priority,
  file.path(
    OUT_DIR,
    "06_priority_QA.csv"
  )
)

if (any(
  qa_priority$PCS_top20 !=
  qa_priority$expected_top20
) ||
any(
  qa_priority$Load_top20 !=
  qa_priority$expected_top20
)) {
  
  stop("Priority classification QA failed.")
  
}


# ============================================================
# 9. Priority composition
# ============================================================

priority_levels <- c(
  "Neither",
  "Load only",
  "PCS only",
  "Both high"
)

summarise_priority <- function(df) {
  
  df %>%
    count(
      nutrient,
      priority_class,
      name = "n"
    ) %>%
    complete(
      nutrient = c("TN", "TP"),
      priority_class = priority_levels,
      fill = list(n = 0)
    ) %>%
    group_by(nutrient) %>%
    mutate(
      total = sum(n),
      fraction = n / total,
      label = paste0(
        n,
        " (",
        round(fraction * 100),
        "%)"
      )
    ) %>%
    ungroup() %>%
    mutate(
      priority_class = factor(
        priority_class,
        levels = priority_levels
      )
    )
  
}

composition_original <- summarise_priority(
  priority_original
)

composition_sensitivity <- summarise_priority(
  priority_sensitivity
)

write_csv(
  composition_original,
  file.path(
    OUT_DIR,
    "07_composition_original.csv"
  )
)

write_csv(
  composition_sensitivity,
  file.path(
    OUT_DIR,
    "08_composition_sensitivity.csv"
  )
)


# ============================================================
# 10. Jaccard calculations
# ============================================================

calculate_jaccard <- function(
    df,
    fractions = TOP_FRACTIONS,
    n_random = N_RANDOM) {
  
  result <- list()
  idx <- 1L
  
  for (nu in c("TN", "TP")) {
    
    x <- df %>%
      filter(nutrient == nu) %>%
      arrange(STAID)
    
    N <- nrow(x)
    
    for (q in fractions) {
      
      k <- max(
        1L,
        ceiling(q * N)
      )
      
      top_pcs <- x %>%
        arrange(
          desc(PCS),
          STAID
        ) %>%
        slice_head(n = k) %>%
        pull(STAID)
      
      top_load <- x %>%
        arrange(
          desc(mean_annual_load),
          STAID
        ) %>%
        slice_head(n = k) %>%
        pull(STAID)
      
      shared <- length(
        intersect(
          top_pcs,
          top_load
        )
      )
      
      observed_jaccard <- shared /
        (2 * k - shared)
      
      random_values <- replicate(
        n_random,
        {
          s1 <- sample.int(
            N,
            size = k,
            replace = FALSE
          )
          
          s2 <- sample.int(
            N,
            size = k,
            replace = FALSE
          )
          
          shared_random <- length(
            intersect(s1, s2)
          )
          
          shared_random /
            (2 * k - shared_random)
        }
      )
      
      result[[idx]] <- data.frame(
        nutrient = nu,
        top_fraction = q,
        N = N,
        k = k,
        shared = shared,
        Jaccard = observed_jaccard,
        random_mean = mean(
          random_values
        ),
        random_low = unname(
          quantile(
            random_values,
            0.025
          )
        ),
        random_high = unname(
          quantile(
            random_values,
            0.975
          )
        ),
        random_p = (
          1 + sum(
            random_values >= observed_jaccard
          )
        ) / (n_random + 1)
      )
      
      idx <- idx + 1L
      
    }
  }
  
  bind_rows(result)
  
}

set.seed(RANDOM_SEED)

jaccard_original <- calculate_jaccard(
  d_original
)

set.seed(RANDOM_SEED)

jaccard_sensitivity <- calculate_jaccard(
  d_sensitivity
)

write_csv(
  jaccard_original,
  file.path(
    OUT_DIR,
    "09_jaccard_original.csv"
  )
)

write_csv(
  jaccard_sensitivity,
  file.path(
    OUT_DIR,
    "10_jaccard_sensitivity.csv"
  )
)

jaccard_comparison <- bind_rows(
  jaccard_original %>%
    mutate(scenario = "Original"),
  jaccard_sensitivity %>%
    mutate(scenario = "TN minimum removed")
)

write_csv(
  jaccard_comparison,
  file.path(
    OUT_DIR,
    "11_jaccard_comparison.csv"
  )
)


# ============================================================
# 11. Plotting theme
# ============================================================

theme_s4 <- theme_bw(
  base_size = 12
) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(
      face = "bold",
      size = 12
    ),
    plot.subtitle = element_text(
      size = 10
    ),
    legend.position = "bottom",
    strip.text = element_text(
      face = "bold"
    )
  )


# ============================================================
# 12. Figure creation function
# ============================================================

make_figure_s4 <- function(
    df,
    correlations,
    composition,
    jaccard,
    scenario_label) {
  
  # ----------------------------------------------------------
  # (a) TN
  # ----------------------------------------------------------
  
  tn_stats <- correlations %>%
    filter(nutrient == "TN")
  
  p_a <- df %>%
    filter(nutrient == "TN") %>%
    ggplot(
      aes(
        x = PCS,
        y = mean_annual_load
      )
    ) +
    geom_point(
      color = COLORS[["TN"]],
      alpha = 0.85,
      size = 2.5
    ) +
    scale_y_log10(
      labels = label_number(
        big.mark = " "
      )
    ) +
    labs(
      title = "(a) TN: network influence and load",
      subtitle = sprintf(
        "Spearman rho = %.3f; n = %d",
        tn_stats$spearman_rho,
        tn_stats$n
      ),
      x = "Global PCS (TN)",
      y = "Mean annual TN load [original unit]"
    ) +
    theme_s4
  
  # ----------------------------------------------------------
  # (b) TP
  # ----------------------------------------------------------
  
  tp_stats <- correlations %>%
    filter(nutrient == "TP")
  
  p_b <- df %>%
    filter(nutrient == "TP") %>%
    ggplot(
      aes(
        x = PCS,
        y = mean_annual_load
      )
    ) +
    geom_point(
      color = COLORS[["TP"]],
      alpha = 0.85,
      size = 2.5
    ) +
    scale_y_log10(
      labels = label_number(
        big.mark = " "
      )
    ) +
    labs(
      title = "(b) TP: network influence and load",
      subtitle = sprintf(
        "Spearman rho = %.3f; n = %d",
        tp_stats$spearman_rho,
        tp_stats$n
      ),
      x = "Global PCS (TP)",
      y = "Mean annual TP load [original unit]"
    ) +
    theme_s4
  
  # ----------------------------------------------------------
  # (c) Priority classes
  # ----------------------------------------------------------
  
  p_c <- ggplot(
    composition,
    aes(
      x = nutrient,
      y = fraction,
      fill = priority_class
    )
  ) +
    geom_col(
      width = 0.62,
      color = "white",
      linewidth = 0.35
    ) +
    geom_text(
      aes(
        label = if_else(
          n > 0,
          label,
          ""
        )
      ),
      position = position_stack(
        vjust = 0.5
      ),
      size = 3.5
    ) +
    scale_fill_manual(
      values = PRIORITY_COLORS,
      drop = FALSE
    ) +
    scale_y_continuous(
      labels = percent_format(
        accuracy = 1
      ),
      limits = c(0, 1),
      expand = expansion(
        mult = c(0, 0.02)
      )
    ) +
    labs(
      title = "(c) Agreement of high-priority stations",
      subtitle = paste0(
        "Top 20% by PCS and annual load"
      ),
      x = NULL,
      y = "Share of matched stations",
      fill = NULL
    ) +
    theme_s4
  
  # ----------------------------------------------------------
  # (d) Jaccard overlap
  # ----------------------------------------------------------
  
  p_d <- ggplot(
    jaccard,
    aes(
      x = top_fraction * 100,
      y = Jaccard,
      color = nutrient,
      group = nutrient
    )
  ) +
    geom_ribbon(
      aes(
        ymin = random_low,
        ymax = random_high,
        fill = nutrient
      ),
      alpha = 0.12,
      color = NA
    ) +
    geom_line(
      aes(
        y = random_mean
      ),
      linetype = "dotted",
      linewidth = 0.8
    ) +
    geom_line(
      linewidth = 1.1
    ) +
    geom_point(
      size = 2.7
    ) +
    scale_color_manual(
      values = COLORS
    ) +
    scale_fill_manual(
      values = COLORS,
      guide = "none"
    ) +
    scale_x_continuous(
      breaks = c(
        10, 20, 30
      )
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(
        0, 1, 0.2
      )
    ) +
    labs(
      title = "(d) Top-ranked station overlap",
      subtitle = paste0(
        "Solid: observed | dotted: random mean | ",
        "shaded: random 95% range"
      ),
      x = "Top station fraction (%)",
      y = "Jaccard index",
      color = NULL
    ) +
    theme_s4
  
  # ----------------------------------------------------------
  # Assemble
  # ----------------------------------------------------------
  
  fig <- (
    (p_a | p_b) /
      (p_c | p_d)
  ) +
    plot_annotation(
      title = paste0(
        "Figure S4 | PCS versus NTN nutrient loads"
      ),
      subtitle = scenario_label,
      caption = paste0(
        "Source: 06_unified_station_metrics.csv. ",
        "Original load units and exact 11-year ",
        "observation period require verification. ",
        "Random bands are reference distributions, ",
        "not confidence intervals."
      )
    )
  
  list(
    figure = fig,
    panels = list(
      a = p_a,
      b = p_b,
      c = p_c,
      d = p_d
    )
  )
  
}


# ============================================================
# 13. Produce original Figure S4
# ============================================================

original_fig <- make_figure_s4(
  df = d_original,
  correlations = corr_original,
  composition = composition_original,
  jaccard = jaccard_original,
  scenario_label = paste0(
    "Fig. 4 unified station cohort | ",
    "original station set"
  )
)

ggsave(
  file.path(
    OUT_DIR,
    "Figure_S4_all_stations.pdf"
  ),
  original_fig$figure,
  width = 15,
  height = 10.5,
  bg = "white"
)

ggsave(
  file.path(
    OUT_DIR,
    "Figure_S4_all_stations.png"
  ),
  original_fig$figure,
  width = 15,
  height = 10.5,
  dpi = 400,
  bg = "white"
)


# ============================================================
# 14. Produce sensitivity Figure S4
# ============================================================

sensitivity_fig <- make_figure_s4(
  df = d_sensitivity,
  correlations = corr_sensitivity,
  composition = composition_sensitivity,
  jaccard = jaccard_sensitivity,
  scenario_label = paste0(
    "Sensitivity analysis | ",
    "one lowest-load TN station excluded"
  )
)

ggsave(
  file.path(
    OUT_DIR,
    "Figure_S4_TN_lowest_load_removed.pdf"
  ),
  sensitivity_fig$figure,
  width = 15,
  height = 10.5,
  bg = "white"
)

ggsave(
  file.path(
    OUT_DIR,
    "Figure_S4_TN_lowest_load_removed.png"
  ),
  sensitivity_fig$figure,
  width = 15,
  height = 10.5,
  dpi = 400,
  bg = "white"
)


# ============================================================
# 15. Export independent sensitivity panels
# ============================================================

for (nm in names(sensitivity_fig$panels)) {
  
  ggsave(
    file.path(
      OUT_DIR,
      paste0(
        "Figure_S4_sensitivity_",
        nm,
        ".pdf"
      )
    ),
    sensitivity_fig$panels[[nm]],
    width = 7.5,
    height = 5.2,
    bg = "white"
  )
  
}


# ============================================================
# 16. Final outputs
# ============================================================

cat("\n========================================\n")
cat("Figure S4 completed.\n")
cat("Output directory:\n", OUT_DIR, "\n")
cat("========================================\n")

cat("\nRemoved TN station:\n")
print(tn_min, width = Inf)

cat("\nSpearman before / after:\n")
# 恢复默认缺失值打印设置
options(na.print = "NA")

# 转换为 tibble，确保使用正确的打印方法
jaccard_comparison <- tibble::as_tibble(
  jaccard_comparison
)

print(jaccard_comparison, n = Inf, width = Inf)

cat("\nPriority QA:\n")
print(
  qa_priority,
  width = Inf
)

cat("\nJaccard before / after:\n")
print(
  jaccard_comparison,
  n = Inf
)

print(sensitivity_fig$figure)