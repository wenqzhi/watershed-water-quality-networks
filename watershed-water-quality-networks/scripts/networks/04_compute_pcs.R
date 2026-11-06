source("config.R")

# ============================================================
# RECOMPUTE MULTI-STEP CASCADE PCS
# Based on existing Fusion-Net weighted matrices
#
# Nutrients: TN / TP
# Period: 2011–2023
#
# IMPORTANT:
# This is a reproducible weighted-walk cascade
# score, not yet verified as the original PCS.
#
# No TE / Phys-Net / Fusion-Net recalculation.
# No arbitrary physical-lag assignment.
# ============================================================


# ============================================================
# 0. PACKAGES
# ============================================================

PKGS <- c(
  "readr",
  "dplyr",
  "tidyr",
  "purrr",
  "ggplot2",
  "tibble"
)

for (p in PKGS) {
  
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

invisible(
  lapply(
    PKGS,
    library,
    character.only = TRUE
  )
)


# ============================================================
# 1. DIRECTORIES
# ============================================================

ROOT <- get_data_root()

FUSION_DIR <- file.path(
  ROOT,
  "Fusion_DualLag_2011_2023"
)

OUT_DIR <- file.path(
  ROOT,
  "Fig4_PCS_recomputed"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. MODEL SETTINGS
# ============================================================

H_MAIN <- 5L

ALPHA_MAIN <- 0.85

H_TEST <- c(3L, 5L, 8L, 12L)

ALPHA_TEST <- c(
  0.50,
  0.70,
  0.85,
  0.95
)

HOTSPOT_SHARE <- 0.20

# Do not use pad/substring to change station IDs.
# Preserve 8-, 9-, and 10-digit IDs exactly.

clean_id <- function(x) {
  
  trimws(
    as.character(x)
  )
}


# ============================================================
# 3. READ FULL FUSION WEIGHT MATRIX
# ============================================================

read_fusion_matrix <- function(nutrient) {
  
  path <- file.path(
    FUSION_DIR,
    nutrient,
    "Fusion_weight_matrix_full.csv"
  )
  
  if (!file.exists(path)) {
    stop("Missing matrix: ", path)
  }
  
  raw <- readr::read_csv(
    path,
    
    col_types = readr::cols(
      STAID = readr::col_character(),
      .default = readr::col_double()
    ),
    
    show_col_types = FALSE
  )
  
  if (!"STAID" %in% names(raw)) {
    stop("Matrix missing STAID column.")
  }
  
  row_ids <- clean_id(raw$STAID)
  
  column_ids <- clean_id(
    setdiff(names(raw), "STAID")
  )
  
  if (anyDuplicated(row_ids)) {
    stop("Duplicate matrix row IDs.")
  }
  
  if (anyDuplicated(column_ids)) {
    stop("Duplicate matrix column IDs.")
  }
  
  if (
    !setequal(
      row_ids,
      column_ids
    )
  ) {
    stop(
      "Matrix row/column stations do not match: ",
      nutrient
    )
  }
  
  # Align columns to the same ordering as rows.
  
  raw <- raw %>%
    select(
      all_of(c("STAID", row_ids))
    )
  
  W <- as.matrix(
    raw[, -1, drop = FALSE]
  )
  
  storage.mode(W) <- "double"
  
  rownames(W) <- row_ids
  colnames(W) <- row_ids
  
  if (
    any(!is.finite(W))
  ) {
    stop(
      "Non-finite Fusion weights in ",
      nutrient
    )
  }
  
  if (
    any(W < 0)
  ) {
    stop(
      "Negative Fusion weights in ",
      nutrient
    )
  }
  
  diag(W) <- 0
  
  cat(
    "\n",
    nutrient,
    "matrix dimension:",
    nrow(W),
    "x",
    ncol(W),
    "\n"
  )
  
  cat(
    "Positive edges:",
    sum(W > 0),
    "\n"
  )
  
  W
}


# ============================================================
# 4. MULTI-STEP CASCADE SCORE
# ============================================================

# Step 1:
# A = W / max(row sums of W)
#
# Step 2:
# P_i = sum_{h=1}^H alpha^(h-1)
#       rowSum(A^h)_i
#
# This is weighted directed-walk influence.
#
# Importantly, do not normalize each row of W
# separately to 1: that would make many nodes
# have identical total outgoing probability.
#
# Uses global normalization instead.


calculate_pcs <- function(
    W,
    H = 5L,
    alpha = 0.85
) {
  
  stopifnot(
    is.matrix(W),
    nrow(W) == ncol(W),
    H >= 1,
    alpha >= 0,
    alpha <= 1
  )
  
  n <- nrow(W)
  
  out_strength <- rowSums(W)
  
  max_out <- max(out_strength)
  
  if (
    !is.finite(max_out) ||
    max_out <= 0
  ) {
    stop(
      "Matrix contains no positive outgoing weights."
    )
  }
  
  A <- W / max_out
  
  # Recursive calculation of:
  #
  # A^h %*% 1
  #
  # Avoids explicitly constructing A^h.
  #
  # v_0 = 1
  # v_h = A %*% v_(h-1)
  
  v <- rep(1, n)
  
  pcs_raw <- rep(0, n)
  
  step_results <- matrix(
    0,
    nrow = n,
    ncol = H
  )
  
  for (h in seq_len(H)) {
    
    v <- as.vector(
      A %*% v
    )
    
    weighted_step <- (
      alpha^(h - 1)
    ) * v
    
    step_results[, h] <- weighted_step
    
    pcs_raw <- pcs_raw + weighted_step
  }
  
  # Normalize scores across nodes.
  #
  # Sum PCS_share = 1.
  #
  # This is a share of computed structural
  # cascade scores, NOT a share of nutrient load.
  
  total_pcs <- sum(pcs_raw)
  
  pcs_share <- if (
    total_pcs > 0
  ) {
    
    pcs_raw / total_pcs
    
  } else {
    
    rep(NA_real_, n)
  }
  
  output <- tibble(
    STAID = rownames(W),
    
    PCS_recomputed = pcs_raw,
    
    PCS_share = pcs_share,
    
    direct_component = step_results[, 1],
    
    indirect_component = (
      pcs_raw - step_results[, 1]
    ),
    
    out_strength = out_strength,
    
    in_strength = colSums(W),
    
    out_degree = rowSums(W > 0),
    
    in_degree = colSums(W > 0)
  )
  
  # Store every step so we can evaluate
  # the contribution of indirect walks.
  
  for (h in seq_len(H)) {
    
    output[[paste0(
      "step_", h
    )]] <- step_results[, h]
  }
  
  output
}


# ============================================================
# 5. READ EXISTING NODE COORDINATES
# ============================================================

read_existing_nodes <- function(nutrient) {
  
  path <- file.path(
    FUSION_DIR,
    nutrient,
    "Fusion_node_metrics.csv"
  )
  
  if (!file.exists(path)) {
    return(
      tibble(
        STAID = character(),
        lon = numeric(),
        lat = numeric()
      )
    )
  }
  
  raw <- readr::read_csv(
    path,
    
    col_types = readr::cols(
      STAID = readr::col_character(),
      .default = readr::col_guess()
    ),
    
    show_col_types = FALSE
  )
  
  raw %>%
    
    transmute(
      STAID = clean_id(STAID),
      
      lon = if (
        "lon" %in% names(raw)
      ) {
        as.numeric(lon)
      } else {
        NA_real_
      },
      
      lat = if (
        "lat" %in% names(raw)
      ) {
        as.numeric(lat)
      } else {
        NA_real_
      }
    ) %>%
    
    distinct(
      STAID,
      .keep_all = TRUE
    )
}


# ============================================================
# 6. CALCULATE TN AND TP PCS
# ============================================================

PCS_RESULTS <- list()

WEIGHT_MATRICES <- list()

for (nu in c("TN", "TP")) {
  
  cat(
    "\n====================================\n"
  )
  
  cat(
    "CALCULATING PCS:",
    nu,
    "\n"
  )
  
  cat(
    "====================================\n"
  )
  
  W <- read_fusion_matrix(nu)
  
  WEIGHT_MATRICES[[nu]] <- W
  
  pcs_nu <- calculate_pcs(
    W,
    H = H_MAIN,
    alpha = ALPHA_MAIN
  )
  
  old_nodes <- read_existing_nodes(nu)
  
  pcs_nu <- pcs_nu %>%
    
    left_join(
      old_nodes,
      by = "STAID"
    ) %>%
    
    mutate(
      nutrient = nu,
      H = H_MAIN,
      alpha = ALPHA_MAIN,
      .before = 1
    ) %>%
    
    arrange(
      desc(PCS_recomputed)
    ) %>%
    
    mutate(
      PCS_rank = row_number(),
      
      PCS_percentile = percent_rank(
        PCS_recomputed
      ),
      
      top20_PCS = (
        PCS_rank <= ceiling(
          n() * HOTSPOT_SHARE
        )
      )
    )
  
  PCS_RESULTS[[nu]] <- pcs_nu
  
  readr::write_csv(
    pcs_nu,
    
    file.path(
      OUT_DIR,
      paste0(
        nu,
        "_Global_PCS_recomputed.csv"
      )
    )
  )
  
  cat("\nTop 10 stations:\n")
  
  print(
    pcs_nu %>%
      select(
        STAID,
        PCS_rank,
        PCS_recomputed,
        PCS_share,
        direct_component,
        indirect_component
      ) %>%
      head(10),
    n = 10
  )
}


# ============================================================
# 7. COMBINED TN / TP PCS TABLE
# ============================================================

pcs_all <- bind_rows(
  PCS_RESULTS
)

readr::write_csv(
  pcs_all,
  
  file.path(
    OUT_DIR,
    "Global_PCS_TN_TP_recomputed.csv"
  )
)


# ============================================================
# 8. CUMULATIVE PCS CONTRIBUTION
# ============================================================

pcs_curve <- pcs_all %>%
  
  group_by(nutrient) %>%
  
  arrange(
    desc(PCS_recomputed),
    .by_group = TRUE
  ) %>%
  
  mutate(
    station_rank = row_number(),
    
    station_fraction = (
      station_rank / n()
    ),
    
    cumulative_PCS = cumsum(
      PCS_share
    )
  ) %>%
  
  ungroup()

readr::write_csv(
  pcs_curve,
  
  file.path(
    OUT_DIR,
    "PCS_cumulative_curves.csv"
  )
)


# ============================================================
# 9. CONCENTRATION STATISTICS
# ============================================================

pcs_summary <- pcs_curve %>%
  
  group_by(nutrient) %>%
  
  summarise(
    n_stations = n(),
    
    n_top20 = ceiling(
      n() * HOTSPOT_SHARE
    ),
    
    top20_PCS_share = sum(
      PCS_share[
        station_rank <=
          ceiling(
            n() * HOTSPOT_SHARE
          )
      ]
    ),
    
    stations_for_80pct = min(
      station_rank[
        cumulative_PCS >= 0.80
      ]
    ),
    
    station_fraction_for_80pct = (
      stations_for_80pct / n()
    ),
    
    mean_PCS = mean(
      PCS_recomputed
    ),
    
    median_PCS = median(
      PCS_recomputed
    ),
    
    max_PCS = max(
      PCS_recomputed
    ),
    
    .groups = "drop"
  )

readr::write_csv(
  pcs_summary,
  
  file.path(
    OUT_DIR,
    "PCS_concentration_summary.csv"
  )
)

cat("\nPCS concentration results:\n")

print(
  pcs_summary,
  n = Inf
)


# ============================================================
# 10. PARAMETER SENSITIVITY
# ============================================================

parameter_grid <- tidyr::crossing(
  H = H_TEST,
  alpha = ALPHA_TEST
)

sensitivity_results <- list()

for (nu in c("TN", "TP")) {
  
  W <- WEIGHT_MATRICES[[nu]]
  
  for (r in seq_len(
    nrow(parameter_grid)
  )) {
    
    h <- parameter_grid$H[r]
    
    a <- parameter_grid$alpha[r]
    
    temp <- calculate_pcs(
      W,
      H = h,
      alpha = a
    )
    
    temp <- temp %>%
      
      transmute(
        nutrient = nu,
        
        STAID,
        
        H = h,
        
        alpha = a,
        
        PCS = PCS_recomputed
      )
    
    sensitivity_results[[
      length(sensitivity_results) + 1
    ]] <- temp
  }
}

pcs_sensitivity <- bind_rows(
  sensitivity_results
)

readr::write_csv(
  pcs_sensitivity,
  
  file.path(
    OUT_DIR,
    "PCS_parameter_sensitivity.csv"
  )
)


# ============================================================
# 11. SENSITIVITY METRICS
# ============================================================

baseline <- pcs_all %>%
  
  select(
    nutrient,
    STAID,
    baseline_PCS = PCS_recomputed,
    baseline_top20 = top20_PCS
  )

sensitivity_comparison <- pcs_sensitivity %>%
  
  left_join(
    baseline,
    
    by = c(
      "nutrient",
      "STAID"
    )
  ) %>%
  
  group_by(
    nutrient,
    H,
    alpha
  ) %>%
  
  mutate(
    scenario_rank = min_rank(
      desc(PCS)
    ),
    
    scenario_top20 = (
      scenario_rank <= ceiling(
        n() * HOTSPOT_SHARE
      )
    )
  ) %>%
  
  summarise(
    n_stations = n(),
    
    rank_spearman = cor(
      PCS,
      baseline_PCS,
      method = "spearman"
    ),
    
    top20_overlap = (
      sum(
        scenario_top20 &
          baseline_top20
      ) /
        sum(baseline_top20)
    ),
    
    .groups = "drop"
  )

readr::write_csv(
  sensitivity_comparison,
  
  file.path(
    OUT_DIR,
    "PCS_sensitivity_metrics.csv"
  )
)

cat("\nPCS sensitivity results:\n")

print(
  sensitivity_comparison,
  n = Inf
)


# ============================================================
# 12. DIRECT VS INDIRECT INFLUENCE
# ============================================================

pcs_components <- pcs_all %>%
  
  mutate(
    indirect_fraction = if_else(
      PCS_recomputed > 0,
      
      indirect_component /
        PCS_recomputed,
      
      NA_real_
    )
  ) %>%
  
  select(
    nutrient,
    STAID,
    PCS_rank,
    PCS_recomputed,
    direct_component,
    indirect_component,
    indirect_fraction
  )

readr::write_csv(
  pcs_components,
  
  file.path(
    OUT_DIR,
    "PCS_direct_indirect_components.csv"
  )
)


# ============================================================
# 13. FIGURE: CUMULATIVE PCS
# ============================================================

p_cumulative <- ggplot(
  pcs_curve,
  
  aes(
    x = station_fraction,
    y = cumulative_PCS,
    color = nutrient
  )
) +
  
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "grey65"
  ) +
  
  geom_line(
    linewidth = 1.2
  ) +
  
  geom_vline(
    xintercept = HOTSPOT_SHARE,
    linetype = "dotted",
    color = "grey50"
  ) +
  
  scale_color_manual(
    values = c(
      TN = "#277F9C",
      TP = "#C7754F"
    )
  ) +
  
  scale_x_continuous(
    labels = scales::label_percent(),
    limits = c(0, 1)
  ) +
  
  scale_y_continuous(
    labels = scales::label_percent(),
    limits = c(0, 1)
  ) +
  
  labs(
    title = paste0(
      "Concentration of Weighted ",
      "Network Cascade Influence"
    ),
    
    subtitle = paste0(
      "H = ",
      H_MAIN,
      "; attenuation = ",
      ALPHA_MAIN
    ),
    
    x = "Cumulative fraction of stations ranked by PCS",
    y = "Cumulative fraction of PCS",
    color = NULL
  ) +
  
  theme_bw(base_size = 11) +
  
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )


# ============================================================
# 14. FIGURE: DIRECT / INDIRECT COMPONENTS
# ============================================================

p_components <- ggplot(
  pcs_components,
  
  aes(
    x = PCS_rank,
    y = indirect_fraction,
    color = nutrient
  )
) +
  
  geom_point(
    alpha = 0.75,
    size = 1.8
  ) +
  
  geom_smooth(
    method = "loess",
    se = FALSE
  ) +
  
  facet_wrap(
    ~nutrient,
    scales = "free_x"
  ) +
  
  scale_color_manual(
    values = c(
      TN = "#277F9C",
      TP = "#C7754F"
    ),
    guide = "none"
  ) +
  
  scale_y_continuous(
    labels = scales::label_percent()
  ) +
  
  labs(
    title = "Indirect contribution to network cascade score",
    
    x = "Station rank by PCS",
    
    y = "Indirect / total PCS"
  ) +
  
  theme_bw(base_size = 11)


# ============================================================
# 15. FIGURE: PCS SENSITIVITY
# ============================================================

p_sensitivity <- ggplot(
  sensitivity_comparison,
  
  aes(
    x = factor(H),
    y = rank_spearman,
    fill = factor(alpha)
  )
) +
  
  geom_col(
    position = position_dodge(
      width = 0.8
    ),
    width = 0.7
  ) +
  
  facet_wrap(
    ~nutrient
  ) +
  
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  
  labs(
    title = "PCS ranking sensitivity",
    
    x = "Maximum walk length (H)",
    
    y = "Spearman correlation with baseline",
    
    fill = "Attenuation"
  ) +
  
  theme_bw(base_size = 11)


# ============================================================
# 16. SAVE FIGURES
# ============================================================

ggsave(
  file.path(
    OUT_DIR,
    "PCS_cumulative_curves.pdf"
  ),
  
  p_cumulative,
  
  width = 8,
  height = 6,
  bg = "white"
)

ggsave(
  file.path(
    OUT_DIR,
    "PCS_direct_indirect_components.pdf"
  ),
  
  p_components,
  
  width = 10,
  height = 5,
  bg = "white"
)

ggsave(
  file.path(
    OUT_DIR,
    "PCS_parameter_sensitivity.pdf"
  ),
  
  p_sensitivity,
  
  width = 9,
  height = 5.5,
  bg = "white"
)


# ============================================================
# 17. FINAL QA
# ============================================================

cat("\n====================================\n")
cat("PCS RECOMPUTATION FINISHED\n")
cat("====================================\n")

cat("\nNumber of nodes:\n")

print(
  pcs_all %>%
    count(nutrient)
)

cat("\nPCS values:\n")

print(
  pcs_all %>%
    group_by(nutrient) %>%
    summarise(
      min_PCS = min(PCS_recomputed),
      median_PCS = median(PCS_recomputed),
      max_PCS = max(PCS_recomputed),
      sum_share = sum(PCS_share),
      .groups = "drop"
    )
)

cat("\nCumulative contribution:\n")

print(pcs_summary)

cat("\nOutput directory:\n")
cat(OUT_DIR, "\n")

print(p_cumulative)
