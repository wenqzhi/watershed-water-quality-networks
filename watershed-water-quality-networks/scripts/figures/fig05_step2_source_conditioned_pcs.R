# ============================================================
# FIG. 5 — STEP 5.2
# SOURCE-CONDITIONED PCS AND TARGET RANKING
#
# Complementary monitoring information:
# F-Net / Global PCS / NTN loading / Random
#
# This is a RETROSPECTIVE structural assessment.
# NOT an independent out-of-sample forecast.
# ============================================================

library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

# ============================================================
# 1. PATHS
# ============================================================

source("config.R")
ROOT <- get_data_root()

STEP51 <- file.path(
  ROOT,
  "Fig5_anomaly_trigger",
  "Step5_1_event_catalog"
)

OUT <- file.path(
  ROOT,
  "Fig5_anomaly_trigger",
  "Step5_2_source_PCS"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

EVENT_FILE <- file.path(
  STEP51,
  "05_all_real_anomaly_triggers.csv"
)

CATALOG_FILE <- file.path(
  STEP51,
  "08_Fig5_candidate_event_catalog.csv"
)

OBS_FILE <- file.path(
  STEP51,
  "06_event_target_observation_matrix.csv"
)

PCS_FILE <- file.path(
  ROOT,
  "Fig4_PCS_recomputed",
  "Global_PCS_TN_TP_recomputed.csv"
)

LOAD_FILE <- file.path(
  ROOT,
  "Fig4_final_independent",
  "03_NTN_mean_annual_load.csv"
)

FUSION_DIR <- file.path(
  ROOT,
  "Fusion_DualLag_2011_2023"
)

# ============================================================
# 2. SETTINGS
# ============================================================

H <- 5L
ALPHA <- 0.85

TOP_K <- c(5L, 10L, 20L)

Z_THRESHOLD <- 1.5

# Require complete 0–3 month observation window.
REQUIRE_FULL_WINDOW <- TRUE

# If Global PCS or NTN load unavailable,
# those strategies rank only valid candidates.

# ============================================================
# 3. HELPERS
# ============================================================

clean_id <- function(x) {
  x <- trimws(as.character(x))
  sub("\\.0+$", "", x)
}

save_csv <- function(x, name) {
  readr::write_csv(
    x,
    file.path(OUT, name),
    na = ""
  )
}

read_csv_ids <- function(path, id_cols) {
  
  if (!file.exists(path)) {
    stop("Missing file: ", path)
  }
  
  spec <- readr::cols(.default = readr::col_guess())
  
  for (nm in id_cols) {
    spec$cols[[nm]] <- readr::col_character()
  }
  
  readr::read_csv(
    path,
    col_types = spec,
    show_col_types = FALSE
  )
}

# ============================================================
# 4. READ EVENTS
# ============================================================

events <- read_csv_ids(
  EVENT_FILE,
  c("source")
) %>%
  mutate(
    source = clean_id(source),
    nutrient = toupper(nutrient),
    trigger_date = as.Date(trigger_date)
  )

catalog <- read_csv_ids(
  CATALOG_FILE,
  c("source")
) %>%
  mutate(
    source = clean_id(source)
  )

obs <- read_csv_ids(
  OBS_FILE,
  c("source", "target")
) %>%
  mutate(
    source = clean_id(source),
    target = clean_id(target),
    nutrient = toupper(nutrient),
    target_date = as.Date(target_date),
    lag_months = as.integer(lag_months)
  )

pcs <- read_csv_ids(
  PCS_FILE,
  c("STAID")
) %>%
  transmute(
    STAID = clean_id(STAID),
    nutrient = toupper(nutrient),
    global_PCS = as.numeric(PCS_recomputed)
  )

# ============================================================
# 5. READ NTN LOAD BASELINE
# ============================================================

if (file.exists(LOAD_FILE)) {
  
  ntn_load <- read_csv_ids(
    LOAD_FILE,
    c("STAID")
  ) %>%
    transmute(
      STAID = clean_id(STAID),
      nutrient = toupper(nutrient),
      mean_annual_load = as.numeric(
        mean_annual_load
      )
    )
  
} else {
  
  warning(
    "NTN load baseline file not found. ",
    "Loading-based comparison will be omitted."
  )
  
  ntn_load <- tibble(
    STAID = character(),
    nutrient = character(),
    mean_annual_load = numeric()
  )
}

# ============================================================
# 6. READ F-NET
# ============================================================

read_fusion <- function(nu) {
  
  path <- file.path(
    FUSION_DIR,
    nu,
    "Fusion_all_edges.csv"
  )
  
  d <- read_csv_ids(
    path,
    c("from", "to")
  )
  
  d %>%
    transmute(
      nutrient = nu,
      from = clean_id(from),
      to = clean_id(to),
      weight = as.numeric(W_fusion)
    ) %>%
    filter(
      is.finite(weight),
      weight > 0,
      !is.na(from),
      !is.na(to),
      from != to
    ) %>%
    group_by(from, to, nutrient) %>%
    summarise(
      weight = max(weight),
      .groups = "drop"
    )
}

edges <- bind_rows(
  read_fusion("TN"),
  read_fusion("TP")
)

save_csv(edges, "01_FNet_edges_used.csv")

# ============================================================
# 7. CONSTRUCT TRANSITION MATRICES
# ============================================================

# P is row-normalized directed F-Net.
#
# Dangling nodes have zero transition rows:
# their walks terminate rather than restart.

make_transition <- function(nu) {
  
  nodes <- pcs %>%
    filter(nutrient == nu) %>%
    arrange(STAID) %>%
    pull(STAID)
  
  d <- edges %>%
    filter(nutrient == nu)
  
  N <- length(nodes)
  
  W <- matrix(
    0,
    nrow = N,
    ncol = N,
    dimnames = list(nodes, nodes)
  )
  
  i <- match(d$from, nodes)
  j <- match(d$to, nodes)
  
  ok <- !is.na(i) & !is.na(j)
  
  for (r in which(ok)) {
    W[i[r], j[r]] <- d$weight[r]
  }
  
  totals <- rowSums(W)
  
  P <- W
  
  active <- totals > 0
  
  P[active, ] <- sweep(
    W[active, , drop = FALSE],
    1,
    totals[active],
    "/"
  )
  
  list(
    nodes = nodes,
    W = W,
    P = P,
    active = active
  )
}

networks <- list(
  TN = make_transition("TN"),
  TP = make_transition("TP")
)

# ============================================================
# 8. SOURCE-CONDITIONED PCS
# ============================================================

# Compute:
# I_Aj = sum_l alpha^(l-1) * (P^l)_Aj
#
# No self-contribution in target ranking.
# This version retains all matrix walks,
# including possible repeated nodes.

source_influence <- function(nu, source) {
  
  net <- networks[[nu]]
  
  nodes <- net$nodes
  P <- net$P
  
  idx <- match(source, nodes)
  
  if (is.na(idx)) {
    return(NULL)
  }
  
  v <- rep(0, length(nodes))
  v[idx] <- 1
  
  influence <- rep(0, length(nodes))
  
  for (l in seq_len(H)) {
    
    v <- as.numeric(v %*% P)
    
    influence <- influence +
      ALPHA^(l - 1) * v
  }
  
  influence[idx] <- NA_real_
  
  tibble(
    nutrient = nu,
    source = source,
    target = nodes,
    source_PCS = influence
  ) %>%
    filter(target != source)
}

# ============================================================
# 9. CACHE BY UNIQUE SOURCE
# ============================================================

unique_sources <- events %>%
  distinct(nutrient, source)

influence_cache <- vector(
  "list",
  nrow(unique_sources)
)

for (k in seq_len(nrow(unique_sources))) {
  
  nu <- unique_sources$nutrient[k]
  so <- unique_sources$source[k]
  
  influence_cache[[k]] <- source_influence(
    nu, so
  )
}

influence_all <- bind_rows(influence_cache)

if (nrow(influence_all) == 0) {
  stop("No valid source-conditioned influence.")
}

save_csv(
  influence_all,
  "02_unique_source_target_influence.csv"
)

# ============================================================
# 10. EVENT-SPECIFIC TARGET RANKING
# ============================================================

rankings <- events %>%
  select(
    event_id,
    nutrient,
    source,
    trigger_date,
    source_Z
  ) %>%
  left_join(
    influence_all,
    by = c("nutrient", "source")
  ) %>%
  left_join(
    pcs %>%
      rename(
        target = STAID
      ),
    by = c("nutrient", "target")
  ) %>%
  left_join(
    ntn_load %>%
      rename(
        target = STAID
      ),
    by = c("nutrient", "target")
  ) %>%
  group_by(event_id) %>%
  mutate(
    rank_FNet = rank(
      -source_PCS,
      ties.method = "first",
      na.last = "keep"
    ),
    
    rank_global = rank(
      -global_PCS,
      ties.method = "first",
      na.last = "keep"
    ),
    
    rank_load = rank(
      -mean_annual_load,
      ties.method = "first",
      na.last = "keep"
    ),
    
    positive_reach = (
      is.finite(source_PCS) &
        source_PCS > 0
    )
  ) %>%
  ungroup()

# Strict QA
if (anyDuplicated(
  rankings[c("event_id", "target")]
)) {
  stop("Duplicate event-target ranking rows.")
}

save_csv(
  rankings,
  "03_all_event_target_rankings.csv"
)

# ============================================================
# 11. JOIN REAL OBSERVED ANOMALIES
# ============================================================

# Evaluate a fixed 0–3 month horizon.
# Lag-aware windows will be assigned in Step 5.3.

observed_window <- obs %>%
  group_by(
    event_id,
    nutrient,
    source,
    target
  ) %>%
  summarise(
    n_lags_present = n_distinct(lag_months),
    
    n_lags_observed = sum(target_observed),
    
    any_anomaly = any(
      target_anomaly_15,
      na.rm = TRUE
    ),
    
    first_anomaly_lag = if (
      any(target_anomaly_15, na.rm = TRUE)
    ) {
      min(lag_months[target_anomaly_15])
    } else {
      NA_integer_
    },
    
    .groups = "drop"
  ) %>%
  mutate(
    complete_window = (
      n_lags_present == 4 &
        n_lags_observed == 4
    )
  )

evaluation <- rankings %>%
  left_join(
    observed_window,
    by = c(
      "event_id",
      "nutrient",
      "source",
      "target"
    )
  )

if (REQUIRE_FULL_WINDOW) {
  
  evaluation <- evaluation %>%
    filter(complete_window)
}

save_csv(
  evaluation,
  "04_rankings_with_observed_outcomes.csv"
)

# ============================================================
# 12. TOP-K PERFORMANCE
# ============================================================

# Important:
# These are retrospective associations,
# not independent predictive performance.
#
# Each event is evaluated within the same
# candidate target set.

score_strategy <- function(d, column, k) {
  
  valid <- d %>%
    filter(
      is.finite(.data[[column]]),
      !is.na(any_anomaly)
    ) %>%
    arrange(.data[[column]], target)
  
  if (nrow(valid) < k) {
    return(tibble(
      n_eligible = nrow(valid),
      hits = NA_integer_,
      precision = NA_real_,
      recall = NA_real_
    ))
  }
  
  top <- head(valid, k)
  
  total_positive <- sum(valid$any_anomaly)
  
  hits <- sum(top$any_anomaly)
  
  tibble(
    n_eligible = nrow(valid),
    hits = hits,
    precision = hits / k,
    recall = if (total_positive > 0) {
      hits / total_positive
    } else {
      NA_real_
    }
  )
}

results <- list()
counter <- 1L

for (ev in unique(evaluation$event_id)) {
  
  d <- evaluation %>%
    filter(event_id == ev)
  
  nu <- unique(d$nutrient)
  
  for (k in TOP_K) {
    
    for (strategy in c(
      "rank_FNet",
      "rank_global",
      "rank_load"
    )) {
      
      r <- score_strategy(
        d, strategy, k
      )
      
      results[[counter]] <- tibble(
        event_id = ev,
        nutrient = nu,
        K = k,
        strategy = strategy
      ) %>%
        bind_cols(r)
      
      counter <- counter + 1L
    }
  }
}

performance <- bind_rows(results)

save_csv(
  performance,
  "05_event_level_TopK_performance.csv"
)

performance_summary <- performance %>%
  group_by(
    nutrient,
    K,
    strategy
  ) %>%
  summarise(
    n_evaluated_events =
      sum(is.finite(precision)),
    
    mean_precision = if (
      any(is.finite(precision))
    ) {
      mean(precision, na.rm = TRUE)
    } else {
      NA_real_
    },
    
    mean_recall = if (
      any(is.finite(recall))
    ) {
      mean(recall, na.rm = TRUE)
    } else {
      NA_real_
    },
    
    .groups = "drop"
  )

save_csv(
  performance_summary,
  "06_TopK_performance_summary.csv"
)

# ============================================================
# 13. SOURCE-INFLUENCE COVERAGE
# ============================================================

coverage <- rankings %>%
  group_by(
    event_id,
    nutrient,
    source
  ) %>%
  summarise(
    n_candidate_targets = n(),
    
    n_positive_influence = sum(
      positive_reach,
      na.rm = TRUE
    ),
    
    fraction_positive_influence =
      n_positive_influence /
      n_candidate_targets,
    
    max_source_PCS = if (
      any(is.finite(source_PCS))
    ) {
      max(source_PCS, na.rm = TRUE)
    } else {
      NA_real_
    },
    
    .groups = "drop"
  )

save_csv(
  coverage,
  "07_source_influence_coverage.csv"
)

# ============================================================
# 14. PLOT PRELIMINARY PERFORMANCE
# ============================================================

p <- performance_summary %>%
  filter(is.finite(mean_precision)) %>%
  mutate(
    strategy = recode(
      strategy,
      rank_FNet = "Source-conditioned F-Net",
      rank_global = "Global PCS",
      rank_load = "NTN loading"
    )
  ) %>%
  ggplot(
    aes(
      x = K,
      y = mean_precision,
      color = strategy,
      group = strategy
    )
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.4) +
  facet_wrap(
    ~nutrient,
    ncol = 2
  ) +
  scale_y_continuous(
    labels = scales::label_percent(),
    limits = c(0, 1)
  ) +
  scale_x_continuous(
    breaks = TOP_K
  ) +
  labs(
    title = "Retrospective anomaly attention: Top-K",
    subtitle = paste0(
      "Any observed target anomaly within ",
      "0–3 months; no lag-aware filtering yet"
    ),
    x = "Number of prioritized target stations",
    y = "Observed anomaly fraction within Top-K",
    color = "Ranking strategy"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  file.path(
    OUT,
    "Fig5S_preliminary_TopK.pdf"
  ),
  p,
  width = 11,
  height = 6,
  bg = "white"
)

# ============================================================
# 15. SUMMARY
# ============================================================

cat("\n==================================\n")
cat("FIG.5 STEP 5.2 COMPLETED\n")
cat("==================================\n")

cat("\nUnique sources:\n")
print(unique_sources %>% count(nutrient))

cat("\nSource influence coverage:\n")
print(
  coverage %>%
    group_by(nutrient) %>%
    summarise(
      median_fraction_positive =
        median(fraction_positive_influence),
      .groups = "drop"
    )
)

cat("\nPreliminary Top-K comparison:\n")
print(
  performance_summary,
  n = Inf
)

cat("\nOutput:\n", OUT, "\n")

print(p)