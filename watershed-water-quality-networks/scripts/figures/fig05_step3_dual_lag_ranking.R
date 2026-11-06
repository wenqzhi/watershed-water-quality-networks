# ============================================================
# FIGURE 5 — STEP 5.3
# TOP-K DUAL-LAG-AWARE MONITORING ATTENTION
#
# Main:
#   Top-5 Source-conditioned PCS
#
# Supplement:
#   Top-10 Source-conditioned PCS
#
# Lag rules:
#   R=1  -> physical travel time (days)
#   R=0 & significant TE -> TE lag (months)
#   Other -> undefined
#
# Multi-step:
#   representative highest-weight valid path,
#   maximum 5 edges, no repeated nodes.
#
# Evaluation:
#   observed target anomaly in suggested month(s)
#
# IMPORTANT:
#   retrospective attention evaluation,
#   not independent prediction or mass transport.
# ============================================================


# ============================================================
# 0. PACKAGES
# ============================================================

PKGS <- c(
  "readr", "dplyr", "tidyr",
  "tibble", "ggplot2", "lubridate"
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

BASE <- file.path(
  ROOT, "Fig5_anomaly_trigger"
)

STEP51 <- file.path(
  BASE, "Step5_1_event_catalog"
)

STEP52 <- file.path(
  BASE, "Step5_2_source_PCS"
)

OUT <- file.path(
  BASE, "Step5_3_dual_lag_TopK"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

RANK_FILE <- file.path(
  STEP52, "03_all_event_target_rankings.csv"
)

OBS_FILE <- file.path(
  STEP51, "06_event_target_observation_matrix.csv"
)

FUSION_DIR <- file.path(
  ROOT, "Fusion_DualLag_2011_2023"
)

PNET_DIR <- file.path(
  ROOT, "PNet_2011_2023_BlockTracing"
)

TE_DIR <- file.path(
  ROOT, "TE_network_results"
)


# ============================================================
# 2. SETTINGS
# ============================================================

TOP_K <- c(5L, 10L)

MAX_PATH_EDGES <- 5L

TE_P <- 0.10

Z_THRESHOLD <- 1.5

MAX_EVAL_LAG <- 3L

# Approximate month length used only to convert
# hydraulic travel days to monthly attention bins.
DAYS_PER_MONTH <- 30.44

# No extra fitted coefficients.


# ============================================================
# 3. HELPERS
# ============================================================

need <- function(path) {
  
  if (!file.exists(path)) {
    stop("Missing file:\n", path)
  }
}

clean_id <- function(x) {
  
  x <- trimws(as.character(x))
  x <- sub("\\.0+$", "", x)
  
  x[x %in% c("", "NA")] <- NA_character_
  
  x
}

as_bool <- function(x) {
  
  as.character(x) %in%
    c("TRUE", "True", "true", "T", "1")
}

save_csv <- function(d, name) {
  
  readr::write_csv(
    d,
    file.path(OUT, name),
    na = ""
  )
}

optional_num <- function(d, names_to_try) {
  
  found <- names_to_try[
    names_to_try %in% names(d)
  ]
  
  if (!length(found)) {
    return(rep(NA_real_, nrow(d)))
  }
  
  suppressWarnings(
    as.numeric(d[[found[1]]])
  )
}

find_input <- function(directory, filename) {
  
  candidate <- file.path(
    directory, filename
  )
  
  if (file.exists(candidate)) return(candidate)
  
  found <- list.files(
    directory,
    pattern = paste0(
      "^", gsub("\\.", "\\\\.", filename), "$"
    ),
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (!length(found)) return(NA_character_)
  
  found[1]
}


# ============================================================
# 4. READ STEP 5.2 RANKINGS
# ============================================================

need(RANK_FILE)
need(OBS_FILE)

rankings <- readr::read_csv(
  RANK_FILE,
  
  col_types = readr::cols(
    source = readr::col_character(),
    target = readr::col_character(),
    .default = readr::col_guess()
  ),
  
  show_col_types = FALSE
) %>%
  mutate(
    source = clean_id(source),
    target = clean_id(target),
    
    nutrient = toupper(nutrient),
    
    trigger_date = as.Date(trigger_date),
    
    rank_FNet = as.integer(rank_FNet),
    
    source_PCS = as.numeric(source_PCS)
  ) %>%
  filter(
    !is.na(rank_FNet),
    is.finite(source_PCS),
    source_PCS > 0
  )

obs <- readr::read_csv(
  OBS_FILE,
  
  col_types = readr::cols(
    source = readr::col_character(),
    target = readr::col_character(),
    .default = readr::col_guess()
  ),
  
  show_col_types = FALSE
) %>%
  mutate(
    source = clean_id(source),
    target = clean_id(target),
    
    nutrient = toupper(nutrient),
    
    lag_months = as.integer(lag_months),
    
    target_observed = as_bool(target_observed),
    
    target_anomaly_15 = as_bool(
      target_anomaly_15
    )
  )

if (anyDuplicated(
  rankings[c("event_id", "target")]
)) {
  stop("Duplicate event-target rankings.")
}

if (anyDuplicated(
  obs[c("event_id", "target", "lag_months")]
)) {
  stop("Duplicate event-target-lag observations.")
}


# ============================================================
# 5. READ ORIGINAL F-NET FULL EDGES
# ============================================================

read_fnet <- function(nu) {
  
  f <- file.path(
    FUSION_DIR,
    nu,
    "Fusion_all_edges.csv"
  )
  
  need(f)
  
  d <- readr::read_csv(
    f,
    
    col_types = readr::cols(
      from = readr::col_character(),
      to = readr::col_character(),
      .default = readr::col_guess()
    ),
    
    show_col_types = FALSE
  )
  
  tibble(
    nutrient = nu,
    
    from = clean_id(d$from),
    to = clean_id(d$to),
    
    W_fusion = as.numeric(d$W_fusion),
    
    R_fusion = optional_num(d, c("R")),
    
    physical_days_fusion = optional_num(
      d,
      c(
        "travel_days_medium",
        "travel_days",
        "tau_phys_days"
      )
    ),
    
    assigned_lag = optional_num(
      d,
      c("assigned_lag")
    ),
    
    lag_source = if (
      "lag_source" %in% names(d)
    ) {
      as.character(d$lag_source)
    } else {
      rep(NA_character_, nrow(d))
    },
    
    lag_unit = if (
      "lag_unit" %in% names(d)
    ) {
      as.character(d$lag_unit)
    } else {
      rep(NA_character_, nrow(d))
    }
  ) %>%
    filter(
      is.finite(W_fusion),
      W_fusion > 0,
      !is.na(from),
      !is.na(to),
      from != to
    )
}

fnet <- bind_rows(
  read_fnet("TN"),
  read_fnet("TP")
)

if (anyDuplicated(
  fnet[c("nutrient", "from", "to")]
)) {
  stop("Duplicate directed Fusion-Net edges.")
}


# ============================================================
# 6. READ P-NET PHYSICAL TRAVEL TIMES
# ============================================================

read_pnet <- function(nu) {
  
  file <- find_input(
    file.path(PNET_DIR, nu),
    "PNet_directed_reachable_pairs.csv"
  )
  
  if (is.na(file)) {
    return(tibble(
      nutrient = character(),
      from = character(),
      to = character(),
      physical_days_pnet = numeric()
    ))
  }
  
  d <- readr::read_csv(
    file,
    col_types = readr::cols(
      from = readr::col_character(),
      to = readr::col_character(),
      .default = readr::col_guess()
    ),
    show_col_types = FALSE
  )
  
  tibble(
    nutrient = nu,
    
    from = clean_id(d$from),
    to = clean_id(d$to),
    
    physical_days_pnet = optional_num(
      d,
      c(
        "travel_days_medium",
        "travel_days",
        "tau_phys_days"
      )
    )
  ) %>%
    distinct(
      nutrient, from, to,
      .keep_all = TRUE
    )
}

pnet <- bind_rows(
  read_pnet("TN"),
  read_pnet("TP")
)


# ============================================================
# 7. READ TE SIGNIFICANCE AND LAGS
# ============================================================

read_te <- function(nu) {
  
  file <- file.path(
    TE_DIR,
    nu,
    "main_2011_2023_all_edges.csv"
  )
  
  need(file)
  
  d <- readr::read_csv(
    file,
    show_col_types = FALSE
  )
  
  from_name <- intersect(
    c("from", "source"), names(d)
  )[1]
  
  to_name <- intersect(
    c("to", "target"), names(d)
  )[1]
  
  if (
    is.na(from_name) ||
    is.na(to_name)
  ) {
    stop("TE source/target columns missing.")
  }
  
  tibble(
    nutrient = nu,
    
    from = clean_id(d[[from_name]]),
    to = clean_id(d[[to_name]]),
    
    TE_p = optional_num(
      d, c("p_value", "p")
    ),
    
    TE_value = optional_num(
      d,
      c("TE_bias_corrected", "TE")
    ),
    
    TE_lag = optional_num(
      d,
      c("lag_months", "lag")
    )
  )
}

te <- bind_rows(
  read_te("TN"),
  read_te("TP")
)

if (anyDuplicated(
  te[c("nutrient", "from", "to")]
)) {
  stop("Duplicate TE directed pairs.")
}


# ============================================================
# 8. ASSIGN DUAL LAGS TO F-NET EDGES
# ============================================================

edge_lags <- fnet %>%
  left_join(
    pnet,
    by = c("nutrient", "from", "to")
  ) %>%
  left_join(
    te,
    by = c("nutrient", "from", "to")
  ) %>%
  mutate(
    # A verified P-Net reachable path can establish R=1.
    # Missing P-Net records never establish R=0.
    
    R = case_when(
      is.finite(R_fusion) ~ R_fusion,
      
      is.finite(physical_days_pnet) ~ 1,
      
      TRUE ~ NA_real_
    ),
    
    # Fusion file time has priority; fall back
    # to standalone P-Net.
    
    physical_days = coalesce(
      physical_days_fusion,
      physical_days_pnet
    ),
    
    # If original dual-lag data explicitly assign
    # a physical lag, use it when needed.
    
    physical_days = if_else(
      R == 1 &
        !is.finite(physical_days) &
        lag_source == "Physical" &
        lag_unit == "days",
      
      assigned_lag,
      
      physical_days
    ),
    
    TE_significant = (
      is.finite(TE_p) &
        TE_p < TE_P &
        is.finite(TE_value) &
        TE_value > 0
    ),
    
    edge_type = case_when(
      R == 1 ~ "Physical",
      
      R == 0 & TE_significant ~ "Statistical",
      
      TRUE ~ "Undefined"
    ),
    
    lag_in_months = case_when(
      edge_type == "Physical" &
        is.finite(physical_days) &
        physical_days >= 0 ~
        physical_days / DAYS_PER_MONTH,
      
      edge_type == "Statistical" &
        is.finite(TE_lag) &
        TE_lag >= 0 ~ TE_lag,
      
      TRUE ~ NA_real_
    )
  ) %>%
  mutate(
    lag_defined = is.finite(lag_in_months)
  )

save_csv(
  edge_lags,
  "01_FNet_dual_lag_edges.csv"
)

cat("\nDual-lag assignment:\n")

print(
  edge_lags %>%
    count(
      nutrient,
      edge_type,
      lag_defined
    ),
  n = Inf
)

if (!any(edge_lags$lag_defined)) {
  stop(
    "No F-Net edge has a valid dual lag. ",
    "Check P-Net travel days, TE lag fields ",
    "and verified R values."
  )
}


# ============================================================
# 9. TRANSITION WEIGHTS
# ============================================================

# Use the original F-Net W_fusion to form
# row-normalized transition probabilities.
#
# Normalize using ALL positive edges,
# including edges without assigned lag.
#
# Otherwise dropping undefined edges would
# artificially increase path probabilities.

edge_lags <- edge_lags %>%
  group_by(
    nutrient, from
  ) %>%
  mutate(
    P = W_fusion /
      sum(W_fusion)
  ) %>%
  ungroup()

save_csv(
  edge_lags,
  "02_FNet_transition_and_lags.csv"
)


# ============================================================
# 10. REPRESENTATIVE VALID PATH — FIXED
# ============================================================

best_paths_for_source <- function(
    nu, source, edge_data
) {
  
  d <- edge_data %>%
    filter(
      nutrient == nu,
      lag_defined,
      is.finite(P),
      P > 0,
      is.finite(lag_in_months),
      lag_in_months >= 0
    ) %>%
    select(
      from,
      to,
      P,
      lag_in_months,
      edge_type
    )
  
  if (nrow(d) == 0) {
    return(tibble())
  }
  
  source <- as.character(source)
  
  # Each path stores its complete visited-node sequence.
  frontier <- tibble(
    node = source,
    probability = 1,
    months = 0,
    n_edges = 0L,
    route = source,
    types = ""
  )
  
  results <- list()
  
  for (step in seq_len(MAX_PATH_EDGES)) {
    
    if (nrow(frontier) == 0) break
    
    # Expand every currently retained path.
    next_paths <- frontier %>%
      inner_join(
        d,
        by = c("node" = "from"),
        relationship = "many-to-many"
      )
    
    if (nrow(next_paths) == 0) break
    
    # Vectorized cycle check.
    # fixed=TRUE avoids regular-expression issues.
    #
    # Check each candidate target against the
    # visited nodes of its OWN route.
    
    is_cycle <- mapply(
      function(route_i, to_i) {
        
        if (
          is.na(route_i) ||
          is.na(to_i) ||
          !nzchar(route_i)
        ) {
          return(TRUE)
        }
        
        visited <- strsplit(
          route_i,
          ">",
          fixed = TRUE
        )[[1]]
        
        to_i %in% visited
      },
      next_paths$route,
      next_paths$to,
      USE.NAMES = FALSE
    )
    
    next_paths <- next_paths[
      !is_cycle,
      ,
      drop = FALSE
    ]
    
    if (nrow(next_paths) == 0) break
    
    next_paths <- next_paths %>%
      mutate(
        probability = probability * P,
        
        months = months + lag_in_months,
        
        n_edges = step,
        
        route = paste0(
          route, ">", to
        ),
        
        types = if_else(
          types == "",
          edge_type,
          paste(
            types,
            edge_type,
            sep = "+"
          )
        ),
        
        node = to
      ) %>%
      select(
        node,
        probability,
        months,
        n_edges,
        route,
        types
      )
    
    # Save all valid paths at the current step.
    results[[length(results) + 1L]] <- next_paths
    
    # Limit excessive growth:
    # retain strongest 5 candidates per
    # destination at each step.
    #
    # This is a computational approximation,
    # not a new physical parameter.
    
    frontier <- next_paths %>%
      group_by(node) %>%
      arrange(
        desc(probability),
        .by_group = TRUE
      ) %>%
      slice_head(n = 5) %>%
      ungroup()
  }
  
  if (length(results) == 0) {
    return(tibble())
  }
  
  all_paths <- bind_rows(results) %>%
    filter(
      node != source,
      is.finite(probability),
      probability > 0
    )
  
  if (nrow(all_paths) == 0) {
    return(tibble())
  }
  
  # Representative maximum-product path
  # for each target station.
  
  best <- all_paths %>%
    group_by(node) %>%
    arrange(
      desc(probability),
      n_edges,
      .by_group = TRUE
    ) %>%
    slice_head(n = 1) %>%
    ungroup()
  
  best %>%
    transmute(
      nutrient = nu,
      source = source,
      target = node,
      
      path_probability = probability,
      
      path_months = months,
      
      path_steps = n_edges,
      
      representative_route = route,
      
      route_type = case_when(
        grepl("Physical", types) &
          grepl("Statistical", types) ~
          "Mixed",
        
        grepl("Physical", types) ~
          "Physical",
        
        grepl("Statistical", types) ~
          "Statistical",
        
        TRUE ~ "Undefined"
      )
    )
}

# ============================================================
# 11. CACHE PATHS BY SOURCE
# ============================================================

sources <- rankings %>%
  distinct(
    nutrient, source
  )

path_list <- vector(
  "list",
  nrow(sources)
)

for (i in seq_len(nrow(sources))) {
  
  path_list[[i]] <- best_paths_for_source(
    nu = sources$nutrient[i],
    
    source = sources$source[i],
    
    edge_data = edge_lags
  )
}

paths <- bind_rows(path_list)

save_csv(
  paths,
  "03_representative_dual_lag_paths.csv"
)


# ============================================================
# 12. SELECT TOP-K TARGETS
# ============================================================

# ============================================================
# PATH QUALITY ASSURANCE
# ============================================================

if (nrow(paths) == 0) {
  stop(
    "No valid dual-lag paths found. ",
    "Check edge_lags and assigned lags."
  )
}

if (anyDuplicated(
  paths[c("nutrient", "source", "target")]
)) {
  stop("Duplicate source-target paths.")
}

if (any(
  !is.finite(paths$path_months) |
  paths$path_months < 0
)) {
  stop("Invalid cumulative path lag.")
}

cat("\nRoute type distribution:\n")

print(
  paths %>%
    count(
      nutrient,
      route_type
    ),
  n = Inf
)



topk <- rankings %>%
  filter(
    rank_FNet <= max(TOP_K)
  ) %>%
  left_join(
    paths,
    by = c(
      "nutrient", "source", "target"
    )
  )

# For a noninteger cumulative time T,
# allow floor(T) and ceiling(T) months.
#
# Examples:
#   T = 0.26 -> months 0 and 1
#   T = 1.00 -> month 1
#   T = 1.60 -> months 1 and 2
#
# This handles monthly timing uncertainty
# without adding a fitted parameter.

topk <- topk %>%
  mutate(
    lag_min = if_else(
      is.finite(path_months),
      as.integer(floor(path_months)),
      NA_integer_
    ),
    
    lag_max = if_else(
      is.finite(path_months),
      as.integer(ceiling(path_months)),
      NA_integer_
    ),
    
    window_defined = (
      is.finite(path_months) &
        lag_min <= MAX_EVAL_LAG &
        lag_max >= 0
    ),
    
    window_label = case_when(
      !is.finite(path_months) ~
        "Unassigned",
      
      lag_min > MAX_EVAL_LAG ~
        ">3 months",
      
      TRUE ~ paste0(
        lag_min, "–", lag_max, " months"
      )
    )
  )

save_csv(
  topk,
  "04_TopK_dual_lag_windows.csv"
)


# ============================================================
# 13. JOIN REAL OBSERVED ANOMALIES
# ============================================================

# Retain observations at ALL available lags,
# then evaluate only the recommended months.

evaluation_months <- topk %>%
  inner_join(
    obs %>%
      select(
        event_id,
        target,
        lag_months,
        target_observed,
        target_anomaly_15,
        target_Z
      ),
    by = c("event_id", "target"),
    relationship = "many-to-many"
  ) %>%
  mutate(
    in_window = (
      window_defined &
        lag_months >= lag_min &
        lag_months <= lag_max
    ),
    
    evaluable_month = (
      in_window &
        target_observed
    ),
    
    anomaly_hit = (
      evaluable_month &
        target_anomaly_15
    )
  )

save_csv(
  evaluation_months,
  "05_TopK_observed_month_comparison.csv"
)


# ============================================================
# 14. TARGET-LEVEL VERIFICATION
# ============================================================

target_verification <- evaluation_months %>%
  group_by(
    event_id,
    nutrient,
    source,
    target,
    rank_FNet,
    source_PCS,
    window_label,
    route_type,
    path_months
  ) %>%
  summarise(
    n_window_months = sum(in_window),
    
    n_observed_months = sum(
      evaluable_month
    ),
    
    observed_anomaly = any(
      anomaly_hit,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    window_fully_observed = (
      n_window_months > 0 &
        n_observed_months == n_window_months
    ),
    
    verification_status = case_when(
      n_window_months == 0 ~
        "Window unavailable",
      
      !window_fully_observed ~
        "Incomplete observation",
      
      observed_anomaly ~
        "Observed anomaly",
      
      TRUE ~ "No observed anomaly"
    )
  )

save_csv(
  target_verification,
  "06_TopK_target_verification.csv"
)


# ============================================================
# 15. EVENT-LEVEL TOP-K EVALUATION
# ============================================================

score_topk <- function(data, k) {
  
  selected <- data %>%
    filter(
      rank_FNet <= k
    )
  
  eligible <- selected %>%
    filter(
      window_fully_observed
    )
  
  tibble(
    K = k,
    
    n_recommended = nrow(selected),
    
    n_evaluable = nrow(eligible),
    
    n_hits = sum(
      eligible$observed_anomaly
    ),
    
    hit_fraction = if (
      nrow(eligible) > 0
    ) {
      mean(eligible$observed_anomaly)
    } else {
      NA_real_
    },
    
    coverage_fraction = if (
      nrow(selected) > 0
    ) {
      nrow(eligible) / nrow(selected)
    } else {
      NA_real_
    }
  )
}

event_results <- list()

ids <- unique(
  target_verification$event_id
)

counter <- 1L

for (id in ids) {
  
  one <- target_verification %>%
    filter(event_id == id)
  
  for (k in TOP_K) {
    
    tmp <- score_topk(one, k)
    
    event_results[[counter]] <- tibble(
      event_id = id,
      
      nutrient = unique(one$nutrient),
      
      source = unique(one$source)
    ) %>%
      bind_cols(tmp)
    
    counter <- counter + 1L
  }
}

event_performance <- bind_rows(
  event_results
)

save_csv(
  event_performance,
  "07_event_TopK_dual_lag_performance.csv"
)


# ============================================================
# 16. OVERALL DESCRIPTIVE PERFORMANCE
# ============================================================

summary_results <- event_performance %>%
  group_by(
    nutrient, K
  ) %>%
  summarise(
    n_events = n(),
    
    n_evaluable_events = sum(
      n_evaluable > 0
    ),
    
    n_total_recommended = sum(
      n_recommended
    ),
    
    n_total_evaluable = sum(
      n_evaluable
    ),
    
    n_total_hits = sum(
      n_hits
    ),
    
    pooled_hit_fraction = if (
      sum(n_evaluable) > 0
    ) {
      sum(n_hits) / sum(n_evaluable)
    } else {
      NA_real_
    },
    
    mean_event_hit_fraction = if (
      any(is.finite(hit_fraction))
    ) {
      mean(
        hit_fraction,
        na.rm = TRUE
      )
    } else {
      NA_real_
    },
    
    coverage_fraction = if (
      sum(n_recommended) > 0
    ) {
      sum(n_evaluable) /
        sum(n_recommended)
    } else {
      NA_real_
    },
    
    .groups = "drop"
  )

save_csv(
  summary_results,
  "08_TopK_dual_lag_summary.csv"
)


# ============================================================
# 17. FIGURE — TOP-K ATTENTION RESULTS
# ============================================================

plot_data <- summary_results %>%
  mutate(
    K = factor(K)
  )

p <- ggplot(
  plot_data,
  
  aes(
    x = K,
    y = mean_event_hit_fraction,
    fill = nutrient
  )
) +
  geom_col(
    width = 0.6
  ) +
  geom_text(
    aes(
      label = scales::percent(
        mean_event_hit_fraction,
        accuracy = 0.1
      )
    ),
    vjust = -0.5,
    size = 3.5
  ) +
  facet_wrap(
    ~nutrient,
    ncol = 2
  ) +
  scale_fill_manual(
    values = c(
      TN = "#287F9C",
      TP = "#CE7953"
    ),
    guide = "none"
  ) +
  scale_y_continuous(
    labels = scales::label_percent(),
    expand = expansion(
      mult = c(0, 0.18)
    )
  ) +
  labs(
    title = "Dual-lag-aware anomaly monitoring attention",
    subtitle = paste0(
      "Source-conditioned F-Net Top-K; ",
      "observed anomalies within assigned windows"
    ),
    x = "Number of prioritized stations",
    y = "Observed anomaly hit fraction"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(
    OUT, "Fig5S_dual_lag_TopK.pdf"
  ),
  p,
  width = 10,
  height = 5.5,
  bg = "white"
)


# ============================================================
# 18. FINAL QA
# ============================================================

cat("\n====================================\n")
cat("STEP 5.3 COMPLETED\n")
cat("====================================\n")

cat("\nDual-lag edge coverage:\n")

print(
  edge_lags %>%
    count(
      nutrient,
      edge_type,
      lag_defined
    ),
  n = Inf
)

cat("\nTop-K window assignment:\n")

print(
  topk %>%
    count(
      nutrient,
      window_label
    ),
  n = Inf
)

cat("\nTop-K verification:\n")

print(
  target_verification %>%
    count(
      nutrient,
      verification_status
    ),
  n = Inf
)

cat("\nPerformance:\n")

print(
  summary_results,
  n = Inf
)

cat("\nOutput:\n", OUT, "\n")

print(p)