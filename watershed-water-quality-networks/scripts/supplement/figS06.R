source("config.R")

# ================================================================
# Figure S6 | Revised monitoring-budget experiment
#
# Training: 2011-10 to 2020-12
# Independent testing: 2021-01 to 2023-12
#
# (a) Annual Recall@10
# (b) Network-method Precision@10
# (c) Fixed-station trigger coverage
# (d) Overall observation-or-flagging coverage
#
# KEY REVISIONS
# 1. One GLOBAL Top-10 recommendation list per month.
# 2. Frozen, within-source percentile ranking of PCS.
# 3. All selected sources may contribute to monthly ranking.
# 4. Station-month anomalies counted only once.
# 5. Precomputed matrices for faster budget replay.
# 6. Report monitoring completeness and direct/extra coverage.
# ================================================================


# ================================================================
# 1. Packages
# ================================================================

packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "ggplot2",
  "patchwork",
  "tibble",
  "scales"
)

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

suppressPackageStartupMessages(
  invisible(
    lapply(
      packages,
      library,
      character.only = TRUE
    )
  )
)


# ================================================================
# 2. Paths
# ================================================================

ROOT <- get_data_root()

TRAIN_ROOT <- file.path(
  ROOT,
  "Fig6_independent_training_2011_2020"
)

PNET_ROOT <- file.path(
  ROOT,
  "PNet_2011_2023_BlockTracing"
)

STEP2_ROOT <- file.path(
  ROOT,
  "Fig6_independent_validation",
  "Step2_event_replay"
)

OLD_OUT <- file.path(
  ROOT,
  "Figure_S6_monitoring_robustness"
)

OUT <- file.path(
  ROOT,
  "Figure_S6_revised_global_top10"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

NTN_FILE <- file.path(
  ROOT,
  "NTN_2023_MonthlyLoadTable.csv"
)


# ================================================================
# 3. Parameters
# ================================================================

NUTRIENTS <- c("TN", "TP")

BUDGETS <- c(
  8L, 12L, 16L, 24L, 28L
)

TOP_K <- 10L

RANDOM_REPS <- 500L

BOOT_REPS <- 1000L

SEED <- 20261101L

TRAIN_START <- as.Date("2011-10-01")
TRAIN_END <- as.Date("2020-12-01")

TEST_START <- as.Date("2021-01-01")
TEST_END <- as.Date("2023-12-01")

Z_THRESHOLD <- 1.5

MIN_TRAIN_Z <- 36L

MIN_MONTH_COVERAGE <- 0.80

MAX_LAG <- 3L

# Use previously computed S6(a-b) results if possible.
REUSE_EVENT_RESULTS <- TRUE

# Common monitoring input population:
# all frozen training-network stations.

# Monthly monitoring selection:
# only observed positive fixed stations can trigger.

# Recommendation:
# one shared global Top 10 per month, NOT per source.

COLORS <- c(
  TN = "#307FA6",
  TP = "#DE8750"
)

STRATEGY_COLORS <- c(
  Random = "#9CA3AF",
  Physical = "#5375A9",
  TE = "#946BA4",
  PCS = "#D78B4D"
)


# ================================================================
# 4. Helper functions
# ================================================================

read_file <- function(
    path,
    ids = character()) {
  
  if (!file.exists(path)) {
    stop("Missing file: ", path)
  }
  
  sp <- readr::cols(
    .default = readr::col_guess()
  )
  
  for (nm in ids) {
    sp$cols[[nm]] <- readr::col_character()
  }
  
  readr::read_csv(
    path,
    col_types = sp,
    show_col_types = FALSE
  )
}


req <- function(
    x,
    required,
    label) {
  
  missing <- setdiff(
    required,
    names(x)
  )
  
  if (length(missing) > 0) {
    
    stop(
      label,
      " missing columns: ",
      paste(missing, collapse = ", "),
      "\nAvailable: ",
      paste(names(x), collapse = ", ")
    )
  }
}


id8 <- function(x) {
  
  x <- trimws(
    sub(
      "\\.0+$",
      "",
      as.character(x)
    )
  )
  
  good <- !is.na(x) &
    grepl("^[0-9]{1,8}$", x)
  
  x[good] <- sprintf(
    "%08.0f",
    as.numeric(x[good])
  )
  
  x
}


to_bool <- function(x) {
  
  if (is.logical(x)) {
    return(x)
  }
  
  z <- tolower(
    trimws(as.character(x))
  )
  
  ifelse(
    z %in% c(
      "true", "t", "1", "yes"
    ),
    TRUE,
    ifelse(
      z %in% c(
        "false", "f", "0", "no"
      ),
      FALSE,
      NA
    )
  )
}


save_csv <- function(
    x,
    filename) {
  
  readr::write_csv(
    x,
    file.path(OUT, filename),
    na = ""
  )
}


tie_value <- function(
    source,
    target) {
  
  x <- paste(
    source,
    target,
    sep = "-"
  )
  
  vapply(
    x,
    function(z) {
      
      h <- 7
      
      for (a in utf8ToInt(z)) {
        h <- (h * 137 + a) %% 1000000007
      }
      
      h / 1000000007
    },
    numeric(1)
  )
}


bootstrap_cluster <- function(
    values,
    source,
    B = BOOT_REPS) {
  
  valid <- is.finite(values) &
    !is.na(source)
  
  values <- values[valid]
  source <- source[valid]
  
  if (
    length(unique(source)) < 2
  ) {
    
    return(c(
      low = NA_real_,
      high = NA_real_
    ))
  }
  
  groups <- split(
    values,
    source
  )
  
  ids <- names(groups)
  
  simulations <- replicate(
    B,
    {
      
      sample_ids <- sample(
        ids,
        length(ids),
        replace = TRUE
      )
      
      mean(
        unlist(
          groups[sample_ids],
          use.names = FALSE
        )
      )
    }
  )
  
  unname(
    quantile(
      simulations,
      c(0.025, 0.975)
    )
  )
}


# ================================================================
# 5. Read frozen training networks
# ================================================================

read_training <- function(nu) {
  
  cat("\nReading:", nu, "\n")
  
  folder <- file.path(
    TRAIN_ROOT,
    nu
  )
  
  F_RAW <- read_file(
    file.path(
      folder,
      "02_training_FNet_all_pairs.csv"
    ),
    c("from", "to")
  )
  
  TE_RAW <- read_file(
    file.path(
      folder,
      "train_2011_2020_all_edges.csv"
    ),
    c("from", "to")
  )
  
  SOURCE_RAW <- read_file(
    file.path(
      folder,
      "04_training_source_conditioned_PCS.csv"
    ),
    c("source", "target")
  )
  
  GLOBAL_RAW <- read_file(
    file.path(
      folder,
      "05_training_global_PCS.csv"
    ),
    "STAID"
  )
  
  req(
    F_RAW,
    c("from", "to", "W_fus"),
    "Training F-Net"
  )
  
  req(
    SOURCE_RAW,
    c("source", "target", "source_PCS"),
    "Source-conditioned PCS"
  )
  
  req(
    GLOBAL_RAW,
    c("STAID", "PCS_in"),
    "Global PCS"
  )
  
  F <- F_RAW %>%
    transmute(
      source = id8(from),
      target = id8(to),
      W = as.numeric(W_fus)
    )
  
  SOURCE <- SOURCE_RAW %>%
    transmute(
      source = id8(source),
      target = id8(target),
      source_PCS = as.numeric(source_PCS)
    )
  
  GLOBAL <- GLOBAL_RAW %>%
    transmute(
      STAID = id8(STAID),
      PCS = as.numeric(PCS_in)
    )
  
  ids <- sort(
    unique(
      c(F$source, F$target)
    )
  )
  
  ids <- ids[
    !is.na(ids)
  ]
  
  if (
    anyDuplicated(
      SOURCE[c("source", "target")]
    ) ||
    anyDuplicated(GLOBAL["STAID"])
  ) {
    
    stop(nu, ": duplicate PCS keys")
  }
  
  if (
    any(!is.finite(SOURCE$source_PCS)) ||
    any(SOURCE$source_PCS < 0)
  ) {
    
    stop(nu, ": invalid source PCS")
  }
  
  
  # ------------------------------------------------
  # 5.1 TE edge identification
  # ------------------------------------------------
  
  filtered_cols <- c(
    "TE_weight",
    "te_weight",
    "TE_sig_weight",
    "TE_significant_weight"
  )
  
  raw_cols <- c(
    "TE_raw",
    "te_raw",
    "TE_value"
  )
  
  p_cols <- c(
    "p_value",
    "TE_p",
    "te_p",
    "p",
    "pval"
  )
  
  sig_cols <- c(
    "TE_significant",
    "te_significant",
    "significant"
  )
  
  te_sources <- list(
    independent = TE_RAW,
    training_fusion = F_RAW
  )
  
  selected_te <- NULL
  
  # Prefer explicitly filtered TE weights.
  
  for (nm in names(te_sources)) {
    
    x <- te_sources[[nm]]
    
    col <- intersect(
      filtered_cols,
      names(x)
    )
    
    if (length(col) > 0) {
      
      selected_te <- list(
        data = x,
        field = col[1],
        mode = "filtered",
        source = nm
      )
      
      break
    }
  }
  
  # Otherwise use raw TE with significance.
  
  if (is.null(selected_te)) {
    
    for (nm in names(te_sources)) {
      
      x <- te_sources[[nm]]
      
      col <- intersect(
        raw_cols,
        names(x)
      )
      
      has_p <- length(
        intersect(
          p_cols,
          names(x)
        )
      ) > 0
      
      has_sig <- length(
        intersect(
          sig_cols,
          names(x)
        )
      ) > 0
      
      if (
        length(col) > 0 &&
        (has_p || has_sig)
      ) {
        
        selected_te <- list(
          data = x,
          field = col[1],
          mode = "raw",
          source = nm
        )
        
        break
      }
    }
  }
  
  if (is.null(selected_te)) {
    
    stop(
      nu,
      ": no valid training TE layer"
    )
  }
  
  x <- selected_te$data
  
  te_raw <- suppressWarnings(
    as.numeric(
      x[[selected_te$field]]
    )
  )
  
  te_score <- rep(
    0,
    length(te_raw)
  )
  
  if (
    selected_te$mode == "filtered"
  ) {
    
    if (
      any(
        is.finite(te_raw) &
        te_raw < 0
      )
    ) {
      stop(nu, ": negative TE weights")
    }
    
    te_score[
      is.finite(te_raw)
    ] <- te_raw[
      is.finite(te_raw)
    ]
    
  } else {
    
    pfield <- intersect(
      p_cols,
      names(x)
    )
    
    sfield <- intersect(
      sig_cols,
      names(x)
    )
    
    p <- rep(
      NA_real_,
      length(te_raw)
    )
    
    sig <- rep(
      NA,
      length(te_raw)
    )
    
    if (length(pfield) > 0) {
      
      p <- suppressWarnings(
        as.numeric(
          x[[pfield[1]]]
        )
      )
    }
    
    if (length(sfield) > 0) {
      
      sig <- to_bool(
        x[[sfield[1]]]
      )
    }
    
    known_p <- is.finite(p) &
      p >= 0 &
      p <= 1
    
    known_sig <- !known_p &
      !is.na(sig)
    
    positive <- is.finite(te_raw) &
      te_raw > 0
    
    unresolved <- positive &
      !known_p &
      !known_sig
    
    if (any(unresolved)) {
      
      save_csv(
        tibble(
          nutrient = nu,
          from = id8(
            x$from[unresolved]
          ),
          to = id8(
            x$to[unresolved]
          ),
          TE_raw = te_raw[unresolved],
          p_value = p[unresolved]
        ),
        paste0(
          "S6_",
          nu,
          "_unresolved_TE.csv"
        )
      )
      
      stop(
        nu,
        ": positive TE without significance"
      )
    }
    
    significant <- rep(
      FALSE,
      length(te_raw)
    )
    
    significant[known_p] <-
      p[known_p] < 0.10
    
    significant[known_sig] <-
      sig[known_sig]
    
    keep <- significant & positive
    
    te_score[keep] <- te_raw[keep]
  }
  
  TE <- tibble(
    source = id8(x$from),
    target = id8(x$to),
    TE = te_score
  ) %>%
    filter(
      !is.na(source),
      !is.na(target),
      source != target
    )
  
  if (
    anyDuplicated(
      TE[c("source", "target")]
    )
  ) {
    
    stop(nu, ": duplicate TE pairs")
  }
  
  cat(
    "TE source:",
    selected_te$source,
    "\nTE field:",
    selected_te$field,
    "\nPositive TE edges:",
    sum(TE$TE > 0),
    "\n"
  )
  
  
  # ------------------------------------------------
  # 5.2 Physical network
  # ------------------------------------------------
  
  phys_path <- file.path(
    PNET_ROOT,
    nu,
    "PNet_directed_reachable_pairs.csv"
  )
  
  PH_RAW <- read_file(
    phys_path,
    c("from", "to")
  )
  
  req(
    PH_RAW,
    c("from", "to"),
    phys_path
  )
  
  PH <- PH_RAW %>%
    transmute(
      source = id8(from),
      target = id8(to)
    ) %>%
    filter(
      !is.na(source),
      !is.na(target),
      source != target
    ) %>%
    distinct()
  
  # Missing PH pairs = unknown physical relation,
  # not explicitly unreachable.
  
  
  # ------------------------------------------------
  # 5.3 Training selection scores
  # ------------------------------------------------
  
  phys_out <- PH %>%
    filter(
      source %in% ids,
      target %in% ids
    ) %>%
    count(
      source,
      name = "Physical"
    )
  
  te_out <- TE %>%
    group_by(source) %>%
    summarise(
      TE = sum(TE),
      .groups = "drop"
    )
  
  SCORES <- tibble(
    STAID = ids
  ) %>%
    left_join(
      phys_out,
      by = c("STAID" = "source")
    ) %>%
    left_join(
      te_out,
      by = c("STAID" = "source")
    ) %>%
    left_join(
      GLOBAL,
      by = "STAID"
    )
  
  for (
    nm in c("Physical", "TE", "PCS")
  ) {
    
    SCORES[[nm]][
      is.na(SCORES[[nm]])
    ] <- 0
  }
  
  if (
    any(
      !is.finite(
        as.matrix(
          SCORES[
            c(
              "Physical",
              "TE",
              "PCS"
            )
          ]
        )
      )
    )
  ) {
    
    stop(nu, ": invalid selection score")
  }
  
  
  list(
    ids = ids,
    F = F,
    SOURCE = SOURCE,
    GLOBAL = GLOBAL,
    TE = TE,
    PH = PH,
    SCORES = SCORES
  )
}


TR <- setNames(
  lapply(
    NUTRIENTS,
    read_training
  ),
  NUTRIENTS
)

save_csv(
  bind_rows(
    lapply(
      NUTRIENTS,
      function(nu) {
        
        TR[[nu]]$SCORES %>%
          mutate(
            nutrient = nu
          )
      }
    )
  ),
  "S6_00_training_selection_scores.csv"
)


# ================================================================
# 6. S6(a-b): Read or compute event metrics
# ================================================================

event_cache <- file.path(
  OLD_OUT,
  "S6_02_event_method_expected_metrics.csv"
)

if (
  REUSE_EVENT_RESULTS &&
  file.exists(event_cache)
) {
  
  cat(
    "\nReusing previously calculated ",
    "event metrics.\n"
  )
  
  EVENT_METRIC <- read_file(
    event_cache,
    "source"
  )
  
  req(
    EVENT_METRIC,
    c(
      "nutrient",
      "event_id",
      "source",
      "date",
      "method",
      "precision",
      "recall"
    ),
    event_cache
  )
  
  EVENT_METRIC <- EVENT_METRIC %>%
    mutate(
      source = id8(source),
      date = as.Date(date)
    )
  
} else {
  
  cat(
    "\nRecalculating S6(a-b)...\n"
  )
  
  EV_RAW <- read_file(
    file.path(
      STEP2_ROOT,
      "03_independent_test_events.csv"
    ),
    "source"
  )
  
  LABEL_RAW <- read_file(
    file.path(
      STEP2_ROOT,
      "05_common_event_target_evaluation.csv"
    ),
    c("source", "target")
  )
  
  req(
    EV_RAW,
    c(
      "event_id",
      "nutrient",
      "source",
      "trigger_date"
    ),
    "Fig6 events"
  )
  
  req(
    LABEL_RAW,
    c(
      "event_id",
      "nutrient",
      "source",
      "target",
      "evaluable",
      "hit_03"
    ),
    "Fig6 labels"
  )
  
  EV <- EV_RAW %>%
    transmute(
      event_id = as.character(event_id),
      nutrient = toupper(
        as.character(nutrient)
      ),
      source = id8(source),
      date = as.Date(trigger_date)
    )
  
  U <- LABEL_RAW %>%
    transmute(
      event_id = as.character(event_id),
      nutrient = toupper(
        as.character(nutrient)
      ),
      source = id8(source),
      target = id8(target),
      evaluable = to_bool(evaluable),
      hit = to_bool(hit_03)
    ) %>%
    inner_join(
      EV,
      by = c(
        "event_id",
        "nutrient",
        "source"
      )
    ) %>%
    filter(
      evaluable %in% TRUE,
      !is.na(hit),
      source != target
    )
  
  U <- bind_rows(
    lapply(
      NUTRIENTS,
      function(nu) {
        
        U %>%
          filter(
            nutrient == nu,
            source %in% TR[[nu]]$ids,
            target %in% TR[[nu]]$ids
          )
      }
    )
  ) %>%
    group_by(event_id) %>%
    filter(n() >= TOP_K) %>%
    ungroup()
  
  if (
    anyDuplicated(
      U[c("event_id", "target")]
    )
  ) {
    
    stop("Duplicated Fig6 event-target pairs")
  }
  
  if (nrow(U) == 0) {
    stop("No evaluable Fig6 events")
  }
  
  
  # ----------------------------------------------
  # Event-level four-method comparison
  # ----------------------------------------------
  
  rank_once <- function(
    dat,
    nu,
    method) {
    
    tr <- TR[[nu]]
    
    if (method == "F-Net") {
      
      SC <- tr$SOURCE %>%
        transmute(
          source,
          target,
          score = source_PCS
        )
    }
    
    if (method == "TE-Net") {
      
      SC <- tr$TE %>%
        transmute(
          source,
          target,
          score = TE
        )
    }
    
    if (method == "Phys-Net") {
      
      SC <- tr$PH %>%
        mutate(score = 1)
    }
    
    if (method == "Random") {
      
      SC <- tibble(
        source = character(),
        target = character(),
        score = numeric()
      )
    }
    
    dat %>%
      left_join(
        SC,
        by = c("source", "target")
      ) %>%
      mutate(
        score = ifelse(
          is.finite(score),
          score,
          0
        ),
        tie = runif(n())
      ) %>%
      group_by(event_id) %>%
      arrange(
        desc(score),
        tie,
        .by_group = TRUE
      ) %>%
      mutate(rank = row_number()) %>%
      summarise(
        source = first(source),
        date = first(date),
        precision = sum(
          hit[rank <= TOP_K]
        ) / TOP_K,
        recall = if (
          sum(hit) > 0
        ) {
          sum(
            hit[rank <= TOP_K]
          ) / sum(hit)
        } else {
          NA_real_
        },
        .groups = "drop"
      ) %>%
      mutate(
        nutrient = nu,
        method = method
      )
  }
  
  set.seed(SEED)
  
  EVENT_REPS <- bind_rows(
    lapply(
      NUTRIENTS,
      function(nu) {
        
        dat <- U %>%
          filter(nutrient == nu)
        
        bind_rows(
          lapply(
            c(
              "Random",
              "Phys-Net",
              "TE-Net",
              "F-Net"
            ),
            function(method) {
              
              bind_rows(
                lapply(
                  seq_len(RANDOM_REPS),
                  function(r) {
                    
                    rank_once(
                      dat,
                      nu,
                      method
                    ) %>%
                      mutate(rep = r)
                  }
                )
              )
            }
          )
        )
      }
    )
  )
  
  save_csv(
    EVENT_REPS,
    "S6_event_random_repetitions.csv"
  )
  
  EVENT_METRIC <- EVENT_REPS %>%
    group_by(
      nutrient,
      event_id,
      source,
      date,
      method
    ) %>%
    summarise(
      precision = mean(precision),
      recall = if (
        all(is.na(recall))
      ) {
        NA_real_
      } else {
        mean(
          recall,
          na.rm = TRUE
        )
      },
      .groups = "drop"
    )
}

save_csv(
  EVENT_METRIC,
  "S6_02_event_method_expected_metrics.csv"
)


# ================================================================
# 7. S6(a): Annual Recall@10
# ================================================================

set.seed(SEED + 1L)

ANNUAL <- EVENT_METRIC %>%
  filter(
    method == "F-Net"
  ) %>%
  mutate(
    year = as.integer(
      format(date, "%Y")
    )
  ) %>%
  group_by(
    nutrient,
    year
  ) %>%
  group_modify(
    ~ {
      
      x <- .x[
        is.finite(.x$recall),
        ,
        drop = FALSE
      ]
      
      ci <- bootstrap_cluster(
        x$recall,
        x$source
      )
      
      tibble(
        n_events = nrow(.x),
        n_defined = nrow(x),
        mean_recall = if (
          nrow(x) > 0
        ) {
          mean(x$recall)
        } else {
          NA_real_
        },
        low = ci[1],
        high = ci[2]
      )
    }
  ) %>%
  ungroup()

save_csv(
  ANNUAL,
  "S6_05_annual_recall.csv"
)


# ================================================================
# 8. S6(b): Method Precision@10
# ================================================================

set.seed(SEED + 2L)

METHOD <- EVENT_METRIC %>%
  group_by(
    nutrient,
    method
  ) %>%
  group_modify(
    ~ {
      
      ci <- bootstrap_cluster(
        .x$precision,
        .x$source
      )
      
      tibble(
        n_events = nrow(.x),
        mean_precision = mean(
          .x$precision
        ),
        low = ci[1],
        high = ci[2]
      )
    }
  ) %>%
  ungroup()

save_csv(
  METHOD,
  "S6_06_method_precision.csv"
)


# ================================================================
# 9. Read NTN monthly concentration
# ================================================================

NTN_RAW <- read_file(
  NTN_FILE,
  c("STAID", "PCODE")
)

req(
  NTN_RAW,
  c(
    "STAID",
    "PCODE",
    "Year",
    "Month",
    "KConc"
  ),
  NTN_FILE
)

NTN <- NTN_RAW %>%
  transmute(
    nutrient = case_when(
      PCODE == "P00600" ~ "TN",
      PCODE == "P00665" ~ "TP",
      TRUE ~ NA_character_
    ),
    STAID = id8(STAID),
    year = as.integer(Year),
    month = as.integer(Month),
    concentration = as.numeric(KConc)
  ) %>%
  filter(
    !is.na(nutrient),
    year >= 2011,
    year <= 2023,
    month %in% 1:12
  ) %>%
  mutate(
    date = as.Date(
      sprintf(
        "%04d-%02d-01",
        year,
        month
      )
    )
  )

if (
  anyDuplicated(
    NTN[c(
      "nutrient",
      "STAID",
      "date"
    )]
  )
) {
  
  stop(
    "Duplicated NTN station-month data"
  )
}


# ================================================================
# 10. Frozen training climatology
# ================================================================

CLIM <- NTN %>%
  filter(
    date >= TRAIN_START,
    date <= TRAIN_END,
    is.finite(concentration)
  ) %>%
  group_by(
    nutrient,
    STAID,
    month
  ) %>%
  summarise(
    mu = mean(concentration),
    sigma = sd(concentration),
    .groups = "drop"
  )

NTN <- NTN %>%
  left_join(
    CLIM,
    by = c(
      "nutrient",
      "STAID",
      "month"
    )
  ) %>%
  mutate(
    Z = if_else(
      is.finite(concentration) &
        is.finite(sigma) &
        sigma > 0,
      (concentration - mu) / sigma,
      NA_real_
    )
  )

TRAIN_COVERAGE <- NTN %>%
  filter(
    date >= TRAIN_START,
    date <= TRAIN_END
  ) %>%
  group_by(
    nutrient,
    STAID
  ) %>%
  summarise(
    n_train_Z = sum(is.finite(Z)),
    .groups = "drop"
  )


# ================================================================
# 11. Independent monthly anomaly matrix
# ================================================================

TEST_ANOM <- bind_rows(
  lapply(
    NUTRIENTS,
    function(nu) {
      
      tidyr::expand_grid(
        nutrient = nu,
        STAID = TR[[nu]]$ids,
        date = seq(
          TEST_START,
          TEST_END,
          by = "month"
        )
      ) %>%
        left_join(
          NTN %>%
            select(
              nutrient,
              STAID,
              date,
              Z
            ),
          by = c(
            "nutrient",
            "STAID",
            "date"
          )
        ) %>%
        left_join(
          TRAIN_COVERAGE,
          by = c(
            "nutrient",
            "STAID"
          )
        ) %>%
        mutate(
          observed = is.finite(Z) &
            is.finite(n_train_Z) &
            n_train_Z >= MIN_TRAIN_Z,
          anomaly = if_else(
            observed,
            Z > Z_THRESHOLD,
            NA
          )
        )
    }
  )
)

save_csv(
  TEST_ANOM,
  "S6_01_test_anomaly_matrix.csv"
)


# ================================================================
# 12. Month-level observation QA
# ================================================================

MONTH_QA <- TEST_ANOM %>%
  group_by(
    nutrient,
    date
  ) %>%
  summarise(
    total_stations = n(),
    observed_stations = sum(observed),
    observation_fraction =
      observed_stations / total_stations,
    positive_stations = sum(
      anomaly %in% TRUE
    ),
    eligible =
      observation_fraction >=
      MIN_MONTH_COVERAGE,
    .groups = "drop"
  )

save_csv(
  MONTH_QA,
  "S6_01_month_coverage_QA.csv"
)


# ================================================================
# 13. Convert to fast matrices
# ================================================================

make_matrices <- function(nu) {
  
  ids <- TR[[nu]]$ids
  
  dates <- seq(
    TEST_START,
    TEST_END,
    by = "month"
  )
  
  n <- length(ids)
  nt <- length(dates)
  
  observed <- matrix(
    FALSE,
    nrow = n,
    ncol = nt,
    dimnames = list(
      ids,
      as.character(dates)
    )
  )
  
  anomaly <- observed
  
  dat <- TEST_ANOM %>%
    filter(nutrient == nu)
  
  ri <- match(
    dat$STAID,
    ids
  )
  
  ci <- match(
    dat$date,
    dates
  )
  
  valid <- !is.na(ri) &
    !is.na(ci)
  
  observed[
    cbind(
      ri[valid],
      ci[valid]
    )
  ] <- dat$observed[valid]
  
  anomaly[
    cbind(
      ri[valid],
      ci[valid]
    )
  ] <- dat$anomaly[valid] %in% TRUE
  
  eligible_month <- MONTH_QA %>%
    filter(
      nutrient == nu
    ) %>%
    arrange(date) %>%
    pull(eligible)
  
  if (
    length(eligible_month) != nt
  ) {
    
    stop(nu, ": month QA length mismatch")
  }
  
  list(
    ids = ids,
    dates = dates,
    observed = observed,
    anomaly = anomaly,
    eligible = eligible_month
  )
}

MATS <- setNames(
  lapply(
    NUTRIENTS,
    make_matrices
  ),
  NUTRIENTS
)


# ================================================================
# 14. Frozen source-conditioned PCS score matrix
#
# Rows: triggering source stations
# Cols: potential hidden target stations
#
# Within-source positive-score percentile ranking:
# strongest target -> 1
# zero-score or unavailable -> 0
#
# This makes source scores comparable for aggregation.
# ================================================================

make_pcs_matrix <- function(nu) {
  
  ids <- TR[[nu]]$ids
  n <- length(ids)
  
  raw <- matrix(
    0,
    nrow = n,
    ncol = n,
    dimnames = list(ids, ids)
  )
  
  s <- TR[[nu]]$SOURCE
  
  ri <- match(
    s$source,
    ids
  )
  
  ci <- match(
    s$target,
    ids
  )
  
  valid <- !is.na(ri) &
    !is.na(ci) &
    is.finite(s$source_PCS)
  
  raw[
    cbind(
      ri[valid],
      ci[valid]
    )
  ] <- s$source_PCS[valid]
  
  diag(raw) <- 0
  
  out <- matrix(
    0,
    nrow = n,
    ncol = n,
    dimnames = list(ids, ids)
  )
  
  for (i in seq_len(n)) {
    
    positive <- which(
      raw[i, ] > 0
    )
    
    if (
      length(positive) == 0
    ) {
      next
    }
    
    v <- raw[i, positive]
    
    # Fixed training-only percentile score.
    # Average ranks are used for ties.
    
    out[i, positive] <- rank(
      v,
      ties.method = "average"
    ) / length(v)
  }
  
  diag(out) <- 0
  
  out
}

PCS_MATS <- setNames(
  lapply(
    NUTRIENTS,
    make_pcs_matrix
  ),
  NUTRIENTS
)


# ================================================================
# 15. Fixed monitoring station sets
# ================================================================

set.seed(SEED + 3L)

make_station_sets <- function(nu) {
  
  scores <- TR[[nu]]$SCORES
  
  ids <- TR[[nu]]$ids
  
  n <- length(ids)
  
  if (
    max(BUDGETS) > n
  ) {
    
    stop(nu, ": budget exceeds network nodes")
  }
  
  out <- list()
  idx <- 1L
  
  for (
    method in c(
      "Physical",
      "TE",
      "PCS"
    )
  ) {
    
    ranking <- scores %>%
      mutate(
        tie = tie_value(
          paste0(nu, "_selection"),
          STAID
        )
      ) %>%
      arrange(
        desc(.data[[method]]),
        tie
      ) %>%
      pull(STAID)
    
    for (b in BUDGETS) {
      
      out[[idx]] <- tibble(
        nutrient = nu,
        method = method,
        budget = b,
        rep = 1L,
        STAID = ranking[
          seq_len(b)
        ]
      )
      
      idx <- idx + 1L
    }
  }
  
  for (
    r in seq_len(RANDOM_REPS)
  ) {
    
    ranking <- sample(
      ids,
      length(ids),
      replace = FALSE
    )
    
    for (b in BUDGETS) {
      
      out[[idx]] <- tibble(
        nutrient = nu,
        method = "Random",
        budget = b,
        rep = r,
        STAID = ranking[
          seq_len(b)
        ]
      )
      
      idx <- idx + 1L
    }
  }
  
  bind_rows(out)
}

SETS <- bind_rows(
  lapply(
    NUTRIENTS,
    make_station_sets
  )
)

save_csv(
  SETS,
  "S6_02_fixed_station_sets.csv"
)


# ================================================================
# 16. Fast budget replay
#
# At each month:
# - use only positive OBSERVED fixed stations as sources
# - aggregate their frozen PCS rankings
# - choose <= Top 10 distinct HIDDEN targets
# - maintain recommendation for months 0-3
# - do not allow hidden outcomes to trigger
#
# Candidate ranking is independent of test outcomes.
# ================================================================

replay_fast <- function(
    nu,
    selected_ids,
    budget,
    method,
    rep_id) {
  
  dat <- MATS[[nu]]
  pcs <- PCS_MATS[[nu]]
  
  ids <- dat$ids
  
  n <- length(ids)
  nt <- length(dat$dates)
  
  selected <- match(
    selected_ids,
    ids
  )
  
  if (
    anyNA(selected)
  ) {
    
    stop(
      nu,
      ": fixed-station ID not found"
    )
  }
  
  hidden <- setdiff(
    seq_len(n),
    selected
  )
  
  obs <- dat$observed
  A <- dat$anomaly
  
  eligible <- dat$eligible
  
  # Exclude non-eligible evaluation months.
  # Keep an intact monthly timeline so
  # recommendation windows retain calendar lags.
  
  A_eval <- A
  
  A_eval[
    ,
    !eligible
  ] <- FALSE
  
  total_positive <- sum(A_eval)
  
  if (
    total_positive == 0
  ) {
    
    return(tibble())
  }
  
  anomalous_months <- colSums(
    A_eval
  ) > 0
  
  n_anomalous_months <- sum(
    anomalous_months
  )
  
  # Direct fixed-station observations.
  
  fixed_positive <- A[
    selected,
    ,
    drop = FALSE
  ]
  
  fixed_positive[
    ,
    !eligible
  ] <- FALSE
  
  direct_flag <- matrix(
    FALSE,
    nrow = n,
    ncol = nt
  )
  
  direct_flag[
    selected,
    
  ] <- fixed_positive
  
  direct_count <- sum(
    direct_flag & A_eval
  )
  
  triggered <- colSums(
    fixed_positive
  ) > 0
  
  trigger_count <- sum(
    triggered & anomalous_months
  )
  
  trigger_coverage <-
    trigger_count / n_anomalous_months
  
  
  # ----------------------------------------------
  # 16.1 Monitoring completeness
  # ----------------------------------------------
  
  monitoring_completeness <- mean(
    obs[
      selected,
      eligible,
      drop = FALSE
    ]
  )
  
  
  # ----------------------------------------------
  # 16.2 Monthly GLOBAL Top-10 recommendations
  # ----------------------------------------------
  
  # Monthly recommendation matrix:
  # columns represent trigger months.
  
  recommended <- matrix(
    FALSE,
    nrow = n,
    ncol = nt
  )
  
  # Training-only deterministic tie priority.
  tie <- tie_value(
    paste0(nu, "_global_recommendation"),
    ids
  )
  
  recommended_counts <- integer(nt)
  active_source_counts <- integer(nt)
  zero_score_counts <- integer(nt)
  
  for (j in which(triggered)) {
    
    sources <- selected[
      fixed_positive[, j]
    ]
    
    active_source_counts[j] <-
      length(sources)
    
    if (
      length(sources) == 0
    ) {
      next
    }
    
    # Aggregation: max normalized PCS rank
    # over all currently abnormal fixed sources.
    
    combined <- apply(
      pcs[
        sources,
        ,
        drop = FALSE
      ],
      2,
      max
    )
    
    candidate_score <- combined[hidden]
    
    candidate_tie <- tie[hidden]
    
    # Rank ALL hidden candidates.
    # Zero-score fallback uses fixed order,
    # not observations from 2021-2023.
    
    ord <- order(
      -candidate_score,
      candidate_tie,
      ids[hidden]
    )
    
    k <- min(
      TOP_K,
      length(ord)
    )
    
    top_targets <- hidden[
      ord[
        seq_len(k)
      ]
    ]
    
    recommended[
      top_targets,
      j
    ] <- TRUE
    
    recommended_counts[j] <-
      length(top_targets)
    
    zero_score_counts[j] <- sum(
      combined[top_targets] <= 0
    )
  }
  
  
  # ----------------------------------------------
  # 16.3 Recommendation windows
  # ----------------------------------------------
  
  # 0-3 months and 1-3 months,
  # calculated from SAME monthly ranking.
  
  make_active <- function(lag_start) {
    
    active <- matrix(
      FALSE,
      nrow = n,
      ncol = nt
    )
    
    for (
      lag in lag_start:MAX_LAG
    ) {
      
      if (
        lag >= nt
      ) {
        next
      }
      
      source_columns <- seq_len(
        nt - lag
      )
      
      target_columns <- source_columns + lag
      
      active[
        ,
        target_columns
      ] <- active[
        ,
        target_columns,
        drop = FALSE
      ] | recommended[
        ,
        source_columns,
        drop = FALSE
      ]
    }
    
    active
  }
  
  active_03 <- make_active(0L)
  active_13 <- make_active(1L)
  
  
  # ----------------------------------------------
  # 16.4 Hidden anomaly coverage
  # ----------------------------------------------
  
  hidden_positive <- A_eval
  
  hidden_positive[
    selected,
    
  ] <- FALSE
  
  additional_03 <- sum(
    hidden_positive & active_03
  )
  
  additional_13 <- sum(
    hidden_positive & active_13
  )
  
  # QA: no overlap between direct/hidden anomalies.
  
  if (
    any(
      direct_flag & hidden_positive
    )
  ) {
    
    stop(
      nu,
      ": direct and hidden anomaly overlap"
    )
  }
  
  if (
    direct_count + additional_03 >
    total_positive
  ) {
    
    stop(
      nu,
      ": coverage exceeds 100%"
    )
  }
  
  
  # ----------------------------------------------
  # 16.5 Outputs
  # ----------------------------------------------
  
  tibble(
    nutrient = nu,
    method = method,
    budget = budget,
    rep = rep_id,
    
    n_eligible_months = sum(eligible),
    n_anomaly_months = n_anomalous_months,
    
    total_anomaly_station_months =
      total_positive,
    
    triggered_months =
      trigger_count,
    
    trigger_coverage =
      trigger_coverage,
    
    monitoring_completeness =
      monitoring_completeness,
    
    direct_count =
      direct_count,
    
    direct_coverage =
      direct_count / total_positive,
    
    flagged_03 =
      additional_03,
    
    flagged_13 =
      additional_13,
    
    additional_coverage_03 =
      additional_03 / total_positive,
    
    additional_coverage_13 =
      additional_13 / total_positive,
    
    overall_coverage_03 =
      (
        direct_count + additional_03
      ) / total_positive,
    
    overall_coverage_13 =
      (
        direct_count + additional_13
      ) / total_positive,
    
    n_alarm_months =
      sum(triggered),
    
    n_recommendations =
      sum(recommended_counts),
    
    n_sources_activated =
      sum(active_source_counts),
    
    n_zero_score_recommendations =
      sum(zero_score_counts),
    
    zero_score_fraction =
      if (
        sum(recommended_counts) > 0
      ) {
        sum(zero_score_counts) /
          sum(recommended_counts)
      } else {
        NA_real_
      },
    
    max_recommendations_per_month =
      max(recommended_counts)
  )
}


# ================================================================
# 17. Run all budget simulations
# ================================================================

cat(
  "\nStarting revised S6(c-d) replay...\n"
)

SET_GROUPS <- split(
  SETS,
  interaction(
    SETS$nutrient,
    SETS$method,
    SETS$budget,
    SETS$rep,
    drop = TRUE
  )
)

n_groups <- length(
  SET_GROUPS
)

RESULTS <- vector(
  "list",
  n_groups
)

for (i in seq_along(SET_GROUPS)) {
  
  x <- SET_GROUPS[[i]]
  
  RESULTS[[i]] <- replay_fast(
    nu = x$nutrient[1],
    selected_ids = x$STAID,
    budget = x$budget[1],
    method = x$method[1],
    rep_id = x$rep[1]
  )
  
  if (
    i %% 100 == 0 ||
    i == n_groups
  ) {
    
    cat(
      sprintf(
        "Completed %d / %d scenarios (%.1f%%)\n",
        i,
        n_groups,
        100 * i / n_groups
      )
    )
  }
}

BUDGET_RESULT <- bind_rows(
  RESULTS
)

if (
  nrow(BUDGET_RESULT) == 0
) {
  stop("No budget results")
}

save_csv(
  BUDGET_RESULT,
  "S6_03_04_budget_replay_replicates.csv"
)


# ================================================================
# 18. Budget summaries
# ================================================================

BUDGET_SUMMARY <- BUDGET_RESULT %>%
  group_by(
    nutrient,
    method,
    budget
  ) %>%
  summarise(
    n_replicates = n(),
    
    trigger_mean = mean(
      trigger_coverage
    ),
    
    trigger_low = quantile(
      trigger_coverage,
      0.025
    ),
    
    trigger_high = quantile(
      trigger_coverage,
      0.975
    ),
    
    direct_mean = mean(
      direct_coverage
    ),
    
    additional_mean_03 = mean(
      additional_coverage_03
    ),
    
    additional_mean_13 = mean(
      additional_coverage_13
    ),
    
    overall_mean_03 = mean(
      overall_coverage_03
    ),
    
    overall_low_03 = quantile(
      overall_coverage_03,
      0.025
    ),
    
    overall_high_03 = quantile(
      overall_coverage_03,
      0.975
    ),
    
    overall_mean_13 = mean(
      overall_coverage_13
    ),
    
    monitoring_completeness = mean(
      monitoring_completeness
    ),
    
    zero_score_fraction = mean(
      zero_score_fraction,
      na.rm = TRUE
    ),
    
    denominator = first(
      total_anomaly_station_months
    ),
    
    .groups = "drop"
  )

save_csv(
  BUDGET_SUMMARY,
  "S6_03_04_budget_summary.csv"
)


# ================================================================
# 19. QA: shared denominator and Top-K limit
# ================================================================

if (
  any(
    BUDGET_RESULT$
    max_recommendations_per_month >
    TOP_K
  )
) {
  
  stop(
    "Global Top-K budget was violated"
  )
}

DENOMINATOR_QA <- BUDGET_RESULT %>%
  group_by(nutrient) %>%
  summarise(
    n_denominators = n_distinct(
      total_anomaly_station_months
    ),
    denominator = first(
      total_anomaly_station_months
    ),
    n_month_counts = n_distinct(
      n_eligible_months
    ),
    .groups = "drop"
  )

if (
  any(
    DENOMINATOR_QA$n_denominators != 1
  ) ||
  any(
    DENOMINATOR_QA$n_month_counts != 1
  )
) {
  
  stop(
    "Evaluation denominators are inconsistent"
  )
}

save_csv(
  DENOMINATOR_QA,
  "S6_denominator_QA.csv"
)

COMPLETENESS_QA <- BUDGET_SUMMARY %>%
  select(
    nutrient,
    method,
    budget,
    monitoring_completeness,
    trigger_mean,
    direct_mean,
    additional_mean_03,
    overall_mean_03,
    zero_score_fraction
  )

save_csv(
  COMPLETENESS_QA,
  "S6_monitoring_completeness_QA.csv"
)


# ================================================================
# 20. Plot theme
# ================================================================

theme_s6 <- theme_bw(
  base_size = 11
) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(
      face = "bold",
      size = 12
    ),
    strip.text = element_text(
      face = "bold"
    ),
    legend.position = "bottom"
  )


# ================================================================
# 21. Panel (a): Annual Recall@10
# ================================================================

p_a <- ANNUAL %>%
  mutate(
    year = factor(year)
  ) %>%
  ggplot(
    aes(
      x = year,
      y = mean_recall,
      color = nutrient,
      group = nutrient
    )
  ) +
  geom_line(
    linewidth = 0.9
  ) +
  geom_point(
    size = 2.7
  ) +
  geom_errorbar(
    aes(
      ymin = low,
      ymax = high
    ),
    width = 0.08
  ) +
  geom_text(
    aes(
      label = paste0(
        "n=",
        n_defined
      )
    ),
    vjust = -1.2,
    size = 3,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = COLORS
  ) +
  scale_y_continuous(
    labels = scales::label_percent(),
    expand = expansion(
      mult = c(0.05, 0.15)
    )
  ) +
  labs(
    title = "(a) Annual F-Net Recall@10",
    x = "Independent test year",
    y = "Mean event Recall@10",
    color = NULL
  ) +
  theme_s6


# ================================================================
# 22. Panel (b): Network-method Precision@10
# ================================================================

p_b <- METHOD %>%
  mutate(
    method = factor(
      method,
      levels = c(
        "Random",
        "Phys-Net",
        "TE-Net",
        "F-Net"
      )
    )
  ) %>%
  ggplot(
    aes(
      x = method,
      y = mean_precision,
      color = nutrient
    )
  ) +
  geom_pointrange(
    aes(
      ymin = low,
      ymax = high
    ),
    position = position_dodge(
      width = 0.5
    ),
    linewidth = 0.7
  ) +
  scale_color_manual(
    values = COLORS
  ) +
  scale_y_continuous(
    labels = scales::label_percent()
  ) +
  labs(
    title = "(b) Independent-test Precision@10",
    subtitle = "Source-cluster bootstrap 95% CI",
    x = NULL,
    y = "Mean event Precision@10",
    color = NULL
  ) +
  theme_s6 +
  theme(
    axis.text.x = element_text(
      angle = 18,
      hjust = 1
    )
  )


# ================================================================
# 23. Budget plotting data
# ================================================================

PLOT_BUDGET <- BUDGET_SUMMARY %>%
  mutate(
    method = factor(
      method,
      levels = c(
        "Random",
        "Physical",
        "TE",
        "PCS"
      )
    )
  )


# ================================================================
# 24. Panel (c): Trigger coverage
# ================================================================

p_c <- ggplot(
  PLOT_BUDGET,
  aes(
    x = budget,
    y = 100 * trigger_mean,
    color = method,
    group = method
  )
) +
  geom_ribbon(
    data = PLOT_BUDGET %>%
      filter(method == "Random"),
    aes(
      ymin = 100 * trigger_low,
      ymax = 100 * trigger_high,
      fill = method
    ),
    alpha = 0.12,
    color = NA
  ) +
  geom_line(
    linewidth = 0.9
  ) +
  geom_point(
    size = 2.2
  ) +
  facet_wrap(
    ~ nutrient,
    nrow = 1
  ) +
  scale_color_manual(
    values = STRATEGY_COLORS
  ) +
  scale_fill_manual(
    values = STRATEGY_COLORS,
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = BUDGETS
  ) +
  coord_cartesian(
    ylim = c(0, 100)
  ) +
  labs(
    title = "(c) Trigger coverage of anomalous months",
    subtitle = "Fixed monitoring inputs; random-selection 95% interval",
    x = "Number of retained monitoring stations",
    y = "Triggered anomaly months (%)",
    color = "Selection strategy"
  ) +
  theme_s6


# ================================================================
# 25. Panel (d): Overall anomaly coverage
# ================================================================

p_d <- ggplot(
  PLOT_BUDGET,
  aes(
    x = budget,
    y = 100 * overall_mean_03,
    color = method,
    group = method
  )
) +
  geom_ribbon(
    data = PLOT_BUDGET %>%
      filter(method == "Random"),
    aes(
      ymin = 100 * overall_low_03,
      ymax = 100 * overall_high_03,
      fill = method
    ),
    alpha = 0.12,
    color = NA
  ) +
  geom_line(
    linewidth = 0.9
  ) +
  geom_point(
    size = 2.2
  ) +
  facet_wrap(
    ~ nutrient,
    nrow = 1
  ) +
  scale_color_manual(
    values = STRATEGY_COLORS
  ) +
  scale_fill_manual(
    values = STRATEGY_COLORS,
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = BUDGETS
  ) +
  coord_cartesian(
    ylim = c(0, 100)
  ) +
  labs(
    title = "(d) Fixed-budget anomaly observation or flagging",
    subtitle = "One global Top 10 per trigger month; 0–3-month window",
    x = "Number of retained monitoring stations",
    y = "Anomalous station-month coverage (%)",
    color = "Selection strategy"
  ) +
  theme_s6


# ================================================================
# 26. Assemble Figure S6
# ================================================================

FIG_S6 <- (
  (p_a | p_b) /
    (p_c | p_d)
) +
  patchwork::plot_annotation(
    title = paste0(
      "Figure S6 | Robustness of Event Screening ",
      "and Monitoring Efficiency"
    ),
    subtitle = paste0(
      "Training: 2011–2020 | ",
      "Independent testing: 2021–2023 | ",
      "Global Top-10 monitoring budget"
    ),
    caption = paste0(
      "Panels (a,b) use independent Fig. 6 events. ",
      "Panels (c,d) use monthly fixed-station alarms. ",
      "One global Top-10 hidden-target list is issued ",
      "per triggered month. ",
      "Random bands represent 95% station-selection ",
      "simulation intervals. ",
      "Flagged anomalies are retrospective opportunities, ",
      "not confirmed new observations."
    )
  )

print(FIG_S6)


# ================================================================
# 27. Export figures
# ================================================================

ggsave(
  file.path(
    OUT,
    "Figure_S6_revised.pdf"
  ),
  FIG_S6,
  width = 16,
  height = 11,
  bg = "white"
)

ggsave(
  file.path(
    OUT,
    "Figure_S6_revised.png"
  ),
  FIG_S6,
  width = 16,
  height = 11,
  dpi = 400,
  bg = "white"
)

PANELS <- list(
  a = p_a,
  b = p_b,
  c = p_c,
  d = p_d
)

for (nm in names(PANELS)) {
  
  ggsave(
    file.path(
      OUT,
      paste0(
        "Figure_S6_",
        nm,
        "_revised.pdf"
      )
    ),
    PANELS[[nm]],
    width = 8,
    height = 5.5,
    bg = "white"
  )
}


# ================================================================
# 28. Print results and checks
# ================================================================

cat(
  "\n====================================\n",
  "Figure S6 revised completed\n",
  "====================================\n"
)

cat("\nOutput:\n", OUT, "\n")

cat("\nAnnual Recall@10:\n")
print(
  as.data.frame(ANNUAL),
  row.names = FALSE
)

cat("\nMethod Precision@10:\n")
print(
  as.data.frame(METHOD),
  row.names = FALSE
)

cat("\nBudget summary:\n")
print(
  as.data.frame(BUDGET_SUMMARY),
  row.names = FALSE
)

cat("\nDenominator QA:\n")
print(
  as.data.frame(DENOMINATOR_QA),
  row.names = FALSE
)

cat(
  "\nMaximum monthly recommendation count:\n"
)

print(
  max(
    BUDGET_RESULT$
      max_recommendations_per_month
  )
)

cat("\nFinished.\n")
