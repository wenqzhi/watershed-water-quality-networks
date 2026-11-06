
# ================================================================
# FIG6 — STEP 1
# TRAINING NETWORK RECONSTRUCTION
#
# Training: 2011-10 to 2020-12
# Testing:  2021-01 to 2023-12
#
# TN / TP:
# 1. Training concentration matrix
# 2. TE-Net
# 3. Physical-network fusion
# 4. Source-conditioned PCS
# 5. Global influence scores
#
# No test-period concentration data used.
# ================================================================


# ================================================================
# 0. PACKAGES
# ================================================================

PKGS <- c(
  "readr",
  "dplyr",
  "tidyr",
  "lubridate",
  "igraph",
  "tibble",
  "stringr"
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

CSV <- file.path(
  ROOT,
  "NTN_2023_MonthlyLoadTable.csv"
)

P_FILE <- file.path(
  ROOT,
  "PNet_2011_2023_BlockTracing",
  "PNet_directed_reachable_pairs.csv"
)

OUT_ROOT <- file.path(
  ROOT,
  "Fig6_independent_training_2011_2020"
)

dir.create(
  OUT_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)


# ================================================================
# 2. PARAMETERS
# ================================================================

TRAIN_START <- as.Date("2011-10-01")
TRAIN_END   <- as.Date("2020-12-01")

TRAIN <- c(
  as.character(TRAIN_START),
  as.character(TRAIN_END)
)

MAX_LAG <- 3L

TE_P_THRESHOLD <- 0.10

# First QA run: 199
# Final formal run: 999
B <- 199L

SEED <- 20261023L

# Use Inf for all eligible stations.
# Use 12 for a small smoke test.
MAX_NODES <- Inf

# Original provisional fusion settings.
# Verify against your manuscript before final use.
OMEGA <- 0.5

BACKBONE_Q <- 0.85

VELOCITY_M_S <- 0.5

TRAVEL_DECAY_DAYS <- 10

PCS_DEPTH <- 5L

PCS_DISCOUNT <- 0.85

# Original code's analyte codes
PCODES <- c(
  TN = "P00600",
  TP = "P00665"
)

VALUE_COL <- "KConc"


# ================================================================
# 3. UTILITIES
# ================================================================

id8 <- function(x) {
  
  x <- trimws(as.character(x))
  
  x <- sub("\\.0+$", "", x)
  
  x[x %in% c(
    "", "NA", "NaN"
  )] <- NA_character_
  
  ix <- !is.na(x) &
    grepl("^[0-9]{1,8}$", x)
  
  x[ix] <- stringr::str_pad(
    x[ix],
    width = 8,
    side = "left",
    pad = "0"
  )
  
  x
}


check_columns <- function(d, required, label) {
  
  missing <- setdiff(
    required,
    names(d)
  )
  
  if (length(missing) > 0) {
    stop(
      label,
      " missing fields: ",
      paste(missing, collapse = ", ")
    )
  }
}


normalize_positive <- function(x) {
  
  x <- as.numeric(x)
  
  y <- numeric(length(x))
  
  ix <- which(
    is.finite(x) & x > 0
  )
  
  if (!length(ix)) {
    return(y)
  }
  
  lo <- min(x[ix])
  hi <- max(x[ix])
  
  if (hi > lo) {
    
    y[ix] <- (
      x[ix] - lo
    ) / (hi - lo)
    
  } else {
    
    y[ix] <- 1
  }
  
  y
}


# ================================================================
# 4. READ RAW CONCENTRATION DATA
# ================================================================

if (!file.exists(CSV)) {
  stop(
    "Raw NTN concentration file not found:\n",
    CSV
  )
}

raw_all <- readr::read_csv(
  CSV,
  col_types = readr::cols(
    STAID = readr::col_character(),
    PCODE = readr::col_character(),
    .default = readr::col_guess()
  ),
  show_col_types = FALSE
)

check_columns(
  raw_all,
  c(
    "STAID",
    "PCODE",
    "Year",
    "Month",
    VALUE_COL
  ),
  "NTN input"
)

raw_all <- raw_all %>%
  dplyr::mutate(
    STAID = id8(STAID),
    PCODE = as.character(PCODE),
    Year = as.integer(Year),
    Month = as.integer(Month),
    Date = as.Date(
      sprintf(
        "%04d-%02d-01",
        Year,
        Month
      )
    ),
    KConc = as.numeric(KConc)
  )

cat("\n====================================\n")
cat("RAW NTN DATA LOADED\n")
cat("====================================\n")

cat(
  "Rows:",
  nrow(raw_all),
  "\n"
)


# ================================================================
# 5. PREPARE()
# ================================================================
#
# IMPORTANT:
# This function was missing in the previous script.
#
# It:
# - reads data already filtered to a nutrient;
# - extracts the requested training months;
# - keeps stations with complete finite observations;
# - does not fill or interpolate missing concentrations;
# - returns dates, station IDs and matrix.
#
# ================================================================

prepare <- function(
    period,
    eligible = NULL
) {
  
  if (!exists("d", inherits = TRUE)) {
    stop(
      "Nutrient-specific data object d does not exist."
    )
  }
  
  start_date <- as.Date(period[1])
  end_date <- as.Date(period[2])
  
  expected_dates <- seq(
    from = start_date,
    to = end_date,
    by = "month"
  )
  
  n_expected <- length(expected_dates)
  
  # --------------------------------
  # 5.1 Filter the training window
  # --------------------------------
  
  dd <- d %>%
    dplyr::filter(
      Date >= start_date,
      Date <= end_date,
      !is.na(STAID)
    ) %>%
    dplyr::select(
      STAID,
      Date,
      KConc
    )
  
  if (!is.null(eligible)) {
    
    eligible <- id8(eligible)
    
    dd <- dd %>%
      dplyr::filter(
        STAID %in% eligible
      )
  }
  
  # --------------------------------
  # 5.2 Duplicate checks
  # --------------------------------
  
  duplicate_rows <- dd %>%
    dplyr::count(
      STAID,
      Date,
      name = "n"
    ) %>%
    dplyr::filter(n > 1)
  
  if (nrow(duplicate_rows) > 0) {
    
    print(
      utils::head(
        duplicate_rows,
        20
      )
    )
    
    stop(
      "Duplicate station/month records ",
      "in training data."
    )
  }
  
  # --------------------------------
  # 5.3 Keep finite, nonnegative data
  # --------------------------------
  
  dd <- dd %>%
    dplyr::mutate(
      valid = (
        is.finite(KConc) &
          KConc >= 0
      )
    )
  
  station_qa <- dd %>%
    dplyr::group_by(STAID) %>%
    dplyr::summarise(
      n_months = dplyr::n_distinct(Date),
      
      n_valid = sum(valid),
      
      complete = (
        n_months == n_expected &
          n_valid == n_expected
      ),
      
      .groups = "drop"
    )
  
  complete_ids <- station_qa %>%
    dplyr::filter(complete) %>%
    dplyr::arrange(STAID) %>%
    dplyr::pull(STAID)
  
  cat("\nTraining window:\n")
  
  cat(
    as.character(start_date),
    "to",
    as.character(end_date),
    "\n"
  )
  
  cat(
    "Expected months:",
    n_expected,
    "\n"
  )
  
  cat(
    "Stations with data:",
    nrow(station_qa),
    "\n"
  )
  
  cat(
    "Complete stations:",
    length(complete_ids),
    "\n"
  )
  
  if (length(complete_ids) < 3) {
    stop(
      "Fewer than three complete stations ",
      "for network calculation."
    )
  }
  
  # --------------------------------
  # 5.4 Optional QA subset
  # --------------------------------
  
  if (is.finite(MAX_NODES)) {
    
    complete_ids <- utils::head(
      complete_ids,
      as.integer(MAX_NODES)
    )
    
    cat(
      "QA subset stations:",
      length(complete_ids),
      "\n"
    )
  }
  
  # --------------------------------
  # 5.5 Construct wide matrix
  # --------------------------------
  
  wide <- dd %>%
    dplyr::filter(
      STAID %in% complete_ids
    ) %>%
    dplyr::select(
      Date,
      STAID,
      KConc
    ) %>%
    tidyr::pivot_wider(
      names_from = STAID,
      values_from = KConc
    ) %>%
    dplyr::arrange(Date)
  
  # Ensure correct month ordering.
  
  if (!identical(
    as.Date(wide$Date),
    expected_dates
  )) {
    
    stop(
      "Training dates are not identical ",
      "to the required monthly sequence."
    )
  }
  
  wide <- wide %>%
    dplyr::select(
      Date,
      dplyr::all_of(complete_ids)
    )
  
  x <- as.matrix(
    wide[, -1, drop = FALSE]
  )
  
  storage.mode(x) <- "double"
  
  if (any(!is.finite(x))) {
    stop(
      "Prepared training matrix contains ",
      "NA / NaN / Inf."
    )
  }
  
  if (any(x < 0)) {
    stop(
      "Negative concentration detected."
    )
  }
  
  cat(
    "Final training matrix:",
    nrow(x),
    "months x",
    ncol(x),
    "stations\n"
  )
  
  list(
    x = x,
    dates = expected_dates,
    ids = colnames(x),
    station_qa = station_qa
  )
}


# ================================================================
# 6. PREPROCESS TIME SERIES
# ================================================================

preprocess <- function(a, dates) {
  
  if (any(a < 0, na.rm = TRUE)) {
    stop(
      "Negative concentrations found."
    )
  }
  
  z <- log1p(a)
  
  mo <- lubridate::month(dates)
  
  tt <- seq_len(nrow(z))
  
  for (j in seq_len(ncol(z))) {
    
    x <- z[, j]
    
    # Remove monthly climatological means.
    x <- x - ave(
      x,
      mo,
      FUN = mean
    )
    
    # Remove linear trend.
    fit <- stats::lm(
      x ~ tt
    )
    
    z[, j] <- stats::residuals(fit)
  }
  
  z
}


# ================================================================
# 7. DISCRETIZE TIME SERIES
# ================================================================

discretize <- function(z) {
  
  out <- apply(
    z,
    2,
    function(x) {
      
      br <- stats::quantile(
        x,
        probs = c(1 / 3, 2 / 3),
        names = FALSE
      )
      
      if (br[1] >= br[2]) {
        
        r <- rank(
          x,
          ties.method = "average"
        )
        
        as.integer(
          cut(
            r,
            breaks = c(
              -Inf,
              length(x) / 3,
              2 * length(x) / 3,
              Inf
            ),
            labels = FALSE
          )
        )
        
      } else {
        
        as.integer(
          cut(
            x,
            breaks = c(
              -Inf,
              br[1],
              br[2],
              Inf
            ),
            labels = FALSE
          )
        )
      }
    }
  )
  
  if (is.null(dim(out))) {
    
    out <- matrix(
      out,
      ncol = 1
    )
  }
  
  colnames(out) <- colnames(z)
  
  storage.mode(out) <- "integer"
  
  if (anyNA(out)) {
    stop(
      "Discretization produced missing values."
    )
  }
  
  out
}


# ================================================================
# 8. CONDITIONAL MUTUAL INFORMATION
# ================================================================

cmi <- function(x, y, z, K = 3L) {
  
  n <- length(x)
  
  if (
    n == 0L ||
    anyNA(x) ||
    anyNA(y) ||
    anyNA(z)
  ) {
    return(NA_real_)
  }
  
  index <- x +
    K * (y - 1L) +
    K * K * (z - 1L)
  
  arr <- array(
    tabulate(
      index,
      nbins = K^3
    ),
    dim = c(K, K, K)
  )
  
  p <- arr / n
  
  pxz <- apply(
    p,
    c(1, 3),
    sum
  )
  
  pyz <- apply(
    p,
    c(2, 3),
    sum
  )
  
  pz <- apply(
    p,
    3,
    sum
  )
  
  val <- 0
  
  for (k in seq_len(K)) {
    
    for (j in seq_len(K)) {
      
      for (i in seq_len(K)) {
        
        if (
          p[i, j, k] > 0 &&
          pxz[i, k] > 0 &&
          pyz[j, k] > 0 &&
          pz[k] > 0
        ) {
          
          val <- val +
            p[i, j, k] *
            log(
              p[i, j, k] * pz[k] /
                (
                  pxz[i, k] *
                    pyz[j, k]
                )
            )
        }
      }
    }
  }
  
  max(
    0,
    val
  )
}


# ================================================================
# 9. TRANSFER ENTROPY AT DIFFERENT LAGS
# ================================================================

lag_te <- function(
    source,
    target,
    maxlag = MAX_LAG
) {
  
  tix <- seq.int(
    maxlag + 1L,
    length(target)
  )
  
  vapply(
    seq_len(maxlag),
    function(lag) {
      
      cmi(
        source[tix - lag],
        target[tix],
        target[tix - 1L]
      )
    },
    numeric(1)
  )
}


# ================================================================
# 10. SURROGATE TEST
# ================================================================

surrogate_shift <- function(n) {
  
  good <- seq.int(
    1L,
    n - 1L
  )
  
  dist <- pmin(
    good,
    n - good
  )
  
  good <- good[
    dist > max(
      12L,
      MAX_LAG + 1L
    )
  ]
  
  if (!length(good)) {
    stop(
      "Series too short for surrogate shift."
    )
  }
  
  sample(
    good,
    1
  )
}


shift_vec <- function(x, k) {
  
  n <- length(x)
  
  x[c(
    (n - k + 1L):n,
    seq_len(n - k)
  )]
}


# ================================================================
# 11. CALCULATE TRAINING TE NETWORK
# ================================================================

calc_network <- function(
    period,
    label,
    eligible = NULL
) {
  
  raw <- prepare(
    period,
    eligible
  )
  
  z <- preprocess(
    raw$x,
    raw$dates
  )
  
  mat <- discretize(z)
  
  n <- ncol(mat)
  
  nt <- nrow(mat)
  
  if (nt <= 2 * MAX_LAG + 28L) {
    stop(
      "Too few months for TE calculation."
    )
  }
  
  n_pairs <- n * (n - 1L)
  
  cat("\n====================================\n")
  cat("CALCULATING TE NETWORK\n")
  cat("====================================\n")
  
  cat(
    label, "\n",
    "Stations:", n, "\n",
    "Months:", nt, "\n",
    "Directed pairs:", n_pairs, "\n",
    "Surrogates:", B, "\n"
  )
  
  pairs <- which(
    row(matrix(0, n, n)) !=
      col(matrix(0, n, n)),
    arr.ind = TRUE
  )
  
  set.seed(SEED)
  
  res <- vector(
    "list",
    nrow(pairs)
  )
  
  for (k in seq_len(nrow(pairs))) {
    
    i <- pairs[k, 1]
    j <- pairs[k, 2]
    
    xs <- mat[, i]
    yt <- mat[, j]
    
    tv <- lag_te(
      xs,
      yt
    )
    
    score <- max(tv)
    
    best <- which.max(tv)
    
    null <- replicate(
      B,
      {
        shift <- surrogate_shift(
          length(xs)
        )
        
        xs_null <- shift_vec(
          xs,
          shift
        )
        
        max(
          lag_te(
            xs_null,
            yt
          )
        )
      }
    )
    
    p <- (
      1 + sum(null >= score)
    ) / (
      B + 1
    )
    
    res[[k]] <- tibble::tibble(
      from = colnames(mat)[i],
      to = colnames(mat)[j],
      lag_months = best,
      TE_nats = score,
      TE_bias_corrected = max(
        0,
        score - mean(null)
      ),
      p_value = p
    )
    
    if (k %% 500 == 0) {
      
      cat(
        label,
        ": ",
        k,
        "/",
        nrow(pairs),
        " pairs tested\n",
        sep = ""
      )
    }
  }
  
  ed <- dplyr::bind_rows(
    res
  ) %>%
    dplyr::mutate(
      q_value = p.adjust(
        p_value,
        method = "BH"
      ),
      significant = (
        p_value < TE_P_THRESHOLD
      )
    )
  
  readr::write_csv(
    ed,
    file.path(
      OUT,
      paste0(
        label,
        "_all_edges.csv"
      )
    )
  )
  
  readr::write_csv(
    dplyr::filter(
      ed,
      significant
    ),
    file.path(
      OUT,
      paste0(
        label,
        "_significant_edges.csv"
      )
    )
  )
  
  # Training matrix
  
  input_matrix <- tibble::as_tibble(
    raw$x,
    .name_repair = "minimal"
  )
  
  input_matrix <- dplyr::bind_cols(
    tibble::tibble(
      Date = raw$dates
    ),
    input_matrix
  )
  
  readr::write_csv(
    input_matrix,
    file.path(
      OUT,
      paste0(
        label,
        "_input_matrix.csv"
      )
    )
  )
  
  readr::write_csv(
    raw$station_qa,
    file.path(
      OUT,
      paste0(
        label,
        "_station_coverage_QA.csv"
      )
    )
  )
  
  readr::write_csv(
    tibble::tibble(
      STAID = raw$ids
    ),
    file.path(
      OUT,
      paste0(
        label,
        "_station_ids.csv"
      )
    )
  )
  
  # Network
  
  ed_sig <- ed %>%
    dplyr::filter(
      significant,
      TE_bias_corrected > 0
    )
  
  g <- igraph::graph_from_data_frame(
    ed_sig,
    directed = TRUE,
    vertices = data.frame(
      name = raw$ids
    )
  )
  
  nodes <- tibble::tibble(
    STAID = raw$ids,
    
    indegree = igraph::degree(
      g,
      mode = "in"
    ),
    
    outdegree = igraph::degree(
      g,
      mode = "out"
    ),
    
    in_strength = igraph::strength(
      g,
      mode = "in",
      weights = igraph::E(g)$TE_bias_corrected
    ),
    
    out_strength = igraph::strength(
      g,
      mode = "out",
      weights = igraph::E(g)$TE_bias_corrected
    )
  )
  
  readr::write_csv(
    nodes,
    file.path(
      OUT,
      paste0(
        label,
        "_node_metrics.csv"
      )
    )
  )
  
  igraph::write_graph(
    g,
    file.path(
      OUT,
      paste0(
        label,
        "_network.graphml"
      )
    ),
    format = "graphml"
  )
  
  metrics <- tibble::tibble(
    nutrient = NU,
    window = label,
    train_start = TRAIN_START,
    train_end = TRAIN_END,
    months = nt,
    stations = n,
    candidate_edges = nrow(ed),
    significant_p_edges = sum(
      ed$significant
    ),
    positive_significant_edges = nrow(ed_sig),
    density = nrow(ed_sig) /
      (n * (n - 1L)),
    surrogates = B,
    max_lag = MAX_LAG
  )
  
  readr::write_csv(
    metrics,
    file.path(
      OUT,
      paste0(
        label,
        "_summary.csv"
      )
    )
  )
  
  list(
    edges = ed,
    nodes = nodes,
    metrics = metrics,
    ids = raw$ids
  )
}


# ================================================================
# 12. FUSION NETWORK
# ================================================================

build_fusion <- function(
    nu,
    te_edges,
    node_ids
) {
  
  out <- file.path(
    OUT_ROOT,
    nu
  )
  
  dir.create(
    out,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  nodes <- sort(
    unique(
      id8(node_ids)
    )
  )
  
  pairs <- tidyr::expand_grid(
    from = nodes,
    to = nodes
  ) %>%
    dplyr::filter(
      from != to
    )
  
  # --------------------------------
  # Training TE layer
  # --------------------------------
  
  te <- te_edges %>%
    dplyr::filter(
      significant %in% TRUE,
      TE_bias_corrected > 0
    ) %>%
    dplyr::transmute(
      from = id8(from),
      to = id8(to),
      T_raw = TE_bias_corrected,
      TE_p = p_value,
      TE_lag_months = lag_months
    )
  
  full <- pairs %>%
    dplyr::left_join(
      te,
      by = c("from", "to")
    ) %>%
    dplyr::mutate(
      T_raw = dplyr::coalesce(
        T_raw,
        0
      )
    )
  
  # --------------------------------
  # Physical layer
  # --------------------------------
  
  if (!file.exists(P_FILE)) {
    
    warning(
      "Physical network file missing:\n",
      P_FILE,
      "\nTraining TE-Net has been saved, ",
      "but F-Net cannot be generated."
    )
    
    readr::write_csv(
      full,
      file.path(
        out,
        "01_training_TE_pair_layer.csv"
      )
    )
    
    return(
      invisible(NULL)
    )
  }
  
  phy_raw <- readr::read_csv(
    P_FILE,
    col_types = readr::cols(
      from = readr::col_character(),
      to = readr::col_character(),
      .default = readr::col_guess()
    ),
    show_col_types = FALSE
  )
  
  check_columns(
    phy_raw,
    c(
      "from",
      "to",
      "river_distance_km"
    ),
    "Physical network"
  )
  
  phy <- phy_raw %>%
    dplyr::transmute(
      from = id8(from),
      to = id8(to),
      river_distance_km = as.numeric(
        river_distance_km
      )
    ) %>%
    dplyr::filter(
      from != to,
      is.finite(river_distance_km),
      river_distance_km > 0
    ) %>%
    dplyr::group_by(
      from,
      to
    ) %>%
    dplyr::summarise(
      river_distance_km = min(
        river_distance_km
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      travel_days =
        river_distance_km * 1000 /
        (
          VELOCITY_M_S * 86400
        ),
      
      P_raw = exp(
        -travel_days /
          TRAVEL_DECAY_DAYS
      )
    )
  
  full <- full %>%
    dplyr::left_join(
      phy,
      by = c("from", "to")
    )
  
  full <- full %>%
    dplyr::mutate(
      physical_verified_path = !is.na(P_raw),
      
      physical_status = ifelse(
        physical_verified_path,
        "traced_downstream",
        "unresolved_or_no_path"
      ),
      
      P_effective = dplyr::coalesce(
        P_raw,
        0
      )
    )
  
  # NA physical reachability is not
  # interpreted as verified nonconnection.
  
  # --------------------------------
  # Normalize and fuse
  # --------------------------------
  
  full$T_norm <- normalize_positive(
    full$T_raw
  )
  
  full$P_norm <- normalize_positive(
    full$P_effective
  )
  
  full$W_fus <-
    OMEGA * full$T_norm +
    (1 - OMEGA) * full$P_norm
  
  positive <- full$W_fus[
    is.finite(full$W_fus) &
      full$W_fus > 0
  ]
  
  if (!length(positive)) {
    stop(
      "No positive fusion edges for ",
      nu
    )
  }
  
  cutoff <- as.numeric(
    stats::quantile(
      positive,
      probs = BACKBONE_Q,
      names = FALSE
    )
  )
  
  full$backbone <- (
    full$W_fus > 0 &
      full$W_fus >= cutoff
  )
  
  readr::write_csv(
    full,
    file.path(
      out,
      "02_training_FNet_all_pairs.csv"
    )
  )
  
  readr::write_csv(
    dplyr::filter(
      full,
      backbone
    ),
    file.path(
      out,
      "03_training_FNet_backbone.csv"
    )
  )
  
  # ============================================================
  # 13. SOURCE-CONDITIONED PCS
  # ============================================================
  
  nn <- length(nodes)
  
  W <- matrix(
    0,
    nrow = nn,
    ncol = nn,
    dimnames = list(
      nodes,
      nodes
    )
  )
  
  ix <- match(
    full$from,
    nodes
  )
  
  jx <- match(
    full$to,
    nodes
  )
  
  W[
    cbind(ix, jx)
  ] <- full$W_fus
  
  rs <- rowSums(W)
  
  P <- W
  
  good <- rs > 0
  
  P[good, ] <- W[
    good,
    ,
    drop = FALSE
  ] / rs[good]
  
  P[!good, ] <- 0
  
  # Discounted walk influence
  S <- matrix(
    0,
    nrow = nn,
    ncol = nn,
    dimnames = list(
      nodes,
      nodes
    )
  )
  
  M <- diag(nn)
  
  for (ell in seq_len(PCS_DEPTH)) {
    
    M <- M %*% P
    
    S <- S +
      PCS_DISCOUNT^(ell - 1L) * M
  }
  
  diag(S) <- 0
  
  source_rank <- tidyr::expand_grid(
    source = nodes,
    target = nodes
  ) %>%
    dplyr::filter(
      source != target
    ) %>%
    dplyr::mutate(
      source_PCS = S[
        cbind(
          match(source, nodes),
          match(target, nodes)
        )
      ]
    ) %>%
    dplyr::group_by(source) %>%
    dplyr::arrange(
      dplyr::desc(source_PCS),
      target,
      .by_group = TRUE
    ) %>%
    dplyr::mutate(
      rank_FNet = dplyr::row_number()
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      nutrient = nu
    )
  
  readr::write_csv(
    source_rank,
    file.path(
      out,
      "04_training_source_conditioned_PCS.csv"
    )
  )
  
  # ============================================================
  # 14. GLOBAL INFLUENCE SCORES
  # ============================================================
  #
  # PCS_in is the discounted incoming influence.
  # PCS_out is the discounted outgoing influence.
  #
  # This must be checked against your manuscript's
  # final Global PCS definition before Fig. 6.
  # ============================================================
  
  global <- tibble::tibble(
    STAID = nodes,
    nutrient = nu,
    
    PCS_in = colSums(S),
    PCS_out = rowSums(S),
    
    FNet_out_strength = rowSums(W),
    FNet_in_strength = colSums(W)
  )
  
  readr::write_csv(
    global,
    file.path(
      out,
      "05_training_global_PCS.csv"
    )
  )
  
  # ============================================================
  # 15. NETWORK QA
  # ============================================================
  
  qa <- tibble::tibble(
    nutrient = nu,
    
    training_start = TRAIN_START,
    training_end = TRAIN_END,
    
    nodes = nn,
    
    candidate_pairs = nrow(full),
    
    TE_positive_significant_edges = nrow(te),
    
    physical_traced = sum(
      full$physical_verified_path
    ),
    
    positive_fused_edges = sum(
      full$W_fus > 0
    ),
    
    backbone_edges = sum(
      full$backbone
    ),
    
    omega = OMEGA,
    
    backbone_quantile = BACKBONE_Q,
    
    PCS_depth = PCS_DEPTH,
    
    PCS_discount = PCS_DISCOUNT,
    
    surrogate_iterations = B
  )
  
  readr::write_csv(
    qa,
    file.path(
      out,
      "06_training_network_QA.csv"
    )
  )
  
  cat("\n====================================\n")
  
  cat(
    nu,
    " TRAINING F-NET COMPLETED\n"
  )
  
  print(qa)
  
  cat(
    "Output:",
    out,
    "\n"
  )
}


# ================================================================
# 16. RUN TN AND TP
# ================================================================

for (NU in c("TN", "TP")) {
  
  cat("\n\n====================================\n")
  cat("START:", NU, "\n")
  cat("====================================\n")
  
  current_pcode <- unname(
    PCODES[[NU]]
  )
  
  OUT <- file.path(
    OUT_ROOT,
    NU
  )
  
  dir.create(
    OUT,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # Important:
  # Use .data$PCODE and current_pcode
  # to avoid data masking conflicts.
  
  d <- raw_all %>%
    dplyr::filter(
      .data$PCODE == .env$current_pcode
    ) %>%
    dplyr::select(
      STAID,
      Date,
      KConc
    )
  
  cat(
    "Nutrient:",
    NU,
    "\n"
  )
  
  cat(
    "PCODE:",
    current_pcode,
    "\n"
  )
  
  cat(
    "Raw records:",
    nrow(d),
    "\n"
  )
  
  if (nrow(d) == 0) {
    stop(
      "No input records for ",
      NU,
      " / ",
      current_pcode
    )
  }
  
  # Duplicate QA before network calculation
  
  dup <- d %>%
    dplyr::count(
      STAID,
      Date
    ) %>%
    dplyr::filter(n > 1)
  
  if (nrow(dup) > 0) {
    
    print(
      utils::head(
        dup,
        20
      )
    )
    
    stop(
      "Duplicate station/month for ",
      NU
    )
  }
  
  # --------------------------------
  # Calculate TRAINING TE-Net
  # --------------------------------
  
  fit <- calc_network(
    period = TRAIN,
    label = "train_2011_2020"
  )
  
  # --------------------------------
  # Calculate TRAINING F-Net
  # --------------------------------
  
  build_fusion(
    nu = NU,
    te_edges = fit$edges,
    node_ids = fit$ids
  )
  
  cat(
    "\nFINISHED:",
    NU,
    "\n"
  )
}


# ================================================================
# 17. FINAL STATUS
# ================================================================

cat("\n\n====================================\n")
cat("TRAINING NETWORK STEP COMPLETED\n")
cat("====================================\n")

cat(
  "Training period:",
  as.character(TRAIN_START),
  "to",
  as.character(TRAIN_END),
  "\n"
)

cat(
  "Test concentrations used: NO\n"
)

cat(
  "Output root:\n",
  OUT_ROOT,
  "\n"
)

for (NU in c("TN", "TP")) {
  
  cat(
    "\n",
    NU,
    " network summary:\n",
    sep = ""
  )
  
  qa_file <- file.path(
    OUT_ROOT,
    NU,
    "train_2011_2020_summary.csv"
  )
  
  if (file.exists(qa_file)) {
    
    print(
      readr::read_csv(
        qa_file,
        show_col_types = FALSE
      )
    )
  }
  
  fusion_file <- file.path(
    OUT_ROOT,
    NU,
    "06_training_network_QA.csv"
  )
  
  if (file.exists(fusion_file)) {
    
    print(
      readr::read_csv(
        fusion_file,
        show_col_types = FALSE
      )
    )
    
  } else {
    
    cat(
      "F-Net QA not available. ",
      "Check P_FILE.\n"
    )
  }
}

cat("\nDONE.\n")
