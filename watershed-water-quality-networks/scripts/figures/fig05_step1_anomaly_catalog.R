# ============================================================
# FIGURE 5 — STEP 5.1
# REAL NTN NUTRIENT ANOMALY EVENT CATALOG
#
# Goal:
# Identify real TN/TP anomaly triggers and
# their observation coverage for future
# source-conditioned F-Net monitoring ranking.
#
# This step does NOT:
# - Recalculate TE or F-Net
# - Predict future nutrient concentrations
# - Infer physical nutrient transport
# - Select successful cases after inspection
#
# All events and observation availability
# are preserved for subsequent validation.
# ============================================================


# ============================================================
# 0. PACKAGES
# ============================================================

PKGS <- c(
  "readr", "dplyr", "tidyr",
  "tibble", "ggplot2", "lubridate",
  "scales", "stringr"
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

OUT <- file.path(
  ROOT,
  "Fig5_anomaly_trigger",
  "Step5_1_event_catalog"
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

STATION_FILE <- file.path(
  ROOT,
  "NTN-station.txt"
)

PCS_FILE <- file.path(
  ROOT,
  "Fig4_PCS_recomputed",
  "Global_PCS_TN_TP_recomputed.csv"
)


# ============================================================
# 2. SETTINGS
# ============================================================

START_DATE <- as.Date("2011-10-01")
END_DATE   <- as.Date("2023-09-01")

TRAIN_END <- as.Date("2019-12-01")

TEST_START <- as.Date("2020-01-01")

# Main anomaly threshold
Z_TRIGGER <- 1.5

# Previous month must be below this level
Z_PREVIOUS_MAX <- 1.0

# Sensitivity analysis
Z_SENSITIVITY <- c(1.0, 1.5, 2.0)

# Require enough baseline observations
MIN_TRAIN_TOTAL <- 60L

# Each calendar-month climatology
# requires a sufficient number of years
MIN_MONTH_SAMPLES <- 5L

# Minimum nonzero within-calendar-month SD
MIN_MONTH_SD <- 1e-8

# Evaluate availability for the current
# month plus the next 3 months.
MAX_FOLLOWUP_MONTHS <- 3L

# Only proceed with events that have
# sufficient target-station coverage.
MIN_TARGETS_OBSERVED <- 10L

# Set TRUE if you want every source
# to have an observed previous month.
REQUIRE_PREVIOUS_OBSERVATION <- TRUE

# Do not interpolate missing months.
ALLOW_INTERPOLATION <- FALSE


# ============================================================
# 3. HELPERS
# ============================================================

clean_id <- function(x) {
  
  x <- trimws(as.character(x))
  
  x <- sub("\\.0+$", "", x)
  
  x[x %in% c("", "NA", "NaN")] <- NA_character_
  
  x
}

save_csv <- function(d, name) {
  
  readr::write_csv(
    d,
    file.path(OUT, name),
    na = ""
  )
}

if (!file.exists(NTN_FILE)) {
  stop("NTN file missing: ", NTN_FILE)
}

if (!file.exists(PCS_FILE)) {
  stop("PCS file missing: ", PCS_FILE)
}


# ============================================================
# 4. LOAD F-NET STATION SET
# ============================================================

pcs <- readr::read_csv(
  PCS_FILE,
  
  col_types = readr::cols(
    STAID = readr::col_character(),
    .default = readr::col_guess()
  ),
  
  show_col_types = FALSE
) %>%
  transmute(
    STAID = clean_id(STAID),
    
    nutrient = toupper(
      as.character(nutrient)
    ),
    
    PCS = as.numeric(PCS_recomputed)
  ) %>%
  filter(
    nutrient %in% c("TN", "TP"),
    !is.na(STAID),
    is.finite(PCS)
  ) %>%
  distinct(STAID, nutrient, .keep_all = TRUE)

cat("\nF-Net network stations:\n")

print(
  pcs %>% count(nutrient)
)


# ============================================================
# 5. READ REAL NTN MONTHLY CONCENTRATION
# ============================================================

ntn_raw <- readr::read_csv(
  NTN_FILE,
  
  col_types = readr::cols(
    STAID = readr::col_character(),
    PCODE = readr::col_character(),
    .default = readr::col_guess()
  ),
  
  show_col_types = FALSE
)

required <- c(
  "STAID", "PCODE",
  "Year", "Month", "KConc"
)

missing <- setdiff(required, names(ntn_raw))

if (length(missing)) {
  stop(
    "Missing NTN columns: ",
    paste(missing, collapse = ", ")
  )
}

ntn <- ntn_raw %>%
  transmute(
    STAID = clean_id(STAID),
    
    PCODE = trimws(
      as.character(PCODE)
    ),
    
    year = as.integer(Year),
    
    month = as.integer(Month),
    
    concentration = suppressWarnings(
      as.numeric(KConc)
    )
  ) %>%
  mutate(
    nutrient = case_when(
      PCODE == "P00600" ~ "TN",
      PCODE == "P00665" ~ "TP",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    nutrient %in% c("TN", "TP"),
    year >= 2011,
    year <= 2023,
    month %in% 1:12
  ) %>%
  mutate(
    date = as.Date(
      sprintf(
        "%04d-%02d-01",
        year, month
      )
    )
  ) %>%
  filter(
    date >= START_DATE,
    date <= END_DATE
  ) %>%
  inner_join(
    pcs,
    by = c("STAID", "nutrient")
  )

# Valid numeric concentrations only.
# Negative concentrations are treated as
# invalid; zero is retained provisionally.

ntn <- ntn %>%
  mutate(
    concentration = if_else(
      is.finite(concentration) &
        concentration >= 0,
      concentration,
      NA_real_
    )
  )

# No duplicate station/nutrient/month rows.
dups <- ntn %>%
  count(
    STAID,
    nutrient,
    date
  ) %>%
  filter(n > 1)

if (nrow(dups) > 0) {
  
  save_csv(
    dups,
    "QA_duplicate_station_months.csv"
  )
  
  stop(
    "Duplicate NTN station-month records. ",
    "Inspect QA_duplicate_station_months.csv."
  )
}

save_csv(
  ntn,
  "01_NTN_real_monthly_concentrations.csv"
)


# ============================================================
# 6. OBSERVATION COVERAGE
# ============================================================

all_months <- tibble(
  date = seq(
    START_DATE,
    END_DATE,
    by = "month"
  )
)

full_grid <- tidyr::crossing(
  pcs %>% select(STAID, nutrient),
  all_months
) %>%
  mutate(
    year = lubridate::year(date),
    
    month = lubridate::month(date)
  )

ntn_full <- full_grid %>%
  left_join(
    ntn %>%
      select(
        STAID,
        nutrient,
        date,
        concentration,
        PCS
      ),
    by = c(
      "STAID",
      "nutrient",
      "date"
    )
  ) %>%
  left_join(
    pcs,
    by = c("STAID", "nutrient"),
    suffix = c("", ".full")
  ) %>%
  select(
    STAID,
    nutrient,
    date,
    year,
    month,
    concentration,
    PCS = PCS.full
  )

coverage <- ntn_full %>%
  group_by(
    nutrient,
    STAID
  ) %>%
  summarise(
    n_total_months = n(),
    
    n_valid = sum(
      is.finite(concentration)
    ),
    
    n_train = sum(
      date <= TRAIN_END &
        is.finite(concentration)
    ),
    
    n_test = sum(
      date >= TEST_START &
        is.finite(concentration)
    ),
    
    valid_fraction =
      n_valid / n_total_months,
    
    .groups = "drop"
  )

save_csv(
  coverage,
  "02_station_coverage.csv"
)


# ============================================================
# 7. TRAINING-ONLY MONTHLY CLIMATOLOGY
# ============================================================

train_data <- ntn_full %>%
  filter(
    date <= TRAIN_END,
    is.finite(concentration)
  )

train_total <- train_data %>%
  count(
    STAID,
    nutrient,
    name = "n_train_total"
  )

climatology <- train_data %>%
  group_by(
    STAID,
    nutrient,
    month
  ) %>%
  summarise(
    n_month_train = n(),
    
    mean_month = mean(concentration),
    
    sd_month = if (n() >= 2) {
      sd(concentration)
    } else {
      NA_real_
    },
    
    .groups = "drop"
  ) %>%
  left_join(
    train_total,
    by = c("STAID", "nutrient")
  ) %>%
  mutate(
    valid_climatology = (
      n_train_total >= MIN_TRAIN_TOTAL &
        n_month_train >= MIN_MONTH_SAMPLES &
        is.finite(sd_month) &
        sd_month > MIN_MONTH_SD
    )
  )

save_csv(
  climatology,
  "03_training_monthly_climatology.csv"
)


# ============================================================
# 8. REAL DESEASONALIZED ANOMALIES
# ============================================================

anomaly <- ntn_full %>%
  left_join(
    climatology,
    by = c(
      "STAID",
      "nutrient",
      "month"
    )
  ) %>%
  mutate(
    Z = if_else(
      !is.na(valid_climatology) &
        valid_climatology &
        is.finite(concentration),
      (concentration - mean_month) /
        sd_month,
      NA_real_
    ),
    
    period = case_when(
      date <= TRAIN_END ~ "Training",
      date >= TEST_START ~ "Evaluation",
      TRUE ~ "Other"
    )
  )

save_csv(
  anomaly,
  "04_all_monthly_concentration_anomalies.csv"
)


# ============================================================
# 9. EVENT ONSET DETECTION
# ============================================================

# Use the previous CALENDAR month,
# never merely the previous nonmissing row.
#
# All calendar months were explicitly created
# in Step 6.

event_series <- anomaly %>%
  arrange(
    nutrient,
    STAID,
    date
  ) %>%
  group_by(
    nutrient,
    STAID
  ) %>%
  mutate(
    Z_previous = lag(Z, 1),
    
    previous_date = lag(date, 1),
    
    concentration_previous =
      lag(concentration, 1)
  ) %>%
  ungroup()

event_series <- event_series %>%
  mutate(
    triggered = (
      period == "Evaluation" &
        is.finite(Z) &
        Z > Z_TRIGGER &
        is.finite(Z_previous) &
        Z_previous <= Z_PREVIOUS_MAX &
        previous_date ==
        (date %m-% months(1))
    )
  )

events <- event_series %>%
  filter(triggered) %>%
  transmute(
    nutrient,
    
    source = STAID,
    
    trigger_date = date,
    
    trigger_year = year,
    
    trigger_month = month,
    
    source_concentration = concentration,
    
    source_Z = Z,
    
    previous_Z = Z_previous,
    
    source_PCS = PCS
  ) %>%
  arrange(
    nutrient,
    trigger_date,
    source
  ) %>%
  mutate(
    event_id = sprintf(
      "EV%04d", row_number()
    ),
    .before = 1
  )

save_csv(
  events,
  "05_all_real_anomaly_triggers.csv"
)

cat("\nAll real anomaly triggers:\n")

print(
  events %>%
    count(nutrient)
)


# ============================================================
# 10. OBSERVATION AVAILABILITY AFTER TRIGGER
# ============================================================

# For each event:
# - exclude source station
# - count valid target stations
#   at each lag k = 0,1,2,3 months
#
# No target anomaly outcomes are used to
# select a source event.

candidate_stations <- pcs %>%
  select(
    target = STAID,
    nutrient
  )

event_targets <- events %>%
  select(
    event_id,
    nutrient,
    source,
    trigger_date
  ) %>%
  inner_join(
    candidate_stations,
    by = "nutrient"
  ) %>%
  filter(
    source != target
  ) %>%
  tidyr::crossing(
    lag_months = 0:MAX_FOLLOWUP_MONTHS
  ) %>%
  mutate(
    target_date = trigger_date %m+%
      months(lag_months)
  ) %>%
  filter(
    target_date <= END_DATE
  ) %>%
  left_join(
    anomaly %>%
      select(
        nutrient,
        STAID,
        date,
        concentration,
        Z
      ) %>%
      rename(
        target = STAID,
        target_date = date,
        target_concentration = concentration,
        target_Z = Z
      ),
    by = c(
      "nutrient",
      "target",
      "target_date"
    )
  ) %>%
  mutate(
    target_observed = is.finite(target_Z),
    
    target_anomaly_15 = (
      target_observed &
        target_Z > Z_TRIGGER
    )
  )

save_csv(
  event_targets,
  "06_event_target_observation_matrix.csv"
)


# ============================================================
# 11. EVENT COVERAGE METRICS
# ============================================================

# Number of distinct candidate targets with
# at least one observation at the trigger
# month or within the subsequent 3 months.

event_coverage <- event_targets %>%
  group_by(
    event_id,
    nutrient,
    source,
    trigger_date
  ) %>%
  summarise(
    n_candidate_targets =
      n_distinct(target),
    
    n_observed_target_months =
      sum(target_observed),
    
    n_targets_observed =
      n_distinct(
        target[target_observed]
      ),
    
    n_targets_with_anomaly =
      n_distinct(
        target[
          target_anomaly_15
        ]
      ),
    
    n_targets_lag0 =
      n_distinct(
        target[
          lag_months == 0 &
            target_observed
        ]
      ),
    
    n_targets_lag1 =
      n_distinct(
        target[
          lag_months == 1 &
            target_observed
        ]
      ),
    
    n_targets_lag2 =
      n_distinct(
        target[
          lag_months == 2 &
            target_observed
        ]
      ),
    
    n_targets_lag3 =
      n_distinct(
        target[
          lag_months == 3 &
            target_observed
        ]
      ),
    
    .groups = "drop"
  ) %>%
  mutate(
    sufficient_observation =
      n_targets_observed >=
      MIN_TARGETS_OBSERVED
  )

save_csv(
  event_coverage,
  "07_event_observation_coverage.csv"
)


# ============================================================
# 12. EVENT CATALOG FOR FIGURE 5
# ============================================================

# Selection criteria here depend ONLY on
# source anomaly and observation availability.
#
# DO NOT select representative figures by
# whether F-Net happened to succeed.

event_catalog <- events %>%
  left_join(
    event_coverage %>%
      select(
        event_id,
        n_candidate_targets,
        n_targets_observed,
        n_observed_target_months,
        sufficient_observation,
        starts_with("n_targets_lag")
      ),
    by = "event_id"
  ) %>%
  arrange(
    nutrient,
    desc(sufficient_observation),
    desc(n_targets_observed),
    desc(source_Z)
  )

save_csv(
  event_catalog,
  "08_Fig5_candidate_event_catalog.csv"
)

cat("\nEligible events:\n")

print(
  event_catalog %>%
    group_by(nutrient) %>%
    summarise(
      n_events = n(),
      n_eligible = sum(
        sufficient_observation,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
)


# ============================================================
# 13. TRIGGER THRESHOLD SENSITIVITY
# ============================================================

threshold_stats <- lapply(
  Z_SENSITIVITY,
  
  function(th) {
    
    event_series %>%
      filter(
        period == "Evaluation",
        is.finite(Z),
        is.finite(Z_previous),
        Z > th,
        Z_previous <= Z_PREVIOUS_MAX
      ) %>%
      count(
        nutrient,
        name = "n_events"
      ) %>%
      mutate(
        Z_threshold = th
      )
  }
) %>%
  bind_rows()

save_csv(
  threshold_stats,
  "09_trigger_threshold_sensitivity.csv"
)


# ============================================================
# 14. SOURCE-STATION TRIGGER FREQUENCY
# ============================================================

source_summary <- events %>%
  group_by(
    nutrient,
    source
  ) %>%
  summarise(
    n_triggers = n(),
    
    max_source_Z = max(source_Z),
    
    mean_source_Z = mean(source_Z),
    
    .groups = "drop"
  ) %>%
  arrange(
    nutrient,
    desc(n_triggers)
  )

save_csv(
  source_summary,
  "10_source_station_trigger_frequency.csv"
)


# ============================================================
# 15. VISUAL QA — MONTHLY ANOMALY DISTRIBUTION
# ============================================================

p_z <- anomaly %>%
  filter(
    period == "Evaluation",
    is.finite(Z)
  ) %>%
  ggplot(
    aes(
      x = Z,
      fill = nutrient
    )
  ) +
  geom_histogram(
    bins = 60,
    alpha = 0.75,
    color = "white"
  ) +
  geom_vline(
    xintercept = Z_TRIGGER,
    linetype = "dashed",
    linewidth = 0.7,
    color = "#B84745"
  ) +
  facet_wrap(
    ~nutrient,
    ncol = 2,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = c(
      TN = "#287F9C",
      TP = "#CE7953"
    ),
    guide = "none"
  ) +
  labs(
    title = "Distribution of observed TN/TP anomalies",
    subtitle = paste0(
      "Evaluation period: 2020–2023; ",
      "reference climatology from training period"
    ),
    x = "Deseasonalized concentration anomaly (Z)",
    y = "Station-month observations"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(
    OUT, "Fig5S_anomaly_distribution.pdf"
  ),
  p_z,
  width = 10,
  height = 5,
  bg = "white"
)


# ============================================================
# 16. VISUAL QA — NUMBER OF ANOMALY TRIGGERS
# ============================================================

events_monthly <- events %>%
  count(
    nutrient,
    trigger_date,
    name = "n_triggers"
  )

p_events <- ggplot(
  events_monthly,
  aes(
    x = trigger_date,
    y = n_triggers,
    color = nutrient
  )
) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.3) +
  facet_wrap(
    ~nutrient,
    ncol = 1,
    scales = "free_y"
  ) +
  scale_color_manual(
    values = c(
      TN = "#287F9C",
      TP = "#CE7953"
    ),
    guide = "none"
  ) +
  labs(
    title = "Observed nutrient anomaly triggers",
    subtitle = "First-onset triggers, 2020–2023",
    x = "Date",
    y = "Number of triggering stations"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(
    OUT, "Fig5S_trigger_timeline.pdf"
  ),
  p_events,
  width = 11,
  height = 6,
  bg = "white"
)


# ============================================================
# 17. FINAL REPORT
# ============================================================

cat("\n========================================\n")
cat("FIGURE 5 — STEP 5.1 COMPLETED\n")
cat("========================================\n")

cat("\nTraining period:\n")
print(c(START_DATE, TRAIN_END))

cat("\nEvaluation period:\n")
print(c(TEST_START, END_DATE))

cat("\nObservation coverage:\n")

print(
  coverage %>%
    group_by(nutrient) %>%
    summarise(
      n_stations = n(),
      
      median_valid_fraction =
        median(valid_fraction),
      
      .groups = "drop"
    )
)

cat("\nActual anomaly triggers:\n")

print(
  events %>% count(nutrient)
)

cat("\nEligible events with target coverage:\n")

print(
  event_catalog %>%
    filter(sufficient_observation) %>%
    count(nutrient)
)

cat("\nThreshold sensitivity:\n")

print(
  threshold_stats,
  n = Inf
)

cat("\nTop eligible cases by observation coverage:\n")

print(
  event_catalog %>%
    filter(sufficient_observation) %>%
    group_by(nutrient) %>%
    slice_head(n = 10) %>%
    ungroup() %>%
    select(
      event_id,
      nutrient,
      source,
      trigger_date,
      source_Z,
      n_targets_observed
    ),
  n = Inf
)

cat("\nOutput directory:\n")
cat(OUT, "\n")

print(p_z)
print(p_events)