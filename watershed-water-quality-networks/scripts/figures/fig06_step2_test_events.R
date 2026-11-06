# ============================================================
# FIG6 STEP 2
# Independent anomaly event replay and evaluation dataset
# ============================================================

library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)

# ============================================================
# 1. PATHS
# ============================================================

source("config.R")
ROOT <- get_data_root()

TRAIN_ROOT <- file.path(
  ROOT, "Fig6_independent_training_2011_2020"
)

RAW_FILE <- file.path(
  ROOT, "NTN_2023_MonthlyLoadTable.csv"
)

OUT <- file.path(
  ROOT, "Fig6_independent_validation",
  "Step2_event_replay"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. SETTINGS
# ============================================================

TRAIN_START <- as.Date("2011-10-01")
TRAIN_END   <- as.Date("2020-12-01")

TEST_START <- as.Date("2021-01-01")
TEST_END   <- as.Date("2023-12-01")

# Require complete 0–3 month follow-up.
MAX_LAG <- 3L

Z_THRESHOLD <- 1.5

# At least six observations for each calendar-month baseline.
MIN_CLIM_N <- 6L

# Separate triggers from the same source by at least 4 months.
MIN_EVENT_GAP <- 4L

PCODES <- c(
  TN = "P00600",
  TP = "P00665"
)

# ============================================================
# 3. HELPERS
# ============================================================

id8 <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("\\.0+$", "", x)
  
  ii <- !is.na(x) & grepl("^[0-9]{1,8}$", x)
  
  x[ii] <- stringr::str_pad(
    x[ii], width = 8, side = "left", pad = "0"
  )
  
  x
}

write_out <- function(d, filename) {
  readr::write_csv(
    d,
    file.path(OUT, filename),
    na = ""
  )
}

# ============================================================
# 4. READ ORIGINAL MONTHLY CONCENTRATIONS
# ============================================================

if (!file.exists(RAW_FILE)) {
  stop("Raw concentration file missing: ", RAW_FILE)
}

raw <- readr::read_csv(
  RAW_FILE,
  col_types = readr::cols(
    STAID = readr::col_character(),
    PCODE = readr::col_character(),
    .default = readr::col_guess()
  ),
  show_col_types = FALSE
)

required <- c("STAID", "PCODE", "Year", "Month", "KConc")

if (length(setdiff(required, names(raw)))) {
  stop("Missing essential columns in raw data.")
}

raw <- raw %>%
  mutate(
    STAID = id8(STAID),
    PCODE = as.character(PCODE),
    Date = as.Date(sprintf(
      "%04d-%02d-01",
      as.integer(Year),
      as.integer(Month)
    )),
    KConc = as.numeric(KConc)
  ) %>%
  filter(
    Date >= TRAIN_START,
    Date <= TEST_END,
    is.finite(KConc),
    KConc >= 0
  )

if (anyDuplicated(raw[c("STAID", "PCODE", "Date")])) {
  stop("Duplicate station/nutrient/month records.")
}

# ============================================================
# 5. READ TRAINING NETWORK STATION IDS
# ============================================================

read_training_ids <- function(nu) {
  
  path <- file.path(
    TRAIN_ROOT,
    nu,
    "train_2011_2020_station_ids.csv"
  )
  
  if (!file.exists(path)) {
    stop("Training station ID file missing: ", path)
  }
  
  readr::read_csv(
    path,
    col_types = readr::cols(
      STAID = readr::col_character()
    ),
    show_col_types = FALSE
  ) %>%
    transmute(
      nutrient = nu,
      STAID = id8(STAID)
    )
}

training_ids <- bind_rows(
  read_training_ids("TN"),
  read_training_ids("TP")
)

# ============================================================
# 6. BUILD TRAINING-ONLY MONTHLY CLIMATOLOGY
# ============================================================

monthly <- bind_rows(
  lapply(names(PCODES), function(nu) {
    
    pcode_now <- unname(PCODES[[nu]])
    
    raw %>%
      filter(.data$PCODE == .env$pcode_now) %>%
      transmute(
        nutrient = nu,
        STAID,
        Date,
        KConc
      )
  })
) %>%
  semi_join(
    training_ids,
    by = c("nutrient", "STAID")
  ) %>%
  mutate(
    calendar_month = month(Date)
  )

climatology <- monthly %>%
  filter(
    Date >= TRAIN_START,
    Date <= TRAIN_END
  ) %>%
  group_by(
    nutrient,
    STAID,
    calendar_month
  ) %>%
  summarise(
    n_train = n(),
    mean_train = mean(KConc),
    sd_train = sd(KConc),
    .groups = "drop"
  ) %>%
  mutate(
    baseline_valid = (
      n_train >= MIN_CLIM_N &
        is.finite(sd_train) &
        sd_train > 0
    )
  )

write_out(
  climatology,
  "01_training_monthly_climatology.csv"
)

# ============================================================
# 7. APPLY FROZEN CLIMATOLOGY TO TEST DATA
# ============================================================

test_monthly <- monthly %>%
  filter(
    Date >= TEST_START,
    Date <= TEST_END
  ) %>%
  left_join(
    climatology,
    by = c(
      "nutrient",
      "STAID",
      "calendar_month"
    )
  ) %>%
  mutate(
    Z = if_else(
      baseline_valid %in% TRUE,
      (KConc - mean_train) / sd_train,
      NA_real_
    ),
    
    anomaly = if_else(
      is.finite(Z),
      Z > Z_THRESHOLD,
      NA
    )
  )

write_out(
  test_monthly,
  "02_test_monthly_anomalies.csv"
)

# ============================================================
# 8. IDENTIFY SOURCE ANOMALY EVENTS
# ============================================================

# Only evaluate triggers with a complete future 0–3 month horizon.
LAST_TRIGGER <- TEST_END %m-% months(MAX_LAG)

event_candidates <- test_monthly %>%
  filter(
    Date <= LAST_TRIGGER,
    anomaly %in% TRUE
  ) %>%
  arrange(nutrient, STAID, Date)

# Non-overlapping triggers from the same source.
# This is a pre-specified event declustering rule.

events_list <- list()

groups <- split(
  event_candidates,
  interaction(
    event_candidates$nutrient,
    event_candidates$STAID,
    drop = TRUE
  )
)

for (g in groups) {
  
  if (!nrow(g)) next
  
  g <- g[order(g$Date), , drop = FALSE]
  
  accepted <- logical(nrow(g))
  last_date <- as.Date(NA)
  
  for (i in seq_len(nrow(g))) {
    
    current_date <- g$Date[i]
    
    current_index <- (
      year(current_date) * 12L +
        month(current_date)
    )
    
    if (is.na(last_date)) {
      
      accepted[i] <- TRUE
      last_date <- current_date
      
    } else {
      
      previous_index <- (
        year(last_date) * 12L +
          month(last_date)
      )
      
      if (
        current_index - previous_index >= MIN_EVENT_GAP
      ) {
        
        accepted[i] <- TRUE
        last_date <- current_date
      }
    }
  }
  
  events_list[[length(events_list) + 1L]] <- g[
    accepted,
    ,
    drop = FALSE
  ]
}

events <- bind_rows(events_list) %>%
  arrange(nutrient, Date, STAID) %>%
  mutate(
    event_id = sprintf(
      "TEST_EV%05d",
      row_number()
    )
  ) %>%
  transmute(
    event_id,
    nutrient,
    source = STAID,
    trigger_date = Date,
    source_Z = Z,
    source_concentration = KConc,
    year = year(trigger_date),
    season = case_when(
      month(trigger_date) %in% c(12, 1, 2) ~ "DJF",
      month(trigger_date) %in% c(3, 4, 5) ~ "MAM",
      month(trigger_date) %in% c(6, 7, 8) ~ "JJA",
      TRUE ~ "SON"
    )
  )

write_out(
  events,
  "03_independent_test_events.csv"
)

# ============================================================
# 9. BUILD SHARED EVENT-TARGET CANDIDATES
# ============================================================

# All eligible training-network nodes are candidates,
# except the source itself.

candidates <- events %>%
  inner_join(
    training_ids %>%
      rename(target = STAID),
    by = "nutrient",
    relationship = "many-to-many"
  ) %>%
  filter(source != target) %>%
  select(
    event_id,
    nutrient,
    source,
    target,
    trigger_date
  )

# ============================================================
# 10. BUILD EVENT-TARGET-MONTH OBSERVATIONS
# ============================================================

event_months <- candidates %>%
  tidyr::crossing(
    lag_months = 0:MAX_LAG
  ) %>%
  mutate(
    target_date = trigger_date %m+%
      months(lag_months)
  ) %>%
  left_join(
    test_monthly %>%
      transmute(
        nutrient,
        target = STAID,
        target_date = Date,
        target_Z = Z,
        target_anomaly = anomaly,
        target_concentration = KConc
      ),
    by = c(
      "nutrient",
      "target",
      "target_date"
    )
  ) %>%
  mutate(
    target_observed = is.finite(target_Z)
  )

if (anyDuplicated(
  event_months[c(
    "event_id",
    "target",
    "lag_months"
  )]
)) {
  stop("Duplicate event-target-month after join.")
}

write_out(
  event_months,
  "04_event_target_month_matrix.csv"
)

# ============================================================
# 11. BUILD SHARED EVALUATION TARGET TABLE
# ============================================================

target_evaluation <- event_months %>%
  group_by(
    event_id,
    nutrient,
    source,
    target,
    trigger_date
  ) %>%
  summarise(
    n_observed_months = sum(
      target_observed
    ),
    
    n_required_months = MAX_LAG + 1L,
    
    n_anomaly_months = sum(
      target_anomaly %in% TRUE
    ),
    
    fully_observed = (
      n_observed_months ==
        n_required_months
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    # Primary evaluation uses fully observed targets.
    evaluable = fully_observed,
    
    hit_03 = case_when(
      evaluable & n_anomaly_months > 0 ~ TRUE,
      evaluable & n_anomaly_months == 0 ~ FALSE,
      TRUE ~ NA
    )
  )

write_out(
  target_evaluation,
  "05_common_event_target_evaluation.csv"
)

# ============================================================
# 12. QA SUMMARIES
# ============================================================

event_qa <- events %>%
  group_by(nutrient, year, season) %>%
  summarise(
    n_events = n(),
    n_sources = n_distinct(source),
    .groups = "drop"
  )

target_qa <- target_evaluation %>%
  group_by(nutrient) %>%
  summarise(
    n_events = n_distinct(event_id),
    n_candidates = n(),
    n_evaluable = sum(evaluable),
    n_not_evaluable = sum(!evaluable),
    n_positive = sum(hit_03 %in% TRUE),
    n_negative = sum(hit_03 %in% FALSE),
    evaluable_fraction = mean(evaluable),
    .groups = "drop"
  )

event_level_qa <- target_evaluation %>%
  group_by(event_id, nutrient) %>%
  summarise(
    n_candidates = n(),
    n_evaluable = sum(evaluable),
    n_positive = sum(hit_03 %in% TRUE),
    .groups = "drop"
  )

write_out(
  event_qa,
  "06_test_event_QA.csv"
)

write_out(
  target_qa,
  "07_common_target_QA.csv"
)

write_out(
  event_level_qa,
  "08_event_level_coverage_QA.csv"
)

cat("\n====================================\n")
cat("FIG6 STEP 2 COMPLETED\n")
cat("====================================\n")

cat("\nIndependent test events:\n")
print(event_qa, n = Inf)

cat("\nShared evaluation target summary:\n")
print(target_qa, n = Inf)

cat("\nOutput directory:\n")
cat(OUT, "\n")