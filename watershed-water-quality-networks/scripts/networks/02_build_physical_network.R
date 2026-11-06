source("config.R")

# ============================================================
# Chesapeake Bay NTN Physical Propagation Network
# P-Net: D8 downstream tracing + FAC snapping + Daily Q
#
# Memory-safe version for very large NHDPlus rasters
# Main period: 2011-10 to 2023-09
#
# Outputs:
# 1. Physical upstream-downstream edges
# 2. Directed reachability matrix
# 3. Station snapping QA
# 4. River tracing QA
# 5. Daily discharge summaries
# 6. Network metrics and GraphML
# 7. Watershed + rivers + physical network maps
# ============================================================


# ============================================================
# 0. Packages
# ============================================================

packages <- c(
  "terra",
  "sf",
  "readr",
  "dplyr",
  "tidyr",
  "stringr",
  "lubridate",
  "ggplot2",
  "igraph",
  "tibble"
)

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

invisible(
  lapply(packages, library, character.only = TRUE)
)

terraOptions(
  memfrac = 0.35,
  progress = 1
)


# ============================================================
# 1. Paths
# ============================================================

ROOT <- get_data_root()

FDRFAC_DIR <- file.path(
  ROOT,
  "NHDPlusMA/NHDPlus02/NHDPlusFdrFac02a"
)

FDRNULL_DIR <- file.path(
  ROOT,
  "NHDPlusMA/NHDPlus02/NHDPlusFdrNull02a"
)

RIVER_SHP <- file.path(
  ROOT,
  "水系/水系.shp"
)

BASIN_SHP <- file.path(
  ROOT,
  "ntnwatersheds/ntnwatersheds.shp"
)

STATION_FILE <- file.path(
  ROOT,
  "NTN-station.txt"
)

NTN_FILE <- file.path(
  ROOT,
  "NTN_2023_MonthlyLoadTable.csv"
)

Q_FILE <- file.path(
  ROOT,
  "USGS_daily_Q.csv"
)

OUT <- file.path(
  ROOT,
  "PNet_2011_2023_BlockTracing"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. Main settings
# ============================================================

START_DATE <- as.Date("2011-10-01")
END_DATE   <- as.Date("2023-09-01")

PCODE <- "P00665"   # P00600=TN; P00665 = TP
VALUE <- "KConc"

MAIN_ONLY <- TRUE

# Station snapping
SNAP_RADIUS_M <- 1000

# Candidate points must be close enough to original station
MAX_ACCEPTED_SNAP_M <- 1000

# FAC-based snap score
# A larger value puts more emphasis on drainage accumulation
FAC_WEIGHT <- 0.20

# Block cache
BLOCK_SIZE <- 256L
MAX_CACHE_BLOCKS <- 12L

# Maximum traced D8 steps for each station
MAX_TRACE_STEPS <- 2000000L

# NHDPlus FDR is assumed to use ArcGIS D8
# Confirm by checking actual raster documentation
D8_ENCODING <- "ESRI"

# Flow velocity scenarios (m/s)
# These are illustrative, NOT discharge-derived velocities
VELOCITIES <- c(
  slow = 0.2,
  medium = 0.5,
  fast = 1.0
)

# Drawing settings
SHOW_ALL_STATIONS <- TRUE


# ============================================================
# 3. Helpers
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


save_table <- function(x, name) {
  
  readr::write_csv(
    x,
    file.path(OUT, name)
  )
}


# ============================================================
# 4. Find NHDPlus FDR and FAC rasters
# ============================================================

find_nhd_grid <- function(base_dir, grid_name) {
  
  all_files <- list.files(
    base_dir,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE
  )
  
  # ESRI binary grid directories
  hdr <- all_files[
    grepl(
      "hdr\\.adf$",
      all_files,
      ignore.case = TRUE
    )
  ]
  
  grid_dirs <- unique(dirname(hdr))
  
  grid_dirs <- grid_dirs[
    tolower(basename(grid_dirs)) ==
      tolower(grid_name)
  ]
  
  # TIFF / IMG alternatives
  other <- all_files[
    grepl(
      "\\.(tif|tiff|img)$",
      all_files,
      ignore.case = TRUE
    )
  ]
  
  other <- other[
    tolower(
      tools::file_path_sans_ext(
        basename(other)
      )
    ) == tolower(grid_name)
  ]
  
  unique(c(grid_dirs, other))
}


fdr_paths <- find_nhd_grid(
  FDRFAC_DIR,
  "fdr"
)

fac_paths <- find_nhd_grid(
  FDRFAC_DIR,
  "fac"
)

if (
  length(fdr_paths) != 1 ||
  length(fac_paths) != 1
) {
  
  cat("\nFDR candidates:\n")
  print(fdr_paths)
  
  cat("\nFAC candidates:\n")
  print(fac_paths)
  
  stop(
    "Could not uniquely identify FDR/FAC. ",
    "Specify their exact paths."
  )
}

fdr <- terra::rast(fdr_paths[1])
fac <- terra::rast(fac_paths[1])

if (
  !terra::compareGeom(
    fdr,
    fac,
    stopOnError = FALSE
  )
) {
  stop("FDR and FAC raster geometries differ.")
}

if (terra::is.lonlat(fdr)) {
  stop(
    "FDR uses geographic coordinates. ",
    "This script requires a projected metric grid."
  )
}

if (
  is.na(terra::crs(fdr)) ||
  !nzchar(terra::crs(fdr))
) {
  stop("FDR has no defined CRS.")
}

if (D8_ENCODING != "ESRI") {
  stop("Only ESRI D8 encoding implemented.")
}

cat("\n=============== RASTER INFO ===============\n")

print(fdr)
print(fac)

cat(
  "Raster cells:",
  format(terra::ncell(fdr), big.mark = ","),
  "\n"
)

cat(
  "Resolution:",
  terra::res(fdr),
  "\n"
)

cat("===========================================\n")


# ============================================================
# 5. Read stations
# ============================================================

stations_all <- readr::read_csv(
  STATION_FILE,
  col_types = cols(
    SITE_NO = col_character(),
    .default = col_guess()
  ),
  show_col_types = FALSE
) %>%
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

stations <- stations_all


# ============================================================
# 6. Select complete stations for 2011-2023
# ============================================================

if (MAIN_ONLY) {
  
  ntn <- readr::read_csv(
    NTN_FILE,
    col_types = cols(
      STAID = col_character(),
      PCODE = col_character(),
      .default = col_guess()
    ),
    show_col_types = FALSE
  )
  
  months <- seq(
    START_DATE,
    END_DATE,
    by = "month"
  )
  
  eligible <- ntn %>%
    filter(
      PCODE == !!PCODE
    ) %>%
    mutate(
      STAID = norm_id(STAID),
      Date = as.Date(
        sprintf(
          "%04d-%02d-01",
          Year,
          Month
        )
      ),
      val = as.numeric(.data[[VALUE]])
    ) %>%
    filter(
      Date %in% months,
      is.finite(val)
    ) %>%
    group_by(STAID) %>%
    summarise(
      n_months = n_distinct(Date),
      .groups = "drop"
    ) %>%
    filter(
      n_months == length(months)
    )
  
  stations <- stations %>%
    semi_join(
      eligible,
      by = "STAID"
    )
}

if (nrow(stations) < 2) {
  stop("Insufficient eligible NTN stations.")
}

cat(
  "\nSelected NTN stations:",
  nrow(stations),
  "\n"
)

save_table(
  stations,
  "PNet_input_stations.csv"
)


# ============================================================
# 7. Transform stations into FDR coordinate system
# ============================================================

stations_sf <- st_as_sf(
  stations,
  coords = c("lon", "lat"),
  crs = 4326,
  remove = FALSE
)

stations_proj <- st_transform(
  stations_sf,
  terra::crs(fdr)
)

coords <- st_coordinates(stations_proj)

stations$x <- coords[, 1]
stations$y <- coords[, 2]


# ============================================================
# 8. Snap stations using local raster windows
#
# This function DOES NOT load full FDR/FAC.
# ============================================================

snap_station <- function(
    x, y, fdr, fac,
    radius_m = 1000,
    fac_weight = 0.20
) {
  
  search_extent <- terra::ext(
    x - radius_m,
    x + radius_m,
    y - radius_m,
    y + radius_m
  )
  
  fdr_local <- tryCatch(
    terra::crop(
      fdr,
      search_extent,
      snap = "near"
    ),
    error = function(e) NULL
  )
  
  if (is.null(fdr_local)) {
    return(NULL)
  }
  
  fac_local <- terra::crop(
    fac,
    terra::ext(fdr_local),
    snap = "near"
  )
  
  if (
    !terra::compareGeom(
      fdr_local,
      fac_local,
      stopOnError = FALSE
    )
  ) {
    stop("Local FDR/FAC geometry mismatch.")
  }
  
  fd <- terra::values(
    fdr_local,
    mat = FALSE
  )
  
  fa <- terra::values(
    fac_local,
    mat = FALSE
  )
  
  valid_codes <- c(
    1, 2, 4, 8,
    16, 32, 64, 128
  )
  
  valid <- which(
    !is.na(fd) &
      fd %in% valid_codes &
      is.finite(fa) &
      fa > 0
  )
  
  if (!length(valid)) {
    return(NULL)
  }
  
  local_xy <- terra::xyFromCell(
    fdr_local,
    valid
  )
  
  distance <- sqrt(
    (local_xy[, 1] - x)^2 +
      (local_xy[, 2] - y)^2
  )
  
  inside <- which(
    distance <= radius_m
  )
  
  if (!length(inside)) {
    return(NULL)
  }
  
  valid <- valid[inside]
  
  local_xy <- local_xy[
    inside,
    ,
    drop = FALSE
  ]
  
  distance <- distance[inside]
  
  accumulation <- fa[valid]
  
  logfac <- log1p(accumulation)
  
  denom <- max(logfac)
  
  fac_score <- if (denom > 0) {
    1 - logfac / denom
  } else {
    rep(0, length(logfac))
  }
  
  score <- distance / radius_m +
    fac_weight * fac_score
  
  best <- which.min(score)
  
  xy_best <- matrix(
    local_xy[best, ],
    nrow = 1
  )
  
  global_cell <- terra::cellFromXY(
    fdr,
    xy_best
  )
  
  tibble(
    cell = as.numeric(global_cell),
    snapped_x = local_xy[best, 1],
    snapped_y = local_xy[best, 2],
    snap_distance_m = distance[best],
    FAC = accumulation[best],
    score = score[best]
  )
}


# ============================================================
# 9. Snap all stations
# ============================================================

cat("\nSnapping stations...\n")

snap_results <- vector(
  "list",
  nrow(stations)
)

for (i in seq_len(nrow(stations))) {
  
  cat(
    sprintf(
      "\rSnapping %d / %d",
      i,
      nrow(stations)
    )
  )
  
  res <- snap_station(
    stations$x[i],
    stations$y[i],
    fdr,
    fac,
    radius_m = SNAP_RADIUS_M,
    fac_weight = FAC_WEIGHT
  )
  
  if (!is.null(res)) {
    
    snap_results[[i]] <- bind_cols(
      stations[i, ],
      res
    )
    
  } else {
    
    snap_results[[i]] <- stations[i, ] %>%
      mutate(
        cell = NA_real_,
        snapped_x = NA_real_,
        snapped_y = NA_real_,
        snap_distance_m = NA_real_,
        FAC = NA_real_,
        score = NA_real_
      )
  }
}

cat("\n")

snap_QA <- bind_rows(snap_results) %>%
  mutate(
    snap_status = case_when(
      is.na(cell) ~ "failed",
      snap_distance_m > MAX_ACCEPTED_SNAP_M ~
        "distance_exceeded",
      TRUE ~ "snapped"
    )
  )

save_table(
  snap_QA,
  "PNet_station_snap_QA.csv"
)

valid_stations <- snap_QA %>%
  filter(
    snap_status == "snapped"
  )

cat(
  "\nSuccessfully snapped:",
  nrow(valid_stations),
  "/",
  nrow(stations),
  "\n"
)

if (nrow(valid_stations) < 2) {
  stop("Fewer than two stations successfully snapped.")
}

# Detect stations sharing one raster cell
shared_cells <- valid_stations %>%
  group_by(cell) %>%
  filter(n() > 1) %>%
  ungroup()

save_table(
  shared_cells,
  "PNet_shared_snap_cells_QA.csv"
)


# ============================================================
# 10. Block-cached raster reader
#
# Only 256 x 256 cells read per block.
# Maximum 12 blocks retained in RAM.
# ============================================================

make_block_reader <- function(
    raster,
    block_size = 256L,
    max_blocks = 12L
) {
  
  cache <- new.env(
    hash = TRUE,
    parent = emptyenv()
  )
  
  order_keys <- character()
  
  nrows <- terra::nrow(raster)
  ncols <- terra::ncol(raster)
  
  get_value <- function(cell) {
    
    if (
      is.na(cell) ||
      cell < 1 ||
      cell > terra::ncell(raster)
    ) {
      return(NA_real_)
    }
    
    row <- floor((cell - 1) / ncols) + 1
    col <- ((cell - 1) %% ncols) + 1
    
    block_r <- floor(
      (row - 1) / block_size
    )
    
    block_c <- floor(
      (col - 1) / block_size
    )
    
    key <- paste(
      block_r,
      block_c,
      sep = "_"
    )
    
    if (!exists(
      key,
      envir = cache,
      inherits = FALSE
    )) {
      
      r1 <- block_r * block_size + 1
      c1 <- block_c * block_size + 1
      
      r2 <- min(
        r1 + block_size - 1,
        nrows
      )
      
      c2 <- min(
        c1 + block_size - 1,
        ncols
      )
      
      corners <- terra::xyFromCell(
        raster,
        c(
          terra::cellFromRowCol(
            raster, r1, c1
          ),
          terra::cellFromRowCol(
            raster, r2, c2
          )
        )
      )
      
      dx <- terra::res(raster)[1] / 2
      dy <- terra::res(raster)[2] / 2
      
      block_ext <- terra::ext(
        min(corners[, 1]) - dx,
        max(corners[, 1]) + dx,
        min(corners[, 2]) - dy,
        max(corners[, 2]) + dy
      )
      
      block <- terra::crop(
        raster,
        block_ext,
        snap = "near"
      )
      
      block_values <- terra::values(
        block,
        mat = FALSE
      )
      
      assign(
        key,
        list(
          raster = block,
          values = block_values
        ),
        envir = cache
      )
      
      order_keys <<- c(
        order_keys,
        key
      )
      
      if (
        length(order_keys) >
        max_blocks
      ) {
        
        oldest <- order_keys[1]
        
        rm(
          list = oldest,
          envir = cache
        )
        
        order_keys <<- order_keys[-1]
      }
      
    } else {
      
      # Keep recently used blocks in cache
      order_keys <<- c(
        order_keys[order_keys != key],
        key
      )
    }
    
    obj <- get(
      key,
      envir = cache,
      inherits = FALSE
    )
    
    xy <- terra::xyFromCell(
      raster,
      cell
    )
    
    local_cell <- terra::cellFromXY(
      obj$raster,
      xy
    )
    
    if (is.na(local_cell)) {
      return(NA_real_)
    }
    
    as.numeric(
      obj$values[local_cell]
    )
  }
  
  get_value
}


cat("\nCreating cached FDR reader...\n")

get_fdr <- make_block_reader(
  fdr,
  block_size = BLOCK_SIZE,
  max_blocks = MAX_CACHE_BLOCKS
)


# ============================================================
# 11. D8 direction rules
# ============================================================

# ArcGIS D8:
# 1   = East
# 2   = Southeast
# 4   = South
# 8   = Southwest
# 16  = West
# 32  = Northwest
# 64  = North
# 128 = Northeast

D8 <- data.frame(
  code = c(
    1, 2, 4, 8,
    16, 32, 64, 128
  ),
  dr = c(
    0, 1, 1, 1,
    0, -1, -1, -1
  ),
  dc = c(
    1, 1, 0, -1,
    -1, -1, 0, 1
  )
)


# ============================================================
# 12. D8 downstream function
# ============================================================

nc <- terra::ncol(fdr)
nr <- terra::nrow(fdr)

cell_res <- terra::res(fdr)

step_x <- cell_res[1]
step_y <- cell_res[2]

next_downstream <- function(cell) {
  
  code <- get_fdr(cell)
  
  if (
    is.na(code) ||
    !(code %in% D8$code)
  ) {
    return(NULL)
  }
  
  ii <- match(
    code,
    D8$code
  )
  
  dr <- D8$dr[ii]
  dc <- D8$dc[ii]
  
  row <- floor(
    (cell - 1) / nc
  ) + 1
  
  col <- (
    (cell - 1) %% nc
  ) + 1
  
  row2 <- row + dr
  col2 <- col + dc
  
  if (
    row2 < 1 ||
    row2 > nr ||
    col2 < 1 ||
    col2 > nc
  ) {
    return(NULL)
  }
  
  next_cell <- (
    (row2 - 1) * nc + col2
  )
  
  distance_m <- sqrt(
    (dc * step_x)^2 +
      (dr * step_y)^2
  )
  
  list(
    cell = next_cell,
    distance_m = distance_m
  )
}


# ============================================================
# 13. Build station lookup
# ============================================================

# A raster cell can correspond to multiple stations.
# All such cases are reported in QA.

station_cell_lookup <- split(
  valid_stations$STAID,
  as.character(valid_stations$cell)
)

# Store all traced edge records
all_edges <- vector(
  "list",
  nrow(valid_stations)
)

trace_QA <- vector(
  "list",
  nrow(valid_stations)
)


# ============================================================
# 14. Trace downstream from every station
# ============================================================

cat("\n====================================\n")
cat("STARTING D8 DOWNSTREAM TRACING\n")
cat("====================================\n")

for (i in seq_len(nrow(valid_stations))) {
  
  station_id <- valid_stations$STAID[i]
  
  current_cell <- valid_stations$cell[i]
  
  cumulative_distance <- 0
  
  visited <- new.env(
    hash = TRUE,
    parent = emptyenv()
  )
  
  edge_i <- list()
  
  n_steps <- 0L
  status <- "max_steps"
  
  cat(
    "\nTracing",
    i,
    "/",
    nrow(valid_stations),
    ":",
    station_id,
    "\n"
  )
  
  for (step in seq_len(MAX_TRACE_STEPS)) {
    
    key <- as.character(current_cell)
    
    if (exists(
      key,
      envir = visited,
      inherits = FALSE
    )) {
      
      status <- "D8_loop"
      break
    }
    
    assign(
      key,
      TRUE,
      envir = visited
    )
    
    next_info <- next_downstream(
      current_cell
    )
    
    if (is.null(next_info)) {
      
      status <- "sink_boundary_or_nodata"
      break
    }
    
    current_cell <- next_info$cell
    
    cumulative_distance <-
      cumulative_distance +
      next_info$distance_m
    
    n_steps <- n_steps + 1L
    
    downstream_ids <-
      station_cell_lookup[[
        as.character(current_cell)
      ]]
    
    if (!is.null(downstream_ids)) {
      
      for (downstream_id in downstream_ids) {
        
        if (
          downstream_id != station_id
        ) {
          
          edge_i[[length(edge_i) + 1L]] <-
            tibble(
              from = station_id,
              to = downstream_id,
              river_distance_km =
                cumulative_distance / 1000,
              D8_steps = n_steps
            )
        }
      }
    }
    
    if (step %% 100000L == 0L) {
      
      cat(
        "  Steps:",
        step,
        "Distance (km):",
        round(
          cumulative_distance / 1000,
          1
        ),
        "\n"
      )
    }
  }
  
  all_edges[[i]] <- bind_rows(edge_i)
  
  trace_QA[[i]] <- tibble(
    STAID = station_id,
    start_cell = valid_stations$cell[i],
    n_steps = n_steps,
    traced_distance_km =
      cumulative_distance / 1000,
    termination = status,
    downstream_stations_found =
      nrow(all_edges[[i]])
  )
}


# ============================================================
# 15. Save tracing QA
# ============================================================

trace_QA_df <- bind_rows(
  trace_QA
)

save_table(
  trace_QA_df,
  "PNet_traces_QA.csv"
)


# ============================================================
# 16. Construct directed P-Net edges
# ============================================================

edges <- bind_rows(
  all_edges
)

if (nrow(edges) > 0) {
  
  edges <- edges %>%
    filter(
      from != to,
      is.finite(river_distance_km),
      river_distance_km > 0
    ) %>%
    group_by(
      from,
      to
    ) %>%
    slice_min(
      river_distance_km,
      n = 1,
      with_ties = FALSE
    ) %>%
    ungroup()
  
} else {
  
  edges <- tibble(
    from = character(),
    to = character(),
    river_distance_km = numeric(),
    D8_steps = integer()
  )
}


# ============================================================
# 17. Add illustrative travel-time scenarios
# ============================================================

for (nm in names(VELOCITIES)) {
  
  velocity <- VELOCITIES[[nm]]
  
  column <- paste0(
    "travel_days_",
    nm
  )
  
  edges[[column]] <-
    edges$river_distance_km *
    1000 /
    (velocity * 86400)
}

save_table(
  edges,
  "PNet_directed_reachable_pairs.csv"
)


# ============================================================
# 18. Construct directed reachability matrix
# ============================================================

station_ids <- valid_stations$STAID

reach_matrix <- matrix(
  0L,
  nrow = length(station_ids),
  ncol = length(station_ids),
  dimnames = list(
    station_ids,
    station_ids
  )
)

if (nrow(edges) > 0) {
  
  for (i in seq_len(nrow(edges))) {
    
    reach_matrix[
      edges$from[i],
      edges$to[i]
    ] <- 1L
  }
}

# A missing connection is NOT necessarily verified absence:
# raster boundary/incorrect snap may prevent finding a path.

reach_df <- tibble(
  STAID = rownames(reach_matrix)
) %>%
  bind_cols(
    as_tibble(
      reach_matrix,
      .name_repair = "minimal"
    )
  )

save_table(
  reach_df,
  "PNet_reachability_matrix.csv"
)


# ============================================================
# 19. Read and summarize USGS daily Q
# ============================================================

if (file.exists(Q_FILE)) {
  
  cat("\nReading daily Q...\n")
  
  qraw <- readr::read_csv(
    Q_FILE,
    show_col_types = FALSE
  )
  
  writeLines(
    names(qraw),
    file.path(
      OUT,
      "PNet_Q_original_field_names.txt"
    )
  )
  
  identify_col <- function(
    df,
    candidates
  ) {
    
    idx <- match(
      tolower(candidates),
      tolower(names(df))
    )
    
    idx <- idx[!is.na(idx)]
    
    if (!length(idx)) {
      return(NA_character_)
    }
    
    names(df)[idx[1]]
  }
  
  id_col <- identify_col(
    qraw,
    c(
      "STAID",
      "SITE_NO",
      "site_no",
      "site",
      "site_id",
      "station",
      "USGS"
    )
  )
  
  date_col <- identify_col(
    qraw,
    c(
      "Date",
      "DATE",
      "datetime",
      "time"
    )
  )
  
  q_col <- identify_col(
    qraw,
    c(
      "Q",
      "Flow",
      "discharge",
      "Q_cms",
      "flow_cms",
      "00060_00003",
      "X_00060_00003"
    )
  )
  
  if (
    !is.na(id_col) &&
    !is.na(date_col) &&
    !is.na(q_col)
  ) {
    
    dates_raw <- as.character(
      qraw[[date_col]]
    )
    
    parsed_dates <- suppressWarnings(
      lubridate::ymd(
        dates_raw
      )
    )
    
    if (all(is.na(parsed_dates))) {
      
      parsed_dates <- suppressWarnings(
        lubridate::mdy(
          dates_raw
        )
      )
    }
    
    qdat <- tibble(
      STAID = norm_id(
        qraw[[id_col]]
      ),
      Date = as.Date(parsed_dates),
      Q_original_units = suppressWarnings(
        as.numeric(
          qraw[[q_col]]
        )
      )
    )
    
    qstats <- qdat %>%
      filter(
        Date >= START_DATE,
        Date <= END_DATE,
        is.finite(Q_original_units),
        Q_original_units >= 0
      ) %>%
      group_by(STAID) %>%
      summarise(
        Q_mean = mean(
          Q_original_units
        ),
        Q_median = median(
          Q_original_units
        ),
        Q_p10 = as.numeric(
          quantile(
            Q_original_units,
            0.10
          )
        ),
        Q_p90 = as.numeric(
          quantile(
            Q_original_units,
            0.90
          )
        ),
        n_valid_days = n(),
        .groups = "drop"
      )
    
    save_table(
      qstats,
      "PNet_daily_Q_summary.csv"
    )
    
    edges_Q <- edges %>%
      left_join(
        qstats %>%
          select(
            from = STAID,
            Q_upstream_median = Q_median
          ),
        by = "from"
      ) %>%
      left_join(
        qstats %>%
          select(
            to = STAID,
            Q_downstream_median = Q_median
          ),
        by = "to"
      )
    
    save_table(
      edges_Q,
      "PNet_edges_with_daily_Q.csv"
    )
    
  } else {
    
    warning(
      "Daily Q fields were not recognized. ",
      "See PNet_Q_original_field_names.txt."
    )
  }
}


# ============================================================
# 20. Build igraph physical network
# ============================================================

vertices <- tibble(
  name = valid_stations$STAID
)

g <- igraph::graph_from_data_frame(
  edges %>%
    select(
      from,
      to,
      river_distance_km
    ),
  directed = TRUE,
  vertices = vertices
)

igraph::write_graph(
  g,
  file.path(
    OUT,
    "PNet_directed.graphml"
  ),
  format = "graphml"
)


# ============================================================
# 21. Network node metrics
# ============================================================

node_metrics <- tibble(
  STAID = igraph::V(g)$name,
  in_degree = as.numeric(
    igraph::degree(
      g,
      mode = "in"
    )
  ),
  out_degree = as.numeric(
    igraph::degree(
      g,
      mode = "out"
    )
  )
) %>%
  left_join(
    valid_stations %>%
      select(
        STAID,
        lon,
        lat,
        FAC,
        snap_distance_m
      ),
    by = "STAID"
  )

save_table(
  node_metrics,
  "PNet_node_metrics.csv"
)


# ============================================================
# 22. Read watershed and river network for mapping
# ============================================================

basin <- st_read(
  BASIN_SHP,
  quiet = TRUE
)

rivers <- st_read(
  RIVER_SHP,
  quiet = TRUE
)

if (is.na(st_crs(basin))) {
  stop(
    "Watershed CRS missing. ",
    "Assign its verified original CRS."
  )
}

if (is.na(st_crs(rivers))) {
  
  bb <- st_bbox(rivers)
  
  geographic <- all(
    is.finite(bb)
  ) &&
    abs(bb["xmin"]) <= 180 &&
    abs(bb["xmax"]) <= 180 &&
    abs(bb["ymin"]) <= 90 &&
    abs(bb["ymax"]) <= 90
  
  if (geographic) {
    
    st_crs(rivers) <- 4326
    
    warning(
      "River CRS provisionally assigned EPSG:4326. ",
      "Verify the original coordinate system."
    )
    
  } else {
    
    stop(
      "River CRS missing and coordinates are not ",
      "recognizable as longitude/latitude."
    )
  }
}

basin <- st_transform(
  basin,
  4326
)

rivers <- st_transform(
  rivers,
  4326
)


# ============================================================
# 23. Station geometry
# ============================================================

pts <- st_as_sf(
  snap_QA,
  coords = c("lon", "lat"),
  crs = 4326,
  remove = FALSE
) %>%
  left_join(
    node_metrics %>%
      select(
        STAID,
        in_degree,
        out_degree
      ),
    by = "STAID"
  )

pts$group <- ifelse(
  pts$snap_status == "snapped",
  "P-Net nodes",
  "Unmatched nodes"
)


# ============================================================
# 24. Build edge geometry for visualization
#
# Straight connecting lines are schematic.
# The distance in edges is calculated along D8 flow paths.
# ============================================================

edge_lines <- list()

if (nrow(edges) > 0) {
  
  station_pos <- valid_stations %>%
    select(
      STAID,
      lon,
      lat
    )
  
  edge_geo <- edges %>%
    left_join(
      station_pos %>%
        rename(
          from = STAID,
          lon_from = lon,
          lat_from = lat
        ),
      by = "from"
    ) %>%
    left_join(
      station_pos %>%
        rename(
          to = STAID,
          lon_to = lon,
          lat_to = lat
        ),
      by = "to"
    )
  
  edge_lines <- lapply(
    seq_len(nrow(edge_geo)),
    function(i) {
      
      st_linestring(
        matrix(
          c(
            edge_geo$lon_from[i],
            edge_geo$lat_from[i],
            edge_geo$lon_to[i],
            edge_geo$lat_to[i]
          ),
          ncol = 2,
          byrow = TRUE
        )
      )
    }
  )
  
  edges_sf <- st_sf(
    edge_geo,
    geometry = st_sfc(
      edge_lines,
      crs = 4326
    )
  )
  
  st_write(
    edges_sf,
    file.path(
      OUT,
      "PNet_schematic_edges.gpkg"
    ),
    delete_dsn = TRUE,
    quiet = TRUE
  )
  
} else {
  
  edges_sf <- NULL
}


# ============================================================
# 25. P-Net map with main river network
# ============================================================

p <- ggplot() +
  
  geom_sf(
    data = basin,
    fill = "#F7F7F7",
    color = "#A9A9A9",
    linewidth = 0.22
  ) +
  
  geom_sf(
    data = rivers,
    color = "#7BB2D3",
    linewidth = 0.30,
    alpha = 0.85
  )

if (
  !is.null(edges_sf) &&
  nrow(edges_sf) > 0
) {
  
  p <- p +
    geom_sf(
      data = edges_sf,
      color = "#D66A5E",
      linewidth = 0.25,
      linetype = "dashed",
      alpha = 0.45
    )
}

p <- p +
  
  geom_sf(
    data = pts,
    aes(
      color = group
    ),
    size = 2.0,
    alpha = 0.95
  ) +
  
  scale_color_manual(
    values = c(
      "P-Net nodes" = "#B63F38",
      "Unmatched nodes" = "grey50"
    ),
    name = NULL
  ) +
  
  coord_sf(
    datum = NA
  ) +
  
  labs(
    title = "Physical propagation network",
    subtitle = paste0(
      "NHDPlus D8 flow routing | ",
      nrow(valid_stations),
      " matched NTN stations"
    ),
    x = NULL,
    y = NULL
  ) +
  
  theme_bw(
    base_size = 11
  ) +
  
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(
      face = "bold"
    )
  )

ggsave(
  file.path(
    OUT,
    "PNet_station_river_map.pdf"
  ),
  p,
  width = 9,
  height = 8
)

ggsave(
  file.path(
    OUT,
    "PNet_station_river_map.png"
  ),
  p,
  width = 9,
  height = 8,
  dpi = 350
)


# ============================================================
# 26. Additional map: physical output degree
# ============================================================

p_degree <- ggplot() +
  
  geom_sf(
    data = basin,
    fill = "grey97",
    color = "grey75",
    linewidth = 0.25
  ) +
  
  geom_sf(
    data = rivers,
    color = "#85B9D5",
    linewidth = 0.28
  ) +
  
  geom_sf(
    data = pts %>%
      filter(
        snap_status == "snapped"
      ),
    aes(
      size = out_degree
    ),
    color = "#B64D44",
    alpha = 0.85
  ) +
  
  scale_size_continuous(
    range = c(1.4, 5)
  ) +
  
  coord_sf(
    datum = NA
  ) +
  
  labs(
    title = "Physical downstream reachability",
    size = "Downstream stations",
    x = NULL,
    y = NULL
  ) +
  
  theme_bw(
    base_size = 11
  ) +
  
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  file.path(
    OUT,
    "PNet_out_degree_map.pdf"
  ),
  p_degree,
  width = 9,
  height = 8
)

ggsave(
  file.path(
    OUT,
    "PNet_out_degree_map.png"
  ),
  p_degree,
  width = 9,
  height = 8,
  dpi = 350
)


# ============================================================
# 27. Overall QA
# ============================================================

summary <- tibble(
  input_stations = nrow(stations),
  snapped_stations = nrow(valid_stations),
  snap_failed = sum(
    snap_QA$snap_status != "snapped"
  ),
  shared_snap_cells = n_distinct(
    shared_cells$cell
  ),
  directed_physical_edges = nrow(edges),
  possible_directed_pairs =
    nrow(valid_stations) *
    (nrow(valid_stations) - 1),
  density = if (
    nrow(valid_stations) > 1
  ) {
    nrow(edges) /
      (
        nrow(valid_stations) *
          (nrow(valid_stations) - 1)
      )
  } else {
    NA_real_
  },
  n_trace_max_steps = sum(
    trace_QA_df$termination == "max_steps"
  ),
  n_trace_loops = sum(
    trace_QA_df$termination == "D8_loop"
  )
)

save_table(
  summary,
  "PNet_summary_QA.csv"
)

cat("\n")
cat("========================================\n")
cat("P-NET COMPLETED\n")
cat("========================================\n")

print(summary)

cat(
  "\nOutput directory:\n",
  OUT,
  "\n"
)
