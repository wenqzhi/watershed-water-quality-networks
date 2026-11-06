source("config.R")

# ================================================================
# Chesapeake Bay Water Quality Propagation Network
#
# Fusion-Net + Connectivity-dependent Dual-Lag System
#
# Nutrients: TN / TP
# Period: 2011-10 to 2023-09
#
# Fusion:
# W_fusion = omega * T_norm + (1 - omega) * P_norm
#
# Dual-Lag:
# R = 1 -> Physical travel time
# R = 0 and TE p < 0.1 -> TE lag
# R = unknown -> Undefined
#
# NOTE:
# p < 0.1 is exploratory and unadjusted.
# Physical travel-time scenarios are not calibrated speeds.
# ================================================================


# ================================================================
# 0. Packages
# ================================================================

pkgs <- c(
  "readr", "dplyr", "tidyr", "stringr",
  "tibble", "ggplot2", "igraph",
  "sf", "scales", "patchwork"
)

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

invisible(lapply(pkgs, library, character.only = TRUE))


# ================================================================
# 1. Paths
# ================================================================

ROOT <- get_data_root()

NUTRIENTS <- c("TN", "TP")

TE_ROOT <- file.path(
  ROOT, "TE_network_results"
)

PHYS_ROOT <- file.path(
  ROOT, "PNet_2011_2023_BlockTracing"
)

OUT_ROOT <- file.path(
  ROOT, "Fusion_DualLag_2011_2023"
)

STATION_FILE <- file.path(
  ROOT, "NTN-station.txt"
)

RIVER_FILE <- file.path(
  ROOT, "水系/水系.shp"
)

BASIN_FILE <- file.path(
  ROOT, "ntnwatersheds/ntnwatersheds.shp"
)

dir.create(
  OUT_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)


# ================================================================
# 2. Main parameters
# ================================================================

# Fusion weights to evaluate
OMEGA_VALUES <- seq(0, 1, 0.1)

# Display setting - not cross-validated optimum
OMEGA_MAIN <- 0.5

# Strong-link percentiles
BACKBONE_Q <- c(
  0.75, 0.80, 0.85, 0.90
)

BACKBONE_MAIN <- 0.85

# TE significance
TE_P_THRESHOLD <- 0.10

# Use raw p values for exploratory TE screening
TE_SIGNIFICANCE <- "raw_p"

# Whether the TE screening is applied to the
# TE layer before fusion:
# FALSE = preserve continuous TE weights;
#         use p<0.1 only for TE-lag assignment.
# TRUE  = non-significant TE links have TE weight 0.
#
# FALSE follows the manuscript's continuous
# TE-weight fusion most closely.
FILTER_TE_BEFORE_FUSION <- FALSE

# Physical time
PHYSICAL_TIME_FIELD <- "travel_days_medium"

# If unavailable, calculate a scenario from
# river distance and an assumed velocity.
ALLOW_PHYSICAL_FALLBACK <- TRUE

VELOCITY_M_S <- 0.50

# Decay for physical connectivity weight
DECAY_DAYS <- 7

DAYS_PER_MONTH <- 30.4375

# IMPORTANT:
# Do not automatically treat missing P-Net
# edges as confirmed physically unreachable.
ASSUME_MISSING_PHYS_AS_ZERO <- FALSE

# Export all parameter combinations
EXPORT_ALL_BACKBONES <- TRUE

# Spatial display
MAX_DRAW_EDGES <- 350


# ================================================================
# 3. Helper functions
# ================================================================

norm_id <- function(x) {
  
  x <- trimws(as.character(x))
  x <- sub("\\.0+$", "", x)
  
  stringr::str_pad(
    x, width = 8, pad = "0"
  )
}


read_edge_file <- function(file) {
  
  if (!file.exists(file)) {
    stop("Missing input file: ", file)
  }
  
  readr::read_csv(
    file,
    col_types = readr::cols(
      from = readr::col_character(),
      to = readr::col_character(),
      .default = readr::col_guess()
    ),
    show_col_types = FALSE
  ) %>%
    mutate(
      from = norm_id(from),
      to = norm_id(to)
    )
}


check_unique_edges <- function(x, label) {
  
  dup <- x %>%
    count(from, to) %>%
    filter(n > 1)
  
  if (nrow(dup) > 0) {
    stop(
      label,
      " contains duplicate directed edges."
    )
  }
}


write_result <- function(x, dir, filename) {
  
  readr::write_csv(
    x,
    file.path(dir, filename),
    na = ""
  )
}


# Min-max normalization over positive edges.
# Zero remains zero.
normalize_layer <- function(x) {
  
  x <- as.numeric(x)
  
  x[!is.finite(x)] <- 0
  x[x < 0] <- 0
  
  out <- numeric(length(x))
  
  pos <- which(x > 0)
  
  if (!length(pos)) {
    return(out)
  }
  
  xmin <- min(x[pos])
  xmax <- max(x[pos])
  
  if (xmax > xmin) {
    
    out[pos] <- (
      x[pos] - xmin
    ) / (
      xmax - xmin
    )
    
  } else {
    
    # Equal positive weights:
    # preserve support without division by zero.
    out[pos] <- 1
  }
  
  out
}


# ================================================================
# 4. Read TE-Net
# ================================================================

read_te <- function(nutrient) {
  
  file <- file.path(
    TE_ROOT,
    nutrient,
    "main_2011_2023_all_edges.csv"
  )
  
  raw <- read_edge_file(file)
  
  candidates <- c(
    "TE_bias_corrected",
    "TE_effective",
    "TE_nats",
    "TE"
  )
  
  weight_field <- candidates[
    candidates %in% names(raw)
  ][1]
  
  if (is.na(weight_field)) {
    stop("TE weight column not found.")
  }
  
  if (!all(
    c("p_value", "lag_months") %in% names(raw)
  )) {
    stop(
      "TE file requires p_value and lag_months."
    )
  }
  
  te <- raw %>%
    transmute(
      from,
      to,
      TE_raw = as.numeric(
        .data[[weight_field]]
      ),
      TE_p = as.numeric(p_value),
      tau_TE_months = as.numeric(lag_months)
    ) %>%
    filter(from != to)
  
  check_unique_edges(te, "TE-Net")
  
  te <- te %>%
    mutate(
      TE_q = p.adjust(
        TE_p,
        method = "BH"
      ),
      
      TE_significant = is.finite(TE_p) &
        TE_p < TE_P_THRESHOLD,
      
      TE_weight = pmax(
        coalesce(TE_raw, 0),
        0
      )
    )
  
  if (FILTER_TE_BEFORE_FUSION) {
    
    te <- te %>%
      mutate(
        TE_weight = ifelse(
          TE_significant,
          TE_weight,
          0
        )
      )
  }
  
  te
}


# ================================================================
# 5. Read P-Net
# ================================================================

read_phys <- function(nutrient) {
  
  file <- file.path(
    PHYS_ROOT,
    nutrient,
    "PNet_directed_reachable_pairs.csv"
  )
  
  raw <- read_edge_file(file)
  
  if (!"river_distance_km" %in% names(raw)) {
    stop("P-Net lacks river_distance_km.")
  }
  
  distance <- as.numeric(
    raw$river_distance_km
  )
  
  if (PHYSICAL_TIME_FIELD %in% names(raw)) {
    
    travel_days <- as.numeric(
      raw[[PHYSICAL_TIME_FIELD]]
    )
    
  } else if (ALLOW_PHYSICAL_FALLBACK) {
    
    travel_days <- distance * 1000 /
      (VELOCITY_M_S * 86400)
    
    warning(
      nutrient,
      ": physical travel time uses ",
      "assumed constant velocity."
    )
    
  } else {
    
    travel_days <- rep(
      NA_real_,
      nrow(raw)
    )
  }
  
  phys <- raw %>%
    mutate(
      tau_Phys_days = travel_days
    ) %>%
    transmute(
      from,
      to,
      
      river_distance_km =
        as.numeric(river_distance_km),
      
      tau_Phys_days =
        as.numeric(tau_Phys_days)
    ) %>%
    filter(
      from != to,
      is.finite(river_distance_km),
      river_distance_km > 0
    ) %>%
    arrange(
      from, to, river_distance_km
    ) %>%
    group_by(from, to) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    mutate(
      R_phys = 1L,
      
      # Physical propagation strength
      P_raw = ifelse(
        is.finite(tau_Phys_days) &
          tau_Phys_days >= 0,
        
        exp(
          -tau_Phys_days / DECAY_DAYS
        ),
        
        NA_real_
      )
    )
  
  check_unique_edges(phys, "Phys-Net")
  
  phys
}


# ================================================================
# 6. Read independently verified connectivity (optional)
#
# Optional CSV:
# PNet_verified_reachability.csv
#
# Columns:
# from,to,R
#
# R=1: physically downstream reachable
# R=0: verified no directed river path
# R=NA: unresolved
# ================================================================

read_verified_R <- function(nutrient) {
  
  file <- file.path(
    PHYS_ROOT,
    nutrient,
    "PNet_verified_reachability.csv"
  )
  
  if (!file.exists(file)) {
    
    return(
      tibble(
        from = character(),
        to = character(),
        R_verified = integer()
      )
    )
  }
  
  x <- read_edge_file(file)
  
  if (!"R" %in% names(x)) {
    stop("Verified connectivity CSV needs R.")
  }
  
  x <- x %>%
    transmute(
      from, to,
      R_verified = as.integer(R)
    )
  
  check_unique_edges(
    x,
    "Verified connectivity"
  )
  
  if (any(
    !is.na(x$R_verified) &
    !x$R_verified %in% c(0L, 1L)
  )) {
    stop("R must be 0, 1, or NA.")
  }
  
  x
}


# ================================================================
# 7. Read station locations
# ================================================================

stations <- readr::read_csv(
  STATION_FILE,
  col_types = readr::cols(
    SITE_NO = readr::col_character(),
    .default = readr::col_guess()
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


# ================================================================
# 8. Read watershed and river shapefiles
# ================================================================

basin <- sf::st_read(
  BASIN_FILE,
  quiet = TRUE
)

rivers <- sf::st_read(
  RIVER_FILE,
  quiet = TRUE
)

if (is.na(sf::st_crs(basin))) {
  stop("Watershed CRS undefined.")
}

if (is.na(sf::st_crs(rivers))) {
  
  bb <- sf::st_bbox(rivers)
  
  if (
    all(is.finite(bb)) &&
    abs(bb["xmin"]) <= 180 &&
    abs(bb["xmax"]) <= 180 &&
    abs(bb["ymin"]) <= 90 &&
    abs(bb["ymax"]) <= 90
  ) {
    
    sf::st_crs(rivers) <- 4326
    
    warning(
      "River CRS assigned EPSG:4326 ",
      "provisionally. Verify coordinates."
    )
    
  } else {
    
    stop(
      "River CRS missing and ",
      "coordinate system unknown."
    )
  }
}

basin <- sf::st_transform(
  basin, 4326
)

rivers <- sf::st_transform(
  rivers, 4326
)


# ================================================================
# 9. Build TE and physical layers
# ================================================================

build_layers <- function(
    te,
    phys,
    verified
) {
  
  nodes <- sort(
    unique(
      c(te$from, te$to)
    )
  )
  
  # All candidate directed pairs
  pairs <- tidyr::expand_grid(
    from = nodes,
    to = nodes
  ) %>%
    filter(from != to)
  
  layers <- pairs %>%
    left_join(
      te,
      by = c("from", "to")
    ) %>%
    left_join(
      phys,
      by = c("from", "to")
    ) %>%
    left_join(
      verified,
      by = c("from", "to")
    )
  
  # Prevent contradictory connectivity
  conflict <- layers %>%
    filter(
      R_phys == 1L &
        R_verified == 0L
    )
  
  if (nrow(conflict)) {
    stop(
      "Contradictory physical connectivity!"
    )
  }
  
  layers <- layers %>%
    mutate(
      R = case_when(
        !is.na(R_verified) ~ R_verified,
        
        R_phys == 1L ~ 1L,
        
        ASSUME_MISSING_PHYS_AS_ZERO ~ 0L,
        
        TRUE ~ NA_integer_
      ),
      
      physical_status = case_when(
        R == 1L ~ "Connected",
        R == 0L ~ "Verified_unconnected",
        TRUE ~ "Unknown"
      ),
      
      TE_weight = coalesce(
        TE_weight, 0
      ),
      
      P_raw = coalesce(
        P_raw, 0
      )
    )
  
  # Normalize TE and Phys separately
  layers <- layers %>%
    mutate(
      T_norm = normalize_layer(
        TE_weight
      ),
      
      P_norm = normalize_layer(
        P_raw
      )
    )
  
  # Two important distinctions:
  # R=1 is topology, not normalized P>0.
  # A minimum positive P may normalize to zero.
  
  layers <- layers %>%
    mutate(
      layer_support = case_when(
        TE_weight > 0 & R == 1L ~
          "TE + Phys",
        
        TE_weight > 0 & R == 0L ~
          "TE / unconnected",
        
        TE_weight > 0 & is.na(R) ~
          "TE / physical unknown",
        
        TE_weight == 0 & R == 1L ~
          "Phys only",
        
        TRUE ~
          "Other"
      )
    )
  
  list(
    layers = layers,
    nodes = nodes
  )
}


# ================================================================
# 10. Fusion weights
# ================================================================

apply_fusion <- function(
    layers,
    omega
) {
  
  layers %>%
    mutate(
      omega = omega,
      
      W_fusion =
        omega * T_norm +
        (1 - omega) * P_norm
    )
}


# ================================================================
# 11. Strong-link backbone
# ================================================================

extract_backbone <- function(
    fusion,
    q
) {
  
  active <- fusion %>%
    filter(
      is.finite(W_fusion),
      W_fusion > 0
    )
  
  if (!nrow(active)) {
    return(active)
  }
  
  threshold <- as.numeric(
    quantile(
      active$W_fusion,
      probs = q,
      na.rm = TRUE
    )
  )
  
  active %>%
    filter(
      W_fusion >= threshold
    ) %>%
    mutate(
      threshold_q = q,
      threshold_W = threshold
    )
}


# ================================================================
# 12. Dual-Lag assignment
# ================================================================

assign_dual_lag <- function(
    backbone
) {
  
  backbone %>%
    mutate(
      tau_Phys_months =
        tau_Phys_days / DAYS_PER_MONTH,
      
      lag_type = case_when(
        
        R == 1L &
          is.finite(tau_Phys_days) &
          tau_Phys_days >= 0 ~
          "Physical",
        
        R == 1L ~
          "Physical_time_missing",
        
        R == 0L &
          TE_significant %in% TRUE &
          is.finite(tau_TE_months) &
          tau_TE_months > 0 ~
          "Information",
        
        is.na(R) ~
          "Physical_status_unknown",
        
        TRUE ~
          "Undefined"
      ),
      
      tau_Fusion_months = case_when(
        
        lag_type == "Physical" ~
          tau_Phys_months,
        
        lag_type == "Information" ~
          tau_TE_months,
        
        TRUE ~
          NA_real_
      ),
      
      lag_source = case_when(
        
        lag_type == "Physical" ~
          "Physical travel-time scenario",
        
        lag_type == "Information" ~
          "Significant TE lag",
        
        TRUE ~
          "Not assigned"
      ),
      
      lag_assigned = lag_type %in% c(
        "Physical",
        "Information"
      )
    )
}


# ================================================================
# 13. Node metrics
# ================================================================

calculate_nodes <- function(
    edges,
    nodes
) {
  
  if (!nrow(edges)) {
    
    return(
      tibble(
        STAID = nodes,
        in_degree = 0L,
        out_degree = 0L,
        in_strength = 0,
        out_strength = 0
      )
    )
  }
  
  graph <- igraph::graph_from_data_frame(
    edges %>%
      select(
        from,
        to,
        W_fusion
      ),
    directed = TRUE,
    vertices = data.frame(
      name = nodes
    )
  )
  
  tibble(
    STAID = igraph::V(graph)$name,
    
    in_degree = as.integer(
      igraph::degree(
        graph,
        mode = "in"
      )
    ),
    
    out_degree = as.integer(
      igraph::degree(
        graph,
        mode = "out"
      )
    ),
    
    in_strength = as.numeric(
      igraph::strength(
        graph,
        mode = "in",
        weights = igraph::E(graph)$W_fusion
      )
    ),
    
    out_strength = as.numeric(
      igraph::strength(
        graph,
        mode = "out",
        weights = igraph::E(graph)$W_fusion
      )
    )
  )
}


# ================================================================
# 14. Matrices
# ================================================================

make_matrix <- function(
    edges,
    nodes,
    field,
    default = NA
) {
  
  m <- matrix(
    default,
    nrow = length(nodes),
    ncol = length(nodes),
    dimnames = list(nodes, nodes)
  )
  
  if (nrow(edges)) {
    
    for (i in seq_len(nrow(edges))) {
      
      m[
        edges$from[i],
        edges$to[i]
      ] <- edges[[field]][i]
    }
  }
  
  m
}


save_matrix <- function(
    m,
    outdir,
    filename
) {
  
  d <- data.frame(
    STAID = rownames(m),
    m,
    check.names = FALSE
  )
  
  write_result(
    d,
    outdir,
    filename
  )
}


# ================================================================
# 15. Spatial network map
# ================================================================

plot_network <- function(
    edges,
    node_metrics,
    nutrient,
    outdir
) {
  
  coords <- stations %>%
    filter(
      STAID %in% node_metrics$STAID
    )
  
  pts <- sf::st_as_sf(
    coords,
    coords = c("lon", "lat"),
    crs = 4326,
    remove = FALSE
  ) %>%
    left_join(
      node_metrics,
      by = "STAID"
    )
  
  p <- ggplot() +
    
    geom_sf(
      data = basin,
      fill = "#F7F7F7",
      color = "grey70",
      linewidth = 0.2
    ) +
    
    geom_sf(
      data = rivers,
      color = "#8DB9D8",
      linewidth = 0.25,
      alpha = 0.8
    )
  
  if (nrow(edges)) {
    
    draw_edges <- edges %>%
      arrange(
        desc(W_fusion)
      ) %>%
      slice_head(
        n = MAX_DRAW_EDGES
      ) %>%
      left_join(
        stations %>%
          rename(
            from = STAID,
            lon1 = lon,
            lat1 = lat
          ),
        by = "from"
      ) %>%
      left_join(
        stations %>%
          rename(
            to = STAID,
            lon2 = lon,
            lat2 = lat
          ),
        by = "to"
      ) %>%
      filter(
        is.finite(lon1),
        is.finite(lon2),
        is.finite(lat1),
        is.finite(lat2)
      )
    
    p <- p +
      
      geom_segment(
        data = draw_edges,
        aes(
          x = lon1,
          y = lat1,
          xend = lon2,
          yend = lat2,
          color = lag_type,
          linewidth = W_fusion
        ),
        alpha = 0.45,
        arrow = grid::arrow(
          length = grid::unit(
            0.04,
            "inches"
          )
        )
      ) +
      
      scale_color_manual(
        values = c(
          "Physical" = "#B55248",
          "Information" = "#557EA8",
          "Physical_time_missing" = "#DB9B57",
          "Physical_status_unknown" = "grey60",
          "Undefined" = "grey80"
        ),
        name = "Lag type"
      ) +
      
      scale_linewidth_continuous(
        range = c(0.2, 1.1),
        name = "Fusion strength"
      )
  }
  
  p <- p +
    
    geom_sf(
      data = pts,
      aes(
        size = out_strength
      ),
      shape = 21,
      fill = "white",
      color = "#214C6A",
      stroke = 0.6,
      alpha = 0.95
    ) +
    
    scale_size_continuous(
      range = c(1.4, 4.5),
      name = "Out-strength"
    ) +
    
    coord_sf(
      datum = NA
    ) +
    
    labs(
      title = paste0(
        nutrient,
        " Fusion-Net with Dual-Lag"
      ),
      subtitle = paste0(
        "2011-2023 | omega=",
        OMEGA_MAIN,
        " | q=",
        BACKBONE_MAIN,
        " | TE p<0.1"
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
      outdir,
      "Fusion_DualLag_map.pdf"
    ),
    p,
    width = 11,
    height = 8.5
  )
  
  ggsave(
    file.path(
      outdir,
      "Fusion_DualLag_map.png"
    ),
    p,
    width = 11,
    height = 8.5,
    dpi = 350
  )
  
  invisible(p)
}


# ================================================================
# 16. Full processing function
# ================================================================

run_nutrient <- function(nutrient) {
  
  cat("\n=====================================\n")
  cat("Processing:", nutrient, "\n")
  cat("=====================================\n")
  
  OUT <- file.path(
    OUT_ROOT,
    nutrient
  )
  
  dir.create(
    OUT,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # ------------------------------------------------------------
  # 16.1 Input
  # ------------------------------------------------------------
  
  te <- read_te(nutrient)
  
  phys <- read_phys(nutrient)
  
  verified <- read_verified_R(nutrient)
  
  prepared <- build_layers(
    te,
    phys,
    verified
  )
  
  layers <- prepared$layers
  nodes <- prepared$nodes
  
  cat(
    "Nodes:", length(nodes), "\n",
    "Candidate pairs:", nrow(layers), "\n"
  )
  
  write_result(
    layers,
    OUT,
    "Fusion_input_layers.csv"
  )
  
  # ------------------------------------------------------------
  # 16.2 Omega and quantile sensitivity
  # ------------------------------------------------------------
  
  sens_list <- list()
  idx <- 1L
  
  for (omega in OMEGA_VALUES) {
    
    fusion <- apply_fusion(
      layers,
      omega
    )
    
    for (q in BACKBONE_Q) {
      
      bb <- extract_backbone(
        fusion, q
      )
      
      dual <- assign_dual_lag(bb)
      
      summary_one <- tibble(
        nutrient = nutrient,
        omega = omega,
        quantile = q,
        
        n_nodes = length(nodes),
        
        n_backbone_edges = nrow(dual),
        
        density = nrow(dual) /
          (
            length(nodes) *
              (length(nodes) - 1)
          ),
        
        n_physical = sum(
          dual$lag_type == "Physical"
        ),
        
        n_information = sum(
          dual$lag_type == "Information"
        ),
        
        n_unknown = sum(
          dual$lag_type ==
            "Physical_status_unknown"
        ),
        
        n_physical_time_missing = sum(
          dual$lag_type ==
            "Physical_time_missing"
        ),
        
        n_undefined = sum(
          dual$lag_type == "Undefined"
        ),
        
        mean_W = if (nrow(dual)) {
          mean(dual$W_fusion)
        } else {
          NA_real_
        }
      )
      
      sens_list[[idx]] <- summary_one
      idx <- idx + 1L
      
      if (EXPORT_ALL_BACKBONES) {
        
        filename <- sprintf(
          "Fusion_DualLag_omega_%0.1f_q%02d.csv",
          omega,
          round(q * 100)
        )
        
        write_result(
          dual,
          OUT,
          filename
        )
      }
    }
  }
  
  sensitivity <- bind_rows(
    sens_list
  )
  
  write_result(
    sensitivity,
    OUT,
    "Fusion_DualLag_sensitivity.csv"
  )
  
  # ------------------------------------------------------------
  # 16.3 Main fusion network
  # ------------------------------------------------------------
  
  fusion_main <- apply_fusion(
    layers,
    OMEGA_MAIN
  )
  
  backbone_main <- extract_backbone(
    fusion_main,
    BACKBONE_MAIN
  )
  
  dual_main <- assign_dual_lag(
    backbone_main
  )
  
  # All candidate edges, including zero weight
  write_result(
    fusion_main,
    OUT,
    "Fusion_all_edges.csv"
  )
  
  # Strong-link network with Dual-Lag labels
  write_result(
    dual_main,
    OUT,
    "Fusion_DualLag_backbone.csv"
  )
  
  assigned <- dual_main %>%
    filter(
      lag_assigned
    )
  
  write_result(
    assigned,
    OUT,
    "Fusion_DualLag_assigned.csv"
  )
  
  # ------------------------------------------------------------
  # 16.4 Node metrics
  # ------------------------------------------------------------
  
  metrics <- calculate_nodes(
    dual_main,
    nodes
  ) %>%
    left_join(
      stations,
      by = "STAID"
    )
  
  write_result(
    metrics,
    OUT,
    "Fusion_node_metrics.csv"
  )
  
  # ------------------------------------------------------------
  # 16.5 Matrices
  # ------------------------------------------------------------
  
  # W: all fusion weights (including weak links)
  W_full <- make_matrix(
    fusion_main,
    nodes,
    "W_fusion",
    default = 0
  )
  
  # W backbone: retained edges only
  W_backbone <- make_matrix(
    dual_main,
    nodes,
    "W_fusion",
    default = 0
  )
  
  # Lag: only assigned backbone edges
  TAU <- make_matrix(
    dual_main,
    nodes,
    "tau_Fusion_months",
    default = NA_real_
  )
  
  TYPE <- make_matrix(
    dual_main,
    nodes,
    "lag_type",
    default = "No_edge"
  )
  
  R <- make_matrix(
    fusion_main,
    nodes,
    "R",
    default = NA_integer_
  )
  
  save_matrix(
    W_full,
    OUT,
    "Fusion_weight_matrix_full.csv"
  )
  
  save_matrix(
    W_backbone,
    OUT,
    "Fusion_weight_matrix_backbone.csv"
  )
  
  save_matrix(
    TAU,
    OUT,
    "Fusion_DualLag_month_matrix.csv"
  )
  
  save_matrix(
    TYPE,
    OUT,
    "Fusion_DualLag_type_matrix.csv"
  )
  
  save_matrix(
    R,
    OUT,
    "Fusion_connectivity_matrix.csv"
  )
  
  # ------------------------------------------------------------
  # 16.6 GraphML
  # ------------------------------------------------------------
  
  g <- igraph::graph_from_data_frame(
    dual_main %>%
      select(
        from,
        to,
        W_fusion,
        T_norm,
        P_norm,
        R,
        lag_type,
        tau_Fusion_months
      ),
    directed = TRUE,
    vertices = data.frame(name = nodes)
  )
  
  igraph::write_graph(
    g,
    file.path(
      OUT,
      "Fusion_DualLag_Network.graphml"
    ),
    format = "graphml"
  )
  
  # ------------------------------------------------------------
  # 16.7 Spatial map
  # ------------------------------------------------------------
  
  plot_network(
    dual_main,
    metrics,
    nutrient,
    OUT
  )
  
  # ------------------------------------------------------------
  # 16.8 Diagnostics
  # ------------------------------------------------------------
  
  lag_summary <- dual_main %>%
    count(
      lag_type,
      name = "n_edges"
    ) %>%
    mutate(
      proportion = n_edges /
        sum(n_edges)
    )
  
  write_result(
    lag_summary,
    OUT,
    "Fusion_DualLag_type_summary.csv"
  )
  
  p_lag <- ggplot(
    lag_summary,
    aes(
      x = lag_type,
      y = proportion,
      fill = lag_type
    )
  ) +
    geom_col(
      width = 0.7
    ) +
    geom_text(
      aes(
        label = scales::percent(
          proportion,
          accuracy = 0.1
        )
      ),
      vjust = -0.5,
      size = 3.5
    ) +
    scale_y_continuous(
      labels = scales::percent_format(),
      expand = expansion(
        mult = c(0, 0.12)
      )
    ) +
    labs(
      title = paste0(
        nutrient,
        " Dual-Lag edge composition"
      ),
      x = NULL,
      y = "Edge proportion"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(
        angle = 20,
        hjust = 1
      )
    )
  
  ggsave(
    file.path(
      OUT,
      "Fusion_DualLag_composition.pdf"
    ),
    p_lag,
    width = 8,
    height = 5
  )
  
  # Lag distributions, if there are assigned edges
  if (nrow(assigned)) {
    
    p_tau <- ggplot(
      assigned,
      aes(
        x = tau_Fusion_months,
        fill = lag_type
      )
    ) +
      geom_histogram(
        bins = 30,
        alpha = 0.8
      ) +
      facet_wrap(
        ~lag_type,
        scales = "free"
      ) +
      labs(
        title = paste0(
          nutrient,
          " Fusion Dual-Lag distributions"
        ),
        x = "Assigned lag (months)",
        y = "Number of edges"
      ) +
      theme_bw(
        base_size = 11
      )
    
    ggsave(
      file.path(
        OUT,
        "Fusion_DualLag_distribution.pdf"
      ),
      p_tau,
      width = 9,
      height = 5
    )
  }
  
  # ------------------------------------------------------------
  # 16.9 Summary
  # ------------------------------------------------------------
  
  summary <- tibble(
    nutrient = nutrient,
    
    n_nodes = length(nodes),
    
    n_candidate_pairs = nrow(layers),
    
    n_TE_p_lt_010 = sum(
      layers$TE_significant %in% TRUE
    ),
    
    n_physical_connected = sum(
      layers$R == 1L,
      na.rm = TRUE
    ),
    
    n_physical_unconnected = sum(
      layers$R == 0L,
      na.rm = TRUE
    ),
    
    n_physical_unknown = sum(
      is.na(layers$R)
    ),
    
    n_backbone = nrow(dual_main),
    
    n_dual_lag_assigned = nrow(assigned),
    
    n_backbone_physical = sum(
      dual_main$lag_type == "Physical"
    ),
    
    n_backbone_information = sum(
      dual_main$lag_type == "Information"
    ),
    
    omega = OMEGA_MAIN,
    
    backbone_quantile = BACKBONE_MAIN
  )
  
  write_result(
    summary,
    OUT,
    "Fusion_DualLag_summary.csv"
  )
  
  print(summary)
  
  invisible(summary)
}


# ================================================================
# 17. Execute TN and TP
# ================================================================

summaries <- list()

for (nutrient in NUTRIENTS) {
  
  summaries[[nutrient]] <-
    run_nutrient(nutrient)
}

final_summary <- bind_rows(
  summaries
)

write_result(
  final_summary,
  OUT_ROOT,
  "Fusion_TN_TP_summary.csv"
)

cat("\n====================================\n")
cat("FUSION + DUAL-LAG COMPLETED\n")
cat("====================================\n")

print(final_summary)

cat(
  "\nResults saved in:\n",
  OUT_ROOT,
  "\n"
)
