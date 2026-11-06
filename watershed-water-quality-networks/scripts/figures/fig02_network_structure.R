
# ==============================================================
# FIGURE 2
# Physical–information network structure and consistency
#
# TN / TP
# TE-Net | Phys-Net | Fusion-Net
# Physical support ratios | Null-model comparison
#
# Edge color:
# TE-Net     = T_norm
# Phys-Net   = P_norm
# Fusion-Net = W_fusion
#
# Common color scale: 0–1
# All monitoring stations have equal symbol size
# ==============================================================

# ==============================================================
# 0. Packages
# ==============================================================

pkgs <- c(
  "readr",
  "dplyr",
  "tidyr",
  "stringr",
  "tibble",
  "ggplot2",
  "sf",
  "igraph",
  "patchwork",
  "scales"
)

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

invisible(
  lapply(pkgs, library, character.only = TRUE)
)

# ==============================================================
# 1. Paths
# ==============================================================

source("config.R")
ROOT <- get_data_root()

TE_ROOT <- file.path(
  ROOT, "TE_network_results"
)

PHYS_ROOT <- file.path(
  ROOT, "PNet_2011_2023_BlockTracing"
)

STATION_FILE <- file.path(
  ROOT, "NTN-station.txt"
)

BASIN_FILE <- file.path(
  ROOT, "ntnwatersheds/ntnwatersheds.shp"
)

RIVER_FILE <- file.path(
  ROOT, "水系/水系.shp"
)

OUT <- file.path(
  ROOT, "Fig2_network_consistency"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

# ==============================================================
# 2. Settings
# ==============================================================

NUTRIENTS <- c("TN", "TP")

TE_P_THRESHOLD <- 0.10

OMEGA_FUSION <- 0.5
FUSION_Q <- 0.85

# TE-Net map only: retain strongest significant edges
TE_MAP_Q <- 0.60

# Phys-Net displays ALL recognized physical connections
SHOW_ALL_PHYSICAL_EDGES <- TRUE

PHYS_TIME_FIELD <- "travel_days_medium"

ALLOW_PHYS_FALLBACK <- TRUE
VELOCITY_M_S <- 0.5
DECAY_DAYS <- 7

N_NULL <- 5000L
set.seed(123)

# Node size
NODE_SIZE <- 1.8

# Edge appearance
EDGE_WIDTH <- 0.45
EDGE_ALPHA <- 0.85

# Common normalized-strength palette
EDGE_COLORS <- c(
  
  "#FDE725",   # 黄色：高强度
  "#90D743",  # 黄绿色
  "#35B779",  # 绿色
  "#21918C",  # 青绿色
  "#31688E",  # 蓝色
  "#482878"  # 紫色
  
)

# Physical missing edges are NOT automatically R=0.
# Verified connectivity is used if available.
REQUIRE_VERIFIED_R0 <- TRUE

# ============================================================== 
# 3. Helpers
# ==============================================================

norm_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("\\.0+$", "", x)
  stringr::str_pad(x, width = 8, pad = "0")
}

normalize_positive <- function(x) {
  
  x <- as.numeric(x)
  x[!is.finite(x)] <- 0
  x[x < 0] <- 0
  
  out <- rep(0, length(x))
  idx <- which(x > 0)
  
  if (!length(idx)) return(out)
  
  xmin <- min(x[idx])
  xmax <- max(x[idx])
  
  if (xmax > xmin) {
    out[idx] <- (
      x[idx] - xmin
    ) / (
      xmax - xmin
    )
  } else {
    out[idx] <- 1
  }
  
  out
}

read_edge_file <- function(file) {
  
  if (!file.exists(file)) {
    stop("Missing file: ", file)
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

check_unique <- function(x, label) {
  
  dup <- x %>%
    count(from, to) %>%
    filter(n > 1)
  
  if (nrow(dup)) {
    stop(label, " contains duplicate directed pairs.")
  }
}

write_csv2 <- function(x, filename) {
  readr::write_csv(x, filename, na = "")
}

safe_ratio <- function(a, b) {
  ifelse(b > 0, a / b, NA_real_)
}

# ==============================================================
# 4. Read station coordinates
# ==============================================================

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

# ==============================================================
# 5. Read watershed and river map
# ==============================================================

basin <- sf::st_read(
  BASIN_FILE,
  quiet = TRUE
)

rivers <- sf::st_read(
  RIVER_FILE,
  quiet = TRUE
)

if (is.na(sf::st_crs(basin))) {
  stop(
    "Watershed CRS is missing. ",
    "Please assign its verified source CRS."
  )
}

if (is.na(sf::st_crs(rivers))) {
  
  bb <- sf::st_bbox(rivers)
  
  geographic <- all(is.finite(bb)) &&
    abs(bb["xmin"]) <= 180 &&
    abs(bb["xmax"]) <= 180 &&
    abs(bb["ymin"]) <= 90 &&
    abs(bb["ymax"]) <= 90
  
  if (geographic) {
    
    sf::st_crs(rivers) <- 4326
    
    warning(
      "River CRS temporarily set to EPSG:4326. ",
      "Verify the original CRS."
    )
    
  } else {
    
    stop(
      "River CRS missing and not recognizable ",
      "as geographic coordinates."
    )
  }
}

basin <- sf::st_transform(basin, 4326)
rivers <- sf::st_transform(rivers, 4326)

# ==============================================================
# 6. Read TE-Net
# ==============================================================

read_te <- function(nutrient) {
  
  file <- file.path(
    TE_ROOT,
    nutrient,
    "main_2011_2023_all_edges.csv"
  )
  
  raw <- read_edge_file(file)
  
  weight_candidates <- c(
    "TE_bias_corrected",
    "TE_effective",
    "TE_nats",
    "TE"
  )
  
  weight_col <- weight_candidates[
    weight_candidates %in% names(raw)
  ][1]
  
  if (is.na(weight_col)) {
    stop("Cannot find TE weight column: ", file)
  }
  
  if (!"p_value" %in% names(raw)) {
    stop("Missing p_value in TE table.")
  }
  
  te <- raw %>%
    transmute(
      from,
      to,
      
      TE_weight = pmax(
        as.numeric(.data[[weight_col]]),
        0
      ),
      
      TE_p = as.numeric(p_value),
      
      tau_TE_months = if (
        "lag_months" %in% names(raw)
      ) {
        as.numeric(lag_months)
      } else {
        NA_real_
      }
    ) %>%
    filter(from != to)
  
  check_unique(te, paste0(nutrient, " TE-Net"))
  
  te %>%
    mutate(
      TE_sig = is.finite(TE_p) &
        TE_p < TE_P_THRESHOLD
    )
}

# ==============================================================
# 7. Read Phys-Net
# ==============================================================

read_phys <- function(nutrient) {
  
  file <- file.path(
    PHYS_ROOT,
    nutrient,
    "PNet_directed_reachable_pairs.csv"
  )
  
  raw <- read_edge_file(file)
  
  if (!"river_distance_km" %in% names(raw)) {
    stop("Phys-Net missing river_distance_km.")
  }
  
  distance <- as.numeric(
    raw$river_distance_km
  )
  
  if (PHYS_TIME_FIELD %in% names(raw)) {
    
    travel_days <- as.numeric(
      raw[[PHYS_TIME_FIELD]]
    )
    
  } else if (ALLOW_PHYS_FALLBACK) {
    
    travel_days <- distance * 1000 /
      (VELOCITY_M_S * 86400)
    
    warning(
      nutrient,
      ": using assumed velocity for ",
      "physical travel-time scenario."
    )
    
  } else {
    
    travel_days <- rep(
      NA_real_,
      nrow(raw)
    )
  }
  
  phys <- raw %>%
    transmute(
      from,
      to,
      river_distance_km = distance,
      tau_Phys_days = travel_days
    ) %>%
    filter(
      from != to,
      is.finite(river_distance_km),
      river_distance_km > 0
    ) %>%
    arrange(from, to, river_distance_km) %>%
    group_by(from, to) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    mutate(
      R_phys = 1L,
      P_raw = ifelse(
        is.finite(tau_Phys_days) &
          tau_Phys_days >= 0,
        exp(-tau_Phys_days / DECAY_DAYS),
        NA_real_
      )
    )
  
  check_unique(
    phys,
    paste0(nutrient, " Phys-Net")
  )
  
  phys
}

# ==============================================================
# 8. Optional verified reachability
# ==============================================================

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
  
  raw <- read_edge_file(file)
  
  if (!"R" %in% names(raw)) {
    stop("Verified reachability file requires R.")
  }
  
  out <- raw %>%
    transmute(
      from,
      to,
      R_verified = as.integer(R)
    )
  
  check_unique(out, "Verified reachability")
  
  if (any(
    !is.na(out$R_verified) &
    !(out$R_verified %in% c(0L, 1L))
  )) {
    stop("Verified R must be 0, 1 or NA.")
  }
  
  out
}

# ==============================================================
# 9. Build three network layers
# ==============================================================

build_networks <- function(nutrient) {
  
  message("\nProcessing ", nutrient, "...")
  
  te <- read_te(nutrient)
  phys <- read_phys(nutrient)
  verified <- read_verified_R(nutrient)
  
  # Keep TE-Net nodes as the master node set.
  te_nodes <- sort(
    unique(c(te$from, te$to))
  )
  
  phys_nodes <- sort(
    unique(c(phys$from, phys$to))
  )
  
  # Only stations represented in BOTH input networks
  # are compared in Fig. 2.
  #
  # This prevents missing physical station coverage
  # from being silently interpreted as missing edges.
  nodes <- sort(
    intersect(te_nodes, phys_nodes)
  )
  
  if (length(nodes) < 2) {
    stop(
      nutrient,
      ": fewer than 2 common TE/Phys nodes."
    )
  }
  
  te <- te %>%
    filter(from %in% nodes, to %in% nodes)
  
  phys <- phys %>%
    filter(from %in% nodes, to %in% nodes)
  
  verified <- verified %>%
    filter(from %in% nodes, to %in% nodes)
  
  all_pairs <- tidyr::expand_grid(
    from = nodes,
    to = nodes
  ) %>%
    filter(from != to)
  
  layers <- all_pairs %>%
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
  
  # Check for contradictory topology
  conflict <- layers %>%
    filter(
      R_phys == 1L,
      R_verified == 0L
    )
  
  if (nrow(conflict) > 0) {
    stop(
      nutrient,
      ": physical reachability contradiction."
    )
  }
  
  layers <- layers %>%
    mutate(
      R = case_when(
        !is.na(R_verified) ~ R_verified,
        R_phys == 1L ~ 1L,
        TRUE ~ NA_integer_
      ),
      
      physical_status = case_when(
        R == 1L ~ "Verified downstream",
        R == 0L ~ "Verified no downstream path",
        TRUE ~ "Unresolved"
      ),
      
      TE_weight = coalesce(
        TE_weight, 0
      ),
      
      TE_sig = coalesce(
        TE_sig, FALSE
      ),
      
      P_raw = coalesce(
        P_raw, 0
      )
    ) %>%
    mutate(
      T_norm = normalize_positive(
        TE_weight
      ),
      
      P_norm = normalize_positive(
        P_raw
      ),
      
      W_fusion =
        OMEGA_FUSION * T_norm +
        (1 - OMEGA_FUSION) * P_norm
    )
  
  # TE-Net: significant edges for the map
  te_sig <- layers %>%
    filter(
      TE_sig,
      TE_weight > 0
    )
  
  if (nrow(te_sig) > 0) {
    
    te_threshold <- as.numeric(
      quantile(
        te_sig$TE_weight,
        TE_MAP_Q,
        na.rm = TRUE
      )
    )
    
    te_map <- te_sig %>%
      filter(TE_weight >= te_threshold) %>%
      mutate(
        plot_weight = T_norm
      )
    
  } else {
    
    te_map <- te_sig %>%
      mutate(
        plot_weight = numeric()
      )
  }
  
  # Phys-Net: show ALL recognized connections
  phys_edges <- layers %>%
    filter(
      R == 1L,
      is.finite(river_distance_km),
      river_distance_km > 0
    )
  
  phys_map <- phys_edges %>%
    mutate(
      plot_weight = P_norm
    )
  
  # Fusion-Net: strong-link backbone
  active_fusion <- layers %>%
    filter(
      is.finite(W_fusion),
      W_fusion > 0
    )
  
  if (nrow(active_fusion) > 0) {
    
    fusion_threshold <- as.numeric(
      quantile(
        active_fusion$W_fusion,
        FUSION_Q,
        na.rm = TRUE
      )
    )
    
    fusion_backbone <- active_fusion %>%
      filter(
        W_fusion >= fusion_threshold
      ) %>%
      mutate(
        phys_supported = R == 1L,
        plot_weight = W_fusion
      )
    
  } else {
    
    fusion_backbone <- active_fusion %>%
      mutate(
        phys_supported = logical(),
        plot_weight = numeric()
      )
  }
  
  # ----------------------------------------------------------
  # Physical support statistics
  # ----------------------------------------------------------
  
  n_total_pairs <- length(nodes) *
    (length(nodes) - 1)
  
  n_reachable <- sum(
    layers$R == 1L,
    na.rm = TRUE
  )
  
  te_support <- te_sig %>%
    mutate(
      phys_supported = R == 1L,
      phys_unknown = is.na(R)
    )
  
  n_te_sig <- nrow(te_support)
  
  n_te_supported <- sum(
    te_support$phys_supported,
    na.rm = TRUE
  )
  
  total_te_strength <- sum(
    te_support$TE_weight,
    na.rm = TRUE
  )
  
  phys_te_strength <- sum(
    te_support$TE_weight[
      te_support$phys_supported %in% TRUE
    ],
    na.rm = TRUE
  )
  
  support_stats <- tibble(
    nutrient = nutrient,
    n_nodes = length(nodes),
    
    total_directed_pairs = n_total_pairs,
    
    reachable_pairs = n_reachable,
    
    reachable_pair_ratio = safe_ratio(
      n_reachable,
      n_total_pairs
    ),
    
    unresolved_pairs = sum(
      is.na(layers$R)
    ),
    
    verified_no_path_pairs = sum(
      layers$R == 0L,
      na.rm = TRUE
    ),
    
    sig_te_edges = n_te_sig,
    
    supported_sig_te_edges = n_te_supported,
    
    supported_sig_te_edge_ratio = safe_ratio(
      n_te_supported,
      n_te_sig
    ),
    
    unresolved_sig_te_edges = sum(
      te_support$phys_unknown
    ),
    
    sig_te_weight_sum = total_te_strength,
    
    supported_sig_te_weight_sum =
      phys_te_strength,
    
    supported_sig_te_weight_ratio = safe_ratio(
      phys_te_strength,
      total_te_strength
    ),
    
    fusion_edges = nrow(fusion_backbone),
    
    fusion_phys_supported = sum(
      fusion_backbone$phys_supported,
      na.rm = TRUE
    ),
    
    fusion_phys_supported_ratio = safe_ratio(
      sum(
        fusion_backbone$phys_supported,
        na.rm = TRUE
      ),
      nrow(fusion_backbone)
    )
  )
  
  # ----------------------------------------------------------
  # Null model
  #
  # Shuffle TE weights among the significant TE edges
  # while keeping physical-supported positions fixed.
  #
  # This tests weight enrichment conditional on the
  # observed TE edge set and its fixed physical support.
  #
  # It does NOT test preferential edge overlap against
  # randomized network topology.
  # ----------------------------------------------------------
  
  null_df <- tibble(
    nutrient = character(),
    null_mean_TE = numeric()
  )
  
  null_summary <- tibble(
    nutrient = nutrient,
    observed_mean_TE = NA_real_,
    null_mean = NA_real_,
    null_sd = NA_real_,
    empirical_p_ge = NA_real_,
    supported_edge_n = 0L
  )
  
  n_support <- n_te_supported
  
  if (
    n_te_sig > 0 &&
    n_support > 0 &&
    n_support < n_te_sig
  ) {
    
    observed_mean <- mean(
      te_support$TE_weight[
        te_support$phys_supported %in% TRUE
      ]
    )
    
    null_means <- replicate(
      N_NULL,
      mean(
        sample(
          te_support$TE_weight,
          size = n_support,
          replace = FALSE
        )
      )
    )
    
    empirical_p <- (
      1 + sum(null_means >= observed_mean)
    ) / (N_NULL + 1)
    
    null_df <- tibble(
      nutrient = nutrient,
      null_mean_TE = null_means
    )
    
    null_summary <- tibble(
      nutrient = nutrient,
      observed_mean_TE = observed_mean,
      null_mean = mean(null_means),
      null_sd = sd(null_means),
      empirical_p_ge = empirical_p,
      supported_edge_n = n_support
    )
  }
  
  # ----------------------------------------------------------
  # Network structural metrics
  # ----------------------------------------------------------
  
  network_metrics <- function(
    edges,
    nodes,
    weight_col
  ) {
    
    if (nrow(edges) == 0) {
      
      return(
        tibble(
          n_nodes = length(nodes),
          n_edges = 0L,
          density = 0,
          reciprocity = NA_real_,
          weak_components = length(nodes),
          largest_weak_component_ratio =
            1 / length(nodes),
          mean_weight = NA_real_
        )
      )
    }
    
    ed <- edges %>%
      transmute(
        from,
        to,
        weight = as.numeric(
          .data[[weight_col]]
        )
      )
    
    g <- igraph::graph_from_data_frame(
      ed,
      directed = TRUE,
      vertices = data.frame(
        name = nodes
      )
    )
    
    comp <- igraph::components(
      g,
      mode = "weak"
    )
    
    tibble(
      n_nodes = length(nodes),
      n_edges = nrow(edges),
      
      density = nrow(edges) /
        (length(nodes) * (length(nodes) - 1)),
      
      reciprocity = suppressWarnings(
        igraph::reciprocity(g)
      ),
      
      weak_components = comp$no,
      
      largest_weak_component_ratio =
        max(comp$csize) / length(nodes),
      
      mean_weight = mean(
        ed$weight,
        na.rm = TRUE
      )
    )
  }
  
  metrics <- dplyr::bind_rows(
    
    network_metrics(
      te_sig, nodes, "TE_weight"
    ) %>%
      mutate(layer = "TE-Net"),
    
    network_metrics(
      phys_edges, nodes, "P_raw"
    ) %>%
      mutate(layer = "Phys-Net"),
    
    network_metrics(
      fusion_backbone, nodes, "W_fusion"
    ) %>%
      mutate(layer = "Fusion-Net")
    
  ) %>%
    mutate(nutrient = nutrient) %>%
    relocate(nutrient, layer)
  
  list(
    nodes = nodes,
    layers = layers,
    te_sig = te_sig,
    te_map = te_map,
    phys_edges = phys_edges,
    phys_map = phys_map,
    fusion_backbone = fusion_backbone,
    support_stats = support_stats,
    null_df = null_df,
    null_summary = null_summary,
    metrics = metrics
  )
}

# ==============================================================
# 10. Network map function
# ==============================================================

make_network_plot <- function(
    edges,
    nodes,
    layer_name,
    nutrient,
    title_tag
) {
  
  coords <- stations %>%
    filter(
      STAID %in% nodes
    )
  
  pts <- sf::st_as_sf(
    coords,
    coords = c("lon", "lat"),
    crs = 4326,
    remove = FALSE
  )
  
  draw_edges <- edges %>%
    left_join(
      coords %>%
        rename(
          from = STAID,
          lon1 = lon,
          lat1 = lat
        ),
      by = "from"
    ) %>%
    left_join(
      coords %>%
        rename(
          to = STAID,
          lon2 = lon,
          lat2 = lat
        ),
      by = "to"
    ) %>%
    filter(
      is.finite(lon1),
      is.finite(lat1),
      is.finite(lon2),
      is.finite(lat2),
      is.finite(plot_weight)
    )
  
  p <- ggplot() +
    
    geom_sf(
      data = basin,
      fill = "#F7F7F7",
      color = "#B4B4B4",
      linewidth = 0.18
    ) +
    
    geom_sf(
      data = rivers,
      color = "#A6CBDD",
      linewidth = 0.23,
      alpha = 0.8
    )
  
  if (nrow(draw_edges) > 0) {
    
    p <- p +
      geom_segment(
        data = draw_edges,
        aes(
          x = lon1,
          y = lat1,
          xend = lon2,
          yend = lat2,
          color = plot_weight
        ),
        
        linewidth = EDGE_WIDTH,
        alpha = EDGE_ALPHA,
        
        arrow = grid::arrow(
          length = grid::unit(
            0.030, "inches"
          ),
          type = "closed"
        )
      )
  }
  
  p +
    
    geom_sf(
      data = pts,
      shape = 21,
      size = NODE_SIZE,
      fill = "white",
      color = "#333333",
      stroke = 0.43,
      alpha = 0.98
    ) +
    
    scale_color_gradientn(
      colours = EDGE_COLORS,
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      oob = scales::squish,
      
      name = "Normalized propagation strength",
      
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = grid::unit(5, "cm"),
        barheight = grid::unit(0.35, "cm"),
        ticks.colour = "black",
        frame.colour = "grey50"
      )
    ) +
    
    coord_sf(
      datum = NA,
      expand = FALSE
    ) +
    
    labs(
      title = paste0(
        title_tag, " ",
        nutrient, " ",
        layer_name
      ),
      x = NULL,
      y = NULL
    ) +
    
    theme_bw(
      base_size = 10
    ) +
    
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(
        face = "bold",
        size = 11
      ),
      legend.position = "bottom"
    )
}

# ==============================================================
# 11. Calculate TN and TP
# ==============================================================

res <- lapply(
  NUTRIENTS,
  build_networks
)

names(res) <- NUTRIENTS

support_df <- bind_rows(
  lapply(res, `[[`, "support_stats")
)

null_df <- bind_rows(
  lapply(res, `[[`, "null_df")
)

null_sum <- bind_rows(
  lapply(res, `[[`, "null_summary")
)

metrics_df <- bind_rows(
  lapply(res, `[[`, "metrics")
)

write_csv2(
  support_df,
  file.path(
    OUT, "Fig2_support_statistics.csv"
  )
)

write_csv2(
  null_sum,
  file.path(
    OUT, "Fig2_null_model_summary.csv"
  )
)

write_csv2(
  metrics_df,
  file.path(
    OUT, "Fig2_network_metrics.csv"
  )
)

# ============================================================== 
# 12. Physical-support ratios (Panel G)
# ==============================================================

support_plot_df <- bind_rows(
  
  support_df %>%
    transmute(
      nutrient,
      metric = "Reachable\ndirected pairs",
      value = reachable_pair_ratio
    ),
  
  support_df %>%
    transmute(
      nutrient,
      metric = "Supported\nTE edges",
      value = supported_sig_te_edge_ratio
    ),
  
  support_df %>%
    transmute(
      nutrient,
      metric = "Supported\nTE weight",
      value = supported_sig_te_weight_ratio
    ),
  
  support_df %>%
    transmute(
      nutrient,
      metric = "Supported\nFusion backbone",
      value = fusion_phys_supported_ratio
    )
)

support_plot_df$metric <- factor(
  support_plot_df$metric,
  levels = c(
    "Reachable\ndirected pairs",
    "Supported\nTE edges",
    "Supported\nTE weight",
    "Supported\nFusion backbone"
  )
)

p_support <- ggplot(
  support_plot_df,
  aes(
    x = metric,
    y = value,
    fill = nutrient
  )
) +
  
  geom_col(
    position = position_dodge(
      width = 0.75
    ),
    width = 0.68
  ) +
  
  geom_text(
    aes(
      label = ifelse(
        is.finite(value),
        scales::percent(
          value,
          accuracy = 0.1
        ),
        "NA"
      )
    ),
    position = position_dodge(
      width = 0.75
    ),
    vjust = -0.45,
    size = 3.0
  ) +
  
  scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    ),
    limits = c(0, 0.3)
  ) +
  
  scale_fill_manual(
    values = c(
      "TN" = "#BD625D",
      "TP" = "#5983AD"
    )
  ) +
  
  labs(
    title = "G Physical support ratios",
    x = NULL,
    y = "Proportion",
    fill = NULL
  ) +
  
  theme_bw(
    base_size = 10
  ) +
  
  theme(
    axis.text.x = element_text(
      angle = 20,
      hjust = 1,
      size = 8
    ),
    legend.position = "bottom",
    plot.title = element_text(
      face = "bold",
      size = 11
    )
  )

# ==============================================================
# 13. Null model (Panel H)
# ==============================================================

null_label <- null_sum %>%
  mutate(
    label = paste0(
      "Observed = ",
      round(observed_mean_TE, 3),
      
      "\nNull mean = ",
      round(null_mean, 3),
      
      "\nEmpirical p = ",
      signif(empirical_p_ge, 3)
    )
  )

p_null <- ggplot(
  null_df %>%
    filter(
      is.finite(null_mean_TE)
    ),
  aes(
    x = null_mean_TE
  )
) +
  
  geom_histogram(
    bins = 35,
    fill = "#ADB9C3",
    color = "white",
    linewidth = 0.2
  ) +
  
  geom_vline(
    data = null_sum %>%
      filter(
        is.finite(observed_mean_TE)
      ),
    
    aes(
      xintercept = observed_mean_TE
    ),
    
    color = "#BB5149",
    linewidth = 0.9,
    linetype = "dashed"
  ) +
  
  facet_wrap(
    ~nutrient,
    scales = "free"
  ) +
  
  geom_text(
    data = null_label,
    
    aes(
      x = Inf,
      y = Inf,
      label = label
    ),
    
    hjust = 1.05,
    vjust = 1.5,
    size = 3,
    inherit.aes = FALSE
  ) +
  
  labs(
    title = paste0(
      "H Null-model comparison of ",
      "TE strength"
    ),
    x = "Randomized mean TE strength",
    y = "Frequency"
  ) +
  
  theme_bw(
    base_size = 10
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 11
    )
  )

# ==============================================================
# 14. Six spatial network panels
# ==============================================================

p_A <- make_network_plot(
  edges = res$TN$te_map,
  nodes = res$TN$nodes,
  layer_name = "TE-Net",
  nutrient = "TN",
  title_tag = "A"
)

p_B <- make_network_plot(
  edges = res$TN$phys_map,
  nodes = res$TN$nodes,
  layer_name = "Phys-Net",
  nutrient = "TN",
  title_tag = "B"
)

p_C <- make_network_plot(
  edges = res$TN$fusion_backbone,
  nodes = res$TN$nodes,
  layer_name = "Fusion-Net",
  nutrient = "TN",
  title_tag = "C"
)

p_D <- make_network_plot(
  edges = res$TP$te_map,
  nodes = res$TP$nodes,
  layer_name = "TE-Net",
  nutrient = "TP",
  title_tag = "D"
)

p_E <- make_network_plot(
  edges = res$TP$phys_map,
  nodes = res$TP$nodes,
  layer_name = "Phys-Net",
  nutrient = "TP",
  title_tag = "E"
)

p_F <- make_network_plot(
  edges = res$TP$fusion_backbone,
  nodes = res$TP$nodes,
  layer_name = "Fusion-Net",
  nutrient = "TP",
  title_tag = "F"
)

# ==============================================================
# 15. Combine all panels — fixed layout
# ==============================================================

library(patchwork)

# 六张网络图：2行 × 3列
network_panels <- patchwork::wrap_plots(
  p_A, p_B, p_C,
  p_D, p_E, p_F,
  ncol = 3,
  nrow = 2,
  guides = "collect"
) &
  theme(
    legend.position = "bottom"
  )

# 底部两张统计图：1行 × 2列
bottom_panels <- patchwork::wrap_plots(
  p_support,
  p_null,
  ncol = 2,
  nrow = 1
)

# 上下组合
fig2 <- patchwork::wrap_plots(
  network_panels,
  bottom_panels,
  ncol = 1,
  heights = c(2.4, 1.0)
) +
  patchwork::plot_annotation(
    title = paste0(
      "Physical–information network ",
      "structure and consistency"
    ),
    
    subtitle = paste0(
      "TN and TP | TE p < ",
      TE_P_THRESHOLD,
      " | Fusion omega = ",
      OMEGA_FUSION,
      " | Backbone q = ",
      FUSION_Q
    ),
    
    caption = paste0(
      "Network lines are schematic station-to-station ",
      "connections, not routed river paths. ",
      "Edge color represents normalized strength (0–1)."
    ),
    
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 15
      ),
      plot.subtitle = element_text(size = 10),
      plot.caption = element_text(
        size = 8,
        color = "grey45"
      )
    )
  )

print(fig2)
# ==============================================================
# 16. Export full figure
# ==============================================================

ggsave(
  filename = file.path(
    OUT,
    "Fig2_multilayer_network_normalized_strength.pdf"
  ),
  plot = fig2,
  width = 16,
  height = 15,
  bg = "white"
)

ggsave(
  filename = file.path(
    OUT,
    "Fig2_multilayer_network_normalized_strength.png"
  ),
  plot = fig2,
  width = 16,
  height = 15,
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = file.path(
    OUT,
    "Fig2_multilayer_network_normalized_strength.tif"
  ),
  plot = fig2,
  width = 16,
  height = 15,
  dpi = 400,
  compression = "lzw",
  bg = "white"
)

# ==============================================================
# 17. Export maps separately (optional)
# ==============================================================

map_list <- list(
  TN_TE = p_A,
  TN_Phys = p_B,
  TN_Fusion = p_C,
  TP_TE = p_D,
  TP_Phys = p_E,
  TP_Fusion = p_F
)

for (nm in names(map_list)) {
  
  ggsave(
    filename = file.path(
      OUT,
      paste0(
        "Fig2_",
        nm,
        "_normalized_strength.pdf"
      )
    ),
    plot = map_list[[nm]],
    width = 5,
    height = 7,
    bg = "white"
  )
}

# ==============================================================
# 18. Output underlying network edge tables
# ==============================================================

for (nutrient in NUTRIENTS) {
  
  write_csv2(
    res[[nutrient]]$layers,
    file.path(
      OUT,
      paste0(
        nutrient,
        "_all_aligned_network_edges.csv"
      )
    )
  )
  
  write_csv2(
    res[[nutrient]]$te_map,
    file.path(
      OUT,
      paste0(
        nutrient,
        "_TE_map_edges.csv"
      )
    )
  )
  
  write_csv2(
    res[[nutrient]]$phys_map,
    file.path(
      OUT,
      paste0(
        nutrient,
        "_Phys_map_edges.csv"
      )
    )
  )
  
  write_csv2(
    res[[nutrient]]$fusion_backbone,
    file.path(
      OUT,
      paste0(
        nutrient,
        "_Fusion_backbone_edges.csv"
      )
    )
  )
}

# ==============================================================
# 19. Print results
# ==============================================================

cat("\n")
cat("=============================================\n")
cat("FIGURE 2 COMPLETED\n")
cat("=============================================\n")

cat("\nPhysical support statistics:\n")
print(support_df)

cat("\nNull-model results:\n")
print(null_sum)

cat("\nNetwork structural metrics:\n")
print(metrics_df)

cat("\nInput/output comparison:\n")

for (nutrient in NUTRIENTS) {
  
  cat(
    "\n", nutrient,
    "\nStations:",
    length(res[[nutrient]]$nodes),
    "\nTE edges plotted:",
    nrow(res[[nutrient]]$te_map),
    "\nPhys edges plotted:",
    nrow(res[[nutrient]]$phys_map),
    "\nFusion edges plotted:",
    nrow(res[[nutrient]]$fusion_backbone),
    "\n"
  )
}

cat("\nResults saved to:\n")
cat(OUT, "\n")

print(fig2)

