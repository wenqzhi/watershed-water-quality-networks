source("config.R")
# NTN Chesapeake Bay: Temporal robustness of transfer-entropy networks
# Separate runs of this script do not overwrite the original main network.
# Requires NTN_2023_MonthlyLoadTable.csv
# All tests use original surrogate p-values, including 0.005 and 1.0.

packages <- c("readr", "dplyr", "tidyr", "ggplot2", "igraph", "lubridate", "patchwork")
for (p in packages) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
invisible(lapply(packages, library, character.only = TRUE))

# ------------ USER CONFIG ------------
ROOT <- get_data_root()
CSV  <- file.path(ROOT, "NTN_2023_MonthlyLoadTable.csv")
OUT  <- file.path(ROOT, "TE_temporal_robustness")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

PCODE <- "P00665"      # P00600=TN; P00665 = TP
VARIABLE <- "KConc"
MAX_LAG <- 3L
B <- 199L                # preliminary; higher B for final inferential work
ALPHA <- 0.05
SEED <- 20261023
MAX_NODES <- Inf       # first test 12; set Inf for complete analysis
RUN_GROUP <- "historical_144months"
TOP_FRAC <- 0.05       # exploratory effect-size-only comparison, not significance

# Three NON-OVERLAPPING equal-duration 144-month historical windows.
# Recompute all windows on the SAME set of nodes and directed pairs.
# This comparison evaluates temporal network stability, not predictive accuracy.
WINDOWS <- list(
  hist_1987_1999 = as.Date(c("1987-10-01", "1999-09-01")),
  hist_1999_2011 = as.Date(c("1999-10-01", "2011-09-01")),
  main_2011_2023 = as.Date(c("2011-10-01", "2023-09-01"))
)
active <- names(WINDOWS)
RUN_GROUP <- "historical_144months"

# ------------ INPUT / VALIDATION ------------
d <- readr::read_csv(CSV, col_types = cols(STAID = col_character(), PCODE = col_character(), .default = col_guess()), show_col_types=FALSE) %>%
  mutate(STAID = stringr::str_pad(STAID, 8, pad = "0"),
         Date = as.Date(sprintf("%04d-%02d-01", Year, Month))) %>%
  filter(PCODE == !!PCODE)
stopifnot(VARIABLE %in% names(d), !anyDuplicated(d[c("STAID", "Date")]))
if (any(!is.finite(d[[VARIABLE]]) & !is.na(d[[VARIABLE]]))) stop("Infinite values in selected parameter")

get_valid_ids <- function(window) {
  dates <- seq(window[1], window[2], by = "month")
  d %>% filter(Date %in% dates, is.finite(.data[[VARIABLE]])) %>%
    count(STAID, name = "n") %>% filter(n == length(dates)) %>% pull(STAID)
}

# Determine FIXED node set across selected windows. This is essential for fair
# network comparisons; equal windows can become small for historical periods.
all_valid <- lapply(WINDOWS[active], get_valid_ids)
node_ids <- sort(Reduce(intersect, all_valid))
if (is.finite(MAX_NODES)) node_ids <- head(node_ids, MAX_NODES)
if (length(node_ids) < 3) stop("Fewer than 3 common complete stations across selected windows. Choose fewer windows.")
write_csv(tibble(STAID = node_ids), file.path(OUT,paste0(RUN_GROUP,"_fixed_common_stations.csv")))
cat("Fixed station count across ", RUN_GROUP, ": ", length(node_ids), "\n", sep="")

prepare <- function(window) {
  dates <- seq(window[1], window[2], by = "month")
  wide <- d %>% filter(Date %in% dates, STAID %in% node_ids) %>%
    select(Date, STAID, all_of(VARIABLE)) %>%
    pivot_wider(names_from = STAID, values_from = all_of(VARIABLE)) %>%
    right_join(tibble(Date=dates), by="Date") %>% arrange(Date)
  wide <- wide[,c("Date",node_ids)]
  if (anyNA(wide)) stop("Unexpected missing values in prepared station matrix")
  list(x = as.matrix(wide[,-1,drop=FALSE]), dates = dates)
}

preprocess <- function(x,dates) {
  if(any(x < 0,na.rm=TRUE)) stop("Concentrations cannot be negative")
  z <- log1p(x)
  mo <- lubridate::month(dates)
  tt <- seq_len(nrow(z))
  for(j in seq_len(ncol(z))) {
    v <- z[,j]
    v <- v - ave(v,mo,FUN=mean)
    z[,j] <- residuals(lm(v ~ tt))
  }
  z
}

discretize <- function(z) {
  ans <- apply(z,2,function(x) {
    br <- unname(quantile(x,c(1/3,2/3)))
    if(br[1]>=br[2]) {
      as.integer(cut(rank(x,ties.method="first"),breaks=c(-Inf,length(x)/3,2*length(x)/3,Inf)))
    } else {
      as.integer(cut(x,c(-Inf,br[1],br[2],Inf)))
    }
  })
  colnames(ans) <- colnames(z)
  ans
}

cmi <- function(x,y,z,K=3L) {
  n <- length(x)
  a <- array(tabulate(x + K*(y-1L) + K*K*(z-1L),nbins=K^3),c(K,K,K))/n
  xz <- apply(a,c(1,3),sum); yz <- apply(a,c(2,3),sum); zz <- apply(a,3,sum)
  val <- 0
  for(k in 1:K) for(j in 1:K) for(i in 1:K) {
    if (a[i,j,k] > 0) val <- val + a[i,j,k]*log(a[i,j,k]*zz[k]/(xz[i,k]*yz[j,k]))
  }
  max(0,val)
}

lag_te <- function(x,y) {
  ix <- seq.int(MAX_LAG+1L,length(y))
  vapply(seq_len(MAX_LAG),function(l) cmi(x[ix-l],y[ix],y[ix-1L]),numeric(1))
}

shift_x <- function(x,k) {n <- length(x); x[c((n-k+1):n,seq_len(n-k))]}
shift_candidates <- function(n) {
  k <- 1:(n-1L)
  k[pmin(k,n-k)>max(12L,MAX_LAG+1L)]
}

calc <- function(name) {
  xy <- prepare(WINDOWS[[name]])
  mat <- discretize(preprocess(xy$x,xy$dates))
  n <- ncol(mat); nt <- nrow(mat)
  choices <- shift_candidates(nt)
  if (length(choices)<5) stop("Not enough valid circular shifts")
  pairs <- which(row(matrix(0,n,n)) != col(matrix(0,n,n)),arr.ind=TRUE)
  set.seed(SEED + match(name,names(WINDOWS)))
  results <- vector("list",nrow(pairs))
  cat("Window ",name,": ",n," nodes, ",nt," months, ",nrow(pairs)," directed pairs\n",sep="")
  for(k in seq_len(nrow(pairs))) {
    i <- pairs[k,1]; j <- pairs[k,2]
    x <- mat[,i]; y <- mat[,j]
    obs <- lag_te(x,y); best <- which.max(obs); score <- max(obs)
    null <- replicate(B, max(lag_te(shift_x(x,sample(choices,1)),y)))
    p <- (1 + sum(null >= score))/(B+1)
    results[[k]] <- data.frame(from=colnames(mat)[i],to=colnames(mat)[j],
                               lag_months=best,TE_nats=score,TE_bias_corrected=max(0,score-mean(null)),p_value=p)
    if(k %% 1000L == 0L) cat("  ",k," / ",nrow(pairs),"\n")
  }
  edges <- bind_rows(results) %>% mutate(q_value=p.adjust(p_value,"BH"),significant=q_value<ALPHA)
  write_csv(edges,file.path(OUT,paste0(name,"_all_edges.csv")))
  write_csv(filter(edges,significant),file.path(OUT,paste0(name,"_FDR_edges.csv")))
  g <- graph_from_data_frame(filter(edges,significant),directed=TRUE,vertices=data.frame(name=node_ids))
  node <- tibble(STAID=node_ids,
                 out_strength=strength(g,mode="out",weights=E(g)$TE_bias_corrected),
                 in_strength=strength(g,mode="in",weights=E(g)$TE_bias_corrected),
                 outdegree=degree(g,mode="out"),indegree=degree(g,mode="in"))
  write_csv(node,file.path(OUT,paste0(name,"_node_metrics.csv")))
  write_graph(g,file.path(OUT,paste0(name,"_FDR_network.graphml")),format="graphml")
  summary <- tibble(window=name,start=as.character(min(xy$dates)),end=as.character(max(xy$dates)),
                    months=nt,stations=n,candidate_edges=nrow(edges),raw_p_lt_005=sum(edges$p_value<ALPHA),raw_p_lt_010=sum(edges$p_value<0.10),
                    FDR_edges=sum(edges$significant),FDR_density=mean(edges$significant),
                    median_TE=median(edges$TE_nats),median_bias_corrected_TE=median(edges$TE_bias_corrected),
                    min_p=min(edges$p_value),min_q=min(edges$q_value))
  list(edges=edges,node=node,summary=summary)
}

runs <- setNames(lapply(active,calc),active)
summaries <- bind_rows(lapply(runs,`[[`,"summary"))
write_csv(summaries,file.path(OUT,paste0(RUN_GROUP,"_window_summary.csv")))

# ------------ PAIRWISE NETWORK ROBUSTNESS ------------
# Pairwise comparisons across windows use the EXACT SAME nodes, ordered pairs,
# estimators, and lag search. Full edge TE-rank agreement works even when no
# edges survive BH-FDR. Significant-edge Jaccard can be NA if both are empty.
comparisons <- combn(active,2,simplify=FALSE)
comparison_table <- bind_rows(lapply(comparisons,function(pair) {
  a <- runs[[pair[1]]]$edges
  b <- runs[[pair[2]]]$edges
  stopifnot(identical(a$from,b$from),identical(a$to,b$to))
  top_n <- max(1,ceiling(nrow(a)*TOP_FRAC))
  at <- order(-a$TE_bias_corrected, a$from,a$to)[seq_len(top_n)]
  bt <- order(-b$TE_bias_corrected, b$from,b$to)[seq_len(top_n)]
  ai <- which(a$significant); bi <- which(b$significant)
  union_n <- length(union(ai,bi))
  tibble(window_a=pair[1],window_b=pair[2],
         n_common_stations=length(node_ids),n_ordered_pairs=nrow(a),
         spearman_TE=suppressWarnings(cor(a$TE_bias_corrected,b$TE_bias_corrected,method="spearman")),
         top_fraction=TOP_FRAC,
         top_TE_jaccard=length(intersect(at,bt))/length(union(at,bt)),
         FDR_jaccard=if(union_n==0) NA_real_ else length(intersect(ai,bi))/union_n,
         FDR_common_edges=length(intersect(ai,bi)),
         FDR_edges_a=length(ai),FDR_edges_b=length(bi))
}))
write_csv(comparison_table,file.path(OUT,paste0(RUN_GROUP,"_pairwise_network_comparison.csv")))

# ------------ FOUR PANEL HISTORICAL ROBUSTNESS FIGURE ------------
# A: distribution of p < .1 among tests; B: TE effect-size distribution
# C: correlation of paired directed-edge TE rankings against the main period
# D: overlap + phi coefficient for top-5% TE pathways (effect-size ranking)
MAIN <- "main_2011_2023"
if (!(MAIN %in% names(runs))) stop("Main period missing")
robust_long <- bind_rows(lapply(names(runs), function(w) {
  runs[[w]]$edges %>% mutate(window=w)
}))
plot_a <- ggplot(filter(robust_long,p_value<0.1),aes(x=p_value,fill=window)) +
  geom_histogram(binwidth=1/(B+1),boundary=0,color="white",linewidth=.12) +
  facet_wrap(~window,ncol=1,scales="free_y") +
  scale_fill_brewer(palette="Set2",guide="none") +
  labs(title="A  Surrogate p values < 0.1",x="Surrogate p value",y="Directed pairs") +
  theme_bw(base_size=11) + theme(panel.grid.minor=element_blank())
plot_b <- ggplot(robust_long,aes(x=TE_nats,fill=window)) +
  geom_density(alpha=.3,adjust=1,na.rm=TRUE) +
  labs(title="B  Transfer-entropy strength distributions",x="Observed TE (nats)",y="Density",fill="Period") +
  theme_bw(base_size=11) + theme(legend.position="bottom")
main_edges <- runs[[MAIN]]$edges %>% select(from,to,TE_main=TE_bias_corrected)
other_windows <- setdiff(names(runs),MAIN)
paired <- bind_rows(lapply(other_windows,function(w){
  runs[[w]]$edges %>% select(from,to,TE_history=TE_bias_corrected) %>%
    inner_join(main_edges,by=c("from","to")) %>% mutate(window=w)
}))
# C compares RANKS, not raw effect sizes. Spearman rho annotated in each facet.
paired <- paired %>% group_by(window) %>%
  mutate(rank_history=rank(-TE_history,ties.method="average"),
         rank_main=rank(-TE_main,ties.method="average")) %>% ungroup()
c_rho <- paired %>% group_by(window) %>%
  summarise(rho=cor(rank_history,rank_main,method="spearman"),.groups="drop") %>%
  mutate(label=sprintf("Spearman rho = %.3f",rho))
plot_c <- ggplot(paired,aes(x=rank_history,y=rank_main)) +
  geom_point(alpha=.15,size=.55,color="#376F93") +
  geom_abline(slope=1,intercept=0,linetype="dashed",color="#B45146") +
  geom_text(data=c_rho,aes(x=Inf,y=-Inf,label=label),inherit.aes=FALSE,
            hjust=1.03,vjust=-.5,size=3.5) +
  facet_wrap(~window) +
  labs(title="C  Directed-edge TE ranking stability",
       x="Historical TE rank (1 = strongest)",
       y="Main-period TE rank (1 = strongest)") +
  theme_bw(base_size=11)

# Main vs historical top-k pathway agreement; phi is the Pearson correlation
# between two binary (top-k or not) edge-status vectors; Jaccard is overlap.
path_comparison <- bind_rows(lapply(other_windows,function(w) {
  a <- runs[[MAIN]]$edges; b <- runs[[w]]$edges
  stopifnot(identical(a$from,b$from),identical(a$to,b$to))
  n <- nrow(a); k <- max(1L, ceiling(TOP_FRAC*n))
  aa <- seq_len(n) %in% order(-a$TE_bias_corrected,a$from,a$to)[seq_len(k)]
  bb <- seq_len(n) %in% order(-b$TE_bias_corrected,b$from,b$to)[seq_len(k)]
  tibble(window=w, n_directed_pairs=n, top_k=k,
         common_top_paths=sum(aa&bb),
         jaccard=sum(aa&bb)/sum(aa|bb),
         phi=as.numeric(cor(as.integer(aa),as.integer(bb))),
         spearman=suppressWarnings(cor(a$TE_bias_corrected,b$TE_bias_corrected,method="spearman")))
}))
write_csv(path_comparison,file.path(OUT,paste0(RUN_GROUP,"_main_vs_history_path_overlap.csv")))
# D: Top-5% pathway overlap; random expectation under independent fixed-size sets.
path_comparison <- path_comparison %>%
  mutate(random_expected_jaccard=top_k/(2*n_directed_pairs-top_k))
write_csv(path_comparison,file.path(OUT,paste0(RUN_GROUP,"_main_vs_history_path_overlap.csv")))
plot_d <- ggplot(path_comparison,aes(x=window,y=jaccard)) +
  geom_col(width=.55,fill="#477F9C") +
  geom_point(aes(y=random_expected_jaccard),color="#B45146",size=3,shape=18) +
  geom_text(aes(label=sprintf("%.3f",jaccard)),vjust=-.5,size=3.6) +
  scale_y_continuous(limits=c(0,max(.1,path_comparison$jaccard*1.25)),expand=expansion(mult=c(0,.05))) +
  labs(title="D  Top 5% TE-path overlap (Jaccard)",
       subtitle="Red diamonds: expected overlap for random fixed-size selections",
       x=NULL,y="Jaccard similarity") +
  theme_bw(base_size=11) +
  theme(axis.text.x=element_text(angle=15,hjust=1))

fig <- (plot_a | plot_b) / (plot_c | plot_d) +
  plot_annotation(title="NTN transfer entropy — historical robustness",
                  subtitle=paste0("Fixed stations: ",length(node_ids)," | common directed pairs: ",length(node_ids)*(length(node_ids)-1),
                                  " | max lag: ",MAX_LAG," months | surrogates: ",B,". Windows are non-overlapping; network edges within windows remain dependent."))
ggsave(file.path(OUT,paste0(RUN_GROUP,"_historical_robustness_4panels.pdf")),fig,width=15,height=10,bg="white")
ggsave(file.path(OUT,paste0(RUN_GROUP,"_historical_robustness_4panels.png")),fig,width=15,height=10,dpi=300,bg="white")
write_csv(robust_long,file.path(OUT,paste0(RUN_GROUP,"_all_windows_long_edges.csv")))
print(summaries)
print(path_comparison)
cat("\nFinished: ",OUT,"\n")



# ========================================================================
# NTN historical TE robustness: A p-values, B TE distribution,
# C node out-strength inequality (Gini/Lorenz), D Top-4 source-node retention.
# Uses the PREVIOUSLY COMPUTED three all_edges files; DOES NOT recalculate TE.
# All three networks must share the same directed station pairs.
# ========================================================================
# ========================================================================
# NTN historical TE robustness: A p-values, B TE distribution,
# C node out-strength inequality (Gini/Lorenz), D optimal TE lag distribution.
# Uses the PREVIOUSLY COMPUTED three all_edges files; DOES NOT recalculate TE.
# All three networks must share the same directed station pairs.
# ========================================================================
pkgs <- c("readr", "dplyr", "ggplot2", "patchwork", "tidyr", "scales")
for (x in pkgs) if (!requireNamespace(x, quietly=TRUE)) install.packages(x)
invisible(lapply(pkgs, library, character.only=TRUE))

# Change only this path if your previous script used a different OUT
ROOT <- get_data_root()
OUT <- file.path(ROOT, "TE_temporal_robustness")
FIGDIR <- file.path(OUT, "historical_gini_lag_robustness")
dir.create(FIGDIR, showWarnings=FALSE, recursive=TRUE)

WINDOWS <- c("hist_1987_1999", "hist_1999_2011", "main_2011_2023")
MAIN <- "main_2011_2023"
HIST <- setdiff(WINDOWS,MAIN)
MAX_LAG <- 3L
LAG_FILTER <- "p_lt_0.10"   # "all" or "p_lt_0.10"; unadjusted selection only
P_CUTOFF <- 0.10
B <- 199L
# All pairs (recommended for functional topology); do not use only FDR edges.
# TE_bias_corrected = max(0, observed TE - mean surrogate TE)
STRENGTH_FIELD <- "TE_bias_corrected"
# For comparison, use "TE_nats" above and rerun; report sensitivity.

# -------- load and validate --------
read_edges <- function(w) {
  f <- file.path(OUT,paste0(w,"_all_edges.csv"))
  if (!file.exists(f)) stop("File missing: ",f)
  x <- readr::read_csv(f,show_col_types=FALSE,col_types=readr::cols(
    from=readr::col_character(),to=readr::col_character(),.default=readr::col_guess()))
  required <- c("from","to","p_value","TE_nats",STRENGTH_FIELD,"lag_months")
  if (!all(required %in% names(x))) stop(w, " lacks fields: ", paste(setdiff(required,names(x)),collapse=", "))
  x <- x %>% mutate(from=stringr::str_pad(from,8,pad="0"),
                    to=stringr::str_pad(to,8,pad="0"),
                    window=w)
  if (anyDuplicated(x[c("from","to")])) stop("Duplicate directed pairs: ",w)
  if (any(x$from==x$to)) stop("Self-pairs exist: ",w)
  if (any(!is.finite(x$p_value)|!is.finite(x$TE_nats)|!is.finite(x[[STRENGTH_FIELD]])))
    stop("Nonfinite TE/p values in ",w)
  if (any(x[[STRENGTH_FIELD]] < 0)) stop("Negative bias-corrected TE values in ",w)
  x$lag_months <- suppressWarnings(as.integer(x$lag_months))
  if (anyNA(x$lag_months) || any(!x$lag_months %in% seq_len(MAX_LAG)))
    stop("Unexpected lag_months in ",w,"; expected 1..",MAX_LAG)
  x %>% arrange(from,to)
}
# stringr only used for preserving leading-zero USGS IDs
if (!requireNamespace("stringr",quietly=TRUE)) install.packages("stringr")
runs <- setNames(lapply(WINDOWS,read_edges),WINDOWS)
keys <- lapply(runs,function(x) paste(x$from,x$to,sep="->"))
if (!all(vapply(keys[-1],identical,logical(1),keys[[1]]))) {
  stop("Different directed pairs between periods. Recompute using the SAME station set.")
}
ids <- sort(unique(c(runs[[MAIN]]$from,runs[[MAIN]]$to)))
if (length(ids)<2L) stop("Too few nodes.")
cat("Fixed stations:",length(ids),"; directed edges:",nrow(runs[[MAIN]]),"\n")

# -------- compute node strength using ALL edges --------
node_stats <- dplyr::bind_rows(lapply(WINDOWS,function(w) {
  x <- runs[[w]]
  strength <- x %>% group_by(from) %>%
    summarise(out_strength=sum(.data[[STRENGTH_FIELD]]),
              mean_out_strength=mean(.data[[STRENGTH_FIELD]]),
              n_out_pairs=n(),.groups="drop") %>%
    rename(STAID=from)
  tibble(STAID=ids) %>% left_join(strength,by="STAID") %>%
    mutate(window=w)
}))
if (anyNA(node_stats$out_strength)) stop("Incomplete station list in node strength table.")
readr::write_csv(node_stats,file.path(FIGDIR,"node_out_strength_by_period.csv"))

# -------- C: Gini and Lorenz by period --------
# Gini summarizes across-node concentration of outbound information strength.
# It is descriptive; equality of Gini does not imply stable node identities.
gini_nonnegative <- function(z) {
  if (any(!is.finite(z)) || any(z < 0)) stop("Gini needs finite, nonnegative strength")
  if (sum(z) == 0) return(NA_real_)  # Undefined when all out-strengths are zero
  z <- sort(as.numeric(z))
  n <- length(z)
  sum((2 * seq_len(n) - n - 1) * z) / (n * sum(z))
}
lorenz_df <- bind_rows(lapply(WINDOWS, function(w) {
  x <- node_stats %>% filter(window == w) %>% arrange(out_strength)
  total <- sum(x$out_strength)
  tibble(window=w,
         cum_station=c(0,seq_len(nrow(x))/nrow(x)),
         cum_output=if (total>0) c(0,cumsum(x$out_strength)/total) else rep(NA_real_,nrow(x)+1))
}))
gini_summary <- node_stats %>% group_by(window) %>%
  summarise(n_stations=n(), gini=gini_nonnegative(out_strength),
            total_out_strength=sum(out_strength),
            median_out_strength=median(out_strength),.groups="drop")
readr::write_csv(gini_summary,file.path(FIGDIR,"node_out_strength_gini_by_period.csv"))
readr::write_csv(lorenz_df,file.path(FIGDIR,"node_out_strength_lorenz_curve.csv"))

# -------- D: distribution of optimal tested lags --------
# 'lag_months' is selected as the maximum observed TE over tested lags;
# it is NOT a physically measured travel time. The p filter is exploratory
# and based on existing surrogate p-values without FDR correction.
long_lags <- bind_rows(runs)
lag_frequency <- tidyr::expand_grid(window=WINDOWS, lag_months=seq_len(MAX_LAG),
                                    subset=c("all","p_lt_0.10")) %>%
  left_join(
    bind_rows(
      long_lags %>% mutate(subset="all"),
      long_lags %>% filter(p_value < P_CUTOFF) %>% mutate(subset="p_lt_0.10")
    ) %>%
      count(window,subset,lag_months,name="n_edges"),
    by=c("window","subset","lag_months")
  ) %>%
  mutate(n_edges=replace_na(n_edges,0L)) %>%
  group_by(window,subset) %>%
  mutate(total=sum(n_edges), proportion=if_else(total>0,n_edges/total,NA_real_)) %>%
  ungroup()
write_csv(lag_frequency,file.path(FIGDIR,"optimal_lag_counts_and_proportions.csv"))

# A simple descriptive chi-square across periods ignores dependence between
# dyads and selection over multiple lags; not a valid formal network test.
cat("\nOptimal lag counts and proportions:\n")
print(lag_frequency,n=Inf)

# -------- plot A B: original distributions --------
long <- bind_rows(runs) %>% mutate(window=factor(window,levels=WINDOWS))
labels <- c("1987–1999","1999–2011","2011–2023")
names(labels) <- WINDOWS
palette <- c(hist_1987_1999="#D98A82",hist_1999_2011="#61B99E",main_2011_2023="#8BA7D1")
base <- theme_bw(base_size=11) + theme(panel.grid.minor=element_blank(),
                                       plot.title=element_text(face="bold",size=12))

pA <- ggplot(filter(long,p_value<P_CUTOFF),aes(p_value,fill=window)) +
  geom_histogram(binwidth=1/(B+1),boundary=0,color="white",linewidth=.12) +
  facet_wrap(~window,ncol=1,scales="free_y",labeller=as_labeller(labels)) +
  scale_fill_manual(values=palette,guide="none") +
  labs(title="A  Surrogate p values < 0.1",x="Surrogate p value",y="Directed pairs")+base

pB <- ggplot(long,aes(TE_nats,fill=window,color=window)) +
  geom_density(alpha=.22,linewidth=.8,adjust=1) +
  scale_fill_manual(values=palette,labels=labels) +
  scale_color_manual(values=palette,labels=labels,guide="none") +
  labs(title="B  Transfer-entropy strength distributions",
       x="Observed TE (nats)",y="Density",fill="Period") +
  base + theme(legend.position="bottom")

# -------- plot C: Lorenz curves with Gini coefficients --------
gini_label <- setNames(
  ifelse(is.finite(gini_summary$gini),
         paste0(labels[gini_summary$window], "  (Gini = ",
                sprintf("%.3f",gini_summary$gini),")"),
         paste0(labels[gini_summary$window],"  (Gini undefined)")),
  gini_summary$window
)
pC <- ggplot(lorenz_df,aes(x=cum_station,y=cum_output,color=window)) +
  geom_abline(slope=1,intercept=0,linetype="dashed",color="grey50") +
  geom_line(linewidth=1.05,na.rm=TRUE) +
  scale_color_manual(values=palette,labels=gini_label,name="Period") +
  scale_x_continuous(labels=scales::percent_format(accuracy=1),limits=c(0,1)) +
  scale_y_continuous(labels=scales::percent_format(accuracy=1),limits=c(0,1)) +
  coord_equal() +
  labs(title="C  Concentration of node information output",
       subtitle=paste0("Lorenz curves | Gini from all edges (",STRENGTH_FIELD,")"),
       x="Cumulative share of stations (ascending out-strength)",
       y="Cumulative share of total out-strength") +
  base + theme(legend.position="bottom")

# -------- plot D: optimal lag composition (100% stacked bars) --------
lag_plot <- lag_frequency %>%
  filter(subset==LAG_FILTER) %>%
  mutate(window=factor(window,levels=WINDOWS),
         lag_label=factor(paste0(lag_months," month",ifelse(lag_months==1,"","s")),
                          levels=paste0(seq_len(MAX_LAG)," month",ifelse(seq_len(MAX_LAG)==1,"","s"))))
if (any(lag_plot$total==0)) {
  stop("No edges after selected lag filter in one or more periods. Set LAG_FILTER <- 'all'.")
}
lag_colors <- c("#4F81A3","#74B7A9","#DEA07B")
if (MAX_LAG > length(lag_colors)) lag_colors <- grDevices::hcl.colors(MAX_LAG,"Set 2")
names(lag_colors) <- levels(lag_plot$lag_label)
pD <- ggplot(lag_plot,aes(x=window,y=proportion,fill=lag_label)) +
  geom_col(width=.70,color="white",linewidth=.45) +
  geom_text(aes(label=ifelse(n_edges>0,paste0(n_edges," (",round(100*proportion),"%)"),"")),
            position=position_stack(vjust=.5),size=3.2,color="#243746") +
  scale_fill_manual(values=lag_colors,name="Optimal lag") +
  scale_x_discrete(labels=labels) +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.2),
                     labels=scales::percent_format(accuracy=1),expand=expansion(mult=c(0,0))) +
  labs(title="D  Optimal TE lag composition (1–3 months)",
       subtitle=if (LAG_FILTER=="all") "All directed pairs (descriptive)" else
         paste0("Exploratory: surrogate p < ",P_CUTOFF," (uncorrected); counts shown"),
       x=NULL,y="Proportion of directed pairs") +
  base + theme(legend.position="bottom")

fig <- (pA|pB)/(pC|pD) +
  patchwork::plot_annotation(
    title="NTN transfer entropy — historical structural comparison",
    subtitle=paste0("Fixed stations: ",length(ids),
                    " | common directed pairs: ",nrow(runs[[MAIN]]),
                    " | equal 144-month windows | existing surrogate p-values (B=",B,")"))
ggsave(file.path(FIGDIR,"historical_gini_lag_4panels.pdf"),fig,
       width=15,height=10,bg="white")
ggsave(file.path(FIGDIR,"historical_gini_lag_4panels.png"),fig,
       width=15,height=10,dpi=350,bg="white")

# Summary for manuscripts and transparency
sum_period <- long %>% group_by(window) %>%
  summarise(n_pairs=n(),p_lt_010=sum(p_value<P_CUTOFF),
            p_lt_005=sum(p_value<.05),
            mean_TE=mean(TE_nats),median_TE=median(TE_nats),.groups="drop")
write_csv(sum_period,file.path(FIGDIR,"per_period_significance_strength_summary.csv"))
print(sum_period)
print(gini_summary)
print(lag_frequency,n=Inf)
print(fig)
cat("\nAll outputs in: ",FIGDIR,"\n")
# CAUTIONS:
# - C measures concentration, not identity stability; D reflects the maximum TE lag over only the tested 1–3 months.
# - The edge selection in D is exploratory, with unadjusted p < 0.1.
# - The pooled dyads are not independent; proportions alone do not establish temporal changes.
# - Raw versus bias-corrected TE and lag range should undergo sensitivity analysis.
# - A is UNADJUSTED p-values; do not call p<0.1 FDR significant.
# - Random top-K reference is descriptive only; stations/networks are dependent.

