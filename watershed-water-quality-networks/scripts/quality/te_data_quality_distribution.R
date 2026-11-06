source("config.R")
# NTN Transfer Entropy Network — significance diagnostics
# Run AFTER NTN_TE_network_full.R has produced main_2011_2023_all_edges.csv
# No re-estimation: only summarizes and plots existing TE/p/q.
pkgs <- c("readr","dplyr","tidyr","ggplot2","scales","patchwork")
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
invisible(lapply(pkgs, library, character.only=TRUE))

ROOT <- get_data_root()
OUT <- file.path(ROOT,"TE_network_results")
INPUT <- file.path(OUT,"main_2011_2023_all_edges.csv")
FIG_OUT <- file.path(OUT,"TE_significance_diagnostics")
dir.create(FIG_OUT,recursive=TRUE,showWarnings=FALSE)
B <- 199L
ALPHA <- 0.05

stopifnot(file.exists(INPUT))
d <- readr::read_csv(INPUT,show_col_types=FALSE)
stopifnot(all(c("from","to","TE_nats","p_value") %in% names(d)))
d <- d %>%
  mutate(
    p_value=as.numeric(p_value),
    TE_nats=as.numeric(TE_nats),
    q_value=p.adjust(p_value,method="BH"),
    raw_sig=is.finite(p_value) & p_value < ALPHA,
    fdr_sig=is.finite(q_value) & q_value < ALPHA,
    lag_months=if ("lag_months" %in% names(d)) as.integer(lag_months) else NA_integer_
  )
stopifnot(all(is.finite(d$p_value)), all(d$p_value>=0 & d$p_value<=1))
if (anyDuplicated(d[c("from","to")])) warning("Repeated directed pairs in edge table!")
if (any(d$from==d$to)) warning("Self-connections found!")
if ("q_value" %in% names(readr::read_csv(INPUT,n_max=1,show_col_types=FALSE))) {
  message("Recomputed q_value from p_value over all rows to ensure global BH correction.")
}
n_nodes <- length(unique(c(d$from,d$to)))
m <- nrow(d)
min_p <- min(d$p_value)
p_resolution <- 1/(B+1)
summary <- tibble::tibble(
  nodes_in_edges=n_nodes,
  tested_directed_pairs=m,
  theoretical_pairs=n_nodes*(n_nodes-1),
  bootstrap_B=B,
  p_resolution=p_resolution,
  nominal_p_lt_005=sum(d$raw_sig),
  nominal_percent=100*mean(d$raw_sig),
  BH_FDR_q_lt_005=sum(d$fdr_sig),
  BH_FDR_percent=100*mean(d$fdr_sig),
  observed_min_p=min_p,
  observed_min_q=min(d$q_value),
  min_p_count=sum(d$p_value==min_p),
  median_TE_nats=median(d$TE_nats,na.rm=TRUE),
  max_TE_nats=max(d$TE_nats,na.rm=TRUE)
)
readr::write_csv(summary,file.path(FIG_OUT,"TE_significance_summary.csv"))
readr::write_csv(d,file.path(FIG_OUT,"TE_all_edges_recomputed_BH.csv"))
cat("\n========== TE significance diagnostics ==========\n")
print(summary,width=Inf)
cat("\nTop 15 edges by q and TE:\n")
print(d %>% arrange(q_value,desc(TE_nats)) %>%
        select(from,to,lag_months,TE_nats,p_value,q_value) %>% head(15),n=15)

theme_paper <- theme_bw(base_size=12) +
  theme(panel.grid.minor=element_blank(),
        plot.title=element_text(face="bold",size=12),
        plot.subtitle=element_text(size=9),
        legend.position="bottom")

# A: observed p distribution; avoid falsely interpreting 0–0.05 peak as independent edges.
p1 <- ggplot(d,aes(x=p_value)) +
  geom_histogram(binwidth=1/(B+1),fill="#456A89",color="white",linewidth=.16,
                 boundary=0) +
  geom_vline(xintercept=ALPHA,linetype="dashed",color="#B74C43") +
  labs(title="A  Surrogate-test p values",
       subtitle=paste0("Discrete resolution 1/(B+1) = ",signif(p_resolution,3)),
       x="Permutation p value",y="Directed station pairs") +
  theme_paper

# B: raw p vs BH-adjusted p (log axis for readability)
p2 <- ggplot(d,aes(x=p_value,y=q_value)) +
  geom_point(alpha=.22,size=.75,color="#456A89") +
  geom_hline(yintercept=ALPHA,linetype="dashed",color="#B74C43") +
  geom_abline(slope=1,intercept=0,linetype="dotted",color="grey45") +
  scale_x_continuous(limits=c(0,1),breaks=seq(0,1,.2)) +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.2)) +
  labs(title="B  Raw p versus BH-adjusted q",
       x="Raw p value",y="BH-adjusted q value") + theme_paper

# C: sorted p-values and Benjamini-Hochberg rejection boundary
ranked <- d %>% arrange(p_value) %>%
  mutate(rank=row_number(),BH_boundary=ALPHA*rank/m)
p3 <- ggplot(ranked,aes(x=rank)) +
  geom_line(aes(y=p_value,color="Observed p"),linewidth=.8) +
  geom_line(aes(y=BH_boundary,color="BH threshold"),linewidth=.85) +
  scale_y_continuous(limits=c(0,min(.1,max(min_p*4,.07))),
                     oob=scales::oob_squish) +
  scale_color_manual(values=c("Observed p"="#456A89","BH threshold"="#B74C43")) +
  labs(title="C  Multiple-testing burden",
       subtitle=paste0("m = ",format(m,big.mark=",")," tested directed pairs"),
       x="Sorted test rank",y="p value (zoomed)",color=NULL) +
  theme_paper

# D: TE distribution split by nominal significance
p4 <- ggplot(d,aes(x=TE_nats,fill=raw_sig)) +
  geom_histogram(bins=55,position="identity",alpha=.58,color=NA) +
  scale_fill_manual(values=c("FALSE"="#9AA6AF","TRUE"="#BD6655"),
                    labels=c("p >= 0.05","p < 0.05"),name="Unadjusted test") +
  labs(title="D  Transfer entropy effect sizes",
       x="Observed TE (nats)",y="Directed station pairs") +
  theme_paper

combined <- (p1 | p2) / (p3 | p4) +
  plot_annotation(title="NTN transfer entropy — statistical significance diagnostics",
                  subtitle="Main network: 2011–2023 | inference based on existing all_edges file")
ggsave(file.path(FIG_OUT,"TE_significance_4panels.pdf"),combined,
       width=13,height=9,bg="white")
ggsave(file.path(FIG_OUT,"TE_significance_4panels.png"),combined,
       width=13,height=9,dpi=350,bg="white")

# Diagnostic: which p values occur, and how many
p_counts <- d %>%
  count(p_value,name="n_pairs") %>%
  arrange(p_value)
readr::write_csv(p_counts,file.path(FIG_OUT,"TE_p_value_frequency.csv"))

# Node counts including isolates need station ID csv; inspect if available
station_ids_path <- file.path(OUT,"main_2011_2023_station_ids.csv")
if(file.exists(station_ids_path)) {
  ids <- readr::read_csv(station_ids_path,show_col_types=FALSE)
  cat("\nStations in station_id file:",nrow(ids),"\n")
}
cat("\nFiles saved in:",FIG_OUT,"\n")
