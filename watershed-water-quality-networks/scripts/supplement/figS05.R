source("config.R")
# ============================================================
# Figure S5 | Sensitivity and limitations of event-triggered priorities
# Training: 2011-10 ... 2020-12; independent test: 2021-01 ... 2023-12
# (a) P80/P90/P95 test anomaly sensitivity (new test events)
# (b) Global PCS vs source-conditioned PCS (frozen Fig6 events)
# (c) Direct F-Net edge vs multi-step PCS (frozen Fig6 events)
# (d) Reproducibly selected poor-recall event + missed targets
# No full-period data, no retraining, no inference of missing outcomes.
# ============================================================
options(stringsAsFactors=FALSE)
pkgs <- c('readr','dplyr','tidyr','ggplot2','patchwork','tibble','scales')
for (p in pkgs) if (!requireNamespace(p,quietly=TRUE)) install.packages(p)
invisible(lapply(pkgs,library,character.only=TRUE))
ROOT <- get_data_root()
TRAIN <- file.path(ROOT,'Fig6_independent_training_2011_2020')
STEP2 <- file.path(ROOT,'Fig6_independent_validation','Step2_event_replay')
OUT <- file.path(ROOT,'Figure_S5_event_triggered_sensitivity_revised')
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
NTN_FILE <- file.path(ROOT,'NTN_2023_MonthlyLoadTable.csv')
K <- 10L; B <- 1000L; SEED <- 20261101L
TRAIN_START <- as.Date('2011-10-01'); TRAIN_END <- as.Date('2020-12-01')
TEST_START <- as.Date('2021-01-01'); TEST_END <- as.Date('2023-12-01')
HORIZON <- 3L                    # 0..3 months, NOT a causal forecast
REFRACTORY_MONTHS <- 4L          # same-source trigger spacing
PCTS <- c(.80,.90,.95)
MIN_TRAIN_MONTHS <- 36L
COLS <- c(TN='#307FA6',TP='#DE8750')

read_file <- function(path,ids=character()) {
  if (!file.exists(path)) stop('Missing input: ',path)
  spec <- readr::cols(.default=readr::col_guess())
  for (z in ids) spec$cols[[z]] <- readr::col_character()
  readr::read_csv(path,col_types=spec,show_col_types=FALSE)
}
require_fields <- function(x,req,where) {
  miss <- setdiff(req,names(x))
  if (length(miss)) stop(where,' missing columns: ',paste(miss,collapse=', '), '\nPresent: ',paste(names(x),collapse=', '))
}
id8 <- function(x) {
  x <- trimws(sub('\\.0+$','',as.character(x)))
  ix <- !is.na(x) & grepl('^[0-9]{1,8}$',x)
  x[ix] <- sprintf('%08.0f',as.numeric(x[ix]))
  x
}
tobool <- function(x) {
  if (is.logical(x)) return(x)
  z <- tolower(trimws(as.character(x)))
  ifelse(z %in% c('true','t','1','yes'),TRUE,ifelse(z %in% c('false','f','0','no'),FALSE,NA))
}
writeout <- function(x,name) readr::write_csv(x,file.path(OUT,name),na='')
month_idx <- function(date) as.integer(format(date,'%Y'))*12L + as.integer(format(date,'%m'))
add_month <- function(date,n) {
  x <- as.POSIXlt(as.Date(date)); x$mday <- 1L; x$mon <- x$mon + as.integer(n)
  as.Date(x)
}
mean_valid <- function(x) if (all(is.na(x))) NA_real_ else mean(x,na.rm=TRUE)

# ------------------------------------------------------------
# 1. Load frozen TRAINING rankings. Do not use full-period FNet.
# ------------------------------------------------------------
read_train <- function(nu) {
  folder <- file.path(TRAIN,nu)
  s <- read_file(file.path(folder,'04_training_source_conditioned_PCS.csv'),c('source','target'))
  g <- read_file(file.path(folder,'05_training_global_PCS.csv'),c('STAID'))
  w <- read_file(file.path(folder,'02_training_FNet_all_pairs.csv'),c('from','to'))
  require_fields(s,c('source','target','source_PCS'),nu)
  require_fields(g,c('STAID','PCS_in'),nu)
  require_fields(w,c('from','to','W_fus'),nu)
  source_scores <- s %>% transmute(nutrient=nu,source=id8(source),target=id8(target),score=as.numeric(source_PCS))
  global_scores <- g %>% transmute(nutrient=nu,target=id8(STAID),score=as.numeric(PCS_in))
  direct_scores <- w %>% transmute(nutrient=nu,source=id8(from),target=id8(to),score=as.numeric(W_fus))
  for (x in list(source_scores,global_scores,direct_scores)) {
    if (anyNA(x$score) || any(x$score < 0)) stop(nu,': missing/negative training score')
  }
  if (anyDuplicated(source_scores[c('source','target')]) || anyDuplicated(global_scores['target']) || anyDuplicated(direct_scores[c('source','target')])) stop(nu,': duplicate training score keys')
  ids <- sort(unique(c(direct_scores$source,direct_scores$target)))
  list(source=source_scores,global=global_scores,direct=direct_scores,ids=ids)
}
tr <- lapply(c('TN','TP'),read_train)
S <- dplyr::bind_rows(lapply(tr,`[[`,'source'))
G <- dplyr::bind_rows(lapply(tr,`[[`,'global'))
D <- dplyr::bind_rows(lapply(tr,`[[`,'direct'))
NETWORK_IDS <- dplyr::bind_rows(lapply(seq_along(tr),function(i) tibble(nutrient=c('TN','TP')[i],STAID=tr[[i]]$ids)))

# ------------------------------------------------------------
# 2. Shared ranking and event metrics; score ties use common pseudo-random priorities.
# ------------------------------------------------------------
# Stable deterministic tie priority, identical for a given event-target across methods.
# This avoids systematically prioritizing low station IDs when scores are zero.
# Multiple tie realizations would be a stronger sensitivity check if zero ties dominate.
tie_key <- function(event_id, target) {
  z <- paste(event_id,target,sep=':')
  vapply(z,function(t) {
    a <- utf8ToInt(t)
    h <- 17
    for (v in a) h <- (h*131 + v) %% 1000000007
    as.numeric(h)/1000000007
  },numeric(1))
}
rank_one <- function(universe,score_data,join_keys,method) {
  base <- universe %>% left_join(score_data,by=join_keys)
  if (nrow(base)!=nrow(universe)) stop(method,': score table created duplicate event-target rows')
  base %>% mutate(score=ifelse(is.finite(score),score,0),method=method,
                  no_positive_score=score<=0,tie_priority=tie_key(event_id,target)) %>%
    group_by(event_id) %>% arrange(desc(score),tie_priority,target,.by_group=TRUE) %>%
    mutate(rank=row_number()) %>% ungroup()
}
rank_methods <- function(U,methods=c('Source PCS','Global PCS','Direct link')) {
  ans <- list()
  for (m in methods) {
    if (m=='Source PCS') ans[[m]] <- rank_one(U,S,c('nutrient','source','target'),m)
    if (m=='Global PCS') ans[[m]] <- rank_one(U,G,c('nutrient','target'),m)
    if (m=='Direct link') ans[[m]] <- rank_one(U,D,c('nutrient','source','target'),m)
  }
  bind_rows(ans)
}
metrics <- function(rank_df,K=10L) {
  rank_df %>% group_by(nutrient,event_id,source,method) %>%
    group_modify(~{
      z <- .x[order(.x$rank),,drop=FALSE]
      N <- nrow(z); P <- sum(z$hit)
      hits <- sum(z$hit[seq_len(K)])
      tibble(N=N,P=P,hits=hits,precision=hits/K,
             recall=if(P>0) hits/P else NA_real_,
             lift=if(P>0) (hits/K)/(P/N) else NA_real_)
    }) %>% ungroup()
}
cluster_boot <- function(event_metric,group_col='source',B=1000L) {
  # resample sources, carry all events of each sampled source together
  set.seed(SEED)
  splits <- split(event_metric,event_metric[[group_col]])
  groups <- names(splits)
  if(length(groups)<2) return(c(low=NA_real_,high=NA_real_))
  vals <- replicate(B,{
    select <- sample(groups,length(groups),replace=TRUE)
    mean_valid(unlist(lapply(select,function(k) splits[[k]]$value)))
  })
  unname(stats::quantile(vals,c(.025,.975),na.rm=TRUE))
}

# ------------------------------------------------------------
# 3. Frozen Fig6 test outcomes for panels b/c/d
# ------------------------------------------------------------
e0 <- read_file(file.path(STEP2,'03_independent_test_events.csv'),c('source'))
o0 <- read_file(file.path(STEP2,'05_common_event_target_evaluation.csv'),c('source','target'))
require_fields(e0,c('event_id','nutrient','source','trigger_date'),'Step2 events')
require_fields(o0,c('event_id','nutrient','source','target','evaluable','hit_03'),'Step2 target outcomes')
events <- e0 %>% transmute(event_id=as.character(event_id),nutrient=toupper(as.character(nutrient)),source=id8(source),trigger_date=as.Date(trigger_date))
if (anyDuplicated(events$event_id)) stop('Step2 event_id duplicated')
if (any(events$trigger_date < TEST_START | events$trigger_date > TEST_END)) stop('Frozen Step2 events include non-test dates')
U0 <- o0 %>% transmute(event_id=as.character(event_id),nutrient=toupper(as.character(nutrient)),source=id8(source),target=id8(target),evaluable=tobool(evaluable),hit=tobool(hit_03)) %>%
  inner_join(events,by=c('event_id','nutrient','source')) %>%
  filter(source!=target,evaluable %in% TRUE,!is.na(hit))
if(anyDuplicated(U0[c('event_id','target')])) stop('Duplicate event-target observations')
U <- U0 %>% group_by(event_id) %>% filter(n()>=K) %>% ungroup()
if(!nrow(U)) stop('No frozen events with >= K evaluable targets')
rank_frozen <- rank_methods(U)
metric_frozen <- metrics(rank_frozen,K)
writeout(metric_frozen,'S5_bc_event_metrics.csv')
# Ranking QA: how often are top-10 positions filled by absent/zero scores?
rank_qa <- rank_frozen %>% group_by(nutrient,method) %>%
  summarise(n_events=n_distinct(event_id),n_top_slots=sum(rank<=K),
            n_zero_top=sum(rank<=K & no_positive_score),
            pct_zero_top=100*n_zero_top/n_top_slots,.groups='drop')
writeout(rank_qa,'S5_bc_zero_score_ranking_QA.csv')


# (b/c) paired data; exclude events with no positives from recall analyses
wide <- metric_frozen %>% select(nutrient,event_id,source,method,recall,P,N) %>%
  tidyr::pivot_wider(names_from=method,values_from=recall)
if(!all(c('Source PCS','Global PCS','Direct link') %in% names(wide))) stop('Missing method in paired comparison')
compare_pair <- function(dat,baseline) {
  dat %>% filter(is.finite(.data[['Source PCS']]),is.finite(.data[[baseline]])) %>%
    mutate(delta=.data[['Source PCS']]-.data[[baseline]],comparison=baseline)
}
paired <- bind_rows(compare_pair(wide,'Global PCS'),compare_pair(wide,'Direct link'))
paired_stats <- paired %>% group_by(nutrient,comparison) %>% group_modify(~{
  b <- cluster_boot(tibble(source=.x$source,value=.x$delta),B=B)
  tibble(n_events=nrow(.x),n_sources=n_distinct(.x$source),mean_delta=mean(.x$delta),median_delta=median(.x$delta),CI_low=b[1],CI_high=b[2])
}) %>% ungroup()
writeout(paired,'S5_bc_paired_events.csv')
writeout(paired_stats,'S5_bc_paired_CI.csv')

# ------------------------------------------------------------
# 4. Reconstruct P80/90/95 solely from training concentrations
#    Z = (KConc - TRAIN month-of-year mean)/TRAIN month SD.
#    This matches the Fig6 seasonal Z approach, NOT an STL refit.
# ------------------------------------------------------------
n0 <- read_file(NTN_FILE,c('STAID','PCODE'))
require_fields(n0,c('STAID','PCODE','Year','Month','KConc'),'NTN')
ntn <- n0 %>% transmute(STAID=id8(STAID),nutrient=case_when(as.character(PCODE)=='P00600'~'TN',as.character(PCODE)=='P00665'~'TP',TRUE~NA_character_),
                        year=as.integer(Year),month=as.integer(Month),conc=as.numeric(KConc)) %>%
  filter(!is.na(nutrient),year>=2011,year<=2023,month %in% 1:12) %>%
  mutate(date=as.Date(sprintf('%04d-%02d-01',year,month))) %>%
  inner_join(NETWORK_IDS,by=c('nutrient','STAID'))
if (anyDuplicated(ntn[c('nutrient','STAID','date')])) stop('NTN duplicated station-month observations')
train_stat <- ntn %>% filter(date>=TRAIN_START,date<=TRAIN_END,is.finite(conc)) %>%
  group_by(nutrient,STAID,month) %>% summarise(mu=mean(conc),sd=sd(conc),n=n(),.groups='drop')
ntn <- ntn %>% left_join(train_stat,by=c('nutrient','STAID','month')) %>%
  mutate(Z=if_else(is.finite(sd) & sd>0 & is.finite(conc),(conc-mu)/sd,NA_real_))
train_n <- ntn %>% filter(date>=TRAIN_START,date<=TRAIN_END) %>% group_by(nutrient,STAID) %>%
  summarise(n_train_z=sum(is.finite(Z)),.groups='drop')
ntn <- ntn %>% inner_join(train_n %>% filter(n_train_z>=MIN_TRAIN_MONTHS),by=c('nutrient','STAID'))
train_z <- ntn %>% filter(date>=TRAIN_START,date<=TRAIN_END,is.finite(Z)) %>%
  group_by(nutrient,STAID) %>% summarise(q80=unname(quantile(Z,.80)),q90=unname(quantile(Z,.90)),q95=unname(quantile(Z,.95)),.groups='drop')
ntn <- ntn %>% left_join(train_z,by=c('nutrient','STAID'))
writeout(train_z,'S5_a_train_station_thresholds.csv')

make_events <- function(qcol,tag) {
  z <- ntn %>% filter(date>=TEST_START,date<=TEST_END) %>%
    mutate(flag=is.finite(Z) & is.finite(.data[[qcol]]) & Z>.data[[qcol]]) %>%
    arrange(nutrient,STAID,date)
  # Only source triggers: keep first trigger of a run and enforce 4-month separation.
  src <- z %>% group_by(nutrient,STAID) %>%
    mutate(previous_flag=lag(flag,default=FALSE),new_onset=flag & !previous_flag) %>% ungroup() %>%
    filter(new_onset, date<=add_month(TEST_END,-HORIZON)) %>%
    group_by(nutrient,STAID) %>% arrange(date,.by_group=TRUE) %>%
    group_modify(~{
      if(!nrow(.x)) return(.x)
      keep <- rep(FALSE,nrow(.x)); prev <- -Inf
      for (i in seq_len(nrow(.x))) {
        m <- month_idx(.x$date[i]); if(m-prev>=REFRACTORY_MONTHS) {keep[i]<-TRUE;prev<-m}
      }
      .x[keep,,drop=FALSE]
    }) %>% ungroup() %>%
    mutate(event_id=paste0(tag,'_',sprintf('%05d',row_number()))) %>%
    transmute(event_id,nutrient,source=STAID,trigger_date=date,source_Z=Z)
  if(!nrow(src)) stop('No test events under ',tag)
  # Score candidate targets under exact same full-observed 0..3 rule as frozen Fig6.
  candidates <- src %>% inner_join(NETWORK_IDS,by='nutrient',relationship='many-to-many') %>%
    rename(target=STAID) %>% filter(source!=target)
  # Build observation label for four future/same-months. Missing != negative.
  target_obs <- z %>% select(nutrient,target=STAID,date,flag,Z)
  windows <- candidates %>% tidyr::crossing(lag_month=0:HORIZON) %>%
    mutate(date=add_month(trigger_date,lag_month)) %>%
    left_join(target_obs,by=c('nutrient','target','date'))
  Uq <- windows %>% group_by(event_id,nutrient,source,target,trigger_date) %>%
    summarise(n_valid=sum(is.finite(Z)),hit=any(flag %in% TRUE),.groups='drop') %>%
    filter(n_valid==HORIZON+1L) %>% group_by(event_id) %>% filter(n()>=K) %>% ungroup()
  if(!nrow(Uq)) stop('No eligible complete-observation events for ',tag)
  writeout(src,paste0('S5_a_events_',tag,'.csv'))
  writeout(Uq,paste0('S5_a_evaluable_targets_',tag,'.csv'))
  rankq <- rank_methods(Uq,methods='Source PCS')
  mq <- metrics(rankq,K) %>% mutate(threshold=tag)
  mq
}
qa <- bind_rows(lapply(seq_along(PCTS),function(i) make_events(c('q80','q90','q95')[i],c('P80','P90','P95')[i])))
writeout(qa,'S5_a_event_metrics.csv')
a_summary <- qa %>% group_by(nutrient,threshold) %>% group_modify(~{
  z <- .x
  zz <- z %>% filter(is.finite(recall))
  b <- cluster_boot(tibble(source=zz$source,value=zz$recall),B=B)
  tibble(n_events=nrow(z),n_recall_defined=nrow(zz),n_sources=n_distinct(z$source),
         recall=mean_valid(z$recall),precision=mean(z$precision),lift=mean_valid(z$lift),
         mean_baseline_prevalence=mean(z$P/z$N),CI_low=b[1],CI_high=b[2])
}) %>% ungroup() %>% mutate(threshold=factor(threshold,levels=c('P80','P90','P95')))
writeout(a_summary,'S5_a_threshold_summary.csv')

# ------------------------------------------------------------
# 5. (d) Select one failed event separately for TN and TP.
#    Show exactly Top 10 plus the five highest-ranked missed positives.
# ------------------------------------------------------------
source_event <- metric_frozen %>% filter(method=='Source PCS',is.finite(recall),P>0) %>%
  left_join(events,by=c('event_id','nutrient','source'))
chosen <- source_event %>% group_by(nutrient) %>%
  arrange(recall,desc(P),trigger_date,event_id,.by_group=TRUE) %>%
  slice(1) %>% ungroup()
if (!setequal(chosen$nutrient,c('TN','TP'))) stop('Could not select both TN and TP failed events')
shown <- bind_rows(lapply(seq_len(nrow(chosen)), function(i) {
  ev <- chosen[i,]
  x <- rank_frozen %>% filter(method=='Source PCS',event_id==ev$event_id) %>% arrange(rank)
  top <- x %>% filter(rank<=K)
  missed <- x %>% filter(rank>K,hit) %>% slice_head(n=5)
  bind_rows(top,missed) %>% mutate(
    plot_group=if_else(rank<=K,'Top 10 recommendations','Five highest-ranked misses'),
    status=case_when(rank<=K & hit~'Top10: anomaly',rank<=K~'Top10: normal',TRUE~'Missed anomaly'),
    station_rank=sprintf('Rank %02d | %s',rank,target),
    n_all_missed=ev$P-ev$hits,
    recall_event=ev$recall,
    all_positives=ev$P,
    trigger_date=ev$trigger_date
  )
}))
writeout(chosen,'S5_d_selected_events_TN_TP.csv')
writeout(shown,'S5_d_top10_and_5_misses.csv')

# ------------------------------------------------------------
# 6. Build figure. Panel d is faceted, still one logical panel.
# ------------------------------------------------------------
theme_s5 <- theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),plot.title=element_text(face='bold',size=12),strip.text=element_text(face='bold'),legend.position='bottom')
p_a <- ggplot(a_summary,aes(x=threshold,y=recall,color=nutrient,group=nutrient))+
  geom_errorbar(aes(ymin=CI_low,ymax=CI_high),width=.10,linewidth=.5)+geom_line(linewidth=1)+geom_point(size=2.6)+
  facet_wrap(~nutrient,nrow=1)+scale_color_manual(values=COLS,guide='none')+
  scale_y_continuous(labels=scales::percent_format(accuracy=1))+
  labs(title='(a) Anomaly threshold sensitivity',subtitle='P80/P90/P95 from training Z; source-cluster 95% bootstrap CI',x='Training-period anomaly percentile',y='Mean event Recall@10')+theme_s5

paired_plot <- function(baseline,title,xlab,ylab) {
  z <- paired %>% filter(comparison==baseline)
  stat <- paired_stats %>% filter(comparison==baseline) %>%
    mutate(label=sprintf('Mean delta %+.3f [%.3f, %.3f]\nn=%d',mean_delta,CI_low,CI_high,n_events))
  ggplot(z,aes(x=.data[[baseline]],y=.data[['Source PCS']],color=nutrient))+
    geom_abline(intercept=0,slope=1,linetype='dashed',color='#888888')+
    geom_point(alpha=.50,size=1.8,position=position_jitter(width=.005,height=.005,seed=SEED))+
    geom_text(data=stat,aes(x=.02,y=.48,label=label),inherit.aes=FALSE,
              hjust=0,vjust=1,size=2.6,color='#323232')+
    facet_wrap(~nutrient,nrow=1)+scale_color_manual(values=COLS,guide='none')+
    coord_equal(xlim=c(0,.5),ylim=c(0,.5),expand=FALSE)+
    scale_x_continuous(breaks=seq(0,.5,.1))+scale_y_continuous(breaks=seq(0,.5,.1))+
    labs(title=title,subtitle='Shared events | source-cluster bootstrap 95% CI',x=xlab,y=ylab)+theme_s5
}
p_b <- paired_plot('Global PCS','(b) Global vs source-conditioned PCS','Global PCS Recall@10','Source-conditioned PCS Recall@10')
p_c <- paired_plot('Direct link','(c) Direct vs multi-step F-Net','Direct-link Recall@10','Multi-step Recall@10')

# Panel (d): one facet per nutrient (TN | TP), 15 bars in each facet.
# Top 10 recommended stations (grey/orange), plus five best-ranked
# missed anomalous stations (red). No separate rows for recommendation groups.
plot_d <- shown %>%
  mutate(
    nutrient = factor(nutrient, levels=c('TN','TP')),
    # Unique station key across TN/TP and across rank groups.
    row_key = paste(nutrient, sprintf('%03d', rank), target, sep='::'),
    station_label = sprintf('Rank %02d | %s', rank, target),
    status = factor(
      status,
      levels=c('Top10: normal','Top10: anomaly','Missed anomaly')
    )
  ) %>%
  arrange(nutrient, desc(rank), target)
# Ordering: higher-ranked recommendations appear at top of each facet;
# missed positives (rank >10) appear below them.
plot_d$row_key <- factor(plot_d$row_key, levels=unique(plot_d$row_key))
station_label_map <- stats::setNames(
  plot_d$station_label, as.character(plot_d$row_key)
)

p_d <- ggplot(plot_d, aes(x=score, y=row_key, fill=status)) +
  geom_col(width=.74) +
  facet_wrap(~nutrient, nrow=1, scales='free_y') +
  scale_y_discrete(labels=function(x) unname(station_label_map[x])) +
  scale_fill_manual(
    values=c('Top10: normal'='#AEB5BC',
             'Top10: anomaly'='#D88743',
             'Missed anomaly'='#C45454'),
    drop=FALSE
  ) +
  labs(
    title='(d) Unsuccessful events: TN and TP',
    subtitle=paste(
      sprintf('%s: %s; positives=%d; Recall@10=%.2f',
              chosen$nutrient, chosen$event_id,
              chosen$P, chosen$recall),
      collapse=' | '
    ),
    x='Source-conditioned PCS score',
    y='Target station (global recommendation rank)',
    fill='Observed 0–3-month outcome'
  ) +
  theme_s5 +
  theme(axis.text.y=element_text(size=6.5),
        panel.spacing.x=grid::unit(1.0,'lines'),
        legend.position='bottom')

# Figure (d) QA: 10 top-ranked plus up to 5 missed anomalies in each facet.
d_display_qa <- plot_d %>% group_by(nutrient) %>%
  summarise(n_recommended=sum(rank<=K),
            n_missed_shown=sum(rank>K & hit),
            n_unique_stations=dplyr::n_distinct(target),
            .groups='drop')
if (any(d_display_qa$n_recommended != K))
  stop('S5(d): expected exactly 10 recommendations per nutrient')
writeout(d_display_qa, 'S5_d_display_QA.csv')

fig <- ((p_a|p_b)/(p_c|p_d))+plot_annotation(
  title='Figure S5 | Sensitivity and limitations of event-triggered monitoring priorities',
  subtitle='Training: 2011–2020 | Independent test: 2021–2023 | Top 10 | 0–3-month retrospective outcomes',
  caption='P80/P90/P95 thresholds are trained on seasonal Z and recreate test events; panels b–d use frozen Fig.6 events (Z > 1.5). Recall undefined for events with no positive targets. Source PCS is a network score, not an anomaly probability. Missing observations are excluded, not coded normal.'
)
print(fig)
ggsave(file.path(OUT,'Figure_S5.pdf'),fig,width=18,height=12,bg='white')
ggsave(file.path(OUT,'Figure_S5.png'),fig,width=18,height=12,dpi=400,bg='white')
for (nm in c('a','b','c','d')) ggsave(file.path(OUT,paste0('Figure_S5_',nm,'.pdf')),get(paste0('p_',nm)),width=8,height=5.5,bg='white')
writeout(tibble(item=c('threshold_event_rule','threshold_training','panel_bcd_trigger_rule','time_window','ranking','GLOBAL_COL','bootstrap_reps'),
                value=c('test Z > train station percentile, de-duplicate source onsets','train monthly climatology for Z; station quantile on training Z','frozen Fig6 Z>1.5 Step2 events','0–3 months inclusive','Top10 from frozen training F-Net; common deterministic pseudo-random ties','PCS_in (verify against manuscript)',as.character(B))),
         'S5_analysis_settings_QA.csv')
cat('\nCompleted:',OUT,'\n')
print(a_summary)
print(paired_stats)
print(chosen)
print(rank_qa)
