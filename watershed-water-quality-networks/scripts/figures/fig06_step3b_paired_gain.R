# ============================================================================
# FIG6 STEP 3b | Paired incremental value, evidence coverage and unique hits
# Inputs: Step3_fair_comparison/02_all_method_ranked_candidates.csv
#         Step3_fair_comparison/05_all_event_metrics_including_random.csv
# Output: Step3b_paired_gain_QA
# No network retraining; no tuning on test outcomes.
# ============================================================================
options(stringsAsFactors = FALSE)
for (p in c('readr','dplyr','tidyr','tibble','ggplot2')) {
  if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
}
suppressPackageStartupMessages({
  library(readr);library(dplyr);library(tidyr);library(tibble);library(ggplot2)
})
source("config.R")
ROOT <- get_data_root()
STEP3 <- file.path(ROOT,'Fig6_independent_validation','Step3_fair_comparison')
OUT <- file.path(ROOT,'Fig6_independent_validation','Step3b_paired_gain_QA')
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
K <- 5L
B_BOOT <- 2000L              # 500 for QA; 2000+ for final estimates
SEED <- 20261029L
BASE_METHOD <- 'Source-conditioned F-Net'
COMPARATORS <- c('Random (expected)','TE-Net','Phys-Net','Global PCS')
METHOD_ORDER <- c(BASE_METHOD,COMPARATORS)

read_required <- function(name, ids=character()) {
  f <- file.path(STEP3,name)
  if (!file.exists(f)) stop('Step3 file not found: ',f)
  spec <- readr::cols(.default=readr::col_guess())
  for (nm in ids) spec$cols[[nm]] <- readr::col_character()
  readr::read_csv(f,col_types=spec,show_col_types=FALSE)
}
require_cols <- function(x,cols,label) {
  bad <- setdiff(cols,names(x))
  if(length(bad)) stop(label,' missing columns: ',paste(bad,collapse=', '))
}
as_bool <- function(x) {
  if (is.logical(x)) return(x)
  z <- tolower(trimws(as.character(x)))
  ifelse(is.na(z),NA,ifelse(z %in% c('true','t','1','yes'),TRUE,
                            ifelse(z %in% c('false','f','0','no'),FALSE,NA)))
}
id8 <- function(x) {
  x <- trimws(sub('\\.0+$','',as.character(x)))
  idx <- !is.na(x) & grepl('^[0-9]{1,8}$',x)
  x[idx] <- sprintf('%08.0f',as.numeric(x[idx]))
  x
}
mean_safe <- function(x) {x <- x[is.finite(x)]; if(length(x)) mean(x) else NA_real_}
write_out <- function(x,nm) readr::write_csv(x,file.path(OUT,nm),na='')

# ---- 1. Import and verify event-level and ranking identities -----------------
ranked <- read_required('02_all_method_ranked_candidates.csv',c('event_id','source','target'))
metric <- read_required('05_all_event_metrics_including_random.csv',c('event_id'))
require_cols(ranked,c('event_id','nutrient','source','target','hit_03','N','P','method',
                      'rank','score','score_available','strong_score'),'Rankings')
require_cols(metric,c('event_id','nutrient','method','K','N','P','hits','precision','recall','lift'),'Metrics')
ranked <- ranked %>% transmute(
  event_id=as.character(event_id), nutrient=toupper(as.character(nutrient)),
  source=id8(source),target=id8(target), hit_03=as_bool(hit_03),
  N=as.integer(N),P=as.integer(P),method=as.character(method),
  rank=as.integer(rank),score=as.numeric(score),
  score_available=as_bool(score_available),strong_score=as_bool(strong_score))
metric <- metric %>% transmute(
  event_id=as.character(event_id),nutrient=toupper(as.character(nutrient)),
  method=as.character(method),K=as.integer(K),N=as.integer(N),P=as.integer(P),
  hits=as.numeric(hits),precision=as.numeric(precision),
  recall=as.numeric(recall),lift=as.numeric(lift))
if (anyNA(ranked$hit_03)) stop('Ranked candidates contain missing outcome labels.')
if (anyDuplicated(ranked[c('nutrient','event_id','method','target')])) stop('Duplicated event-method-target rows.')
if (anyDuplicated(ranked[c('nutrient','event_id','method','rank')])) stop('Duplicated ranks.')
if (anyDuplicated(metric[c('nutrient','event_id','method','K')])) stop('Duplicated event-method-K metrics.')
if (!all(c(BASE_METHOD,COMPARATORS) %in% unique(metric$method)))
  stop('Missing method(s): ',paste(setdiff(METHOD_ORDER,unique(metric$method)),collapse=', '))
metrics5 <- metric %>% filter(K==!!K,method %in% METHOD_ORDER)
if (!nrow(metrics5)) stop('No K=5 metrics.')
source_keys <- ranked %>% distinct(nutrient,event_id,source)
if (anyDuplicated(source_keys[c('nutrient','event_id')])) stop('Event ID maps to multiple sources.')
metrics5 <- metrics5 %>% left_join(source_keys,by=c('nutrient','event_id'))
if (anyNA(metrics5$source)) stop('Cannot link event metrics to source station.')
event_methods <- metrics5 %>% count(nutrient,event_id,name='n_methods')
if (any(event_methods$n_methods != length(METHOD_ORDER)))
  stop('Some events do not have all five methods; cannot do matched analysis.')
counts <- metrics5 %>% group_by(nutrient,event_id) %>%
  summarise(n_N=n_distinct(N),n_P=n_distinct(P),.groups='drop')
if(any(counts$n_N!=1 | counts$n_P!=1)) stop('N/P denominator differs between methods.')
write_out(metrics5,'01_focalK5_event_metrics.csv')

# ---- 2. Event-paired delta and source-cluster bootstrap ----------------------
wide <- metrics5 %>% select(nutrient,event_id,source,method,K,N,P,hits,precision,recall,lift) %>%
  tidyr::pivot_wider(names_from=method,values_from=c(hits,precision,recall,lift),
                     names_sep='___')
paired <- dplyr::bind_rows(lapply(COMPARATORS,function(ref) {
  inner_join(
    metrics5 %>% filter(method==BASE_METHOD) %>%
      select(nutrient,event_id,source,N,P,
             fnet_hits=hits,fnet_precision=precision,fnet_recall=recall,fnet_lift=lift),
    metrics5 %>% filter(method==ref) %>%
      select(nutrient,event_id,ref_hits=hits,ref_precision=precision,
             ref_recall=recall,ref_lift=lift),
    by=c('nutrient','event_id')) %>%
    mutate(comparator=ref,delta_hits=fnet_hits-ref_hits,
           delta_precision=fnet_precision-ref_precision,
           delta_recall=fnet_recall-ref_recall,
           delta_lift=fnet_lift-ref_lift)
}))
write_out(paired,'02_event_paired_differences.csv')
set.seed(SEED)
boot_one <- function(d) {
  # Resample source IDs with replacement. Keep ALL events attached to each
  # selected source; repeated selections are kept as repeated source clusters.
  sources <- unique(d$source)
  if(length(sources)<2) return(NULL)
  groups <- split(d,d$source)
  selected <- sample(sources,length(sources),replace=TRUE)
  sampled <- dplyr::bind_rows(lapply(selected,function(s) groups[[s]]))
  c(delta_hits=mean_safe(sampled$delta_hits),
    delta_precision=mean_safe(sampled$delta_precision),
    delta_recall=mean_safe(sampled$delta_recall),
    delta_lift=mean_safe(sampled$delta_lift))
}
boot_results <- list()
ix <- 0L
for (nu in sort(unique(paired$nutrient))) {
  for (ref in COMPARATORS) {
    d <- paired %>% filter(nutrient==nu,comparator==ref)
    nsrc <- n_distinct(d$source)
    if(nsrc < 2) stop('Too few independent source clusters: ',nu,' ',ref)
    reps <- t(replicate(B_BOOT,boot_one(d)))
    ix <- ix+1L
    boot_results[[ix]] <- as_tibble(reps) %>%
      mutate(nutrient=nu,comparator=ref,rep=seq_len(n()))
  }
}
bootstrap <- bind_rows(boot_results)
write_out(bootstrap,'03_cluster_bootstrap_draws.csv')
metric_cols <- c('delta_hits','delta_precision','delta_recall','delta_lift')
point_est <- paired %>% group_by(nutrient,comparator) %>%
  summarise(n_events=n(),n_sources=n_distinct(source),
            across(all_of(metric_cols),mean_safe),.groups='drop') %>%
  pivot_longer(cols=all_of(metric_cols),names_to='metric',values_to='estimate')
ci <- bootstrap %>% pivot_longer(cols=all_of(metric_cols),names_to='metric',values_to='estimate_boot') %>%
  group_by(nutrient,comparator,metric) %>%
  summarise(ci_low=if(all(is.na(estimate_boot))) NA_real_ else unname(quantile(estimate_boot,.025,na.rm=TRUE)),
            ci_high=if(all(is.na(estimate_boot))) NA_real_ else unname(quantile(estimate_boot,.975,na.rm=TRUE)),
            .groups='drop')
paired_summary <- left_join(point_est,ci,by=c('nutrient','comparator','metric'))
write_out(paired_summary,'04_paired_gain_95CI.csv')

# ---- 3. Top-K support and tie diagnostics -----------------------------------
# Score-qualified target: positive score for FNet/TE/Global, present verified
# travel-time score for Phys-Net (negative time is fine).
# Ranking a missing score at the end is a deterministic *fallback*, not network evidence.
rank5 <- ranked %>% filter(rank<=K,method %in% setdiff(METHOD_ORDER,'Random (expected)'))
rank5 <- rank5 %>% mutate(
  evidenced=case_when(method=='Phys-Net' ~ score_available %in% TRUE,
                      TRUE ~ strong_score %in% TRUE),
  fallback=!evidenced,
  score_tie_position=(!evidenced) | is.na(score))
coverage_event <- rank5 %>% group_by(nutrient,event_id,source,method) %>%
  summarise(n_top5=n(),n_evidenced=sum(evidenced),n_fallback=sum(fallback),
            hit_evidenced=sum(hit_03 & evidenced),
            hit_fallback=sum(hit_03 & fallback),
            top5_hits=sum(hit_03),.groups='drop')
if(any(coverage_event$n_top5!=K)) stop('Rankings do not have five rows per event.')
coverage_summary <- coverage_event %>% group_by(nutrient,method) %>%
  summarise(n_events=n(),mean_evidenced=mean(n_evidenced),
            fraction_evidenced=sum(n_evidenced)/sum(n_top5),
            total_fallback=sum(n_fallback),fallback_hit_count=sum(hit_fallback),
            evidence_hit_count=sum(hit_evidenced),.groups='drop')
write_out(coverage_event,'05_top5_evidence_event.csv')
write_out(coverage_summary,'06_top5_evidence_summary.csv')

# Additional QA: ties at Kth and unverified target proportions across ALL ranks.
tie_qa <- ranked %>% filter(method %in% setdiff(METHOD_ORDER,'Random (expected)')) %>%
  group_by(nutrient,event_id,method) %>%
  summarise(kth_score=score[rank==K][1],
            kth_score_available=score_available[rank==K][1],
            n_equal_to_kth=if(is.finite(kth_score)) sum(is.finite(score) & score==kth_score) else NA_integer_,
            n_all_scored=sum(score_available,na.rm=TRUE),.groups='drop')
write_out(tie_qa,'07_rank_boundary_ties_QA.csv')

# ---- 4. FNet-only and baseline-only recommended stations --------------------
# Method-paired set comparisons are computed from exactly the same event/target
# universe and the exact same K=5 budget.
top5 <- ranked %>% filter(rank<=K,method %in% setdiff(METHOD_ORDER,'Random (expected)')) %>%
  select(nutrient,event_id,source,method,target,hit_03,rank,score,score_available,strong_score)
novelty_detail <- bind_rows(lapply(setdiff(COMPARATORS,'Random (expected)'),function(ref) {
  f <- top5 %>% filter(method==BASE_METHOD) %>% select(nutrient,event_id,source,target,hit_03)
  r <- top5 %>% filter(method==ref) %>% select(nutrient,event_id,source,target,hit_03)
  f <- f %>% mutate(in_fnet=TRUE)
  r <- r %>% mutate(in_ref=TRUE)
  full_join(f,r,by=c('nutrient','event_id','source','target','hit_03')) %>%
    mutate(in_fnet=coalesce(in_fnet,FALSE),in_ref=coalesce(in_ref,FALSE),
           comparator=ref,
           category=case_when(in_fnet & in_ref ~ 'Shared',
                              in_fnet & !in_ref ~ 'F-Net only',
                              !in_fnet & in_ref ~ 'Baseline only',TRUE ~ 'Neither'))
}))
write_out(novelty_detail,'08_novelty_target_level.csv')
novelty_event <- novelty_detail %>% group_by(nutrient,event_id,source,comparator) %>%
  summarise(n_shared=sum(category=='Shared'),
            n_fnet_only=sum(category=='F-Net only'),
            n_ref_only=sum(category=='Baseline only'),
            fnet_only_hits=sum(category=='F-Net only' & hit_03),
            ref_only_hits=sum(category=='Baseline only' & hit_03),
            shared_hits=sum(category=='Shared' & hit_03),
            net_extra_hits=fnet_only_hits-ref_only_hits,
            .groups='drop')
if(any(novelty_event$n_shared+novelty_event$n_fnet_only!=K |
       novelty_event$n_shared+novelty_event$n_ref_only!=K))
  stop('Set comparison did not preserve the same K.')
novelty_summary <- novelty_event %>% group_by(nutrient,comparator) %>%
  summarise(n_events=n(),mean_overlap=mean(n_shared),
            total_fnet_only=sum(n_fnet_only),total_ref_only=sum(n_ref_only),
            fnet_only_hits=sum(fnet_only_hits),ref_only_hits=sum(ref_only_hits),
            shared_hits=sum(shared_hits),net_extra_hits=sum(net_extra_hits),
            mean_net_extra_hits=mean(net_extra_hits),.groups='drop')
write_out(novelty_event,'09_novelty_event.csv')
write_out(novelty_summary,'10_novelty_summary.csv')

# ---- 5. Three publication-style QA plots ------------------------------------
colors <- c('Random (expected)'='#92989e','TE-Net'='#8564b4',
            'Phys-Net'='#347db1','Global PCS'='#35906d')
theme_sci <- theme_bw(base_size=11)+theme(
  panel.grid.minor=element_blank(),strip.background=element_rect(fill='#f1f3f4'),
  strip.text=element_text(face='bold'),plot.title=element_text(face='bold'),
  legend.position='bottom')
plot_gain <- paired_summary %>% filter(metric=='delta_precision') %>%
  mutate(comparator=factor(comparator,levels=COMPARATORS)) %>%
  ggplot(aes(x=estimate,y=comparator,color=comparator))+
  geom_vline(xintercept=0,linetype='dashed',color='#7e8790')+
  geom_errorbarh(aes(xmin=ci_low,xmax=ci_high),height=.18,linewidth=.7)+
  geom_point(size=2.8)+facet_wrap(~nutrient,nrow=1)+
  scale_color_manual(values=colors,guide='none')+
  scale_x_continuous(labels=function(x) paste0(round(x*100,1),' pp'))+
  labs(title='(a) Paired Precision@5 gain: F-Net minus baseline',
       subtitle='95% source-cluster bootstrap intervals; zero indicates no incremental gain',
       x='Difference in Precision@5 (percentage points)',y=NULL)+theme_sci
plot_coverage <- coverage_summary %>%
  mutate(method=factor(method,levels=setdiff(METHOD_ORDER,'Random (expected)'))) %>%
  ggplot(aes(method,fraction_evidenced,fill=method))+
  geom_col(width=.68,color='white')+
  geom_text(aes(label=sprintf('%.1f%%',100*fraction_evidenced)),vjust=-.35,size=3)+
  facet_wrap(~nutrient,nrow=1)+
  scale_y_continuous(labels=scales::label_percent(accuracy=1),limits=c(0,1.10))+
  labs(title='(b) Top-5 candidates supported by nonzero/verified network evidence',
       subtitle='Unscored target padding is reported, not treated as a supported connection',
       x=NULL,y='Share of Top-5 placements with network evidence')+
  theme_sci+theme(axis.text.x=element_text(angle=25,hjust=1),legend.position='none')
plot_novelty <- novelty_summary %>%
  mutate(comparator=factor(comparator,levels=setdiff(COMPARATORS,'Random (expected)'))) %>%
  ggplot(aes(comparator,net_extra_hits,fill=comparator))+
  geom_hline(yintercept=0,linetype='dashed',color='#7e8790')+
  geom_col(width=.65)+
  geom_text(aes(label=ifelse(net_extra_hits>=0,paste0('+',net_extra_hits),as.character(net_extra_hits))),
            vjust=ifelse(novelty_summary$net_extra_hits>=0,-.35,1.35),size=3.5)+
  facet_wrap(~nutrient,nrow=1,scales='free_y')+
  labs(title='(c) Added anomaly hits from F-Net-only recommended targets',
       subtitle='Net extra hits = F-Net-only hits minus baseline-only hits; shared hits excluded',
       x=NULL,y='Net additional hit placements across test events')+
  theme_sci+theme(axis.text.x=element_text(angle=18,hjust=1),legend.position='none')
for (nm in c('gain','coverage','novelty')) {
  p <- get(paste0('plot_',nm))
  ggsave(file.path(OUT,paste0('Fig6_Step3b_',nm,'.png')),p,width=10.5,height=5,dpi=350,bg='white')
  ggsave(file.path(OUT,paste0('Fig6_Step3b_',nm,'.pdf')),p,width=10.5,height=5,bg='white')
}
cat('\n========== STEP 3b DONE ==========\n')
cat('Output: ',OUT,'\n',sep='')
cat('\n1. Paired F-Net gain at K=5 (clustered 95% CI):\n')
print(paired_summary %>% filter(metric=='delta_precision'),n=Inf)
cat('\n2. Evidence coverage:\n')
print(coverage_summary,n=Inf)
cat('\n3. Unique hit complementarity:\n')
print(novelty_summary,n=Inf)
cat('\nIMPORTANT: The fitted fusion weights and Global PCS definitions are provisional.\n')
cat('Interpretation requires method-consistency checks; physical unknown != R=0.\n')
