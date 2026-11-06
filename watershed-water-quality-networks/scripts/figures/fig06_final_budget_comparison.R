# Fig.6 final: same-budget comparison of Random/Phys/TE/F-Net follow-up; budget decomposition & QA
# TRAINING network is not refit. Requires existing Step2/Step3/Step3b results.
options(stringsAsFactors=FALSE)
packages <- c('readr','dplyr','tidyr','ggplot2','patchwork','tibble','lubridate','scales')
for(p in packages) if(!requireNamespace(p,quietly=TRUE)) install.packages(p)
invisible(lapply(packages,library,character.only=TRUE))
source("config.R")
ROOT <- get_data_root()
V <- file.path(ROOT,'Fig6_independent_validation')
S2 <- file.path(V,'Step2_event_replay'); S3 <- file.path(V,'Step3_fair_comparison')
S3B <- file.path(V,'Step3b_paired_gain_QA')
TRAIN <- file.path(ROOT,'Fig6_independent_training_2011_2020')
OUT <- file.path(V,'Fig6_final_equal_budget_methods')
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
KSET <- 1:10; K0 <- 5L; B_PHYS <- 500L; B_BUDGET <- 500L; B_BOOT <- 2000L
BUDGETS <- c(8,12,16,24,28); GREEDY_NEIGHBORS <- 5L
GLOBAL_SCORE <- 'PCS_in' # provisional: verify manuscript formula
set.seed(20261029)
readx <- function(p,ids=character()) {
  if(!file.exists(p))stop('Missing: ',p)
  spec<-readr::cols(.default=readr::col_guess());for(n in ids)spec$cols[[n]]<-readr::col_character()
  readr::read_csv(p,col_types=spec,show_col_types=FALSE)
}
need <- function(x,cols,label){m<-setdiff(cols,names(x));if(length(m))stop(label,' missing: ',paste(m,collapse=', '),' | found: ',paste(names(x),collapse=', '))}
id8<-function(x){x<-trimws(sub('\\.0+$','',as.character(x)));j<-!is.na(x)&grepl('^[0-9]{1,8}$',x);x[j]<-sprintf('%08.0f',as.numeric(x[j]));x}
bool<-function(x){if(is.logical(x))return(x);v<-tolower(as.character(x));ifelse(v%in%c('true','t','1'),TRUE,ifelse(v%in%c('false','f','0'),FALSE,NA))}
save<-function(x,f)readr::write_csv(x,file.path(OUT,f),na='')
mean_na<-function(x)if(all(is.na(x)))NA_real_ else mean(x,na.rm=TRUE)

# 1. Locked test universe and previously computed rankings
E<-readx(file.path(S2,'03_independent_test_events.csv'),c('source'))%>%
  transmute(event_id=as.character(event_id),nutrient=toupper(nutrient),source=id8(source),date=as.Date(trigger_date))
L<-readx(file.path(S2,'05_common_event_target_evaluation.csv'),c('source','target'))%>%
  transmute(event_id=as.character(event_id),nutrient=toupper(nutrient),source=id8(source),target=id8(target),evaluable=bool(evaluable),hit=bool(hit_03))
R<-readx(file.path(S3,'02_all_method_ranked_candidates.csv'),c('source','target'))%>%
  mutate(event_id=as.character(event_id),nutrient=toupper(nutrient),source=id8(source),target=id8(target),rank=as.integer(rank),hit_03=bool(hit_03))
M<-readx(file.path(S3,'05_all_event_metrics_including_random.csv'))%>%
  mutate(event_id=as.character(event_id),nutrient=toupper(nutrient),K=as.integer(K))
C<-readx(file.path(S3B,'04_paired_gain_95CI.csv'))
need(R,c('method','event_id','nutrient','source','target','rank','score','strong_score','hit_03'),'R')
need(M,c('method','event_id','nutrient','K','hits','precision','recall','lift','P','N'),'M')
if(anyDuplicated(R[c('method','event_id','target')]))stop('Duplicate ranked candidates')
if(anyDuplicated(E$event_id))stop('Duplicate events')
# Only identical evaluable events previously selected in Step3
eligible<-R%>%distinct(event_id,nutrient,source)
base<-L%>%inner_join(eligible,by=c('event_id','nutrient','source'))%>%filter(evaluable%in%TRUE,!is.na(hit))%>%
  group_by(event_id)%>%mutate(N=n(),P=sum(hit))%>%ungroup()
if(anyDuplicated(base[c('event_id','target')]))stop('Duplicate targets')
qa<-base%>%distinct(event_id,nutrient,source,N,P)
if(any(qa$N<max(KSET)))stop('One or more events have <10 evaluable targets')
if(any(!eligible$event_id%in%qa$event_id))stop('Some Step3 events missing from Step2')
save(qa,'01_event_universe_QA.csv')

# 2. Recompute Phys-Net top K with verified physical paths FIRST,
#    uniform randomized completion among nonverified candidates.
#    These are UNKNOWN reachability, NOT confirmed physical nonconnection.
phys<-R%>%filter(method=='Phys-Net')%>%select(event_id,nutrient,source,target,rank,score,strong_score,hit_03,N,P)
if(nrow(phys)!=nrow(base))stop('Phys candidate universe differs from common universe')
phys_chunks<-split(phys,phys$event_id)
physMC<-vector('list',length(phys_chunks)*B_PHYS)
iout<-0L
for(dat in phys_chunks){
  dat<-dat[order(dat$target),,drop=FALSE]
  verified<-which(bool(dat$strong_score)%in%TRUE)
  verified<-verified[order(-dat$score[verified],dat$target[verified])]
  other<-setdiff(seq_len(nrow(dat)),verified)
  for(b in seq_len(B_PHYS)){
    ord<-c(verified,sample(other,length(other)))
    h<-cumsum(as.integer(dat$hit_03[ord]))
    kk<-KSET
    iout<-iout+1L
    physMC[[iout]]<-tibble(event_id=dat$event_id[1],nutrient=dat$nutrient[1],source=dat$source[1],replicate=b,
                           K=kk,N=dat$N[1],P=dat$P[1],hits=h[kk],precision=h[kk]/kk,
                           recall=if(dat$P[1]>0)h[kk]/dat$P[1] else NA_real_,
                           lift=if(dat$P[1]>0)(h[kk]/kk)/(dat$P[1]/dat$N[1]) else NA_real_,
                           verified_topK=vapply(kk,function(k)sum(seq_len(k)<=length(verified)),integer(1)))
  }
}
physMC<-bind_rows(physMC);save(physMC,'02_phys_random_fallback_replicates.csv')
phys_event<-physMC%>%group_by(event_id,nutrient,K)%>%summarise(
  hits=mean(hits),precision=mean(precision),recall=mean_na(recall),lift=mean_na(lift),
  N=first(N),P=first(P),verified_topK=mean(verified_topK),.groups='drop')%>%mutate(method='Phys-Net (random fallback)')
# Replace old deterministic-ID Phys only for reporting; keep it in QA for reference.
base_metrics<-M%>%filter(method!='Phys-Net')
metrics_all<-bind_rows(base_metrics,phys_event%>%select(all_of(names(base_metrics))))
summary<-metrics_all%>%group_by(nutrient,method,K)%>%summarise(n_events=n(),
                                                               precision=mean_na(precision),recall=mean_na(recall),lift=mean_na(lift),
                                                               pooled_hits=sum(hits),n_recommended=sum(K),.groups='drop')
save(summary,'03_revised_method_curves.csv')
save(summary%>%filter(K==K0),'04_revised_K5_summary.csv')
# A deterministic Phy fallback tie-breaking is now diagnosed, not used as a benchmark.
phys_evidence<-phys%>%filter(rank<=K0)%>%group_by(nutrient)%>%summarise(
  top5_positions=n(),verified=sum(bool(strong_score)%in%TRUE),
  unverified=n()-verified,.groups='drop')
save(phys_evidence,'05_phys_evidence_QA.csv')

# 3. Paired gains: recompute Phys comparison with randomized completion.
#    Source-cluster bootstrap over events; at each event use mean randomized
#    Phys result. Does not include Monte Carlo uncertainty in CI, but B=500
#    reduces it; raw replicates are saved separately.
fnet5<-M%>%filter(method=='Source-conditioned F-Net',K==K0)%>%select(event_id,nutrient,precision_fnet=precision)
phy5<-phys_event%>%filter(K==K0)%>%select(event_id,nutrient,precision_phys=precision)
D<-fnet5%>%inner_join(phy5,by=c('event_id','nutrient'))%>%left_join(E%>%select(event_id,source),by='event_id')%>%
  mutate(delta=precision_fnet-precision_phys)
boot_cluster<-function(dat,B){
  groups<-split(dat$delta,dat$source);n<-length(groups)
  if(n<2)return(c(NA_real_,NA_real_))
  vals<-replicate(B,{s<-sample(seq_len(n),n,replace=TRUE);mean(unlist(groups[s]))})
  unname(quantile(vals,c(.025,.975),na.rm=TRUE))
}
physC<-D%>%group_by(nutrient)%>%group_modify(~{
  ci<-boot_cluster(.x,B_BOOT)
  tibble(comparator='Phys-Net (random fallback)',metric='delta_precision',estimate=mean(.x$delta),ci_low=ci[1],ci_high=ci[2],n_events=nrow(.x),n_sources=n_distinct(.x$source))
})%>%ungroup()
# Other Step3b CIs unchanged. Do not mix OLD Phys-Net CI with revised new Phys-Net.
paired<-C%>%filter(metric=='delta_precision',comparator!='Phys-Net')%>%
  select(any_of(names(physC)))%>%bind_rows(physC)
save(paired,'06_revised_paired_CI.csv')

# 4. Training-only fixed station selection (budget experiment B)
readtrain<-function(nu){
  p<-file.path(TRAIN,nu)
  nodes<-readx(file.path(p,'05_training_global_PCS.csv'),c('STAID'))
  links<-readx(file.path(p,'04_training_source_conditioned_PCS.csv'),c('source','target'))
  need(nodes,c('STAID',GLOBAL_SCORE),'Training Global PCS')
  need(links,c('source','target','source_PCS'),'Training source PCS')
  nd<-nodes%>%transmute(nutrient=nu,station=id8(STAID),value=as.numeric(.data[[GLOBAL_SCORE]]))%>%distinct(station,.keep_all=TRUE)
  lk<-links%>%transmute(source=id8(source),target=id8(target),score=as.numeric(source_PCS))%>%
    filter(source!=target,source%in%nd$station,target%in%nd$station,is.finite(score),score>0)%>%
    arrange(source,desc(score),target)%>%group_by(source)%>%slice_head(n=GREEDY_NEIGHBORS)%>%ungroup()
  list(nodes=nd,links=lk)
}
tr<-lapply(c('TN','TP'),readtrain)
greedy_order<-function(nodes,links){
  bysrc<-split(links$target,links$source);chosen<-character();covered<-character()
  while(length(chosen)<nrow(nodes)){
    rem<-setdiff(nodes$station,chosen)
    gains<-vapply(rem,function(s)length(setdiff(bysrc[[s]],covered)),integer(1))
    val<-nodes$value[match(rem,nodes$station)];val[!is.finite(val)]<--Inf
    winner<-rem[order(-gains,-val,rem)][1]
    chosen<-c(chosen,winner);covered<-union(covered,bysrc[[winner]])
  }
  chosen
}
selections<-list();m<-0L
for(nu in c('TN','TP')){
  tt<-tr[[match(nu,c('TN','TP'))]]
  nodes<-tt$nodes; nmax<-nrow(nodes)
  if(max(BUDGETS)>nmax)stop('Budget larger than training nodes')
  pcs<-nodes%>%arrange(desc(value),station)%>%pull(station)
  greedy<-greedy_order(nodes,tt$links)
  rand<-replicate(B_BUDGET,sample(nodes$station),simplify=FALSE)
  for(budget in BUDGETS)for(strategy in c('PCS-topN','Greedy influence coverage','Random')){
    samples<-switch(strategy,'PCS-topN'=list(pcs),'Greedy influence coverage'=list(greedy),'Random'=rand)
    for(b in seq_along(samples)){
      m<-m+1L;selections[[m]]<-tibble(nutrient=nu,budget=budget,strategy=strategy,replicate=b,station=head(samples[[b]],budget))
    }
  }
}
selected<-bind_rows(selections);save(selected,'07_budget_selected_stations.csv')
# Full test-event universe, same across all station strategies.
# Event-triggered Top5 hits are retrospective follow-up opportunities.
F5<-M%>%filter(method=='Source-conditioned F-Net',K==K0)%>%select(event_id,nutrient,hits,precision,recall,P)
events<-E%>%inner_join(F5,by=c('event_id','nutrient'))
if(anyDuplicated(events$event_id))stop('Budget event IDs duplicated')
all_tot<-events%>%group_by(nutrient)%>%summarise(all_events=n(),all_positive_pairs=sum(P),.groups='drop')
keys<-selected%>%distinct(nutrient,budget,strategy,replicate,station)
triggers<-events%>%inner_join(keys,by=c('nutrient','source'='station'),relationship='many-to-many')
bytrial<-triggers%>%group_by(nutrient,budget,strategy,replicate)%>%summarise(
  triggered_events=n(),hit_pairs=sum(hits),
  conditional_precision=mean(precision),
  conditional_recall=if(sum(P)>0)sum(hits)/sum(P) else NA_real_,
  positive_pairs_in_triggered=sum(P),.groups='drop')
budget<-keys%>%distinct(nutrient,budget,strategy,replicate)%>%
  left_join(bytrial,by=c('nutrient','budget','strategy','replicate'))%>%
  left_join(all_tot,by='nutrient')%>%
  mutate(triggered_events=coalesce(triggered_events,0L),hit_pairs=coalesce(hit_pairs,0),
         positive_pairs_in_triggered=coalesce(positive_pairs_in_triggered,0),
         trigger_coverage=triggered_events/all_events,
         conditional_precision=if_else(triggered_events>0,conditional_precision,NA_real_),
         conditional_recall=if_else(triggered_events>0,conditional_recall,NA_real_),
         overall_opportunity=if_else(all_positive_pairs>0,hit_pairs/all_positive_pairs,NA_real_),
         identity_ok=abs(overall_opportunity-if_else(all_positive_pairs>0,
                                                     (positive_pairs_in_triggered/all_positive_pairs)*coalesce(conditional_recall,0),NA_real_))<1e-10)
if(any(!budget$identity_ok,na.rm=TRUE))stop('Budget decomposition identity failed')
save(budget,'08_budget_replicates_decomposed.csv')
summaryB<-budget%>%group_by(nutrient,budget,strategy)%>%summarise(
  trigger=mean(trigger_coverage),trigger_lo=quantile(trigger_coverage,.025),trigger_hi=quantile(trigger_coverage,.975),
  conditional_precision=mean_na(conditional_precision),conditional_recall=mean_na(conditional_recall),
  opportunity=mean(overall_opportunity),opportunity_lo=quantile(overall_opportunity,.025),opportunity_hi=quantile(overall_opportunity,.975),
  n_rep=n(),.groups='drop')
save(summaryB,'09_budget_summary_decomposed.csv')

# 4B. Same fixed monitoring stations, DIFFERENT follow-up ranking methods
#     Fig 6(e): each scheme/replicate has exactly the same trigger events
#     for Random, Phys-Net, TE-Net and Source-conditioned F-Net.
#     Score: sum of event-target anomaly hits in Top-5 across triggered
#     events divided by all event-target positive outcomes across the
#     entire common evaluation universe (includes untriggered events).
#     Historical retrospective follow-up opportunity, NOT observed detection.
METHODS_E <- c('Random (expected)', 'Phys-Net (random fallback)',
               'TE-Net', 'Source-conditioned F-Net')
BUDGET_STRATEGY_MAIN <- 'PCS-topN'  # fixed for the main Fig.6(e)
# Other two fixed-station strategies are retained as supplementary comparisons.

method5 <- metrics_all %>%
  dplyr::filter(K == K0, method %in% METHODS_E) %>%
  dplyr::select(method, event_id, nutrient, hits, precision, recall, P, N)
if (anyDuplicated(method5[c('method','event_id')])) stop('Duplicated method-event at K=5')
if (any(!is.finite(method5$hits))) stop('Missing Top-5 hits in method comparisons')
if (any(method5$hits < 0 | method5$hits > K0)) stop('Impossible hit count')
if (any(method5$P < 0 | method5$P > method5$N)) stop('Impossible target positives')

# Common event universe must be identical between all methods.
method_qa <- method5 %>%
  dplyr::group_by(nutrient, method) %>%
  dplyr::summarise(n_events=dplyr::n_distinct(event_id), .groups='drop')
for (nu in c('TN','TP')) {
  rows <- method_qa %>% dplyr::filter(nutrient == nu)
  if (nrow(rows) != length(METHODS_E)) stop('Missing follow-up method for ',nu)
  if (length(unique(rows$n_events)) != 1L) stop('Method event universe mismatch for ',nu)
  expected_ids <- sort(unique(events$event_id[events$nutrient==nu]))
  for (meth in METHODS_E) {
    actual_ids <- sort(unique(method5$event_id[method5$nutrient==nu & method5$method==meth]))
    if (!identical(actual_ids,expected_ids)) stop('Event IDs mismatch: ',nu,' / ',meth)
  }
}

# Assign a unique key to each fixed-station selection trial.
trial_keys <- selected %>%
  dplyr::distinct(nutrient,budget,strategy,replicate)
# Crucial: join triggers ONCE before joining any method ranking.
trigger_trials <- events %>%
  dplyr::select(event_id,nutrient,source) %>%
  dplyr::inner_join(selected, by=c('nutrient','source'='station'),relationship='many-to-many') %>%
  dplyr::select(event_id,nutrient,budget,strategy,replicate)
if (anyDuplicated(trigger_trials[c('event_id','nutrient','budget','strategy','replicate')])) {
  stop('Repeated trigger event within the same fixed-station trial')
}

# Each method sees exactly the same triggered event IDs in each trial.
trial_rankings <- trigger_trials %>%
  dplyr::inner_join(method5, by=c('event_id','nutrient'),relationship='many-to-many')

# Aggregate at one selection-trial x one method.
trial_perf <- trial_rankings %>%
  dplyr::group_by(nutrient,budget,strategy,replicate,method) %>%
  dplyr::summarise(n_triggered=dplyr::n_distinct(event_id),
                   hits_top5=sum(hits),
                   total_pos_triggered=sum(P),
                   recommended_positions=sum(rep(K0,dplyr::n())),
                   .groups='drop')

# Add zero-trigger trials by explicitly expanding method labels.
all_trials_methods <- tidyr::crossing(trial_keys,method=METHODS_E)
trial_perf <- all_trials_methods %>%
  dplyr::left_join(trial_perf,
                   by=c('nutrient','budget','strategy','replicate','method')) %>%
  dplyr::left_join(all_tot,by='nutrient') %>%
  dplyr::mutate(n_triggered=dplyr::coalesce(n_triggered,0L),
                hits_top5=dplyr::coalesce(hits_top5,0),
                total_pos_triggered=dplyr::coalesce(total_pos_triggered,0),
                recommended_positions=dplyr::coalesce(recommended_positions,0L),
                trigger_coverage=n_triggered/all_events,
                conditional_recall=dplyr::if_else(total_pos_triggered>0,
                                                  hits_top5/total_pos_triggered,NA_real_),
                overall_opportunity=dplyr::if_else(all_positive_pairs>0,
                                                   hits_top5/all_positive_pairs,NA_real_))

# Same fixed station set => trigger counts and denominators must match
# between ALL methods in each trial.
trial_common_QA <- trial_perf %>%
  dplyr::group_by(nutrient,budget,strategy,replicate) %>%
  dplyr::summarise(n_methods=dplyr::n(),
                   same_trigger=(dplyr::n_distinct(n_triggered)==1L),
                   same_denominator=(dplyr::n_distinct(all_positive_pairs)==1L),
                   same_triggered_positives=(dplyr::n_distinct(total_pos_triggered)==1L),
                   .groups='drop')
if (any(trial_common_QA$n_methods!=length(METHODS_E) |
        !trial_common_QA$same_trigger |
        !trial_common_QA$same_denominator |
        !trial_common_QA$same_triggered_positives))
  stop('Budget methods DO NOT share identical fixed-station triggers / outcome denominator')
save(trial_common_QA,'11_equal_budget_comparability_QA.csv')
save(trial_perf,'12_equal_budget_followup_replicates.csv')

# Summarize each method under each fixed-station selection strategy.
method_budget <- trial_perf %>%
  dplyr::group_by(nutrient,budget,strategy,method) %>%
  dplyr::summarise(n_trials=dplyr::n(),
                   overall_opportunity=mean(overall_opportunity),
                   overall_lo=stats::quantile(overall_opportunity,.025,na.rm=TRUE),
                   overall_hi=stats::quantile(overall_opportunity,.975,na.rm=TRUE),
                   conditional_recall=mean_na(conditional_recall),
                   trigger_coverage=mean(trigger_coverage),
                   .groups='drop')
save(method_budget,'13_equal_budget_followup_summary.csv')

# Paired difference against random follow-up, using the SAME trial keys.
rand_budget <- trial_perf %>%
  dplyr::filter(method=='Random (expected)') %>%
  dplyr::select(nutrient,budget,strategy,replicate,
                random_opportunity=overall_opportunity)
paired_budget <- trial_perf %>%
  dplyr::left_join(rand_budget,
                   by=c('nutrient','budget','strategy','replicate')) %>%
  dplyr::mutate(delta_vs_random=overall_opportunity-random_opportunity)
save(paired_budget,'14_equal_budget_paired_vs_random.csv')
cat('\nMethod comparison under identical fixed-station triggers: QA PASSED\n')

# 5. Event stability (NO interpolation of empty year-season cells)
stability<-M%>%filter(method=='Source-conditioned F-Net',K==K0)%>%left_join(E%>%select(event_id,date),by='event_id')%>%
  mutate(year=lubridate::year(date),season=case_when(lubridate::month(date)%in%c(12,1,2)~'DJF',
                                                     lubridate::month(date)%in%3:5~'MAM',lubridate::month(date)%in%6:8~'JJA',TRUE~'SON'))%>%
  group_by(nutrient,year,season)%>%summarise(n_events=n(),precision=mean(precision),recall=mean_na(recall),.groups='drop')%>%
  pivot_longer(c(precision,recall),names_to='metric',values_to='value')%>%
  mutate(metric=factor(metric,levels=c('precision','recall'),labels=c('Precision@5','Recall@5')),
         season=factor(season,levels=c('DJF','MAM','JJA','SON')))
save(stability,'10_year_season_stability.csv')

# 6. Plot
method_cols<-c('Source-conditioned F-Net'='#D7733B','Phys-Net (random fallback)'='#3475AB',
               'TE-Net'='#9366B2','Global PCS'='#3A9369','Random (expected)'='#8A9499')
budget_cols<-c('PCS-topN'='#3475AB','Greedy influence coverage'='#D7733B','Random'='#8A9499')
th<-theme_bw(base_size=9)+theme(panel.grid.minor=element_blank(),strip.background=element_rect(fill='#F1F3F4',color='#D4DADF'),
                                strip.text=element_text(face='bold',size=8),plot.title=element_text(face='bold',size=10),
                                plot.subtitle=element_text(size=8,color='#66737D'),legend.title=element_blank(),legend.position='bottom',
                                axis.text=element_text(size=7),plot.margin=margin(6,6,6,6))
summary$method<-factor(summary$method,levels=names(method_cols))
p_a<-ggplot(summary,aes(K,precision,color=method,group=method))+geom_line(linewidth=.8)+geom_point(size=1.1)+
  facet_wrap(~nutrient,nrow=1)+scale_color_manual(values=method_cols)+scale_y_continuous(labels=scales::label_percent())+
  labs(title='(a) Independent-test Precision@K',subtitle='Phys-Net completion randomized 500 times',x='Recommended targets (K)',y='Precision@K')+th
p_b<-summary%>%select(nutrient,method,K,recall,lift)%>%pivot_longer(c(recall,lift),names_to='metric',values_to='value')%>%
  mutate(metric=factor(metric,levels=c('recall','lift'),labels=c('Recall@K','Lift@K')))%>%
  ggplot(aes(K,value,color=method,group=method))+geom_line(linewidth=.8)+geom_point(size=1.1)+
  facet_grid(metric~nutrient,scales='free_y')+scale_color_manual(values=method_cols)+
  labs(title='(b) Recall@K and Lift@K',x='K',y=NULL)+th
paired$nutrient<-factor(paired$nutrient,levels=c('TN','TP'))
p_c<-ggplot(paired,aes(y=comparator,x=100*estimate))+geom_vline(xintercept=0,linetype='dashed',color='#A0ADB2')+
  geom_segment(aes(x=100*ci_low,xend=100*ci_high,y=comparator,yend=comparator),color='#637E8B',linewidth=.7)+
  geom_point(color='#D7733B',size=2.4)+facet_wrap(~nutrient,nrow=1)+
  labs(title='(c) Paired gain in F-Net Precision@5',subtitle='Source-cluster bootstrap; revised Phys-Net fallback',
       x='F-Net minus baseline (percentage points)',y=NULL)+th
p_d<-ggplot(summaryB,aes(budget,trigger,color=strategy,group=strategy))+
  geom_ribbon(data=summaryB%>%filter(strategy=='Random'),aes(ymin=trigger_lo,ymax=trigger_hi,fill=strategy),alpha=.12,color=NA)+
  geom_line(linewidth=.8)+geom_point(size=1.7)+facet_wrap(~nutrient,nrow=1)+
  scale_color_manual(values=budget_cols)+scale_fill_manual(values=budget_cols)+scale_x_continuous(breaks=BUDGETS)+
  scale_y_continuous(labels=scales::label_percent())+
  labs(title='(d) Trigger coverage by fixed station budget',subtitle='Random band: sampling variation, not a 95% CI of model performance',
       x='Fixed stations (N)',y='Triggered source-event share')+th
# Final panel (e): follow-up METHOD comparison with EXACT SAME
# fixed stations and trigger events within each budget.
# Main uses PCS-topN fixed stations; Greedy/Random station schemes
# are saved in supplementary tables (13/14) for robustness.
Eplot <- method_budget %>% dplyr::filter(strategy==BUDGET_STRATEGY_MAIN)
if(!nrow(Eplot))stop('No budget simulations for main fixed-station strategy')
Eplot$method <- factor(Eplot$method,levels=METHODS_E)
p_e <- ggplot(Eplot,aes(x=budget,y=overall_opportunity,color=method,group=method)) +
  geom_line(linewidth=.9) + geom_point(size=1.7) +
  facet_wrap(~nutrient,nrow=1) +
  scale_color_manual(values=method_cols,drop=TRUE) +
  scale_x_continuous(breaks=BUDGETS) +
  scale_y_continuous(labels=scales::label_percent(accuracy=.1)) +
  labs(title='(e) Equal-budget follow-up methods',
       subtitle=paste0('Identical PCS-topN sensors and events; K=',K0,
                       ' | Greedy/Random sensor layouts in supplementary tables'),
       x='Fixed stations (N)',
       y='Overall retrospective follow-up opportunity') + th
p_f<-ggplot(stability,aes(season,factor(year),fill=value))+geom_tile(color='white')+
  geom_text(aes(label=paste0(sprintf('%.2f',value),'\nn=',n_events)),size=2.3)+
  facet_grid(metric~nutrient)+scale_fill_viridis_c(option='C',na.value='#EDF0F1')+
  labs(title='(f) Year–season stability of F-Net',subtitle='Small event counts imply high uncertainty',x='Season',y='Test year',fill='Value')+th+
  theme(legend.position='right')
fig<-((p_a|p_b)/(p_c|p_d)/(p_e|p_f))+patchwork::plot_layout(heights=c(1,1,1),widths=c(1,1.2))+
  plot_annotation(title='Fig. 6 | Independent-test anomaly screening and budget-constrained follow-up',
                  subtitle='Training: 2011–2020; test replay: 2021–2023 | 0–3-month outcomes; fixed K=5 for budget simulation',
                  caption=paste('Phys-Net: verified directed paths first, otherwise random completion; unknown paths are not deemed disconnected.',
                                'Budget results describe retrospective event-target opportunities, not actual detection without additional sampling.',
                                'F-Net fusion and Global PCS definitions remain provisional until matched to manuscript.'),
                  theme=theme(plot.title=element_text(face='bold',size=14),plot.subtitle=element_text(size=9),
                              plot.caption=element_text(size=8,hjust=0,color='#66737D')))
ggsave(file.path(OUT,'Fig6_final_equal_budget.pdf'),fig,width=18,height=18,bg='white')
ggsave(file.path(OUT,'Fig6_final_equal_budget.png'),fig,width=18,height=18,dpi=350,bg='white')
for(nm in letters[1:6])ggsave(file.path(OUT,paste0('Fig6_',nm,'.pdf')),get(paste0('p_',nm)),width=10,height=6,bg='white')
cat('\nRevised K=5:\n');print(summary%>%filter(K==5),n=Inf)
cat('\nPhys evidence:\n');print(phys_evidence,n=Inf)
cat('\nBudget summary:\n');print(summaryB,n=Inf)
cat('\nEqual-budget method summary (main panel e):\n');print(Eplot,n=Inf)
cat('\nOutput: ',OUT,'\n',sep='')
print(fig)
