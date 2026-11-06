# ================================================================
# Fig.6 Step 3 | Fair independent anomaly-replay evaluation
# Training networks: 2011-10--2020-12; test events: 2021--2023
# Same evaluable targets and same positive set for every method.
# ================================================================
options(stringsAsFactors = FALSE)
pkgs <- c('readr','dplyr','tidyr','ggplot2','tibble')
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
invisible(lapply(pkgs, library, character.only=TRUE))
source("config.R")
ROOT <- get_data_root()
TRAIN <- file.path(ROOT,'Fig6_independent_training_2011_2020')
STEP2 <- file.path(ROOT,'Fig6_independent_validation','Step2_event_replay')
PNET <- file.path(ROOT,'PNet_2011_2023_BlockTracing')
OUT <- file.path(ROOT,'Fig6_independent_validation','Step3_fair_comparison')
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
FOCAL_K <- 5L
K_SEQ <- 1:10
# No random sampling required: exact expected metrics for uniform random ranking.
# If Phys-Net lacks a verified path, assign zero priority, not evidence of R=0.
# Global PCS choice is provisional. Must match paper definition before formal comparison.
GLOBAL_COL <- 'PCS_in'
TE_PMAX <- 0.10

id8 <- function(x) {
  x <- trimws(sub('\\.0+$','',as.character(x)))
  ii <- !is.na(x) & grepl('^[0-9]{1,8}$',x)
  x[ii] <- sprintf('%08.0f',as.numeric(x[ii]))
  x
}
getcsv <- function(path, ids=character()) {
  if (!file.exists(path)) stop('Missing input: ',path)
  spec <- readr::cols(.default=readr::col_guess())
  for (nm in ids) spec$cols[[nm]] <- readr::col_character()
  readr::read_csv(path,col_types=spec,show_col_types=FALSE)
}
needs <- function(d,nms,label) {
  miss <- setdiff(nms,names(d))
  if (length(miss)) stop(label,' missing columns: ',paste(miss,collapse=', '),
                         '\nActual: ',paste(names(d),collapse=', '))
}
writeout <- function(x,filename) readr::write_csv(x,file.path(OUT,filename),na='')

# ---------------------------------------------------------------
# 1. Frozen TEST outcome table and selected test events
# ---------------------------------------------------------------
events <- getcsv(file.path(STEP2,'03_independent_test_events.csv'),c('source')) %>%
  transmute(event_id=as.character(event_id),
            nutrient=toupper(as.character(nutrient)),
            source=id8(source), trigger_date=as.Date(trigger_date))
if (anyDuplicated(events$event_id)) stop('Duplicated test event IDs.')
obs0 <- getcsv(file.path(STEP2,'05_common_event_target_evaluation.csv'),c('source','target'))
needs(obs0,c('event_id','nutrient','source','target','evaluable','hit_03'), 'Step2 outcomes')
tobool <- function(z) {
  if (is.logical(z)) return(z)
  zz <- tolower(trimws(as.character(z)))
  ifelse(is.na(zz),NA, ifelse(zz %in% c('true','t','1','yes'),TRUE,
                              ifelse(zz %in% c('false','f','0','no'),FALSE,NA)))
}
obs <- obs0 %>%
  transmute(event_id=as.character(event_id), nutrient=toupper(as.character(nutrient)),
            source=id8(source),target=id8(target),
            evaluable=tobool(evaluable), hit_03=tobool(hit_03)) %>%
  inner_join(events,by=c('event_id','nutrient','source')) %>%
  filter(source != target)
if (anyDuplicated(obs[c('event_id','target')])) stop('Duplicate event-target rows in Step2.')
if (any(obs$evaluable & is.na(obs$hit_03),na.rm=TRUE)) stop('Evaluable targets have missing labels.')
# Primary metric: complete observations only. All methods share this universe.
U <- obs %>% filter(evaluable %in% TRUE, !is.na(hit_03)) %>%
  group_by(event_id) %>% mutate(N=n(),P=sum(hit_03)) %>% ungroup()
qa_events <- U %>% distinct(event_id,nutrient,source,trigger_date,N,P) %>%
  mutate(eligible = N >= max(K_SEQ))
writeout(qa_events,'01_event_eligibility_QA.csv')
eligible_ids <- qa_events %>% filter(eligible) %>% pull(event_id)
U <- U %>% filter(event_id %in% eligible_ids)
if (!nrow(U)) stop('No eligible test events with >= ',max(K_SEQ),' evaluable targets.')
message('Eligible events: ',dplyr::n_distinct(U$event_id),
        '; complete evaluable event-target pairs: ',nrow(U))

# ---------------------------------------------------------------
# 2. Read TRAINING scores with strict key/column checks
# ---------------------------------------------------------------
read_training <- function(nu) {
  dir <- file.path(TRAIN,nu)
  fnet <- getcsv(file.path(dir,'04_training_source_conditioned_PCS.csv'),c('source','target'))
  needs(fnet,c('source','target','source_PCS'),'Training F-Net')
  fnet <- fnet %>% transmute(nutrient=nu,source=id8(source),target=id8(target),
                             score=as.numeric(source_PCS)) %>%
    distinct(nutrient,source,target,.keep_all=TRUE)
  glob <- getcsv(file.path(dir,'05_training_global_PCS.csv'),c('STAID'))
  needs(glob,c('STAID',GLOBAL_COL),'Training Global PCS')
  global <- tibble(nutrient=nu,target=id8(glob$STAID),
                   score=as.numeric(glob[[GLOBAL_COL]])) %>%
    distinct(nutrient,target,.keep_all=TRUE)
  te <- getcsv(file.path(dir,'train_2011_2020_all_edges.csv'),c('from','to'))
  needs(te,c('from','to','p_value','TE_bias_corrected'),'Training TE')
  te <- te %>% transmute(nutrient=nu,source=id8(from),target=id8(to),
                         score=ifelse(is.finite(p_value) & p_value < TE_PMAX,
                                      pmax(as.numeric(TE_bias_corrected),0),0)) %>%
    group_by(nutrient,source,target) %>% summarise(score=max(score,na.rm=TRUE),.groups='drop')
  phy <- getcsv(file.path(PNET,nu,'PNet_directed_reachable_pairs.csv'),c('from','to'))
  needs(phy,c('from','to'),'Physical reachability')
  # Explicitly read slow / medium / fast if present. Do not recompute from arbitrary speed.
  tc <- c('travel_days_medium','travel_days_slow','travel_days_fast')
  col <- tc[tc %in% names(phy)][1]
  if (length(col)==0 || is.na(col)) stop(nu,': no measured/derived PNet travel-days column. ',paste(names(phy),collapse=', '))
  physical <- tibble(nutrient=nu,source=id8(phy$from),target=id8(phy$to),
                     travel_days=as.numeric(phy[[col]])) %>%
    filter(is.finite(travel_days),travel_days>=0) %>%
    group_by(nutrient,source,target) %>%
    summarise(score=-min(travel_days),verified_path=TRUE,.groups='drop')
  list(fnet=fnet,global=global,te=te,physical=physical)
}
tr <- lapply(c('TN','TP'),read_training)
F <- bind_rows(lapply(tr,`[[`,'fnet'))
G <- bind_rows(lapply(tr,`[[`,'global'))
T <- bind_rows(lapply(tr,`[[`,'te'))
PH <- bind_rows(lapply(tr,`[[`,'physical'))

# ---------------------------------------------------------------
# 3. Full rankings. Every method has exactly N candidates/event.
#    Unscored/absent edges rank LAST with a deterministic ID tie-break.
#    They are marked explicitly; missing physical edges are NOT R=0.
# ---------------------------------------------------------------
base <- U %>% select(event_id,nutrient,source,target,hit_03,N,P)
mk <- function(data,method,keys,positive_method=TRUE) {
  ans <- base %>% left_join(data,by=keys) %>%
    mutate(method=method,
           score_available=is.finite(score),
           rank_score=ifelse(score_available,score,-Inf),
           # 'strong' is a diagnostic, not part of the common denominator.
           strong_score=score_available & (if (positive_method) score>0 else TRUE)) %>%
    group_by(event_id) %>% arrange(desc(rank_score),target,.by_group=TRUE) %>%
    mutate(rank=row_number()) %>% ungroup() %>%
    select(event_id,nutrient,source,target,hit_03,N,P,method,
           score,score_available,strong_score,rank)
  if (nrow(ans)!=nrow(base)) stop('Unexpected join multiplicity: ',method)
  ans
}
rankings <- bind_rows(
  mk(F,'Source-conditioned F-Net',c('nutrient','source','target')),
  mk(T,'TE-Net',c('nutrient','source','target')),
  mk(G,'Global PCS',c('nutrient','target')),
  mk(PH %>% select(-verified_path),'Phys-Net',c('nutrient','source','target'),positive_method=FALSE)
)
writeout(rankings,'02_all_method_ranked_candidates.csv')
coverage <- rankings %>% group_by(method,nutrient,event_id) %>%
  summarise(N=first(N),n_scores=sum(score_available),n_strong=sum(strong_score),
            n_positive_targets=sum(hit_03),.groups='drop')
writeout(coverage,'03_ranking_coverage_QA.csv')

# ---------------------------------------------------------------
# 4. Calculate event-level metrics with COMMON denominator P.
#    P=0 => Recall and Lift undefined; precision remains defined.
#    Lift = Precision / base rate (P/N) for that event.
# ---------------------------------------------------------------
by_event <- rankings %>% group_by(method,nutrient,event_id) %>%
  group_modify(~ {
    z <- .x[order(.x$rank),,drop=FALSE]
    N <- z$N[1]; P <- z$P[1]
    bind_rows(lapply(K_SEQ,function(k) {
      h <- sum(z$hit_03[seq_len(k)])
      precision <- h/k
      tibble(K=k,N=N,P=P,hits=h,precision=precision,
             recall=if(P>0) h/P else NA_real_,
             lift=if(P>0) precision/(P/N) else NA_real_)
    }))
  }) %>% ungroup()
writeout(by_event,'04_event_level_metrics.csv')

# Uniform random expected values, exact conditional expectation.
# A random Top-K sample has E[hits]=K*(P/N), hence E[Precision]=P/N,
# E[Recall]=K/N (P>0), E[Lift]=1 (P>0).
random_events <- qa_events %>% filter(eligible) %>%
  tidyr::crossing(K=K_SEQ) %>%
  transmute(method='Random (expected)', nutrient,event_id,K,N,P,
            hits=K*P/N,precision=P/N,
            recall=ifelse(P>0,K/N,NA_real_),
            lift=ifelse(P>0,1,NA_real_))
metrics_all <- bind_rows(by_event,random_events)
writeout(metrics_all,'05_all_event_metrics_including_random.csv')

# Macro-average across events. Include P=0 in Precision, but not Recall/Lift.
mean_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x,na.rm=TRUE)
summary <- metrics_all %>% group_by(nutrient,method,K) %>%
  summarise(n_events=n(),n_events_with_anomaly=sum(P>0),
            precision=mean_na(precision),recall=mean_na(recall),
            lift=mean_na(lift),
            pooled_hits=sum(hits),pooled_recommendations=sum(K),
            micro_precision=pooled_hits/pooled_recommendations,
            .groups='drop')
writeout(summary,'06_topK_method_summary.csv')
focal <- summary %>% filter(K==FOCAL_K) %>%
  arrange(nutrient,desc(precision),method)
writeout(focal,'07_focalK5_summary.csv')

# ---------------------------------------------------------------
# 5. Plot panel a and b (no false CIs).
# ---------------------------------------------------------------
cols <- c('Random (expected)'='#999999','Phys-Net'='#3578B3',
          'TE-Net'='#9467BD','Global PCS'='#3B9264',
          'Source-conditioned F-Net'='#E07B39')
th <- theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),
                                   strip.background=element_rect(fill='#F2F4F5',color='#CCD3D8'),
                                   legend.position='bottom')
p_a <- ggplot(summary,aes(K,precision,color=method,group=method))+
  geom_line(linewidth=.9)+geom_point(size=1.6)+
  facet_wrap(~nutrient,nrow=1)+scale_color_manual(values=cols)+
  scale_x_continuous(breaks=K_SEQ)+scale_y_continuous(labels=scales::label_percent())+
  labs(title='(a) Independent-test Precision@K',x='Number of recommended stations (K)',
       y='Mean Precision@K',color=NULL)+th
p_b <- summary %>% select(nutrient,method,K,recall,lift) %>%
  pivot_longer(c(recall,lift),names_to='metric',values_to='value') %>%
  ggplot(aes(K,value,color=method,group=method))+
  geom_line(linewidth=.9)+geom_point(size=1.6)+
  facet_grid(metric~nutrient,scales='free_y')+scale_color_manual(values=cols)+
  scale_x_continuous(breaks=K_SEQ)+
  labs(title='(b) Independent-test Recall@K and Lift@K',x='K',y=NULL,color=NULL)+th
for (nm in c('a','b')) {
  p <- get(paste0('p_',nm))
  ggsave(file.path(OUT,paste0('Fig6_panel_',nm,'.png')),p,
         width=11,height=6,dpi=350,bg='white')
  ggsave(file.path(OUT,paste0('Fig6_panel_',nm,'.pdf')),p,
         width=11,height=6,bg='white')
}

# ---------------------------------------------------------------
# 6. Print clear QA and focal result
# ---------------------------------------------------------------
cat('\n================ DATA QA ================\n')
print(qa_events %>% group_by(nutrient) %>% summarise(
  total_events=n(),eligible_events=sum(eligible),
  mean_evaluable_targets=mean(N),events_with_positive=sum(P>0),.groups='drop'))
cat('\n============= K=5 RESULTS ==============\n')
print(focal %>% select(nutrient,method,n_events,n_events_with_anomaly,
                       precision,recall,lift,micro_precision),n=Inf,width=Inf)
cat('\n========== RANKING COVERAGE ============\n')
print(coverage %>% group_by(nutrient,method) %>% summarise(
  avg_scored=mean(n_scores),avg_strong=mean(n_strong),.groups='drop'),n=Inf)
cat('\nOutput directory: ',OUT,'\n',sep='')
print(p_a)
print(p_b)
