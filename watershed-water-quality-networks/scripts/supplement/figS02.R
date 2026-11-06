source("config.R")
# Figure S2 | Upstream-downstream information and dual-lag characteristics
# 2011-2023 full-period TN and TP. Run in RStudio on the original I: drive.
# Inputs: full-pair TE from Fusion_Net_2011_2023, PNet reachability,
#         PNet directed reachable pairs, monthly NTN concentrations.
# No TE re-estimation, no FDR, no silent imputation of unknown physical status.

packages <- c('readr','dplyr','tidyr','ggplot2','patchwork','stringr','tibble','scales')
for(p in packages) if(!requireNamespace(p,quietly=TRUE)) stop('Install package: ',p)
invisible(lapply(packages,library,character.only=TRUE))
ROOT <- get_data_root()
FR <- file.path(ROOT,'Fusion_Net_2011_2023')
PR <- file.path(ROOT,'PNet_2011_2023_BlockTracing')
DATA_FILE <- file.path(ROOT,'NTN_2023_MonthlyLoadTable.csv')
OUT <- file.path(ROOT,'Figure_S2_fullperiod')
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
ALPHA <- .10
WINDOW <- 6L # months before and after source event
COL <- c(TN='#327FA6',TP='#E18652')
EXCLUDED_FIG3 <- c('01529500>01578310','01562000>01578310',
                   '01599000>01646580','01631000>01646580',
                   '02034000>02035000','02024752>02035000')

check_file <- function(x) {if(!file.exists(x))stop('Input missing: ',x,call.=FALSE);x}
id8 <- function(x) {x<-trimws(as.character(x));x<-sub('\\.0+$','',x);x[x %in% c('','NA','NaN')]<-NA_character_;good<-!is.na(x)&grepl('^[0-9]{1,8}$',x);x[good]<-stringr::str_pad(x[good],8,pad='0');x}
colpick <- function(d,opts,required=TRUE){z<-opts[opts %in% names(d)];if(length(z))return(z[1]);if(required)stop('Missing columns ',paste(opts,collapse='/'),' available: ',paste(names(d),collapse=', '));NA_character_}
checkpair <- function(d,lab){if(anyNA(d$from)|anyNA(d$to)|any(d$from==d$to)|anyDuplicated(d[c('from','to')]))stop('Invalid/duplicate pairs in ',lab)}

read_reach <- function(nu,ids){
  p<-file.path(PR,nu,'PNet_reachability_matrix.csv')
  d<-readr::read_csv(check_file(p),col_types=readr::cols(.default=readr::col_character()),show_col_types=FALSE)
  f<-colpick(d,c('from','source','upstream'),FALSE)
  t<-colpick(d,c('to','target','downstream'),FALSE)
  r<-colpick(d,c('R','reachable','reachability','is_reachable','R_phys'),FALSE)
  if(!is.na(f)&&!is.na(t)&&!is.na(r)) {
    z<-tibble::tibble(from=id8(d[[f]]),to=id8(d[[t]]),raw=d[[r]])
  } else {
    z<-d |> dplyr::rename(from=1) |> tidyr::pivot_longer(-from,names_to='to',values_to='raw') |>
      dplyr::mutate(from=id8(from),to=id8(to))
  }
  v<-tolower(trimws(as.character(z$raw)))
  z$R<-dplyr::case_when(v %in% c('1','1.0','true','yes','reachable')~1L,
                        v %in% c('0','0.0','false','no','unreachable','non-reachable')~0L,
                        TRUE~NA_integer_)
  z<-z |> dplyr::filter(from %in% ids,to %in% ids,from!=to) |> dplyr::select(from,to,R)
  if(anyDuplicated(z[c('from','to')]))stop(nu,' reachability has duplicate pairs')
  z
}

read_nu <- function(nu){
  p<-file.path(FR,nu,'Fusion_all_pairs_layers.csv')
  d<-readr::read_csv(check_file(p),col_types=readr::cols(from=readr::col_character(),to=readr::col_character(),.default=readr::col_guess()),show_col_types=FALSE)
  needed<-c('from','to','TE_raw','p_value','lag_months')
  if(length(setdiff(needed,names(d))))stop(nu,' missing: ',paste(setdiff(needed,names(d)),collapse=', '))
  z<-d |> dplyr::transmute(from=id8(from),to=id8(to),TE=as.numeric(TE_raw),p=as.numeric(p_value),lag=as.numeric(lag_months))
  checkpair(z,paste(nu,'TE'))
  ids<-sort(unique(c(z$from,z$to)))
  if(nrow(z)!=length(ids)*(length(ids)-1L)) stop(nu,' TE is not full directed candidate-pair table; reverse-direction S2(a) unsafe')
  if(any(!is.finite(z$TE))|anyNA(z$p)|any(z$p<0|z$p>1))stop(nu,' nonfinite TE / invalid p-values')
  z$TE<-pmax(z$TE,0)
  r<-read_reach(nu,ids)
  z<-dplyr::left_join(z,r,by=c('from','to')) |> dplyr::mutate(nutrient=nu)
  cat(nu,': ',length(ids),' nodes; ',nrow(z),' ordered pairs; R=1 ',sum(z$R==1,na.rm=TRUE),'; R=? ',sum(is.na(z$R)),'\n',sep='')
  z
}
all<-dplyr::bind_rows(lapply(c('TN','TP'),read_nu))
readr::write_csv(all |> dplyr::count(nutrient,R,name='n'),file.path(OUT,'00_physical_status_QA.csv'))

# (a) Physical upstream -> downstream, original TE in both orientations.
# TE_raw can be zero because of estimation or prior filtering; diagnose before interpretation.
forward <- all |> dplyr::filter(R==1) |>
  dplyr::transmute(nutrient, up=from, down=to, TE_forward=TE, p_forward=p)
reverse <- all |>
  dplyr::transmute(nutrient, up=to, down=from, TE_reverse=TE, p_reverse=p)
a <- dplyr::inner_join(forward,reverse,by=c('nutrient','up','down')) |>
  dplyr::filter(up!=down) |>
  dplyr::mutate(delta=TE_forward-TE_reverse,
                comparison=dplyr::case_when(delta>1e-10~'Forward stronger',
                                            delta< -1e-10~'Reverse stronger',
                                            TRUE~'Tie'),
                zero_pattern=dplyr::case_when(TE_forward==0 & TE_reverse==0~'Both zero',
                                              TE_forward==0~'Forward zero',
                                              TE_reverse==0~'Reverse zero',
                                              TRUE~'Both positive'))
if(!nrow(a))stop('S2(a): no confirmed R=1 directional pairs with both TE orientations')
if(anyDuplicated(a[c('nutrient','up','down')]))stop('S2(a): duplicated physical direction pair')
a_summary <- a |> dplyr::group_by(nutrient) |>
  dplyr::summarise(n_pairs=dplyr::n(),median_delta=median(delta),
                   n_forward=sum(comparison=='Forward stronger'),
                   n_reverse=sum(comparison=='Reverse stronger'),
                   n_tie=sum(comparison=='Tie'),
                   pct_forward=100*n_forward/n_pairs,
                   pct_reverse=100*n_reverse/n_pairs,
                   pct_tie=100*n_tie/n_pairs,
                   zero_either=sum(TE_forward==0 | TE_reverse==0),
                   n_forward_zero_sig=sum(TE_forward==0 & p_forward<ALPHA),
                   n_reverse_zero_sig=sum(TE_reverse==0 & p_reverse<ALPHA),
                   .groups='drop')
readr::write_csv(a,file.path(OUT,'s2a_bidirectional_TE_pairs.csv'))
readr::write_csv(a_summary,file.path(OUT,'s2a_directionality_summary.csv'))
cat('S2(a) zero-value QA (inspect if zeros arise from filtering rather than raw estimation):\n')
print(a_summary)
labels_a <- a_summary |> dplyr::mutate(
  label=sprintf('n=%d | median delta=%.4f\nF %.1f%% / R %.1f%% / tie %.1f%%',
                n_pairs,median_delta,pct_forward,pct_reverse,pct_tie))
p_a <- ggplot(a,aes(TE_forward,TE_reverse))+
  geom_abline(slope=1,intercept=0,linetype='dashed',colour='#838B91')+
  geom_point(aes(colour=nutrient,shape=zero_pattern),alpha=.72,size=1.8)+
  facet_wrap(~nutrient,scales='free')+
  scale_colour_manual(values=COL)+
  scale_shape_manual(values=c('Both positive'=16,'Both zero'=4,
                              'Forward zero'=17,'Reverse zero'=15))+
  geom_label(data=labels_a,aes(x=Inf,y=Inf,label=label),inherit.aes=FALSE,
             hjust=1.02,vjust=1.1,size=2.65,label.size=0,fill='white')+
  labs(title='(a) Information directionality of connected pairs',
       subtitle='All original TE values; zero-value source must be verified',
       x='Upstream-to-downstream TE',y='Downstream-to-upstream TE',shape='TE zero pattern')+
  theme(legend.position='bottom')

# (b) Discrete lag frequency: exact integer 1-3 (or actual available lags), not boxplots.
b_all <- all |> dplyr::filter(p<ALPHA,TE>0,is.finite(lag))
b <- b_all |> dplyr::filter(!is.na(R)) |>
  dplyr::mutate(status=factor(ifelse(R==1,'Reachable','Non-reachable'),
                              levels=c('Reachable','Non-reachable')))
if(!nrow(b))stop('S2(b): no significant TE edges with confirmed physical status')
if(any(abs(b$lag-round(b$lag))>1e-8))stop('S2(b): TE lags are not integer months')
b$lag <- as.integer(round(b$lag))
lag_levels <- sort(unique(b$lag))
b_stats <- b |> dplyr::group_by(nutrient,status) |>
  dplyr::summarise(n=dplyr::n(),median_lag=median(lag),q25=quantile(lag,.25),
                   q75=quantile(lag,.75),.groups='drop')
b_unknown <- b_all |> dplyr::group_by(nutrient) |>
  dplyr::summarise(n_unknown=sum(is.na(R)),.groups='drop')
b_freq <- b |> dplyr::count(nutrient,status,lag,name='n_lag') |>
  tidyr::complete(nutrient,status,lag=lag_levels,fill=list(n_lag=0)) |>
  dplyr::left_join(b_stats |> dplyr::select(nutrient,status,n),by=c('nutrient','status')) |>
  dplyr::mutate(percentage=100*n_lag/n,lag=factor(lag,levels=lag_levels))
readr::write_csv(b_stats,file.path(OUT,'s2b_lag_group_summary.csv'))
readr::write_csv(b_unknown,file.path(OUT,'s2b_unknown_status.csv'))
readr::write_csv(b_freq,file.path(OUT,'s2b_lag_frequency.csv'))
labels_b <- b_stats |> dplyr::mutate(label=sprintf('n=%d; med. %.1f',n,median_lag))
p_b <- ggplot(b_freq,aes(x=lag,y=percentage,fill=status))+
  geom_col(width=.68,colour='white',linewidth=.25)+
  facet_grid(nutrient~status)+
  geom_text(aes(label=ifelse(n_lag>0,sprintf('%.0f%%',percentage),'')),
            vjust=-.4,size=2.6)+
  scale_fill_manual(values=c(Reachable='#69A5A4','Non-reachable'='#B1BCC5'))+
  scale_y_continuous(limits=c(0,110),breaks=c(0,25,50,75,100),labels=function(x)paste0(x,'%'))+
  labs(title='(b) Distribution of significant TE lags',
       subtitle='Raw p < 0.10; known physical statuses only',
       x='TE lag (months)',y='Fraction of edges (%)')+
  theme(legend.position='none')
# (c) Physical travel time: preferred verified path-list, one observation per directed pair.
read_travel<-function(nu){
  path<-file.path(PR,nu,'PNet_directed_reachable_pairs.csv')
  d<-readr::read_csv(check_file(path),col_types=readr::cols(.default=readr::col_character()),show_col_types=FALSE)
  fc<-colpick(d,c('from','source','upstream'));tc<-colpick(d,c('to','target','downstream'))
  tcol<-colpick(d,c('travel_days_medium','travel_days','travel_time_days','travel_days_fast'))
  z<-tibble::tibble(nutrient=nu,from=id8(d[[fc]]),to=id8(d[[tc]]),travel_days=as.numeric(d[[tcol]]),source_col=tcol) |>
    dplyr::filter(is.finite(travel_days),travel_days>=0,!is.na(from),!is.na(to),from!=to)
  if(anyDuplicated(z[c('from','to')]))stop(nu,': duplicate physical travel time directed pairs')
  z
}
travel<-dplyr::bind_rows(lapply(c('TN','TP'),read_travel))
# Restrict to R=1 in same node universe as full-period TE; report records lacking this confirmation.
known<-all |> dplyr::select(nutrient,from,to,R)
travel<-travel |> dplyr::left_join(known,by=c('nutrient','from','to'))
travel_qa<-travel |> dplyr::group_by(nutrient,source_col) |> dplyr::summarise(n_available=dplyr::n(),n_confirmed=sum(R==1,na.rm=TRUE),n_unresolved=sum(is.na(R)),.groups='drop')
readr::write_csv(travel_qa,file.path(OUT,'s2c_physical_travel_QA.csv'))
travel<-travel |> dplyr::filter(R==1)
if(!nrow(travel))stop('No confirmed R=1 physical travel times in full-period TE node universe')
travel<-travel |> dplyr::mutate(time_class=cut(travel_days,breaks=c(-Inf,1,7,30,Inf),labels=c('<=1 d','1-7 d','7-30 d','>30 d'),right=TRUE))
travel_stats<-travel |> dplyr::group_by(nutrient) |> dplyr::summarise(n=dplyr::n(),median_days=median(travel_days),pct_under30=100*mean(travel_days<30),.groups='drop')
travel_class<-travel |> dplyr::count(nutrient,time_class,name='n')
readr::write_csv(travel_stats,file.path(OUT,'s2c_travel_summary.csv'))
readr::write_csv(travel_class,file.path(OUT,'s2c_travel_categories.csv'))
p_c<-ggplot(travel,aes(travel_days,fill=nutrient))+
  geom_histogram(bins=22,colour='white',linewidth=.3)+
  geom_vline(data=travel_stats,aes(xintercept=median_days),linetype='dashed',colour='#3E4C59')+
  facet_wrap(~nutrient,scales='free_y')+
  geom_label(data=travel_stats,aes(x=Inf,y=Inf,label=sprintf('n=%d | median %.1f d | <30 d: %.1f%%',n,median_days,pct_under30)),inherit.aes=FALSE,hjust=1.04,vjust=1.25,size=2.4,label.size=0)+
  scale_fill_manual(values=COL)+
  labs(title='(c) Hydraulic travel time of verified downstream paths',x='Hydraulic travel time (days)',y='Number of directed station pairs')+
  theme(legend.position='none')

# (d) Choose an additional example, explicitly excluding pairs already in Fig. 3.
obs<-readr::read_csv(check_file(DATA_FILE),col_types=readr::cols(STAID=readr::col_character(),PCODE=readr::col_character(),.default=readr::col_guess()),show_col_types=FALSE)
req<-c('STAID','PCODE','Year','Month','KConc')
if(length(setdiff(req,names(obs))))stop('NTN file missing ',paste(setdiff(req,names(obs)),collapse=', '))
obs<-obs |> dplyr::transmute(station=id8(STAID),nutrient=dplyr::case_when(PCODE=='P00600'~'TN',PCODE=='P00665'~'TP',TRUE~NA_character_),
                             date=as.Date(sprintf('%04d-%02d-01',as.integer(Year),as.integer(Month))),month=as.integer(Month),value=as.numeric(KConc)) |>
  dplyr::filter(!is.na(nutrient),date>=as.Date('2011-01-01'),date<=as.Date('2023-12-01'),is.finite(value))
if(anyDuplicated(obs[c('station','nutrient','date')]))stop('Duplicate NTN monthly concentration station/nutrient/date rows. Resolve before S2(d).')
# Month-of-year standardized anomalies; descriptive statistics from 2011-2023.
climatology<-obs |> dplyr::group_by(nutrient,station,month) |>
  dplyr::summarise(mu=mean(value),sd=sd(value),n=dplyr::n(),.groups='drop')
obs<-obs |> dplyr::left_join(climatology,by=c('nutrient','station','month')) |>
  dplyr::mutate(Z=ifelse(is.finite(sd)&sd>0&n>=3,(value-mu)/sd,NA_real_))
# Eligible pairs: confirmed upstream -> downstream, not the designated Fig.3 examples.
physical_pairs<-a |> dplyr::filter(up!=down) |> dplyr::transmute(nutrient,up,down,key=paste(up,down,sep='>')) |>
  dplyr::filter(!key %in% EXCLUDED_FIG3)
# Deterministic selection: highest total shared months with standardized observations,
# then largest upstream standardized event. No selection based on downstream response.
upobs<-obs |> dplyr::select(nutrient,up=station,date,Z_up=Z)
dnobs<-obs |> dplyr::select(nutrient,down=station,date,Z_down=Z)
pair_month<-physical_pairs |> dplyr::inner_join(upobs,by=c('nutrient','up')) |>
  dplyr::inner_join(dnobs,by=c('nutrient','down','date')) |>
  dplyr::filter(is.finite(Z_up),is.finite(Z_down))
pair_rank<-pair_month |> dplyr::group_by(nutrient,up,down) |>
  dplyr::summarise(n_shared=dplyr::n(),max_up=max(Z_up),.groups='drop') |>
  dplyr::filter(n_shared>=36) |> dplyr::arrange(dplyr::desc(n_shared),dplyr::desc(max_up),up,down)
readr::write_csv(pair_rank,file.path(OUT,'s2d_available_pairs_ranked.csv'))
if(!nrow(pair_rank)) {
  warning('S2(d): no additional physical pair has >=36 shared standardized months. Panel will show data-availability explanation.')
  p_d<-ggplot()+annotate('text',x=0,y=0,label='No eligible additional physical pair with >=36 common months')+
    theme_void()+labs(title='(d) Observed monthly anomaly example')
} else {
  sel<-pair_rank |> dplyr::group_by(nutrient) |> dplyr::slice_head(n=1) |> dplyr::ungroup()
  selected<-pair_month |> dplyr::inner_join(sel |> dplyr::select(nutrient,up,down),by=c('nutrient','up','down')) |>
    dplyr::filter(date>=as.Date('2012-01-01'),date<=as.Date('2023-06-01')) |>
    dplyr::group_by(nutrient,up,down) |> dplyr::slice_max(Z_up,n=1,with_ties=FALSE) |> dplyr::ungroup()
  if(any(selected$up==selected$down))stop('S2(d): self-pair selected; check IDs')
  readr::write_csv(selected,file.path(OUT,'s2d_selected_events.csv'))
  series<-dplyr::bind_rows(lapply(seq_len(nrow(selected)),function(i){
    s<-selected[i,]
    if(is.na(s$up) || is.na(s$down) || s$up==s$down)stop('S2(d): source and downstream must be different station IDs')
    low<-seq(s$date,by='-1 month',length.out=WINDOW+1)[WINDOW+1]
    high<-seq(s$date,by='1 month',length.out=WINDOW+1)[WINDOW+1]
    obs |> dplyr::filter(nutrient==s$nutrient,station %in% c(s$up,s$down),date>=low,date<=high) |>
      dplyr::mutate(role=ifelse(station==s$up,'Upstream','Downstream'),
                    display=paste0(nutrient,' | ',s$up,' -> ',s$down,' | source event: ',format(s$date,'%Y-%m')))
  }))
  readr::write_csv(series,file.path(OUT,'s2d_observed_standardized_series.csv'))
  if(anyNA(series$role) || !all(unique(series$role) %in% c('Upstream','Downstream')))stop('S2(d): unexpected station role')
  p_d<-ggplot(series,aes(date,Z,colour=role,linetype=role,group=role))+
    geom_hline(yintercept=0,colour='#AAB2B8',linewidth=.3)+
    geom_hline(yintercept=1.5,linetype='dotted',colour='#B36647')+
    geom_line(linewidth=.75,na.rm=TRUE)+geom_point(size=1.25,na.rm=TRUE)+
    facet_wrap(~display,ncol=1,scales='free_x')+
    scale_colour_manual(values=c(Upstream='#327FA6',Downstream='#E18652'))+
    scale_linetype_manual(values=c(Upstream='solid',Downstream='dashed'))+
    scale_x_date(date_breaks='3 months',date_labels='%Y-%m')+
    labs(title='(d) Observed TN/TP anomalies in additional connected pairs',subtitle='Illustrative; selected by shared coverage then largest upstream anomaly',x=NULL,y='Standardized concentration anomaly (Z)',colour=NULL,linetype=NULL)+
    theme(axis.text.x=element_text(angle=35,hjust=1),legend.position='bottom')
}

base<-theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),plot.title=element_text(face='bold',size=10),strip.background=element_rect(fill='#F1F3F4',colour='#D0D5DA'))
panels<-list(p_a=p_a+base,p_b=p_b+base,p_c=p_c+base,p_d=p_d+base)
figure<-(panels$p_a+panels$p_b)/(panels$p_c+panels$p_d)+
  patchwork::plot_annotation(title='Figure S2 | Upstream-Downstream Information and Dual-Lag Characteristics',subtitle='Full period: 2011-2023 | TN and TP | physical reachability confirmed only',theme=theme(plot.title=element_text(face='bold',size=14)))
for(nm in names(panels))ggsave(file.path(OUT,paste0('Figure_S2_',nm,'.pdf')),panels[[nm]],width=8.5,height=5.0,device=cairo_pdf)
ggsave(file.path(OUT,'Figure_S2_fullperiod.pdf'),figure,width=15.8,height=11.8,device=cairo_pdf)
ggsave(file.path(OUT,'Figure_S2_fullperiod.png'),figure,width=15.8,height=11.8,dpi=400,bg='white')
writeLines(c('S2(a) uses TE_raw for both directions. Zero-value patterns are audited; TE_raw must be verified to be unfiltered before making claims of original information directionality.','S2(b) shows discrete monthly TE lag proportions for p<0.10 positive TE; unknown R excluded. See s2b_lag_frequency.csv.','S2(c) uses hydraulic travel_days_medium from verified directed paths; 30-day fraction is approximate month-scale comparison.','S2(d) standardization uses the 2011-2023 calendar-month climatology and is descriptive, not out-of-sample.','S2(d) pair selection based on data coverage and upstream event, never downstream response.'),file.path(OUT,'README_methods.txt'))
cat('\nDONE: ',OUT,'\n')
