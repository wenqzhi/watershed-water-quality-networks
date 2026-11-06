source("config.R")
# ================================================================
# Figure S3 | Sensitivity and Stability of Network Influence
# Revised panel (d): marginal share by propagation step
# 2011-2023 TN/TP | explicit provenance; NO FDR / NO TE recomputation
# IMPORTANT: PCS_FUNCTION below is a PROXY until Eq.(11)/Eq.(S10)
# and Fig.4 code are checked. Do not call results manuscript PCS yet.
# ================================================================
need <- c('readr','dplyr','ggplot2','patchwork','tidyr')
missing <- need[!vapply(need, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) stop('Install packages first: ',paste(missing,collapse=', '))

ROOT <- get_data_root()
FROOT <- file.path(ROOT,'Fusion_Net_2011_2023')
OUT <- file.path(ROOT,'Figure_S3_PCS_sensitivity')
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

# USE ONLY AFTER confirmed against manuscript PCS definition:
PCS_DEFINITION <- 'discounted_walk_proxy' # explicit method name, NOT final manuscript PCS
OMEGA0 <- 0.5
L0 <- 6L  # replace with actual manuscript baseline if different
LS <- c(2L,3L,4L,5L,6L,8L,10L)
OMEGAS <- c(0,0.2,0.4,0.5,0.6,0.8,1)
DISCOUNT <- 0.85 # PROXY ONLY; replace using the manuscript if different
BACKBONE_Q <- 0.85 # from existing network display script
USE_BACKBONE <- TRUE # TRUE: positive edge weight >= qth quantile; zero edges excluded
COL <- c(TN='#327FA6',TP='#E18652')

stopifnot(L0 %in% LS,OMEGA0 %in% OMEGAS,DISCOUNT>=0,DISCOUNT<=1)

read_layers <- function(nu) {
  path <- file.path(FROOT,nu,'Fusion_all_pairs_layers.csv')
  if(!file.exists(path)) stop(nu,': cannot find ',path)
  d <- readr::read_csv(path,show_col_types=FALSE,col_types=readr::cols(.default=readr::col_guess(),from=readr::col_character(),to=readr::col_character()))
  req <- c('from','to','T_norm','P_norm')
  if(!all(req %in% names(d))) stop(nu,': missing columns: ',paste(setdiff(req,names(d)),collapse=', '))
  d <- dplyr::transmute(d,from=trimws(as.character(.data$from)),to=trimws(as.character(.data$to)),T_norm=as.numeric(.data$T_norm),P_norm=as.numeric(.data$P_norm))
  if(anyNA(d) || any(!is.finite(d$T_norm)) || any(!is.finite(d$P_norm))) stop(nu,': invalid layer values')
  if(any(d$T_norm<0|d$P_norm<0)) stop(nu,': negative layer weights not allowed')
  if(any(d$from==d$to)) stop(nu,': contains self-links')
  if(anyDuplicated(paste(d$from,d$to,sep='->'))) stop(nu,': duplicate directed pairs')
  ids <- sort(unique(c(d$from,d$to)))
  if(nrow(d)!=length(ids)*(length(ids)-1L)) stop(nu,': require complete candidate-pair edge layer: ',nrow(d),' vs ',length(ids)*(length(ids)-1))
  list(edges=d,ids=ids)
}

make_W <- function(obj,omega) {
  d <- obj$edges; ids <- obj$ids
  w <- omega*d$T_norm + (1-omega)*d$P_norm
  if(USE_BACKBONE) {
    pos <- w[is.finite(w)&w>0]
    if(!length(pos)) stop('Fusion network has zero positive weights')
    cut <- as.numeric(stats::quantile(pos,BACKBONE_Q,names=FALSE,type=7))
    w[w<cut] <- 0
  }
  W <- matrix(0,length(ids),length(ids),dimnames=list(ids,ids))
  W[cbind(match(d$from,ids),match(d$to,ids))] <- w
  diag(W) <- 0
  W
}

# ---------------------------------------------------------------
# REPLACE ONLY THIS FUNCTION after verifying the REAL PCS formula.
# PROXY: sum of discounted powers of globally scaled adjacency.
# A = W / max(rowSums(W)); S = A + discount*A^2 + ...
# Global output influence is the ROW sum of S.
# Global scaling retains row-strength heterogeneity, unlike row-normalization.
# If the manuscript uses other normalization, path restrictions,
# gates, attenuation or directional convention, adapt this function.
# ---------------------------------------------------------------
PCS_FUNCTION <- function(W,L,discount=DISCOUNT) {
  if(any(!is.finite(W)) || any(W<0)) stop('W must be finite and nonnegative')
  rs <- rowSums(W)
  if(max(rs)<=0) stop('No positive network edges')
  P <- W / max(rs)
  S <- matrix(0,nrow(W),ncol(W),dimnames=dimnames(W))
  power <- diag(nrow(W))
  for(h in seq_len(L)) {
    power <- power %*% P
    S <- S + discount^(h-1L)*power
  }
  setNames(rowSums(S),rownames(W))
}

# S3(d): step-level components of EXACTLY THE SAME propagation operator.
# Denominator is the combined contribution of steps 1..MARGINAL_MAX_L.
# This is mathematical decomposition, not empirical nutrient travel time.
MARGINAL_MAX_L <- 10L
step_contributions <- function(W, max_L=MARGINAL_MAX_L, discount=DISCOUNT) {
  if(any(!is.finite(W)) || any(W<0)) stop('W must be finite and nonnegative')
  rs <- rowSums(W)
  if(max(rs)<=0) stop('No positive network edges')
  P <- W / max(rs)
  power <- diag(nrow(W))
  raw <- numeric(max_L)
  for(h in seq_len(max_L)) {
    power <- power %*% P
    raw[h] <- sum(power) * discount^(h-1L)
  }
  if(any(!is.finite(raw)) || any(raw<0) || sum(raw)<=0)
    stop('Invalid step-wise contribution')
  data.frame(step=seq_len(max_L), raw_contribution=raw,
             share_pct=100*raw/sum(raw),
             cumulative_pct=100*cumsum(raw)/sum(raw))
}

# If the proxy ranking is almost flat, interpret no concentration/importance.
# Most importantly, this proxy is NOT assumed to equal the manuscript PCS.
get_pcs <- function(obj,omega,L) {
  W <- make_W(obj,omega)
  pcs <- PCS_FUNCTION(W,L)
  if(any(!is.finite(pcs)) || any(pcs<0)) stop('Invalid PCS values')
  if(sum(pcs)<=0) stop('All PCS are zero')
  pcs
}

top_ids <- function(x) {
  n <- length(x); k <- ceiling(.2*n)
  names(x)[order(-as.numeric(x),names(x))][seq_len(k)]
}
contrib20 <- function(x) 100*sum(x[top_ids(x)])/sum(x)
jacc <- function(a,b) length(intersect(a,b))/length(union(a,b))

# ---------------------------------------------------------------
# Parameter analysis. Fix omega=OMEGA0 for all L comparisons.
# S3(c): fixed baseline top-20% stations; the contribution share is
# PCS_i(L)/sum_j PCS_j(L), not unnormalized absolute PCS.
# Delta is in PERCENTAGE POINTS relative to baseline L0.
# S3(d) separately decomposes the sum of steps 1 to MARGINAL_MAX_L.
# ---------------------------------------------------------------
res <- lapply(c('TN','TP'), function(nu) {
  obj <- read_layers(nu)
  cat('\n',nu,': ',length(obj$ids),'nodes / ',nrow(obj$edges),' pairs\n',sep='')
  base <- get_pcs(obj,OMEGA0,L0)
  base_ids <- top_ids(base)
  base_share <- base / sum(base)
  
  byL <- dplyr::bind_rows(lapply(LS,function(L) {
    pcs <- get_pcs(obj,OMEGA0,L)
    data.frame(nutrient=nu,L=L,node=names(pcs),PCS=unname(pcs),
               share_pct=100*unname(pcs/sum(pcs)),
               delta_share_pp=100*unname(pcs/sum(pcs)-base_share[names(pcs)]),
               baseline_top20=names(pcs)%in%base_ids,
               current_top20=names(pcs)%in%top_ids(pcs),
               stringsAsFactors=FALSE)
  }))
  byOmega <- dplyr::bind_rows(lapply(OMEGAS,function(om) {
    pcs <- get_pcs(obj,om,L0)
    data.frame(nutrient=nu,omega=om,node=names(pcs),PCS=unname(pcs),
               top20=names(pcs)%in%top_ids(pcs),stringsAsFactors=FALSE)
  }))
  
  A <- byL |>
    dplyr::group_by(.data$nutrient,.data$L) |>
    dplyr::summarise(contribution20=100*sum(.data$PCS[.data$current_top20])/sum(.data$PCS),.groups='drop')
  B <- dplyr::bind_rows(lapply(OMEGAS,function(om) {
    pcs <- get_pcs(obj,om,L0)
    data.frame(nutrient=nu,omega=om,
               spearman=if(stats::sd(pcs)==0||stats::sd(base)==0) NA_real_ else
                 stats::cor(pcs[names(base)],base,method='spearman'),
               jaccard20=jacc(top_ids(pcs),base_ids))
  }))
  C_nodes <- byL |>
    dplyr::filter(.data$baseline_top20) |>
    dplyr::mutate(baseline_L=L0)
  C <- C_nodes |>
    dplyr::group_by(.data$nutrient,.data$L) |>
    dplyr::summarise(median_delta=stats::median(.data$delta_share_pp),
                     q25=as.numeric(stats::quantile(.data$delta_share_pp,0.25)),
                     q75=as.numeric(stats::quantile(.data$delta_share_pp,0.75)),
                     min_delta=min(.data$delta_share_pp),
                     max_delta=max(.data$delta_share_pp),
                     n_nodes=dplyr::n(),.groups='drop')
  # Actual step contributions at baseline omega, on a FIXED W.
  # Also independently reconstruct PCS at max step to check exact equivalence.
  W0 <- make_W(obj,OMEGA0)
  D <- step_contributions(W0) |>
    dplyr::mutate(nutrient=nu, .before=1)
  pcs_max <- PCS_FUNCTION(W0,MARGINAL_MAX_L)
  # This QA ensures the step-wise denominator equals the sum of all node PCS.
  if(!isTRUE(all.equal(sum(D$raw_contribution),sum(pcs_max),
                       tolerance=1e-9)))
    stop(nu,': step decomposition does not match PCS_FUNCTION')
  # Original hub Jaccard across L stays in QA table, but no longer plotted.
  D_hubs <- dplyr::bind_rows(lapply(LS,function(L) {
    pcs <- get_pcs(obj,OMEGA0,L)
    data.frame(nutrient=nu,L=L,
               spearman=if(stats::sd(pcs)==0||stats::sd(base)==0) NA_real_ else
                 stats::cor(pcs[names(base)],base,method='spearman'),
               jaccard20=jacc(top_ids(pcs),base_ids))
  }))
  QA <- data.frame(nutrient=nu,n_nodes=length(base),top20_n=length(base_ids),
                   positive_edges=sum(make_W(obj,OMEGA0)>0),
                   pcs_min=min(base),pcs_max=max(base),pcs_sd=stats::sd(base),
                   concentration20=contrib20(base),method=PCS_DEFINITION,
                   omega0=OMEGA0,L0=L0,discount=DISCOUNT,
                   backbone_q=if(USE_BACKBONE) BACKBONE_Q else NA_real_)
  list(a=A,b=B,c=C,c_nodes=C_nodes,d=D,d_hubs=D_hubs,byL=byL,byOmega=byOmega,qa=QA,
       baseline_hubs=data.frame(nutrient=nu,node=base_ids,baseline_PCS=unname(base[base_ids]),
                                baseline_share_pct=100*unname(base_share[base_ids])))
})

bind <- function(x) dplyr::bind_rows(lapply(res,`[[`,x))
a<-bind('a');b<-bind('b');cc<-bind('c');dd<-bind('d');qa<-bind('qa')
for(nm in c('a','b','c','c_nodes','d','d_hubs','byL','byOmega','qa','baseline_hubs'))
  readr::write_csv(bind(nm),file.path(OUT,paste0('S3_',nm,'.csv')))

if(any(qa$pcs_sd < 1e-10)) warning('Near-constant PCS: proxy rankings degenerate.')
if(anyNA(b$spearman)) warning('NA Spearman: inspect constant node scores.')

th <- ggplot2::theme_bw(base_size=11)+ggplot2::theme(
  plot.title=ggplot2::element_text(face='bold',size=11),
  legend.position='bottom',panel.grid.minor=ggplot2::element_blank(),
  plot.margin=ggplot2::margin(7,7,7,7))

pa <- ggplot2::ggplot(a,ggplot2::aes(x=L,y=contribution20,color=nutrient))+
  ggplot2::geom_line(linewidth=.8)+ggplot2::geom_point(size=2)+
  ggplot2::scale_color_manual(values=COL)+ggplot2::scale_x_continuous(breaks=LS)+
  ggplot2::labs(title='(a) Top-20% influence concentration',x='Maximum network steps',
                y='Top-20% contribution (%)',color=NULL)+th

pb <- ggplot2::ggplot(b,ggplot2::aes(x=omega,y=spearman,color=nutrient))+
  ggplot2::geom_line(linewidth=.8)+ggplot2::geom_point(size=2)+
  ggplot2::geom_vline(xintercept=OMEGA0,linetype=2,color='grey60')+
  ggplot2::scale_color_manual(values=COL)+ggplot2::scale_x_continuous(breaks=OMEGAS)+
  ggplot2::scale_y_continuous(limits=c(0,1))+
  ggplot2::labs(title='(b) Rank stability across fusion weights',x='TE coefficient (omega)',
                y='Spearman rank correlation',color=NULL)+th

# S3(c): every thin line is the SAME baseline hub over L. The bold line is
# the median across hubs. Shaded region is the interquartile RANGE, not a CI.
cn <- bind('c_nodes')
pc <- ggplot2::ggplot() +
  ggplot2::geom_hline(yintercept=0,linetype=2,color='grey60')+
  ggplot2::geom_line(data=cn,ggplot2::aes(x=L,y=delta_share_pp,group=node),
                     color='grey55',alpha=.22,linewidth=.35)+
  ggplot2::geom_ribbon(data=cc,ggplot2::aes(x=L,ymin=q25,ymax=q75),
                       fill='#B6BFC4',alpha=.3)+
  ggplot2::geom_line(data=cc,ggplot2::aes(x=L,y=median_delta,color=nutrient),linewidth=.85)+
  ggplot2::geom_point(data=cc,ggplot2::aes(x=L,y=median_delta,color=nutrient),size=1.8)+
  ggplot2::geom_vline(xintercept=L0,linetype=3,color='grey65')+
  ggplot2::facet_wrap(~nutrient,nrow=1,scales='free_y')+
  ggplot2::scale_color_manual(values=COL,guide='none')+
  ggplot2::scale_x_continuous(breaks=LS)+
  ggplot2::labs(title='(c) Relative contribution change of baseline hubs',
                subtitle='Fixed baseline top 20%; median and IQR across hubs (not CI)',
                x='Maximum network steps',y='Change in PCS share (percentage points)')+th

pd <- ggplot2::ggplot(dd,ggplot2::aes(x=step,y=share_pct,color=nutrient))+
  ggplot2::geom_line(linewidth=.85)+ggplot2::geom_point(size=2.1)+
  ggplot2::geom_vline(xintercept=L0,linetype=2,color='grey60')+
  ggplot2::scale_color_manual(values=COL)+
  ggplot2::scale_x_continuous(breaks=seq_len(MARGINAL_MAX_L))+
  ggplot2::scale_y_continuous(labels=function(x) paste0(x,'%'),expand=ggplot2::expansion(mult=c(.02,.08)))+
  ggplot2::labs(title='(d) Marginal contribution of each propagation step',
                subtitle='Step share of cumulative influence (steps 1–10)',
                x='Propagation step',y='Share of total influence (%)',color=NULL)+th

p <- ((pa+pb)/(pc+pd)) + patchwork::plot_annotation(
  title='Figure S3 | Sensitivity and Stability of Network Influence',
  subtitle=paste0('2011–2023 | PROVISIONAL ',PCS_DEFINITION,
                  '; omega0=',OMEGA0,', L0=',L0,'; verify manuscript PCS before publication'),
  theme=ggplot2::theme(plot.title=ggplot2::element_text(face='bold',size=15)))

ggplot2::ggsave(file.path(OUT,'Figure_S3_PROVISIONAL.pdf'),p,width=13,height=9.3)
ggplot2::ggsave(file.path(OUT,'Figure_S3_PROVISIONAL.png'),p,width=13,height=9.3,dpi=400)
for(z in c('a','b','c','d'))
  ggplot2::ggsave(file.path(OUT,paste0('Figure_S3_',z,'.pdf')),get(paste0('p',z)),width=6.5,height=4.7)
cat('\nOutputs in: ',OUT,'\n',sep='')
print(qa)
cat('\nBaseline hub overlaps across L (Jaccard, QA only):\n')
print(bind('d_hubs')[,c('nutrient','L','jaccard20')])
cat('\nStep decomposition (sum of shares should be 100% for each nutrient):\n')
print(dd)
cat('\nIMPORTANT: PCS_FUNCTION is a proxy; match Eq.(11), Eq.(S10), and Fig.4 before publication.\n')
