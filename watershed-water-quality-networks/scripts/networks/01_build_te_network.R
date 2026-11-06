source("config.R")
# NTN Chesapeake Bay: directed transfer entropy networks
# MAIN 2011-10 -- 2023-09; historical non-overlapping controls; spatial extension
# Requires: readr, dplyr, tidyr, lubridate, ggplot2, sf, igraph
# NOTE: TE = predictive information dependence, NOT proof of transport.

pkgs <- c('readr','dplyr','tidyr','lubridate','ggplot2','sf','igraph')
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
invisible(lapply(pkgs, library, character.only=TRUE))

# =========== CONFIGURATION: EDIT THESE PATHS =============
ROOT <- get_data_root()
CSV <- file.path(ROOT, 'NTN_2023_MonthlyLoadTable.csv')
STATIONS <- file.path(ROOT, 'NTN-station.txt')
WATERSHED <- file.path(ROOT, 'ntnwatersheds/ntnwatersheds.shp')
RIVERS <- file.path(ROOT, '水系/水系.shp')
OUT <- file.path(ROOT, 'TE_network_results')
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

PCODE <- 'P00665'      # total nitrogen. Change to P00665 for TP
VALUE <- 'KConc'       # use KConc, FNConc, KLoad, FNLoad separately
MAX_LAG <- 3L         # sensitivity: 1, 3, 6 months
B <- 199L             # increase to 499/999 for final inference
FDR <- 0.05
SEED <- 20261023
WINDOW_MODE <- 'main' # 'main', 'early', 'older', 'spatial'
# For entire 107-node network use all ordered station pairs, ~11k pairs.
# Computationally expensive. Start by setting MAX_NODES=12; set Inf for full run.
MAX_NODES <- Inf
# Estimate lag on training series only, with a circular-shift surrogate max over lags.
# Use identical effective time support for all lags.
# B surrogates => smallest attainable p is 1/(B+1)

WINDOWS <- list(
  main    = c('2011-10-01','2023-09-01'),
  early   = c('2004-10-01','2011-09-01'),
  older   = c('1984-10-01','2004-09-01'),
  spatial = c('2014-10-01','2023-09-01')
)

# =========== DATA =============
d <- readr::read_csv(CSV, col_types=readr::cols(STAID='c', PCODE='c', .default=readr::col_guess()), show_col_types=FALSE) %>%
  dplyr::filter(PCODE==!!PCODE) %>%
  dplyr::mutate(STAID=sprintf('%08s',STAID), Date=as.Date(sprintf('%04d-%02d-01',Year,Month)))
d$STAID <- gsub(' ','0',d$STAID,fixed=TRUE)
stopifnot(VALUE %in% names(d))
if (anyDuplicated(d[c('STAID','Date')])) stop('Duplicate station/month records for selected PCODE')

# Screening: exactly one valid selected parameter for every month in selected period
prepare <- function(period, eligible=NULL) {
  dates <- seq(as.Date(period[1]), as.Date(period[2]), by='month')
  dd <- d %>% filter(Date %in% dates, is.finite(.data[[VALUE]]))
  ids <- dd %>% count(STAID, name='n') %>% filter(n==length(dates)) %>% pull(STAID)
  if(!is.null(eligible)) ids <- intersect(ids,eligible)
  ids <- sort(ids)
  if (is.finite(MAX_NODES)) ids <- head(ids,MAX_NODES)
  wide <- dd %>% filter(STAID %in% ids) %>% select(Date,STAID,all_of(VALUE)) %>%
    tidyr::pivot_wider(names_from=STAID,values_from=all_of(VALUE)) %>%
    right_join(tibble(Date=dates),by='Date') %>% arrange(Date)
  if (ncol(wide)<3L || anyNA(wide)) stop('Insufficient complete stations or unexpected missing values')
  list(x=as.matrix(wide[,-1,drop=FALSE]),dates=dates,ids=ids)
}

# =========== PREPROCESS =============
# log1p then remove within-window mean month-of-year and simple linear trend.
# Treat 0 as valid. Negative concentrations cause halt.
preprocess <- function(a,dates) {
  if(any(a<0,na.rm=TRUE)) stop('Negative concentrations found')
  z <- log1p(a)
  mo <- lubridate::month(dates)
  tt <- seq_len(nrow(z))
  for(j in seq_len(ncol(z))) {
    x <- z[,j]
    x <- x - ave(x,mo,FUN=mean)
    x <- residuals(lm(x~tt))
    z[,j] <- x
  }
  z
}

# Discretize each series into terciles, fixed across all lags and surrogates.
discretize <- function(z) {
  apply(z,2,function(x) {
    br <- quantile(x,c(1/3,2/3),names=FALSE)
    if(br[1]>=br[2]) {
      # Falls back to rank terciles in presence of ties
      as.integer(cut(rank(x,ties.method='average'), breaks=c(-Inf,length(x)/3,2*length(x)/3,Inf),labels=FALSE))
    } else as.integer(cut(x,breaks=c(-Inf,br[1],br[2],Inf),labels=FALSE))
  })
}

# Conditional MI I(X;Y|Z) in nats, plug-in estimate.
# X=source past, Y=target current, Z=target own previous month.
cmi <- function(x,y,z,K=3L) {
  n <- length(x)
  if(n==0L || anyNA(x) || anyNA(y) || anyNA(z)) return(NA_real_)
  arr <- array(tabulate(x + K*(y-1L) + K*K*(z-1L), nbins=K^3),c(K,K,K))
  p <- arr / n
  pxz <- apply(p,c(1,3),sum)
  pyz <- apply(p,c(2,3),sum)
  pz <- apply(p,3,sum)
  val <- 0
  for(k in seq_len(K)) for(j in seq_len(K)) for(i in seq_len(K)) {
    if(p[i,j,k]>0 && pxz[i,k]>0 && pyz[j,k]>0 && pz[k]>0) {
      val <- val + p[i,j,k]*log(p[i,j,k]*pz[k]/(pxz[i,k]*pyz[j,k]))
    }
  }
  max(0,val)
}

# Use common target indices, so no inflation from different sample sizes across lags.
lag_te <- function(source,target,maxlag=MAX_LAG) {
  tix <- seq.int(maxlag+1L,length(target))
  vapply(seq_len(maxlag),function(lag) {
    cmi(source[tix-lag],target[tix],target[tix-1L])
  },numeric(1))
}

# Circular shifts of source destroy particular alignment while preserving its
# distribution and serial pattern. Exclude near-zero shifts (±12 months).
# Seasonal structure has already been removed.
surrogate_shift <- function(n) {
  good <- seq.int(1L,n-1L)
  dist <- pmin(good,n-good)
  good <- good[dist>max(12L,MAX_LAG+1L)]
  if(!length(good)) stop('Series too short for circular-shift surrogate')
  sample(good,1)
}
shift_vec <- function(x,k) {
  n <- length(x)
  x[c((n-k+1L):n,seq_len(n-k))]
}

calc_network <- function(period,label,eligible=NULL) {
  raw <- prepare(period,eligible)
  mat <- discretize(preprocess(raw$x,raw$dates))
  n <- ncol(mat)
  if(nrow(mat)<=2*MAX_LAG+28L) stop('Too few months for TE test')
  cat(sprintf('\n%s: %s stations x %s months; %s directed pairs\n',label,n,nrow(mat),n*(n-1L)))
  pairs <- which(row(matrix(0,n,n))!=col(matrix(0,n,n)),arr.ind=TRUE)
  # pairs indices are from row-to-column; source=column row number, target=column column number
  set.seed(SEED)
  res <- vector('list',nrow(pairs))
  for(k in seq_len(nrow(pairs))) {
    i <- pairs[k,1]; j <- pairs[k,2]
    xs <- mat[,i]; yt <- mat[,j]
    tv <- lag_te(xs,yt)
    score <- max(tv); best <- which.max(tv)
    # Max-T correction across lag search within each pair
    null <- replicate(B, max(lag_te(shift_vec(xs,surrogate_shift(length(xs))),yt)))
    p <- (1+sum(null>=score))/(B+1)
    res[[k]] <- data.frame(from=colnames(mat)[i],to=colnames(mat)[j],
                           lag_months=best,TE_nats=score,
                           TE_bias_corrected=max(0,score-mean(null)),
                           p_value=p, stringsAsFactors=FALSE)
    if(k %% 500 == 0) cat(label,': ',k,'/',nrow(pairs),' edges tested\n',sep='')
  }
  ed <- bind_rows(res) %>% mutate(q_value=p.adjust(p_value,method='BH'),
                                  significant=q_value<FDR)
  # Important: if B=199 and large number of pairs, BH may return no edges.
  # Report full evidence; do NOT relax q after viewing results.
  write_csv(ed,file.path(OUT,paste0(label,'_all_edges.csv')))
  write_csv(filter(ed,significant),file.path(OUT,paste0(label,'_significant_edges.csv')))
  write_csv(tibble(Date=raw$dates,as.data.frame(raw$x,check.names=FALSE)),
            file.path(OUT,paste0(label,'_input_matrix.csv')))
  write_csv(tibble(STAID=raw$ids),file.path(OUT,paste0(label,'_station_ids.csv')))
  # Keep all nodes, including isolated nodes. Network: weighted by bias-corrected TE.
  g <- igraph::graph_from_data_frame(filter(ed,significant),directed=TRUE,
                                     vertices=data.frame(name=raw$ids))
  nodes <- tibble(STAID=raw$ids,
                  indegree=degree(g,mode='in'),outdegree=degree(g,mode='out'),
                  in_strength=strength(g,mode='in',weights=E(g)$TE_bias_corrected),
                  out_strength=strength(g,mode='out',weights=E(g)$TE_bias_corrected))
  write_csv(nodes,file.path(OUT,paste0(label,'_node_metrics.csv')))
  write_graph(g,file.path(OUT,paste0(label,'_network.graphml')),format='graphml')
  metrics <- tibble(window=label,month_count=nrow(mat),stations=n,
                    candidate_edges=nrow(ed),significant_edges=sum(ed$significant),
                    density=sum(ed$significant)/(n*(n-1)))
  write_csv(metrics,file.path(OUT,paste0(label,'_summary.csv')))
  list(edges=ed,nodes=nodes,metrics=metrics,ids=raw$ids)
}

# =========== ESTIMATE THE MAIN NETWORK =============
main <- calc_network(WINDOWS$main,'main_2011_2023')

# OPTIONAL: un-comment ONE of the following lines to run a separate analysis.
# NOTE historical networks have different nodes and shorter samples; NOT independent
# validation of the 2011-2023 edges if you use future data for tuning.
# early <- calc_network(WINDOWS$early,'early_2004_2011',eligible=main$ids)
# older <- calc_network(WINDOWS$older,'older_1984_2004',eligible=main$ids)
# spatial <- calc_network(WINDOWS$spatial,'spatial_2014_2023')

# =========== STATION MAP =============
# Station coordinates must use real lon/lat in NTN-station.txt.
if(file.exists(STATIONS)) {
  ss <- read_csv(STATIONS,col_types=cols(SITE_NO='c',.default=col_guess()),show_col_types=FALSE) %>%
    transmute(STAID=sprintf('%08s',SITE_NO),lon=as.numeric(lon),lat=as.numeric(lat))
  ss$STAID <- gsub(' ','0',ss$STAID,fixed=TRUE)
  ss <- left_join(ss,main$nodes,by='STAID') %>%
    mutate(group=if_else(STAID %in% main$ids,'Main TE nodes','Other NTN sites'))
  if(all(c('lon','lat')%in%names(ss))) {
    pts <- st_as_sf(filter(ss,is.finite(lon),is.finite(lat)),coords=c('lon','lat'),crs=4326,remove=FALSE)
    p <- ggplot()
    if(file.exists(WATERSHED)) {
      ws <- st_read(WATERSHED,quiet=TRUE)
      if(is.na(st_crs(ws))) stop('Watershed CRS missing; assign the verified original CRS first')
      p <- p + geom_sf(data=st_transform(ws,4326),fill='grey96',color='grey75',linewidth=.2)
    }
    if(file.exists(RIVERS)) {
      rv <- st_read(RIVERS,quiet=TRUE)
      if(is.na(st_crs(rv))) {
        bb <- st_bbox(rv)
        if(all(is.finite(bb)) && bb['xmin'] >= -180 && bb['xmax'] <= 180 &&
           bb['ymin'] >= -90 && bb['ymax'] <= 90) {
          warning('River shapefile missing CRS; coordinates look geographic. Assigning EPSG:4326; verify this.')
          st_crs(rv) <- 4326
        } else stop('River CRS missing and coordinates not geographic. Supply its real CRS.')
      }
      p <- p + geom_sf(data=st_transform(rv,4326),color='#8BBBD6',linewidth=.18,alpha=.7)
    }
    p <- p + geom_sf(data=pts,aes(color=group,size=out_strength),alpha=.9) +
      scale_color_manual(values=c('Main TE nodes'='#B44D43','Other NTN sites'='grey55')) +
      scale_size_continuous(range=c(1.3,4.2)) +
      theme_bw(base_size=11) + labs(color=NULL,size='TE out-strength',
                                    title='NTN monitoring stations and main-period information output')
    ggsave(file.path(OUT,'main_TE_station_map.pdf'),p,width=9,height=8)
    ggsave(file.path(OUT,'main_TE_station_map.png'),p,width=9,height=8,dpi=350)
  }
}

cat('\nDONE. Results folder: ',OUT,'\n',sep='')
