# Sequential benchmark worker. Run from the source root with a chosen library:
# Rscript tests/native/iid-benchmark.R LIB OUT prepare
# Rscript tests/native/iid-benchmark.R LIB OUT run CASE MODE CHAIN ITERATIONS [BLOCKING]
# MODE is separate or collapsed. Run fresh processes under /usr/bin/time -v;
# use one BLAS/OpenMP thread and the same CPU affinity for every worker.
args <- commandArgs(TRUE)
stopifnot(length(args) >= 3L)
.libPaths(c(args[1], .libPaths()))
suppressPackageStartupMessages(library(stLMM))
out <- args[2]; dir.create(out, recursive=TRUE, showWarnings=FALSE)
cases <- data.frame(
  case=c('intercept100','independent100','correlated100','correlated300',
    'correlated3_300','crossed300','nested100','by100','nngp_local100',
    'nngp_dispersed100','probit100','nb30'),
  groups=c(100L,100L,100L,300L,300L,300L,100L,100L,100L,100L,100L,30L),
  k=c(1L,2L,2L,2L,3L,2L,2L,2L,2L,2L,2L,2L),
  family=c(rep('gaussian',10),'probit','negative_binomial'),
  stringsAsFactors=FALSE)
make_data <- function(i) {
  spec <- cases[i,]; name <- spec$case; ng <- spec$groups; k <- spec$k
  spatial <- grepl('^nngp_',name)
  # The two spatial layouts share coordinates, covariates, process and noise.
  seed <- if(spatial) 88009L else 88000L+i
  n <- ng*10L
  set.seed(seed)
  d <- data.frame(group=factor(rep(seq_len(ng),each=10L)), x=rnorm(n), z=rnorm(n))
  if(spatial) {
    coords <- expand.grid(lon=seq(0,1,length.out=25),lat=seq(0,1,length.out=20))
    d$lon <- rep(coords$lon,each=2); d$lat <- rep(coords$lat,each=2)
    if(name=='nngp_dispersed100') {
      rng <- .Random.seed; set.seed(88999)
      d$group <- factor(rep(sample(rep(seq_len(ng),each=5)),each=2), levels=seq_len(ng))
      .Random.seed <<- rng
    }
  }
  S <- matrix(c(.8,.25,-.1,.25,.4,.08,-.1,.08,.3),3)[seq_len(k),seq_len(k),drop=FALSE]
  if(name=='independent100') S <- diag(diag(S))
  X <- cbind(1,d$x,d$z)[,seq_len(k),drop=FALSE]
  beta <- c(1,.4,-.25)[seq_len(k)]
  if(name=='probit100' || name=='nb30') beta[1] <- .1
  d$parent <- factor(rep(seq_len(if(name=='by100') 4L else 10L),each=n/(if(name=='by100') 4L else 10L)))
  matrices <- if(name=='by100') lapply(c(.5,1,1.5,2),function(a) S*a) else list(S)
  b <- matrix(rnorm(ng*k),ng,k) %*% chol(S)
  if(name=='by100') for(h in 1:4) {
    jj <- ((h-1)*ng/4+1):(h*ng/4)
    b[jj,] <- b[jj,]*sqrt(c(.5,1,1.5,2)[h])
  }
  effects <- b[as.integer(d$group),,drop=FALSE]
  formula <- switch(k,y~iid(group),y~x+iid(group,~1+x),y~x+z+iid(group,~1+x+z))
  pr <- list(iid_1=if(k==1L || name=='independent100') list(sigma_sq=ig(3,1)) else
    list(covariance=iw(k+2,diag(k))))
  covariance_truth <- list(iid_1=S)
  if(name=='independent100') formula <- y~x+iid(group,~1+x,cov_model='independent')
  if(name=='by100') {
    formula <- y~x+iid(group,~1+x,by=parent)
    covariance_truth$iid_1 <- setNames(matrices,levels(d$parent))
  }
  if(name=='crossed300') {
    d$crossed <- factor(sample(rep(seq_len(100L),length.out=n)))
    S2 <- matrix(c(.5,-.15,-.15,.3),2)
    b2 <- matrix(rnorm(200),100,2)%*%chol(S2)
    effects <- effects+b2[as.integer(d$crossed),]
    formula <- y~x+iid(group,~1+x)+iid(crossed,~1+x)
    pr$iid_2 <- list(covariance=iw(4,diag(2)));covariance_truth$iid_2 <- S2
  }
  if(name=='nested100') {
    effects[,1] <- effects[,1]+rnorm(10,sd=sqrt(.35))[as.integer(d$parent)]
    formula <- y~x+iid(parent)+iid(group,~1+x)
    pr <- list(iid_1=list(sigma_sq=ig(3,1)),iid_2=list(covariance=iw(4,diag(2))))
    covariance_truth <- list(iid_1=matrix(.35),iid_2=S)
  }
  fixed_iid <- drop(X%*%beta)+rowSums(X*effects)
  process <- rep(0,n)
  if(spatial) {
    formula <- update(formula,.~.+nngp(lon,lat,m=10,ordering='maxmin',cov_model='exp'))
    pr$nngp_1 <- list(sigma_sq=ig(3,.6),phi=uniform(1,20))
    d$y <- rep(0,n)
    setup <- suppressWarnings(stLMM(formula,d,priors=pr,n_samples=1,warmup=FALSE,verbose=FALSE))
    gr <- setup$backend$graphs[[1]]; term <- setup$backend$process_terms[[1]]
    m <- nrow(gr$coords_ord); C <- .3*exp(-6*as.matrix(dist(gr$coords_ord)))
    set.seed(seed+100L);w <- numeric(m)
    for(j in seq_len(m)) {
      count <- gr$nnIndxLU[m+j]
      nn <- if(count) gr$nnIndx[gr$nnIndxLU[j]+seq_len(count)]+1L else integer()
      weights <- if(count) solve(C[nn,nn,drop=FALSE],C[nn,j]) else numeric()
      variance <- C[j,j]-sum(weights*C[nn,j]);stopifnot(variance>0)
      w[j] <- sum(weights*w[nn])+sqrt(variance)*rnorm(1)
    }
    process <- w[term$map]
  }
  eta <- fixed_iid+process
  set.seed(seed+200L)
  d$y <- switch(spec$family,gaussian=eta+rnorm(n,sd=sqrt(.3)),
    probit=rbinom(n,1,pnorm(eta)),negative_binomial=rnbinom(n,size=2,mu=exp(eta)))
  list(spec=spec,data=d,formula=formula,priors=pr,beta=beta,
    covariance_truth=covariance_truth,fixed_iid=fixed_iid,predictor=eta,
    combined=sweep(effects,2,beta,'+'),seed=seed)
}
if(args[3]=='prepare') {
  for(i in seq_len(nrow(cases))) saveRDS(make_data(i),file.path(out,paste0(cases$case[i],'-data.rds')))
  write.csv(cases,file.path(out,'cases.csv'),row.names=FALSE)
  capture.output(sessionInfo(),file=file.path(out,'session-info.txt'))
  quit(status=0)
}
stopifnot(args[3]=='run',length(args)>=7L)
case <- args[4];mode <- args[5];chain <- as.integer(args[6]);ns <- as.integer(args[7])
blocking <- if(length(args)>=8L) args[8] else 'joint'
stopifnot(mode %in% c('separate','collapsed'),chain %in% 1:3,ns>=100L)
x <- readRDS(file.path(out,paste0(case,'-data.rds')))
burn <- ns%/%4L;keep <- seq.int(burn+1L,ns);spatial <- grepl('^nngp_',case)
starts <- list();multiplier <- c(.5,1,2)[chain]
for(nm in names(x$covariance_truth)) {
  S <- x$covariance_truth[[nm]]
  if(is.list(S)) S <- S[[2L]]
  starts[[nm]] <- if('sigma_sq' %in% names(x$priors[[nm]]))
    list(sigma_sq=if(nrow(S)==1L) S[1,1]*multiplier else
      setNames(diag(S)*multiplier,c('(Intercept)','x'))) else list(covariance=S*multiplier)
}
if(x$spec$family=='gaussian') starts$resid <- list(tau_sq=.3*multiplier)
if(spatial) starts$nngp_1 <- list(sigma_sq=.3*multiplier,phi=c(3,6,12)[chain])
set.seed(x$seed+1000L+chain+if(mode=='collapsed') 100000L else 0L)
warnings <- character()
gc()
fit_seconds <- system.time(fit <- withCallingHandlers(stLMM(x$formula,x$data,
  priors=x$priors,starting=starts,family=x$spec$family,
  size=if(x$spec$family=='negative_binomial') 2 else NULL,
  n_samples=ns,n_omp_threads=1L,collapse_iid=mode=='collapsed',
  metropolis=blocking,warmup=TRUE,verbose=FALSE),warning=function(w) {
    warnings <<- c(warnings,conditionMessage(w));invokeRestart('muffleWarning')
  }))['elapsed']
# Native probes use these small backends separately; no full fitted objects saved.
if(chain==1L) saveRDS(fit$backend,file.path(out,paste(case,mode,blocking,'backend.rds',sep='-')))
meta <- fit$term_description$global
fit_mib <- as.numeric(object.size(fit))/1024^2
recovery_seconds <- 0
if(spatial || isTRUE(fit$backend$defer_iid_recovery))
  recovery_seconds <- system.time(fit <- recover(fit,sub_sample=list(start=burn+1L)))['elapsed']
alpha_keep <- if(isTRUE(fit$backend$defer_iid_recovery)) match(keep,fit$recover_iter) else keep
rows <- unique(as.integer(seq(1,nrow(x$data),length.out=6)))
Z <- as.matrix(fit$backend$Z[rows,,drop=FALSE]);X <- fit$backend$X[rows,,drop=FALSE]
fixed_iid <- fit$beta_samples[keep,,drop=FALSE]%*%t(X)+fit$alpha_samples[alpha_keep,,drop=FALSE]%*%t(Z)
predictor <- fixed_iid
if(spatial) predictor <- predictor+fit$w_samples_ordered$nngp_1[,fit$backend$process_terms[[1]]$map[rows],drop=FALSE]
coefs <- unlist(lapply(fit$term_description$random_effects$terms,function(b) rep(b$coefficients,b$n_levels)))
combined <- list()
for(cf in colnames(X)) {
  zz <- Z; zz[,coefs!=cf] <- 0
  if(cf!='(Intercept)') zz <- zz/x$data[[cf]][rows]
  combined[[cf]] <- fit$beta_samples[keep,cf]+fit$alpha_samples[alpha_keep,,drop=FALSE]%*%t(zz)
  colnames(combined[[cf]]) <- paste0('combined_',cf,'[',rows,']')
}
parameters <- cbind(fit$beta_samples[keep,,drop=FALSE],
  tau_sq=if(!is.null(fit$tau_sq_samples)) fit$tau_sq_samples[keep] else NULL,
  fit$iid_sigma_sq_samples[keep,,drop=FALSE],fit$iid_covariance_samples[keep,,drop=FALSE],
  fit$sigma_sq_samples[keep,,drop=FALSE],fit$theta_samples[keep,,drop=FALSE])
truth <- x$beta
if(x$spec$family=='gaussian') truth <- c(truth,.3)
scalar_truth <- covariance_truth <- numeric()
for(b in fit$term_description$random_effects$terms) {
  nm <- if(is.null(b$formula_term)) b$name else b$formula_term
  S <- x$covariance_truth[[nm]]
  if(!is.null(b$by)) S <- S[[b$by_level]]
  if(identical(b$cov_model,'unstructured')) covariance_truth <- c(covariance_truth,S[lower.tri(S,diag=TRUE)]) else {
    j <- match(b$coefficients,c('(Intercept)','x','z'))
    scalar_truth <- c(scalar_truth,if(nrow(S)==1L) S[1,1] else S[j,j])
  }
}
truth <- c(truth,scalar_truth,covariance_truth,if(spatial) c(.3,6))
stopifnot(length(truth)==ncol(parameters))
names(truth) <- colnames(parameters)
colnames(fixed_iid) <- paste0('fixed_iid[',rows,']');colnames(predictor) <- paste0('predictor[',rows,']')
selected <- which(colSums(abs(Z))>0)
alpha <- fit$alpha_samples[alpha_keep,selected,drop=FALSE]
colnames(alpha) <- paste0('raw_alpha[',selected,']')
draws <- cbind(parameters,alpha,do.call(cbind,combined),fixed_iid,predictor)
truth <- c(truth,setNames(rep(NA_real_,ncol(alpha)),colnames(alpha)),
  setNames(as.numeric(x$combined[rows,,drop=FALSE]),colnames(do.call(cbind,combined))),
  setNames(x$fixed_iid[rows],colnames(fixed_iid)),setNames(x$predictor[rows],colnames(predictor)))
stopifnot(identical(names(truth),colnames(draws)),all(is.finite(draws)))
metrics <- data.frame(case=case,mode=mode,chain=chain,blocking=blocking,n=nrow(x$data),
  groups=x$spec$groups,q=fit$backend$q,n_samples=ns,burn=burn,
  fit_seconds=unname(fit_seconds),recovery_seconds=unname(recovery_seconds),
  initialization_seconds=unname(fit$timing$phases['initialization']),
  warmup_seconds=unname(fit$timing$phases['warmup']),
  sampling_seconds=unname(fit$timing$phases['sampling']),
  fit_object_mib=fit_mib,recovered_object_mib=as.numeric(object.size(fit))/1024^2,
  M_nnz=meta$M_nnz,chol_nnz=meta$cholmod_lnz,chol_flops=meta$cholmod_flops,
  ordering=meta$cholmod_ordering,acceptance=fit$covariance_acceptance)
file <- file.path(out,paste(case,mode,blocking,chain,ns,sep='-'))
saveRDS(list(metrics=metrics,draws=draws,truth=truth,adaptive=fit$adaptive_metropolis,
  warnings=warnings,environment=c(extSoftVersion()['BLAS'],Sys.getenv(c('OPENBLAS_NUM_THREADS',
    'OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')))),paste0(file,'.rds'))
write.csv(metrics,paste0(file,'.csv'),row.names=FALSE)
cat(case,mode,blocking,chain,ns,'fit',round(fit_seconds,3),'recovery',round(recovery_seconds,3),'\n')
