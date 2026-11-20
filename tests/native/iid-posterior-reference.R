# Gaussian posterior comparisons under matched IG/IW priors.
# Uses the backend control until the public collapsing interface is enabled.
Sys.setenv(STLMM_REFERENCE_ONLY='1')
source('tests/native/latent-reference.R')
Sys.unsetenv('STLMM_REFERENCE_ONLY')
stopifnot(requireNamespace('posterior',quietly=TRUE))
set.seed(76701)
n <- 360L
d <- data.frame(y=rnorm(n),x=rnorm(n),group=factor(rep(seq_len(30),each=12)),
  parent=rep(c('A','B','C'),each=120),child=rep(seq_len(10),each=12,length.out=n),
  crossed=rep(seq_len(12),30),time=rep(seq_len(60),each=6),
  lon=rep(runif(60),each=6),lat=rep(runif(60),each=6))
d$nested <- interaction(d$parent,d$child,drop=TRUE)
cases <- expand.grid(kind=c('scalar','independent','covariance','crossed','nested','by'),
  structure=c('none','ar1','nngp'),stringsAsFactors=FALSE)
results <- list()
for(cc in seq_len(nrow(cases))) {
  kind <- cases$kind[cc];structure <- cases$structure[cc]
  form <- as.formula(paste('y~x+',forms[[kind]],switch(structure,none='',
    ar1='+ar1(time)',nngp='+nngp(lon,lat,m=8,cov_model="exp")')))
  pr <- get_priors(kind)
  st <- list(tau_sq=.25)
  if(structure=='ar1') {
    pr$ar1_1 <- list(sigma_sq=ig(3,.6),phi=uniform(-.9,.9))
    st$ar1_1 <- list(sigma_sq=.3,phi=.5)
  }
  if(structure=='nngp') {
    pr$nngp_1 <- list(sigma_sq=ig(3,.6),phi=uniform(.5,10))
    st$nngp_1 <- list(sigma_sq=.3,phi=3)
  }
  f <- suppressWarnings(stLMM(form,d,priors=pr,starting=st,n_samples=1,
    warmup=FALSE,verbose=FALSE))
  be <- f$backend
  for(b in seq_along(be$iid_covariance)) {
    if(is.null(be$iid_covariance[[b]])) be$sigma_sq_re_starting[b] <- .5 else
      be$iid_covariance[[b]]$starting <- matrix(c(.7,.25,.25,.35),2)
  }
  ref <- reference(be)
  set.seed(76701+cc)
  latent <- drop(t(chol(ref$K))%*%rnorm(ncol(ref$K)))
  be$y_model_obs <- drop(be$X_obs%*%c(1,.4)+ref$H%*%latent+rnorm(n,sd=.5))
  be$y_obs <- be$y_model_obs
  be$n_samples <- 12000L
  be$warmup <- list(enabled=TRUE,batch_length=50L,min_batches=3L,max_batches=20L,
    target=c(.15,.4),near_zero=.01)
  be$recover_process <- structure!='none';be$recover_start <- 3001L;be$recover_thin <- 1L
  keep <- 3001:12000
  scalar <- which(vapply(be$iid_covariance,is.null,logical(1)))
  covariance_truth <- unlist(lapply(be$iid_covariance,function(b) {
    if(is.null(b)) return(numeric())
    b$starting[lower.tri(b$starting,diag=TRUE)]
  }),use.names=FALSE)
  covariance_names <- unlist(lapply(seq_along(be$iid_covariance),function(b) {
    cb <- be$iid_covariance[[b]]
    if(is.null(cb)) return(character())
    ij <- which(lower.tri(cb$starting,diag=TRUE),arr.ind=TRUE)
    paste0(names(be$iid_covariance)[b],'_cov[',ij[,1],',',ij[,2],']')
  }),use.names=FALSE)
  parameter_names <- c(colnames(be$X_obs),'tau_sq',
    if(length(scalar)) paste0(names(be$iid_covariance)[scalar],'_sigma_sq'),
    covariance_names,if(structure!='none') paste0(structure,c('_sigma_sq','_phi')),
    'fixed_plus_iid_first','fixed_plus_iid_last','predictor_first','predictor_last')
  nw <- ncol(ref$H)-be$q
  truth <- c(1,.4,.25,be$sigma_sq_re_starting[scalar],covariance_truth,
    if(structure!='none') c(.3,if(structure=='ar1') .5 else 3),
    drop(be$X_obs[c(1,n),,drop=FALSE]%*%c(1,.4)+
         be$Z_obs[c(1,n),,drop=FALSE]%*%latent[nw+seq_len(be$q)]),
    drop(be$X_obs[c(1,n),,drop=FALSE]%*%c(1,.4)+ref$H[c(1,n),,drop=FALSE]%*%latent))
  modes <- list()
  for(joint in c(FALSE,TRUE)) {
    draws <- NULL;elapsed <- numeric(3);fill <- numeric(3)
    be$collapse_iid <- joint
    for(chain in 1:3) {
      set.seed(76801+cc*10+chain+1000*joint)
      elapsed[chain] <- system.time(out <- .Call('stLMM_collapsed_sampler',be,PACKAGE='stLMM'))['elapsed']
      base <- out$beta_samples[keep,,drop=FALSE]%*%t(be$X_obs[c(1,n),,drop=FALSE])+
        out$alpha_samples[keep,,drop=FALSE]%*%t(as.matrix(be$Z_obs[c(1,n),,drop=FALSE]))
      predictor <- base
      if(nw) predictor <- predictor+out$w_samples%*%t(ref$H[c(1,n),seq_len(nw),drop=FALSE])
      dr <- cbind(out$beta_samples[keep,,drop=FALSE],tau_sq=out$tau_sq_samples[keep],
        out$sigma_sq_re_samples[keep,scalar,drop=FALSE],out$iid_covariance_samples[keep,,drop=FALSE],
        out$sigma_sq_samples[keep,,drop=FALSE],out$theta_samples[keep,,drop=FALSE],base,predictor)
      if(is.null(draws)) draws <- array(NA_real_,c(nrow(dr),3,ncol(dr)),
        dimnames=list(NULL,NULL,parameter_names))
      draws[,chain,] <- dr
      fill[chain] <- out$term_description$global$cholmod_lnz
    }
    summary <- posterior::summarise_draws(posterior::as_draws_array(draws),
      'mean','mcse_mean','rhat','ess_bulk','ess_tail',
      q025=function(x) unname(quantile(x,.025)),q975=function(x) unname(quantile(x,.975)))
    summary$truth <- truth;summary$seconds <- sum(elapsed)
    summary$ess_per_second <- summary$ess_bulk/sum(elapsed)
    name <- if(joint) 'collapsed' else 'separate'
    modes[[name]] <- list(summary=summary,seconds=elapsed,chol_nnz=fill)
    cat(kind,structure,name,'seconds',round(sum(elapsed),2),'max R-hat',round(max(summary$rhat),4),
      'min ESS',round(min(summary$ess_bulk)),'\n');flush.console()
  }
  a <- modes$separate$summary;b <- modes$collapsed$summary
  err <- abs(a$mean-b$mean)/sqrt(a$mcse_mean^2+b$mcse_mean^2)
  modes$mean_difference_mcse <- err
  results[[paste(kind,structure,sep='/')]] <- modes
  saveRDS(results,args[2])
  stopifnot(all(is.finite(err)),max(err)<6,max(b$rhat)<1.05)
}
cat(length(results),'estimated-covariance posterior comparisons passed\n')
