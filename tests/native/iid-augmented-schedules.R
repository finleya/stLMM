# Exercise blocking, warmup, missing responses, offsets and saved-draw alignment.
Sys.setenv(STLMM_REFERENCE_ONLY='1')
source('tests/native/latent-reference.R')
Sys.unsetenv('STLMM_REFERENCE_ONLY')
set.seed(77201)
d <- data.frame(y=rep(c(0,1,2),12),x=rnorm(36),off=rep(c(-.2,.3),18),
  group=factor(rep(1:6,each=6),levels=1:7),parent=rep(c('A','B'),each=18),
  child=rep(1:3,each=6,length.out=36),crossed=rep(1:9,4),time=rep(1:12,each=3),
  lon=rep(runif(12),each=3),lat=rep(runif(12),each=3))
d$nested <- interaction(d$parent,d$child,drop=TRUE)
families <- c('binomial','probit','negative_binomial','estimated_size')
blocking <- c('joint','by_term','residual_process','variance_theta','process_variance','scalar')
fields <- c('beta_samples','alpha_samples','sigma_sq_re_samples','iid_covariance_samples',
  'theta_samples','sigma_sq_samples','size_samples','adaptive_metropolis','size_metropolis')
cases <- 0L
for(family in families) for(kind in c('scalar','independent','covariance','crossed','nested','by'))
for(structure in c('none','ar1','nngp')) {
  d$y <- if(family=='probit') rep(c(0,1),18) else rep(c(0,1,2),12)
  d$y[2] <- NA_real_
  form <- as.formula(paste('y~x+offset(off)+',forms[[kind]],switch(structure,
    none='',ar1='+ar1(time)',nngp='+nngp(lon,lat,m=3,cov_model="exp")')))
  pr <- get_priors(kind)
  if(structure=='ar1') pr$ar1_1 <- list(sigma_sq=ig(3,.6),phi=uniform(-.8,.8))
  if(structure=='nngp') pr$nngp_1 <- list(sigma_sq=ig(3,.6),phi=uniform(.5,8))
  if(family=='estimated_size') pr$size <- gamma_dist(2,1)
  fit <- suppressWarnings(stLMM(form,d,priors=pr,
    family=if(family=='estimated_size') 'negative_binomial' else family,
    trials=if(family=='binomial') rep(3,36) else NULL,
    size=if(family=='negative_binomial') 2 else NULL,
    n_samples=1,warmup=FALSE,verbose=FALSE))
  be <- fit$backend;be$collapse_iid <- TRUE;be$n_samples <- 32L
  for(block in seq_along(blocking)) for(warm in c(FALSE,TRUE)) {
    cases <- cases+1L
    be$metropolis_blocking <- as.integer(block-1L)
    be$metropolis_target_accept <- if(blocking[block]=='scalar') .44 else .234
    be$warmup <- list(enabled=warm,batch_length=10L,min_batches=2L,max_batches=2L,
      target=c(.15,.4),near_zero=.01)
    be$recover_process <- structure!='none';be$recover_start <- 3L;be$recover_thin <- 5L
    original <- serialize(be,NULL)
    set.seed(77201+cases)
    out <- withCallingHandlers(.Call('stLMM_collapsed_sampler',be,PACKAGE='stLMM'),
      warning=function(w) stop(conditionMessage(w),call.=FALSE))
    rng <- .Random.seed
    stopifnot(identical(original,serialize(be,NULL)),all(is.finite(out$alpha_samples)),
      all(is.finite(out$beta_samples)),all(is.finite(out$iid_covariance_samples)),
      all(is.finite(out$w_samples)),all(is.finite(out$size_samples)))
    if(warm) stopifnot(out$adaptive_metropolis$warmup$n_batches==2L)
    if(family=='estimated_size') {
      stopifnot(all(out$size_samples>0),out$size_metropolis$warmup$n_attempted==if(warm) 20L else 0L)
    }
    if(structure!='none') {
      stopifnot(identical(out$recover_iter,seq.int(3L,32L,by=5L)))
      be$recover_process <- FALSE
      set.seed(77201+cases)
      unsaved <- .Call('stLMM_collapsed_sampler',be,PACKAGE='stLMM')
      stopifnot(identical(out[intersect(fields,names(out))],unsaved[intersect(fields,names(out))]),
        identical(rng,.Random.seed))
      be$recover_process <- TRUE;be$recover_start <- 1L;be$recover_thin <- 1L
      set.seed(77201+cases)
      full <- .Call('stLMM_collapsed_sampler',be,PACKAGE='stLMM')
      stopifnot(identical(out$w_samples,full$w_samples[out$recover_iter,,drop=FALSE]),
        identical(rng,.Random.seed))
    }
  }
  cat(family,kind,structure,'schedules passed\n');flush.console()
}
# Size-only warmup must still run when every IID covariance is held fixed.
f <- suppressWarnings(stLMM(y~0+iid(group),data.frame(y=c(0,1,2,1,3,0),group='A'),
  family='negative_binomial',priors=list(size=gamma_dist(2,1),iid_1=list(sigma_sq=ig(3,.6))),
  n_samples=1,warmup=FALSE,verbose=FALSE))
b <- f$backend;b$collapse_iid <- TRUE;b$n_samples <- 40L
b$sigma_sq_re_starting <- .4;b$sigma_sq_re_tuning <- 0
b$warmup <- list(enabled=TRUE,batch_length=10L,min_batches=2L,max_batches=2L,
  target=c(.15,.4),near_zero=.01)
out <- .Call('stLMM_collapsed_sampler',b,PACKAGE='stLMM')
stopifnot(length(out$adaptive_metropolis$parameter_labels)==0,
  out$size_metropolis$warmup$n_attempted==20L,all(out$sigma_sq_re_samples==.4),
  all(is.finite(out$alpha_samples)),var(out$size_samples)>0)
saveRDS(list(cases=cases,size_only=TRUE),args[2])
cat(cases,'augmented blocking/warmup schedules passed\n')
