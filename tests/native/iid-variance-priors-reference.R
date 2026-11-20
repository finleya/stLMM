# Independent density/Jacobian, proposal-state and Gaussian quadrature checks.
# Rscript tests/native/iid-variance-priors-reference.R PROBE_LIBRARY OUTPUT.rds
args <- commandArgs(TRUE);stopifnot(length(args)==2L)
.libPaths(c(args[1],.libPaths()));suppressPackageStartupMessages(library(stLMM))
stopifnot(requireNamespace('posterior',quietly=TRUE))
priors <- list(ig=ig(4,2),log_normal=log_normal(-.4,.6),gamma=gamma_dist(2.5,1.5),
  half_normal=half_normal(.8),half_t=half_t(7,.7))
# Normalized R densities evaluated independently on the variance scale.
log_prior <- function(v,name) switch(name,
  ig=dgamma(1/v,4,rate=2,log=TRUE)-2*log(v),
  log_normal=dlnorm(v,-.4,.6,log=TRUE),gamma=dgamma(v,2.5,rate=1.5,log=TRUE),
  half_normal=dnorm(sqrt(v),0,.8,log=TRUE)-.5*log(v),
  half_t=dt(sqrt(v)/.7,df=7,log=TRUE)-log(.7)-.5*log(v))
set.seed(89301)
d <- data.frame(g=factor(rep(letters[1:6],each=4),levels=letters[1:7]),
  time=rep(1:8,each=3),x=0)
a <- rnorm(7,sd=sqrt(.7));d$y <- a[as.integer(d$g)]+rnorm(24,sd=sqrt(.3))
results <- list();density_error <- numeric();proposals <- list()
for(name in names(priors)) {
  f <- stLMM(y~0+iid(g),d,priors=list(iid_1=list(sigma_sq=priors[[name]])),
    starting=list(tau_sq=fixed(.3)),n_samples=1,warmup=FALSE,verbose=FALSE)
  be <- f$backend;be$probe_alpha <- as.double(seq_len(be$q)/9)
  eta <- seq(-12,6,length.out=81);v <- exp(eta)
  actual <- .Call('stLMM_iid_density_probe',be,eta,PACKAGE='stLMM')
  expected <- cbind(log_prior(v,name)+eta,log_prior(v,name)+eta+
    vapply(v,function(x) sum(dnorm(be$probe_alpha,0,sqrt(x),log=TRUE)),numeric(1)))
  error <- abs(sweep(actual,2,actual[41,])-sweep(expected,2,expected[41,]))/(1+abs(expected))
  stopifnot(max(error)<1e-12)
  density_error[name] <- max(error)
  extreme <- .Call('stLMM_iid_density_probe',be,c(-2000,2000),PACKAGE='stLMM')
  stopifnot(all(extreme == -Inf))
  if(name=='ig') next
  for(collapse in c(FALSE,TRUE)) for(block in c('joint','by_term','residual_process',
      'variance_theta','process_variance','scalar')) {
    f <- stLMM(y~iid(g)+ar1(time),d,
      priors=list(iid_1=list(sigma_sq=priors[[name]]),
        ar1_1=list(sigma_sq=ig(3,1),phi=uniform(-.8,.8))),
      collapse_iid=collapse,metropolis=block,n_samples=1,warmup=FALSE,verbose=FALSE)
    set.seed(89310+length(proposals))
    z <- .Call('stLMM_iid_mh_probe',f$backend,PACKAGE='stLMM')
    stopifnot(all(z>0))
    proposals[[paste(name,collapse,block)]] <- z
  }
}
for(kind in c('iid','ar1','zero')) for(name in names(priors)) {
  formula <- switch(kind,iid=y~0+iid(g),ar1=y~0+iid(g)+ar1(time),zero=y~0+x:iid(g))
  start <- list(tau_sq=fixed(.3))
  if(kind=='ar1') start$ar1_1 <- list(sigma_sq=fixed(.4),phi=fixed(.6))
  f <- stLMM(formula,d,starting=start,priors=list(iid_1=list(sigma_sq=priors[[name]])),
    n_samples=1,warmup=FALSE,verbose=FALSE)
  be <- f$backend;Z <- as.matrix(be$Z_obs);R <- diag(.3,nrow(Z))
  if(kind=='ar1') R <- R+.4*.6^abs(outer(d$time,d$time,'-'))
  L <- t(chol(R));W <- forwardsolve(L,Z);yy <- forwardsolve(L,be$y_model_obs)
  e <- svd(W,nv=0);lam <- e$d^2
  proj <- drop(crossprod(e$u,yy))
  target <- function(eta) vapply(eta,function(z) {
    v<-exp(z);-.5*sum(log1p(v*lam)+proj^2/(1+v*lam))+log_prior(v,name)+z
  },numeric(1))
  shift <- optimize(target,c(-15,8),maximum=TRUE)$objective
  den <- function(z) exp(target(z)-shift)
  normalizer <- integrate(den,-60,20,rel.tol=1e-10)$value
  exact <- c(variance=integrate(function(z)exp(z)*den(z),-60,20,rel.tol=1e-9)$value/normalizer,
    log_variance=integrate(function(z)z*den(z),-60,20,rel.tol=1e-9)$value/normalizer)
  for(collapse in c(FALSE,TRUE)) {
    be$collapse_iid <- collapse;be$n_samples <- 24000L
    be$warmup <- list(enabled=TRUE,batch_length=25L,min_batches=3L,max_batches=20L,
      target=c(.15,.4),near_zero=.01)
    samples <- matrix(NA_real_,18000,3)
    for(ch in 1:3) {
      set.seed(89400+100*match(kind,c('iid','ar1','zero'))+10*match(name,names(priors))+ch+1000*collapse)
      be$sigma_sq_re_starting[] <- c(.1,.8,2)[ch]
      fit <- .Call('stLMM_collapsed_sampler',be,PACKAGE='stLMM')
      samples[,ch] <- fit$sigma_sq_re_samples[6001:24000,1]
    }
    dr <- list(variance=samples,log_variance=log(samples))
    estimate <- vapply(dr,mean,numeric(1));se <- vapply(dr,posterior::mcse_mean,numeric(1))
    error <- abs(estimate-exact)/se;rhat <- vapply(dr,posterior::rhat,numeric(1))
    id <- paste(kind,name,collapse)
    results[[id]] <- list(exact=exact,estimate=estimate,mcse=se,error=error,rhat=rhat,
      ess=vapply(dr,posterior::ess_bulk,numeric(1)),interval=quantile(samples,c(.025,.975)))
    saveRDS(list(density_error=density_error,proposals=proposals,posterior=results),args[2])
    cat(id,'max error/MCSE',max(error),'R-hat',max(rhat),'\n');flush.console()
    stopifnot(max(error)<6,max(rhat)<1.02)
  }
}
cat(length(density_error),'prior densities,',length(proposals),'proposal-state cases and',
    length(results),'Gaussian posterior references passed\n')
