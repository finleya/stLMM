# Marginal variance posterior checked by independent one-dimensional quadrature.
args <- commandArgs(TRUE)
stopifnot(length(args)==2L)
.libPaths(c(args[1],.libPaths()))
suppressPackageStartupMessages(library(stLMM))
stopifnot(requireNamespace('posterior',quietly=TRUE))
set.seed(76401)
d <- data.frame(g=factor(rep(letters[1:6],each=4),levels=letters[1:7]),
                time=rep(1:8,each=3),x=0,y=rnorm(24))
results <- list()
for(kind in c('iid','ar1','zero')) {
  formula <- switch(kind,iid=y~0+iid(g),ar1=y~0+iid(g)+ar1(time),zero=y~0+x:iid(g))
  start <- list(tau_sq=fixed(.3))
  if(kind=='ar1') start$ar1_1 <- list(sigma_sq=fixed(.4),phi=fixed(.6))
  f <- suppressWarnings(stLMM(formula,d,starting=start,
    priors=list(iid_1=list(sigma_sq=ig(4,2))),n_samples=1,warmup=FALSE,verbose=FALSE))
  b <- f$backend;b$collapse_iid <- TRUE;b$n_samples <- 18000L
  Z <- as.matrix(b$Z_obs); R <- diag(.3,nrow(Z))
  if(kind=='ar1') R <- R + .4*.6^abs(outer(d$time,d$time,'-'))
  Rchol <- t(chol(R)); W <- forwardsolve(Rchol,Z)
  yy <- forwardsolve(Rchol,b$y_model_obs)
  ev <- eigen(tcrossprod(W),symmetric=TRUE)
  lam <- pmax(0,ev$values);proj <- drop(crossprod(ev$vectors,yy))
  target <- function(eta) vapply(eta,function(e) {
    v<-exp(e)
    -.5*sum(log1p(v*lam)+proj^2/(1+v*lam))-4*e-2/v
  },numeric(1))
  mode <- optimize(target,c(-10,10),maximum=TRUE)$objective
  den <- function(e) exp(target(e)-mode)
  norm <- integrate(den,-20,20,rel.tol=1e-11)$value
  exact <- c(variance=integrate(function(e) exp(e)*den(e),-20,20,rel.tol=1e-10)$value/norm,
    log_variance=integrate(function(e) e*den(e),-20,20,rel.tol=1e-10)$value/norm)
  samples <- matrix(NA_real_,16000,3)
  for(chain in 1:3) {
    set.seed(76401+10*match(kind,c('iid','ar1','zero'))+chain)
    b$sigma_sq_re_starting[] <- c(.2,1,2)[chain]
    out <- .Call('stLMM_collapsed_sampler',b,PACKAGE='stLMM')
    samples[,chain] <- out$sigma_sq_re_samples[2001:18000,1]
    stopifnot(all(is.finite(samples[,chain])),all(samples[,chain]>0))
  }
  means <- c(mean(samples),mean(log(samples)))
  se <- c(posterior::mcse_mean(samples),posterior::mcse_mean(log(samples)))
  error <- abs(means-exact)/se
  stopifnot(max(error)<6,posterior::rhat(samples)<1.02)
  results[[kind]] <- c(exact,estimate=means,mcse=se,standardized_error=error,
    rhat=posterior::rhat(samples),ess=posterior::ess_bulk(samples))
  b$sigma_sq_re_starting[] <- .73
  b$sigma_sq_re_tuning <- 0
  b$n_samples <- 5L
  held <- .Call('stLMM_collapsed_sampler',b,PACKAGE='stLMM')
  stopifnot(all(held$sigma_sq_re_samples==.73))
}
saveRDS(results,args[2]);print(do.call(rbind,results))
cat(length(results),'scalar covariance quadrature cases passed\n')
