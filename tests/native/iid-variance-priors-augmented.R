# Ordinary-likelihood quadrature for nonconjugate IID variance priors.
# Integrate standardized alpha = sqrt(variance) * z, with z ~ N(0,1).
args <- commandArgs(TRUE);stopifnot(length(args)==2L)
.libPaths(c(args[1],.libPaths()));suppressPackageStartupMessages(library(stLMM))
stopifnot(requireNamespace('posterior',quietly=TRUE))
priors <- list(log_normal=log_normal(-.4,.6),gamma=gamma_dist(2.5,1.5),
  half_normal=half_normal(.8),half_t=half_t(7,.7))
log_prior <- function(v,name) switch(name,
  log_normal=dlnorm(v,-.4,.6,log=TRUE),gamma=dgamma(v,2.5,rate=1.5,log=TRUE),
  half_normal=dnorm(sqrt(v),0,.8,log=TRUE)-.5*log(v),
  half_t=dt(sqrt(v)/.7,df=7,log=TRUE)-log(.7)-.5*log(v))
quadrature <- function(n,lo,hi) {
  j <- seq_len(n-1);J <- matrix(0,n,n)
  J[cbind(j,j+1)] <- J[cbind(j+1,j)] <- j/sqrt(4*j*j-1)
  e <- eigen(J,symmetric=TRUE)
  list(x=(lo+hi)/2+(hi-lo)/2*e$values,w=(hi-lo)*e$vectors[1,]^2)
}
results <- list()
for(family in c('binomial','probit','negative_binomial')) for(name in names(priors)) {
  d <- data.frame(g=factor(rep('A',12)),off=rep(c(-.3,.2,.4),4),
    y=c(0,1,0,2,1,3,0,2,1,3,2,1))
  if(family=='probit') d$y <- as.numeric(d$y>1)
  reference <- function(nq) {
    z <- quadrature(nq,-9,9);eta <- quadrature(nq,-50,8)
    grid <- expand.grid(z=z$x,eta=eta$x);v <- exp(grid$eta);a <- sqrt(v)*grid$z
    ll <- dnorm(grid$z,log=TRUE)+log_prior(v,name)+grid$eta
    for(i in seq_len(nrow(d))) {
      mu <- a+d$off[i]
      ll <- ll+switch(family,binomial=dbinom(d$y[i],3,plogis(mu),log=TRUE),
        probit=pnorm(if(d$y[i]) mu else -mu,log.p=TRUE),
        negative_binomial=dnbinom(d$y[i],size=2,mu=exp(mu),log=TRUE))
    }
    weight <- as.vector(outer(z$w,eta$w))*exp(ll-max(ll));weight <- weight/sum(weight)
    drop(crossprod(weight,cbind(alpha=a,alpha_sq=a^2,variance=v,log_variance=grid$eta)))
  }
  exact <- reference(600L);fine <- reference(840L)
  stopifnot(max(abs(exact-fine))<2e-6)
  f <- stLMM(y~0+iid(g)+offset(off),d,family=family,
    priors=list(iid_1=list(sigma_sq=priors[[name]])),
    trials=if(family=='binomial') rep(3,12) else NULL,
    size=if(family=='negative_binomial') 2 else NULL,
    n_samples=1,warmup=FALSE,verbose=FALSE)
  for(collapse in c(FALSE,TRUE)) {
    be <- f$backend;be$collapse_iid <- collapse;be$n_samples <- 24000L
    be$warmup <- list(enabled=TRUE,batch_length=25L,min_batches=3L,max_batches=20L,
      target=c(.15,.4),near_zero=.01)
    draws <- array(NA_real_,c(18000,3,4))
    for(ch in 1:3) {
      set.seed(89700+100*match(family,c('binomial','probit','negative_binomial'))+
        10*match(name,names(priors))+ch+1000*collapse)
      be$sigma_sq_re_starting[] <- c(.1,.8,2)[ch]
      out <- .Call('stLMM_collapsed_sampler',be,PACKAGE='stLMM');keep <- 6001:24000
      a <- out$alpha_samples[keep,1];v <- out$sigma_sq_re_samples[keep,1]
      draws[,ch,] <- cbind(a,a^2,v,log(v))
    }
    estimate <- apply(draws,3,mean);se <- apply(draws,3,posterior::mcse_mean)
    error <- abs(estimate-fine)/se;rhat <- apply(draws,3,posterior::rhat)
    id <- paste(family,name,collapse)
    results[[id]] <- list(exact=fine,estimate=estimate,mcse=se,error=error,rhat=rhat,
      ess=apply(draws,3,posterior::ess_bulk))
    saveRDS(results,args[2]);cat(id,'max error/MCSE',max(error),'R-hat',max(rhat),'\n');flush.console()
    stopifnot(max(error)<6,max(rhat)<1.02)
  }
}
cat(length(results),'augmented posterior references passed\n')
