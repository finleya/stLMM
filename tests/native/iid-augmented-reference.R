# Compare the actual augmented sampler with ordinary-likelihood quadrature.
args <- commandArgs(TRUE)
stopifnot(length(args)==2L)
.libPaths(c(args[1],.libPaths()))
suppressPackageStartupMessages(library(stLMM))
stopifnot(requireNamespace('posterior',quietly=TRUE))
quadrature <- function(n,lo,hi) {
  j <- seq_len(n-1);J <- matrix(0,n,n)
  J[cbind(j,j+1)] <- J[cbind(j+1,j)] <- j/sqrt(4*j*j-1)
  e <- eigen(J,symmetric=TRUE)
  list(x=(lo+hi)/2+(hi-lo)/2*e$values,w=(hi-lo)*e$vectors[1,]^2)
}
results <- list()
for(formulation in c('scalar','covariance','fixed_beta'))
for(family in c('binomial','probit','negative_binomial','estimated_size')) {
  correlated <- formulation=='covariance';fixed_beta <- formulation=='fixed_beta'
  if((correlated || fixed_beta) && family=='estimated_size') next
  d <- data.frame(g=factor(rep('A',12)),x=seq(-1.5,1.5,length.out=12),
    off=rep(c(-.3,.2,.4),4),y=c(0,1,0,2,1,3,0,2,1,3,2,1))
  if(family=='probit') d$y <- as.numeric(d$y>1)
  form <- if(correlated) y~0+iid(g,~1+x)+offset(off) else y~0+iid(g)+offset(off)
  if(fixed_beta) form <- y~0+x+iid(g)+offset(off)
  Psi <- matrix(c(.6,.15,.15,.4),2);df <- 6
  pr <- list(iid_1=if(correlated) list(covariance=iw(df,Psi)) else list(sigma_sq=ig(3,.6)))
  if(fixed_beta) pr$beta <- normal(.2,.7)
  if(family=='estimated_size') pr$size <- gamma_dist(2,1)
  f <- suppressWarnings(stLMM(form,d,priors=pr,
    family=if(family=='estimated_size') 'negative_binomial' else family,
    trials=if(family=='binomial') rep(3,12) else NULL,
    size=if(family=='negative_binomial') 2 else NULL,n_samples=1,warmup=FALSE,verbose=FALSE))
  reference <- function(nq) {
    a <- quadrature(nq,-12,12)
    two <- correlated || fixed_beta || family=='estimated_size'
    b <- if(correlated || fixed_beta) a else if(two) quadrature(nq,-10,5) else list(x=0,w=1)
    grid <- expand.grid(a=a$x,b=b$x)
    weight <- as.vector(outer(a$w,b$w))
    P <- solve(Psi)
    ll <- if(correlated) -(df+1)/2*log1p(P[1,1]*grid$a^2+2*P[1,2]*grid$a*grid$b+P[2,2]*grid$b^2) else
      dt(grid$a/sqrt(.6/3),df=6,log=TRUE)
    r <- if(family=='estimated_size') exp(grid$b) else 2
    if(family=='estimated_size') ll <- ll+dgamma(r,2,1,log=TRUE)+grid$b
    if(fixed_beta) ll <- ll+dnorm(grid$b,.2,.7,log=TRUE)
    for(i in seq_len(nrow(d))) {
      eta <- grid$a+d$off[i]+if(correlated || fixed_beta) d$x[i]*grid$b else 0
      ll <- ll+switch(family,binomial=dbinom(d$y[i],3,plogis(eta),log=TRUE),
        probit=pnorm(if(d$y[i]==1) eta else -eta,log.p=TRUE),
        dnbinom(d$y[i],size=r,mu=exp(eta),log=TRUE))
    }
    weight <- weight*exp(ll-max(ll));weight <- weight/sum(weight)
    moments <- if(correlated) cbind(grid$a,grid$b,
      (Psi[1,1]+grid$a^2)/(df-2),(Psi[2,1]+grid$a*grid$b)/(df-2),
      (Psi[2,2]+grid$b^2)/(df-2)) else
      cbind(grid$a,grid$a^2,(.6+grid$a^2/2)/2.5)
    if(fixed_beta) moments[,2] <- grid$b
    if(family=='estimated_size') moments <- cbind(moments,exp(grid$b),grid$b)
    drop(crossprod(weight,moments))
  }
  exact <- reference(360L);fine <- reference(480L)
  stopifnot(max(abs(exact-fine))<1e-7)
  be <- f$backend;be$collapse_iid <- TRUE;be$n_samples <- 24000L
  be$warmup <- list(enabled=TRUE,batch_length=25L,min_batches=3L,max_batches=20L,
    target=c(.15,.4),near_zero=.01)
  draws <- array(NA_real_,c(20000,3,length(exact)))
  for(chain in 1:3) {
    set.seed(77100+10*match(family,c('binomial','probit','negative_binomial','estimated_size'))+chain+100*correlated)
    if(correlated) be$iid_covariance[[1]]$starting <- Psi*c(.5,1,2)[chain] else
      be$sigma_sq_re_starting[] <- c(.1,.4,1)[chain]
    if(family=='estimated_size') be$nb_size <- c(.5,2,5)[chain]
    out <- .Call('stLMM_collapsed_sampler',be,PACKAGE='stLMM');keep <- 4001:24000
    dr <- if(correlated) cbind(out$alpha_samples[keep,],out$iid_covariance_samples[keep,]) else
      cbind(out$alpha_samples[keep,1],out$alpha_samples[keep,1]^2,out$sigma_sq_re_samples[keep,1])
    if(fixed_beta) dr[,2] <- out$beta_samples[keep,1]
    if(family=='estimated_size') dr <- cbind(dr,out$size_samples[keep],log(out$size_samples[keep]))
    draws[,chain,] <- dr
  }
  estimate <- apply(draws,3,mean);se <- apply(draws,3,posterior::mcse_mean)
  err <- abs(estimate-fine)/se;rhat <- apply(draws,3,posterior::rhat)
  results[[paste(family,formulation)]] <- list(exact=fine,estimate=estimate,mcse=se,
    error=err,rhat=rhat,ess=apply(draws,3,posterior::ess_bulk))
  saveRDS(results,args[2])
  cat(family,formulation,'maximum error/MCSE',max(err),'maximum R-hat',max(rhat),'\n')
  flush.console()
  stopifnot(max(err)<6,max(rhat)<1.02)
}
cat(length(results),'augmented ordinary-likelihood quadrature cases passed\n')
