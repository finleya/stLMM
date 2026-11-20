# Verify accept/reject state consistency and all proposal blocking schedules.
args <- commandArgs(TRUE)
stopifnot(length(args)==2L)
.libPaths(c(args[1],.libPaths()))
suppressPackageStartupMessages(library(stLMM))
set.seed(76501)
d <- data.frame(y=rnorm(36),x=rnorm(36),g=factor(rep(letters[1:6],each=6)),
  parent=rep(c('A','B'),each=18),child=rep(1:3,12),time=rep(1:12,each=3))
blocking <- c('joint','by_term','residual_process','variance_theta','process_variance','scalar')
results <- list()
for(kind in c('scalar','independent','by')) for(structure in c(FALSE,TRUE)) {
  form <- switch(kind,scalar='iid(g)',independent='iid(g,~1+x,cov_model="independent")',by='iid(child,by=parent)')
  pr <- if(kind=='independent') list(iid_1=list(sigma_sq=list('(Intercept)'=ig(3,1),x=ig(3,1)))) else
    list(iid_1=list(sigma_sq=ig(3,1)))
  if(structure) pr$ar1_1 <- list(sigma_sq=ig(3,1),phi=uniform(-.8,.8))
  formula <- as.formula(paste('y~x+',form,if(structure) '+ar1(time)' else ''))
  for(block in blocking) {
    f <- suppressWarnings(stLMM(formula,d,priors=pr,n_samples=1,verbose=FALSE,warmup=FALSE,
      metropolis=list(blocking=block,batch_length=10)))
    b <- f$backend;b$collapse_iid <- TRUE
    set.seed(76501+length(results))
    accepted <- .Call('stLMM_iid_mh_probe',b,PACKAGE='stLMM')
    stopifnot(all(accepted>0L))
    for(warm in c(FALSE,TRUE)) {
      b$n_samples <- 120L
      b$warmup <- list(enabled=warm,batch_length=20L,min_batches=2L,max_batches=2L,
                      target=c(.15,.4),near_zero=.01)
      fit <- .Call('stLMM_collapsed_sampler',b,PACKAGE='stLMM')
      stopifnot(all(is.finite(fit$beta_samples)),all(is.finite(fit$alpha_samples)),
        all(is.finite(fit$sigma_sq_re_samples)),all(fit$sigma_sq_re_samples>0),
        length(fit$adaptive_metropolis$parameter_labels)==1+length(b$sigma_sq_re_starting)+2*structure)
      if(warm) stopifnot(fit$adaptive_metropolis$warmup$n_batches==2L)
    }
    results[[paste(kind,structure,block)]] <- accepted
  }
}
saveRDS(results,args[2]);print(colSums(do.call(rbind,results)))
cat(length(results),'proposal-state cases and',2*length(results),'sampling schedules passed\n')
