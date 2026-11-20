# Verify accept/reject state consistency and all proposal blocking schedules.
args <- commandArgs(TRUE)
stopifnot(length(args)==2L)
.libPaths(c(args[1],.libPaths()))
suppressPackageStartupMessages(library(stLMM))
set.seed(76501)
d <- data.frame(y=rnorm(36),x=rnorm(36),g=factor(rep(letters[1:6],each=6)),
  parent=rep(c('A','B'),each=18),child=rep(1:3,12),time=rep(1:12,each=3))
for(j in 1:4) d[[paste0('z',j)]] <- rnorm(nrow(d))
blocking <- c('joint','by_term','residual_process','variance_theta','process_variance','scalar')
results <- list()
for(kind in c('scalar','independent','by','covariance','covariance_by','mixed','dimension6')) for(structure in c(FALSE,TRUE)) {
  form <- switch(kind,scalar='iid(g)',independent='iid(g,~1+x,cov_model="independent")',by='iid(child,by=parent)',
    covariance='iid(g,~1+x)',covariance_by='iid(child,~1+x,by=parent)',
    mixed='iid(g,~1+x)+iid(parent)',dimension6='iid(g,~1+x+z1+z2+z3+z4)')
  pr <- if(kind=='independent') list(iid_1=list(sigma_sq=list('(Intercept)'=ig(3,1),x=ig(3,1)))) else
    list(iid_1=list(sigma_sq=ig(3,1)))
  if(kind %in% c('covariance','covariance_by','mixed')) pr <- list(iid_1=list(covariance=iw(5,diag(2))))
  if(kind=='mixed') pr$iid_2 <- list(sigma_sq=ig(3,1))
  if(kind=='dimension6') pr <- list(iid_1=list(covariance=iw(9,diag(6))))
  if(structure) pr$ar1_1 <- list(sigma_sq=ig(3,1),phi=uniform(-.8,.8))
  formula <- as.formula(paste('y~x+',form,if(structure) '+ar1(time)' else ''))
  for(block in blocking) {
    f <- suppressWarnings(stLMM(formula,d,priors=pr,n_samples=1,verbose=FALSE,warmup=FALSE,
      metropolis=list(blocking=block,batch_length=10)))
    b <- f$backend;b$collapse_iid <- TRUE
    set.seed(76501+length(results))
    accepted <- withCallingHandlers(.Call('stLMM_iid_mh_probe',b,PACKAGE='stLMM'),
      warning=function(w) stop(conditionMessage(w),call.=FALSE))
    stopifnot(all(accepted>0L))
    for(warm in c(FALSE,TRUE)) {
      b$n_samples <- 120L
      b$warmup <- list(enabled=warm,batch_length=20L,min_batches=2L,max_batches=2L,
                      target=c(.15,.4),near_zero=.01)
      fit <- withCallingHandlers(.Call('stLMM_collapsed_sampler',b,PACKAGE='stLMM'),
        warning=function(w) stop(conditionMessage(w),call.=FALSE))
      stopifnot(all(is.finite(fit$beta_samples)),all(is.finite(fit$alpha_samples)),
        all(is.finite(fit$sigma_sq_re_samples)),all(fit$sigma_sq_re_samples>0),
        length(fit$adaptive_metropolis$parameter_labels)==1+sum(vapply(b$iid_covariance,function(cb) if(is.null(cb)) 1 else cb$dim*(cb$dim+1)/2,numeric(1)))+2*structure)
      if(warm) stopifnot(fit$adaptive_metropolis$warmup$n_batches==2L)
      stopifnot(all(vapply(fit$term_description$random_effects$terms,
        function(term) identical(term$sampler$status,'collapsed Metropolis'),logical(1))))
    }
    results[[paste(kind,structure,block)]] <- accepted
  }
}
saveRDS(results,args[2]);print(colSums(do.call(rbind,results)))
cat(length(results),'proposal-state cases and',2*length(results),'sampling schedules passed\n')
