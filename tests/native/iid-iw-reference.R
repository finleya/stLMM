# IW density in log-Cholesky coordinates and prior-only distribution checks.
args <- commandArgs(TRUE)
stopifnot(length(args)==2L)
.libPaths(c(args[1],.libPaths()))
suppressPackageStartupMessages(library(stLMM))
stopifnot(requireNamespace('posterior',quietly=TRUE))
set.seed(76601)
results <- list();density_error <- jacobian_error <- 0
for(k in c(1L,2L,3L,6L)) {
  d <- data.frame(y=rnorm(12),g=factor(rep(letters[1:4],each=3)))
  for(j in seq_len(k)) d[[paste0('x',j)]] <- 0
  form <- as.formula(paste0('y~0+iid(g,~0+',paste0('x',seq_len(k),collapse='+'),')'))
  scale <- diag(seq(1,2,length.out=k));scale[lower.tri(scale)] <- .1
  scale[upper.tri(scale)] <- t(scale)[upper.tri(scale)]
  df <- k+7
  f <- stLMM(form,d,starting=list(tau_sq=fixed(.3)),priors=list(iid_1=list(covariance=iw(df,scale))),
    n_samples=1,warmup=FALSE,verbose=FALSE)
  b <- f$backend
  ij <- which(lower.tri(matrix(0,k,k),diag=TRUE),arr.ind=TRUE)
  # Singular and indefinite inputs must error, without repairing the matrices.
  for(field in c('scale','starting')) for(diagonal in c(0,-1)) {
    bad <- b;bad$collapse_iid <- TRUE
    bad$iid_covariance[[1]][[field]] <- diag(k)
    bad$iid_covariance[[1]][[field]][k,k] <- diagonal
    err <- tryCatch(.Call('stLMM_collapsed_sampler',bad,PACKAGE='stLMM'),error=identity)
    stopifnot(inherits(err,'error'),grepl('positive definite',conditionMessage(err)))
  }
  if(k>1L) {
    # An invalid current target must error rather than freeze the chain in place.
    bad <- b;bad$collapse_iid <- TRUE
    bad$iid_covariance[[1]]$scale <- diag(1e308,k)
    bad$iid_covariance[[1]]$starting <- diag(k)
    err <- tryCatch(.Call('stLMM_collapsed_sampler',bad,PACKAGE='stLMM'),error=identity)
    stopifnot(inherits(err,'error'),grepl('non-finite.*current state',conditionMessage(err)))
  }
  # Zero tuning holds the supplied covariance exactly, including off-diagonals.
  held <- b;held$collapse_iid <- TRUE;held$n_samples <- 30L
  held$iid_covariance[[1]]$tuning <- 0
  held$iid_covariance[[1]]$starting <- scale
  out <- .Call('stLMM_collapsed_sampler',held,PACKAGE='stLMM')
  stopifnot(identical(unname(out$iid_covariance_samples),
    matrix(rep(scale[ij],each=30),30,nrow(ij))),
    identical(out$term_description$random_effects$terms[[1]]$sampler$status,'fixed'))
  unpack <- function(eta) {
    L <- matrix(0,k,k);L[ij] <- eta;diag(L) <- exp(diag(L));L
  }
  for(point in 1:8) {
    eta <- rnorm(nrow(ij),sd=.25);L <- unpack(eta);Sigma <- tcrossprod(L);P <- solve(Sigma)
    logdet <- as.numeric(determinant(Sigma,logarithm=TRUE)$modulus)
    jac <- k*log(2)+sum((k-seq_len(k)+2)*log(diag(L)))
    expected <- -.5*((df+k+1)*logdet+sum(scale*P))+jac
    out <- .Call('stLMM_iid_density_probe',b,as.double(eta),PACKAGE='stLMM')
    stopifnot(isTRUE(out$valid))
    error <- max(abs(out$covariance-Sigma),abs(out$precision-P),abs(out$cholesky-L),
      abs(out$logdet-logdet),abs(out$target-expected))
    density_error <- max(density_error,error)
    jac_num <- vapply(seq_along(eta),function(j) {
      lo <- hi <- eta;lo[j]<-lo[j]-1e-6;hi[j]<-hi[j]+1e-6
      (tcrossprod(unpack(hi))[ij]-tcrossprod(unpack(lo))[ij])/2e-6
    },numeric(length(eta)))
    jac_num <- matrix(jac_num,length(eta),length(eta))
    jacobian_error <- max(jacobian_error,abs(as.numeric(determinant(jac_num,logarithm=TRUE)$modulus)-jac))
    stopifnot(error<1e-9,jacobian_error<1e-7)
  }
  for(value in c(-1000,1000)) {
    eta <- rep(0,nrow(ij));eta[1]<-value
    stopifnot(!.Call('stLMM_iid_density_probe',b,eta,PACKAGE='stLMM')$valid)
  }
  if(k==6L) next
  b$collapse_iid <- TRUE;b$n_samples <- 24000L
  # In the zero-design model, integrating over coefficients leaves the IW prior.
  nc <- nrow(ij);stats <- array(NA_real_,c(22000,3,2*nc+1))
  for(chain in 1:3) {
    set.seed(76601+10*k+chain)
    out <- .Call('stLMM_collapsed_sampler',b,PACKAGE='stLMM')
    dr <- out$iid_covariance_samples[2001:24000,,drop=FALSE]
    stats[,chain,seq_len(nc)] <- dr
    for(i in seq_len(nrow(dr))) {
      S <- matrix(0,k,k);S[ij]<-dr[i,];S[upper.tri(S)]<-t(S)[upper.tri(S)]
      U <- chol(S)
      stats[i,chain,nc+seq_len(nc)] <- chol2inv(U)[ij]
      stats[i,chain,2*nc+1] <- 2*sum(log(diag(U)))
    }
  }
  exact <- c((scale/(df-k-1))[ij],(df*solve(scale))[ij],
    as.numeric(determinant(scale,logarithm=TRUE)$modulus)-sum(digamma((df+1-seq_len(k))/2))-k*log(2))
  estimate <- apply(stats,3,mean);se <- apply(stats,3,posterior::mcse_mean)
  errors <- abs(estimate-exact)/se
  stopifnot(max(errors)<6,max(apply(stats,3,posterior::rhat))<1.02)
  results[[paste0('dimension',k)]] <- list(exact=exact,estimate=estimate,mcse=se,error=errors,
    rhat=apply(stats,3,posterior::rhat),ess=apply(stats,3,posterior::ess_bulk))
  cat('IW dimension',k,'maximum error/MCSE',max(errors),'maximum R-hat',max(results[[length(results)]]$rhat),'\n')
  flush.console()
}
saveRDS(list(density_error=density_error,jacobian_error=jacobian_error,prior=results),args[2])
cat('IW target and Jacobian checks passed:',density_error,jacobian_error,'\n')
