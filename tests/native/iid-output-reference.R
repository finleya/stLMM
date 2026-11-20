# Public fitting/recovery draws against the joint Gaussian posterior, with
# covariance parameters held fixed. Run with library and output-RDS arguments.
Sys.setenv(STLMM_REFERENCE_ONLY = '1')
source('tests/native/latent-reference.R')
Sys.unsetenv('STLMM_REFERENCE_ONLY')
d$off <- seq(-.3, .4, length.out = nrow(d))
results <- list(); case <- 0L
for(kind in names(forms)) for(weighted in c(FALSE, TRUE)) {
  case <- case + 1L
  structure <- c('none', 'ar1', 'nngp', 'gp')[(case-1L) %% 4L + 1L]
  suffix <- switch(structure, none = '', ar1 = '+ ar1(time)',
    nngp = '+ nngp(lon, lat, m=3, cov_model="exp")', gp = '+ gp(lon, lat, cov_model="exp")')
  form <- as.formula(paste('y ~ x + offset(off) +', forms[[kind]], suffix,
    if(weighted) '+ resid(model="fixed", variance=v)' else ''))
  pr <- get_priors(kind); pr$beta <- normal(c(.4, -.3), 1.3)
  tune <- lapply(pr[grepl('^iid_', names(pr))], function(b) setNames(list(0), names(b)))
  st <- switch(structure, none = list(),
    ar1 = list(ar1_1 = list(sigma_sq=fixed(.4), phi=fixed(.6))),
    nngp = list(nngp_1 = list(sigma_sq=fixed(.4), phi=fixed(2))),
    gp = list(gp_1 = list(sigma_sq=fixed(.4), phi=fixed(2))))
  if(!weighted) st$resid <- list(tau_sq = fixed(.7))
  set.seed(87300 + case)
  fit <- suppressWarnings(stLMM(form, d, priors=pr, starting=st, tuning=tune,
    collapse_iid=TRUE, n_samples=5000, warmup=FALSE, verbose=FALSE))
  be <- fit$backend; ref <- reference(be)
  # Independent dense reference includes the proper fixed-effect prior and offset.
  design <- cbind(be$X_obs, ref$H)
  J <- crossprod(design, design * ref$W)
  ii <- seq_len(be$p); jj <- be$p + seq_len(ncol(ref$H))
  J[cbind(ii,ii)] <- J[cbind(ii,ii)] + be$beta_prior_precision
  J[jj,jj] <- J[jj,jj] + solve(ref$K)
  rhs <- drop(crossprod(design, be$y_model_obs * ref$W))
  rhs[ii] <- rhs[ii] + be$beta_prior_precision * be$beta_prior_mean
  mu <- drop(solve(J, rhs))
  original <- fit$alpha_samples
  if(structure != 'none') fit <- recover(fit)
  stopifnot(identical(fit$alpha_samples, original))
  draws <- cbind(fit$beta_samples, fit$w_samples_stacked, fit$alpha_samples)
  # Joint collapsing with fixed covariances gives independent posterior draws.
  z <- sweep(draws, 2, mu) %*% t(chol(J))
  mean_error <- max(abs(colMeans(z))) * sqrt(nrow(z))
  covariance_error <- abs(cov(z)-diag(ncol(z))) * sqrt(nrow(z)-1)
  diag(covariance_error) <- diag(covariance_error)/sqrt(2)
  stopifnot(mean_error < 6.5, max(covariance_error) < 6.5)
  obs <- be$observed_index
  eta <- fitted(fit, summary=FALSE, scale='link')
  expected <- sweep(draws %*% t(design), 2, be$offset_obs, '+')
  stopifnot(max(abs(eta[,obs]-expected)) < 1e-10)
  id <- paste(kind, structure, weighted, sep='/')
  results[[id]] <- c(mean_error=mean_error, covariance_error=max(covariance_error))
}
saveRDS(results, args[2])
print(apply(do.call(rbind, results), 2, max))
cat(length(results), 'public fitting/recovery Gaussian posterior checks passed\n')
