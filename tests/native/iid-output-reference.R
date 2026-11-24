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
  stopifnot(is.null(fit$alpha_samples))
  fit <- recover(fit)
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
# Vary every covariance draw, including parent-specific blocks and grouped
# residual variances. Conditional whitening tests restoration and joint pairing,
# independently of how well the short parameter chain mixes.
for(case in 1:4){
  kind <- c('independent', 'crossed', 'by', 'collinear')[case]
  structure <- c('none', 'ar1', 'nngp', 'gp')[case]
  suffix <- switch(structure, none = '', ar1 = '+ ar1(time)',
    nngp = '+ nngp(lon, lat, m=3, cov_model="exp")', gp = '+ gp(lon, lat, cov_model="exp")')
  grouped <- case %% 2L == 0L
  form <- as.formula(paste('y ~ x + offset(off) +', forms[[kind]], suffix,
    if(grouped) '+ resid(model="group", group=parent)' else ''))
  pr <- get_priors(kind)
  pr$resid <- list(tau_sq=ig(3,1))
  if(structure != 'none') pr[[paste0(structure, '_1')]] <-
    c(list(sigma_sq=ig(3, .6)), if(structure=='ar1') list(phi=uniform(.1,.8)) else list(phi=uniform(.5,3)))
  set.seed(87400 + case)
  f <- suppressWarnings(stLMM(form, d, priors=pr, collapse_iid=TRUE,
    n_samples=1500, warmup=FALSE, verbose=FALSE))
  r <- recover(f, sub_sample=list(start=301, thin=2))
  idx <- r$recover_iter
  latent <- cbind(r$w_samples_stacked, r$alpha_samples)
  z <- matrix(NA_real_, nrow(latent), ncol(latent))
  for(j in seq_along(idx)){
    be <- f$backend; row <- idx[j]; scalar_col <- cov_col <- 0L
    for(b in seq_along(be$iid_covariance)){
      cb <- be$iid_covariance[[b]]
      if(is.null(cb)){
        scalar_col <- scalar_col + 1L
        be$sigma_sq_re_starting[b] <- f$iid_sigma_sq_samples[row, scalar_col]
      } else {
        k <- cb$dim; columns <- seq_len(k*(k+1L)/2L) + cov_col
        S <- matrix(0,k,k); S[lower.tri(S,diag=TRUE)] <- f$iid_covariance_samples[row,columns]
        S <- S+t(S)-diag(diag(S),k)
        be$iid_covariance[[b]]$starting <- S
        cov_col <- max(columns)
      }
    }
    for(t in seq_along(be$process_terms_obs)){
      be$process_terms_obs[[t]]$sigma_sq_starting <- f$sigma_sq_samples[row,t]
      be$process_terms_obs[[t]]$theta_starting <- f$theta_samples[row,t]
    }
    if(grouped){
      groups <- match(d$parent[be$observed_index], be$residual_model$groups)
      be$residual_model$type <- 'fixed_variance'
      be$residual_model$obs_precision_obs <- 1/f$residual_variance_samples[row,groups]
    } else be$tau_sq_starting <- f$tau_sq_samples[row]
    ref <- reference(be)
    M <- solve(ref$K)+crossprod(ref$H,ref$H*ref$W)
    rhs <- crossprod(ref$H,ref$W*(be$y_model_obs-drop(be$X_obs%*%f$beta_samples[row,])))
    mu <- drop(solve(M,rhs))
    z[j,] <- drop(chol(M)%*%(latent[j,]-mu))
  }
  mean_error <- max(abs(colMeans(z)))*sqrt(nrow(z))
  covariance_error <- abs(cov(z)-diag(ncol(z)))*sqrt(nrow(z)-1)
  diag(covariance_error) <- diag(covariance_error)/sqrt(2)
  stopifnot(mean_error < 6.5, max(covariance_error) < 6.5)
  results[[paste('varying',kind,structure,sep='/')]] <-
    c(mean_error=mean_error,covariance_error=max(covariance_error))
}
saveRDS(results, args[2])
print(apply(do.call(rbind, results), 2, max))
cat(length(results), 'public fitting/recovery Gaussian posterior checks passed\n')
