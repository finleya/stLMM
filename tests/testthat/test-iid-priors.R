iid_prior_data <- function(){
  set.seed(89210)
  data.frame(y = rnorm(24), x = rnorm(24), group = rep(letters[1:6], each = 4),
             parent = rep(c('A', 'B'), each = 12), time = 1:24,
             off = rep(c(-.2, .3), 12))
}

test_that('scalar IID priors retain their family and scale in both strategies', {
  d <- iid_prior_data()
  priors <- list(ig(3, 1), log_normal(-.4, .7), gamma_dist(2, 3), half_normal(.8), half_t(4, .8))
  for(collapse in c(FALSE, TRUE)) for(i in seq_along(priors)){
    pr <- priors[[i]]
    # Term controls are normalized without changing coefficient names.
    f <- stLMM(y ~ x + iid(group), d, priors = list(iid_1 = list(SIGMA_SQ = pr)),
      collapse_iid = collapse, n_samples = 20, warmup = FALSE, verbose = FALSE)
    expect_equal(unname(f$backend$sigma_sq_re_prior[1, 1]), c(0, 2, 3, 4, 5)[i])
    term <- f$term_description$random_effects$terms[[1]]
    expect_identical(term$sampler$status, if(collapse) 'collapsed Metropolis' else
      if(i == 1) 'gibbs' else 'conditional Metropolis')
    expect_equal(term$sigma_sq$prior_hyperparameters, pr$parameters)
    if(i > 1) expect_identical(term$sigma_sq$prior_scale,
      if(i >= 4) 'standard deviation' else 'variance')
    expect_true(all(is.finite(f$iid_sigma_sq_samples) & f$iid_sigma_sq_samples > 0))
  }
})

test_that('mixed IG, nonconjugate and IW terms share the existing blocking machinery', {
  d <- iid_prior_data()
  form <- y ~ x + iid(group, ~ 1 + x, cov_model = 'independent', by = parent) +
    iid(parent, ~ 1 + x) + ar1(time)
  pr <- list(iid_1 = list(sigma_sq = list('(Intercept)' = ig(3, 1), x = half_t(4, .8))),
    iid_2 = list(covariance = iw(4, diag(2))))
  for(collapse in c(FALSE, TRUE)) for(blocking in c('joint', 'by_term', 'residual_process',
      'variance_theta', 'process_variance', 'scalar')){
    f <- stLMM(form, d, priors = pr, starting = list(ar1_1 = list(sigma_sq = fixed(.3), phi = fixed(.4))),
      tuning = list(iid_1 = list(sigma_sq = c('(Intercept)' = .1, x = .3))),
      collapse_iid = collapse, metropolis = blocking, n_samples = 20, verbose = FALSE,
      warmup = list(batch_length = 5, min_batches = 2, max_batches = 2))
    expect_equal(f$adaptive_metropolis$dimension, if(collapse) 8L else 3L)
    expect_equal(f$adaptive_metropolis$warmup$n_batches, 2L)
    status <- vapply(f$term_description$random_effects$terms, function(t) t$sampler$status, character(1))
    expect_equal(unname(status), if(collapse) rep('collapsed Metropolis', 5) else
      c('gibbs', 'gibbs', 'conditional Metropolis', 'conditional Metropolis', 'gibbs'))
    expect_equal(unname(f$backend$sigma_sq_re_prior[, 1]), c(0, 0, 5, 5, 0))
    expect_true(all(is.finite(f$iid_covariance_samples)))
  }
})

test_that('zero tuning fixes only Metropolis-updated IID coefficients across chains', {
  d <- iid_prior_data()
  form <- y ~ x + iid(group, ~ 1 + x, cov_model = 'independent', by = parent)
  for(collapse in c(FALSE, TRUE)){
    f <- stLMM(form, d,
      priors = list(iid_1 = list(sigma_sq = list('(Intercept)' = ig(3, 1), x = log_normal(-.2, .5)))),
      starting = list(iid_1 = list(sigma_sq = .7)),
      tuning = list(iid_1 = list(sigma_sq = 0)), collapse_iid = collapse,
      n_samples = 15, chains = 2, chain_control = list(seed = 89211), warmup = FALSE, verbose = FALSE)
    for(ch in f$chains){
      expect_true(all(ch$iid_sigma_sq_samples[, 3:4] == .7))
      expect_identical(ch$backend$sigma_sq_re_starting[3:4], c(.7, .7))
      if(collapse) expect_true(all(ch$iid_sigma_sq_samples == .7)) else
        expect_true(all(apply(ch$iid_sigma_sq_samples[, 1:2], 2, sd) > 0))
    }
  }
  f <- stLMM(y ~ 0 + iid(group), d, priors = list(iid_1 = list(sigma_sq = half_normal(1))),
    starting = list(tau_sq = fixed(.3)), n_samples = 12, verbose = FALSE,
    warmup = list(batch_length = 5, min_batches = 2, max_batches = 2))
  expect_equal(f$adaptive_metropolis$dimension, 1L)
  expect_equal(f$adaptive_metropolis$warmup$n_batches, 2L)
})

test_that('both independent IID syntaxes share nonconjugate updates', {
  d <- iid_prior_data()
  for(collapse in c(FALSE, TRUE)){
    set.seed(89212)
    a <- stLMM(y ~ x + iid(group) + x:iid(group), d,
      priors = list(iid_1 = list(sigma_sq = half_normal(.8)), iid_2 = list(sigma_sq = log_normal(-.4, .7))),
      n_samples = 30, collapse_iid = collapse, warmup = FALSE, verbose = FALSE)
    rng <- .Random.seed
    set.seed(89212)
    b <- stLMM(y ~ x + iid(group, ~ 1 + x, cov_model = 'independent'), d,
      priors = list(iid_1 = list(sigma_sq = list('(Intercept)' = half_normal(.8), x = log_normal(-.4, .7)))),
      n_samples = 30, collapse_iid = collapse, warmup = FALSE, verbose = FALSE)
    expect_identical(.Random.seed, rng)
    expect_equal(a$beta_samples, b$beta_samples)
    expect_equal(a$alpha_samples, b$alpha_samples, ignore_attr = TRUE)
    expect_equal(a$iid_sigma_sq_samples, b$iid_sigma_sq_samples, ignore_attr = TRUE)
  }
})

test_that('IID prior validation rejects unsupported families and malformed controls', {
  d <- iid_prior_data()
  fit <- function(pr) stLMM(y ~ iid(group), d, priors = list(iid_1 = list(sigma_sq = pr)),
    n_samples = 2, warmup = FALSE, verbose = FALSE)
  for(pr in list(uniform(0, 1), beta_dist(2, 2), flat(), normal()))
    expect_error(fit(pr), 'unsupported variance prior')
  for(pr in list(log_normal(0, 1, support = c(.1, 2)), gamma_dist(2, 1, support = c(.1, 2))))
    expect_error(fit(pr), 'do not currently accept finite support')
  for(pr in list(log_normal(0, 0), gamma_dist(-1, 1), half_normal(-1), half_t(0, 1)))
    expect_error(fit(pr), 'requires')
})

test_that('legacy IG backends without generic prior metadata remain readable', {
  d <- iid_prior_data()
  f <- stLMM(y ~ x + iid(group), d, priors = list(iid_1 = list(sigma_sq = ig(3, 1))),
    n_samples = 12, warmup = FALSE, verbose = FALSE)
  for(collapse in c(FALSE, TRUE)){
    b <- f$backend; b$collapse_iid <- collapse
    set.seed(89213); a <- .Call('stLMM_collapsed_sampler', b, PACKAGE = 'stLMM'); rng <- .Random.seed
    b$sigma_sq_re_prior <- NULL
    set.seed(89213); z <- .Call('stLMM_collapsed_sampler', b, PACKAGE = 'stLMM')
    expect_identical(.Random.seed, rng)
    fields <- grep('samples$|adaptive_metropolis|term_description', names(a), value = TRUE)
    expect_identical(a[fields], z[fields])
  }
})

test_that('new IID priors preserve paired prediction with structured effects and missing responses', {
  for(family in c('gaussian', 'binomial', 'probit', 'negative_binomial', 'estimated_size')){
    d <- iid_prior_data()
    if(family != 'gaussian') d$y <- as.numeric(d$y > 0)
    d$y[c(2, 9)] <- NA_real_
    d$lon <- d$time / nrow(d); d$lat <- d$x
    for(collapse in c(FALSE, TRUE)){
      pr <- list(iid_1 = list(sigma_sq = list('(Intercept)' = half_t(4, .8), x = ig(3, 1))),
        iid_2 = list(covariance = iw(4, diag(2))))
      if(family == 'estimated_size') pr$size <- gamma_dist(2, 1)
      args <- list(formula = y ~ x + offset(off) +
        iid(group, ~ 1 + x, cov_model = 'independent', by = parent) +
        iid(parent, ~ 1 + x) + nngp(lon, lat, m = 5, cov_model = 'exp'), data = d,
        priors = pr, starting = list(nngp_1 = list(sigma_sq = fixed(.3), phi = fixed(2))),
        family = if(family == 'estimated_size') 'negative_binomial' else family,
        size = if(family == 'negative_binomial') 2 else NULL, collapse_iid = collapse,
        n_samples = 40, warmup = FALSE, verbose = FALSE,
        save_process = if(family != 'gaussian') list(start = 2, thin = 2) else NULL)
      if(family == 'estimated_size') expect_warning(f <- do.call(stLMM, args), 'weakly identified')
      else f <- do.call(stLMM, args)
      ss <- list(start = 4, thin = 4)
      f <- recover(f, sub_sample = ss)
      eta <- fitted(f, summary = FALSE, scale = 'link', sub_sample = ss)
      idx <- attr(eta, 'draw_index')
      expected <- f$beta_samples[idx, , drop = FALSE] %*% t(f$backend$X) +
        f$alpha_samples[idx, , drop = FALSE] %*% t(as.matrix(f$backend$Z)) + f$w_samples$nngp_1
      expected <- sweep(expected, 2, d$off, '+')
      expect_equal(eta, expected, ignore_attr = TRUE)
      p <- predict(f, d[, setdiff(names(d), 'y')], scale = 'link', sub_sample = ss)
      expect_equal(p$draw_index, idx)
      expect_equal(p$mu_samples, expected, ignore_attr = TRUE)
      obs <- which(!is.na(d$y))
      ll <- t(vapply(seq_along(idx), function(i) switch(family,
        gaussian = dnorm(d$y[obs], expected[i, obs], sqrt(f$tau_sq_samples[idx[i]]), log = TRUE),
        binomial = dbinom(d$y[obs], 1, plogis(expected[i, obs]), log = TRUE),
        probit = dbinom(d$y[obs], 1, pnorm(expected[i, obs]), log = TRUE),
        dnbinom(d$y[obs], size = if(family == 'estimated_size') f$size_samples[idx[i]] else 2,
          mu = exp(expected[i, obs]), log = TRUE)), numeric(length(obs))))
      expect_equal(log_lik(f, sub_sample = ss), ll, ignore_attr = TRUE)
    }
  }
})
