collapsed_iid_data <- function(){
  set.seed(87201)
  data.frame(y = rnorm(24), x = rnorm(24), group = rep(letters[1:6], each = 4),
    parent = rep(c("A", "B"), each = 12), time = 1:24, off = rep(c(-.2, .3), 12))
}

test_that("IID strategy selection is explicit and preserves the default", {
  d <- collapsed_iid_data()
  fit <- function(formula = y ~ x + iid(group), ...){
    set.seed(87202)
    stLMM(formula, d, priors = list(iid_1 = list(sigma_sq = ig(3, 1)), resid = list(tau_sq = ig(3, 1))),
      n_samples = 12, warmup = FALSE, verbose = FALSE, ...)
  }
  a <- fit(); rng <- .Random.seed
  b <- fit(collapse_iid = FALSE)
  expect_identical(a$samples, b$samples)
  expect_identical(.Random.seed, rng)
  expect_false(a$backend$collapse_iid)
  expect_identical(a$backend$collapse_iid_requested, FALSE)
  a <- fit(collapse_iid = NULL); rng <- .Random.seed
  b <- fit(collapse_iid = TRUE)
  expect_identical(a$samples, b$samples)
  expect_identical(.Random.seed, rng)
  expect_true(a$backend$collapse_iid)
  expect_null(a$backend$collapse_iid_requested)
  expect_output(print(b), "IID sampler: collapsed")
  expect_identical(b$term_description$random_effects$terms[[1]]$sampler$status,
                   "collapsed Metropolis")
  for(formula in list(y ~ x, y ~ x + iid(group) + ar1(time))){
    st <- list(ar1_1 = list(sigma_sq = fixed(.4), phi = fixed(.5)))
    a <- fit(formula, collapse_iid = NULL, starting = st); rng <- .Random.seed
    b <- fit(formula, collapse_iid = FALSE, starting = st)
    expect_false(a$backend$collapse_iid)
    expect_identical(a$samples, b$samples)
    expect_identical(.Random.seed, rng)
  }
  # Grouped residual variances do not select the structured-process strategy.
  a <- fit(y ~ x + iid(group) + resid(model = "group", group = parent),
           collapse_iid = NULL)
  expect_true(a$backend$collapse_iid)
  expect_error(fit(y ~ x, collapse_iid = TRUE), "requires at least one iid")
  for(value in list(NA, 1, "auto", logical(), c(TRUE, FALSE), list(TRUE)))
    expect_error(fit(collapse_iid = value), "collapse_iid must be")
})

test_that("collapsed tuning follows coefficient and parent controls", {
  d <- collapsed_iid_data()
  form <- y ~ x + iid(group, ~ 1 + x, cov_model = "independent", by = parent) +
    iid(parent, ~ 1 + x)
  make <- function(tuning, ...)
    stLMM(form, d, priors = list(iid_1 = list(sigma_sq = ig(3, 1)),
      iid_2 = list(covariance = iw(4, diag(2)))), tuning = tuning,
      n_samples = 30, chains = 2, chain_control = list(seed = 87203),
      collapse_iid = TRUE, warmup = FALSE, verbose = FALSE, ...)
  f <- make(list(iid_1 = list(sigma_sq = c(x = .3, "(Intercept)" = 0)),
                 iid_2 = list(covariance = 0)))
  expect_output(print(f), "IID sampler: collapsed")
  for(ch in f$chains){
    expect_identical(ch$backend$collapse_iid_requested, TRUE)
    expect_true(ch$backend$collapse_iid)
    expect_equal(ch$backend$sigma_sq_re_tuning, c(0, 0, .3, .3, .2))
    expect_equal(ch$iid_sigma_sq_samples[, 1:2], matrix(1, 30, 2), ignore_attr = TRUE)
    expect_equal(ch$iid_covariance_samples, matrix(rep(c(1/7, 0, 1/7), each = 30), 30),
                 ignore_attr = TRUE)
    expect_equal(ch$adaptive_metropolis$dimension, 3L) # residual and two slopes
  }
  for(value in list(-1, Inf, NA_real_, numeric(), c(.1, .2), "0.1", fixed(.1))){
    expect_error(make(list(iid_1 = list(sigma_sq = value))), "sigma_sq|coefficient")
    expect_error(make(list(iid_2 = list(covariance = value))), "covariance")
  }
  expect_error(make(list(iid_2 = list(sigma_sq = .1))), "must contain covariance")
  expect_error(make(list(iid_1 = list(sigma_sq = list(x = .1)))), "coefficient names")
})

test_that("collapsed public fits retain paired predictors, likelihoods, and size draws", {
  for(family in c("gaussian", "binomial", "probit", "negative_binomial", "estimated_size")){
    d <- collapsed_iid_data()
    if(family != "gaussian") d$y <- as.numeric(d$y > 0)
    d$y[c(2, 9)] <- NA_real_
    for(process in c(FALSE, TRUE)){
      formula <- y ~ x + offset(off) + iid(group, ~ 1 + x) + iid(parent)
      if(process) formula <- update(formula, . ~ . + ar1(time))
      pr <- list(iid_1 = list(covariance = iw(4, diag(2))),
                 iid_2 = list(sigma_sq = ig(3, 1)))
      if(family == "estimated_size") pr$size <- gamma_dist(2, 1)
      args <- list(formula = formula, data = d, priors = pr,
        family = if(family == "estimated_size") "negative_binomial" else family,
        size = if(family == "negative_binomial") 2 else NULL,
        starting = list(ar1_1 = list(sigma_sq = fixed(.4), phi = fixed(.5))),
        n_samples = if(family == "gaussian") 320 else 80,
        chains = 2, chain_control = list(seed = 87204),
        warmup = FALSE, verbose = FALSE, collapse_iid = TRUE,
        save_process = if(process && family != "gaussian") list(start = 2, thin = 2) else NULL)
      if(family == "estimated_size")
        expect_warning(f <- do.call(stLMM, args), "weakly identified")
      else
        f <- do.call(stLMM, args)
      if(process || family == "gaussian") f <- recover(f, sub_sample = list(start = 4, thin = 4))
      for(ch in f$chains){
        ss <- list(start = 4, thin = 4)
        eta <- fitted(ch, summary = FALSE, scale = "link", sub_sample = ss)
        idx <- attr(eta, "draw_index")
        expected <- ch$beta_samples[idx, , drop = FALSE] %*% t(ch$backend$X) +
          ch$alpha_samples[if(family == "gaussian") match(idx, ch$recover_iter) else idx, , drop = FALSE] %*% t(as.matrix(ch$backend$Z))
        expected <- sweep(expected, 2, d$off, "+")
        if(process) expected <- expected + ch$w_samples$ar1_1
        expect_equal(eta, expected, ignore_attr = TRUE)
        p <- predict(ch, d[, setdiff(names(d), "y")], scale = "link", sub_sample = ss)
        expect_equal(p$draw_index, idx)
        expect_equal(p$mu_samples, expected, ignore_attr = TRUE)
        obs <- which(!is.na(d$y))
        ll <- t(vapply(seq_along(idx), function(i){
          switch(family,
            gaussian = dnorm(d$y[obs], expected[i, obs], sqrt(ch$tau_sq_samples[idx[i]]), log = TRUE),
            binomial = dbinom(d$y[obs], 1, plogis(expected[i, obs]), log = TRUE),
            probit = dbinom(d$y[obs], 1, pnorm(expected[i, obs]), log = TRUE),
            dnbinom(d$y[obs], size = if(family == "estimated_size") ch$size_samples[idx[i]] else 2,
                    mu = exp(expected[i, obs]), log = TRUE))
        }, numeric(length(obs))))
        expect_equal(log_lik(ch, sub_sample = ss), ll, ignore_attr = TRUE)
        if(family == "estimated_size"){
          set.seed(87205)
          p <- predict(ch, d[, setdiff(names(d), "y")], y_samples = TRUE, sub_sample = ss)
          set.seed(87205)
          ys <- t(vapply(seq_along(idx), function(i)
            rnbinom(nrow(d), size = ch$size_samples[idx[i]], mu = p$mu_samples[i, ]), numeric(nrow(d))))
          expect_equal(p$y_samples, ys, ignore_attr = TRUE)
        }
      }
      if(requireNamespace("posterior", quietly = TRUE)){
        s <- summary(f, type = "combined", terms = c("iid_1", "iid_2"), diagnostic_method = "rank")
        expect_true(all(is.finite(as.matrix(s$diagnostics[, c("rhat", "ess_bulk", "ess_tail")]))))
      }
    }
  }
})

test_that("public IID covariance controls enter every Metropolis blocking scheme", {
  d <- collapsed_iid_data()
  dimensions <- list(joint = 11L, by_term = c(1L, 2L, 1L, 1L, 3L, 3L),
    residual_process = c(1L, 10L), variance_theta = c(10L, 1L),
    process_variance = c(1L, 9L, 1L), scalar = rep(1L, 11))
  for(blocking in names(dimensions)){
    f <- stLMM(y ~ x + iid(group, ~ 1 + x, cov_model = "independent") +
      iid(group, ~ 1 + x, by = parent) + ar1(time), d,
      priors = list(iid_1 = list(sigma_sq = ig(3, 1)),
        iid_2 = list(covariance = iw(4, diag(2))),
        ar1_1 = list(sigma_sq = ig(3, 1), phi = uniform(-.8, .8))),
      tuning = list(iid_1 = list(sigma_sq = c(x = .3, "(Intercept)" = .1)),
        iid_2 = list(covariance = .15)), collapse_iid = TRUE,
      n_samples = 12, metropolis = blocking, verbose = FALSE,
      warmup = list(batch_length = 5, min_batches = 2, max_batches = 2))
    a <- f$adaptive_metropolis
    expect_equal(a$dimension, 11L, info = blocking)
    expect_equal(unname(vapply(a$blocks, `[[`, integer(1), "dimension")),
                 dimensions[[blocking]], info = blocking)
    expect_equal(a$warmup$n_batches, 2L, info = blocking)
    expect_true(all(is.finite(f$iid_covariance_samples)), info = blocking)
  }
})
