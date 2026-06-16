test_that("as_samples returns posterior draw data frames for stLMM fits", {
  set.seed(6001)
  dat <- data.frame(y = rnorm(8), x = rnorm(8))

  fit <- stLMM(
    y ~ x,
    data = dat,
    starting = list(tau_sq = 1),
    priors = list(tau_sq = ig(2, 1)),
    n_samples = 6,
    verbose = FALSE
  )

  draws <- as_samples(fit, burn = 1, thin = 2)

  expect_s3_class(draws, "data.frame")
  expect_equal(draws$.chain, rep(1L, 3))
  expect_equal(draws$.iteration, c(2L, 4L, 6L))
  expect_true(all(c("(Intercept)", "x", "tau_sq") %in% names(draws)))
  expect_equal(draws[["(Intercept)"]], fit$beta_samples[c(2, 4, 6), "(Intercept)"])
  expect_equal(draws$tau_sq, fit$tau_sq_samples[c(2, 4, 6)])

  no_meta <- as_samples(fit, burn = 1, thin = 2, metadata = FALSE)
  expect_false(".chain" %in% names(no_meta))
  expect_false(".iteration" %in% names(no_meta))
  expect_equal(no_meta$x, fit$beta_samples[c(2, 4, 6), "x"])
})

test_that("as_samples combines multi-chain fit draws and preserves chain ids", {
  set.seed(6002)
  dat <- data.frame(y = rnorm(10), x = rnorm(10))

  fit <- stLMM(
    y ~ x,
    data = dat,
    n_samples = 6,
    chains = 2,
    chain_control = list(seed = 6002),
    verbose = FALSE
  )

  draws <- as_samples(fit, burn = 1, thin = 3)

  expect_equal(draws$.chain, c(1L, 1L, 2L, 2L))
  expect_equal(draws$.iteration, c(2L, 5L, 2L, 5L))
  expect_equal(draws[draws$.chain == 1L, "tau_sq"], fit$chains[[1]]$tau_sq_samples[c(2, 5)])
  expect_equal(draws[draws$.chain == 2L, "tau_sq"], fit$chains[[2]]$tau_sq_samples[c(2, 5)])

  by_chain <- as_samples(fit, burn = 1, thin = 3, combine_chains = FALSE)
  expect_type(by_chain, "list")
  expect_length(by_chain, 2)
  expect_equal(by_chain[[2]]$.chain, c(2L, 2L))
})

test_that("as_samples aligns recovered parameters with recovered process draws", {
  set.seed(6003)
  dat <- data.frame(
    y = rnorm(4),
    time = seq_len(4)
  )

  fit <- stLMM(
    y ~ ar1(time),
    data = dat,
    starting = list(tau_sq = 1, ar1_1 = c(sigma_sq = 1, phi = 0.4)),
    tuning = list(tau_sq = 0, ar1_1 = c(sigma_sq = 0, phi = 0)),
    priors = list(
      tau_sq = ig(2, 1),
      ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.9, 0.9))
    ),
    n_samples = 8,
    verbose = FALSE
  )

  rec <- recover(fit, sub_sample = list(start = 2, thin = 2))
  draws <- as_samples(rec, include_w = TRUE)
  m <- as.matrix(as_mcmc(rec, include_w = TRUE))

  expect_equal(draws$.iteration, rec$recover_iter)
  expect_equal(draws[["(Intercept)"]], rec$beta_samples[rec$recover_iter, "(Intercept)"])
  expect_equal(draws$tau_sq, rec$tau_sq_samples[rec$recover_iter])
  expect_equal(draws$w_ar1_1_1, rec$w_samples$ar1_1[, 1])
  expect_equal(nrow(m), length(rec$recover_iter))
  expect_equal(m[, "(Intercept)"], rec$beta_samples[rec$recover_iter, "(Intercept)"])
  expect_equal(m[, "w_ar1_1_1"], rec$w_samples$ar1_1[, 1])

  thinned <- as_samples(rec, burn = 3, thin = 2, include_w = TRUE)
  m_thinned <- as.matrix(as_mcmc(rec, burn = 3, thin = 2, include_w = TRUE))
  expect_equal(thinned$.iteration, c(4L, 8L))
  expect_equal(thinned$w_ar1_1_1, rec$w_samples$ar1_1[c(2, 4), 1])
  expect_equal(nrow(m_thinned), 2L)
  expect_equal(m_thinned[, "w_ar1_1_1"], rec$w_samples$ar1_1[c(2, 4), 1])
})

test_that("as_samples handles prediction and fitted sample outputs", {
  set.seed(6004)
  dat <- data.frame(y = rnorm(8), x = seq(-1, 1, length.out = 8))

  fit <- stLMM(
    y ~ x,
    data = dat,
    starting = list(tau_sq = 1),
    priors = list(tau_sq = ig(2, 1)),
    n_samples = 6,
    verbose = FALSE
  )

  pred <- predict(
    fit,
    newdata = data.frame(x = c(-0.5, 0.5)),
    sub_sample = list(start = 2, thin = 2),
    y_samples = TRUE
  )
  pred_draws <- as_samples(pred, sample = "all")

  expect_equal(pred_draws$.iteration, pred$draw_index)
  expect_equal(pred_draws$mu_1, pred$mu_samples[, 1])
  expect_equal(pred_draws$y_2, pred$y_samples[, 2])

  pred_mu_only <- predict(
    fit,
    newdata = data.frame(x = c(-0.5, 0.5)),
    sub_sample = list(start = 2, thin = 2)
  )
  pred_mu_only_draws <- as_samples(pred_mu_only, sample = "all", metadata = FALSE)
  expect_named(pred_mu_only_draws, c("mu_1", "mu_2"))
  expect_error(as_samples(pred_mu_only, sample = "y"), "y samples are not available")

  mu <- fitted(fit, summary = FALSE, sub_sample = list(start = 2, thin = 2))
  fitted_draws <- as_samples(mu)

  expect_equal(fitted_draws$.iteration, attr(mu, "draw_index"))
  expect_equal(fitted_draws[[3]], mu[, 1])
})

test_that("as_samples handles fitted sample lists from recovered chains", {
  set.seed(6005)
  dat <- data.frame(
    y = rnorm(4),
    time = seq_len(4)
  )

  fit <- stLMM(
    y ~ 0 + ar1(time),
    data = dat,
    n_samples = 6,
    chains = 2,
    priors = list(
      tau_sq = ig(2, 1),
      ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.8, 0.8))
    ),
    verbose = FALSE
  )

  rec <- recover(fit, sub_sample = list(start = 2, thin = 2))
  mu <- fitted(rec, summary = FALSE)
  draws <- as_samples(mu)

  expect_s3_class(mu, "stLMM_fitted_chains")
  expect_equal(draws$.chain, c(rep(1L, 3), rep(2L, 3)))
  expect_equal(draws$.iteration, rep(c(2L, 4L, 6L), 2))
  expect_equal(draws[draws$.chain == 2L, 3], mu[[2]][, 1])
})
