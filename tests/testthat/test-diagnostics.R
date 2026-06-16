test_that("log_lik computes Gaussian pointwise log likelihood", {
  set.seed(301)
  dat <- data.frame(y = rnorm(6), x = rnorm(6))

  fit <- stLMM(
    y ~ x,
    data = dat,
    starting = list(tau_sq = 1),
    priors = list(tau_sq = ig(2, 1)),
    n_samples = 6,
    verbose = FALSE
  )

  sub <- list(start = 2, thin = 2)
  ll <- log_lik(fit, sub_sample = sub)
  eta <- fitted(fit, summary = FALSE, sub_sample = sub, scale = "link")
  draw_index <- attr(eta, "draw_index")
  expected <- matrix(
    stats::dnorm(
      rep(dat$y, each = nrow(eta)),
      mean = as.vector(eta),
      sd = rep(sqrt(fit$tau_sq_samples[draw_index]), nrow(dat)),
      log = TRUE
    ),
    nrow = nrow(eta)
  )

  expect_equal(dim(ll), c(3L, 6L))
  expect_equal(attr(ll, "draw_index"), draw_index)
  expect_equal(as.vector(ll), as.vector(expected))
})

test_that("log_lik excludes missing responses and works after recovery", {
  set.seed(302)
  dat <- data.frame(
    y = c(0.2, NA, -0.1, 0.4, NA),
    time = seq_len(5)
  )

  fit <- stLMM(
    y ~ 0 + ar1(time),
    data = dat,
    starting = list(tau_sq = 1, ar1_1 = c(sigma_sq = 1, phi = 0.4)),
    tuning = list(tau_sq = 0, ar1_1 = c(sigma_sq = 0, phi = 0)),
    priors = list(
      tau_sq = ig(2, 1),
      ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.9, 0.9))
    ),
    n_samples = 6,
    verbose = FALSE
  )

  expect_error(log_lik(fit), "requires saved or recovered latent process samples")

  rec <- recover(fit, sub_sample = list(start = 2, thin = 2))
  ll <- log_lik(rec)
  eta <- fitted(rec, summary = FALSE, scale = "link")
  obs <- rec$backend$observed_index
  expected <- matrix(
    stats::dnorm(
      rep(dat$y[obs], each = nrow(eta)),
      mean = as.vector(eta[, obs, drop = FALSE]),
      sd = rep(sqrt(rec$tau_sq_samples[rec$recover_iter]), length(obs)),
      log = TRUE
    ),
    nrow = nrow(eta)
  )

  expect_equal(dim(ll), c(length(rec$recover_iter), length(obs)))
  expect_equal(attr(ll, "observed_index"), as.integer(obs))
  expect_equal(as.vector(ll), as.vector(expected))
})

test_that("log_lik computes binomial and negative-binomial likelihoods", {
  set.seed(303)
  dat_bin <- data.frame(y = c(0L, 2L, 1L, 3L), x = c(-1, 0, 1, 2),
                        trials = c(1L, 3L, 2L, 4L))
  fit_bin <- stLMM(
    y ~ x,
    data = dat_bin,
    family = "binomial",
    trials = "trials",
    n_samples = 5,
    verbose = FALSE
  )
  eta_bin <- fitted(fit_bin, summary = FALSE, scale = "link")
  ll_bin <- log_lik(fit_bin)
  expected_bin <- matrix(
    stats::dbinom(
      rep(dat_bin$y, each = nrow(eta_bin)),
      size = rep(dat_bin$trials, each = nrow(eta_bin)),
      prob = stats::plogis(as.vector(eta_bin)),
      log = TRUE
    ),
    nrow = nrow(eta_bin)
  )

  size <- 4
  dat_nb <- data.frame(y = c(0L, 1L, 4L, 2L), x = c(-0.5, 0, 0.5, 1))
  fit_nb <- stLMM(
    y ~ x,
    data = dat_nb,
    family = "negative_binomial",
    size = size,
    n_samples = 5,
    verbose = FALSE
  )
  eta_nb <- fitted(fit_nb, summary = FALSE, scale = "link")
  ll_nb <- log_lik(fit_nb)
  expected_nb <- matrix(
    stats::dnbinom(
      rep(dat_nb$y, each = nrow(eta_nb)),
      size = size,
      mu = exp(as.vector(eta_nb)),
      log = TRUE
    ),
    nrow = nrow(eta_nb)
  )

  expect_equal(as.vector(ll_bin), as.vector(expected_bin))
  expect_equal(as.vector(ll_nb), as.vector(expected_nb))
})

test_that("log_lik uses saved Polya-Gamma process draws", {
  set.seed(304)
  dat <- data.frame(
    y = c(0L, 1L, 1L, 0L, 1L, 0L),
    x = stats::rnorm(6),
    time = seq_len(6)
  )

  fit <- stLMM(
    y ~ x + ar1(time),
    data = dat,
    family = "binomial",
    n_samples = 6,
    warmup = FALSE,
    save_process = list(start = 2, thin = 2),
    starting = list(ar1_1 = c(sigma_sq = 0.5, phi = 0.2)),
    tuning = list(ar1_1 = c(sigma_sq = 0, phi = 0)),
    priors = list(ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.9, 0.9))),
    verbose = FALSE
  )

  ll <- log_lik(fit)
  eta <- fitted(fit, summary = FALSE, scale = "link")
  expected <- matrix(
    stats::dbinom(
      rep(dat$y, each = nrow(eta)),
      size = 1L,
      prob = stats::plogis(as.vector(eta)),
      log = TRUE
    ),
    nrow = nrow(eta)
  )

  expect_equal(attr(ll, "draw_index"), as.integer(c(2, 4, 6)))
  expect_equal(dim(ll), c(3L, nrow(dat)))
  expect_equal(as.vector(ll), as.vector(expected))
})

test_that("log_lik combines multi-chain likelihood draws", {
  set.seed(305)
  dat <- data.frame(y = stats::rnorm(4), x = stats::rnorm(4))
  fit <- stLMM(
    y ~ x,
    data = dat,
    starting = list(tau_sq = 1),
    priors = list(tau_sq = ig(2, 1)),
    chains = 2,
    n_samples = 4,
    verbose = FALSE
  )

  ll <- log_lik(fit, sub_sample = list(start = 2, thin = 2))

  expect_equal(dim(ll), c(4L, nrow(dat)))
  expect_equal(attr(ll, "chain"), c(1L, 1L, 2L, 2L))
  expect_equal(attr(ll, "draw_index"), c(2L, 4L, 2L, 4L))
  expect_equal(attr(ll, "observed_index"), seq_len(nrow(dat)))
})

test_that("waic delegates to loo when available", {
  testthat::skip_if_not_installed("loo")

  set.seed(306)
  dat <- data.frame(y = rnorm(5), x = rnorm(5))
  fit <- stLMM(
    y ~ x,
    data = dat,
    starting = list(tau_sq = 1),
    priors = list(tau_sq = ig(2, 1)),
    n_samples = 5,
    verbose = FALSE
  )

  got <- suppressWarnings(waic(fit))
  expected <- suppressWarnings(loo::waic(log_lik(fit)))

  expect_s3_class(got, "waic")
  expect_equal(got$estimates, expected$estimates)
})
