test_that("fixed-effect sampler uses observed rows when response is missing", {
  set.seed(101)
  dat <- data.frame(y = rnorm(8), x = rnorm(8))
  dat_mis <- dat
  dat_mis$y[c(3, 7)] <- NA_real_
  dat_obs <- dat[-c(3, 7), , drop = FALSE]

  args <- list(
    formula = y ~ x,
    starting = list(tau_sq = 1),
    tuning = list(tau_sq = 0),
    priors = list(tau_sq = ig(2, 1)),
    n_samples = 6,
    verbose = FALSE
  )

  set.seed(102)
  fit_mis <- do.call(stLMM, c(args, list(data = dat_mis)))
  set.seed(102)
  fit_obs <- do.call(stLMM, c(args, list(data = dat_obs)))

  expect_equal(fit_mis$backend$n, 8L)
  expect_equal(fit_mis$backend$n_obs, 6L)
  expect_equal(fit_mis$backend$n_missing_response, 2L)
  expect_equal(fit_mis$backend$observed_index, as.integer(c(1, 2, 4, 5, 6, 8)))
  expect_equal(unname(fit_mis$beta_samples), unname(fit_obs$beta_samples))
  expect_equal(unname(fit_mis$tau_sq_samples), unname(fit_obs$tau_sq_samples))
})

test_that("iid random-effect sampler uses observed rows when response is missing", {
  set.seed(103)
  dat <- data.frame(
    y = rnorm(10),
    x = rnorm(10),
    group = factor(rep(letters[1:5], each = 2))
  )
  dat_mis <- dat
  dat_mis$y[c(2, 9)] <- NA_real_
  dat_obs <- droplevels(dat[-c(2, 9), , drop = FALSE])

  args <- list(
    formula = y ~ x + iid(group),
     starting = list(tau_sq = 1, iid_1 = list(sigma_sq = 1)),
    tuning = list(tau_sq = 0),
    priors = list(tau_sq = ig(2, 1), iid_1 = list(sigma_sq = ig(2, 1))),
    n_samples = 6,
    verbose = FALSE
  )

  set.seed(104)
  fit_mis <- do.call(stLMM, c(args, list(data = dat_mis)))
  set.seed(104)
  fit_obs <- do.call(stLMM, c(args, list(data = dat_obs)))

  expect_equal(fit_mis$backend$n_obs, 8L)
  expect_equal(unname(fit_mis$beta_samples), unname(fit_obs$beta_samples))
  expect_equal(unname(fit_mis$alpha_samples), unname(fit_obs$alpha_samples))
  expect_equal(unname(fit_mis$iid_sigma_sq_samples), unname(fit_obs$iid_sigma_sq_samples))
})

test_that("process recovery keeps full support with missing response rows", {
  set.seed(105)
  dat <- data.frame(
    y = rnorm(6),
    time = seq_len(6)
  )
  dat$y[4] <- NA_real_

  fit <- stLMM(
    y ~ 0 + ar1(time),
    data = dat,
    starting = list(tau_sq = 1, ar1_1 = c(sigma_sq = 1, phi = 0.3)),
    tuning = list(tau_sq = 0, ar1_1 = c(sigma_sq = 0, phi = 0)),
    priors = list(
      tau_sq = ig(2, 1),
      ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.9, 0.9))
    ),
    n_samples = 6,
    verbose = FALSE
  )

  rec <- recover(fit, sub_sample = list(start = 2, thin = 2))
  mu <- fitted(rec, summary = FALSE)

  expect_equal(fit$backend$n, 6L)
  expect_equal(fit$backend$n_obs, 5L)
  expect_equal(length(fitted(rec)), 6L)
  expect_equal(ncol(mu), 6L)
  expect_equal(ncol(rec$w_samples$ar1_1), 6L)
})

test_that("missing predictors still error and summary reports missing response count", {
  dat <- data.frame(y = rnorm(5), x = rnorm(5))
  dat$y[2] <- NA_real_

  fit <- stLMM(
    y ~ x,
    data = dat,
    starting = list(tau_sq = 1),
    tuning = list(tau_sq = 0),
    priors = list(tau_sq = ig(2, 1)),
    n_samples = 4,
    verbose = FALSE
  )

  sm <- summary(fit)
  expect_equal(sm$n, 5L)
  expect_equal(sm$n_obs, 4L)
  expect_equal(sm$n_missing_response, 1L)

  dat_bad <- dat
  dat_bad$x[1] <- NA_real_
  expect_error(
    stLMM(
      y ~ x,
      data = dat_bad,
      starting = list(tau_sq = 1),
      priors = list(tau_sq = ig(2, 1)),
      n_samples = 4,
      verbose = FALSE
    ),
    "missing data in model predictors"
  )
})
