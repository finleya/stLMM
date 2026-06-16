library(stLMM)

test_that("beta prior defaults follow likelihood family", {
  dat_g <- data.frame(y = c(1, 2, 3), x = c(-1, 0, 1))
  fit_g <- stLMM(
    y ~ x,
    data = dat_g,
    n_samples = 2,
    starting = list(tau_sq = 1),
    tuning = list(tau_sq = 0),
    priors = list(tau_sq = ig(2, 1)),
    verbose = FALSE
  )
  expect_equal(fit_g$backend$beta_prior_type, 0L)
  expect_equal(fit_g$backend$beta_prior_precision, c(0, 0))

  dat_b <- data.frame(y = c(0, 1, 1), x = c(-1, 0, 1))
  fit_b <- stLMM(
    y ~ x,
    data = dat_b,
    family = "binomial",
    n_samples = 2,
    verbose = FALSE
  )
  expect_equal(fit_b$backend$beta_prior_type, 1L)
  expect_equal(fit_b$backend$beta_prior_mean, c(0, 0))
  expect_equal(fit_b$backend$beta_prior_precision, c(0.01, 0.01))
})

test_that("family accepts documented strings and stats family objects", {
  dat <- data.frame(y = c(0, 1, 1), x = c(-1, 0, 1))

  fit_string <- stLMM(
    y ~ x,
    data = dat,
    family = "binomial",
    n_samples = 2,
    verbose = FALSE
  )

  fit_function <- stLMM(
    y ~ x,
    data = dat,
    family = stats::binomial(),
    n_samples = 2,
    verbose = FALSE
  )

  expect_equal(fit_string$backend$family, "binomial")
  expect_equal(fit_function$backend$family, "binomial")
})

test_that("normal beta prior accepts scalar and vector controls", {
  dat <- data.frame(y = c(1, 2, 3), x = c(-1, 0, 1))

  fit_scalar <- stLMM(
    y ~ x,
    data = dat,
    n_samples = 2,
    starting = list(tau_sq = 1),
    tuning = list(tau_sq = 0),
    priors = list(tau_sq = ig(2, 1), beta = normal(mean = 1, sd = 2)),
    verbose = FALSE
  )
  expect_equal(fit_scalar$backend$beta_prior_type, 1L)
  expect_equal(fit_scalar$backend$beta_prior_mean, c(1, 1))
  expect_equal(fit_scalar$backend$beta_prior_precision, c(0.25, 0.25))

  fit_vector <- stLMM(
    y ~ x,
    data = dat,
    n_samples = 2,
    starting = list(tau_sq = 1),
    tuning = list(tau_sq = 0),
    priors = list(
      tau_sq = ig(2, 1),
      beta = normal(mean = c(1, -1), sd = c(2, 4))
    ),
    verbose = FALSE
  )
  expect_equal(fit_vector$backend$beta_prior_mean, c(1, -1))
  expect_equal(fit_vector$backend$beta_prior_precision, c(0.25, 0.0625))
})

test_that("flat and normal priors are rejected outside beta", {
  dat <- data.frame(y = c(1, 2, 3), x = c(-1, 0, 1))

  expect_error(
    stLMM(
      y ~ x,
      data = dat,
      n_samples = 2,
      priors = list(tau_sq = flat()),
      verbose = FALSE
    ),
    "unsupported variance prior flat"
  )

  expect_error(
    stLMM(
      y ~ x,
      data = dat,
      n_samples = 2,
      priors = list(tau_sq = normal()),
      verbose = FALSE
    ),
    "unsupported variance prior normal"
  )

  expect_error(
    stLMM(
      y ~ x + nngp(x, m = 2),
      data = dat,
      n_samples = 2,
      priors = list(
        tau_sq = ig(2, 1),
        nngp_1 = list(sigma_sq = ig(2, 1), phi = flat())
      ),
      verbose = FALSE
    ),
    "not supported"
  )
})

test_that("normal beta prior contributes to beta update", {
  set.seed(24)
  dat <- data.frame(y = 0)
  fit <- stLMM(
    y ~ 1,
    data = dat,
    n_samples = 200,
    starting = list(tau_sq = 1),
    tuning = list(tau_sq = 0),
    priors = list(
      tau_sq = ig(2, 1),
      beta = normal(mean = 3, sd = 0.01)
    ),
    verbose = FALSE
  )

  expect_gt(mean(fit$beta_samples[, 1]), 2.95)
  expect_lt(abs(mean(fit$beta_samples[, 1]) - 3), 0.03)
})
