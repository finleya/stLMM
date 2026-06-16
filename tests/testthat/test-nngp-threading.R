test_that("NNGP B/F updates are invariant to n_omp_threads", {
  set.seed(3001)
  n <- 180
  dat <- data.frame(
    lon = runif(n),
    lat = runif(n),
    x = rnorm(n)
  )
  dat$y <- 1 + 0.5 * dat$x + rnorm(n, sd = 0.4)

  pri_exp <- list(
    resid = list(tau_sq = half_t(df = 3, scale = 0.5)),
    nngp_1 = list(sigma_sq = half_t(df = 3, scale = 1), phi = uniform(0.1, 15))
  )

  set.seed(3002)
  fit_exp_1 <- stLMM(
    y ~ x + nngp(lon, lat, m = 10, cov_model = "exp", ordering = "coord"),
    data = dat,
    priors = pri_exp,
    n_samples = 20,
    n_omp_threads = 1,
    verbose = FALSE
  )

  set.seed(3002)
  fit_exp_2 <- suppressWarnings(stLMM(
    y ~ x + nngp(lon, lat, m = 10, cov_model = "exp", ordering = "coord"),
    data = dat,
    priors = pri_exp,
    n_samples = 20,
    n_omp_threads = 2,
    verbose = FALSE
  ))

  expect_identical(fit_exp_1$beta_samples, fit_exp_2$beta_samples)
  expect_identical(fit_exp_1$theta_samples, fit_exp_2$theta_samples)
  expect_identical(fit_exp_1$sigma_sq_samples, fit_exp_2$sigma_sq_samples)
  expect_identical(fit_exp_1$tau_sq_samples, fit_exp_2$tau_sq_samples)

  pri_matern <- list(
    resid = list(tau_sq = half_t(df = 3, scale = 0.5)),
    nngp_1 = list(
      sigma_sq = half_t(df = 3, scale = 1),
      phi = uniform(0.1, 15),
      nu = uniform(0.1, 2)
    )
  )

  set.seed(3003)
  fit_matern_1 <- stLMM(
    y ~ x + nngp(lon, lat, m = 10, cov_model = "matern", ordering = "coord"),
    data = dat,
    priors = pri_matern,
    n_samples = 20,
    n_omp_threads = 1,
    verbose = FALSE
  )

  set.seed(3003)
  fit_matern_2 <- suppressWarnings(stLMM(
    y ~ x + nngp(lon, lat, m = 10, cov_model = "matern", ordering = "coord"),
    data = dat,
    priors = pri_matern,
    n_samples = 20,
    n_omp_threads = 2,
    verbose = FALSE
  ))

  expect_identical(fit_matern_1$beta_samples, fit_matern_2$beta_samples)
  expect_identical(fit_matern_1$theta_samples, fit_matern_2$theta_samples)
  expect_identical(fit_matern_1$sigma_sq_samples, fit_matern_2$sigma_sq_samples)
  expect_identical(fit_matern_1$tau_sq_samples, fit_matern_2$tau_sq_samples)
})
