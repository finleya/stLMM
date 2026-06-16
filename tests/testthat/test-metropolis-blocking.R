test_that("metropolis blocking keeps joint as the default", {
  set.seed(910)
  dat <- data.frame(
    y = rnorm(18),
    x = rnorm(18),
    lon = rep(seq_len(6), 3),
    lat = rep(c(1, 2, 3), each = 6)
  )

  fit <- suppressWarnings(stLMM(
    y ~ x + nngp(lon, lat, m = 2),
    data = dat,
    priors = list(
      tau_sq = ig(2, 1),
      nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5))
    ),
    n_samples = 30,
    warmup = FALSE,
    verbose = FALSE
  ))

  expect_equal(fit$backend$metropolis$blocking, "joint")
  expect_equal(fit$backend$metropolis$batch_length, 25L)
  expect_equal(fit$adaptive_metropolis$batch_length, 25L)
  expect_equal(fit$adaptive_metropolis$dimension, 3L)
  expect_length(fit$adaptive_metropolis$blocks, 1L)
  expect_equal(fit$adaptive_metropolis$blocks[[1]]$dimension, 3L)
  expect_equal(fit$adaptive_metropolis$blocks[[1]]$parameter_labels,
               fit$adaptive_metropolis$parameter_labels)
  expect_equal(fit$covariance_acceptance, fit$adaptive_metropolis$acceptance)
})

test_that("built-in metropolis blocking partitions active parameters", {
  set.seed(911)
  n <- 20L
  dat <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    lon = runif(n),
    lat = runif(n),
    day = rep(seq_len(5), length.out = n)
  )

  pri <- list(
    tau_sq = ig(2, 1),
    ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.8, 0.8)),
    nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5))
  )

  by_term <- suppressWarnings(stLMM(
    y ~ x + ar1(day) + nngp(lon, lat, m = 3),
    data = dat,
    priors = pri,
    n_samples = 20,
    warmup = FALSE,
    metropolis = list(blocking = "by_term"),
    verbose = FALSE
  ))

  expect_equal(length(by_term$adaptive_metropolis$blocks), 3L)
  expect_equal(vapply(by_term$adaptive_metropolis$blocks, `[[`, integer(1), "dimension"),
               c(1L, 2L, 2L))
  expect_equal(by_term$adaptive_metropolis$blocks[[1]]$parameter_labels, "log(tau_sq)")
  expect_true(all(grepl("^log\\(ar1_1_sigma_sq\\)|^ar1_1_theta_", by_term$adaptive_metropolis$blocks[[2]]$parameter_labels)))
  expect_true(all(grepl("^log\\(nngp_1_sigma_sq\\)|^nngp_1_theta_", by_term$adaptive_metropolis$blocks[[3]]$parameter_labels)))

  variance_theta <- suppressWarnings(stLMM(
    y ~ x + ar1(day) + nngp(lon, lat, m = 3),
    data = dat,
    priors = pri,
    n_samples = 20,
    warmup = FALSE,
    metropolis = "variance_theta",
    verbose = FALSE
  ))

  expect_equal(length(variance_theta$adaptive_metropolis$blocks), 2L)
  expect_equal(vapply(variance_theta$adaptive_metropolis$blocks, `[[`, integer(1), "dimension"),
               c(3L, 2L))
  expect_true(all(grepl("tau_sq|sigma_sq", variance_theta$adaptive_metropolis$blocks[[1]]$parameter_labels)))
  expect_true(all(grepl("_theta_", variance_theta$adaptive_metropolis$blocks[[2]]$parameter_labels)))

  process_variance <- suppressWarnings(stLMM(
    y ~ x + ar1(day) + nngp(lon, lat, m = 3),
    data = dat,
    priors = pri,
    n_samples = 20,
    warmup = FALSE,
    metropolis = list(blocking = "process_variance"),
    verbose = FALSE
  ))

  expect_equal(length(process_variance$adaptive_metropolis$blocks), 3L)
  expect_equal(vapply(process_variance$adaptive_metropolis$blocks, `[[`, integer(1), "dimension"),
               c(1L, 2L, 2L))
})

test_that("residual_process blocking separates sampled residual variance", {
  set.seed(912)
  n <- 8L
  dat <- data.frame(
    y = rnorm(n),
    day = seq_len(n),
    group = paste0("g", seq_len(n)),
    vhat = seq(0.2, 0.7, length.out = n)
  )

  fit <- stLMM(
    y ~ 0 + ar1(day) + resid(model = "group", group = group, variance = vhat, prior = "ig", shape = 6, tuning = 0.02),
    data = dat,
    priors = list(ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.8, 0.8))),
    n_samples = 20,
    warmup = FALSE,
    metropolis = list(blocking = "residual_process"),
    verbose = FALSE
  )

  expect_null(fit$tau_sq_samples)
  expect_equal(length(fit$adaptive_metropolis$blocks), 2L)
  expect_equal(fit$adaptive_metropolis$blocks[[1]]$dimension, n)
  expect_equal(fit$adaptive_metropolis$blocks[[2]]$dimension, 2L)
  expect_true(all(grepl("residual_variance", fit$adaptive_metropolis$blocks[[1]]$parameter_labels)))
})

test_that("scalar metropolis blocking uses one-dimensional adaptive blocks", {
  set.seed(913)
  n <- 16L
  dat <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    lon = runif(n),
    lat = runif(n)
  )

  fit <- suppressWarnings(stLMM(
    y ~ x + nngp(lon, lat, m = 3),
    data = dat,
    priors = list(
      tau_sq = ig(2, 1),
      nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5))
    ),
    n_samples = 80,
    warmup = list(enabled = TRUE, max_batches = 2),
    metropolis = list(blocking = "scalar", batch_length = 10),
    verbose = FALSE
  ))

  expect_equal(fit$backend$metropolis$blocking, "scalar")
  expect_equal(fit$backend$metropolis$target_accept, 0.44)
  expect_equal(fit$backend$metropolis$batch_length, 10L)
  expect_equal(length(fit$adaptive_metropolis$blocks), fit$adaptive_metropolis$dimension)
  expect_true(all(vapply(fit$adaptive_metropolis$blocks, `[[`, integer(1), "dimension") == 1L))
  expect_equal(fit$adaptive_metropolis$target_accept, 0.44)
  expect_true(is.finite(fit$covariance_acceptance))
})

test_that("metropolis controls expose retained adaptation settings", {
  set.seed(914)
  dat <- data.frame(y = rnorm(12), x = rnorm(12))

  fit <- stLMM(
    y ~ x,
    data = dat,
    priors = list(tau_sq = ig(2, 1)),
    n_samples = 60,
    warmup = FALSE,
    metropolis = list(target_accept = 0.31, batch_length = 12),
    verbose = FALSE
  )

  expect_equal(fit$backend$metropolis$target_accept, 0.31)
  expect_equal(fit$backend$metropolis$batch_length, 12L)
  expect_equal(fit$adaptive_metropolis$target_accept, 0.31)
  expect_equal(fit$adaptive_metropolis$batch_length, 12L)
  expect_equal(length(fit$adaptive_metropolis$batch_acceptance_history), 5L)
})
