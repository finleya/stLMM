test_that("explicit prior constructors are required", {
  dat <- data.frame(y = rnorm(8), x = rnorm(8))

  expect_error(
    stLMM(y ~ x, data = dat, priors = list(tau_sq = c(2, 1)),
          n_samples = 4, verbose = FALSE),
    "must use a prior constructor"
  )

  fit <- stLMM(
    y ~ x,
    data = dat,
    priors = list(tau_sq = ig(2, 1)),
    n_samples = 4,
    verbose = FALSE
  )
  expect_equal(fit$backend$tau_sq_prior[1], 0)
})

test_that("variance prior families are accepted for Metropolis-updated variances", {
  dat <- data.frame(y = rnorm(9), x = rnorm(9), s = seq_len(9))

  fit_tau <- stLMM(
    y ~ x,
    data = dat,
    priors = list(tau_sq = half_t(df = 4, scale = 1)),
    n_samples = 4,
    verbose = FALSE
  )
  expect_equal(fit_tau$backend$tau_sq_prior[1], 5)

  fit_proc <- stLMM(
    y ~ x + ar1(s),
    data = dat,
    priors = list(
      ar1_1 = list(
        sigma_sq = half_normal(scale = 1),
        phi = uniform(-0.8, 0.8)
      )
    ),
    n_samples = 4,
    verbose = FALSE
  )
  expect_equal(fit_proc$backend$process_terms[[1]]$sigma_sq_prior[1], 4)
})

test_that("theta priors are checked against parameter domains", {
  dat <- data.frame(
    y = rnorm(10),
    x = rnorm(10),
    lon = runif(10),
    lat = runif(10)
  )

  fit <- stLMM(
    y ~ x + nngp(lon, lat, m = 3, cov_model = "matern"),
    data = dat,
    priors = list(
      nngp_1 = list(
        sigma_sq = ig(2, 1),
        phi = log_normal(log(2), 0.5, support = c(0.1, 8)),
        nu = gamma_dist(2, 2, support = c(0.1, 3))
      )
    ),
    n_samples = 4,
    verbose = FALSE
  )
  expect_equal(fit$backend$process_terms[[1]]$theta_prior[, 1], c(2, 3))

  expect_error(
    stLMM(
      y ~ x + nngp(lon, lat, m = 3, cov_model = "matern"),
      data = dat,
      priors = list(
        nngp_1 = list(
          sigma_sq = ig(2, 1),
          phi = log_normal(log(2), 0.5),
          nu = uniform(0.1, 3)
        )
      ),
      n_samples = 4,
      verbose = FALSE
    ),
    "must declare finite theta support"
  )

  expect_error(
    stLMM(
      y ~ x + nngp(lon, lat, m = 3, cov_model = "matern"),
      data = dat,
      priors = list(
        nngp_1 = list(
          sigma_sq = ig(2, 1),
          phi = beta_dist(2, 2),
          nu = uniform(0.1, 3)
        )
      ),
      n_samples = 4,
      verbose = FALSE
    ),
    "not supported for positive theta"
  )
})

test_that("grouped random-effect variances remain conjugate IG", {
  dat <- data.frame(y = rnorm(10), x = rnorm(10), group = rep(letters[1:5], each = 2))

  expect_error(
    stLMM(
      y ~ x + iid(group),
      data = dat,
      priors = list(tau_sq = ig(2, 1), iid_1 = list(sigma_sq = half_t(4, 1))),
      n_samples = 4,
      verbose = FALSE
    ),
    "must use ig"
  )
})

test_that("beta priors are accepted for unit-interval theta parameters", {
  ids <- paste0("a", 1:6)
  adj <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 3, 3, 4, 4, 5, 5, 6),
    j = c(2, 1, 3, 2, 4, 3, 5, 4, 6, 5),
    x = 1,
    dims = c(6, 6),
    dimnames = list(ids, ids)
  )
  g <- car_graph(adj)
  dat <- data.frame(y = rnorm(6), area = ids)

  fit_car <- stLMM(
    y ~ car(area, graph = g),
    data = dat,
    priors = list(car_1 = list(sigma_sq = ig(2, 1), rho = beta_dist(2, 2))),
    n_samples = 4,
    verbose = FALSE
  )
  expect_equal(fit_car$backend$process_terms[[1]]$theta_prior[1, 1], 6)
  expect_equal(unname(fit_car$backend$process_terms[[1]]$theta_bounds[1, ]), c(0, 1))

  dat_time <- data.frame(
    y = rnorm(12),
    area = rep(ids, 2),
    time = rep(1:2, each = length(ids))
  )
  fit_car_time <- stLMM(
    y ~ car_time(area, time, graph = g),
    data = dat_time,
    priors = list(car_time_1 = list(
      sigma_sq = ig(2, 1),
      rho = beta_dist(2, 2),
      phi = uniform(-0.8, 0.8)
    )),
    n_samples = 4,
    verbose = FALSE
  )
  expect_equal(fit_car_time$backend$process_terms[[1]]$theta_prior[1, 1], 6)
})

test_that("fixed covariance parameters use starting values without priors", {
  set.seed(101)
  dat <- data.frame(
    y = rnorm(14),
    x = rnorm(14),
    lon = runif(14),
    lat = runif(14),
    time = rep(1:2, each = 7)
  )

  fit_tau <- stLMM(
    y ~ x,
    data = dat,
    starting = list(tau_sq = fixed(0.45)),
    n_samples = 6,
    verbose = FALSE,
    warmup = FALSE
  )
  expect_equal(unique(fit_tau$tau_sq_samples), 0.45)
  expect_equal(fit_tau$adaptive_metropolis$dim, 0)

  fit_process <- stLMM(
    y ~ x + nngp(lon, lat, time, m = 4, cov_model = "gneiting"),
    data = dat,
    starting = list(
      tau_sq = fixed(0.6),
      nngp_1 = list(
        sigma_sq = fixed(1.2),
        a = fixed(0.8),
        c = fixed(2.5),
        alpha = fixed(1),
        beta = fixed(0),
        gamma = fixed(0.5),
        delta = fixed(0)
      )
    ),
    n_samples = 6,
    verbose = FALSE,
    warmup = FALSE
  )
  expect_equal(unique(fit_process$tau_sq_samples), 0.6)
  expect_equal(unique(fit_process$sigma_sq_samples[, "nngp_1_sigma_sq"]), 1.2)
  expect_equal(unique(fit_process$theta_samples[, "nngp_1_a"]), 0.8)
  expect_equal(unique(fit_process$theta_samples[, "nngp_1_c"]), 2.5)
  expect_equal(unique(fit_process$theta_samples[, "nngp_1_alpha"]), 1)
  expect_equal(unique(fit_process$theta_samples[, "nngp_1_beta"]), 0)
  expect_equal(unique(fit_process$theta_samples[, "nngp_1_gamma"]), 0.5)
  expect_equal(unique(fit_process$theta_samples[, "nngp_1_delta"]), 0)
  expect_equal(fit_process$adaptive_metropolis$dim, 0)
  expect_length(fit_process$adaptive_metropolis$blocks, 0)
})

test_that("fixed process parameters can be mixed with free parameters", {
  set.seed(102)
  dat <- data.frame(
    y = rnorm(16),
    x = rnorm(16),
    lon = runif(16),
    lat = runif(16)
  )

  fit <- stLMM(
    y ~ x + nngp(lon, lat, m = 4),
    data = dat,
    starting = list(
      tau_sq = fixed(0.5),
      nngp_1 = list(sigma_sq = 1, phi = fixed(3))
    ),
    priors = list(nngp_1 = list(sigma_sq = ig(2, 1))),
    n_samples = 10,
    verbose = FALSE,
    warmup = FALSE
  )

  expect_equal(fit$adaptive_metropolis$parameter_labels, "log(nngp_1_sigma_sq)")
  expect_equal(unique(fit$tau_sq_samples), 0.5)
  expect_equal(unique(fit$theta_samples[, "nngp_1_phi"]), 3)
  expect_gt(length(unique(fit$sigma_sq_samples[, "nngp_1_sigma_sq"])), 1)

  expect_error(
    stLMM(
      y ~ x + nngp(lon, lat, m = 4),
      data = dat,
      starting = list(nngp_1 = list(sigma_sq = fixed(1), phi = fixed(3))),
      tuning = list(nngp_1 = list(sigma_sq = 0.1, phi = 0)),
      n_samples = 4,
      verbose = FALSE
    ),
    "cannot have positive tuning"
  )
})

test_that("fixed residual group variances do not require priors", {
  set.seed(103)
  dat <- data.frame(
    y = rnorm(18),
    x = rnorm(18),
    group = rep(letters[1:3], each = 6)
  )

  fit <- stLMM(
    y ~ x + resid(model = "group", group = group),
    data = dat,
    starting = list(resid = list(tau_sq = list(
      a = fixed(0.4),
      b = 0.8,
      c = fixed(1.2)
    ))),
    priors = list(resid = list(b = ig(3, 1))),
    n_samples = 10,
    verbose = FALSE,
    warmup = FALSE,
    metropolis = "scalar"
  )

  expect_equal(fit$adaptive_metropolis$parameter_labels, "log(residual_variance_2)")
  expect_equal(unique(fit$residual_variance_samples[, "a"]), 0.4)
  expect_equal(unique(fit$residual_variance_samples[, "c"]), 1.2)
  expect_gt(length(unique(fit$residual_variance_samples[, "b"])), 1)
})
