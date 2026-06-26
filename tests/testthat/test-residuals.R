test_that("fixed residual variances remove tau_sq and recover weighted regression", {
  set.seed(320)
  n <- 24L
  dat <- data.frame(
    x = seq(-1, 1, length.out = n),
    vhat = seq(0.2, 1.1, length.out = n)
  )
  X <- model.matrix(~ x, dat)
  beta <- c(1.2, -0.7)
  dat$y <- as.numeric(X %*% beta + rnorm(n, sd = sqrt(dat$vhat)))

  fit <- stLMM(
    y ~ x + resid(model = "fixed", variance = vhat),
    data = dat,
    n_samples = 3000,
    verbose = FALSE
  )

  w <- 1 / dat$vhat
  beta_wls <- solve(crossprod(X, X * w), crossprod(X, dat$y * w))

  expect_null(fit$tau_sq_samples)
  expect_null(fit$samples$tau_sq)
  expect_false("log(tau_sq)" %in% fit$adaptive_metropolis$parameter_labels)
  expect_lt(max(abs(colMeans(fit$beta_samples) - as.vector(beta_wls))), 0.08)
})

test_that("resid formula term supports default and fixed residual variance models", {
  set.seed(327)
  n <- 18L
  dat <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    vhat = seq(0.2, 0.8, length.out = n)
  )

  fit_default <- stLMM(
    y ~ x + resid(),
    data = dat,
    starting = list(resid = list(tau_sq = 1)),
    tuning = list(resid = list(tau_sq = 0)),
    priors = list(resid = list(tau_sq = ig(2, 1))),
    n_samples = 6,
    verbose = FALSE
  )
  expect_identical(fit_default$backend$residual_model$type, "global_tau")
  expect_false(any(grepl("resid", colnames(fit_default$backend$X), fixed = TRUE)))

  fit_fixed <- stLMM(
    y ~ x + resid(model = "fixed", variance = vhat),
    data = dat,
    n_samples = 6,
    verbose = FALSE
  )
  expect_identical(fit_fixed$backend$residual_model$type, "fixed_variance")
  expect_null(fit_fixed$tau_sq_samples)
})

test_that("resid formula term supports grouped and scaled residual variance models", {
  set.seed(328)
  dat <- data.frame(
    y = rnorm(24),
    x = rnorm(24),
    group = rep(letters[1:3], each = 8),
    vhat = rep(c(0.2, 0.4, 0.6), each = 8),
    n_eff = rep(c(5, 8, 11), each = 8)
  )

  fit_group <- stLMM(
    y ~ x + resid(model = "group", group = group),
    data = dat,
    priors = list(resid = list(tau_sq = half_t(df = 3, scale = 1))),
    n_samples = 8,
    verbose = FALSE
  )
  expect_identical(fit_group$backend$residual_model$type, "group_ig_variance")
  expect_identical(fit_group$backend$residual_model$method, "prior")
  expect_equal(colnames(fit_group$residual_variance_samples), paste0("tau_sq_", letters[1:3]))

  fit_shannon <- stLMM(
    y ~ x + resid(model = "group", group = group, variance = vhat, n = n_eff, prior = "shannon"),
    data = dat,
    n_samples = 8,
    verbose = FALSE
  )
  expect_identical(fit_shannon$backend$residual_model$method, "shannon")
  expect_equal(colnames(fit_shannon$residual_variance_samples), paste0("tau_sq_", letters[1:3]))

  fit_scaled <- stLMM(
    y ~ x + resid(model = "scaled", variance = vhat, n = n_eff),
    data = dat,
    n_samples = 8,
    verbose = FALSE
  )
  expect_identical(fit_scaled$backend$residual_model$type, "scaled_variance")
  expect_equal(colnames(fit_scaled$residual_variance_samples), c("kappa", "tau0_sq"))
})

test_that("resid formula term is restricted to one Gaussian residual model", {
  dat <- data.frame(y = rbinom(8, 1, 0.5), x = rnorm(8), vhat = rep(0.2, 8))

  expect_error(
    stLMM(y ~ x + resid(), data = dat, family = "binomial", n_samples = 4, verbose = FALSE),
    "not used with family = \"binomial\"",
    fixed = TRUE
  )
  expect_error(
    stLMM(y ~ x + resid() + resid(model = "fixed", variance = vhat),
          data = dat, n_samples = 4, verbose = FALSE),
    "only one resid"
  )
  expect_error(
    stLMM(y ~ x + x:resid(), data = dat, n_samples = 4, verbose = FALSE),
    "cannot be used in interactions"
  )
  expect_error(
    stLMM(
      y ~ x + resid(model = "group", group = group),
      data = data.frame(y = rnorm(9), x = rnorm(9), group = rep(letters[1:3], each = 3)),
      priors = list(tau_sq = ig(2, 1)),
      n_samples = 4,
      verbose = FALSE
    ),
    "priors\\$tau_sq is not used"
  )
})

test_that("fixed residual variances validate observed rows only", {
  dat <- data.frame(
    y = c(1, NA, 2),
    x = c(0, 1, 2),
    vhat = c(0.4, NA, 0.8)
  )

  fit <- stLMM(
    y ~ x + resid(model = "fixed", variance = vhat),
    data = dat,
    n_samples = 4,
    verbose = FALSE
  )

  expect_identical(fit$backend$residual_model$type, "fixed_variance")
  expect_error(
    stLMM(
      y ~ x + resid(model = "fixed", variance = vhat),
      data = transform(dat, vhat = c(0.4, 0.5, -0.8)),
      n_samples = 4,
      verbose = FALSE
    ),
    "finite and positive"
  )
  expect_error(
    stLMM(
      y ~ x + resid(model = "fixed", variance = vhat),
      data = dat,
      priors = list(tau_sq = ig(2, 1)),
      n_samples = 4,
      verbose = FALSE
    ),
    "priors\\$tau_sq is not used"
  )
})

test_that("fixed residual variances are used in prediction draws", {
  set.seed(321)
  dat <- data.frame(
    y = rnorm(18),
    x = rnorm(18),
    vhat = rep(c(0.1, 0.9), each = 9)
  )

  fit <- stLMM(
    y ~ x + resid(model = "fixed", variance = vhat),
    data = dat,
    n_samples = 2500,
    verbose = FALSE
  )
  pred <- predict(fit, y_samples = TRUE)

  residual_draws <- pred$y_samples - pred$mu_samples
  empirical_var <- apply(residual_draws, 2, var)
  expect_lt(max(abs(empirical_var - dat$vhat)), 0.12)

  expect_error(
    predict(fit, newdata = data.frame(x = 0, vhat = NA_real_), y_samples = TRUE),
    "finite and positive"
  )
})

test_that("fixed residual variances enter AR1 process recovery", {
  set.seed(322)
  n <- 4L
  phi <- 0.35
  dat <- data.frame(
    y = c(0, NA, 0, 0),
    time = seq_len(n),
    vhat = c(0.4, NA, 0.8, 1.3)
  )

  fit <- stLMM(
    y ~ 0 + ar1(time) + resid(model = "fixed", variance = vhat),
    data = dat,
    starting = list(ar1_1 = c(sigma_sq = 1, phi = phi)),
    tuning = list(ar1_1 = c(sigma_sq = 0, phi = 0)),
    priors = list(ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.9, 0.9))),
    n_samples = 16000,
    verbose = FALSE
  )
  rec <- recover(fit)

  Q <- matrix(0, n, n)
  Q[1, 1] <- 1 / (1 - phi^2)
  Q[n, n] <- 1 / (1 - phi^2)
  for(i in 2:(n - 1))
    Q[i, i] <- (1 + phi^2) / (1 - phi^2)
  for(i in 1:(n - 1)){
    Q[i, i + 1] <- -phi / (1 - phi^2)
    Q[i + 1, i] <- -phi / (1 - phi^2)
  }

  obs_precision <- c(1 / dat$vhat[1], 0, 1 / dat$vhat[3], 1 / dat$vhat[4])
  target_cov <- solve(Q + diag(obs_precision, n))
  recovered_cov <- cov(rec$w_samples$ar1_1)

  expect_null(rec$tau_sq_samples)
  expect_lt(max(abs(recovered_cov - target_cov)), 0.025)
})

test_that("group IG residual variance can be fixed and matches fixed residual variance", {
  set.seed(323)
  dat <- data.frame(
    y = rnorm(12),
    x = rnorm(12),
    group = rep(letters[1:3], each = 4),
    vhat = rep(c(0.25, 0.50, 0.90), each = 4)
  )

  set.seed(900)
  fit_fixed <- stLMM(
    y ~ x + resid(model = "fixed", variance = vhat),
    data = dat,
    n_samples = 20,
    verbose = FALSE
  )

  set.seed(900)
  fit_group <- stLMM(
    y ~ x + resid(model = "group", group = group, variance = vhat, prior = "ig", shape = 5, tuning = 0),
    data = dat,
    n_samples = 20,
    verbose = FALSE
  )

  expect_null(fit_group$tau_sq_samples)
  expect_equal(colnames(fit_group$residual_variance_samples), paste0("tau_sq_", letters[1:3]))
  expect_equal(unname(fit_group$residual_variance_samples[1, ]), c(0.25, 0.50, 0.90))
  expect_equal(fit_group$beta_samples, fit_fixed$beta_samples, tolerance = 1e-12)
})

test_that("generic group residual variance can be fixed and matches fixed residual variance", {
  set.seed(325)
  group_var <- c(a = 0.25, b = 0.50, c = 0.90)
  dat <- data.frame(
    y = rnorm(12),
    x = rnorm(12),
    group = rep(letters[1:3], each = 4)
  )
  dat$vhat <- group_var[dat$group]

  set.seed(901)
  fit_fixed <- stLMM(
    y ~ x + resid(model = "fixed", variance = vhat),
    data = dat,
    n_samples = 20,
    verbose = FALSE
  )

  set.seed(901)
  fit_group <- stLMM(
    y ~ x + resid(model = "group", group = group),
    data = dat,
    starting = list(resid = list(tau_sq = group_var)),
    tuning = list(resid = list(tau_sq = 0)),
    priors = list(resid = list(tau_sq = ig(5, 1))),
    n_samples = 20,
    verbose = FALSE
  )

  expect_null(fit_group$tau_sq_samples)
  expect_equal(colnames(fit_group$residual_variance_samples), paste0("tau_sq_", letters[1:3]))
  expect_equal(unname(fit_group$residual_variance_samples[1, ]), unname(group_var))
  expect_equal(fit_group$beta_samples, fit_fixed$beta_samples, tolerance = 1e-12)
})

test_that("generic group residual variance accepts half-t and named priors", {
  set.seed(326)
  dat <- data.frame(
    y = rnorm(24),
    x = rnorm(24),
    group = rep(letters[1:3], each = 8)
  )

  fit <- stLMM(
    y ~ x + resid(model = "group", group = group),
    data = dat,
    starting = list(resid = list(tau_sq = c(a = 0.3, b = 0.6, c = 0.9))),
    tuning = list(resid = list(tau_sq = 0.02)),
    priors = list(resid = list(tau_sq = half_t(df = 3, scale = 1))),
    n_samples = 30,
    verbose = FALSE
  )

  expect_equal(dim(fit$residual_variance_samples), c(30L, 3L))
  expect_true(all(fit$residual_variance_samples > 0))
  expect_true(all(grepl("residual_variance", fit$adaptive_metropolis$parameter_labels)))

  fit_named <- stLMM(
    y ~ x + resid(model = "group", group = group),
    data = dat,
    starting = list(resid = 0.5),
    tuning = list(resid = 0),
    priors = list(resid = list(
      a = ig(3, 0.5),
      b = half_t(df = 3, scale = 1),
      c = gamma_dist(shape = 3, rate = 2)
    )),
    n_samples = 5,
    verbose = FALSE
  )

  expect_equal(colnames(fit_named$residual_variance_samples), paste0("tau_sq_", letters[1:3]))

  fit_default_controls <- stLMM(
    y ~ x + resid(model = "group", group = group),
    data = dat,
    priors = list(resid = list(tau_sq = half_t(df = 3, scale = 1))),
    n_samples = 5,
    verbose = FALSE
  )

  expect_equal(colnames(fit_default_controls$residual_variance_samples), paste0("tau_sq_", letters[1:3]))
  expect_true(all(fit_default_controls$residual_variance_samples > 0))

  expect_error(
    stLMM(
      y ~ x + resid(model = "group", group = group),
      data = dat,
      priors = list(resid = list(a = ig(3, 1), b = ig(3, 1))),
      n_samples = 5,
      verbose = FALSE
    ),
    "missing residual group prior"
  )
})

test_that("group IG residual variance is sampled and used by prediction and recovery", {
  set.seed(324)
  n <- 6L
  dat <- data.frame(
    y = rnorm(n),
    time = seq_len(n),
    group = paste0("g", seq_len(n)),
    vhat = seq(0.2, 0.7, length.out = n)
  )

  fit <- stLMM(
    y ~ 0 + ar1(time) + resid(model = "group", group = group, variance = vhat, prior = "ig", shape = 6, tuning = 0.02),
    data = dat,
    starting = list(ar1_1 = c(sigma_sq = 1, phi = 0.3)),
    tuning = list(ar1_1 = c(sigma_sq = 0, phi = 0)),
    priors = list(ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.9, 0.9))),
    n_samples = 60,
    verbose = FALSE
  )
  rec <- recover(fit)
  pred <- predict(rec, y_samples = TRUE)

  expect_null(fit$tau_sq_samples)
  expect_equal(dim(fit$residual_variance_samples), c(60L, n))
  expect_true(all(grepl("residual_variance", fit$adaptive_metropolis$parameter_labels)))
  expect_equal(dim(rec$w_samples$ar1_1), c(60L, n))
  expect_equal(dim(pred$y_samples), c(60L, n))
})

test_that("Shannon residual prior requires effective sample sizes above one", {
  dat <- data.frame(
    y = rnorm(4),
    group = c("a", "a", "b", "b"),
    vhat = c(0.2, 0.2, 0.4, 0.4),
    n_eff = c(5, 5, 1, 1)
  )

  expect_error(
    stLMM(
      y ~ 1 + resid(model = "group", group = group, variance = vhat, n = n_eff, prior = "shannon"),
      data = dat,
      n_samples = 4,
      verbose = FALSE
    ),
    "requires n > 1"
  )
  expect_error(
    stLMM(
      y ~ 1 + resid(model = "group", group = group, variance = vhat, n = n_eff, prior = "shannon", shape = 4),
      data = dat,
      n_samples = 4,
      verbose = FALSE
    ),
    "shape is only used"
  )
})

test_that("scaled residual variance can be fixed and matches fixed residual variance", {
  set.seed(325)
  dat <- data.frame(
    y = rnorm(16),
    x = rnorm(16),
    vhat = rep(c(0.2, 0.8), each = 8)
  )

  set.seed(901)
  fit_fixed <- stLMM(
    y ~ x + resid(model = "fixed", variance = vhat),
    data = dat,
    n_samples = 20,
    verbose = FALSE
  )

  set.seed(901)
  fit_scaled <- stLMM(
    y ~ x + resid(model = "scaled", variance = vhat, starting = c(kappa = 1), tuning = 0),
    data = dat,
    n_samples = 20,
    verbose = FALSE
  )

  expect_null(fit_scaled$tau_sq_samples)
  expect_equal(colnames(fit_scaled$residual_variance_samples), "kappa")
  expect_equal(unname(fit_scaled$residual_variance_samples[, 1]), rep(1, 20))
  expect_equal(fit_scaled$beta_samples, fit_fixed$beta_samples, tolerance = 1e-12)
})

test_that("sample-size aware scaled residual variance uses kappa and tau0_sq safely", {
  set.seed(326)
  n <- 10L
  dat <- data.frame(
    y = rnorm(n),
    time = seq_len(n),
    vhat = seq(0.15, 0.9, length.out = n),
    n_eff = rep(c(2, 8), length.out = n)
  )

  fit <- stLMM(
    y ~ 0 + ar1(time) + resid(model = "scaled", variance = vhat, n = n_eff, shrinkage = 5,
                              starting = c(kappa = 1.1, tau0_sq = 0.4),
                              tuning = c(kappa = 0.02, tau0_sq = 0.02)),
    data = dat,
    starting = list(ar1_1 = c(sigma_sq = 1, phi = 0.3)),
    tuning = list(ar1_1 = c(sigma_sq = 0, phi = 0)),
    priors = list(ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.9, 0.9))),
    n_samples = 60,
    verbose = FALSE
  )
  rec <- recover(fit)
  pred <- predict(rec, y_samples = TRUE)

  expect_null(fit$tau_sq_samples)
  expect_equal(colnames(fit$residual_variance_samples), c("kappa", "tau0_sq"))
  expect_true(all(fit$residual_variance_samples > 0))
  expect_true(all(c("log(kappa)", "log(tau0_sq)") %in% fit$adaptive_metropolis$parameter_labels))
  expect_equal(dim(rec$w_samples$ar1_1), c(60L, n))
  expect_equal(dim(pred$y_samples), c(60L, n))

  expect_error(
    stLMM(
      y ~ 1 + resid(model = "scaled", variance = vhat),
      data = transform(dat, vhat = replace(vhat, 1, 0)),
      n_samples = 4,
      verbose = FALSE
    ),
    "finite and positive"
  )
  expect_error(
    stLMM(
      y ~ 1 + resid(model = "scaled", variance = vhat, n = n_eff),
      data = transform(dat, n_eff = replace(n_eff, 1, 1)),
      n_samples = 4,
      verbose = FALSE
    ),
    "greater than 1"
  )
})
