test_that("fitted.stLMM returns fixed-effect fitted samples", {
  set.seed(10)
  dat <- data.frame(y = rnorm(8), x = rnorm(8))

  fit <- stLMM(
    y ~ x,
    data = dat,
    starting = list(tau_sq = 1),
    priors = list(tau_sq = ig(2, 1)),
    n_samples = 6,
    verbose = FALSE
  )

  mu <- fitted(fit, summary = FALSE, sub_sample = list(start = 2, thin = 2))
  draw_index <- attr(mu, "draw_index")
  expected <- fit$beta_samples[draw_index, , drop = FALSE] %*% t(fit$backend$X)
  attr(mu, "draw_index") <- NULL

  expect_s3_class(fit, "stLMM")
  expect_equal(draw_index, c(2L, 4L, 6L))
  expect_equal(unname(mu), unname(expected))
  expect_equal(fitted(fit, sub_sample = list(start = 2, thin = 2)), as.numeric(colMeans(mu)))
})

test_that("fitted.stLMM includes explicit grouped random effects", {
  set.seed(11)
  dat <- data.frame(
    y = rnorm(12),
    x = rnorm(12),
    group = factor(rep(seq_len(4), each = 3))
  )

  fit <- stLMM(
    y ~ x + iid(group),
    data = dat,
     starting = list(tau_sq = 1, iid_1 = list(sigma_sq = 1)),
    priors = list(tau_sq = ig(2, 1), iid_1 = list(sigma_sq = ig(2, 1))),
    n_samples = 5,
    verbose = FALSE
  )

  mu <- fitted(fit, summary = FALSE)
  expected <- fit$beta_samples %*% t(fit$backend$X) +
    as.matrix(fit$alpha_samples %*% Matrix::t(fit$backend$Z))
  attr(mu, "draw_index") <- NULL

  expect_equal(unname(mu), unname(expected))
})

test_that("fitted.stLMM includes recovered process effects with draw alignment", {
  set.seed(12)
  dat <- data.frame(
    y = rep(0, 6),
    time = rep(seq_len(3), each = 2)
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
    n_samples = 8,
    verbose = FALSE
  )

  rec <- recover(fit, sub_sample = list(start = 2, thin = 2))
  mu <- fitted(rec, summary = FALSE)
  map <- rec$backend$process_terms[[1]]$map
  expected <- rec$w_samples$ar1_1[, map, drop = FALSE]

  expect_equal(attr(mu, "draw_index"), rec$recover_iter)
  attr(mu, "draw_index") <- NULL
  expect_equal(unname(mu), unname(expected))
})

test_that("fitted.stLMM requires recovered process samples for process models", {
  set.seed(13)
  dat <- data.frame(y = rnorm(5), time = seq_len(5))

  fit <- stLMM(
    y ~ 0 + ar1(time),
    data = dat,
    starting = list(tau_sq = 1, ar1_1 = c(sigma_sq = 1, phi = 0.4)),
    tuning = list(tau_sq = 0, ar1_1 = c(sigma_sq = 0, phi = 0)),
    priors = list(
      tau_sq = ig(2, 1),
      ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.9, 0.9))
    ),
    n_samples = 4,
    verbose = FALSE
  )

  expect_error(fitted(fit), "require saved or recovered latent process samples")
})

test_that("standalone recover adds process draws usable by fitted.stLMM", {
  set.seed(131)
  dat <- data.frame(
    y = rep(0, 5),
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
    n_samples = 8,
    verbose = FALSE
  )

  rec <- recover(fit, sub_sample = list(start = 2, thin = 3))
  mu <- fitted(rec, summary = FALSE)
  expected <- rec$w_samples$ar1_1[, rec$backend$process_terms[[1]]$map, drop = FALSE]

  expect_s3_class(rec, "stLMM_recovery")
  expect_equal(rec$recover_iter, c(2L, 5L, 8L))
  expect_equal(attr(mu, "draw_index"), rec$recover_iter)
  attr(mu, "draw_index") <- NULL
  expect_equal(unname(mu), unname(expected))
})

test_that("standalone recover matches collapsed AR1 posterior covariance", {
  set.seed(132)
  n <- 5L
  phi <- 0.4
  dat <- data.frame(y = rep(0, n), time = seq_len(n))

  fit <- stLMM(
    y ~ 0 + ar1(time),
    data = dat,
    starting = list(tau_sq = 1, ar1_1 = c(sigma_sq = 1, phi = phi)),
    tuning = list(tau_sq = 0, ar1_1 = c(sigma_sq = 0, phi = 0)),
    priors = list(
      tau_sq = ig(2, 1),
      ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.9, 0.9))
    ),
    n_samples = 12000,
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

  target_cov <- solve(Q + diag(n))
  recovered_cov <- cov(rec$w_samples$ar1_1)

  expect_lt(max(abs(recovered_cov - target_cov)), 0.04)
})

test_that("fitted.stLMM tracks observed data and simulated grouped effects", {
  set.seed(14)
  n_group <- 10L
  n_rep <- 6L
  group <- factor(rep(seq_len(n_group), each = n_rep))
  x <- rnorm(n_group * n_rep)
  beta_true <- c("(Intercept)" = 1.1, x = -0.6)
  alpha_true <- rnorm(n_group, sd = 0.7)
  alpha_true <- alpha_true - mean(alpha_true)
  y <- beta_true[1] + beta_true[2] * x + alpha_true[as.integer(group)] +
    rnorm(length(x), sd = 0.05)
  dat <- data.frame(y = y, x = x, group = group)

  fit <- stLMM(
    y ~ x + iid(group),
    data = dat,
     starting = list(tau_sq = 0.05, iid_1 = list(sigma_sq = 0.5)),
    priors = list(tau_sq = ig(2, 0.1), iid_1 = list(sigma_sq = ig(2, 0.5))),
    n_samples = 2500,
    verbose = FALSE
  )

  sub_sample <- list(start = 1001, thin = 3)
  draw_index <- seq.int(sub_sample$start, fit$backend$n_samples, by = sub_sample$thin)
  beta_hat <- colMeans(fit$beta_samples[draw_index, , drop = FALSE])
  alpha_hat <- colMeans(fit$alpha_samples[draw_index, , drop = FALSE])
  mu_hat <- fitted(fit, sub_sample = sub_sample)

  expect_lt(sqrt(mean((mu_hat - y)^2)), 0.12)
  expect_lt(max(abs(beta_hat - beta_true[names(beta_hat)])), 0.15)
  expect_lt(sqrt(mean((as.numeric(alpha_hat - mean(alpha_hat)) - alpha_true)^2)), 0.2)
})

test_that("recovered AR1 effects and fitted values track simulated process data", {
  set.seed(15)
  n_time <- 8L
  n_rep <- 5L
  phi <- 0.55
  sigma_sq <- 0.7
  tau_sq <- 0.0025
  w_true <- numeric(n_time)
  w_true[1] <- rnorm(1, sd = sqrt(sigma_sq))
  for(i in 2:n_time)
    w_true[i] <- phi * w_true[i - 1] + rnorm(1, sd = sqrt(sigma_sq * (1 - phi^2)))

  time <- rep(seq_len(n_time), each = n_rep)
  y <- w_true[time] + rnorm(length(time), sd = sqrt(tau_sq))
  dat <- data.frame(y = y, time = time)

  fit <- stLMM(
    y ~ 0 + ar1(time),
    data = dat,
    starting = list(tau_sq = tau_sq, ar1_1 = c(sigma_sq = sigma_sq, phi = phi)),
    tuning = list(tau_sq = 0, ar1_1 = c(sigma_sq = 0, phi = 0)),
    priors = list(
      tau_sq = ig(2, 1),
      ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.9, 0.9))
    ),
    n_samples = 1500,
    verbose = FALSE
  )

  rec <- recover(fit, sub_sample = list(start = 501, thin = 2))
  w_hat <- colMeans(rec$w_samples$ar1_1)
  mu_hat <- fitted(rec)

  expect_lt(sqrt(mean((w_hat - w_true)^2)), 0.08)
  expect_lt(sqrt(mean((mu_hat - y)^2)), 0.08)
})

test_that("fitted.stLMM applies SVC scaling and recovers simulated SVC process values", {
  set.seed(16)
  n_time <- 7L
  n_rep <- 6L
  phi <- 0.45
  sigma_sq <- 0.6
  tau_sq <- 0.0025
  w_true <- numeric(n_time)
  w_true[1] <- rnorm(1, sd = sqrt(sigma_sq))
  for(i in 2:n_time)
    w_true[i] <- phi * w_true[i - 1] + rnorm(1, sd = sqrt(sigma_sq * (1 - phi^2)))

  time <- rep(seq_len(n_time), each = n_rep)
  x <- runif(length(time), 0.7, 1.8)
  beta_true <- 0.4
  y <- beta_true + x * w_true[time] + rnorm(length(time), sd = sqrt(tau_sq))
  dat <- data.frame(y = y, x = x, time = time)

  fit <- stLMM(
    y ~ 1 + x:ar1(time),
    data = dat,
    starting = list(tau_sq = tau_sq, ar1_1 = c(sigma_sq = sigma_sq, phi = phi)),
    tuning = list(tau_sq = 0, ar1_1 = c(sigma_sq = 0, phi = 0)),
    priors = list(
      tau_sq = ig(2, 1),
      ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.9, 0.9))
    ),
    n_samples = 1800,
    verbose = FALSE
  )

  rec <- recover(fit, sub_sample = list(start = 601, thin = 2))
  mu <- fitted(rec, summary = FALSE)
  map <- rec$backend$process_terms[[1]]$map
  beta_draws <- rec$beta_samples[attr(mu, "draw_index"), , drop = FALSE]
  expected <- beta_draws %*% t(rec$backend$X) +
    sweep(rec$w_samples$ar1_1[, map, drop = FALSE], 2L, x, `*`)
  attr(mu, "draw_index") <- NULL
  w_hat <- colMeans(rec$w_samples$ar1_1)
  mu_hat <- as.numeric(colMeans(mu))

  expect_equal(unname(mu), unname(expected))
  expect_lt(sqrt(mean((w_hat - w_true)^2)), 0.12)
  expect_lt(sqrt(mean((mu_hat - y)^2)), 0.1)
})

test_that("recovered NNGP effects with repeated observations fit simulated process data", {
  set.seed(17)
  coords <- cbind(
    lon = c(0, 1, 2, 0, 1, 2, 0.5, 1.5),
    lat = c(0, 0, 0, 1, 1, 1, 2, 2)
  )
  n_node <- nrow(coords)
  n_rep <- 4L
  phi <- 0.7
  sigma_sq <- 0.8
  tau_sq <- 0.0025
  dist_mat <- as.matrix(stats::dist(coords))
  cov_mat <- sigma_sq * exp(-phi * dist_mat)
  w_true <- as.vector(t(chol(cov_mat)) %*% rnorm(n_node))
  idx <- rep(seq_len(n_node), each = n_rep)
  y <- w_true[idx] + rnorm(length(idx), sd = sqrt(tau_sq))
  dat <- data.frame(y = y, lon = coords[idx, "lon"], lat = coords[idx, "lat"])

  fit <- suppressWarnings(stLMM(
    y ~ 0 + nngp(lon, lat, m = 3, ordering = "maxmin"),
    data = dat,
    starting = list(tau_sq = tau_sq, nngp_1 = c(sigma_sq = sigma_sq, phi = phi)),
    tuning = list(tau_sq = 0, nngp_1 = c(sigma_sq = 0, phi = 0)),
    priors = list(
      tau_sq = ig(2, 1),
      nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5))
    ),
    n_samples = 1500,
    verbose = FALSE
  ))

  rec <- recover(fit, sub_sample = list(start = 501, thin = 2))
  w_hat <- colMeans(rec$w_samples$nngp_1)
  w_hat_ordered <- colMeans(rec$w_samples_ordered$nngp_1)
  mu_hat <- fitted(rec)

  expect_lt(sqrt(mean((w_hat - w_true)^2)), 0.1)
  expect_equal(unname(w_hat_ordered[rec$backend$graphs[[1]]$ord_inv]), unname(w_hat))
  expect_lt(sqrt(mean((mu_hat - y)^2)), 0.08)
})
