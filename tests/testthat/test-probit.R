test_that("fixed-effect probit sampler matches probit GLM", {
  set.seed(601)
  n <- 500
  x <- stats::rnorm(n)
  beta <- c(-0.25, 0.8)
  pr <- stats::pnorm(beta[1] + beta[2] * x)
  dat <- data.frame(y = stats::rbinom(n, 1, pr), x = x)

  fit <- stLMM(
    y ~ x,
    data = dat,
    family = "probit",
    n_samples = 1200,
    verbose = FALSE
  )

  keep <- 401:1200
  glm_coef <- stats::coef(stats::glm(y ~ x, data = dat,
                                     family = stats::binomial(link = "probit")))
  post_mean <- colMeans(fit$beta_samples[keep, , drop = FALSE])
  fitted_prob <- fitted(fit, summary = FALSE, sub_sample = list(start = 401, thin = 100))
  fitted_link <- fitted(fit, summary = FALSE, sub_sample = list(start = 401, thin = 100),
                        scale = "link")
  ll <- log_lik(fit, sub_sample = list(start = 401, thin = 200))

  expect_equal(fit$backend$family, "probit")
  expect_null(fit$tau_sq_samples)
  expect_true(all(abs(post_mean - glm_coef) < 0.22))
  expect_true(all(fitted_prob > 0 & fitted_prob < 1))
  expect_true(max(abs(stats::pnorm(fitted_link) - fitted_prob)) < 1e-10)
  expect_equal(ncol(ll), n)
})

test_that("fixed-effect probit sampler matches independent Albert-Chib reference", {
  sample_truncated_latent <- function(mu, y){
    p0 <- stats::pnorm(0, mean = mu, sd = 1)
    u <- stats::runif(length(mu))
    p <- ifelse(y == 1L, p0 + u * (1 - p0), u * p0)
    stats::qnorm(p, mean = mu, sd = 1)
  }

  reference_ac_probit <- function(dat, n_draw, beta_prior_sd, seed){
    set.seed(seed)
    X <- stats::model.matrix(y ~ x, dat)
    y <- as.integer(dat$y)
    p <- ncol(X)
    prior_precision <- diag(1 / beta_prior_sd^2, p)
    precision <- crossprod(X) + prior_precision
    covariance <- solve(precision)
    covariance_chol <- chol(covariance)
    beta <- rep(0, p)
    out <- matrix(NA_real_, nrow = n_draw, ncol = p)
    colnames(out) <- colnames(X)

    for(i in seq_len(n_draw)){
      z <- sample_truncated_latent(as.numeric(X %*% beta), y)
      mean <- covariance %*% crossprod(X, z)
      beta <- as.numeric(mean + t(covariance_chol) %*% stats::rnorm(p))
      out[i, ] <- beta
    }
    out
  }

  set.seed(605)
  n <- 350
  x <- stats::rnorm(n)
  beta <- c(-0.45, 0.75)
  pr <- stats::pnorm(beta[1] + beta[2] * x)
  dat <- data.frame(y = stats::rbinom(n, 1, pr), x = x)

  fit <- stLMM(
    y ~ x,
    data = dat,
    family = "probit",
    n_samples = 1800,
    warmup = FALSE,
    verbose = FALSE
  )
  stlmm_draws <- fit$beta_samples[601:1800, , drop = FALSE]

  ref_draws <- reference_ac_probit(
    dat = dat,
    n_draw = 5000,
    beta_prior_sd = 10,
    seed = 606
  )
  ref_draws <- ref_draws[2001:5000, , drop = FALSE]

  mean_diff <- abs(colMeans(stlmm_draws) - colMeans(ref_draws))
  sd_ratio <- apply(stlmm_draws, 2L, stats::sd) / apply(ref_draws, 2L, stats::sd)

  expect_true(all(mean_diff < 0.08))
  expect_true(all(sd_ratio > 0.85 & sd_ratio < 1.15))
})

test_that("stats binomial probit family is normalized to probit", {
  set.seed(602)
  dat <- data.frame(y = stats::rbinom(20, 1, 0.5), x = stats::rnorm(20))

  fit <- stLMM(
    y ~ x,
    data = dat,
    family = stats::binomial(link = "probit"),
    n_samples = 12,
    verbose = FALSE
  )

  expect_equal(fit$backend$family, "probit")
  expect_true(all(is.finite(fit$beta_samples)))
})

test_that("unsupported stats binomial links are rejected", {
  dat <- data.frame(y = c(0, 1, 0, 1, 0, 1), x = seq(-1, 1, length.out = 6))

  expect_error(
    stLMM(y ~ x, data = dat, family = stats::binomial(link = "cloglog"),
          n_samples = 4, verbose = FALSE),
    "unsupported binomial link"
  )
  expect_error(
    stLMM(y ~ x, data = dat, family = stats::binomial(link = "cauchit"),
          n_samples = 4, verbose = FALSE),
    "unsupported binomial link"
  )
})

test_that("probit sampler supports iid random effects", {
  set.seed(603)
  n_group <- 12
  n_per_group <- 10
  group <- factor(rep(seq_len(n_group), each = n_per_group))
  x <- stats::rnorm(n_group * n_per_group)
  alpha <- stats::rnorm(n_group, 0, 0.35)
  beta <- c(-0.2, 0.6)
  pr <- stats::pnorm(beta[1] + beta[2] * x + alpha[group])
  dat <- data.frame(y = stats::rbinom(length(x), 1, pr), x = x, group = group)

  fit <- stLMM(
    y ~ x + iid(group),
    data = dat,
    family = "probit",
    n_samples = 250,
    warmup = FALSE,
    verbose = FALSE,
    priors = list(iid_1 = list(sigma_sq = ig(2, 1)))
  )

  expect_equal(ncol(fit$alpha_samples), n_group)
  expect_equal(ncol(fit$iid_sigma_sq_samples), 1L)
  expect_true(all(is.finite(fit$beta_samples)))
  expect_true(all(fit$iid_sigma_sq_samples > 0))
})

test_that("probit sampler supports structured process saving and prediction", {
  set.seed(604)
  n <- 28
  time <- seq_len(n)
  x <- stats::rnorm(n)
  w <- as.numeric(stats::filter(stats::rnorm(n, 0, 0.3), filter = 0.45, method = "recursive"))
  beta <- c(-0.1, 0.7)
  pr <- stats::pnorm(beta[1] + beta[2] * x + w)
  dat <- data.frame(y = stats::rbinom(n, 1, pr), x = x, time = time)

  fit <- stLMM(
    y ~ x + ar1(time),
    data = dat,
    family = "probit",
    n_samples = 120,
    warmup = FALSE,
    verbose = FALSE,
    priors = list(ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.9, 0.9))),
    tuning = list(ar1_1 = list(sigma_sq = 0, phi = 0))
  )

  expect_false(is.null(fit$w_samples))
  expect_equal(fit$recover_iter, seq_len(120))

  rec <- recover(fit, sub_sample = list(start = 41, thin = 20))
  fitted_prob <- fitted(rec, summary = FALSE)
  fitted_link <- fitted(rec, summary = FALSE, scale = "link")
  pred <- predict(rec, newdata = data.frame(x = c(-1, 0, 1), time = c(3, 12, 24)),
                  y_samples = TRUE)

  expect_equal(dim(rec$w_samples$ar1_1), c(length(rec$recover_iter), n))
  expect_true(all(is.finite(rec$w_samples$ar1_1)))
  expect_true(max(abs(stats::pnorm(fitted_link) - fitted_prob)) < 1e-10)
  expect_true(all(pred$mu_samples > 0 & pred$mu_samples < 1))
  expect_true(all(pred$y_samples %in% c(0, 1)))
})

test_that("structured probit sampler matches independent AR1 Albert-Chib reference", {
  sample_truncated_latent <- function(mu, y){
    p0 <- stats::pnorm(0, mean = mu, sd = 1)
    u <- stats::runif(length(mu))
    p <- ifelse(y == 1L, p0 + u * (1 - p0), u * p0)
    stats::qnorm(p, mean = mu, sd = 1)
  }

  ar1_precision <- function(n, sigma_sq, phi){
    Q <- matrix(0, n, n)
    if(n == 1L){
      Q[1, 1] <- 1 / sigma_sq
      return(Q)
    }
    den <- sigma_sq * (1 - phi^2)
    Q[1, 1] <- 1 / den
    Q[n, n] <- 1 / den
    for(i in 2:(n - 1L))
      Q[i, i] <- (1 + phi^2) / den
    for(i in 1:(n - 1L)){
      Q[i, i + 1L] <- -phi / den
      Q[i + 1L, i] <- -phi / den
    }
    Q
  }

  draw_gaussian_precision <- function(precision, rhs){
    R <- chol(precision)
    mean <- backsolve(R, forwardsolve(t(R), rhs))
    as.numeric(mean + backsolve(R, stats::rnorm(length(rhs))))
  }

  reference_ac_ar1_probit <- function(dat, n_draw, sigma_sq, phi, seed){
    set.seed(seed)
    X <- stats::model.matrix(y ~ x, dat)
    y <- as.integer(dat$y)
    n <- nrow(X)
    p <- ncol(X)
    beta_prior_precision <- diag(0.01, p)
    Qw <- ar1_precision(n, sigma_sq, phi)
    precision <- rbind(
      cbind(crossprod(X) + beta_prior_precision, t(X)),
      cbind(X, diag(n) + Qw)
    )
    beta <- rep(0, p)
    w <- rep(0, n)
    beta_out <- matrix(NA_real_, nrow = n_draw, ncol = p)
    w_out <- matrix(NA_real_, nrow = n_draw, ncol = n)
    colnames(beta_out) <- colnames(X)

    for(i in seq_len(n_draw)){
      z <- sample_truncated_latent(as.numeric(X %*% beta) + w, y)
      draw <- draw_gaussian_precision(
        precision,
        c(as.numeric(crossprod(X, z)), z)
      )
      beta <- draw[seq_len(p)]
      w <- draw[p + seq_len(n)]
      beta_out[i, ] <- beta
      w_out[i, ] <- w
    }
    list(beta = beta_out, w = w_out)
  }

  set.seed(606)
  n <- 18
  sigma_sq <- 0.35
  phi <- 0.45
  time <- seq_len(n)
  x <- stats::rnorm(n)
  w <- as.numeric(stats::filter(
    stats::rnorm(n, sd = sqrt(sigma_sq * (1 - phi^2))),
    filter = phi,
    method = "recursive",
    init = stats::rnorm(1, sd = sqrt(sigma_sq))
  ))
  beta <- c(-0.25, 0.7)
  pr <- stats::pnorm(beta[1] + beta[2] * x + w)
  dat <- data.frame(y = stats::rbinom(n, 1, pr), x = x, time = time)

  fit <- stLMM(
    y ~ x + ar1(time),
    data = dat,
    family = "probit",
    n_samples = 1800,
    warmup = FALSE,
    verbose = FALSE,
    starting = list(ar1_1 = c(sigma_sq = sigma_sq, phi = phi)),
    tuning = list(ar1_1 = c(sigma_sq = 0, phi = 0)),
    priors = list(ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.9, 0.9)))
  )
  rec <- recover(fit, sub_sample = list(start = 601, thin = 1))

  ref <- reference_ac_ar1_probit(
    dat = dat,
    n_draw = 5000,
    sigma_sq = sigma_sq,
    phi = phi,
    seed = 607
  )
  ref_beta <- ref$beta[2001:5000, , drop = FALSE]
  ref_w <- ref$w[2001:5000, , drop = FALSE]

  beta_mean_diff <- abs(colMeans(rec$beta_samples[rec$recover_iter, , drop = FALSE]) -
                          colMeans(ref_beta))
  beta_sd_ratio <- apply(rec$beta_samples[rec$recover_iter, , drop = FALSE], 2L, stats::sd) /
    apply(ref_beta, 2L, stats::sd)
  w_mean_diff <- abs(colMeans(rec$w_samples$ar1_1) - colMeans(ref_w))

  expect_true(all(beta_mean_diff < 0.12))
  expect_true(all(beta_sd_ratio > 0.75 & beta_sd_ratio < 1.25))
  expect_true(mean(w_mean_diff) < 0.16)
})

test_that("probit offsets and missing responses are handled in estimation and prediction", {
  set.seed(607)
  n <- 80
  x <- stats::rnorm(n)
  off <- stats::runif(n, -0.4, 0.4)
  eta <- -0.2 + 0.9 * x + off
  y <- stats::rbinom(n, 1, stats::pnorm(eta))
  y[c(3, 19, 77)] <- NA
  dat <- data.frame(y = y, x = x, off = off)

  fit <- stLMM(
    y ~ x + offset(off),
    data = dat,
    family = "probit",
    n_samples = 300,
    warmup = FALSE,
    verbose = FALSE
  )

  fitted_link <- fitted(fit, summary = FALSE, sub_sample = list(start = 101, thin = 50),
                        scale = "link")
  fitted_prob <- fitted(fit, summary = FALSE, sub_sample = list(start = 101, thin = 50))
  ll <- log_lik(fit, sub_sample = list(start = 101, thin = 50))
  pred <- predict(
    fit,
    newdata = data.frame(x = c(-1, 0, 1), off = c(0.2, 0, -0.2)),
    y_samples = TRUE
  )

  expect_equal(length(fit$backend$observed_index), n - 3L)
  expect_equal(ncol(ll), n - 3L)
  expect_true(max(abs(stats::pnorm(fitted_link) - fitted_prob)) < 1e-10)
  expect_true(all(is.finite(ll)))
  expect_true(all(pred$mu_samples > 0 & pred$mu_samples < 1))
  expect_true(all(pred$y_samples %in% c(0, 1)))
})

test_that("probit input checks reject unsupported controls and responses", {
  dat <- data.frame(y = c(0, 1, 0, 1, 0, 1), x = seq(-1, 1, length.out = 6))

  expect_error(
    stLMM(y ~ x, data = transform(dat, y = 2), family = "probit",
          n_samples = 4, verbose = FALSE),
    "binary 0/1"
  )
  expect_error(
    stLMM(y ~ x, data = dat, family = "probit", trials = rep(1L, 6),
          n_samples = 4, verbose = FALSE),
    "trials is used only"
  )
  expect_error(
    stLMM(y ~ x + resid(model = "fixed", variance = 0.2), data = dat, family = "probit",
          n_samples = 4, verbose = FALSE),
    "resid\\(\\) is not used"
  )
  expect_error(
    stLMM(y ~ x, data = dat, family = "probit", save_process = TRUE,
          n_samples = 4, verbose = FALSE),
    "save_process can only"
  )
})
