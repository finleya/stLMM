test_that("fixed-effect negative-binomial PG sampler targets log-mean model", {
  set.seed(201)
  n <- 450
  x <- stats::rnorm(n)
  size <- 4.5
  beta <- c(0.15, -0.55)
  mu <- exp(beta[1] + beta[2] * x)
  dat <- data.frame(y = stats::rnbinom(n, size = size, mu = mu), x = x)

  fit <- stLMM(
    y ~ x,
    data = dat,
    family = "negative_binomial",
    size = size,
    n_samples = 1200,
    verbose = FALSE
  )

  keep <- 401:1200
  post_mean <- colMeans(fit$beta_samples[keep, , drop = FALSE])

  expect_equal(fit$backend$family, "negative_binomial")
  expect_equal(fit$backend$nb_size, size)
  expect_null(fit$tau_sq_samples)
  expect_true(all(abs(post_mean - beta) < c(0.18, 0.18)))
})

test_that("negative-binomial prediction returns means and count samples", {
  set.seed(202)
  n <- 180
  x <- stats::rnorm(n)
  size <- 3
  beta <- c(0.05, 0.4)
  dat <- data.frame(
    y = stats::rnbinom(n, size = size, mu = exp(beta[1] + beta[2] * x)),
    x = x
  )

  fit <- stLMM(
    y ~ x,
    data = dat,
    family = "negbin",
    size = size,
    n_samples = 300,
    verbose = FALSE
  )

  newdata <- data.frame(x = c(-0.5, 0.5, 1.0))
  pred_mean <- predict(fit, newdata = newdata, y_samples = TRUE,
                       sub_sample = list(start = 101, thin = 50))
  pred_link <- predict(fit, newdata = newdata, scale = "link",
                       sub_sample = list(start = 101, thin = 50))

  expect_equal(dim(pred_mean$mu_samples), c(4L, 3L))
  expect_true(all(pred_mean$mu_samples > 0))
  expect_true(max(abs(exp(pred_link$mu_samples) - pred_mean$mu_samples)) < 1e-10)
  expect_equal(dim(pred_mean$y_samples), c(4L, 3L))
  expect_true(all(pred_mean$y_samples >= 0))
  expect_true(all(abs(pred_mean$y_samples - round(pred_mean$y_samples)) < sqrt(.Machine$double.eps)))
})

test_that("large-size negative-binomial PG fit approximates Poisson regression", {
  set.seed(207)
  n <- 220
  x <- stats::rnorm(n)
  beta <- c(0.2, -0.45)
  mu <- exp(beta[1] + beta[2] * x)
  dat <- data.frame(y = stats::rpois(n, mu), x = x)
  size <- 1000

  fit <- stLMM(
    y ~ x,
    data = dat,
    family = "negative_binomial",
    size = size,
    n_samples = 800,
    verbose = FALSE
  )

  keep <- 301:800
  post_mean <- colMeans(fit$beta_samples[keep, , drop = FALSE])
  poisson_coef <- stats::coef(stats::glm(y ~ x, data = dat, family = stats::poisson()))

  expect_equal(fit$backend$nb_size, size)
  expect_true(all(abs(post_mean - poisson_coef) < 0.08))
})

test_that("negative-binomial sampler supports iid random effects", {
  set.seed(203)
  n_group <- 14
  n_per_group <- 10
  group <- factor(rep(seq_len(n_group), each = n_per_group))
  x <- stats::rnorm(n_group * n_per_group)
  alpha <- stats::rnorm(n_group, 0, 0.25)
  beta <- c(0.1, -0.35)
  size <- 5
  mu <- exp(beta[1] + beta[2] * x + alpha[group])
  dat <- data.frame(y = stats::rnbinom(length(x), size = size, mu = mu),
                    x = x, group = group)

  fit <- stLMM(
    y ~ x + iid(group),
    data = dat,
    family = "negative_binomial",
    size = size,
    n_samples = 350,
    warmup = FALSE,
    verbose = FALSE,
    priors = list(iid_1 = list(sigma_sq = ig(2, 1)))
  )

  expect_equal(ncol(fit$alpha_samples), n_group)
  expect_equal(ncol(fit$iid_sigma_sq_samples), 1L)
  expect_true(all(is.finite(fit$beta_samples)))
  expect_true(all(fit$iid_sigma_sq_samples > 0))
})

test_that("negative-binomial sampler supports structured process recovery", {
  set.seed(204)
  n <- 32
  time <- seq_len(n)
  x <- stats::rnorm(n)
  w <- as.numeric(stats::filter(stats::rnorm(n, 0, 0.25), filter = 0.45, method = "recursive"))
  beta <- c(0.05, 0.35)
  size <- 6
  mu <- exp(beta[1] + beta[2] * x + w)
  dat <- data.frame(y = stats::rnbinom(n, size = size, mu = mu), x = x, time = time)

  fit <- stLMM(
    y ~ x + ar1(time),
    data = dat,
    family = "negative_binomial",
    size = size,
    n_samples = 300,
    warmup = FALSE,
    verbose = FALSE,
    priors = list(ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.9, 0.9))),
    tuning = list(ar1_1 = list(sigma_sq = 0.05, phi = 0.05))
  )

  expect_true(all(is.finite(fit$beta_samples)))
  expect_true(all(fit$sigma_sq_samples > 0))
  expect_true(all(fit$theta_samples > -0.9 & fit$theta_samples < 0.9))

  expect_false(is.null(fit$w_samples))
  expect_equal(fit$recover_iter, seq_len(300))

  rec <- recover(fit, sub_sample = list(start = 101, thin = 50))
  expect_equal(rec$recover_iter, seq.int(101, 300, by = 50))
  expect_equal(dim(rec$w_samples$ar1_1), c(length(rec$recover_iter), n))
  expect_true(all(is.finite(rec$w_samples$ar1_1)))

  fitted_mean <- fitted(rec, summary = FALSE)
  fitted_link <- fitted(rec, summary = FALSE, scale = "link")
  expect_equal(dim(fitted_mean), c(length(rec$recover_iter), n))
  expect_true(all(fitted_mean > 0))
  expect_true(all(is.finite(fitted_link)))
  expect_true(max(abs(exp(fitted_link) - fitted_mean)) < 1e-10)

  pred <- predict(rec, y_samples = TRUE)
  expect_equal(dim(pred$mu_samples), c(length(rec$recover_iter), n))
  expect_true(all(pred$mu_samples > 0))
  expect_true(all(pred$y_samples >= 0))
})

test_that("negative-binomial PG process saving supports thinning and prediction", {
  set.seed(753)
  dat <- data.frame(
    y = c(0, 2, 1, 4, 3, 0, 5, 2),
    x = stats::rnorm(8),
    s = seq_len(8)
  )

  fit <- stLMM(
    y ~ x + ar1(s),
    data = dat,
    family = "negative_binomial",
    size = 20,
    save_process = list(start = 2, thin = 5),
    n_samples = 12,
    warmup = FALSE,
    verbose = FALSE,
    priors = list(ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.9, 0.9))),
    tuning = list(ar1_1 = list(sigma_sq = 0, phi = 0))
  )

  expect_equal(fit$recover_iter, as.integer(c(2, 7, 12)))
  pred <- predict(fit, y_samples = TRUE)
  expect_equal(pred$draw_index, as.integer(c(2, 7, 12)))
  expect_true(all(pred$mu_samples > 0))
  expect_true(all(pred$y_samples >= 0))
})

test_that("negative-binomial input checks describe unsupported controls", {
  dat <- data.frame(
    y = c(0, 1, 2, 0, 3, 1),
    x = stats::rnorm(6),
    group = rep(letters[1:3], each = 2),
    s = seq_len(6)
  )

  expect_error(
    stLMM(y ~ x, data = dat, family = "negative_binomial", n_samples = 4, verbose = FALSE),
    "size must be supplied"
  )

  expect_error(
    stLMM(y ~ x, data = dat, family = "negative_binomial", size = 0, n_samples = 4, verbose = FALSE),
    "size must be a finite positive scalar"
  )

  expect_error(
    stLMM(y ~ x, data = transform(dat, y = 1.5), family = "negative_binomial", size = 2, n_samples = 4, verbose = FALSE),
    "integer-valued"
  )

  expect_error(
    stLMM(y ~ x, data = dat, family = "negative_binomial", size = 2, trials = 1, n_samples = 4, verbose = FALSE),
    "trials is used only"
  )

  expect_error(
    stLMM(y ~ x + resid(), data = dat, family = "negative_binomial", size = 2, n_samples = 4, verbose = FALSE),
    "resid\\(\\) is not used"
  )

  expect_error(
    stLMM(y ~ x, data = dat, family = "binomial", size = 2, n_samples = 4, verbose = FALSE),
    "size is used only"
  )
})
