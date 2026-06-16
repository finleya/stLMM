test_that("fixed-effect Bernoulli PG sampler matches logistic likelihood", {
  set.seed(101)
  n <- 500
  x <- stats::rnorm(n)
  beta <- c(-0.35, 0.85)
  pr <- stats::plogis(beta[1] + beta[2] * x)
  dat <- data.frame(y = stats::rbinom(n, 1, pr), x = x)

  fit <- stLMM(
    y ~ x,
    data = dat,
    family = "binomial",
    n_samples = 1200,
    verbose = FALSE
  )

  keep <- 401:1200
  glm_coef <- stats::coef(stats::glm(y ~ x, data = dat, family = stats::binomial()))
  post_mean <- colMeans(fit$beta_samples[keep, , drop = FALSE])

  expect_equal(fit$backend$family, "binomial")
  expect_null(fit$tau_sq_samples)
  expect_true(all(abs(post_mean - glm_coef) < 0.18))
})

test_that("binomial PG sampler uses supplied trial counts", {
  set.seed(102)
  n <- 350
  x <- stats::rnorm(n)
  trials <- sample(2:7, n, replace = TRUE)
  beta <- c(0.25, -0.7)
  pr <- stats::plogis(beta[1] + beta[2] * x)
  dat <- data.frame(y = stats::rbinom(n, trials, pr), x = x, trials = trials)

  fit <- stLMM(
    y ~ x,
    data = dat,
    family = "binomial",
    trials = "trials",
    n_samples = 1200,
    verbose = FALSE
  )

  keep <- 401:1200
  glm_coef <- stats::coef(stats::glm(cbind(y, trials - y) ~ x, data = dat, family = stats::binomial()))
  post_mean <- colMeans(fit$beta_samples[keep, , drop = FALSE])

  expect_equal(fit$backend$trials, trials)
  expect_true(all(abs(post_mean - glm_coef) < 0.15))
})

test_that("binomial PG sampler handles larger trial counts", {
  set.seed(107)
  n <- 90
  x <- stats::rnorm(n)
  trials <- sample(c(1L, 2L, 5L, 20L, 50L, 200L), n, replace = TRUE)
  beta <- c(-0.2, 0.45)
  pr <- stats::plogis(beta[1] + beta[2] * x)
  dat <- data.frame(y = stats::rbinom(n, trials, pr), x = x, trials = trials)

  fit <- stLMM(
    y ~ x,
    data = dat,
    family = "binomial",
    trials = "trials",
    n_samples = 350,
    verbose = FALSE
  )

  newdata <- data.frame(x = c(-1, 0, 1), trials = c(1L, 50L, 200L))
  pred <- predict(fit, newdata = newdata, y_samples = TRUE,
                  sub_sample = list(start = 151, thin = 50))

  expect_equal(fit$backend$trials, trials)
  expect_true(all(is.finite(fit$beta_samples)))
  expect_true(all(pred$mu_samples > 0 & pred$mu_samples < 1))
  expect_true(all(pred$y_samples >= 0))
  expect_true(all(sweep(pred$y_samples, 2L, newdata$trials, `<=`)))
})

test_that("binomial sampler supports iid random effects", {
  set.seed(103)
  n_group <- 18
  n_per_group <- 12
  group <- factor(rep(seq_len(n_group), each = n_per_group))
  x <- stats::rnorm(n_group * n_per_group)
  alpha <- stats::rnorm(n_group, 0, 0.45)
  beta <- c(-0.25, 0.65)
  pr <- stats::plogis(beta[1] + beta[2] * x + alpha[group])
  dat <- data.frame(y = stats::rbinom(length(x), 1, pr), x = x, group = group)

  fit <- stLMM(
    y ~ x + iid(group),
    data = dat,
    family = "binomial",
    n_samples = 500,
    warmup = FALSE,
    verbose = FALSE,
    priors = list(iid_1 = list(sigma_sq = ig(2, 1)))
  )

  expect_equal(ncol(fit$alpha_samples), n_group)
  expect_equal(ncol(fit$iid_sigma_sq_samples), 1L)
  expect_true(all(is.finite(fit$beta_samples)))
  expect_true(all(fit$iid_sigma_sq_samples > 0))
})

test_that("binomial sampler supports structured process fitting", {
  set.seed(104)
  n <- 36
  time <- seq_len(n)
  x <- stats::rnorm(n)
  w <- as.numeric(stats::filter(stats::rnorm(n, 0, 0.35), filter = 0.5, method = "recursive"))
  beta <- c(-0.15, 0.75)
  pr <- stats::plogis(beta[1] + beta[2] * x + w)
  dat <- data.frame(y = stats::rbinom(n, 1, pr), x = x, time = time)

  fit <- stLMM(
    y ~ x + ar1(time),
    data = dat,
    family = "binomial",
    n_samples = 350,
    warmup = FALSE,
    verbose = FALSE,
    priors = list(ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.9, 0.9))),
    tuning = list(ar1_1 = list(sigma_sq = 0.05, phi = 0.05))
  )

  expect_true(all(is.finite(fit$beta_samples)))
  expect_true(all(fit$sigma_sq_samples > 0))
  expect_true(all(fit$theta_samples > -0.9 & fit$theta_samples < 0.9))

  expect_false(is.null(fit$w_samples))
  expect_equal(fit$recover_iter, seq_len(350))

  rec <- recover(fit, sub_sample = list(start = 101, thin = 50))
  expect_equal(rec$recover_iter, seq.int(101, 350, by = 50))
  expect_equal(dim(rec$w_samples$ar1_1), c(length(rec$recover_iter), n))
  expect_true(all(is.finite(rec$w_samples$ar1_1)))

  fitted_prob <- fitted(rec, summary = FALSE)
  fitted_link <- fitted(rec, summary = FALSE, scale = "link")
  expect_equal(dim(fitted_prob), c(length(rec$recover_iter), n))
  expect_true(all(fitted_prob > 0 & fitted_prob < 1))
  expect_true(all(is.finite(fitted_link)))
  expect_true(max(abs(stats::plogis(fitted_link) - fitted_prob)) < 1e-10)

  expect_error(
    recover(fit, sub_sample = list(start = 101, thin = 50, pg_iter = 0)),
    "no longer supported"
  )

  fit_no_save <- stLMM(
    y ~ x + ar1(time),
    data = dat,
    family = "binomial",
    n_samples = 20,
    warmup = FALSE,
    save_process = FALSE,
    verbose = FALSE,
    priors = list(ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.9, 0.9))),
    tuning = list(ar1_1 = list(sigma_sq = 0, phi = 0))
  )
  expect_error(
    recover(fit_no_save),
    "requires process draws saved during stLMM"
  )
})

test_that("binomial prediction returns probabilities and count samples", {
  set.seed(105)
  n <- 160
  x <- stats::rnorm(n)
  trials <- sample(2:6, n, replace = TRUE)
  beta <- c(0.2, -0.55)
  pr <- stats::plogis(beta[1] + beta[2] * x)
  dat <- data.frame(y = stats::rbinom(n, trials, pr), x = x, trials = trials)

  fit <- stLMM(
    y ~ x,
    data = dat,
    family = "binomial",
    trials = "trials",
    n_samples = 250,
    verbose = FALSE
  )

  newdata <- data.frame(x = c(-0.5, 0.5, 1.0), trials = c(2, 4, 6))
  pred_prob <- predict(fit, newdata = newdata, y_samples = TRUE,
                       sub_sample = list(start = 101, thin = 50))
  pred_link <- predict(fit, newdata = newdata, scale = "link",
                       sub_sample = list(start = 101, thin = 50))

  expect_equal(dim(pred_prob$mu_samples), c(3L, 3L))
  expect_equal(attr(pred_prob$mu_samples, "dim"), attr(pred_link$mu_samples, "dim"))
  expect_true(all(pred_prob$mu_samples > 0 & pred_prob$mu_samples < 1))
  expect_true(max(abs(stats::plogis(pred_link$mu_samples) - pred_prob$mu_samples)) < 1e-10)
  expect_equal(dim(pred_prob$y_samples), c(3L, 3L))
  expect_true(all(pred_prob$y_samples >= 0))
  expect_true(all(sweep(pred_prob$y_samples, 2L, newdata$trials, `<=`)))

  pred_default_trials <- predict(fit, newdata = data.frame(x = 0),
                                 y_samples = TRUE, sub_sample = list(start = 101, thin = 50))
  expect_true(all(pred_default_trials$y_samples %in% c(0, 1)))
})

test_that("binomial prediction works after process recovery", {
  set.seed(106)
  n <- 30
  time <- seq_len(n)
  x <- stats::rnorm(n)
  w <- as.numeric(stats::filter(stats::rnorm(n, 0, 0.3), filter = 0.4, method = "recursive"))
  beta <- c(-0.1, 0.5)
  pr <- stats::plogis(beta[1] + beta[2] * x + w)
  dat <- data.frame(y = stats::rbinom(n, 1, pr), x = x, time = time)

  fit <- stLMM(
    y ~ x + ar1(time),
    data = dat,
    family = "binomial",
    n_samples = 300,
    warmup = FALSE,
    verbose = FALSE,
    priors = list(ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.9, 0.9))),
    tuning = list(ar1_1 = list(sigma_sq = 0.05, phi = 0.05))
  )
  rec <- recover(fit, sub_sample = list(start = 101, thin = 50))

  pred_prob <- predict(rec, y_samples = TRUE)
  pred_link <- predict(rec, scale = "link")
  fitted_link <- fitted(rec, summary = FALSE, scale = "link")

  expect_equal(dim(pred_prob$mu_samples), c(length(rec$recover_iter), n))
  expect_true(all(pred_prob$mu_samples > 0 & pred_prob$mu_samples < 1))
  expect_true(max(abs(stats::plogis(pred_link$mu_samples) - pred_prob$mu_samples)) < 1e-10)
  attr(fitted_link, "draw_index") <- NULL
  attr(fitted_link, "scale") <- NULL
  expect_equal(unname(pred_link$mu_samples), unname(fitted_link))
  expect_true(all(pred_prob$y_samples %in% c(0, 1)))

  newdata <- data.frame(x = c(0, 0.5), time = c(2.5, 2.5), trials = c(3, 5))
  pred_new <- predict(rec, newdata = newdata, y_samples = TRUE)
  expect_equal(dim(pred_new$mu_samples), c(length(rec$recover_iter), 2L))
  expect_true(all(pred_new$mu_samples > 0 & pred_new$mu_samples < 1))
  expect_true(all(sweep(pred_new$y_samples, 2L, newdata$trials, `<=`)))
})

test_that("binomial PG process saving supports default, thinning, prediction, and opt-out", {
  set.seed(246)
  dat <- data.frame(
    y = c(0, 1, 1, 0, 1, 0, 1, 0),
    x = stats::rnorm(8),
    s = seq_len(8)
  )
  common_priors <- list(ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.9, 0.9)))
  common_tuning <- list(ar1_1 = list(sigma_sq = 0, phi = 0))

  fit_default <- stLMM(
    y ~ x + ar1(s),
    data = dat,
    family = "binomial",
    n_samples = 12,
    warmup = FALSE,
    verbose = FALSE,
    priors = common_priors,
    tuning = common_tuning
  )
  expect_equal(fit_default$recover_iter, seq_len(12))
  pred_default <- predict(fit_default, sub_sample = list(start = 3, thin = 3))
  expect_equal(pred_default$draw_index, as.integer(c(3, 6, 9, 12)))
  rec_default <- recover(fit_default, sub_sample = list(start = 4, thin = 4))
  expect_equal(rec_default$recover_iter, as.integer(c(4, 8, 12)))

  fit_thin <- stLMM(
    y ~ x + ar1(s),
    data = dat,
    family = "binomial",
    save_process = list(start = 3, thin = 4),
    n_samples = 12,
    warmup = FALSE,
    verbose = FALSE,
    priors = common_priors,
    tuning = common_tuning
  )
  expect_equal(fit_thin$recover_iter, as.integer(c(3, 7, 11)))
  rec_thin <- recover(fit_thin, sub_sample = list(start = 1, thin = 2))
  expect_equal(rec_thin$recover_iter, as.integer(c(3, 7, 11)))
  pred_thin <- predict(fit_thin)
  expect_equal(pred_thin$draw_index, as.integer(c(3, 7, 11)))

  fit_off <- stLMM(
    y ~ x + ar1(s),
    data = dat,
    family = "binomial",
    save_process = FALSE,
    n_samples = 12,
    warmup = FALSE,
    verbose = FALSE,
    priors = common_priors,
    tuning = common_tuning
  )
  expect_null(fit_off$recover_iter)
  expect_error(recover(fit_off), "requires process draws saved during stLMM")
  expect_error(predict(fit_off), "requires saved or recovered latent process samples")
})

test_that("binomial input checks describe unsupported controls", {
  dat <- data.frame(
    y = c(0, 1, 1, 0, 1, 0),
    x = stats::rnorm(6),
    group = rep(letters[1:3], each = 2),
    s = seq_len(6)
  )

  expect_error(
    stLMM(y ~ x, data = dat, family = "binomial", starting = list(tau_sq = 1), n_samples = 4, verbose = FALSE),
    "starting\\$tau_sq"
  )

  expect_error(
    stLMM(y ~ x, data = transform(dat, y = 2), family = "binomial", n_samples = 4, verbose = FALSE),
    "less than or equal to trials"
  )

  expect_error(
    stLMM(y ~ x, data = dat, family = "binomial", trials = c(1, 1), n_samples = 4, verbose = FALSE),
    "length matching"
  )

  expect_error(
    stLMM(y ~ x, data = dat, family = "binomial", save_process = TRUE, n_samples = 4, verbose = FALSE),
    "only be enabled for Polya-Gamma models with structured process terms"
  )

  expect_error(
    stLMM(
      y ~ x + ar1(s),
      data = transform(dat, y = stats::rnorm(6)),
      save_process = TRUE,
      n_samples = 4,
      verbose = FALSE,
      priors = list(ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.9, 0.9))),
      tuning = list(ar1_1 = list(sigma_sq = 0, phi = 0))
    ),
    "only be enabled for Polya-Gamma models with structured process terms"
  )

  fit <- stLMM(
    y ~ x + ar1(s),
    data = dat,
    family = "binomial",
    save_process = TRUE,
    n_samples = 4,
    warmup = FALSE,
    verbose = FALSE,
    priors = list(ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.9, 0.9))),
    tuning = list(ar1_1 = list(sigma_sq = 0, phi = 0))
  )
  expect_equal(fit$backend$save_process$start, 1L)
  expect_equal(fit$backend$save_process$thin, 1L)
  expect_equal(fit$recover_iter, seq_len(4))
})
