test_that("print summary and plot methods work for stLMM fits", {
  set.seed(3001)
  dat <- data.frame(y = rnorm(12), x = rnorm(12))

  fit <- stLMM(
    y ~ x,
    data = dat,
    starting = list(tau_sq = 1),
    priors = list(tau_sq = ig(2, 1)),
    n_samples = 8,
    verbose = FALSE
  )

  expect_output(print(fit), "stLMM fit")
  s <- summary(fit)
  expect_s3_class(s, "summary_stLMM")
  expect_named(s$parameters, c("beta", "tau_sq"))
  expect_true(all(c("mean", "sd", "q2.5", "q50.0", "q97.5") %in% colnames(s$parameters$beta)))
  expect_output(print(s), "stLMM summary")

  s_sel <- summary(fit, parameters = c("x", "tau_sq"))
  expect_named(s_sel$parameters, c("beta", "tau_sq"))
  expect_equal(rownames(s_sel$parameters$beta), "x")
  expect_equal(rownames(s_sel$parameters$tau_sq), "tau_sq")
  expect_error(summary(fit, parameters = "not_a_parameter"),
               "unknown parameter")

  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  expect_invisible(plot(fit, parameters = "tau_sq", burnin = 1, thin = 2))
  grDevices::dev.off()
  expect_true(file.exists(f))
})

test_that("summary can select iid random-effect parameters", {
  set.seed(3006)
  dat <- data.frame(
    y = rnorm(12),
    x = rnorm(12),
    group = factor(rep(letters[1:3], length.out = 12))
  )

  fit <- stLMM(
    y ~ x + iid(group),
    data = dat,
    starting = list(tau_sq = 1, iid_1 = list(sigma_sq = 1)),
    priors = list(tau_sq = ig(2, 1), iid_1 = list(sigma_sq = ig(2, 1))),
    n_samples = 8,
    verbose = FALSE
  )

  all_summary <- summary(fit)
  expect_true(nrow(all_summary$parameters$alpha) > 1L)

  s <- summary(fit, parameters = c("x", "iid_1_a", "tau_sq"))
  expect_named(s$parameters, c("beta", "alpha", "tau_sq"))
  expect_equal(rownames(s$parameters$beta), "x")
  expect_equal(rownames(s$parameters$alpha), "iid_1_a")
  expect_equal(rownames(s$parameters$tau_sq), "tau_sq")
  expect_null(s$parameters$iid_sigma_sq)
})

test_that("print summary and plot methods work for recovered objects", {
  set.seed(3002)
  dat <- data.frame(y = rnorm(6), time = rep(1:3, each = 2))

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

  expect_output(print(rec), "stLMM recovery")
  expect_output(print(rec), "recovered draws")
  s <- summary(rec, include_w = TRUE, max_w = 2)
  expect_s3_class(s, "summary_stLMM_recovery")
  expect_named(s$w, "ar1_1")
  expect_equal(nrow(s$w$ar1_1), 2L)
  expect_output(print(s), "Recovery")

  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  expect_invisible(plot(rec, nodes = 1:2, burnin = 1))
  expect_invisible(plot(rec, type = "fitted", thin = 2))
  grDevices::dev.off()
  expect_true(file.exists(f))
})

test_that("printed formulas do not escape quotes in character arguments", {
  set.seed(3005)
  dat <- data.frame(
    y = rnorm(12),
    x = rnorm(12),
    lon = rep(seq_len(4), 3),
    lat = rep(seq_len(3), each = 4)
  )

  fit <- suppressWarnings(stLMM(
    y ~ x + nngp(lon, lat, m = 3, cov_model = "exp", ordering = "maxmin"),
    data = dat,
    starting = list(
      tau_sq = 1,
      nngp_1 = c(sigma_sq = 1, phi = 1)
    ),
    tuning = list(
      tau_sq = 0,
      nngp_1 = c(sigma_sq = 0, phi = 0)
    ),
    priors = list(
      tau_sq = ig(2, 1),
      nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5))
    ),
    n_samples = 6,
    verbose = FALSE
  ))
  rec <- recover(fit, sub_sample = list(start = 2, thin = 2))

  fit_print <- paste(capture.output(print(fit)), collapse = "\n")
  rec_print <- paste(capture.output(print(rec)), collapse = "\n")

  expect_match(fit_print, 'cov_model = "exp"', fixed = TRUE)
  expect_match(fit_print, 'ordering = "maxmin"', fixed = TRUE)
  expect_false(grepl('\\\\\"', fit_print))
  expect_match(rec_print, 'cov_model = "exp"', fixed = TRUE)
  expect_match(rec_print, 'ordering = "maxmin"', fixed = TRUE)
  expect_false(grepl('\\\\\"', rec_print))
})

test_that("print summary and plot methods work for prediction objects", {
  set.seed(3003)
  dat <- data.frame(y = rnorm(10), x = rnorm(10))

  fit <- stLMM(
    y ~ x,
    data = dat,
    starting = list(tau_sq = 1),
    priors = list(tau_sq = ig(2, 1)),
    n_samples = 8,
    verbose = FALSE
  )

  pred <- predict(fit, newdata = data.frame(x = c(-1, 0, 1)), y_samples = TRUE)

  expect_output(print(pred), "stLMM prediction")
  s <- summary(pred)
  expect_s3_class(s, "summary_stLMM_prediction")
  expect_equal(nrow(s$mu), 3L)
  expect_equal(nrow(s$y), 3L)
  expect_output(print(s), "prediction summary")

  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  expect_invisible(plot(pred, burnin = 1))
  expect_invisible(plot(pred, type = "scatter", observed = c(-0.5, 0.25, 0.75)))
  expect_error(plot(pred, type = "scatter"), "observed values are required")
  grDevices::dev.off()
  expect_true(file.exists(f))
})
