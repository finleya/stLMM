test_that("Gaussian offset is equivalent to centering the response", {
  set.seed(11)
  n <- 35
  dat <- data.frame(
    y = numeric(n),
    x = rnorm(n),
    off = seq(-1, 1, length.out = n)
  )
  dat$y <- dat$off + 0.4 - 0.7 * dat$x + rnorm(n, sd = 0.15)
  dat$y_centered <- dat$y - dat$off

  set.seed(22)
  fit_offset <- stLMM(
    y ~ x + offset(off),
    data = dat,
    n_samples = 20,
    verbose = FALSE
  )
  set.seed(22)
  fit_centered <- stLMM(
    y_centered ~ x,
    data = dat,
    n_samples = 20,
    verbose = FALSE
  )

  expect_equal(fit_offset$beta_samples, fit_centered$beta_samples)
  expect_equal(fit_offset$tau_sq_samples, fit_centered$tau_sq_samples)

  eta_offset <- fitted(fit_offset, summary = FALSE)
  eta_centered <- fitted(fit_centered, summary = FALSE)
  expect_equal(sweep(eta_centered, 2, dat$off, `+`), eta_offset)
})

test_that("binomial offset enters the PG linear predictor", {
  set.seed(12)
  n <- 90
  dat <- data.frame(
    y = integer(n),
    x = rnorm(n),
    off = seq(-2, 2, length.out = n),
    trials = 12L
  )
  eta <- dat$off + 0.25 + 0.15 * dat$x
  dat$y <- rbinom(n, size = dat$trials, prob = plogis(eta))

  fit <- stLMM(
    y ~ x + offset(off),
    data = dat,
    family = "binomial",
    trials = dat$trials,
    priors = list(beta = normal(0, 2)),
    n_samples = 350,
    verbose = FALSE
  )

  beta_hat <- colMeans(fit$beta_samples[151:350, , drop = FALSE])
  expect_lt(abs(beta_hat[["(Intercept)"]] - 0.25), 0.45)
  expect_lt(abs(beta_hat[["x"]] - 0.15), 0.45)

  eta_draws <- fitted(fit, summary = FALSE, scale = "link",
                      sub_sample = list(start = 300))
  manual <- fit$beta_samples[300:350, , drop = FALSE] %*% t(fit$backend$X)
  manual <- sweep(manual, 2, dat$off, `+`)
  attributes(eta_draws) <- attributes(manual) <- NULL
  expect_equal(eta_draws, manual)
})

test_that("negative-binomial offset enters the PG linear predictor", {
  set.seed(13)
  n <- 90
  dat <- data.frame(
    y = integer(n),
    x = rnorm(n),
    off = seq(log(0.5), log(2.5), length.out = n)
  )
  eta <- dat$off + 0.35 - 0.2 * dat$x
  dat$y <- rnbinom(n, size = 8, mu = exp(eta))

  fit <- stLMM(
    y ~ x + offset(off),
    data = dat,
    family = "negative_binomial",
    size = 8,
    priors = list(beta = normal(0, 2)),
    n_samples = 350,
    verbose = FALSE
  )

  beta_hat <- colMeans(fit$beta_samples[151:350, , drop = FALSE])
  expect_lt(abs(beta_hat[["(Intercept)"]] - 0.35), 0.45)
  expect_lt(abs(beta_hat[["x"]] + 0.2), 0.45)

  mu_hat <- fitted(fit, scale = "response")
  expect_equal(length(mu_hat), nrow(dat))
  expect_true(all(is.finite(mu_hat)))
  expect_true(all(mu_hat > 0))
})

test_that("prediction evaluates offsets in newdata", {
  set.seed(14)
  n <- 35
  dat <- data.frame(
    y = numeric(n),
    x = rnorm(n),
    off = rnorm(n)
  )
  dat$y <- dat$off + 1 - 0.5 * dat$x + rnorm(n, sd = 0.1)

  fit <- stLMM(
    y ~ x + offset(off),
    data = dat,
    n_samples = 30,
    verbose = FALSE
  )

  nd1 <- data.frame(x = c(-1, 1), off = c(0, 0))
  nd2 <- data.frame(x = c(-1, 1), off = c(2, 2))
  p1 <- predict(fit, newdata = nd1, sub_sample = list(start = 20))
  p2 <- predict(fit, newdata = nd2, sub_sample = list(start = 20))
  expect_equal(as.numeric(p2$mu_samples - p1$mu_samples), rep(2, length(p1$mu_samples)))
  expect_error(predict(fit, newdata = data.frame(x = 0)), "off")
})

test_that("NB offset works with recovered CAR effects", {
  set.seed(15)
  area <- rep(letters[1:6], each = 3)
  adj <- matrix(0, 6, 6, dimnames = list(letters[1:6], letters[1:6]))
  for(i in 1:5){
    adj[i, i + 1] <- 1
    adj[i + 1, i] <- 1
  }
  g <- car_graph(adj)
  dat <- data.frame(
    area = area,
    x = rnorm(length(area)),
    off = rep(seq(log(0.8), log(1.5), length.out = 6), each = 3)
  )
  w <- rnorm(6, sd = 0.25)
  names(w) <- letters[1:6]
  dat$y <- rnbinom(nrow(dat), size = 10,
                   mu = exp(dat$off + 0.2 + 0.1 * dat$x + w[dat$area]))

  fit <- stLMM(
    y ~ x + offset(off) + car(area, graph = g),
    data = dat,
    family = "negative_binomial",
    size = 10,
    priors = list(
      beta = normal(0, 2),
      car_1 = list(sigma_sq = half_t(3, 0.5), rho = uniform(0.05, 0.95))
    ),
    starting = list(car_1 = list(sigma_sq = 0.2, rho = 0.5)),
    tuning = list(car_1 = list(sigma_sq = 0.03, rho = 0.05)),
    n_samples = 80,
    verbose = FALSE
  )
  rec <- recover(fit, sub_sample = list(start = 41, thin = 4))
  mu <- fitted(rec, scale = "response")

  expect_equal(length(mu), nrow(dat))
  expect_true(all(is.finite(mu)))
  expect_true(all(mu > 0))
})
