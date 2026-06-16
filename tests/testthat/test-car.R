simulate_car <- function(adj, sigma_sq, rho){
  adj <- as.matrix(adj)
  degree <- rowSums(adj != 0)
  Q <- (diag(degree) - rho * adj) / sigma_sq
  L <- chol(Q)
  z <- rnorm(nrow(adj))
  as.numeric(backsolve(L, z))
}

leroux_precision <- function(adj, sigma_sq, rho){
  adj <- as.matrix(adj)
  degree <- rowSums(adj != 0)
  ((1 - rho) * diag(nrow(adj)) + rho * (diag(degree) - adj)) / sigma_sq
}

car_time_exp_precision <- function(time, lambda){
  time <- sort(unique(as.numeric(time)))
  n_time <- length(time)
  Q <- matrix(0, n_time, n_time)
  if(n_time == 1L){
    Q[1, 1] <- 1
    return(Q)
  }

  phi <- exp(-lambda * diff(time))
  den <- 1 - phi^2
  Q[1, 1] <- 1 / den[1]
  Q[n_time, n_time] <- 1 / den[n_time - 1L]
  if(n_time > 2L){
    for(i in 2:(n_time - 1L))
      Q[i, i] <- 1 / den[i - 1L] + phi[i]^2 / den[i]
  }
  for(i in seq_len(n_time - 1L)){
    Q[i, i + 1L] <- -phi[i] / den[i]
    Q[i + 1L, i] <- Q[i, i + 1L]
  }
  Q
}

test_that("car_graph normalizes symmetric adjacency matrices", {
  adj <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 3, 3, 4),
    j = c(2, 1, 3, 2, 4, 3),
    x = 1,
    dims = c(4, 4),
    dimnames = list(letters[1:4], letters[1:4])
  )

  g <- car_graph(adj)

  expect_s3_class(g, "stLMM_car_graph")
  expect_equal(g$ids, letters[1:4])
  expect_equal(g$degree, c(1, 2, 2, 1))
  expect_equal(g$n, 4L)
})

test_that("CAR process fits, recovers full support, and predicts existing areas", {
  set.seed(301)
  n_area <- 8L
  ids <- paste0("a", seq_len(n_area))
  adj <- Matrix::sparseMatrix(
    i = c(seq_len(n_area - 1L), seq_len(n_area - 1L) + 1L),
    j = c(seq_len(n_area - 1L) + 1L, seq_len(n_area - 1L)),
    x = 1,
    dims = c(n_area, n_area),
    dimnames = list(ids, ids)
  )
  g <- car_graph(adj)

  sigma_sq <- 0.7
  rho <- 0.45
  tau_sq <- 0.05
  beta <- c("(Intercept)" = 0.4, x = -0.3)
  w <- simulate_car(adj, sigma_sq, rho)

  dat <- data.frame(
    area = ids,
    x = rnorm(n_area)
  )
  mu <- beta["(Intercept)"] + beta["x"] * dat$x + w
  dat$y <- mu + rnorm(n_area, sd = sqrt(tau_sq))
  dat$y[c(3, 7)] <- NA_real_

  fit <- stLMM(
    y ~ x + car(area, graph = g),
    data = dat,
    starting = list(tau_sq = tau_sq, car_1 = c(sigma_sq = sigma_sq, rho = rho)),
    tuning = list(tau_sq = 0, car_1 = c(sigma_sq = 0, rho = 0)),
    priors = list(
      tau_sq = ig(2, tau_sq),
      car_1 = list(sigma_sq = ig(2, sigma_sq), rho = uniform(0.05, 0.95))
    ),
    n_samples = 10,
    verbose = FALSE
  )

  expect_equal(fit$backend$n, n_area)
  expect_equal(fit$backend$n_obs, n_area - 2L)
  expect_equal(fit$backend$n_missing_response, 2L)
  expect_equal(fit$backend$graphs[[1]]$graph_type, "car")
  expect_equal(fit$backend$process_terms[[1]]$map, seq_len(n_area))

  rec <- recover(fit, sub_sample = list(start = 2, thin = 2))
  expect_equal(ncol(rec$w_samples$car_1), n_area)
  expect_equal(length(fitted(rec)), n_area)

  pred_full <- predict(rec)
  expect_equal(ncol(pred_full$mu_samples), n_area)

  pred_existing <- predict(rec, newdata = dat[c(1, 4), c("area", "x"), drop = FALSE])
  expect_equal(ncol(pred_existing$mu_samples), 2L)

  bad <- dat[1, c("area", "x"), drop = FALSE]
  bad$area <- "missing_area"
  expect_error(predict(rec, newdata = bad), "new CAR area value")
})

test_that("CAR recovery covariance matches analytic observed-row target", {
  set.seed(303)
  n_area <- 7L
  ids <- paste0("a", seq_len(n_area))
  adj <- Matrix::sparseMatrix(
    i = c(seq_len(n_area - 1L), seq_len(n_area - 1L) + 1L),
    j = c(seq_len(n_area - 1L) + 1L, seq_len(n_area - 1L)),
    x = 1,
    dims = c(n_area, n_area),
    dimnames = list(ids, ids)
  )
  g <- car_graph(adj)

  observed_count <- c(3L, 0L, 2L, 1L, 0L, 4L, 1L)
  rows_per_area <- pmax(observed_count, 1L)
  area <- rep(ids, rows_per_area)
  y <- numeric(length(area))
  y[area %in% ids[observed_count == 0L]] <- NA_real_

  dat <- data.frame(y = y, area = area)
  sigma_sq <- 1
  rho <- 0.4
  tau_sq <- 0.7

  fit <- stLMM(
    y ~ 0 + car(area, graph = g),
    data = dat,
    starting = list(tau_sq = tau_sq, car_1 = c(sigma_sq = sigma_sq, rho = rho)),
    tuning = list(tau_sq = 0, car_1 = c(sigma_sq = 0, rho = 0)),
    priors = list(
      tau_sq = ig(2, tau_sq),
      car_1 = list(sigma_sq = ig(2, sigma_sq), rho = uniform(0.05, 0.95))
    ),
    n_samples = 20000,
    verbose = FALSE
  )
  rec <- recover(fit)

  adj_dense <- as.matrix(adj)
  degree <- rowSums(adj_dense != 0)
  Q <- (diag(degree) - rho * adj_dense) / sigma_sq
  target_cov <- solve(Q + diag(observed_count / tau_sq, n_area))
  recovered_cov <- cov(rec$w_samples$car_1)

  expect_equal(fit$backend$n_missing_response, sum(observed_count == 0L))
  expect_lt(max(abs(recovered_cov - target_cov)), 0.03)
})

test_that("CAR terms reuse graph structure and SVC scaling with missing nodes", {
  set.seed(306)
  n_area <- 6L
  ids <- paste0("a", seq_len(n_area))
  adj <- Matrix::sparseMatrix(
    i = c(seq_len(n_area - 1L), seq_len(n_area - 1L) + 1L),
    j = c(seq_len(n_area - 1L) + 1L, seq_len(n_area - 1L)),
    x = 1,
    dims = c(n_area, n_area),
    dimnames = list(ids, ids)
  )
  g <- car_graph(adj)

  dat <- data.frame(
    y = rnorm(n_area),
    area = ids,
    x_svc = seq(0.5, 1.5, length.out = n_area)
  )
  dat$y[c(2, 5)] <- NA_real_

  fit <- stLMM(
    y ~ 0 + car(area, graph = g) + x_svc:car(area, graph = g),
    data = dat,
    starting = list(
      tau_sq = 0.3,
      car_1 = c(sigma_sq = 0.8, rho = 0.4),
      car_2 = c(sigma_sq = 0.6, rho = 0.3)
    ),
    tuning = list(
      tau_sq = 0,
      car_1 = c(sigma_sq = 0, rho = 0),
      car_2 = c(sigma_sq = 0, rho = 0)
    ),
    priors = list(
      tau_sq = ig(2, 0.3),
      car_1 = list(sigma_sq = ig(2, 0.8), rho = uniform(0.05, 0.95)),
      car_2 = list(sigma_sq = ig(2, 0.6), rho = uniform(0.05, 0.95))
    ),
    n_samples = 12,
    verbose = FALSE
  )

  expect_length(fit$backend$graphs, 1L)
  expect_equal(unname(vapply(fit$backend$process_terms, `[[`, integer(1), "graph_index")),
               c(1L, 1L))
  expect_false(is.null(fit$backend$process_terms[[2]]$x))
  expect_equal(fit$backend$process_terms[[2]]$x, dat$x_svc)
  expect_equal(fit$backend$process_terms_obs[[2]]$x, dat$x_svc[!is.na(dat$y)])
  expect_equal(fit$backend$n_missing_response, 2L)

  rec <- recover(fit, sub_sample = list(start = 2, thin = 2))
  expect_equal(ncol(rec$w_samples$car_1), n_area)
  expect_equal(ncol(rec$w_samples$car_2), n_area)

  fitted_draws <- fitted(rec, summary = FALSE)
  manual <- rec$w_samples$car_1 + sweep(rec$w_samples$car_2, 2L, dat$x_svc, `*`)
  expect_equal(as.vector(fitted_draws), as.vector(manual),
               tolerance = 1e-10)

  newdata <- dat[c(1, 4), c("area", "x_svc"), drop = FALSE]
  pred <- predict(rec, newdata = newdata)
  manual_pred <- rec$w_samples$car_1[, c(1, 4), drop = FALSE] +
    sweep(rec$w_samples$car_2[, c(1, 4), drop = FALSE], 2L, newdata$x_svc, `*`)
  expect_equal(as.vector(pred$mu_samples), as.vector(manual_pred),
               tolerance = 1e-10)
})

test_that("Leroux CAR uses iid-to-ICAR precision with existing graph output", {
  set.seed(309)
  n_area <- 7L
  ids <- paste0("a", seq_len(n_area))
  adj <- Matrix::sparseMatrix(
    i = c(seq_len(n_area - 1L), seq_len(n_area - 1L) + 1L),
    j = c(seq_len(n_area - 1L) + 1L, seq_len(n_area - 1L)),
    x = 1,
    dims = c(n_area, n_area),
    dimnames = list(ids, ids)
  )
  g <- car_graph(adj)

  observed_count <- c(3L, 0L, 2L, 1L, 0L, 4L, 1L)
  rows_per_area <- pmax(observed_count, 1L)
  area <- rep(ids, rows_per_area)
  y <- numeric(length(area))
  y[area %in% ids[observed_count == 0L]] <- NA_real_

  dat <- data.frame(y = y, area = area)
  sigma_sq <- 0.9
  rho <- 0.55
  tau_sq <- 0.6

  fit <- stLMM(
    y ~ 0 + car(area, graph = g, car_model = "leroux"),
    data = dat,
    starting = list(tau_sq = tau_sq, car_1 = c(sigma_sq = sigma_sq, rho = rho)),
    tuning = list(tau_sq = 0, car_1 = c(sigma_sq = 0, rho = 0)),
    priors = list(
      tau_sq = ig(2, tau_sq),
      car_1 = list(sigma_sq = ig(2, sigma_sq), rho = uniform(0.01, 0.99))
    ),
    n_samples = 20000,
    verbose = FALSE
  )
  rec <- recover(fit)

  target_cov <- solve(leroux_precision(adj, sigma_sq, rho) +
                        diag(observed_count / tau_sq, n_area))
  recovered_cov <- cov(rec$w_samples$car_1)

  expect_equal(fit$backend$graphs[[1]]$car_model, "leroux")
  expect_equal(fit$backend$process_terms[[1]]$params$car_model, "leroux")
  expect_lt(max(abs(recovered_cov - target_cov)), 0.03)
})

test_that("proper and Leroux CAR terms with the same adjacency use separate graphs", {
  set.seed(310)
  ids <- paste0("a", 1:5)
  adj <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 3, 3, 4, 4, 5),
    j = c(2, 1, 3, 2, 4, 3, 5, 4),
    x = 1,
    dims = c(5, 5),
    dimnames = list(ids, ids)
  )
  g <- car_graph(adj)
  dat <- data.frame(y = rnorm(5), area = ids, x_svc = seq(0.8, 1.2, length.out = 5))

  fit <- stLMM(
    y ~ 0 + car(area, graph = g) + x_svc:car(area, graph = g, car_model = "leroux"),
    data = dat,
    starting = list(
      tau_sq = 0.2,
      car_1 = c(sigma_sq = 0.8, rho = 0.4),
      car_2 = c(sigma_sq = 0.6, rho = 0.3)
    ),
    tuning = list(
      tau_sq = 0,
      car_1 = c(sigma_sq = 0, rho = 0),
      car_2 = c(sigma_sq = 0, rho = 0)
    ),
    priors = list(
      tau_sq = ig(2, 0.2),
      car_1 = list(sigma_sq = ig(2, 0.8), rho = uniform(0.05, 0.95)),
      car_2 = list(sigma_sq = ig(2, 0.6), rho = uniform(0.05, 0.95))
    ),
    n_samples = 4,
    verbose = FALSE
  )

  expect_equal(fit$term_description$global$n_graphs, 2L)
  expect_equal(unname(vapply(fit$backend$graphs, `[[`, character(1), "car_model")),
               c("proper", "leroux"))
})

test_that("CAR graph reuse is based on adjacency structure", {
  set.seed(308)
  ids <- paste0("a", 1:5)
  adj_path <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 3, 3, 4, 4, 5),
    j = c(2, 1, 3, 2, 4, 3, 5, 4),
    x = 1,
    dims = c(5, 5),
    dimnames = list(ids, ids)
  )
  adj_path_copy <- adj_path
  adj_cycle <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 3, 3, 4, 4, 5, 5, 1),
    j = c(2, 1, 3, 2, 4, 3, 5, 4, 1, 5),
    x = 1,
    dims = c(5, 5),
    dimnames = list(ids, ids)
  )
  g1 <- car_graph(adj_path)
  g2 <- car_graph(adj_path_copy)
  g3 <- car_graph(adj_cycle)

  dat <- data.frame(
    y = rnorm(5),
    area = ids,
    x_svc = seq(0.8, 1.2, length.out = 5)
  )

  fit_same <- stLMM(
    y ~ 0 + car(area, graph = g1) + x_svc:car(area, graph = g2),
    data = dat,
    starting = list(
      tau_sq = 0.2,
      car_1 = c(sigma_sq = 0.8, rho = 0.4),
      car_2 = c(sigma_sq = 0.6, rho = 0.3)
    ),
    tuning = list(
      tau_sq = 0,
      car_1 = c(sigma_sq = 0, rho = 0),
      car_2 = c(sigma_sq = 0, rho = 0)
    ),
    priors = list(
      tau_sq = ig(2, 0.2),
      car_1 = list(sigma_sq = ig(2, 0.8), rho = uniform(0.05, 0.95)),
      car_2 = list(sigma_sq = ig(2, 0.6), rho = uniform(0.05, 0.95))
    ),
    n_samples = 4,
    verbose = FALSE
  )

  fit_different <- stLMM(
    y ~ 0 + car(area, graph = g1) + x_svc:car(area, graph = g3),
    data = dat,
    starting = list(
      tau_sq = 0.2,
      car_1 = c(sigma_sq = 0.8, rho = 0.4),
      car_2 = c(sigma_sq = 0.6, rho = 0.3)
    ),
    tuning = list(
      tau_sq = 0,
      car_1 = c(sigma_sq = 0, rho = 0),
      car_2 = c(sigma_sq = 0, rho = 0)
    ),
    priors = list(
      tau_sq = ig(2, 0.2),
      car_1 = list(sigma_sq = ig(2, 0.8), rho = uniform(0.05, 0.95)),
      car_2 = list(sigma_sq = ig(2, 0.6), rho = uniform(0.05, 0.95))
    ),
    n_samples = 4,
    verbose = FALSE
  )

  expect_equal(fit_same$term_description$global$n_graphs, 1L)
  expect_equal(fit_different$term_description$global$n_graphs, 2L)
})

test_that("CAR-time process uses location-major support and predicts existing nodes", {
  set.seed(304)
  n_area <- 5L
  ids <- paste0("a", seq_len(n_area))
  adj <- Matrix::sparseMatrix(
    i = c(seq_len(n_area - 1L), seq_len(n_area - 1L) + 1L),
    j = c(seq_len(n_area - 1L) + 1L, seq_len(n_area - 1L)),
    x = 1,
    dims = c(n_area, n_area),
    dimnames = list(ids, ids)
  )
  g <- car_graph(adj)
  dat <- expand.grid(time = 1:3, area = ids)
  dat <- dat[order(dat$area, dat$time), ]
  dat$x <- rnorm(nrow(dat))
  dat$y <- rnorm(nrow(dat))
  dat$y[c(4, 11)] <- NA_real_

  fit <- stLMM(
    y ~ x + car_time(area, time, graph = g),
    data = dat,
    starting = list(tau_sq = 0.4, car_time_1 = c(sigma_sq = 0.8, rho = 0.35, phi = 0.25)),
    tuning = list(tau_sq = 0, car_time_1 = c(sigma_sq = 0, rho = 0, phi = 0)),
    priors = list(
      tau_sq = ig(2, 0.4),
      car_time_1 = list(sigma_sq = ig(2, 0.8), rho = uniform(0.05, 0.95), phi = uniform(-0.9, 0.9))
    ),
    n_samples = 12,
    verbose = FALSE
  )

  expect_equal(fit$backend$graphs[[1]]$graph_type, "car_time")
  expect_equal(fit$backend$graphs[[1]]$n_space, n_area)
  expect_equal(fit$backend$graphs[[1]]$n_time, 3L)
  expect_equal(fit$backend$process_terms[[1]]$map, seq_len(nrow(dat)))
  expect_equal(colnames(fit$theta_samples), c("car_time_1_rho", "car_time_1_phi"))

  rec <- recover(fit, sub_sample = list(start = 2, thin = 2))
  expect_equal(ncol(rec$w_samples$car_time_1), n_area * 3L)
  expect_equal(length(fitted(rec)), nrow(dat))

  pred_full <- predict(rec)
  expect_equal(ncol(pred_full$mu_samples), nrow(dat))

  pred_existing <- predict(rec, newdata = dat[c(1, 7), c("area", "time", "x"), drop = FALSE])
  expect_equal(ncol(pred_existing$mu_samples), 2L)

  bad <- dat[1, c("area", "time", "x"), drop = FALSE]
  bad$time <- 99
  expect_error(predict(rec, newdata = bad), "new CAR-time value")
})

test_that("CAR-time exponential time model uses lambda and ordered time gaps", {
  set.seed(308)
  ids <- paste0("a", 1:4)
  adj <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 3, 3, 4),
    j = c(2, 1, 3, 2, 4, 3),
    x = 1,
    dims = c(4, 4),
    dimnames = list(ids, ids)
  )
  g <- car_graph(adj)
  time <- c(0, 1, 3, 6)
  dat <- expand.grid(time = time, area = ids)
  dat <- dat[order(dat$area, dat$time), ]
  dat$y <- rnorm(nrow(dat))

  lambda <- 0.45
  Q_time <- car_time_exp_precision(time, lambda)
  C_time <- exp(-lambda * abs(outer(time, time, "-")))
  expect_equal(Q_time, solve(C_time), tolerance = 1e-10)

  fit <- stLMM(
    y ~ 0 + car_time(area, time, graph = g, time_model = "exp"),
    data = dat,
    starting = list(tau_sq = 0.5, car_time_1 = c(sigma_sq = 0.9, rho = 0.4, lambda = lambda)),
    tuning = list(tau_sq = 0, car_time_1 = c(sigma_sq = 0, rho = 0, lambda = 0.04)),
    priors = list(
      tau_sq = ig(2, 0.5),
      car_time_1 = list(sigma_sq = ig(2, 0.9), rho = uniform(0.05, 0.95), lambda = uniform(0.05, 2.5))
    ),
    n_samples = 12,
    verbose = FALSE
  )

  expect_equal(fit$backend$graphs[[1]]$time_model, "exp")
  expect_equal(fit$backend$graphs[[1]]$time_delta, diff(time))
  expect_equal(fit$backend$process_terms[[1]]$params$time_model, "exp")
  expect_equal(colnames(fit$theta_samples), c("car_time_1_rho", "car_time_1_lambda"))
  expect_true(all(is.finite(fit$theta_samples[, "car_time_1_lambda"])))
  expect_true(all(fit$theta_samples[, "car_time_1_lambda"] > 0))

  rec <- recover(fit, sub_sample = list(start = 2, thin = 2))
  expect_equal(ncol(rec$w_samples$car_time_1), length(ids) * length(time))
})

test_that("CAR-time recovery covariance matches separable observed-row target", {
  set.seed(305)
  ids <- paste0("a", 1:4)
  adj <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 3, 3, 4),
    j = c(2, 1, 3, 2, 4, 3),
    x = 1,
    dims = c(4, 4),
    dimnames = list(ids, ids)
  )
  g <- car_graph(adj)
  support <- expand.grid(time = 1:3, area = ids)
  support <- support[order(support$area, support$time), ]
  dat <- support[rep(seq_len(nrow(support)), c(2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1)), ]
  dat$y <- 0
  dat$y[dat$area == "a2" & dat$time == 2] <- NA_real_
  dat$y[dat$area == "a4" & dat$time == 1] <- NA_real_

  sigma_sq <- 1
  rho <- 0.45
  phi <- 0.30
  tau_sq <- 0.8

  fit <- suppressWarnings(stLMM(
    y ~ 0 + car_time(area, time, graph = g),
    data = dat,
    starting = list(tau_sq = tau_sq, car_time_1 = c(sigma_sq = sigma_sq, rho = rho, phi = phi)),
    tuning = list(tau_sq = 0, car_time_1 = c(sigma_sq = 0, rho = 0, phi = 0)),
    priors = list(
      tau_sq = ig(2, tau_sq),
      car_time_1 = list(sigma_sq = ig(2, sigma_sq), rho = uniform(0.05, 0.95), phi = uniform(-0.9, 0.9))
    ),
    n_samples = 20000,
    verbose = FALSE
  ))
  rec <- recover(fit)

  adj_dense <- as.matrix(adj)
  degree <- rowSums(adj_dense != 0)
  Q_space <- diag(degree) - rho * adj_dense
  den <- 1 - phi^2
  Q_time <- matrix(0, 3, 3)
  Q_time[1, 1] <- 1 / den
  Q_time[3, 3] <- 1 / den
  Q_time[2, 2] <- (1 + phi^2) / den
  Q_time[1, 2] <- Q_time[2, 1] <- -phi / den
  Q_time[2, 3] <- Q_time[3, 2] <- -phi / den
  Q_prior <- kronecker(Q_space, Q_time) / sigma_sq

  node_count <- tabulate(fit$backend$process_terms_obs[[1]]$map, nbins = nrow(Q_prior))
  target_cov <- solve(Q_prior + diag(node_count / tau_sq, nrow(Q_prior)))
  recovered_cov <- cov(rec$w_samples$car_time_1)

  expect_lt(max(abs(recovered_cov - target_cov)), 0.035)
})

test_that("Leroux CAR-time recovery covariance matches separable observed-row target", {
  set.seed(311)
  ids <- paste0("a", 1:4)
  adj <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 3, 3, 4),
    j = c(2, 1, 3, 2, 4, 3),
    x = 1,
    dims = c(4, 4),
    dimnames = list(ids, ids)
  )
  g <- car_graph(adj)
  support <- expand.grid(time = 1:3, area = ids)
  support <- support[order(support$area, support$time), ]
  dat <- support[rep(seq_len(nrow(support)), c(2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1)), ]
  dat$y <- 0
  dat$y[dat$area == "a2" & dat$time == 2] <- NA_real_
  dat$y[dat$area == "a4" & dat$time == 1] <- NA_real_

  sigma_sq <- 1.1
  rho <- 0.5
  phi <- 0.25
  tau_sq <- 0.7

  fit <- suppressWarnings(stLMM(
    y ~ 0 + car_time(area, time, graph = g, car_model = "leroux"),
    data = dat,
    starting = list(tau_sq = tau_sq, car_time_1 = c(sigma_sq = sigma_sq, rho = rho, phi = phi)),
    tuning = list(tau_sq = 0, car_time_1 = c(sigma_sq = 0, rho = 0, phi = 0)),
    priors = list(
      tau_sq = ig(2, tau_sq),
      car_time_1 = list(sigma_sq = ig(2, sigma_sq), rho = uniform(0.01, 0.99), phi = uniform(-0.9, 0.9))
    ),
    n_samples = 20000,
    verbose = FALSE
  ))
  rec <- recover(fit)

  den <- 1 - phi^2
  Q_time <- matrix(0, 3, 3)
  Q_time[1, 1] <- 1 / den
  Q_time[3, 3] <- 1 / den
  Q_time[2, 2] <- (1 + phi^2) / den
  Q_time[1, 2] <- Q_time[2, 1] <- -phi / den
  Q_time[2, 3] <- Q_time[3, 2] <- -phi / den
  Q_prior <- kronecker(leroux_precision(adj, sigma_sq = 1, rho = rho), Q_time) / sigma_sq

  node_count <- tabulate(fit$backend$process_terms_obs[[1]]$map, nbins = nrow(Q_prior))
  target_cov <- solve(Q_prior + diag(node_count / tau_sq, nrow(Q_prior)))
  recovered_cov <- cov(rec$w_samples$car_time_1)

  expect_equal(fit$backend$graphs[[1]]$car_model, "leroux")
  expect_equal(fit$backend$process_terms[[1]]$params$car_model, "leroux")
  expect_lt(max(abs(recovered_cov - target_cov)), 0.035)
})

test_that("CAR-time terms reuse graph structure and SVC scaling with missing nodes", {
  set.seed(307)
  ids <- paste0("a", 1:4)
  adj <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 3, 3, 4),
    j = c(2, 1, 3, 2, 4, 3),
    x = 1,
    dims = c(4, 4),
    dimnames = list(ids, ids)
  )
  g <- car_graph(adj)

  dat <- expand.grid(time = 1:3, area = ids)
  dat <- dat[order(dat$area, dat$time), ]
  dat$x_svc <- seq(0.75, 1.75, length.out = nrow(dat))
  dat$y <- rnorm(nrow(dat))
  dat$y[c(3, 8)] <- NA_real_

  fit <- stLMM(
    y ~ 0 + car_time(area, time, graph = g) + x_svc:car_time(area, time, graph = g),
    data = dat,
    starting = list(
      tau_sq = 0.25,
      car_time_1 = c(sigma_sq = 0.8, rho = 0.35, phi = 0.25),
      car_time_2 = c(sigma_sq = 0.6, rho = 0.30, phi = 0.15)
    ),
    tuning = list(
      tau_sq = 0,
      car_time_1 = c(sigma_sq = 0, rho = 0, phi = 0),
      car_time_2 = c(sigma_sq = 0, rho = 0, phi = 0)
    ),
    priors = list(
      tau_sq = ig(2, 0.25),
      car_time_1 = list(sigma_sq = ig(2, 0.8), rho = uniform(0.05, 0.95), phi = uniform(-0.9, 0.9)),
      car_time_2 = list(sigma_sq = ig(2, 0.6), rho = uniform(0.05, 0.95), phi = uniform(-0.9, 0.9))
    ),
    n_samples = 12,
    verbose = FALSE
  )

  expect_length(fit$backend$graphs, 1L)
  expect_equal(unname(vapply(fit$backend$process_terms, `[[`, integer(1), "graph_index")),
               c(1L, 1L))
  expect_false(is.null(fit$backend$process_terms[[2]]$x))
  expect_equal(fit$backend$process_terms[[2]]$x, dat$x_svc)
  expect_equal(fit$backend$process_terms_obs[[2]]$x, dat$x_svc[!is.na(dat$y)])
  expect_equal(fit$backend$n_missing_response, 2L)

  rec <- recover(fit, sub_sample = list(start = 2, thin = 2))
  expect_equal(ncol(rec$w_samples$car_time_1), length(ids) * 3L)
  expect_equal(ncol(rec$w_samples$car_time_2), length(ids) * 3L)

  map <- fit$backend$process_terms[[1]]$map
  fitted_draws <- fitted(rec, summary = FALSE)
  manual <- rec$w_samples$car_time_1[, map, drop = FALSE] +
    sweep(rec$w_samples$car_time_2[, map, drop = FALSE], 2L, dat$x_svc, `*`)
  expect_equal(as.vector(fitted_draws), as.vector(manual),
               tolerance = 1e-10)

  newdata <- dat[c(1, 7), c("area", "time", "x_svc"), drop = FALSE]
  pred <- predict(rec, newdata = newdata)
  pred_map <- map[c(1, 7)]
  manual_pred <- rec$w_samples$car_time_1[, pred_map, drop = FALSE] +
    sweep(rec$w_samples$car_time_2[, pred_map, drop = FALSE], 2L, newdata$x_svc, `*`)
  expect_equal(as.vector(pred$mu_samples), as.vector(manual_pred),
               tolerance = 1e-10)
})

test_that("Leroux CAR works with binomial Polya-Gamma fitting path", {
  set.seed(312)
  ids <- paste0("a", 1:6)
  adj <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 3, 3, 4, 4, 5, 5, 6),
    j = c(2, 1, 3, 2, 4, 3, 5, 4, 6, 5),
    x = 1,
    dims = c(6, 6),
    dimnames = list(ids, ids)
  )
  g <- car_graph(adj)
  dat <- data.frame(
    y = rbinom(6, size = 10, prob = 0.45),
    n = 10L,
    area = ids,
    x = rnorm(6)
  )

  fit <- stLMM(
    y ~ x + car(area, graph = g, car_model = "leroux"),
    data = dat,
    family = "binomial",
    trials = "n",
    starting = list(car_1 = c(sigma_sq = 0.7, rho = 0.35)),
    tuning = list(car_1 = c(sigma_sq = 0, rho = 0)),
    priors = list(
      car_1 = list(sigma_sq = ig(2, 0.7), rho = uniform(0.01, 0.99))
    ),
    n_samples = 10,
    verbose = FALSE
  )

  expect_equal(fit$backend$family, "binomial")
  expect_equal(fit$backend$graphs[[1]]$car_model, "leroux")
  rec <- recover(fit, sub_sample = list(start = 2, thin = 2))
  expect_equal(ncol(rec$w_samples$car_1), length(ids))
})

test_that("CAR graph rejects isolated areas", {
  adj <- diag(0, 3)
  rownames(adj) <- colnames(adj) <- letters[1:3]
  expect_error(car_graph(adj), "isolated area")
})

test_that("CAR graph can bridge sf polygon islands with nearest edges", {
  testthat::skip_if_not_installed("sf")

  square <- function(x0, y0){
    sf::st_polygon(list(matrix(
      c(
        x0, y0,
        x0 + 1, y0,
        x0 + 1, y0 + 1,
        x0, y0 + 1,
        x0, y0
      ),
      ncol = 2,
      byrow = TRUE
    )))
  }

  poly <- sf::st_sf(
    id = c("a", "b", "c"),
    geometry = sf::st_sfc(
      square(0, 0),
      square(1, 0),
      square(10, 0),
      crs = 3857
    )
  )

  expect_error(car_graph(poly, id = "id"), "isolated area")

  g <- car_graph(poly, id = "id", island = "nearest", island_k = 1)
  expect_s3_class(g, "stLMM_car_graph")
  expect_equal(g$island_policy, "nearest")
  expect_equal(g$island_components_initial, 2L)
  expect_equal(nrow(g$island_added_edges), 1L)
  expect_true(all(g$degree > 0))
  expect_equal(max(stLMM:::car_connected_components(g$adjacency)), 1L)
})

test_that("bundled stunitco layer builds a bridged CAR graph", {
  testthat::skip_if_not_installed("sf")

  data(stunitco, package = "stLMM")

  expect_error(car_graph(stunitco, id = "COUNTYFIPS"), "isolated area")

  g <- car_graph(stunitco, id = "COUNTYFIPS", island = "nearest")
  expect_s3_class(g, "stLMM_car_graph")
  expect_equal(g$n, 3108L)
  expect_equal(g$island_policy, "nearest")
  expect_equal(g$island_components_initial, 3L)
  expect_equal(nrow(g$island_added_edges), 2L)
  expect_equal(g$island_added_edges$from, c("25019", "53055"))
  expect_equal(g$island_added_edges$to, c("25001", "53029"))
  expect_equal(max(stLMM:::car_connected_components(g$adjacency)), 1L)
})

test_that("CAR graph rejects disconnected non-isolated components by default", {
  ids <- paste0("a", 1:4)
  adj <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 4),
    j = c(2, 1, 4, 3),
    x = 1,
    dims = c(4, 4),
    dimnames = list(ids, ids)
  )

  expect_error(car_graph(adj), "disconnected components")
  expect_error(car_graph(adj, island = "nearest"), "requires sf polygon input")
})

test_that("CAR rho updates exercise sparse prior determinant path", {
  set.seed(302)
  n_area <- 6L
  ids <- paste0("a", seq_len(n_area))
  adj <- Matrix::sparseMatrix(
    i = c(seq_len(n_area - 1L), seq_len(n_area - 1L) + 1L),
    j = c(seq_len(n_area - 1L) + 1L, seq_len(n_area - 1L)),
    x = 1,
    dims = c(n_area, n_area),
    dimnames = list(ids, ids)
  )
  g <- car_graph(adj)
  dat <- data.frame(
    y = rnorm(n_area),
    area = ids
  )

  fit <- stLMM(
    y ~ 0 + car(area, graph = g),
    data = dat,
    starting = list(tau_sq = 0.2, car_1 = c(sigma_sq = 0.8, rho = 0.4)),
    tuning = list(tau_sq = 0, car_1 = c(sigma_sq = 0, rho = 0.05)),
    priors = list(
      tau_sq = ig(2, 0.2),
      car_1 = list(sigma_sq = ig(2, 0.8), rho = uniform(0.05, 0.95))
    ),
    n_samples = 8,
    verbose = FALSE
  )

  expect_equal(ncol(fit$theta_samples), 1L)
  expect_true(all(is.finite(fit$theta_samples[, "car_1_rho"])))
  expect_true(is.finite(fit$covariance_acceptance) || is.na(fit$covariance_acceptance))
})
