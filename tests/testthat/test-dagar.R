dagar_precision <- function(graph, sigma_sq, rho){
  n <- graph$n
  Q <- matrix(0, n, n)
  for(i in seq_len(n)){
    parents <- if(graph$parent_count[i] > 0L) {
      graph$parent_index[seq.int(
        graph$parent_start[i] + 1L,
        graph$parent_start[i] + graph$parent_count[i]
      )] + 1L
    } else {
      integer(0)
    }
    m <- length(parents)
    denom <- 1 + (m - 1) * rho^2
    b <- if(m > 0L) rho / denom else 0
    f <- denom / (1 - rho^2)
    l <- numeric(n)
    l[i] <- 1
    if(m > 0L)
      l[parents] <- -b
    Q <- Q + f * tcrossprod(l)
  }
  Q / sigma_sq
}

simulate_from_precision <- function(Q){
  L <- chol(Q)
  z <- rnorm(nrow(Q))
  as.numeric(backsolve(L, z))
}

dagar_time_ar1_precision <- function(phi, n_time){
  Q <- matrix(0, n_time, n_time)
  if(n_time == 1L){
    Q[1, 1] <- 1
    return(Q)
  }
  den <- 1 - phi^2
  diag(Q) <- c(1, rep(1 + phi^2, n_time - 2L), 1) / den
  off <- -phi / den
  Q[cbind(seq_len(n_time - 1L), seq_len(n_time - 1L) + 1L)] <- off
  Q[cbind(seq_len(n_time - 1L) + 1L, seq_len(n_time - 1L))] <- off
  Q
}

dagar_time_exp_precision <- function(time, lambda){
  time <- sort(unique(time))
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
  off <- -phi / den
  Q[cbind(seq_len(n_time - 1L), seq_len(n_time - 1L) + 1L)] <- off
  Q[cbind(seq_len(n_time - 1L) + 1L, seq_len(n_time - 1L))] <- off
  Q
}

dagar_time_precision <- function(graph, sigma_sq, rho, Q_time){
  graph_space <- graph
  graph_space$n <- graph$n_space
  graph_space$q <- graph$n_space
  Q_space <- dagar_precision(
    graph_space,
    sigma_sq = 1,
    rho = rho
  )
  kronecker(Q_space, Q_time) / sigma_sq
}

test_that("make_dagar_graph orients adjacency by ordering", {
  ids <- paste0("a", 1:5)
  adj <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 3, 3, 4, 4, 5),
    j = c(2, 1, 3, 2, 4, 3, 5, 4),
    x = 1,
    dims = c(5, 5),
    dimnames = list(ids, ids)
  )
  g <- car_graph(adj)
  dg <- stLMM:::make_dagar_graph(g, ordering = c(3, 1, 5, 2, 4))

  expect_equal(dg$graph_type, "dagar")
  expect_equal(dg$ids, ids[c(3, 1, 5, 2, 4)])
  expect_equal(dg$ordering, "user")
  expect_true(all(dg$parent_count >= 0L))
  for(i in seq_len(dg$n)){
    if(dg$parent_count[i] > 0L){
      idx <- dg$parent_index[seq.int(dg$parent_start[i] + 1L, dg$parent_start[i] + dg$parent_count[i])] + 1L
      expect_true(all(idx < i))
    }
  }
  expect_equal(sum(dg$parent_count), 4L)
})

test_that("DAGAR precision and log determinant match dense construction", {
  ids <- paste0("a", 1:6)
  adj <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 3, 2, 5, 3, 4, 4, 5, 5, 6),
    j = c(2, 1, 3, 2, 5, 2, 4, 3, 5, 4, 6, 5),
    x = 1,
    dims = c(6, 6),
    dimnames = list(ids, ids)
  )
  g <- car_graph(adj)
  dg <- stLMM:::make_dagar_graph(g, ordering = c(2, 5, 1, 3, 6, 4))
  sigma_sq <- 0.8
  rho <- 0.55
  tau_sq <- 0.2
  Q <- dagar_precision(dg, sigma_sq = sigma_sq, rho = rho)

  dat <- data.frame(y = rep(0, length(ids)), area = ids[dg$ord])
  fit <- stLMM(
    y ~ 0 + dagar(area, graph = g, ordering = c(2, 5, 1, 3, 6, 4)),
    data = dat,
    starting = list(tau_sq = fixed(tau_sq), dagar_1 = list(sigma_sq = fixed(sigma_sq), rho = fixed(rho))),
    n_samples = 5,
    verbose = FALSE
  )

  logdet_dense <- as.numeric(determinant(Q, logarithm = TRUE)$modulus)
  logdet_backend <- -dg$n * log(sigma_sq) +
    sum(log(1 + (dg$parent_count - 1) * rho^2) - log(1 - rho^2))
  expect_equal(logdet_backend, logdet_dense, tolerance = 1e-10)
  expect_equal(fit$backend$graphs[[1]]$parent_count, dg$parent_count)
  expect_equal(fit$term_description$process_terms$dagar_1$diagnostics$zero_parent_nodes, sum(dg$parent_count == 0L))
  expect_equal(fit$theta_samples[, "dagar_1_rho"], rep(rho, 5))
})

test_that("DAGAR recovery covariance matches analytic observed-row target", {
  set.seed(812)
  ids <- paste0("a", 1:7)
  adj <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 3, 3, 4, 4, 5, 3, 6, 6, 7),
    j = c(2, 1, 3, 2, 4, 3, 5, 4, 6, 3, 7, 6),
    x = 1,
    dims = c(7, 7),
    dimnames = list(ids, ids)
  )
  g <- car_graph(adj)
  ordering <- c(4, 3, 5, 2, 6, 1, 7)
  dg <- stLMM:::make_dagar_graph(g, ordering = ordering)

  sigma_sq <- 0.9
  rho <- 0.45
  tau_sq <- 0.35
  observed_count <- c(2L, 1L, 0L, 3L, 1L, 0L, 2L)
  rows_per_area <- pmax(observed_count, 1L)
  area <- rep(ids, rows_per_area)
  y <- numeric(length(area))
  y[area %in% ids[observed_count == 0L]] <- NA
  dat <- data.frame(y = y, area = area)

  fit <- stLMM(
    y ~ 0 + dagar(area, graph = g, ordering = ordering),
    data = dat,
    starting = list(tau_sq = fixed(tau_sq), dagar_1 = list(sigma_sq = fixed(sigma_sq), rho = fixed(rho))),
    n_samples = 700,
    verbose = FALSE
  )
  rec <- recover(fit, sub_sample = list(start = 101, thin = 2))

  Q_prior <- dagar_precision(dg, sigma_sq = sigma_sq, rho = rho)
  obs_count_ordered <- observed_count[dg$ord]
  Q_post <- Q_prior + diag(obs_count_ordered / tau_sq)
  target_cov <- solve(Q_post)
  recovered_cov <- cov(rec$w_samples_ordered$dagar_1)

  expect_equal(unname(diag(recovered_cov)), unname(diag(target_cov)), tolerance = 0.08)
  expect_equal(unname(recovered_cov[lower.tri(recovered_cov)]), unname(target_cov[lower.tri(target_cov)]), tolerance = 0.08)
  expect_equal(ncol(rec$w_samples$dagar_1), length(ids))
  expect_equal(colnames(rec$w_samples$dagar_1), paste0("dagar_1_", seq_along(ids)))
})

test_that("DAGAR fits and predicts existing graph areas", {
  set.seed(813)
  ids <- paste0("a", 1:8)
  adj <- Matrix::sparseMatrix(
    i = c(seq_len(7), seq_len(7) + 1L),
    j = c(seq_len(7) + 1L, seq_len(7)),
    x = 1,
    dims = c(8, 8),
    dimnames = list(ids, ids)
  )
  g <- car_graph(adj)
  dg <- stLMM:::make_dagar_graph(g, ordering = "default")
  sigma_sq <- 0.6
  rho <- 0.5
  tau_sq <- 0.08
  w <- simulate_from_precision(dagar_precision(dg, sigma_sq, rho))

  dat <- data.frame(area = rep(ids, each = 3L), x = rnorm(24))
  dat$y <- 0.3 + 0.4 * dat$x + rep(w, each = 3L) + rnorm(nrow(dat), sd = sqrt(tau_sq))
  dat$y[which(dat$area %in% ids[c(2, 7)])[1]] <- NA

  fit <- stLMM(
    y ~ x + dagar(area, graph = g),
    data = dat,
    starting = list(tau_sq = tau_sq, dagar_1 = c(sigma_sq = sigma_sq, rho = rho)),
    tuning = list(tau_sq = 0, dagar_1 = c(sigma_sq = 0, rho = 0)),
    priors = list(
      tau_sq = ig(2, tau_sq),
      dagar_1 = list(sigma_sq = ig(2, sigma_sq), rho = uniform(0.05, 0.95))
    ),
    n_samples = 20,
    verbose = FALSE
  )

  expect_equal(fit$backend$graphs[[1]]$graph_type, "dagar")
  expect_equal(colnames(fit$theta_samples), "dagar_1_rho")
  rec <- recover(fit, sub_sample = list(start = 5, thin = 3))
  pred <- predict(rec, newdata = dat[1:3, c("area", "x"), drop = FALSE])
  expect_equal(ncol(pred$mu_samples), 3L)

  bad <- dat[1, c("area", "x"), drop = FALSE]
  bad$area <- "missing_area"
  expect_error(predict(rec, newdata = bad), "new DAGAR area value")
})

test_that("DAGAR-time precision and recovery covariance match analytic target", {
  set.seed(814)
  ids <- paste0("a", 1:6)
  adj <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 3, 2, 5, 3, 4, 4, 5, 5, 6),
    j = c(2, 1, 3, 2, 5, 2, 4, 3, 5, 4, 6, 5),
    x = 1,
    dims = c(6, 6),
    dimnames = list(ids, ids)
  )
  g <- car_graph(adj)
  ordering <- c(2, 5, 1, 3, 6, 4)
  time <- 1:3
  dg <- stLMM:::make_dagar_time_graph(g, time, ordering = ordering)

  sigma_sq <- 0.7
  rho <- 0.45
  phi <- 0.35
  tau_sq <- 0.25
  Q_prior <- dagar_time_precision(dg, sigma_sq, rho, dagar_time_ar1_precision(phi, length(time)))

  dat <- expand.grid(area = ids, time = time)
  dat$y <- 0
  dat$y[c(3, 11)] <- NA

  fit <- stLMM(
    y ~ 0 + dagar_time(area, time, graph = g, ordering = ordering),
    data = dat,
    starting = list(
      tau_sq = fixed(tau_sq),
      dagar_time_1 = list(sigma_sq = fixed(sigma_sq), rho = fixed(rho), phi = fixed(phi))
    ),
    n_samples = 900,
    verbose = FALSE
  )
  rec <- recover(fit, sub_sample = list(start = 101, thin = 2))

  obs_count <- as.integer(!is.na(dat$y))
  map_support <- match(dat$area, ids)
  map_ordered <- (dg$ord_inv[map_support] - 1L) * length(time) + match(dat$time, time)
  obs_count_ordered <- tabulate(map_ordered[!is.na(dat$y)], nbins = nrow(Q_prior))
  Q_post <- Q_prior + diag(obs_count_ordered / tau_sq)
  target_cov <- solve(Q_post)
  recovered_cov <- cov(rec$w_samples_ordered$dagar_time_1)

  expect_equal(fit$backend$graphs[[1]]$graph_type, "dagar_time")
  expect_equal(colnames(fit$theta_samples), c("dagar_time_1_rho", "dagar_time_1_phi"))
  expect_equal(unname(diag(recovered_cov)), unname(diag(target_cov)), tolerance = 0.08)
  expect_equal(unname(recovered_cov[lower.tri(recovered_cov)]), unname(target_cov[lower.tri(target_cov)]), tolerance = 0.08)
  expect_equal(ncol(rec$w_samples$dagar_time_1), length(ids) * length(time))
})

test_that("DAGAR-time supports uneven exponential time and existing-support prediction", {
  set.seed(815)
  ids <- paste0("a", 1:5)
  adj <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 3, 3, 4, 4, 5),
    j = c(2, 1, 3, 2, 4, 3, 5, 4),
    x = 1,
    dims = c(5, 5),
    dimnames = list(ids, ids)
  )
  g <- car_graph(adj)
  time <- c(2000, 2001.5, 2004)
  sigma_sq <- 0.6
  rho <- 0.4
  lambda <- 0.7
  tau_sq <- 0.12

  dat <- expand.grid(area = ids, time = time)
  dat$x <- rnorm(nrow(dat))
  dat$y <- 0.2 + 0.5 * dat$x + rnorm(nrow(dat), sd = sqrt(tau_sq))

  fit <- stLMM(
    y ~ x + dagar_time(area, time, graph = g, time_model = "exp", ordering = "default"),
    data = dat,
    starting = list(
      tau_sq = fixed(tau_sq),
      dagar_time_1 = list(sigma_sq = fixed(sigma_sq), rho = fixed(rho), lambda = fixed(lambda))
    ),
    n_samples = 20,
    verbose = FALSE
  )

  expect_equal(fit$backend$graphs[[1]]$time_model, "exp")
  expect_equal(colnames(fit$theta_samples), c("dagar_time_1_rho", "dagar_time_1_lambda"))

  rec <- recover(fit, sub_sample = list(start = 5, thin = 3))
  pred <- predict(rec, newdata = dat[1:4, c("area", "time", "x"), drop = FALSE])
  expect_equal(ncol(pred$mu_samples), 4L)

  bad <- dat[1, c("area", "time", "x"), drop = FALSE]
  bad$time <- 2002
  expect_error(predict(rec, newdata = bad), "fitted time support")
})

test_that("CAR and DAGAR graph plot methods return ggplot objects", {
  testthat::skip_if_not_installed("ggplot2")

  ids <- paste0("a", 1:5)
  adj <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 3, 3, 4, 4, 5),
    j = c(2, 1, 3, 2, 4, 3, 5, 4),
    x = 1,
    dims = c(5, 5),
    dimnames = list(ids, ids)
  )
  g <- car_graph(adj)
  dg <- stLMM:::make_dagar_graph(g, ordering = c(3, 1, 5, 2, 4))
  dtg <- stLMM:::make_dagar_time_graph(g, time_support = 1:3, ordering = c(3, 1, 5, 2, 4))

  expect_s3_class(plot(g, show_ids = TRUE, color_by_degree = TRUE), "ggplot")
  expect_s3_class(plot(dg, show_ids = TRUE, color_by_degree = TRUE), "ggplot")
  expect_s3_class(plot(dtg, show_ids = TRUE, color_by_degree = TRUE), "ggplot")
})
