test_that("collapsed AR1 process recovery uses the Cholesky permutation", {
  set.seed(124)
  n <- 5L
  phi <- 0.4
  tau_sq <- 0.7
  dat <- data.frame(y = rep(0, n), time = seq_len(n))

  fit <- stLMM(
    y ~ 0 + ar1(time),
    data = dat,
    starting = list(tau_sq = tau_sq, ar1_1 = c(sigma_sq = 1, phi = phi)),
    tuning = list(tau_sq = 0, ar1_1 = c(sigma_sq = 0, phi = 0)),
    priors = list(
      tau_sq = ig(2, tau_sq),
      ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.9, 0.9))
    ),
    n_samples = 20000,
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

  target_cov <- solve(Q + diag(1 / tau_sq, n))
  recovered_cov <- cov(rec$w_samples$ar1_1)

  expect_lt(max(abs(recovered_cov - target_cov)), 0.025)
})

test_that("collapsed one-node AR1 recovery is phi independent", {
  set.seed(126)
  n_obs <- 4L
  sigma_sq <- 1
  tau_sq <- 0.6
  dat <- data.frame(y = rep(0, n_obs), time = rep(1, n_obs))

  fit <- suppressWarnings(stLMM(
    y ~ 0 + ar1(time),
    data = dat,
    starting = list(tau_sq = tau_sq, ar1_1 = c(sigma_sq = sigma_sq, phi = 0.8)),
    tuning = list(tau_sq = 0, ar1_1 = c(sigma_sq = 0, phi = 0)),
    priors = list(
      tau_sq = ig(2, 1),
      ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.9, 0.9))
    ),
    n_samples = 20000,
    verbose = FALSE
  ))
  rec <- recover(fit)

  target_var <- 1 / (1 / sigma_sq + n_obs / tau_sq)
  recovered_var <- var(as.numeric(rec$w_samples$ar1_1[, 1]))

  expect_lt(abs(recovered_var - target_var), 0.01)
})

test_that("collapsed NNGP process recovery matches the NNGP posterior covariance", {
  set.seed(125)
  coords <- cbind(
    lon = c(0, 1, 2, 0, 1, 2, 0.5, 1.5),
    lat = c(0, 0, 0, 1, 1, 1, 2, 2)
  )
  n <- nrow(coords)
  phi <- 0.7
  tau_sq <- 0.8
  dat <- data.frame(y = rep(0, n), lon = coords[, "lon"], lat = coords[, "lat"])

  fit <- suppressWarnings(stLMM(
    y ~ 0 + nngp(lon, lat, m = 3),
    data = dat,
    starting = list(tau_sq = tau_sq, nngp_1 = c(sigma_sq = 1, phi = phi)),
    tuning = list(tau_sq = 0, nngp_1 = c(sigma_sq = 0, phi = 0)),
    priors = list(
      tau_sq = ig(2, tau_sq),
      nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5))
    ),
    n_samples = 20000,
    verbose = FALSE
  ))
  rec <- recover(fit)

  graph <- fit$backend$graphs[[1]]
  coords_ord <- graph$coords_ord
  nn_indx <- graph$nnIndx
  nn_lu <- graph$nnIndxLU
  Q <- matrix(0, n, n)

  for(i0 in 0:(n - 1)){
    start <- nn_lu[i0 + 1L]
    m <- nn_lu[n + i0 + 1L]

    if(m == 0L){
      nodes <- i0 + 1L
      weights <- 1
      F <- 1
    } else {
      nbrs <- nn_indx[(start + 1L):(start + m)] + 1L
      dist_i_nbrs <- sqrt(rowSums(
        (coords_ord[rep(i0 + 1L, m), , drop = FALSE] -
           coords_ord[nbrs, , drop = FALSE])^2
      ))
      c_i_nbrs <- exp(-phi * dist_i_nbrs)
      C_nbrs <- matrix(1, m, m)

      for(a in seq_len(m)){
        for(b in seq_len(m)){
          C_nbrs[a, b] <- exp(-phi * sqrt(sum(
            (coords_ord[nbrs[a], ] - coords_ord[nbrs[b], ])^2
          )))
        }
      }

      B <- as.vector(solve(C_nbrs, c_i_nbrs))
      F <- 1 / (1 - sum(B * c_i_nbrs))
      nodes <- c(i0 + 1L, nbrs)
      weights <- c(1, -B)
    }

    Q[nodes, nodes] <- Q[nodes, nodes] + F * tcrossprod(weights)
  }

  target_cov <- solve(Q + diag(1 / tau_sq, n))
  recovered_cov <- cov(rec$w_samples_ordered$nngp_1)
  recovered_cov_user <- cov(rec$w_samples$nngp_1)

  expect_lt(max(abs(recovered_cov - target_cov)), 0.025)
  expect_equal(unname(recovered_cov_user), unname(recovered_cov[graph$ord_inv, graph$ord_inv]))

  set.seed(1251)
  rec_threads_1 <- recover(fit, n_omp_threads = 1, sub_sample = list(start = 1, thin = 20))
  set.seed(1251)
  rec_threads_2 <- recover(fit, n_omp_threads = 2, sub_sample = list(start = 1, thin = 20))
  expect_equal(
    unname(rec_threads_1$w_samples_ordered$nngp_1),
    unname(rec_threads_2$w_samples_ordered$nngp_1)
  )
  expect_error(recover(fit, n_omp_threads = 0), "n_omp_threads must be a positive integer")
  progress_messages <- character()
  withCallingHandlers(
    recover(fit, n_omp_threads = 1, verbose = TRUE, sub_sample = list(start = 1, thin = 20)),
    message = function(m){
      progress_messages <<- c(progress_messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_true(any(grepl("recover: sampling latent process draws", progress_messages)))
})

test_that("collapsed NNGP recovery handles repeated observations and orderings", {
  n_unique <- 8L
  coords_unique <- cbind(
    lon = c(0, 1, 2, 0, 1, 2, 0.5, 1.5),
    lat = c(0, 0, 0, 1, 1, 1, 2, 2)
  )
  rep_count <- c(3L, 1L, 2L, 4L, 1L, 3L, 2L, 1L)
  phi <- 0.7
  tau_sq <- 0.75
  orderings <- list(
    coord = "coord",
    maxmin = "maxmin",
    hilbert = "hilbert",
    random = "random",
    user = c(4L, 1L, 7L, 2L, 5L, 8L, 3L, 6L)
  )

  for(ord_name in names(orderings)){
    set.seed(130 + match(ord_name, names(orderings)))
    idx <- rep(seq_len(n_unique), rep_count)
    dat <- data.frame(
      y = rep(0, length(idx)),
      lon = coords_unique[idx, "lon"],
      lat = coords_unique[idx, "lat"]
    )
    form <- if(is.numeric(orderings[[ord_name]])){
      as.formula(paste0(
        "y ~ 0 + nngp(lon, lat, m = 3, ordering = c(",
        paste(as.integer(orderings[[ord_name]]), collapse = ","),
        "))"
      ))
    } else {
      as.formula(paste0(
        "y ~ 0 + nngp(lon, lat, m = 3, ordering = \"",
        orderings[[ord_name]],
        "\")"
      ))
    }

    fit <- suppressWarnings(stLMM(
      form,
      data = dat,
      starting = list(tau_sq = tau_sq, nngp_1 = c(sigma_sq = 1, phi = phi)),
      tuning = list(tau_sq = 0, nngp_1 = c(sigma_sq = 0, phi = 0)),
      priors = list(
        tau_sq = ig(2, tau_sq),
        nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5))
      ),
      n_samples = 12000,
      verbose = FALSE
    ))
    rec <- recover(fit)

    graph <- fit$backend$graphs[[1]]
    coords_ord <- graph$coords_ord
    nn_indx <- graph$nnIndx
    nn_lu <- graph$nnIndxLU
    Q <- matrix(0, n_unique, n_unique)

    for(i0 in 0:(n_unique - 1)){
      start <- nn_lu[i0 + 1L]
      m <- nn_lu[n_unique + i0 + 1L]

      if(m == 0L){
        nodes <- i0 + 1L
        weights <- 1
        F <- 1
      } else {
        nbrs <- nn_indx[(start + 1L):(start + m)] + 1L
        dist_i_nbrs <- sqrt(rowSums(
          (coords_ord[rep(i0 + 1L, m), , drop = FALSE] -
             coords_ord[nbrs, , drop = FALSE])^2
        ))
        c_i_nbrs <- exp(-phi * dist_i_nbrs)
        C_nbrs <- matrix(1, m, m)

        for(a in seq_len(m)){
          for(b in seq_len(m)){
            C_nbrs[a, b] <- exp(-phi * sqrt(sum(
              (coords_ord[nbrs[a], ] - coords_ord[nbrs[b], ])^2
            )))
          }
        }

        B <- as.vector(solve(C_nbrs, c_i_nbrs))
        F <- 1 / (1 - sum(B * c_i_nbrs))
        nodes <- c(i0 + 1L, nbrs)
        weights <- c(1, -B)
      }

      Q[nodes, nodes] <- Q[nodes, nodes] + F * tcrossprod(weights)
    }

    node_nobs <- tabulate(fit$backend$process_terms[[1]]$map, nbins = n_unique)
    target_cov <- solve(Q + diag(node_nobs / tau_sq, n_unique))
    recovered_cov <- cov(rec$w_samples_ordered$nngp_1)
    recovered_cov_user <- cov(rec$w_samples$nngp_1)

    expect_lt(max(abs(recovered_cov - target_cov)), 0.04)
    expect_equal(unname(recovered_cov_user), unname(recovered_cov[graph$ord_inv, graph$ord_inv]))
  }
})
