joint_recovery_fit <- function(process = FALSE, chains = 1L){
  set.seed(91701)
  d <- data.frame(y = rnorm(32), x = rnorm(32), group = rep(letters[1:8], each = 4),
                  parent = rep(c("A", "B"), each = 16), time = rep(1:16, each = 2),
                  off = rep(c(-.2, .3), 16))
  d$y[c(3, 12)] <- NA_real_
  f <- y ~ x + offset(off) + iid(group, ~ 1 + x, by = parent) + iid(parent)
  if(process) f <- update(f, . ~ . + ar1(time))
  stLMM(f, d, priors = list(iid_1 = list(covariance = iw(4, diag(2))),
    iid_2 = list(sigma_sq = ig(3, 1))),
    starting = list(ar1_1 = list(sigma_sq = fixed(.4), phi = fixed(.5))),
    collapse_iid = TRUE, n_samples = 60, chains = chains,
    chain_control = list(seed = 91702), warmup = FALSE, verbose = FALSE)
}

test_that("collapsed Gaussian effects are recovered jointly on a selected grid", {
  for(process in c(FALSE, TRUE)){
    if(process) expect_warning(f <- joint_recovery_fit(process), "duplicated time values")
    else f <- joint_recovery_fit(process)
    expect_null(f$alpha_samples)
    expect_null(f$samples$alpha)
    expect_equal(nrow(as_samples(f)), 60L)
    expect_equal(summary(f)$n_used, 60L)
    expect_error(predict(f), "require recover")
    expect_error(fitted(f), "require recover")
    expect_error(log_lik(f), "recover")
    expect_error(summary(f, type = "combined", terms = "iid_1"), "require recover")
    expect_error(recover(f, sub_sample = list(start = 61)), "no posterior draws")

    set.seed(91703)
    r <- recover(f, sub_sample = list(start = 5, thin = 3))
    idx <- seq.int(5L, 60L, 3L)
    expect_identical(r$recover_iter, idx)
    expect_equal(dim(r$alpha_samples), c(length(idx), f$backend$q))
    expect_identical(colnames(r$alpha_samples), colnames(f$backend$Z))
    expect_identical(r$beta_samples, f$beta_samples)
    expect_identical(r$iid_covariance_samples, f$iid_covariance_samples)
    eta <- r$beta_samples[idx, , drop = FALSE] %*% t(f$backend$X) +
      r$alpha_samples %*% t(as.matrix(f$backend$Z))
    eta <- sweep(eta, 2, f$backend$offset, "+")
    if(process) eta <- eta + r$w_samples$ar1_1[, f$backend$process_terms[[1]]$map]
    expect_equal(fitted(r, summary = FALSE), eta, ignore_attr = TRUE)
    expect_equal(predict(r)$mu_samples, eta, ignore_attr = TRUE)
    expect_identical(attr(log_lik(r), "draw_index"), idx)
    # Parameter selection uses original iteration IDs; alpha uses recovery rows.
    a <- as_samples(r, burn = 10, thin = 2)
    selected <- idx[idx > 10][seq(1, sum(idx > 10), 2)]
    expect_identical(a$.iteration, selected)
    expect_equal(as.matrix(a[, colnames(r$alpha_samples)]),
      r$alpha_samples[match(selected, idx), ], ignore_attr = TRUE)
    expect_equal(as_mcmc(r, burn = 10, thin = 2),
      as.matrix(a[, -(1:2)]), ignore_attr = TRUE)
    s <- summary(r, burn = 10)
    expect_equal(s$n_used, sum(idx > 10))
    expect_equal(s$parameters$alpha[, "mean"], colMeans(r$alpha_samples[idx > 10, ]))
    s <- summary(r, type = "combined", terms = "iid_1", burn = 10)
    expect_identical(s$iteration, idx[idx > 10])
    expect_equal(s$n_used, sum(idx > 10))
    # A second recovery call can select different original draws.
    r2 <- recover(r, sub_sample = list(start = 7, thin = 5))
    expect_identical(r2$recover_iter, seq.int(7L, 60L, 5L))
    expect_equal(nrow(fitted(r2, summary = FALSE)), length(r2$recover_iter))
    bad <- r; bad$recover_iter <- bad$recover_iter[-1]
    expect_error(fitted(bad), "do not align")
  }
})

test_that("recovered IID-only models support new levels and no fixed effects", {
  set.seed(91704)
  d <- data.frame(y = rnorm(12), g = factor(rep(c("a", "b", "c"), each = 4),
                                          levels = c("a", "b", "c", "unused")))
  f <- stLMM(y ~ 0 + iid(g), d, collapse_iid = TRUE, n_samples = 40,
    priors = list(iid_1 = list(sigma_sq = ig(3, 1))), warmup = FALSE, verbose = FALSE)
  r <- recover(f, sub_sample = list(start = 4, thin = 3))
  expect_equal(dim(r$alpha_samples), c(13L, 4L))
  p <- predict(r, data.frame(g = c("a", "a", "unused")))
  expect_equal(p$mu_samples[, 1], p$mu_samples[, 2])
  expect_equal(p$mu_samples[, 3], r$alpha_samples[, 4])
  expect_error(predict(r, data.frame(g = "new")), "allow_new_levels")
  p <- predict(r, data.frame(g = c("new", "new")), allow_new_levels = TRUE)
  expect_equal(p$mu_samples[, 1], p$mu_samples[, 2])
  expect_true(all(is.finite(p$mu_samples)))
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  expect_invisible(plot(r, parameters = colnames(r$alpha_samples)[1], burnin = 10))
  expect_invisible(plot(r, max_parameters = 1, burnin = 3, thin = 2))
  expect_invisible(plot(r, type = "density", max_parameters = 1, burnin = 3, thin = 2))
  expect_invisible(plot(r, type = "fitted"))
  # All single-chain recovery plots count recovered rows for burnin.
  for(type in c("trace", "density", "fitted"))
    expect_error(plot(r, type = type, burnin = nrow(r$alpha_samples)), "burnin removes all")
  grDevices::dev.off()
  expect_output(print(summary(r)), "40 \\(13 used on recovery grid\\)")
  r3 <- recover(f, sub_sample = list(start = 3, thin = 3))
  expect_error(summary(r3, type = "combined", terms = "all", burn = 39), "removes all")
  f$iid_sigma_sq_samples[1, 1] <- 0
  expect_error(recover(f), "finite and positive")
})

test_that("joint recovery checks covariance definiteness and chain alignment", {
  f <- joint_recovery_fit(chains = 2L)
  r <- recover(f, sub_sample = list(start = 5, thin = 3))
  expect_equal(summary(r)$n_used, 19L)
  expect_output(print(summary(r, diagnostics = FALSE)),
                "60 \\(19 used on recovery grid\\)")
  set.seed(91705)
  r2 <- recover(r, sub_sample = list(start = 7, thin = 5))
  rng <- .Random.seed
  set.seed(91705)
  fresh <- recover(f, sub_sample = list(start = 7, thin = 5))
  expect_identical(.Random.seed, rng)
  for(k in seq_along(r2$chains)){
    expect_identical(r2$chains[[k]]$recover_iter, seq.int(7L, 60L, 5L))
    expect_identical(r2$chains[[k]]$alpha_samples, fresh$chains[[k]]$alpha_samples)
    expect_identical(r2$chains[[k]]$beta_samples, f$chains[[k]]$beta_samples)
  }
  s <- summary(r, type = "combined", terms = "iid_1", burn = 10)
  expect_identical(s$iteration, seq.int(11L, 59L, 3L))
  r$chains[[2]] <- recover(f$chains[[2]], sub_sample = list(start = 6, thin = 3))
  expect_error(summary(r, type = "combined", terms = "iid_1"), "matching recovered iterations")
  expect_error(summary(r), "matching recovered iterations")
  bad <- f$chains[[1]]
  bad$iid_covariance_samples[1, 1:3] <- c(1, 2, 1)
  expect_error(recover(bad), "Cholesky failed")
  bad$iid_covariance_samples[1, 1] <- NA_real_
  expect_error(recover(bad), "must be finite")
})

test_that("failed recovery saves RNG progress and permits a clean retry", {
  f <- joint_recovery_fit()
  bad <- f
  bad$iid_covariance_samples <- bad$iid_covariance_samples[, -1, drop = FALSE]
  set.seed(91706)
  rng <- .Random.seed
  expect_error(recover(bad), "incompatible dimensions")
  expect_identical(.Random.seed, rng)
  set.seed(91706)
  first <- recover(f, sub_sample = list(start = 1, thin = 60))
  rng_after_first <- .Random.seed
  # Fail after one completed draw, through both explicit checks and LAPACK.
  for(failure in c("variance", "nonfinite", "nonpositive_definite")){
    bad <- f
    if(failure == "variance"){
      bad$iid_sigma_sq_samples[2, 1] <- 0
      message <- "finite and positive"
    } else if(failure == "nonfinite"){
      bad$iid_covariance_samples[2, 1] <- NA_real_
      message <- "must be finite"
    } else {
      bad$iid_covariance_samples[2, 1:3] <- c(1, 2, 1)
      message <- "Cholesky failed"
    }
    set.seed(91706)
    expect_error(recover(bad), message)
    expect_identical(.Random.seed, rng_after_first)
    set.seed(91706)
    retry <- recover(f, sub_sample = list(start = 1, thin = 60))
    expect_identical(retry$alpha_samples, first$alpha_samples)
    expect_identical(.Random.seed, rng_after_first)
  }
})

test_that("repeated augmented multi-chain recovery only subsets saved draws", {
  d <- data.frame(y = rep(c(0, 1), 12), g = rep(letters[1:6], each = 4), t = 1:24)
  for(collapse in c(FALSE, TRUE)){
    f <- stLMM(y ~ iid(g) + ar1(t), d, family = "binomial", collapse_iid = collapse,
      priors = list(iid_1 = list(sigma_sq = ig(3, 1))),
      starting = list(ar1_1 = list(sigma_sq = fixed(.4), phi = fixed(.5))),
      n_samples = 30, chains = 2, chain_control = list(seed = 91707),
      warmup = FALSE, verbose = FALSE)
    r <- recover(f, sub_sample = list(start = 3, thin = 2))
    rng <- .Random.seed
    r2 <- recover(r, sub_sample = list(start = 7, thin = 4))
    expect_identical(.Random.seed, rng)
    for(k in seq_along(r2$chains)){
      idx <- match(r2$chains[[k]]$recover_iter, r$chains[[k]]$recover_iter)
      expect_identical(r2$chains[[k]]$recover_iter, seq.int(7L, 30L, 4L))
      expect_equal(r2$chains[[k]]$w_samples$ar1_1,
                   r$chains[[k]]$w_samples$ar1_1[idx, , drop = FALSE], ignore_attr = TRUE)
      expect_identical(r2$chains[[k]]$alpha_samples, f$chains[[k]]$alpha_samples)
    }
  }
})
