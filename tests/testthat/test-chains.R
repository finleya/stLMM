test_that("multi-chain fits keep ordinary fit objects and support coda", {
  set.seed(501)
  dat <- data.frame(y = rnorm(18), x = rnorm(18))

  fit <- stLMM(
    y ~ x,
    data = dat,
    n_samples = 12,
    chains = 3,
    chain_control = list(seed = 99),
    verbose = FALSE
  )

  expect_s3_class(fit, "stLMM_chains")
  expect_length(fit$chains, 3)
  expect_true(all(vapply(fit$chains, inherits, logical(1), "stLMM")))
  expect_equal(dim(fit$chains[[1]]$beta_samples), c(12L, 2L))
  expect_true(all(vapply(fit$chains, function(z) "elapsed" %in% names(z$timing$sampler), logical(1))))
  expect_equal(dim(fit$timing$sampler_by_chain)[1], 3L)
  expect_true(all(c("user.self", "sys.self", "elapsed") %in% colnames(fit$timing$sampler_by_chain)))
  expect_equal(unname(fit$timing$sampler_total[["elapsed"]]),
               sum(fit$timing$sampler_by_chain[, "elapsed"]))
  expect_s3_class(as_mcmc(fit), "mcmc.list")

  s <- summary(fit)
  expect_s3_class(s, "summary_stLMM_chains")
  expect_true(all(c("parameter", "rhat", "effective_size") %in% names(s$diagnostics)))
  s_sel <- summary(fit, parameters = c("x", "tau_sq"))
  expect_equal(rownames(s_sel$parameters), c("x", "tau_sq"))
  expect_equal(s_sel$diagnostics$parameter, c("x", "tau_sq"))
  expect_error(summary(fit, parameters = "not_a_parameter"),
               "unknown parameter")

  s1_burn <- summary(fit$chains[[1]], burn = 2)
  m1 <- as.matrix(as_mcmc(fit$chains[[1]]))
  expect_equal(s1_burn$n_used, 10L)
  expect_equal(
    unname(s1_burn$parameters$beta["(Intercept)", "mean"]),
    mean(m1[3:12, "(Intercept)"])
  )

  s_burn <- summary(fit, burn = 2)
  chain_mats <- lapply(as_mcmc(fit), as.matrix)
  manual_burn <- do.call(rbind, lapply(chain_mats, function(z) z[3:12, , drop = FALSE]))
  expect_equal(s_burn$n_used, 10L)
  expect_equal(
    unname(s_burn$parameters["tau_sq", "mean"]),
    mean(manual_burn[, "tau_sq"])
  )

  chains_burn <- do.call(
    coda::mcmc.list,
    lapply(chain_mats, function(z) coda::mcmc(z[3:12, , drop = FALSE]))
  )
  ess_burn <- coda::effectiveSize(chains_burn)
  expect_equal(
    unname(s_burn$diagnostics["tau_sq", "effective_size"]),
    unname(ess_burn["tau_sq"])
  )
  expect_error(summary(fit$chains[[1]], burn = 12), "burn removes all posterior draws")
  expect_error(summary(fit, burn = 12), "burn removes all posterior draws")

  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  expect_invisible(plot(fit, parameters = c("tau_sq", "(Intercept)"),
                        n_col = 2, burnin = 2, thin = 2))
  expect_invisible(plot(fit, type = "density", parameters = "tau_sq",
                        burnin = 1))
  grDevices::dev.off()
  expect_true(file.exists(f))

  expect_error(plot(fit, burnin = 12), "burnin removes all posterior draws")
})

test_that("as_mcmc drops empty sample blocks in no-fixed-effect chain models", {
  set.seed(504)
  dat <- data.frame(
    y = rnorm(12),
    lon = rep(seq_len(4), 3),
    lat = rep(seq_len(3), each = 4)
  )

  fit <- suppressWarnings(stLMM(
    y ~ 0 + nngp(lon, lat, m = 3),
    data = dat,
    n_samples = 8,
    chains = 3,
    priors = list(
      tau_sq = ig(2, 1),
      nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5))
    ),
    verbose = FALSE
  ))

  m <- as_mcmc(fit)
  expect_s3_class(m, "mcmc.list")
  expect_false("(Intercept)" %in% colnames(as.matrix(m[[1]])))
  expect_true(all(c("tau_sq", "nngp_1_sigma_sq", "nngp_1_phi") %in% colnames(as.matrix(m[[1]]))))
})

test_that("multi-chain starting values recycle or map by chain", {
  set.seed(502)
  dat <- data.frame(
    y = rnorm(20),
    x = rnorm(20),
    time = seq_len(20)
  )

  fit <- stLMM(
    y ~ x + ar1(time),
    data = dat,
    n_samples = 8,
    chains = 3,
    starting = list(
      tau_sq = c(0.5, 1, 2),
      ar1_1 = list(sigma_sq = 1, phi = c(-0.2, 0, 0.2))
    ),
    priors = list(
      tau_sq = ig(2, 1),
      ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.8, 0.8))
    ),
    verbose = FALSE
  )

  expect_equal(vapply(fit$chains, function(z) z$backend$tau_sq_starting, numeric(1)),
               c(chain_1 = 0.5, chain_2 = 1, chain_3 = 2))
  expect_equal(vapply(fit$chains, function(z) z$backend$process_terms[[1]]$theta_starting, numeric(1)),
               c(chain_1 = -0.2, chain_2 = 0, chain_3 = 0.2))
  expect_equal(vapply(fit$chains, function(z) z$backend$process_terms[[1]]$sigma_sq_starting, numeric(1)),
               c(chain_1 = 1, chain_2 = 1, chain_3 = 1))

  expect_error(
    stLMM(
      y ~ x,
      data = dat,
      n_samples = 4,
      chains = 3,
      starting = list(tau_sq = c(1, 2)),
      priors = list(tau_sq = ig(2, 1)),
      verbose = FALSE
    ),
    "must have length 1 or chains"
  )
})

test_that("multi-chain dispersed starts are paired across blocking schemes", {
  set.seed(505)
  dat <- data.frame(
    y = rnorm(24),
    x = rnorm(24),
    lon = rep(seq_len(6), 4),
    lat = rep(seq_len(4), each = 6)
  )

  fit_joint <- suppressWarnings(stLMM(
    y ~ x + nngp(lon, lat, m = 4),
    data = dat,
    n_samples = 6,
    chains = 3,
    chain_control = list(seed = 123, dispersion = 1.5),
    priors = list(
      tau_sq = ig(2, 1),
      nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5))
    ),
    metropolis = "joint",
    verbose = FALSE
  ))

  fit_split <- suppressWarnings(stLMM(
    y ~ x + nngp(lon, lat, m = 4),
    data = dat,
    n_samples = 6,
    chains = 3,
    chain_control = list(seed = 123, dispersion = 1.5),
    priors = list(
      tau_sq = ig(2, 1),
      nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5))
    ),
    metropolis = "variance_theta",
    verbose = FALSE
  ))

  starts <- function(fit){
    t(vapply(fit$chains, function(chain){
      term <- chain$backend$process_terms[[1]]
      c(
        tau_sq = chain$backend$tau_sq_starting,
        sigma_sq = term$sigma_sq_starting,
        stats::setNames(term$theta_starting, term$theta_names)
      )
    }, numeric(3)))
  }

  expect_equal(starts(fit_joint), starts(fit_split))
  expect_gt(length(unique(starts(fit_joint)[, "tau_sq"])), 1L)
  expect_gt(length(unique(starts(fit_joint)[, "sigma_sq"])), 1L)
})

test_that("multi-chain recovery and prediction thin within each chain", {
  set.seed(503)
  dat <- data.frame(
    y = rnorm(16),
    x = rnorm(16),
    time = seq_len(16)
  )

  fit <- stLMM(
    y ~ x + ar1(time),
    data = dat,
    n_samples = 20,
    chains = 4,
    priors = list(
      tau_sq = ig(2, 1),
      ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.8, 0.8))
    ),
    verbose = FALSE
  )

  rec <- recover(fit, sub_sample = list(start = 5, thin = 4))
  expect_s3_class(rec, "stLMM_recovery_chains")
  expect_equal(unname(vapply(rec$chains, function(z) length(z$recover_iter), integer(1))), rep(4L, 4))
  expect_s3_class(as_mcmc(rec), "mcmc.list")

  rec_mcmc <- as_mcmc(rec, include_w = TRUE)
  expect_s3_class(rec_mcmc, "mcmc.list")
  expect_equal(length(rec_mcmc), 4L)
  expect_equal(nrow(as.matrix(rec_mcmc[[1L]])), length(rec$chains[[1L]]$recover_iter))
  expect_equal(
    as.matrix(rec_mcmc[[1L]])[, "(Intercept)"],
    rec$chains[[1L]]$beta_samples[rec$chains[[1L]]$recover_iter, "(Intercept)"]
  )
  expect_equal(
    as.matrix(rec_mcmc[[1L]])[, "w_ar1_1_1"],
    rec$chains[[1L]]$w_samples$ar1_1[, 1]
  )

  rec_mcmc_burn <- as_mcmc(rec, include_w = TRUE, burn = 8)
  expect_equal(nrow(as.matrix(rec_mcmc_burn[[1L]])), 3L)
  expect_equal(
    as.matrix(rec_mcmc_burn[[1L]])[, "w_ar1_1_1"],
    rec$chains[[1L]]$w_samples$ar1_1[2:4, 1]
  )
  expect_s3_class(summary(rec, include_w = TRUE, burn = 8), "summary_stLMM_recovery_chains")

  pred <- predict(rec, sub_sample = list(start = 9, thin = 4))
  expect_s3_class(pred, "stLMM_prediction_chains")
  expect_equal(unname(vapply(pred$chains, function(z) nrow(z$mu_samples), integer(1))), rep(3L, 4))
  expect_equal(ncol(pred$chains[[1]]$mu_samples), nrow(dat))
})
