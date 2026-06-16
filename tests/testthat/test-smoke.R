test_that("fixed-effect model fits with underscore API", {
  set.seed(1)
  dat <- data.frame(y = rnorm(12), x = rnorm(12))

  fit <- stLMM(
    y ~ x,
    data = dat,
    starting = list(tau_sq = 1),
    priors = list(tau_sq = ig(2, 1)),
    n_samples = 4,
    verbose = FALSE
  )

  expect_equal(dim(fit$beta_samples), c(4L, 2L))
  expect_null(fit$alpha_samples)
  expect_null(fit$samples$alpha)
  expect_equal(ncol(fit$tau_sq_samples), NULL)
  expect_true(all(c("user.self", "sys.self", "elapsed") %in% names(fit$timing$sampler)))
  expect_true(is.finite(fit$timing$sampler[["elapsed"]]))
  expect_gte(fit$timing$sampler[["elapsed"]], 0)
})

test_that("starting values are optional", {
  set.seed(4)
  dat <- data.frame(
    y = rnorm(16),
    x = rnorm(16),
    group = factor(rep(seq_len(4), each = 4)),
    lon = rep(seq_len(4), 4),
    lat = rep(seq_len(4), each = 4)
  )

  fit_fixed <- stLMM(
    y ~ x,
    data = dat,
    priors = list(tau_sq = ig(2, 1)),
    n_samples = 4,
    verbose = FALSE
  )

  fit_re <- stLMM(
    y ~ x + iid(group),
    data = dat,
    priors = list(tau_sq = ig(2, 1), iid_1 = list(sigma_sq = ig(2, 1))),
    n_samples = 4,
    verbose = FALSE
  )

  fit_proc <- suppressWarnings(stLMM(
    y ~ x + nngp(lon, lat, m = 2),
    data = dat,
    priors = list(
      tau_sq = ig(2, 1),
      nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5))
    ),
    n_samples = 4,
    verbose = FALSE
  ))

  fit_partial_proc <- suppressWarnings(stLMM(
    y ~ x + nngp(lon, lat, m = 2),
    data = dat,
    starting = list(nngp_1 = c(phi = 0.7)),
    priors = list(
      tau_sq = ig(2, 1),
      nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5))
    ),
    n_samples = 4,
    verbose = FALSE
  ))

  expect_equal(dim(fit_fixed$beta_samples), c(4L, 2L))
  expect_equal(dim(fit_re$alpha_samples), c(4L, 4L))
  expect_equal(unname(colnames(fit_proc$sigma_sq_samples)), "nngp_1_sigma_sq")
  expect_true(is.finite(fit_fixed$backend$tau_sq_starting))
  expect_gt(fit_fixed$backend$tau_sq_starting, 0)
  expect_equal(fit_re$backend$sigma_sq_re_starting, 1)
  expect_equal(fit_proc$backend$process_terms[[1]]$sigma_sq_starting, 1)
  expect_equal(fit_proc$backend$process_terms[[1]]$theta_starting, 2.55)
  expect_equal(fit_partial_proc$backend$process_terms[[1]]$sigma_sq_starting, 1)
  expect_equal(fit_partial_proc$backend$process_terms[[1]]$theta_starting, 0.7)
})

test_that("Matern covariance model is available for dense GP and NNGP terms", {
  set.seed(713)
  dat <- data.frame(
    y = rnorm(12),
    x = rnorm(12),
    lon = rep(seq_len(4), 3),
    lat = rep(seq_len(3), each = 4)
  )

  models <- get_cor_models()
  expect_true("matern" %in% names(models))
  expect_equal(models$matern$names, c("phi", "nu"))
  expect_equal(models$matern$types, c(1L, 1L))
  expect_equal(models$matern$distance_mode, 1L)
  expect_true("gneiting" %in% names(models))
  expect_equal(models$gneiting$names, c("a", "c", "alpha", "beta", "gamma", "delta"))
  expect_equal(models$gneiting$types, c(1L, 1L, 2L, 2L, 2L, 1L))
  expect_equal(models$gneiting$distance_mode, 2L)
  expect_false("gneiting_exp_cauchy" %in% names(models))
  expect_false("nonsep_gneiting" %in% names(models))

  gp_fit <- suppressWarnings(stLMM(
    y ~ x + gp(lon, lat, cov_model = "matern"),
    data = dat,
    starting = list(tau_sq = 0.5, gp_1 = c(sigma_sq = 1, phi = 0.8, nu = 1.2)),
    tuning = list(tau_sq = 0, gp_1 = c(sigma_sq = 0, phi = 0, nu = 0)),
    priors = list(
      tau_sq = ig(2, 0.5),
      gp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5), nu = uniform(0.2, 3))
    ),
    n_samples = 4,
    verbose = FALSE
  ))

  nngp_fit <- suppressWarnings(stLMM(
    y ~ x + nngp(lon, lat, m = 3, cov_model = "matern", ordering = "maxmin"),
    data = dat,
    starting = list(tau_sq = 0.5, nngp_1 = c(sigma_sq = 1, phi = 0.8, nu = 1.2)),
    tuning = list(tau_sq = 0, nngp_1 = c(sigma_sq = 0, phi = 0, nu = 0)),
    priors = list(
      tau_sq = ig(2, 0.5),
      nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5), nu = uniform(0.2, 3))
    ),
    n_samples = 4,
    verbose = FALSE
  ))

  expect_equal(gp_fit$backend$process_terms[[1]]$cov_model, "matern")
  expect_equal(nngp_fit$backend$process_terms[[1]]$cov_model, "matern")
  expect_equal(colnames(gp_fit$theta_samples), c("gp_1_phi", "gp_1_nu"))
  expect_equal(colnames(nngp_fit$theta_samples), c("nngp_1_phi", "nngp_1_nu"))
  expect_equal(unname(gp_fit$theta_samples[1, ]), c(0.8, 1.2))
  expect_equal(unname(nngp_fit$theta_samples[1, ]), c(0.8, 1.2))

  gp_rec <- recover(gp_fit, sub_sample = list(start = 2, thin = 2))
  nngp_rec <- recover(nngp_fit, sub_sample = list(start = 2, thin = 2))
  expect_equal(ncol(gp_rec$w_samples$gp_1), nrow(unique(dat[, c("lon", "lat")])))
  expect_equal(ncol(nngp_rec$w_samples$nngp_1), nrow(unique(dat[, c("lon", "lat")])))
})

test_that("Matern covariance with nu 0.5 matches exponential covariance", {
  set.seed(715)
  dat <- data.frame(
    y = rnorm(10),
    x = rnorm(10),
    lon = seq(0, 1, length.out = 10),
    lat = seq(0, 0.5, length.out = 10)
  )

  fixed_resid <- list(tau_sq = fixed(0.2))
  exp_start <- list(
    resid = fixed_resid,
    gp_1 = list(sigma_sq = fixed(1.1), phi = fixed(2.3))
  )
  matern_start <- list(
    resid = fixed_resid,
    gp_1 = list(sigma_sq = fixed(1.1), phi = fixed(2.3), nu = fixed(0.5))
  )

  set.seed(716)
  exp_fit <- suppressWarnings(stLMM(
    y ~ x + gp(lon, lat, cov_model = "exp"),
    data = dat,
    starting = exp_start,
    n_samples = 12,
    verbose = FALSE
  ))

  set.seed(716)
  matern_fit <- suppressWarnings(stLMM(
    y ~ x + gp(lon, lat, cov_model = "matern"),
    data = dat,
    starting = matern_start,
    n_samples = 12,
    verbose = FALSE
  ))

  expect_equal(matern_fit$beta_samples, exp_fit$beta_samples, tolerance = 1e-12)
})

test_that("CHOLMOD ordering controls are passed to sparse factorization", {
  set.seed(714)
  dat <- data.frame(
    y = rnorm(12),
    x = rnorm(12),
    lon = rep(seq_len(4), 3),
    lat = rep(seq_len(3), each = 4)
  )

  fit_natural <- suppressWarnings(stLMM(
    y ~ x + nngp(lon, lat, m = 3),
    data = dat,
    starting = list(
      tau_sq = fixed(0.5),
      nngp_1 = list(sigma_sq = fixed(1), phi = fixed(1))
    ),
    n_samples = 4,
    cholmod_control = list(ordering = "natural", postorder = FALSE),
    verbose = FALSE
  ))

  fit_amd <- suppressWarnings(stLMM(
    y ~ x + nngp(lon, lat, m = 3),
    data = dat,
    starting = list(
      tau_sq = fixed(0.5),
      nngp_1 = list(sigma_sq = fixed(1), phi = fixed(1))
    ),
    n_samples = 4,
    cholmod_control = "amd",
    verbose = FALSE
  ))

  expect_equal(fit_natural$backend$cholmod_control$ordering, "natural")
  expect_false(fit_natural$backend$cholmod_control$postorder)
  expect_equal(fit_natural$term_description$global$cholmod_requested_ordering, "natural")
  expect_false(fit_natural$term_description$global$cholmod_postorder)
  expect_equal(fit_natural$term_description$global$cholmod_ordering, "natural")

  expect_equal(fit_amd$backend$cholmod_control$ordering, "amd")
  expect_true(fit_amd$backend$cholmod_control$postorder)
  expect_equal(fit_amd$term_description$global$cholmod_requested_ordering, "amd")
  expect_true(fit_amd$term_description$global$cholmod_postorder)
  expect_equal(fit_amd$term_description$global$cholmod_ordering, "amd")

  expect_error(
    stLMM(
      y ~ x,
      data = dat,
      n_samples = 4,
      cholmod_control = list(ordering = "not_an_ordering"),
      verbose = FALSE
    ),
    "'arg' should be one of"
  )
})

test_that("tau_sq controls are accepted and missing process priors are clear", {
  set.seed(5)
  dat <- data.frame(
    y = rnorm(8),
    x = rnorm(8),
    lon = rep(seq_len(4), 2),
    lat = rep(c(1, 2), 4)
  )

  fit <- stLMM(
    y ~ x,
    data = dat,
    starting = list(tau_sq = 0.5),
    tuning = list(tau_sq = 0),
    priors = list(tau_sq = ig(3, 2)),
    n_samples = 4,
    verbose = FALSE
  )

  expect_equal(fit$backend$tau_sq_starting, 0.5)
  expect_equal(fit$backend$tau_sq_tuning, 0)
  expect_equal(fit$backend$tau_sq_IG, c(3, 2))
  expect_error(
    suppressWarnings(stLMM(y ~ x + nngp(lon, lat, m = 2), data = dat, n_samples = 4, verbose = FALSE)),
    "missing prior for free process parameter sigma_sq"
  )
})

test_that("missing grouped random-effect priors are clear", {
  set.seed(8)
  dat <- data.frame(
    y = rnorm(8),
    x = rnorm(8),
    group = factor(rep(1:4, each = 2))
  )

  expect_error(
    stLMM(
      y ~ x + iid(group),
      data = dat,
      n_samples = 4,
      priors = list(tau_sq = ig(2, 1)),
      verbose = FALSE
    ),
    "missing prior\\(s\\) for grouped random-effect variance iid_1"
  )
})

test_that("random-effect pipe syntax points to iid", {
  set.seed(18)
  dat <- data.frame(
    y = rnorm(8),
    x = rnorm(8),
    group = factor(rep(1:4, each = 2))
  )

  expect_error(
    stLMM(
      y ~ x + (1 | group),
      data = dat,
      n_samples = 4,
      priors = list(tau_sq = ig(2, 1), iid_1 = list(sigma_sq = ig(2, 1))),
      verbose = FALSE
    ),
    "use iid\\(group\\)"
  )
})

test_that("process-term pipe syntax points to current term forms", {
  set.seed(19)
  dat <- data.frame(
    y = rnorm(8),
    group = factor(rep(1:4, each = 2)),
    time = rep(1:4, each = 2),
    lon = rep(1:4, each = 2),
    lat = rep(c(1, 2), 4)
  )

  expect_error(
    stLMM(y ~ ar1(time | group), data = dat, n_samples = 4, verbose = FALSE),
    "use ar1\\(\\.\\.\\.\\) without \\|"
  )
  expect_error(
    stLMM(y ~ gp(lon, lat | group), data = dat, n_samples = 4, verbose = FALSE),
    "use gp\\(\\.\\.\\.\\) without \\|"
  )
  expect_error(
    stLMM(y ~ nngp(lon, lat | group), data = dat, n_samples = 4, verbose = FALSE),
    "use nngp\\(\\.\\.\\.\\) without \\|"
  )
})

test_that("recover requires structured process terms", {
  set.seed(9)
  dat <- data.frame(y = rnorm(8), x = rnorm(8))
  fit <- stLMM(
    y ~ x,
    data = dat,
    n_samples = 4,
    priors = list(tau_sq = ig(2, 1)),
    verbose = FALSE
  )

  expect_error(recover(fit), "recover\\(\\) requires at least one structured process term")
})

test_that("tuning controls and validation are explicit", {
  set.seed(6)
  dat <- data.frame(
    y = rnorm(12),
    x = rnorm(12),
    group = factor(rep(seq_len(3), each = 4)),
    lon = rep(seq_len(4), 3),
    lat = rep(c(1, 2, 3), each = 4)
  )

  fit_re_nested <- stLMM(
    y ~ x + iid(group),
    data = dat,
    starting = list(iid_1 = list(sigma_sq = 0.7)),
    priors = list(tau_sq = ig(2, 1), iid_1 = list(sigma_sq = ig(3, 2))),
    n_samples = 4,
    verbose = FALSE
  )

  expect_error(
    stLMM(
      y ~ x + iid(group),
      data = dat,
      starting = list(sigma_sq_re = list(iid_1 = 0.7)),
      priors = list(tau_sq = ig(2, 1), iid_1 = list(sigma_sq = ig(3, 2))),
      n_samples = 4,
      verbose = FALSE
    ),
    "starting\\$sigma_sq_re is no longer supported"
  )

  expect_error(
    stLMM(
      y ~ x + iid(group),
      data = dat,
      priors = list(tau_sq = ig(2, 1), sigma_sq_re = list(iid_1 = ig(3, 2))),
      n_samples = 4,
      verbose = FALSE
    ),
    "priors\\$sigma_sq_re is no longer supported"
  )

  expect_error(
    stLMM(
      y ~ x + iid(group),
      data = dat,
      priors = list(tau_sq = ig(2, 1), iid_1 = ig(3, 2)),
      n_samples = 4,
      verbose = FALSE
    ),
    "priors\\$iid_1 must be a list with entry sigma_sq"
  )

  fit_proc_default_tuning <- suppressWarnings(stLMM(
    y ~ x + nngp(lon, lat, m = 2),
    data = dat,
    priors = list(
      tau_sq = ig(2, 1),
      nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5))
    ),
    n_samples = 4,
    verbose = FALSE
  ))

  fit_proc_fixed <- suppressWarnings(stLMM(
    y ~ x + nngp(lon, lat, m = 2),
    data = dat,
    tuning = list(tau_sq = 0, nngp_1 = c(sigma_sq = 0, phi = 0)),
    priors = list(
      tau_sq = ig(2, 1),
      nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5))
    ),
    n_samples = 4,
    verbose = FALSE
  ))

  expect_equal(fit_re_nested$backend$sigma_sq_re_starting, 0.7)
  expect_equal(unname(fit_re_nested$backend$sigma_sq_re_IG["iid_1", ]), c(3, 2))
  expect_equal(fit_proc_default_tuning$backend$process_terms[[1]]$sigma_sq_tuning, 0.1)
  expect_equal(fit_proc_default_tuning$backend$process_terms[[1]]$theta_tuning, 0.1)
  expect_equal(fit_proc_fixed$backend$tau_sq_tuning, 0)
  expect_equal(fit_proc_fixed$backend$process_terms[[1]]$sigma_sq_tuning, 0)
  expect_equal(fit_proc_fixed$backend$process_terms[[1]]$theta_tuning, 0)

  expect_error(
    suppressWarnings(stLMM(
      y ~ x + nngp(lon, lat, m = 2),
      data = dat,
      starting = list(nngp_1 = c(foo = 1)),
      priors = list(tau_sq = ig(2, 1), nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5))),
      n_samples = 4,
      verbose = FALSE
    )),
    "unknown parameter\\(s\\) foo in starting\\$nngp_1"
  )

  expect_error(
    suppressWarnings(stLMM(
      y ~ x + nngp(lon, lat, m = 2),
      data = dat,
      tuning = list(nngp_1 = c(sigma_sq = 0)),
      priors = list(tau_sq = ig(2, 1), nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5))),
      n_samples = 4,
      verbose = FALSE
    )),
    "missing parameter\\(s\\) phi in tuning\\$nngp_1"
  )

  expect_error(
    suppressWarnings(stLMM(
      y ~ x + nngp(lon, lat, m = 2),
      data = dat,
      priors = list(tau_sq = ig(2, 1), nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5), extra = uniform(1, 2))),
      n_samples = 4,
      verbose = FALSE
    )),
    "unknown parameter\\(s\\) extra in priors\\$nngp_1"
  )
})

test_that("covariance Metropolis diagnostics report the global block", {
  set.seed(10)
  dat <- data.frame(
    y = rnorm(18),
    x = rnorm(18),
    lon = rep(seq_len(6), 3),
    lat = rep(c(1, 2, 3), each = 6)
  )

  fit <- suppressWarnings(stLMM(
    y ~ x + nngp(lon, lat, m = 2),
    data = dat,
    priors = list(
      tau_sq = ig(2, 1),
      nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5))
    ),
    n_samples = 60,
    verbose = FALSE
  ))

  am <- fit$adaptive_metropolis
  expect_true(is.numeric(fit$covariance_acceptance))
  expect_equal(fit$covariance_acceptance, am$acceptance)
  expect_true(is.na(fit$term_param_accept["nngp_1"]))
  expect_equal(am$dimension, 3L)
  expect_equal(am$parameter_labels, c("log(tau_sq)", "log(nngp_1_sigma_sq)", "nngp_1_theta_1"))
  expect_equal(names(am$current_eta), am$parameter_labels)
  expect_equal(am$batch_length, 25L)
  expect_equal(length(am$batch_acceptance_history), 2L)
  expect_equal(length(am$proposal_scale_history), 2L)
  expect_true(is.list(am$warmup))
  expect_true(am$warmup$enabled)
  expect_equal(names(am$warmup$starting_transformed), am$parameter_labels)
  expect_equal(names(am$warmup$ending_transformed), am$parameter_labels)
  expect_equal(dim(am$warmup$starting_proposal_cov), c(3L, 3L))
  expect_equal(dim(am$warmup$ending_proposal_cov), c(3L, 3L))
  expect_equal(length(am$warmup$batch_acceptance), am$warmup$n_batches)
  expect_equal(length(am$warmup$proposal_scale_history), am$warmup$n_batches)
  expect_true(am$warmup$n_attempted >= am$warmup$n_batches)
  expect_equal(nrow(fit$tau_sq_samples), NULL)
  expect_equal(length(fit$tau_sq_samples), 60L)
  expect_equal(nrow(fit$samples$beta), 60L)
  expect_true(is.na(fit$term_description$process_terms$nngp_1$sampler$block_acceptance))
  expect_equal(fit$term_description$process_terms$nngp_1$sampler$status, "global covariance block")

  fit_fixed <- suppressWarnings(stLMM(
    y ~ x + nngp(lon, lat, m = 2),
    data = dat,
    tuning = list(tau_sq = 0, nngp_1 = c(sigma_sq = 0, phi = 0)),
    priors = list(
      tau_sq = ig(2, 1),
      nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5))
    ),
    n_samples = 8,
    verbose = FALSE
  ))

  expect_equal(fit_fixed$adaptive_metropolis$dimension, 0L)
  expect_true(is.na(fit_fixed$covariance_acceptance))
  expect_length(fit_fixed$adaptive_metropolis$parameter_labels, 0L)
  expect_false(fit_fixed$adaptive_metropolis$warmup$enabled)
})

test_that("covariance warmup can be disabled without changing sample counts", {
  set.seed(11)
  dat <- data.frame(
    y = rnorm(18),
    x = rnorm(18),
    lon = rep(seq_len(6), 3),
    lat = rep(c(1, 2, 3), each = 6)
  )

  fit <- suppressWarnings(stLMM(
    y ~ x + nngp(lon, lat, m = 2),
    data = dat,
    priors = list(
      tau_sq = ig(2, 1),
      nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5))
    ),
    n_samples = 20,
    warmup = FALSE,
    verbose = FALSE
  ))

  expect_equal(nrow(fit$samples$beta), 20L)
  expect_false(fit$adaptive_metropolis$warmup$enabled)
  expect_equal(fit$adaptive_metropolis$warmup$n_attempted, 0L)
  expect_equal(length(fit$adaptive_metropolis$warmup$batch_acceptance), 0L)
})

test_that("covariance warmup honors minimum batches", {
  set.seed(12)
  dat <- data.frame(y = rnorm(12), x = rnorm(12))

  fit <- stLMM(
    y ~ x,
    data = dat,
    priors = list(tau_sq = ig(2, 1)),
    n_samples = 6,
    warmup = list(batch_length = 2, min_batches = 2, max_batches = 2),
    verbose = FALSE
  )

  expect_true(fit$adaptive_metropolis$warmup$enabled)
  expect_equal(fit$adaptive_metropolis$warmup$min_batches, 2L)
  expect_equal(fit$adaptive_metropolis$warmup$max_batches, 2L)
  expect_equal(fit$adaptive_metropolis$warmup$n_batches, 2L)
  expect_equal(length(fit$adaptive_metropolis$warmup$batch_acceptance), 2L)

  expect_error(
    stLMM(
      y ~ x,
      data = dat,
      priors = list(tau_sq = ig(2, 1)),
      n_samples = 4,
      warmup = list(min_batches = 3, max_batches = 2),
      verbose = FALSE
    ),
    "warmup\\$min_batches"
  )
})

test_that("grouped random-effect model fits", {
  set.seed(2)
  dat <- data.frame(
    y = rnorm(20),
    x = rnorm(20),
    group = factor(rep(seq_len(5), each = 4))
  )

  fit <- stLMM(
    y ~ x + iid(group),
    data = dat,
    starting = list(tau_sq = 1, iid_1 = list(sigma_sq = 1)),
    priors = list(tau_sq = ig(2, 1), iid_1 = list(sigma_sq = ig(2, 1))),
    n_samples = 4,
    verbose = FALSE
  )

  expect_equal(dim(fit$beta_samples), c(4L, 2L))
  expect_equal(dim(fit$alpha_samples), c(4L, 5L))
  expect_null(fit$sigma_sq_re_samples)
  expect_equal(colnames(fit$iid_sigma_sq_samples), "iid_1_sigma_sq")
  expect_equal(colnames(fit$samples$iid_sigma_sq), "iid_1_sigma_sq")
})

test_that("NNGP process model fits and reports term metadata", {
  set.seed(3)
  dat <- data.frame(
    y = rnorm(8),
    x = rnorm(8),
    lon = rep(seq_len(4), 2),
    lat = rep(c(1, 2), 4)
  )

  fit <- suppressWarnings(stLMM(
    y ~ x + nngp(lon, lat, m = 2),
    data = dat,
    starting = list(tau_sq = 1, nngp_1 = c(sigma_sq = 1, phi = 1)),
    priors = list(
      tau_sq = ig(2, 1),
      nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5))
    ),
    n_samples = 4,
    verbose = FALSE
  ))
  rec <- recover(fit, sub_sample = list(start = 1, thin = 2))

  expect_equal(dim(fit$beta_samples), c(4L, 2L))
  expect_equal(unname(colnames(fit$sigma_sq_samples)), "nngp_1_sigma_sq")
  expect_equal(fit$term_description$global$n_process_terms, 1L)
  expect_null(fit$w_samples)
  expect_named(rec$w_samples, "nngp_1")
  expect_named(rec$w_samples_ordered, "nngp_1")
  expect_equal(nrow(rec$w_samples$nngp_1), 2L)
  expect_equal(attr(rec$w_samples$nngp_1, "node_order"), "support")
  expect_equal(attr(rec$w_samples_ordered$nngp_1, "node_order"), "internal")
})

test_that("stLMM nngp_search selects neighbor search backend", {
  set.seed(712)
  dat <- data.frame(
    y = rnorm(8),
    x = rnorm(8),
    lon = runif(8),
    lat = runif(8)
  )

  fit <- stLMM(
    y ~ x + nngp(lon, lat, m = 2),
    data = dat,
    starting = list(tau_sq = 1, nngp_1 = c(sigma_sq = 1, phi = 1)),
    priors = list(
      tau_sq = ig(2, 1),
      nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5))
    ),
    n_samples = 4,
    nngp_search = "brute",
    verbose = FALSE
  )

  expect_equal(fit$backend$nngp_search, "brute")
  expect_equal(fit$backend$graphs[[1]]$nngp_search, "brute")

  expect_error(
    stLMM(
      y ~ x + nngp(lon, lat, m = 2),
      data = dat,
      n_samples = 4,
      nngp_search = "approx",
      verbose = FALSE
    ),
    "'arg' should be one of"
  )
})
