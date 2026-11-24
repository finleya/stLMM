independent_iid_data <- function(){
  set.seed(85101)
  d <- expand.grid(group = c("a", "b", "c"), parent = c("South", "North"), rep = 1:4)
  d$group <- factor(d$group, levels = c("b", "a", "c", "unused"))
  d$parent <- factor(d$parent, levels = c("North", "South", "unused_parent"))
  d$x <- rnorm(nrow(d)); d$x[c(1, 4)] <- 0
  d$X <- rnorm(nrow(d))
  d$time <- seq_len(nrow(d))
  d$lon <- d$time / 3; d$lat <- (d$time %% 5) / 4
  d$y <- rnorm(nrow(d)); d$y[1] <- NA_real_
  d
}

test_that("independent coefficient formulas share the scalar sampler and prediction path", {
  d <- independent_iid_data()
  prior <- list("(Intercept)" = ig(3, 1), x = ig(4, .4), X = ig(5, 2))
  starts <- c("(Intercept)" = .6, x = .3, X = .8)
  for(collapse_iid in c(FALSE, TRUE)) for(family in c("gaussian", "binomial", "probit", "negative_binomial")){
    dat <- d
    if(family != "gaussian") dat$y <- as.numeric(dat$y > 0)
    for(by in c(FALSE, TRUE)) for(process in c("none", "ar1", "nngp")){
      label <- paste(family, by, process, collapse_iid)
      parent <- if(by) ", by = parent" else ""
      scalar <- paste0("iid(group", parent, ")")
      suffix <- switch(process, none = "", ar1 = " + ar1(time)",
                       nngp = ' + x:nngp(lon, lat, m = 3, ordering = "maxmin")')
      expanded <- as.formula(paste0("y ~ x + X + ", scalar, " + x:", scalar, " + X:", scalar, suffix))
      compact <- as.formula(paste0('y ~ x + X + iid(group, ~ 1 + x + X, cov_model = "independent"', parent, ")", suffix))
      common_start <- if(family == "gaussian") list(resid = list(tau_sq = fixed(.4))) else list()
      if(process == "ar1") common_start$ar1_1 <- list(sigma_sq = fixed(.3), phi = fixed(.4))
      if(process == "nngp") common_start$nngp_1 <- list(sigma_sq = fixed(.3), phi = fixed(.7))
      fit_one <- function(form, priors, starting){
        set.seed(85102)
        f <- stLMM(form, dat, family = family,
          trials = if(family == "binomial") rep(1, nrow(dat)) else NULL,
          size = if(family == "negative_binomial") 2 else NULL,
          priors = priors, starting = starting, n_samples = 20, chains = 2,
          chain_control = list(seed = 85103), warmup = FALSE, verbose = FALSE, collapse_iid = collapse_iid)
        list(fit = f, rng = .Random.seed)
      }
      old_start <- common_start; new_start <- common_start
      # Both supplied starts and automatic dispersed starts must remain aligned.
      if(by){
        old_start <- c(old_start, setNames(lapply(starts, function(x) list(sigma_sq = unname(x))), paste0("iid_", 1:3)))
        new_start$iid_1 <- list(sigma_sq = starts[c("X", "(Intercept)", "x")])
      }
      a <- fit_one(expanded, setNames(lapply(prior, function(x) list(sigma_sq = x)), paste0("iid_", 1:3)), old_start)
      b <- fit_one(compact, list(iid_1 = list(sigma_sq = prior[c("X", "(Intercept)", "x")])), new_start)
      expect_identical(a$rng, b$rng, info = label)
      for(ch in 1:2){
        old <- a$fit$chains[[ch]]; new <- b$fit$chains[[ch]]
        for(nm in c("beta_samples", "alpha_samples", "iid_sigma_sq_samples", "tau_sq_samples",
                    "sigma_sq_samples", "theta_samples"))
          expect_identical(unname(old[[nm]]), unname(new[[nm]]), info = paste(label, nm))
        for(slot in c("p", "i", "x", "Dim"))
          expect_identical(methods::slot(old$backend$Z, slot), methods::slot(new$backend$Z, slot), info = label)
        expect_identical(old$backend$re, new$backend$re, info = label)
        expect_identical(unname(old$backend$sigma_sq_re_IG), unname(new$backend$sigma_sq_re_IG), info = label)
        expect_null(new$iid_covariance_samples)
        expect_true(all(vapply(new$backend$iid_covariance, is.null, logical(1))))
        expect_length(new$backend$Z@x, 3 * nrow(dat))
        expect_true(all(vapply(new$term_description$random_effects$terms,
          function(term) identical(term$formula_term, "iid_1") && identical(term$cov_model, "independent"), logical(1))))
      }
      if(process != "none" || (collapse_iid && family == "gaussian")){
        set.seed(85104); a$fit <- recover(a$fit, sub_sample = list(start = 2, thin = 2)); rng <- .Random.seed
        set.seed(85104); b$fit <- recover(b$fit, sub_sample = list(start = 2, thin = 2))
        expect_identical(.Random.seed, rng)
      }
      sa <- summary(a$fit, type = "combined", terms = paste0("iid_", 1:3), diagnostics = FALSE)
      sb <- summary(b$fit, type = "combined", terms = "iid_1", diagnostics = FALSE)
      expect_identical(sa$parameters, sb$parameters, info = label)
      expect_identical(sa$quantities$coefficient, sb$quantities$coefficient)
      expect_identical(sa$quantities[c("group_group", if(by) "group_parent")],
                       sb$quantities[c("group_group", if(by) "group_parent")])
      nd <- dat[c(2, 2, 5), , drop = FALSE]
      for(new_group in c(FALSE, TRUE)){
        if(new_group){
          nd$group <- "new_group"
          if(process == "ar1") nd$time <- c(30, 30, 31)
          if(process == "nngp") { nd$lon <- nd$lon + .13; nd$lat <- nd$lat + .17 }
          expect_error(predict(b$fit, nd), "allow_new_levels", info = label)
        }
        set.seed(85105); pa <- predict(a$fit, nd, allow_new_levels = new_group, y_samples = TRUE); rng <- .Random.seed
        set.seed(85105); pb <- predict(b$fit, nd, allow_new_levels = new_group, y_samples = TRUE)
        expect_identical(.Random.seed, rng, info = label)
        expect_identical(pa$chains, pb$chains, info = label)
      }
    }
  }
})

test_that("shared controls mean separate variances with repeated hyperparameters", {
  d <- independent_iid_data()
  for(by in c(FALSE, TRUE)){
    form <- if(by) y ~ x + iid(group, ~ 1 + x, cov_model = "independent", by = parent) else
      y ~ x + iid(group, ~ 1 + x, cov_model = "independent")
    make <- function(pr, start){
      set.seed(85106)
      stLMM(form, d, priors = list(iid_1 = list(sigma_sq = pr)),
        starting = list(iid_1 = list(sigma_sq = start)), n_samples = 20, verbose = FALSE)
    }
    a <- make(ig(3, 2), .7)
    b <- make(list("(Intercept)" = ig(3, 2), x = ig(3, 2)), list(x = .7, "(Intercept)" = .7))
    expect_identical(a$samples, b$samples)
    expect_equal(ncol(a$iid_sigma_sq_samples), if(by) 4L else 2L)
    expect_false(identical(a$iid_sigma_sq_samples[, 1], a$iid_sigma_sq_samples[, ncol(a$iid_sigma_sq_samples)]))
    expect_true(all(a$backend$sigma_sq_re_IG[, "shape"] == 3))
    expect_true(all(a$backend$sigma_sq_re_IG[, "scale"] == 2))
    expect_true(all(colnames(a$iid_sigma_sq_samples) %in% names(as_samples(a))))
    expect_output(print(summary(a)), "iid_sigma_sq")
  }
})

test_that("independent formulas support slope-only terms and mixed covariance models", {
  d <- independent_iid_data()
  for(coefficients in c("1", "0 + x", "0 + x + X")){
    form <- as.formula(paste0('y ~ x + X + iid(group, ~ ', coefficients, ', cov_model = "independent")'))
    fit <- stLMM(form, d, priors = list(iid_1 = list(sigma_sq = ig(3, 1))), n_samples = 10, verbose = FALSE)
    expect_null(fit$iid_covariance_samples)
    expect_s3_class(predict(fit, d), "stLMM_prediction")
    expect_s3_class(summary(fit, type = "combined", terms = "all"), "summary_stLMM_combined")
  }
  fit <- stLMM(y ~ x + iid(group, ~ 1 + x, cov_model = "independent") + iid(parent, ~ 1 + x) + iid(group), d,
    priors = list(iid_1 = list(sigma_sq = ig(3, 1)), iid_2 = list(covariance = iw(4, diag(2))),
                  iid_3 = list(sigma_sq = ig(2, 1))), n_samples = 20, chains = 2,
    chain_control = list(seed = 85107), verbose = FALSE)
  f <- fit$chains[[1L]]
  expect_equal(ncol(f$iid_sigma_sq_samples), 3L)
  expect_equal(ncol(f$iid_covariance_samples), 3L)
  expect_identical(tail(colnames(f$iid_sigma_sq_samples), 1L), "iid_3_sigma_sq")
  expect_match(colnames(f$iid_sigma_sq_samples)[1], 'iid_1[coefficient="(Intercept)"]_sigma_sq', fixed = TRUE)
  expect_match(colnames(f$iid_covariance_samples)[1], "iid_2_cov", fixed = TRUE)
  expect_s3_class(summary(fit, type = "combined", terms = "all"), "summary_stLMM_combined")
  expect_s3_class(predict(fit, d), "stLMM_prediction_chains")
  if(requireNamespace("posterior", quietly = TRUE)){
    s <- summary(fit, type = "combined", terms = "iid_1", diagnostic_method = "rank")
    expect_named(s$diagnostics, c("parameter", "rhat", "ess_bulk", "ess_tail"))
  }
})

test_that("independent coefficient controls reject ambiguous or invalid inputs", {
  d <- independent_iid_data()
  fit <- function(prior = list(sigma_sq = ig(3, 1)), ...)
    stLMM(y ~ x + iid(group, ~ 1 + x, cov_model = "independent"), d,
      priors = list(iid_1 = prior), n_samples = 2, verbose = FALSE, ...)
  expect_error(fit(list(covariance = iw(4, diag(2)))), "sigma_sq")
  expect_error(fit(list(sigma_sq = list(ig(3, 1), ig(3, 1)))), "coefficient names")
  expect_error(fit(list(sigma_sq = list(x = ig(3, 1)))), "coefficient names")
  expect_error(fit(list(sigma_sq = list("(Intercept)" = ig(3, 1), X = ig(3, 1)))), "coefficient names")
  expect_error(fit(list(sigma_sq = setNames(list(ig(3, 1), ig(3, 1)), c("x", "x")))), "coefficient names")
  expect_error(fit(list(sigma_sq = list("(Intercept)" = ig(3, 1), x = normal()))), "variance prior")
  expect_error(fit(list(sigma_sq = ig(-1, 1))), "positive")
  expect_error(fit(starting = list(iid_1 = list(sigma_sq = c(.2, .3)))), "coefficient names")
  expect_error(fit(starting = list(iid_1 = list(sigma_sq = c(x = .2)))), "coefficient names")
  expect_error(fit(starting = list(iid_1 = list(sigma_sq = c("(Intercept)" = .2, x = -1)))), "positive")
  expect_error(fit(tuning = list(iid_1 = list(sigma_sq = .1))), "Gibbs")
  bad <- d; bad$x[2] <- Inf
  expect_error(stLMM(y ~ iid(group, ~ x, cov_model = "independent"), bad,
    priors = list(iid_1 = list(sigma_sq = ig(3, 1))), n_samples = 2, verbose = FALSE), "finite numeric")
  expect_error(stLMM:::parse_iid_call('x:iid(group, ~ x, cov_model = "independent")'), "external slope")
})
