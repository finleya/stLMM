nb_size_data <- function(){
  set.seed(59101)
  d <- data.frame(x = rnorm(36), g = rep(letters[1:6], each = 6),
    parent = rep(c("A", "B"), each = 18), time = 1:36,
    lon = runif(36), lat = runif(36), off = rep(c(-.2, .3), 18))
  d$y <- rnbinom(36, size = 3, mu = exp(.7 + .2 * d$x + d$off))
  d
}

test_that("NB size estimation requires explicit compatible controls", {
  d <- nb_size_data()
  fit <- function(...) stLMM(y ~ x, d, family = "nb", n_samples = 5, verbose = FALSE, ...)
  expect_error(fit(), "size must be supplied")
  expect_error(fit(size = 2, priors = list(size = gamma_dist(2, 1))), "either")
  expect_error(fit(size = c(1, 2)), "scalar")
  expect_error(fit(size = 2, starting = list(size = 1)), "require priors")
  expect_error(fit(size = fixed(2), tuning = list(size = .1)), "require priors")
  for(pr in list(ig(2, 1), log_normal(0, 1), gamma_dist(2, 1, support = c(1, 5)),
                 gamma_dist(-1, 1), gamma_dist(2, Inf), c(2, 1)))
    expect_error(fit(priors = list(size = pr)), "prior|gamma_dist")
  for(start in list(0, Inf, NA_real_, fixed(2), c(1, 2), "2"))
    expect_error(fit(priors = list(size = gamma_dist(2, 1)), starting = list(size = start)), "starting\\$size")
  for(tune in list(0, Inf, NA_real_, c(.1, .2), "1"))
    expect_error(fit(priors = list(size = gamma_dist(2, 1)), tuning = list(size = tune)), "tuning\\$size")
  expect_error(stLMM(y ~ x, d, priors = list(size = gamma_dist(2, 1)),
                    n_samples = 5, verbose = FALSE), "only with family")
  expect_no_warning(fit(priors = list(size = gamma_dist(2, 1))))
  set.seed(59102); a <- fit(size = 3)
  set.seed(59102); b <- fit(size = fixed(3))
  expect_identical(a$samples, b$samples)
  expect_null(a$size_samples)
  expect_null(a$adaptive_metropolis$size)
})

test_that("size has scalar warmup and retained adaptation without covariance blocks", {
  d <- nb_size_data()
  f <- stLMM(y ~ x, d, family = "nb", priors = list(size = gamma_dist(2, .5)),
    starting = list(size = 2), tuning = list(size = .3), n_samples = 53,
    metropolis = list(batch_length = 10),
    warmup = list(batch_length = 8, min_batches = 2, max_batches = 2), verbose = FALSE)
  a <- f$adaptive_metropolis$size
  expect_equal(a$parameter_labels, "log(size)")
  expect_equal(a$target_accept, .44)
  expect_equal(a$batch_index, 5L)
  expect_length(a$batch_acceptance_history, 5L)
  expect_length(a$proposal_scale_history, 5L)
  expect_equal(a$warmup$n_attempted, 16L)
  expect_equal(a$warmup$n_batches, 2L)
  expect_equal(a$warmup$target, c(.3, .6))
  expect_equal(a$warmup$starting_proposal_cov, .09)
  expect_false(f$adaptive_metropolis$warmup$enabled)
  expect_equal(f$adaptive_metropolis$warmup$n_batches, 0L)
  expect_equal(f$adaptive_metropolis$warmup$stopped_reason, "skipped")
  expect_true(is.finite(a$proposal_cov) && a$proposal_cov > 0)
  expect_true(a$acceptance >= 0 && a$acceptance <= 1)
  expect_true(all(is.finite(f$size_samples) & f$size_samples > 0))
  expect_gt(length(unique(f$size_samples)), 1L)
  expect_equal(colnames(f$samples$size), "size")
  expect_true("size" %in% names(as_samples(f)))
  expect_true("size" %in% colnames(as_mcmc(f)))
  for(warmup in list(FALSE, list(max_batches = 0))){
    g <- stLMM(y ~ x, d, family = "nb", priors = list(size = gamma_dist(2, 1)),
      n_samples = 13, warmup = warmup, metropolis = list(target_accept = .37), verbose = FALSE)
    expect_equal(g$adaptive_metropolis$size$target_accept, .37)
    expect_equal(g$adaptive_metropolis$size$warmup$stopped_reason, "skipped")
    expect_equal(g$adaptive_metropolis$size$warmup$n_attempted, 0L)
  }
})

test_that("size starts and diagnostics retain chain identity", {
  d <- nb_size_data()
  f <- stLMM(y ~ x, d, family = "nb", priors = list(size = gamma_dist(2, 1)),
    starting = list(size = c(1, 5)), chains = 2, n_samples = 100,
    chain_control = list(seed = 59103), warmup = FALSE, verbose = FALSE)
  expect_equal(vapply(f$chains, function(z) z$backend$nb_size, numeric(1)), c(chain_1 = 1, chain_2 = 5))
  expect_equal(vapply(f$chains, function(z) z$adaptive_metropolis$size$warmup$starting_transformed,
                      numeric(1)), log(c(chain_1 = 1, chain_2 = 5)))
  expect_equal(rownames(summary(f, parameters = "size")$parameters), "size")
  s <- summary(f, parameters = "size", diagnostic_method = "rank")
  expect_true(all(is.finite(as.matrix(s$diagnostics[, c("rhat", "ess_bulk", "ess_tail")]))))
  a <- as_samples(f, burn = 10, thin = 3)
  for(ch in 1:2)
    expect_equal(a$size[a$.chain == ch], f$chains[[ch]]$size_samples[seq(11, 100, 3)])
})

test_that("prediction and log likelihood use matched size draws", {
  d <- nb_size_data(); d$y[c(2, 7)] <- NA_real_
  f <- stLMM(y ~ x + offset(off), d, family = "nb",
    priors = list(size = gamma_dist(2, .5)), n_samples = 80, verbose = FALSE)
  idx <- seq(11, 80, 7)
  nd <- d[1:4, c("x", "off")]
  set.seed(59104)
  p <- predict(f, nd, y_samples = TRUE, sub_sample = list(start = 11, thin = 7))
  set.seed(59104)
  expected <- t(vapply(seq_along(idx), function(i)
    rnbinom(nrow(nd), size = f$size_samples[idx[i]], mu = p$mu_samples[i, ]), numeric(nrow(nd))))
  expect_equal(unname(p$y_samples), unname(expected))
  eta <- fitted(f, summary = FALSE, sub_sample = list(start = 11, thin = 7))
  obs <- which(!is.na(d$y))
  expected_ll <- t(vapply(seq_along(idx), function(i)
    dnbinom(d$y[obs], size = f$size_samples[idx[i]], mu = eta[i, obs], log = TRUE), numeric(length(obs))))
  expect_equal(unname(as.matrix(log_lik(f, sub_sample = list(start = 11, thin = 7)))),
               unname(expected_ll), ignore_attr = TRUE)
  broken <- f; broken$size_samples <- NULL
  expect_error(predict(broken, nd, y_samples = TRUE), "size draws")
  expect_error(log_lik(broken), "size draws")
})

test_that("estimated size works across IID and structured terms with saved draw selection", {
  d <- nb_size_data(); d$y[2] <- NA_real_
  forms <- list(y ~ x + iid(g), y ~ x + iid(g, ~ 1 + x),
    y ~ x + iid(g, ~ 1 + x, cov_model = "independent", by = parent),
    y ~ x + iid(g, ~ 1 + x) + ar1(time),
    y ~ x + iid(g) + nngp(lon, lat, m = 3),
    y ~ x + iid(g) + gp(lon, lat))
  for(i in seq_along(forms)){
    priors <- list(size = gamma_dist(2, .5),
      iid_1 = if(i %in% c(2, 4)) list(covariance = iw(4, diag(2))) else list(sigma_sq = ig(3, 1)),
      ar1_1 = list(sigma_sq = ig(3, .2), phi = uniform(-.8, .8)),
      nngp_1 = list(sigma_sq = ig(3, .2), phi = uniform(.2, 3)),
      gp_1 = list(sigma_sq = ig(3, .2), phi = uniform(.2, 3)))
    warnings <- character()
    f <- withCallingHandlers(stLMM(forms[[i]], d, family = "nb", priors = priors,
      n_samples = 30, chains = 2, chain_control = list(seed = 59105), verbose = FALSE,
      save_process = if(i >= 4) list(start = 2, thin = 2) else NULL),
      warning = function(w){warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")})
    expect_length(warnings, 1L)
    expect_match(warnings, "weakly identified")
    if(i >= 4) f <- recover(f, sub_sample = list(start = 4, thin = 4))
    for(ch in f$chains){
      p <- predict(ch, y_samples = TRUE, sub_sample = list(start = 4, thin = 4))
      idx <- p$draw_index
      set.seed(59106)
      p <- predict(ch, y_samples = TRUE, sub_sample = list(start = 4, thin = 4))
      set.seed(59106)
      expected <- t(vapply(seq_along(idx), function(j)
        rnbinom(nrow(d), size = ch$size_samples[idx[j]], mu = p$mu_samples[j, ]), numeric(nrow(d))))
      expect_equal(unname(p$y_samples), unname(expected))
      eta <- fitted(ch, summary = FALSE, sub_sample = list(start = 4, thin = 4))
      obs <- which(!is.na(d$y))
      expect_equal(unname(as.matrix(log_lik(ch, sub_sample = list(start = 4, thin = 4)))),
        unname(t(vapply(seq_along(idx), function(j)
          dnbinom(d$y[obs], size = ch$size_samples[idx[j]], mu = eta[j, obs], log = TRUE),
          numeric(length(obs))))), ignore_attr = TRUE)
    }
  }
})

test_that("the first size proposal agrees with the independent log-scale NB target", {
  d <- data.frame(y = c(0, 1, 7, 2, 0, 4), eta = c(-2, 0, 1, .2, -.5, 2))
  target <- function(r) sum(dnbinom(d$y, size = r, mu = exp(d$eta), log = TRUE)) +
    dgamma(r, shape = 2, rate = .2, log = TRUE) + log(r)
  for(r in c(.4, 4, 40, 1000)) for(seed in 59110:59124){
    set.seed(seed)
    proposal <- exp(log(r) + .8 * rnorm(1))
    expected <- if(log(runif(1)) < target(proposal) - target(r)) proposal else r
    set.seed(seed)
    f <- stLMM(y ~ 0 + offset(eta), d, family = "nb",
      priors = list(size = gamma_dist(2, .2)), starting = list(size = r),
      tuning = list(size = .8), warmup = FALSE, n_samples = 1, verbose = FALSE)
    expect_equal(f$size_samples, expected, tolerance = 1e-13)
  }
})

test_that("estimated size supports areal terms and the random effect warning is once per fit", {
  d <- nb_size_data(); d$area <- rep(letters[1:6], 6)
  adj <- matrix(0, 6, 6, dimnames = list(letters[1:6], letters[1:6]))
  adj[cbind(1:5, 2:6)] <- 1; adj <- adj + t(adj); graph <- car_graph(adj)
  forms <- list(y ~ x + car(area, graph = graph), y ~ x + dagar(area, graph = graph))
  for(form in forms){
    warnings <- character()
    f <- withCallingHandlers(stLMM(form, d, family = "nb", n_samples = 20, chains = 2,
      priors = list(size = gamma_dist(2, .5),
        car_1 = list(sigma_sq = ig(3, .5), rho = uniform(.05, .8)),
        dagar_1 = list(sigma_sq = ig(3, .5), rho = uniform(.05, .8))), verbose = FALSE),
      warning = function(w){warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")})
    expect_length(warnings, 1L)
    expect_match(warnings, "weakly identified")
    expect_true(all(is.finite(log_lik(f))))
    expect_true(all(is.finite(predict(f)$chains[[1]]$mu_samples)))
  }
})

test_that("new-group prediction uses the selected size draws after simulating group effects", {
  d <- nb_size_data()
  f <- suppressWarnings(stLMM(y ~ x + iid(g, ~ 1 + x, by = parent), d, family = "nb",
    priors = list(size = gamma_dist(2, .5), iid_1 = list(covariance = iw(4, diag(2)))),
    n_samples = 40, verbose = FALSE))
  nd <- d[1:4, ]; nd$g <- "new"; nd$y <- NULL
  expect_error(predict(f, nd), "allow_new_levels")
  set.seed(59130)
  mean <- predict(f, nd, allow_new_levels = TRUE, sub_sample = list(start = 5, thin = 4))
  expected <- t(vapply(seq_along(mean$draw_index), function(i)
    rnbinom(nrow(nd), size = f$size_samples[mean$draw_index[i]], mu = mean$mu_samples[i, ]), numeric(nrow(nd))))
  set.seed(59130)
  actual <- predict(f, nd, allow_new_levels = TRUE, y_samples = TRUE,
    sub_sample = list(start = 5, thin = 4))
  expect_equal(unname(actual$y_samples), unname(expected))
})

test_that("extreme PG underflow is reported without clamping size or its weights", {
  expect_error(stLMM(y ~ 0 + offset(eta), data.frame(y = 0, eta = 1), family = "nb",
    priors = list(size = gamma_dist(2, 1)), starting = list(size = 1e-300),
    warmup = FALSE, n_samples = 1, verbose = FALSE), "Polya-Gamma.*shape=.*tilt=")
})
