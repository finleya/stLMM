combined_summary_data <- function(){
  set.seed(46201)
  d <- data.frame(parent = rep(c("B", "A"), each = 12),
                  plot = rep(rep(c("1", "2", "3"), each = 4), 2),
                  crossed = rep(c("east", "west"), 12),
                  x = rep(c(0, -1, 0, 2), 6), y = rnorm(24))
  d$parent <- factor(d$parent, levels = c("A", "B", "unused"))
  d$plot <- factor(d$plot, levels = c("3", "1", "2", "unused"))
  d[sample(nrow(d)), ]
}

combined_summary_fit <- function(d = combined_summary_data(), chains = 3L, ...){
  stLMM(y ~ x + iid(parent) + iid(plot, ~ 1 + x, by = parent), d,
         priors = list(iid_1 = list(sigma_sq = ig(2, 1)),
                       iid_2 = list(covariance = iw(4, diag(2)))),
         n_samples = 40, chains = chains, chain_control = list(seed = 46202),
         verbose = FALSE, ...)
}

test_that("combined iid coefficients reconstruct fitted predictors without weighting slopes", {
  d <- combined_summary_data()
  f <- combined_summary_fit(d)
  rng <- .Random.seed
  s <- summary(f, type = "combined", terms = "all", burn = 7,
               diagnostics = FALSE)
  expect_s3_class(s, "summary_stLMM_combined")
  expect_identical(.Random.seed, rng)
  expect_identical(s$iteration, 8:40)
  expect_identical(s$chain_ids, 1:3)
  expect_equal(s$n_used, 33)
  expect_equal(nrow(s$parameters), 12L)
  expect_true(all(s$quantities$fixed_effect))
  expect_false(any(s$quantities$group_parent == "unused"))
  expect_equal(unique(s$quantities$terms[s$quantities$coefficient == "x"]), "iid_2")
  for(i in seq_len(nrow(s$quantities))){
    info <- s$quantities[i, ]
    # Independent reference via predictions at x = 0 and 1; differencing gives
    # the coefficient of x. All quantities use the original paired draws.
    nd <- d[info$data_row, , drop = FALSE]
    nd$x <- 0
    zero <- predict(f, nd, sub_sample = list(start = 8))
    nd$x <- 1
    one <- predict(f, nd, sub_sample = list(start = 8))
    draws <- unlist(lapply(seq_along(f$chains), function(k){
      if(info$coefficient == "(Intercept)") zero$chains[[k]]$mu_samples[, 1] else
        one$chains[[k]]$mu_samples[, 1] - zero$chains[[k]]$mu_samples[, 1]
    }))
    expected <- c(mean(draws), sd(draws), quantile(draws, c(.025, .5, .975), names = FALSE))
    expect_equal(unname(s$parameters[i, ]), expected, tolerance = 1e-12)
    expect_identical(as.character(d$parent[info$data_row]), info$group_parent)
    expect_identical(as.character(d$plot[info$data_row]), info$group_plot)
  }
  parent <- summary(f, type = "combined", coefficients = "(Intercept)",
                     terms = "iid_1", diagnostics = FALSE)
  expect_equal(nrow(parent$parameters), 2L)
  expect_false("group_plot" %in% names(parent$quantities))
  expect_output(print(s), "no averaging")
})

test_that("crossed terms include only represented combinations and handle label punctuation", {
  d <- combined_summary_data()
  d$parent <- ifelse(d$parent == "A", 'A:[]"', 'A:[]\\"')
  d$crossed[d$parent == 'A:[]"'] <- "east"
  f <- stLMM(y ~ iid(parent) + iid(crossed), d,
    priors = list(iid_1 = list(sigma_sq = ig(2, 1)), iid_2 = list(sigma_sq = ig(2, 1))),
    n_samples = 15, verbose = FALSE)
  s <- summary(f, type = "combined", terms = "all")
  expect_equal(nrow(s$parameters), 3L)
  expect_equal(nrow(unique(s$quantities[c("group_parent", "group_crossed")])), 3L)
  expected <- fitted(f, summary = FALSE)
  expect_equal(unname(s$parameters[, "mean"]), unname(colMeans(expected)[s$quantities$data_row]))
  expect_equal(nrow(summary(f, type = "combined", terms = "iid_2")$parameters), 2L)
})

test_that("zero slope covariates retain membership and absent fixed coefficients contribute zero", {
  d <- combined_summary_data()
  d$x <- 0
  f <- stLMM(y ~ 1 + x:iid(plot, by = parent), d,
    priors = list(iid_1 = list(sigma_sq = ig(2, 1))), n_samples = 20, verbose = FALSE)
  s <- summary(f, type = "combined", terms = "iid_1", coefficients = "x")
  expect_equal(nrow(s$parameters), 6L)
  expect_false(any(s$quantities$fixed_effect))
  expect_true(all(grepl("^0 \\+ alpha", s$quantities$expression)))
  for(i in seq_len(nrow(s$parameters))){
    nd <- d[s$quantities$data_row[i], , drop = FALSE]
    nd$x <- 1
    draws <- predict(f, nd)$mu_samples[, 1] - f$beta_samples[, "(Intercept)"]
    expect_equal(unname(s$parameters[i, "mean"]), mean(draws), tolerance = 1e-12)
  }
  # Damaging the zero entries destroys membership information; do not guess.
  f$backend$Z <- Matrix::drop0(f$backend$Z)
  expect_error(summary(f, type = "combined", terms = "all"), "missing group membership")
})

test_that("rank diagnostics match posterior on transformed draws and retain chain identity", {
  skip_if_not_installed("posterior")
  f <- combined_summary_fit()
  s <- summary(f, type = "combined", coefficients = "(Intercept)",
               terms = "iid_1", burn = 5, diagnostic_method = "rank")
  for(i in seq_len(nrow(s$parameters))){
    level <- match(s$quantities$group_parent[i],
                    f$chains[[1]]$term_description$random_effects$terms[[1]]$levels)
    draws <- vapply(f$chains, function(ch)
      ch$beta_samples[6:40, "(Intercept)"] + ch$alpha_samples[6:40, level], numeric(35))
    expect_equal(s$diagnostics$rhat[i], posterior::rhat(draws))
    expect_equal(s$diagnostics$ess_bulk[i], posterior::ess_bulk(draws))
    expect_equal(s$diagnostics$ess_tail[i], posterior::ess_tail(draws))
  }
  ordinary <- summary(f, parameters = c("(Intercept)", "x"), burn = 5, diagnostic_method = "rank")
  draws <- vapply(f$chains, function(ch) ch$beta_samples[6:40, "(Intercept)"], numeric(35))
  expect_equal(ordinary$diagnostics$rhat[1], posterior::rhat(draws))
  expect_equal(ordinary$diagnostics$ess_bulk[1], posterior::ess_bulk(draws))
  expect_identical(ordinary$diagnostic_method, "rank")
  one <- summary(f$chains[[1]], type = "combined", terms = "iid_1",
                 diagnostics = TRUE, diagnostic_method = "rank")
  expect_true(all(is.na(one$diagnostics$rhat)))
  expect_true(all(is.finite(one$diagnostics$ess_bulk)))
  short <- summary(f, type = "combined", terms = "iid_1", burn = 39, diagnostic_method = "rank")
  expect_true(all(is.na(short$diagnostics$rhat)))
  constant <- matrix(1, 20, 3)
  expect_true(all(is.na(stLMM:::summary_draw_diagnostics(constant, "rank"))))
})

test_that("modern diagnostics detect scale differences and nonstationarity", {
  skip_if_not_installed("posterior")
  set.seed(46203)
  draws <- matrix(rnorm(8000), 2000, 4)
  scale <- sweep(draws, 2, c(1, 1, 4, 4), "*")
  drift <- draws / 4 + seq(-3, 3, length.out = 2000)
  expect_gt(stLMM:::summary_draw_diagnostics(scale, "rank")[["rhat"]], 1.15)
  expect_gt(stLMM:::summary_draw_diagnostics(drift, "rank")[["rhat"]], 1.5)
  draws[, 1L] <- 0
  expect_true(all(is.na(stLMM:::summary_draw_diagnostics(draws, "rank"))))
  draws[1L, 2L] <- Inf
  expect_true(all(is.na(stLMM:::summary_draw_diagnostics(draws, "rank"))))
})

test_that("combined coefficients exclude offsets and retain groups at missing responses", {
  d <- combined_summary_data()
  d$y[1L] <- NA
  d$o <- seq_len(nrow(d)) / 10
  f <- stLMM(y ~ 0 + x + iid(parent, ~ 1 + x) + offset(o), d,
    priors = list(iid_1 = list(covariance = iw(4, diag(2)))),
    n_samples = 30, verbose = FALSE)
  s <- summary(f, type = "combined", terms = "all", probs = .5)
  expect_equal(nrow(s$parameters), 4L)
  expect_equal(ncol(s$parameters), 3L)
  expect_false(any(s$quantities$fixed_effect[s$quantities$coefficient == "(Intercept)"]))
  for(i in seq_len(nrow(s$quantities))){
    nd <- d[s$quantities$data_row[i], , drop = FALSE]
    nd$x <- 0
    zero <- predict(f, nd)$mu_samples[, 1L]
    nd$x <- 1
    one <- predict(f, nd)$mu_samples[, 1L]
    draws <- if(s$quantities$coefficient[i] == "(Intercept)") zero - nd$o else one - zero
    expect_equal(unname(s$parameters[i, ]), c(mean(draws), sd(draws), median(draws)), tolerance = 1e-12)
  }
})

test_that("ordinary rank summaries handle single-chain and recovered process draws", {
  skip_if_not_installed("posterior")
  d <- combined_summary_data()
  d$time <- seq_len(nrow(d))
  f <- stLMM(y ~ x + ar1(time), d, n_samples = 400, chains = 2,
    starting = list(ar1_1 = list(sigma_sq = fixed(.2), phi = fixed(.4))), verbose = FALSE)
  r <- recover(f, sub_sample = list(start = 3, thin = 2))
  for(object in list(r, r$chains[[1L]])){
    s <- summary(object, include_w = TRUE, diagnostics = TRUE, diagnostic_method = "rank", burn = 4)
    expect_true(all(c("rhat", "ess_bulk", "ess_tail") %in% names(s$diagnostics)))
    expect_identical(s$diagnostic_method, "rank")
    if(!inherits(object, "stLMM_recovery_chains")) expect_true(all(is.na(s$diagnostics$rhat)))
  }
  s <- summary(r, include_w = TRUE, diagnostic_method = "rank", burn = 4)
  draws <- vapply(r$chains, function(ch) ch$beta_samples[ch$recover_iter[ch$recover_iter > 4], "x"], numeric(198))
  expect_equal(s$diagnostics["x", "rhat"], posterior::rhat(draws))
  expect_equal(s$diagnostics["x", "ess_bulk"], posterior::ess_bulk(draws))
})

test_that("combined summaries reject ambiguous selectors and unsupported structured contributions", {
  f <- combined_summary_fit(chains = 1L)
  expect_error(summary(f, type = "combined"), "explicit terms")
  expect_error(summary(f, terms = "iid_1"), 'require type = "combined"', fixed = TRUE)
  expect_error(summary(f, type = "combined", terms = c("iid_1", "iid_1")), "without duplicates")
  expect_error(summary(f, type = "combined", terms = "iid_7"), "unknown iid")
  expect_error(summary(f, type = "combined", terms = "iid_1", coefficients = "x"), "no selected iid")
  expect_error(summary(f, type = "combined", terms = "all", parameters = "x"), "use coefficients")
  expect_error(summary(f, type = "combined", terms = "all", burn = 40), "removes all")
  d <- combined_summary_data()
  d$time <- seq_len(nrow(d))
  for(family in c("gaussian", "binomial", "probit", "negative_binomial")){
    if(family != "gaussian") d$y <- as.numeric(d$y > 0)
    fit <- stLMM(y ~ x + iid(parent) + ar1(time), d, family = family,
      priors = list(iid_1 = list(sigma_sq = ig(2, 1))),
      trials = if(family == "binomial") rep(1, nrow(d)) else NULL,
      size = if(family == "negative_binomial") 2 else NULL,
      starting = list(ar1_1 = list(sigma_sq = fixed(.2), phi = fixed(.4))),
      n_samples = 12, chains = 2, verbose = FALSE)
    expect_error(summary(fit, type = "combined", terms = "all"), "structured")
    expect_error(summary(fit, type = "combined", terms = "ar1_1"), "structured")
    a <- summary(fit, type = "combined", terms = "iid_1", diagnostics = FALSE, burn = 2)
    rec <- recover(fit, sub_sample = list(start = 4, thin = 3))
    b <- summary(rec, type = "combined", terms = "iid_1", diagnostics = FALSE, burn = 2)
    expect_identical(a, b) # iid sums use all retained coefficient draws, not the recovery subset
    expect_equal(a$scale, "link coefficient")
    expect_error(summary(rec, type = "combined", terms = "iid_1", include_w = TRUE), "include_w")
    expect_s3_class(summary(rec$chains[[1]], type = "combined", terms = "iid_1"), "summary_stLMM_combined")
  }
})
