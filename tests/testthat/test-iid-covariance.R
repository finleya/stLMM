iid_covariance_data <- function(){
  set.seed(9401)
  data.frame(y = rnorm(36), x = rnorm(36), z = rnorm(36),
             g = factor(rep(letters[1:6], each = 6)),
             h = factor(rep(LETTERS[1:3], 12)), time = rep(1:6, 6))
}

test_that("iid coefficient formulas build sparse group-major designs without scaling", {
  d <- iid_covariance_data()
  f <- stLMM(y ~ x + iid(g, ~ 1 + x + z), d,
             priors = list(iid_1 = list(covariance = iw(5, diag(3)))),
             n_samples = 8, verbose = FALSE)
  Z <- matrix(0, 36, 18)
  for(i in seq_len(36)) Z[i, (as.integer(d$g[i]) - 1L) * 3 + 1:3] <- c(1, d$x[i], d$z[i])
  expect_s4_class(f$backend$Z, "dgCMatrix")
  expect_equal(unname(as.matrix(f$backend$Z)), Z)
  expect_equal(dim(f$iid_covariance_samples), c(8L, 6L))
  expect_null(f$iid_sigma_sq_samples)
  expect_true(all(colnames(f$iid_covariance_samples) %in% names(as_samples(f))))
  expected <- f$beta_samples %*% t(f$backend$X) + f$alpha_samples %*% t(Z)
  expect_equal(unname(predict(f, newdata = d)$mu_samples), unname(as.matrix(expected)))
  expect_s3_class(summary(f), "summary_stLMM")
})

test_that("mixed scalar and multiple covariance terms retain separate distributions", {
  d <- iid_covariance_data()
  d$nested <- interaction(d$g, d$h, drop = TRUE)
  f <- stLMM(y ~ x + iid(g) + iid(h, coefficients = ~ 1 + x, cov_model = "unstructured") +
               iid(nested, ~ 0 + x + z), d,
             priors = list(iid_1 = list(sigma_sq = ig(2, 1)),
                           iid_2 = list(covariance = iw(4, diag(2))),
                           iid_3 = list(covariance = iw(5, diag(c(2, 1))))),
             n_samples = 12, chains = 2, chain_control = list(seed = 947), verbose = FALSE)
  expect_equal(colnames(f$chains[[1]]$iid_sigma_sq_samples), "iid_1_sigma_sq")
  expect_equal(ncol(f$chains[[1]]$iid_covariance_samples), 6L)
  for(ch in f$chains){
    for(term in ch$term_description$random_effects$terms){
      if(is.null(term$covariance)) next
      expect_true(all(eigen(term$covariance$current, symmetric = TRUE)$values > 0))
    }
  }
  expect_s3_class(predict(f, newdata = d), "stLMM_prediction_chains")
  expect_equal(nrow(as_samples(f)), 24L)
})

test_that("iid covariance inputs reject ambiguous or invalid specifications", {
  d <- iid_covariance_data()
  expect_error(iw(1, diag(2)), "greater than")
  expect_error(iw(4, matrix(c(1, 2, 2, 1), 2)), "positive definite")
  expect_error(iw(4, matrix(c(1, 0, 1, 1), 2)), "symmetric")
  fit <- function(formula, prior = iw(4, diag(2)), ...)
    stLMM(formula, d, priors = list(iid_1 = list(covariance = prior)), n_samples = 2, verbose = FALSE, ...)
  expect_error(fit(y ~ iid(g, ~ 1 + x), iw(5, diag(3))), "2 by 2")
  expect_error(fit(y ~ iid(g, ~ x * z)), "additive numeric")
  expect_error(fit(y ~ iid(g, ~ h)), "finite numeric")
  expect_error(fit(y ~ iid(g, ~ 0)), "at least one")
  expect_error(fit(y ~ x:iid(g, ~ 1 + z)), "external slope")
  expect_error(fit(y ~ iid(g, ~ 1 + x, cov_model = "diagonal")), "unstructured")
  expect_error(fit(y ~ iid(g, ~ 1 + x), starting = list(iid_1 = list(covariance = diag(-1, 2)))), "positive definite")
  expect_error(fit(y ~ iid(g, ~ 1 + x), tuning = list(iid_1 = 1)), "Gibbs")
  scale <- diag(2)
  dimnames(scale) <- list(c("x", "(Intercept)"), c("x", "(Intercept)"))
  expect_error(fit(y ~ iid(g, ~ 1 + x), iw(4, scale)), "coefficient order")
})

test_that("new iid levels require opt-in and share coefficient draws across rows", {
  d <- iid_covariance_data()
  for(correlated in c(FALSE, TRUE)){
    f <- if(correlated)
      stLMM(y ~ x + iid(g, ~ 1 + x), d,
            priors = list(iid_1 = list(covariance = iw(4, diag(2)))), n_samples = 200, verbose = FALSE) else
      stLMM(y ~ x + iid(g), d, priors = list(iid_1 = list(sigma_sq = ig(2, 1))), n_samples = 200, verbose = FALSE)
    nd <- data.frame(x = c(0, 0, 1, 0, 0), g = c("new", "new", "new", "other", "a"))
    expect_error(predict(f, nd), "allow_new_levels = TRUE")
    expect_error(predict(f, nd, allow_new_levels = NA), "TRUE or FALSE")
    expect_error(predict(f, transform(nd, g = NA_character_), allow_new_levels = TRUE), "missing values")
    set.seed(9472)
    pr <- predict(f, nd, allow_new_levels = TRUE)
    expect_equal(pr$mu_samples[, 1], pr$mu_samples[, 2])
    expect_gt(sd(pr$mu_samples[, 1] - pr$mu_samples[, 4]), 0)
    expect_equal(unname(pr$mu_samples[, 5]), unname(predict(f, nd[5, , drop = FALSE])$mu_samples[, 1]))
    set.seed(9472)
    py <- predict(f, nd, allow_new_levels = TRUE, y_samples = TRUE)
    expect_equal(py$mu_samples, pr$mu_samples)
    expect_false(identical(py$y_samples[, 1], py$y_samples[, 2]))
    if(!correlated)
      expect_equal(unname(pr$mu_samples[, 3] - pr$mu_samples[, 1]), unname(f$beta_samples[, "x"]))
  }
})

test_that("covariance terms work with recovery, missing responses, and augmented likelihoods", {
  d <- iid_covariance_data()
  for(family in c("gaussian", "binomial", "probit", "negative_binomial")){
    dd <- d
    dd$time <- seq_len(nrow(dd))
    if(family != "gaussian") dd$y <- as.numeric(dd$y > 0)
    dd$y[1] <- NA
    f <- stLMM(y ~ x + iid(g, ~ 1 + x) + ar1(time), dd, family = family,
                trials = if(family == "binomial") rep(1, nrow(dd)) else NULL,
                size = if(family == "negative_binomial") 2 else NULL,
                starting = list(ar1_1 = list(sigma_sq = fixed(.3), phi = fixed(.4))),
                priors = list(iid_1 = list(covariance = iw(4, diag(2)))),
                n_samples = 10, verbose = FALSE)
    r <- recover(f, sub_sample = list(start = 3, thin = 2))
    nd <- data.frame(x = c(0, 0), g = c("new", "new"), time = c(2, 2))
    expect_error(predict(r, nd), "allow_new_levels")
    pr <- predict(r, nd, allow_new_levels = TRUE, scale = "link")
    expect_equal(pr$mu_samples[, 1], pr$mu_samples[, 2])
    expect_true(all(is.finite(pr$mu_samples)))
  }
})

test_that("retained unused levels and collinear coefficient columns have proper priors", {
  d <- iid_covariance_data()
  d$g <- factor(d$g, levels = c(levels(d$g), "unused"))
  d$z <- d$x
  f <- stLMM(y ~ x + iid(g, ~ 1 + x + z), d,
    priors = list(iid_1 = list(covariance = iw(5, diag(3)))),
    n_samples = 20, verbose = FALSE)
  expect_equal(ncol(f$alpha_samples), 21L)
  expect_true(all(is.finite(f$alpha_samples)))
  nd <- data.frame(x = .5, z = .5, g = "unused")
  # Retained factor levels already have coefficient draws, even without observations.
  pr <- predict(f, nd, sub_sample = list(start = 4, thin = 3))
  idx <- seq(4, 20, 3)
  expected <- f$beta_samples[idx, 1] + .5 * f$beta_samples[idx, 2] +
    drop(f$alpha_samples[idx, 19:21, drop = FALSE] %*% c(1, .5, .5))
  expect_equal(unname(pr$mu_samples[, 1]), unname(expected))
  expect_error(predict(f, transform(nd, x = Inf)), "finite")
})

test_that("new-level prediction forwards through chain and saved-process methods", {
  d <- iid_covariance_data()
  d$y <- as.numeric(d$y > 0)
  d$time <- seq_len(nrow(d))
  f <- stLMM(y ~ x + iid(g, ~ 1 + x) + ar1(time), d, family = "probit",
    starting = list(ar1_1 = list(sigma_sq = fixed(.2), phi = fixed(.3))),
    priors = list(iid_1 = list(covariance = iw(4, diag(2)))),
    n_samples = 8, chains = 2, verbose = FALSE)
  nd <- data.frame(x = c(0, 0), g = c("new", "new"), time = c(1, 1))
  for(object in list(f, recover(f))){
    expect_error(predict(object, nd), "allow_new_levels")
    pr <- predict(object, nd, allow_new_levels = TRUE, sub_sample = list(start = 2, thin = 2))
    expect_length(pr$chains, 2)
    for(ch in pr$chains){
      expect_equal(ch$draw_index, c(2L, 4L, 6L, 8L))
      expect_equal(ch$mu_samples[, 1], ch$mu_samples[, 2])
    }
  }
})
