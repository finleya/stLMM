test_that("prediction ignores response columns across likelihoods and chain types", {
  set.seed(84101)
  d <- data.frame(y = rnorm(12), x = seq(-1, 1, length.out = 12),
    category = factor(rep(c("a", "b"), 6)), off = seq(-.2, .2, length.out = 12),
    group = rep(c("g1", "g2", "g3"), each = 4), time = 1:12, trials = 3L)
  prediction_draws <- function(object, newdata, ...){
    set.seed(84102)
    p <- predict(object, newdata, y_samples = TRUE, ...)
    list(draws = if(is.null(p$chains)) list(p) else p$chains, rng = .Random.seed)
  }
  for(family in c("gaussian", "binomial", "probit", "negative_binomial")){
    dat <- d
    if(family != "gaussian") dat$y <- rep(c(0, 1), 6)
    dat$y[3L] <- NA_real_
    for(structured in c(FALSE, TRUE)){
      formula <- y ~ x * category + offset(off) + iid(group, ~ 1 + x)
      starting <- if(family == "gaussian") list(resid = list(tau_sq = fixed(.4))) else list()
      if(structured){
        formula <- update(formula, . ~ . + ar1(time))
        starting$ar1_1 <- list(sigma_sq = fixed(.3), phi = fixed(.4))
      }
      fit <- stLMM(formula, dat, family = family,
        trials = if(family == "binomial") dat$trials else NULL,
        size = if(family == "negative_binomial") 2 else NULL,
        priors = list(beta = normal(0, 2), iid_1 = list(covariance = iw(4, diag(2)))),
        starting = starting, n_samples = 12, chains = 2, warmup = FALSE, verbose = FALSE)
      if(structured) fit <- recover(fit, sub_sample = list(start = 2, thin = 2))
      nd <- dat[c(2, 5, 9), setdiff(names(dat), "y"), drop = FALSE]
      for(object in list(fit$chains[[1L]], fit)){
        expected <- prediction_draws(object, nd)
        for(response in list(c(-1e6, 0, 1e6), c(NA_real_, 1, NA_real_), rep(NA_real_, 3))){
          supplied <- nd
          supplied$y <- response
          expect_identical(prediction_draws(object, supplied), expected,
            info = paste(family, structured, class(object)[1L]))
        }
        # New iid groups and new AR1 times retain the same predictive draws.
        new <- nd
        new$group <- "new_group"
        if(structured) new$time <- c(14L, 16L, 14L)
        expected_new <- prediction_draws(object, new, allow_new_levels = TRUE)
        new$y <- NA_real_
        expect_identical(prediction_draws(object, new, allow_new_levels = TRUE), expected_new)
      }
    }
  }
})

test_that("prediction does not evaluate transformed responses", {
  d <- data.frame(y = exp(seq(-.3, .6, length.out = 8)), x = 1:8,
                  off = seq(-.2, .2, length.out = 8))
  fit <- stLMM(log(y) ~ x + offset(off), d, n_samples = 10, verbose = FALSE)
  nd <- d[c(2, 4, 7), c("x", "off"), drop = FALSE]
  expected <- predict(fit, nd)
  expect_identical(predict(fit, transform(nd, y = NA_real_)), expected)
  expect_identical(predict(fit, transform(nd, y = -1)), expected)
  manual <- sweep(fit$beta_samples %*% t(cbind(1, nd$x)), 2L, nd$off, `+`)
  expect_equal(as.numeric(expected$mu_samples), as.numeric(manual))
})

test_that("response-free prediction preserves rows with empty fixed-effect designs", {
  d <- data.frame(y = c(.2, -.1, .3, .4, -.2, .1), group = rep(c("a", "b"), 3))
  for(formula in list(y ~ 1, y ~ 0 + iid(group))){
    fit <- stLMM(formula, d, n_samples = 10, verbose = FALSE,
      priors = if(length(attr(terms(formula), "term.labels")))
        list(iid_1 = list(sigma_sq = ig(2, 1))) else list())
    nd <- data.frame(group = c("b", "a", "b"))
    expected <- predict(fit, nd)
    expect_equal(dim(expected$mu_samples), c(10L, 3L))
    expect_identical(predict(fit, transform(nd, y = NA_real_)), expected)
  }
})

test_that("ignoring responses preserves required predictor and offset checks", {
  d <- data.frame(y = seq(-.4, .7, length.out = 8), x = 1:8,
                  off = seq(-.2, .2, length.out = 8), group = rep(c("a", "b"), 4))
  fit <- stLMM(y ~ x + offset(off) + iid(group), d, n_samples = 10, verbose = FALSE,
    priors = list(iid_1 = list(sigma_sq = ig(2, 1))))
  nd <- transform(d[1:3, ], y = NA_real_)
  expect_error(predict(fit, transform(nd, x = NA_real_)), "missing values")
  expect_error(predict(fit, nd[, setdiff(names(nd), "x"), drop = FALSE]), "x")
  expect_error(predict(fit, transform(nd, off = NA_real_)), "missing values")
  expect_error(predict(fit, transform(nd, off = Inf)), "offset values must be finite")
  expect_error(predict(fit, nd[, setdiff(names(nd), "off"), drop = FALSE]), "off")
  expect_error(predict(fit, transform(nd, group = NA_character_)), "group.*missing")
})

test_that("a response variable explicitly used in an offset remains an input", {
  d <- data.frame(y = seq(1, 2, length.out = 8), x = 1:8)
  fit <- stLMM(y ~ x + offset(log(y)), d, n_samples = 10, verbose = FALSE)
  nd <- data.frame(x = c(2, 4), y = c(1, 2))
  a <- predict(fit, nd)
  b <- predict(fit, transform(nd, y = 2 * y))
  expect_equal(as.numeric(b$mu_samples - a$mu_samples), rep(log(2), length(a$mu_samples)))
  expect_error(predict(fit, transform(nd, y = NA_real_)), "missing values")
  expect_error(predict(fit, nd["x"]), "y")
})
