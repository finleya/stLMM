test_that("adaptive Metropolis output has exactly its documented fields", {
  set.seed(86101)
  d <- data.frame(y = rnorm(24), x = rnorm(24), g = rep(letters[1:4], each = 6), time = 1:24)
  expected <- c("dimension", "batch_length", "batch_index", "target_accept",
    "last_batch_accept", "sigma_sq_m", "Sigma0", "proposal_cov", "current_eta",
    "parameter_labels", "acceptance", "batch_acceptance_history",
    "proposal_scale_history", "warmup", "blocks")
  forms <- list(y ~ x, y ~ x + iid(g, ~ 1 + x), y ~ x + iid(g, ~ 1 + x) + ar1(time))
  for(form in forms){
    f <- stLMM(form, d, n_samples = 8, verbose = FALSE,
      priors = list(iid_1 = list(covariance = iw(4, diag(2))),
                    ar1_1 = list(sigma_sq = ig(3, 1), phi = uniform(-.9, .9))))
    expect_identical(names(f$adaptive_metropolis), expected)
    expect_length(f$adaptive_metropolis, length(expected))
  }
})

test_that("fitted polynomial and spline bases survive prediction subsets and recovery", {
  set.seed(86102)
  d <- data.frame(y = rnorm(36), x = rnorm(36), z = rnorm(36),
    f = factor(rep(letters[1:3], 12)), off = seq(-.2, .2, length.out = 36),
    g = rep(letters[1:6], each = 6), time = 1:36)
  for(family in c("gaussian", "binomial", "probit", "negative_binomial")){
    dat <- d
    if(family != "gaussian") dat$y <- as.numeric(dat$y > 0)
    dat$y[2] <- NA_real_
    for(structured in c(FALSE, TRUE)){
      form <- y ~ poly(x, 2) + splines::ns(z, df = 3) + f + offset(off) + iid(g, ~ 1 + x)
      starting <- list()
      if(structured){
        form <- update(form, . ~ . + ar1(time))
        starting$ar1_1 <- list(sigma_sq = fixed(.3), phi = fixed(.4))
      }
      fit <- stLMM(form, dat, family = family, size = if(family == "negative_binomial") 2 else NULL,
        priors = list(beta = normal(0, 2), iid_1 = list(covariance = iw(4, diag(2)))),
        starting = starting, n_samples = 12, chains = 2, warmup = FALSE, verbose = FALSE)
      if(structured) fit <- recover(fit, sub_sample = list(start = 2, thin = 2))
      for(object in list(fit$chains[[1]], fit)){
        all <- predict(object, scale = "link")
        for(rows in list(c(1L, 2L, 9L, 9L), 2L)){
          nd <- dat[rows, setdiff(names(dat), "y"), drop = FALSE]
          nd$f <- as.character(nd$f)
          old <- options(contrasts = c("contr.sum", "contr.poly"))
          pr <- tryCatch(predict(object, nd, scale = "link"), finally = options(old))
          actual <- if(is.null(pr$chains)) list(pr) else pr$chains
          expected <- if(is.null(all$chains)) list(all) else all$chains
          for(ch in seq_along(actual))
            expect_equal(unname(actual[[ch]]$mu_samples),
              unname(expected[[ch]]$mu_samples[, rows, drop = FALSE]), tolerance = 1e-12)
        }
      }
    }
  }
})

test_that("basis prediction evaluates new values in the fitted basis", {
  set.seed(86103)
  d <- data.frame(y = rnorm(24), x = seq(-1, 1, length.out = 24), z = rnorm(24))
  f <- stLMM(y ~ poly(x, 2) + splines::ns(z, df = 3), d, n_samples = 10, verbose = FALSE)
  nd <- data.frame(x = c(-2, .25, 1.5), z = c(-1.5, 0, 2))
  X <- cbind(1, predict(poly(d$x, 2), nd$x), predict(splines::ns(d$z, df = 3), nd$z))
  expect_equal(unname(predict(f, nd)$mu_samples), unname(f$beta_samples %*% t(X)), tolerance = 1e-12)
})

test_that("fixed factor levels and custom contrasts are retained", {
  set.seed(86104)
  d <- data.frame(y = rnorm(18), f = factor(rep(letters[1:3], 6)))
  contrasts(d$f) <- contr.sum(3)
  f <- stLMM(y ~ f, d, n_samples = 10, verbose = FALSE)
  nd <- data.frame(f = "b")
  expect_equal(unname(predict(f, nd)$mu_samples),
    unname(f$beta_samples %*% t(f$backend$X[2, , drop = FALSE])))
  expect_error(predict(f, data.frame(f = "new"), allow_new_levels = TRUE), "new level")
  d$f <- as.character(d$f)
  f <- stLMM(y ~ f, d, n_samples = 10, verbose = FALSE)
  expect_equal(unname(predict(f, nd)$mu_samples),
    unname(f$beta_samples %*% t(f$backend$X[2, , drop = FALSE])))
})

test_that("older saved fits retain ordinary numeric prediction", {
  d <- data.frame(y = seq(-.4, .7, length.out = 12), x = 1:12)
  f <- stLMM(y ~ x, d, n_samples = 10, verbose = FALSE)
  expected <- predict(f, d[2:4, ])
  f$backend[c("fixed_terms", "fixed_xlevels", "fixed_contrasts")] <- NULL
  expect_identical(predict(f, d[2:4, ]), expected)
})

test_that("all iid forms accept one numeric column and reject multicolumn covariates", {
  set.seed(86105)
  d <- data.frame(y = rnorm(24), x = rnorm(24), g = rep(letters[1:4], 6),
                  parent = rep(c("A", "B"), each = 12))
  d$x <- as.numeric(scale(d$x))
  forms <- list(y ~ iid(g) + x:iid(g), y ~ iid(g, ~ 1 + x),
    y ~ iid(g, ~ 1 + x, cov_model = "independent"),
    y ~ iid(g, by = parent) + x:iid(g, by = parent),
    y ~ iid(g, ~ 1 + x, by = parent),
    y ~ iid(g, ~ 1 + x, cov_model = "independent", by = parent))
  for(i in seq_along(forms)){
    priors <- if(i %% 3L == 2L) list(iid_1 = list(covariance = iw(4, diag(2)))) else
      list(iid_1 = list(sigma_sq = ig(3, 1)), iid_2 = list(sigma_sq = ig(3, 1)))
    fit <- function(dat){
      set.seed(86106)
      stLMM(forms[[i]], dat, priors = priors, n_samples = 12, verbose = FALSE)
    }
    vector <- fit(d)
    matrix_data <- d; matrix_data$x <- matrix(d$x, ncol = 1)
    one_column <- fit(matrix_data)
    expect_identical(vector$samples, one_column$samples)
    expect_identical(vector$backend$Z, one_column$backend$Z)
    for(new in c(FALSE, TRUE)){
      nd <- d[1:4, ]; if(new) nd$g <- "new"
      set.seed(86107)
      expected <- predict(vector, nd, allow_new_levels = new, y_samples = TRUE)
      nd$x <- matrix(nd$x, ncol = 1)
      set.seed(86107)
      expect_identical(predict(one_column, nd, allow_new_levels = new, y_samples = TRUE), expected)
      nd$x <- cbind(nd$x, nd$x)
      expect_error(predict(vector, nd, allow_new_levels = new), "finite numeric vector or one-column")
      nd$x <- c(1, NA, 2, 3)
      expect_error(predict(vector, nd, allow_new_levels = new), "missing values")
    }
    bad <- d; bad$x <- cbind(d$x, d$x)
    expect_error(fit(bad), "finite numeric vector or one-column")
    bad$x <- d$x; bad$x[1] <- Inf
    expect_error(fit(bad), "finite numeric")
  }
})

test_that("independent iid starting lists separate coefficient and chain indices", {
  set.seed(86108)
  d <- data.frame(y = rnorm(24), x = rnorm(24), g = rep(letters[1:4], each = 6))
  f <- stLMM(y ~ x + iid(g, ~ 1 + x, cov_model = "independent"), d,
    priors = list(iid_1 = list(sigma_sq = ig(3, 1))),
    starting = list(iid_1 = list(sigma_sq = list("(Intercept)" = c(.4, .8), x = c(.2, .6)))),
    chains = 2, n_samples = 8, verbose = FALSE)
  expect_equal(f$chains[[1]]$backend$sigma_sq_re_starting, c(.4, .2))
  expect_equal(f$chains[[2]]$backend$sigma_sq_re_starting, c(.8, .6))
})
