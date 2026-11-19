iid_by_data <- function(){
  set.seed(81230)
  d <- data.frame(parent = rep(c("A", "A", "B", "B", "B"), each = 5),
                  plot = rep(c("1", "2", "1", "3", "4"), each = 5),
                  x = rnorm(25), y = rnorm(25))
  d$parent <- factor(d$parent, levels = c("B", "A", "unused_parent"))
  d$plot <- factor(d$plot, levels = c("1", "2", "3", "4", "unused_child"))
  d[sample(nrow(d)), ]
}

test_that("by terms create sparse local-child blocks with separate covariances", {
  d <- iid_by_data()
  f <- stLMM(y ~ x + iid(parent) + iid(plot, ~ 1 + x, by = parent), d,
    priors = list(iid_1 = list(sigma_sq = ig(2, 1)),
                  iid_2 = list(covariance = iw(4, diag(2)))),
    n_samples = 30, verbose = FALSE)
  # Independent reference columns: parent intercepts, then child intercept/slope
  # pairs within B and A. Unobserved parent/child combinations are not latent groups.
  Z <- sapply(levels(d$parent), function(p) as.numeric(d$parent == p))
  for(p in c("B", "A")) for(g in if(p == "B") c("1", "3", "4") else c("1", "2")){
    active <- as.numeric(d$parent == p & d$plot == g)
    Z <- cbind(Z, active, active * d$x)
  }
  expect_s4_class(f$backend$Z, "dgCMatrix")
  expect_equal(unname(as.matrix(f$backend$Z)), unname(Z))
  expect_length(f$backend$Z@x, 3 * nrow(d))
  expect_equal(ncol(f$alpha_samples), 13L)
  terms <- f$term_description$random_effects$terms
  expect_equal(vapply(terms, `[[`, integer(1), "n_levels"), c(3L, 3L, 2L), ignore_attr = TRUE)
  expect_equal(terms[[2]]$by_levels, c("B", "A"))
  expect_equal(terms[[2]]$formula_term, "iid_2")
  expect_equal(terms[[3]]$levels, c("1", "2"))
  expect_equal(dim(f$iid_covariance_samples), c(30L, 6L))
  expect_equal(colnames(f$iid_sigma_sq_samples), "iid_1_sigma_sq")
  expect_match(colnames(f$iid_covariance_samples)[1], 'iid_2[parent="B"]', fixed = TRUE)
  expect_match(colnames(f$iid_covariance_samples)[4], 'iid_2[parent="A"]', fixed = TRUE)
  expected <- f$beta_samples %*% t(f$backend$X) + f$alpha_samples %*% t(Z)
  expect_equal(unname(predict(f, d)$mu_samples), unname(as.matrix(expected)))
  expect_true(all(colnames(f$iid_covariance_samples) %in% names(as_samples(f))))
  expect_s3_class(summary(f), "summary_stLMM")
})

test_that("scalar by terms and external slopes use independent IG blocks", {
  d <- iid_by_data()
  f <- stLMM(y ~ x + iid(plot, by = parent) + x:iid(plot, by = parent), d,
    priors = list(iid_1 = list(sigma_sq = ig(3, 2)), iid_2 = list(sigma_sq = ig(4, 1))),
    starting = list(iid_1 = list(sigma_sq = .8), iid_2 = list(sigma_sq = .4)),
    n_samples = 15, chains = 2, chain_control = list(seed = 8123), verbose = FALSE)
  for(ch in f$chains){
    expect_equal(ncol(ch$alpha_samples), 10L)
    expect_equal(ncol(ch$iid_sigma_sq_samples), 4L)
    expect_null(ch$iid_covariance_samples)
    expect_equal(ch$backend$sigma_sq_re_starting, c(.8, .8, .4, .4))
    expect_equal(unname(ch$backend$sigma_sq_re_IG), rbind(c(3, 2), c(3, 2), c(4, 1), c(4, 1)))
  }
  expect_s3_class(predict(f, d), "stLMM_prediction_chains")
})

test_that("a single parent reduces exactly to the existing iid sampler", {
  d <- droplevels(subset(iid_by_data(), parent == "A"))
  for(vector in c(FALSE, TRUE)){
    prior <- if(vector) list(covariance = iw(4, diag(2))) else list(sigma_sq = ig(2, 1))
    forms <- if(vector) list(y ~ x + iid(plot, ~ 1 + x), y ~ x + iid(plot, ~ 1 + x, by = parent)) else
      list(y ~ x + iid(plot), y ~ x + iid(plot, by = parent))
    fits <- lapply(forms, function(form){
      set.seed(81231)
      f <- stLMM(form, d, priors = list(iid_1 = prior), n_samples = 40, verbose = FALSE)
      list(f = f, rng = .Random.seed)
    })
    a <- fits[[1]]$f; b <- fits[[2]]$f
    expect_identical(a$beta_samples, b$beta_samples)
    expect_equal(unname(a$alpha_samples), unname(b$alpha_samples), tolerance = 0)
    expect_equal(unname(a$iid_sigma_sq_samples), unname(b$iid_sigma_sq_samples), tolerance = 0)
    expect_equal(unname(a$iid_covariance_samples), unname(b$iid_covariance_samples), tolerance = 0)
    expect_identical(a$tau_sq_samples, b$tau_sq_samples)
    expect_identical(fits[[1]]$rng, fits[[2]]$rng)
  }
})

test_that("prediction matches parent-child pairs and rejects new parents", {
  d <- iid_by_data()
  for(vector in c(FALSE, TRUE)){
    f <- if(vector) stLMM(y ~ x + iid(plot, ~ 1 + x, by = parent), d,
      priors = list(iid_1 = list(covariance = iw(4, diag(2)))), n_samples = 100, verbose = FALSE) else
      stLMM(y ~ x + iid(plot, by = parent), d,
      priors = list(iid_1 = list(sigma_sq = ig(2, 1))), n_samples = 100, verbose = FALSE)
    nd <- data.frame(parent = c("A", "A", "A", "B", "B"),
                     plot = c("3", "3", "3", "3", "3"), x = c(0, 0, 1, 0, 0))
    # Plot 3 is new in A, but fitted in B.
    expect_error(predict(f, nd), "allow_new_levels = TRUE", fixed = TRUE)
    set.seed(81232)
    p <- predict(f, nd, allow_new_levels = TRUE)
    expect_equal(p$mu_samples[, 1], p$mu_samples[, 2])
    expect_equal(p$mu_samples[, 4], p$mu_samples[, 5])
    expect_equal(unname(p$mu_samples[, 4]), unname(predict(f, nd[4, ])$mu_samples[, 1]))
    expect_gt(sd(p$mu_samples[, 1] - p$mu_samples[, 4]), 0)
    both_new <- transform(nd, plot = "new")
    pn <- predict(f, both_new, allow_new_levels = TRUE)
    expect_equal(pn$mu_samples[, 1], pn$mu_samples[, 2])
    expect_gt(sd(pn$mu_samples[, 1] - pn$mu_samples[, 4]), 0)
    expect_error(predict(f, transform(nd, parent = "new")), "allow_new_levels|new parent")
    expect_error(predict(f, transform(nd, parent = "new"), allow_new_levels = TRUE), "new parent")
    expect_error(predict(f, transform(nd, parent = "unused_parent"), allow_new_levels = TRUE), "new parent")
    expect_error(predict(f, nd[, c("plot", "x")], allow_new_levels = TRUE), "parent.*missing")
    expect_error(predict(f, transform(nd, parent = NA_character_), allow_new_levels = TRUE), "parent.*missing")
    expect_error(predict(f, transform(nd, plot = NA_character_), allow_new_levels = TRUE), "missing values")
    if(vector) expect_error(predict(f, transform(nd, x = Inf), allow_new_levels = TRUE), "finite")
  }
})

test_that("by terms validate arguments and keep unambiguous parent labels", {
  d <- iid_by_data()
  expect_error(stLMM:::parse_iid_call('iid(plot, by = "parent")'), "simple variable")
  expect_error(stLMM:::parse_iid_call("iid(plot, by = interaction(parent, plot))"), "simple variable")
  expect_error(stLMM:::parse_iid_call("iid(plot, by = parent, by = plot)"), "duplicate")
  expect_error(stLMM:::parse_iid_call('iid(plot, by = parent, cov_model = "unstructured")'), "requires")
  expect_equal(stLMM:::parse_iid_call("iid(plot, by = parent, ~ 1 + x)")$coefficient,
                c("(Intercept)", "x"))
  expect_equal(stLMM:::parse_iid_call("iid(coefficients = ~ 1 + x, by = parent, plot)")$group_name, "plot")
  d$parent <- ifelse(d$parent == "A", 'A.[]"', 'A.[]\\"')
  f <- stLMM(y ~ x + iid(plot, ~ 1 + x, by = parent), d,
    priors = list(iid_1 = list(covariance = iw(4, diag(2)))), n_samples = 5, verbose = FALSE)
  expect_false(anyDuplicated(colnames(f$alpha_samples)) > 0)
  expect_false(anyDuplicated(colnames(f$iid_covariance_samples)) > 0)
  expected <- f$beta_samples %*% t(f$backend$X) + f$alpha_samples %*% t(as.matrix(f$backend$Z))
  expect_equal(unname(predict(f, d)$mu_samples), unname(as.matrix(expected)))
  bad <- d
  bad$parent <- I(cbind(d$parent, d$parent))
  expect_error(predict(f, bad), "must be a vector")
})

test_that("by terms work across chains, likelihoods, recovery, and missing responses", {
  for(family in c("gaussian", "binomial", "probit", "negative_binomial")){
    d <- iid_by_data()
    if(family != "gaussian") d$y <- as.numeric(d$y > 0)
    d$y[1] <- NA_real_
    d$time <- seq_len(nrow(d))
    f <- stLMM(y ~ x + iid(plot, ~ 1 + x, by = parent) + ar1(time), d,
      family = family, trials = if(family == "binomial") rep(1, nrow(d)) else NULL,
      size = if(family == "negative_binomial") 2 else NULL,
      priors = list(iid_1 = list(covariance = iw(4, diag(2)))),
      starting = list(iid_1 = list(covariance = diag(2)),
                      ar1_1 = list(sigma_sq = fixed(.2), phi = fixed(.4))),
      n_samples = 12, chains = 2, chain_control = list(seed = 81234), verbose = FALSE)
    for(ch in f$chains){
      expect_equal(ch$backend$iid_covariance[[1]]$starting, diag(2))
      expect_equal(ch$backend$iid_covariance[[2]]$starting, diag(2))
      expect_true(all(is.finite(ch$iid_covariance_samples)))
    }
    nd <- data.frame(parent = c("A", "A", "B", "B"), plot = "new", x = 0, time = 2)
    objects <- if(family == "gaussian") list(recover(f)) else list(f, recover(f))
    for(object in objects){
      p <- predict(object, nd, allow_new_levels = TRUE, scale = "link",
                   sub_sample = list(start = 2, thin = 2))
      expect_length(p$chains, 2L)
      for(ch in p$chains){
        expect_equal(ch$draw_index, seq(2, 12, 2))
        expect_equal(ch$mu_samples[, 1], ch$mu_samples[, 2])
        expect_equal(ch$mu_samples[, 3], ch$mu_samples[, 4])
      }
    }
  }
})
