test_that("AR1 new-node prediction matches the conditional moment", {
  n_draw <- 8000L
  sigma_sq <- 2
  phi <- 0.5
  w_fit <- matrix(rep(c(1, 3), each = n_draw), nrow = n_draw)

  object <- list(
    sigma_sq_samples = matrix(sigma_sq, nrow = n_draw, ncol = 1,
                              dimnames = list(NULL, "ar1_1")),
    theta_samples = matrix(phi, nrow = n_draw, ncol = 1,
                           dimnames = list(NULL, "ar1_1_phi")),
    w_samples_ordered = list(ar1_1 = w_fit)
  )
  term <- list(name = "ar1_1")
  graph <- list(support = c(0, 2))

  set.seed(4101)
  draw <- stLMM:::simulate_ar1_prediction_nodes(
    object = object,
    term = term,
    graph = graph,
    recover_row = seq_len(n_draw),
    draw_index = seq_len(n_draw),
    new_values = 1
  )

  expected_mean <- phi * (1 + 3) / (1 + phi^2)
  expected_var <- sigma_sq * (1 - phi^2) / (1 + phi^2)

  expect_equal(mean(draw[, 1]), expected_mean, tolerance = 0.06)
  expect_equal(stats::var(draw[, 1]), expected_var, tolerance = 0.08)
})

test_that("dense GP new-node prediction matches the conditional moment", {
  n_draw <- 8000L
  sigma_sq <- 2
  phi <- 1.1
  support <- matrix(c(0, 1), ncol = 1, dimnames = list(NULL, "s"))
  new_coord <- matrix(0.5, ncol = 1, dimnames = list(NULL, "s"))
  w_fit <- matrix(rep(c(1, 3), each = n_draw), nrow = n_draw)

  object <- list(
    sigma_sq_samples = matrix(sigma_sq, nrow = n_draw, ncol = 1,
                              dimnames = list(NULL, "gp_1")),
    theta_samples = matrix(phi, nrow = n_draw, ncol = 1,
                           dimnames = list(NULL, "gp_1_phi")),
    w_samples_ordered = list(gp_1 = w_fit)
  )
  term <- list(name = "gp_1", theta_names = "phi", cov_model = "exp")
  graph <- list(coords_ord = support)

  C <- stLMM:::gp_covariance_matrix(support, support, "exp", sigma_sq, c(phi = phi))
  c0 <- stLMM:::gp_covariance_matrix(new_coord, support, "exp", sigma_sq, c(phi = phi))
  C_inv_c0 <- solve(C, as.numeric(c0))
  expected_mean <- sum(C_inv_c0 * w_fit[1, ])
  expected_var <- sigma_sq - sum(as.numeric(c0) * C_inv_c0)

  set.seed(4102)
  draw <- stLMM:::simulate_gp_prediction_nodes(
    object = object,
    term = term,
    graph = graph,
    recover_row = seq_len(n_draw),
    draw_index = seq_len(n_draw),
    new_coords = new_coord
  )

  expect_equal(mean(draw[, 1]), expected_mean, tolerance = 0.06)
  expect_equal(stats::var(draw[, 1]), expected_var, tolerance = 0.08)
})

test_that("joint dense GP new-node prediction matches conditional moments", {
  n_draw <- 12000L
  sigma_sq <- 1.5
  phi <- 0.9
  support <- matrix(c(0, 1, 2), ncol = 1, dimnames = list(NULL, "s"))
  new_coord <- matrix(c(0.4, 1.4), ncol = 1, dimnames = list(NULL, "s"))
  w_fit <- matrix(rep(c(0.5, -0.2, 1.1), each = n_draw), nrow = n_draw)

  object <- list(
    sigma_sq_samples = matrix(sigma_sq, nrow = n_draw, ncol = 1,
                              dimnames = list(NULL, "gp_1")),
    theta_samples = matrix(phi, nrow = n_draw, ncol = 1,
                           dimnames = list(NULL, "gp_1_phi")),
    w_samples_ordered = list(gp_1 = w_fit)
  )
  term <- list(name = "gp_1", theta_names = "phi", cov_model = "exp")
  graph <- list(coords_ord = support)

  C_oo <- stLMM:::gp_covariance_matrix(support, support, "exp", sigma_sq, c(phi = phi))
  C_no <- stLMM:::gp_covariance_matrix(new_coord, support, "exp", sigma_sq, c(phi = phi))
  C_nn <- stLMM:::gp_covariance_matrix(new_coord, new_coord, "exp", sigma_sq, c(phi = phi))
  expected_mean <- as.numeric(C_no %*% solve(C_oo, w_fit[1, ]))
  expected_var <- C_nn - C_no %*% solve(C_oo, t(C_no))

  set.seed(4103)
  draw <- stLMM:::simulate_gp_prediction_nodes(
    object = object,
    term = term,
    graph = graph,
    recover_row = seq_len(n_draw),
    draw_index = seq_len(n_draw),
    new_coords = new_coord,
    joint = TRUE
  )

  expect_equal(colMeans(draw), expected_mean, tolerance = 0.04)
  expect_equal(stats::cov(draw), expected_var, tolerance = 0.05)
})

test_that("joint and independent dense GP predictions have highly correlated marginal summaries", {
  n_draw <- 12000L
  sigma_sq <- 1.5
  phi <- 0.9
  support <- matrix(c(0, 1, 2), ncol = 1, dimnames = list(NULL, "s"))
  new_coord <- matrix(c(0.4, 1.4), ncol = 1, dimnames = list(NULL, "s"))
  w_fit <- matrix(rep(c(0.5, -0.2, 1.1), each = n_draw), nrow = n_draw)

  object <- list(
    sigma_sq_samples = matrix(sigma_sq, nrow = n_draw, ncol = 1,
                              dimnames = list(NULL, "gp_1")),
    theta_samples = matrix(phi, nrow = n_draw, ncol = 1,
                           dimnames = list(NULL, "gp_1_phi")),
    w_samples_ordered = list(gp_1 = w_fit)
  )
  term <- list(name = "gp_1", theta_names = "phi", cov_model = "exp")
  graph <- list(coords_ord = support)

  set.seed(4104)
  draw_ind <- stLMM:::simulate_gp_prediction_nodes(
    object = object,
    term = term,
    graph = graph,
    recover_row = seq_len(n_draw),
    draw_index = seq_len(n_draw),
    new_coords = new_coord,
    joint = FALSE
  )

  set.seed(4105)
  draw_joint <- stLMM:::simulate_gp_prediction_nodes(
    object = object,
    term = term,
    graph = graph,
    recover_row = seq_len(n_draw),
    draw_index = seq_len(n_draw),
    new_coords = new_coord,
    joint = TRUE
  )

  probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  summaries_ind <- c(colMeans(draw_ind), apply(draw_ind, 2L, stats::quantile, probs = probs))
  summaries_joint <- c(colMeans(draw_joint), apply(draw_joint, 2L, stats::quantile, probs = probs))

  expect_gt(stats::cor(summaries_ind, summaries_joint), 0.995)
})

test_that("Gneiting covariance has the expected Datta/Gneiting special-case form", {
  coords <- matrix(
    c(0.0, 0.0,
      0.3, 0.2,
      0.8, 0.6),
    ncol = 2,
    byrow = TRUE
  )
  time <- c(1, 2.5, 4)
  coords_time <- cbind(coords, time = time)
  theta_gneiting <- c(a = 0.4, c = 1.3, alpha = 1, beta = 0, gamma = 0.5, delta = 0)

  C_gneiting <- stLMM:::gp_covariance_matrix(
    coords_time,
    coords_time,
    "gneiting",
    sigma_sq = 2,
    theta = theta_gneiting
  )
  h <- as.matrix(stats::dist(coords))
  u <- abs(outer(time, time, "-"))
  spatial_dim <- ncol(coords)
  t <- 1 + theta_gneiting["a"] * u^(2 * theta_gneiting["alpha"])
  C_expected <- 2 * t^(-(theta_gneiting["delta"] + spatial_dim / 2)) *
    exp(-theta_gneiting["c"] * h^(2 * theta_gneiting["gamma"]) /
          t^(theta_gneiting["beta"] * theta_gneiting["gamma"]))
  dimnames(C_expected) <- NULL

  expect_equal(C_gneiting, C_expected, tolerance = 1e-12)

  C_separable <- 2 * exp(-theta_gneiting["c"] * h) *
    (1 + theta_gneiting["a"] * u^2)^(-1)
  dimnames(C_separable) <- NULL
  expect_equal(C_gneiting, C_separable, tolerance = 1e-12)
})

test_that("compiled NNGP new-node prediction matches the neighbor conditional moment", {
  n_draw <- 8000L
  sigma_sq <- 2
  phi <- 1.1
  support <- matrix(c(0, 1), ncol = 1, dimnames = list(NULL, "s"))
  new_coord <- matrix(0.5, ncol = 1, dimnames = list(NULL, "s"))
  w_fit <- matrix(rep(c(1, 3), each = n_draw), nrow = n_draw)

  object <- list(
    sigma_sq_samples = matrix(sigma_sq, nrow = n_draw, ncol = 1,
                              dimnames = list(NULL, "nngp_1")),
    theta_samples = matrix(phi, nrow = n_draw, ncol = 1,
                           dimnames = list(NULL, "nngp_1_phi")),
    w_samples_ordered = list(nngp_1 = w_fit),
    backend = list(n_omp_threads = 2L)
  )
  term <- list(name = "nngp_1", theta_names = "phi", cov_model = "exp")
  graph <- list(coords_ord = support)

  C <- stLMM:::gp_covariance_matrix(support, support, "exp", sigma_sq, c(phi = phi))
  c0 <- stLMM:::gp_covariance_matrix(new_coord, support, "exp", sigma_sq, c(phi = phi))
  C_inv_c0 <- solve(C, as.numeric(c0))
  expected_mean <- sum(C_inv_c0 * w_fit[1, ])
  expected_var <- sigma_sq - sum(as.numeric(c0) * C_inv_c0)

  set.seed(4103)
  draw <- stLMM:::simulate_nngp_prediction_nodes(
    object = object,
    term = term,
    graph = graph,
    recover_row = seq_len(n_draw),
    draw_index = seq_len(n_draw),
    new_coords = new_coord,
    neighbor_index = matrix(c(1L, 2L), nrow = 1)
  )

  expect_equal(mean(draw[, 1]), expected_mean, tolerance = 0.06)
  expect_equal(stats::var(draw[, 1]), expected_var, tolerance = 0.08)
})
