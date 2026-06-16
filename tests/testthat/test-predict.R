test_that("predict.stLMM works for fixed-effect fits without recovery", {
  set.seed(2001)
  dat <- data.frame(y = rnorm(8), x = seq(-1, 1, length.out = 8))

  fit <- stLMM(
    y ~ x,
    data = dat,
    starting = list(tau_sq = 1),
    priors = list(tau_sq = ig(2, 1)),
    n_samples = 6,
    verbose = FALSE
  )
  pred <- predict(fit, sub_sample = list(start = 2, thin = 2))
  newdata <- data.frame(x = c(-0.5, 0.5))
  pred_new <- predict(fit, newdata = newdata, sub_sample = list(start = 2, thin = 2), y_samples = TRUE)

  expect_s3_class(pred, "stLMM_prediction")
  expect_equal(pred$draw_index, c(2L, 4L, 6L))
  expect_equal(unname(pred$mu_samples), unname(fit$beta_samples[c(2, 4, 6), , drop = FALSE] %*% t(fit$backend$X)))
  expect_equal(dim(pred_new$mu_samples), c(3L, 2L))
  expect_equal(dim(pred_new$y_samples), c(3L, 2L))
})

test_that("predict.stLMM works for grouped random-effect fits without recovery", {
  set.seed(2002)
  dat <- data.frame(
    y = rnorm(12),
    x = rep(c(-1, 0, 1), 4),
    group = factor(rep(letters[1:4], each = 3))
  )

  fit <- stLMM(
    y ~ x + iid(group) + x:iid(group),
    data = dat,
     starting = list(tau_sq = 1, iid_1 = list(sigma_sq = 0.5), iid_2 = list(sigma_sq = 0.5)),
    priors = list(tau_sq = ig(2, 1), iid_1 = list(sigma_sq = ig(2, 1)), iid_2 = list(sigma_sq = ig(2, 1))),
    n_samples = 6,
    verbose = FALSE
  )
  newdata <- data.frame(x = c(2, -0.5), group = factor(c("b", "d"), levels = levels(dat$group)))
  pred <- predict(fit, newdata = newdata, sub_sample = list(start = 2, thin = 2))
  pred_backend <- stLMM:::build_existing_support_prediction_backend(fit, newdata)
  draw_index <- c(2L, 4L, 6L)
  manual <- fit$beta_samples[draw_index, , drop = FALSE] %*% t(pred_backend$X)
  manual <- manual + fit$alpha_samples[draw_index, , drop = FALSE] %*% Matrix::t(pred_backend$Z)

  expect_equal(unname(pred$mu_samples), unname(as.matrix(manual)))

  bad <- newdata
  bad$group <- as.character(bad$group)
  bad$group[1] <- "new_group"
  expect_error(predict(fit, newdata = bad), "new grouping levels")

  slope_only_fit <- stLMM(
    y ~ 1 + x:iid(group),
    data = dat,
     starting = list(tau_sq = 1, iid_1 = list(sigma_sq = 0.5)),
    priors = list(tau_sq = ig(2, 1), iid_1 = list(sigma_sq = ig(2, 1))),
    n_samples = 6,
    verbose = FALSE
  )
  missing_slope <- data.frame(x = c(NA_real_, 0.5), group = factor(c("b", "d"), levels = levels(dat$group)))
  missing_slope$x[1] <- NA_real_
  expect_error(predict(slope_only_fit, newdata = missing_slope), "iid slope covariate x in newdata contains missing values")
})

test_that("predict.stLMM requires recover for process terms", {
  set.seed(2003)
  dat <- data.frame(y = rep(0, 4), time = seq_len(4))

  fit <- stLMM(
    y ~ 0 + ar1(time),
    data = dat,
    starting = list(tau_sq = 1, ar1_1 = c(sigma_sq = 1, phi = 0.4)),
    tuning = list(tau_sq = 0, ar1_1 = c(sigma_sq = 0, phi = 0)),
    priors = list(
      tau_sq = ig(2, 1),
      ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.9, 0.9))
    ),
    n_samples = 4,
    verbose = FALSE
  )

  expect_error(predict(fit), "requires saved or recovered latent process samples")
})

test_that("predict.stLMM_recovery with newdata NULL matches fitted samples", {
  set.seed(201)
  dat <- data.frame(
    y = rep(0, 6),
    time = rep(seq_len(3), each = 2)
  )

  fit <- stLMM(
    y ~ 0 + ar1(time),
    data = dat,
    starting = list(tau_sq = 1, ar1_1 = c(sigma_sq = 1, phi = 0.4)),
    tuning = list(tau_sq = 0, ar1_1 = c(sigma_sq = 0, phi = 0)),
    priors = list(
      tau_sq = ig(2, 1),
      ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.9, 0.9))
    ),
    n_samples = 10,
    verbose = FALSE
  )
  rec <- recover(fit, sub_sample = list(start = 2, thin = 2))

  pred <- predict(rec)
  mu <- fitted(rec, summary = FALSE)
  attr(mu, "draw_index") <- NULL

  expect_s3_class(pred, "stLMM_prediction")
  expect_equal(pred$draw_index, rec$recover_iter)
  expect_equal(unname(pred$mu_samples), unname(mu))
})

test_that("predict.stLMM_recovery maps newdata to existing random effects and process nodes", {
  set.seed(202)
  coords <- data.frame(
    lon = c(0, 1, 2, 0, 1, 2),
    lat = c(0, 0, 0, 1, 1, 1)
  )
  idx_node <- rep(seq_len(nrow(coords)), each = 2)
  dat <- data.frame(
    y = rnorm(length(idx_node)),
    x = rnorm(length(idx_node)),
    x_svc = runif(length(idx_node), 0.5, 1.5),
    group = factor(rep(seq_len(4), length.out = length(idx_node))),
    time = rep(seq_len(4), length.out = length(idx_node)),
    lon = coords$lon[idx_node],
    lat = coords$lat[idx_node]
  )

  fit <- suppressWarnings(stLMM(
    y ~ x + iid(group) + x:iid(group) + ar1(time) + x_svc:nngp(lon, lat, m = 3, ordering = "maxmin"),
    data = dat,
    starting = list(
      tau_sq = 1,
      iid_1 = list(sigma_sq = 0.5),
      iid_2 = list(sigma_sq = 0.5),
      ar1_1 = c(sigma_sq = 0.8, phi = 0.3),
      nngp_1 = c(sigma_sq = 0.7, phi = 0.6)
    ),
    tuning = list(
      tau_sq = 0,
      ar1_1 = c(sigma_sq = 0, phi = 0),
      nngp_1 = c(sigma_sq = 0, phi = 0)
    ),
    priors = list(
      tau_sq = ig(2, 1),
      iid_1 = list(sigma_sq = ig(2, 1)),
      iid_2 = list(sigma_sq = ig(2, 1)),
      ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.9, 0.9)),
      nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5))
    ),
    n_samples = 12,
    verbose = FALSE
  ))
  rec <- recover(fit, sub_sample = list(start = 3, thin = 3))

  row_idx <- c(9L, 2L, 7L, 4L)
  newdata <- dat[row_idx, setdiff(names(dat), "y"), drop = FALSE]
  pred <- predict(rec, newdata = newdata)
  full_mu <- fitted(rec, summary = FALSE)
  draw_index <- attr(full_mu, "draw_index")
  attr(full_mu, "draw_index") <- NULL

  expect_equal(pred$draw_index, draw_index)
  expect_equal(unname(pred$mu_samples), unname(full_mu[, row_idx, drop = FALSE]))

  new_group <- newdata
  new_group$group <- as.character(new_group$group)
  new_group$group[1] <- "not_fit"
  expect_error(predict(rec, newdata = new_group), "new grouping levels")

  new_coord <- newdata
  new_coord$lon[1] <- 999
  pred_new <- predict(rec, newdata = new_coord)
  expect_equal(dim(pred_new$mu_samples), c(length(rec$recover_iter), nrow(new_coord)))
  expect_equal(pred_new$draw_index, rec$recover_iter)
  expect_false(isTRUE(all.equal(pred_new$mu_samples[, 1], pred$mu_samples[, 1])))
})

test_that("predict.stLMM_recovery combines AR1 and NNGP new nodes term-specifically", {
  set.seed(2041)
  coords <- data.frame(
    lon = c(0, 1, 2, 0, 1, 2, 0.5, 1.5),
    lat = c(0, 0, 0, 1, 1, 1, 2, 2)
  )
  dat <- data.frame(
    y = rep(0, nrow(coords)),
    day = c(1, 2, 3, 4, 5, 1, 2, 3),
    lon = coords$lon,
    lat = coords$lat
  )

  fit <- suppressWarnings(stLMM(
    y ~ 0 + ar1(day) + nngp(lon, lat, m = 3, ordering = "maxmin"),
    data = dat,
    starting = list(
      tau_sq = 1,
      ar1_1 = c(sigma_sq = 1, phi = 0.4),
      nngp_1 = c(sigma_sq = 1, phi = 0.7)
    ),
    tuning = list(
      tau_sq = 0,
      ar1_1 = c(sigma_sq = 0, phi = 0),
      nngp_1 = c(sigma_sq = 0, phi = 0)
    ),
    priors = list(
      tau_sq = ig(2, 1),
      ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.9, 0.9)),
      nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5))
    ),
    n_samples = 8,
    verbose = FALSE
  ))
  rec <- recover(fit, sub_sample = list(start = 2, thin = 2))
  newdata <- data.frame(
    day = c(2.5, 2.5, 2.5, 6, 3, 2.5),
    lon = c(0.5, 0.5, 10, 0, 1, 0.5),
    lat = c(0.5, 0.5, 10, 0, 0, 0.5)
  )

  set.seed(20411)
  pred <- predict(rec, newdata = newdata)

  pred_backend <- stLMM:::build_existing_support_prediction_backend(rec, newdata)
  recover_row <- seq_along(rec$recover_iter)
  draw_index <- rec$recover_iter
  manual <- matrix(0, nrow = length(draw_index), ncol = nrow(newdata))

  ar1_map <- pred_backend$process_maps[[1]]
  ar1_mu <- matrix(0, nrow = length(draw_index), ncol = nrow(newdata))
  ar1_existing <- !is.na(ar1_map$map)
  ar1_mu[, ar1_existing] <- rec$w_samples$ar1_1[, ar1_map$map[ar1_existing], drop = FALSE]

  nngp_map <- pred_backend$process_maps[[2]]
  nngp_mu <- matrix(0, nrow = length(draw_index), ncol = nrow(newdata))
  nngp_existing <- !is.na(nngp_map$map)
  nngp_mu[, nngp_existing] <- rec$w_samples_ordered$nngp_1[, nngp_map$map[nngp_existing], drop = FALSE]

  set.seed(20411)
  ar1_new <- stLMM:::simulate_ar1_prediction_nodes(
    object = rec,
    term = rec$backend$process_terms[[1]],
    graph = rec$backend$graphs[[rec$backend$process_terms[[1]]$graph_index]],
    w_samples_ordered = rec$w_samples_ordered,
    recover_row = recover_row,
    draw_index = draw_index,
    new_values = ar1_map$new_values
  )
  nngp_new <- stLMM:::simulate_nngp_prediction_nodes(
    object = rec,
    term = rec$backend$process_terms[[2]],
    graph = rec$backend$graphs[[rec$backend$process_terms[[2]]$graph_index]],
    w_samples_ordered = rec$w_samples_ordered,
    recover_row = recover_row,
    draw_index = draw_index,
    new_coords = nngp_map$new_coords,
    neighbor_index = nngp_map$neighbor_index
  )
  ar1_mu[, !ar1_existing] <- ar1_new[, ar1_map$new_id[!ar1_existing], drop = FALSE]
  nngp_mu[, !nngp_existing] <- nngp_new[, nngp_map$new_id[!nngp_existing], drop = FALSE]
  manual <- manual + ar1_mu + nngp_mu

  expect_equal(dim(pred$mu_samples), c(length(rec$recover_iter), nrow(newdata)))
  expect_equal(pred$draw_index, rec$recover_iter)
  expect_named(pred$w_samples, c("ar1_1", "nngp_1"))
  expect_equal(dim(pred$w_samples$ar1_1), dim(pred$mu_samples))
  expect_equal(dim(pred$w_samples$nngp_1), dim(pred$mu_samples))
  expect_equal(unname(pred$w_samples$ar1_1), unname(ar1_mu))
  expect_equal(unname(pred$w_samples$nngp_1), unname(nngp_mu))
  expect_equal(pred$mu_samples[, 1], pred$mu_samples[, 2])
  expect_equal(pred$mu_samples[, 1], pred$mu_samples[, 6])
  expect_equal(unname(pred$mu_samples), unname(manual))

  pred_no_w <- predict(rec, newdata = newdata, return_w_samples = FALSE)
  expect_null(pred_no_w$w_samples)
})

test_that("predict.stLMM_recovery combines grouped random effects with SVC NNGP new nodes", {
  set.seed(2042)
  coords <- data.frame(
    lon = c(0, 1, 2, 0, 1, 2, 0.5, 1.5, 0, 1, 2, 1),
    lat = c(0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 3, 1.5)
  )
  dat <- data.frame(
    y = rnorm(nrow(coords)),
    x = rep(c(-1, 0, 1), 4),
    x_svc = seq(0.8, 1.5, length.out = nrow(coords)),
    group = factor(rep(letters[1:4], each = 3)),
    lon = coords$lon,
    lat = coords$lat
  )

  fit <- suppressWarnings(stLMM(
    y ~ x + iid(group) + x:iid(group) + x_svc:nngp(lon, lat, m = 3, ordering = "maxmin"),
    data = dat,
    starting = list(
      tau_sq = 1,
      iid_1 = list(sigma_sq = 0.5),
      iid_2 = list(sigma_sq = 0.5),
      nngp_1 = c(sigma_sq = 1, phi = 0.7)
    ),
    tuning = list(
      tau_sq = 0,
      nngp_1 = c(sigma_sq = 0, phi = 0)
    ),
    priors = list(
      tau_sq = ig(2, 1),
      iid_1 = list(sigma_sq = ig(2, 1)),
      iid_2 = list(sigma_sq = ig(2, 1)),
      nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5))
    ),
    n_samples = 8,
    verbose = FALSE
  ))
  rec <- recover(fit, sub_sample = list(start = 2, thin = 2))
  newdata <- data.frame(
    x = c(0.2, -0.7, 0.2, 0.5, 0.2),
    x_svc = c(1, 1, 2, 1.5, 1),
    group = factor(c("b", "d", "b", "c", "b"), levels = levels(dat$group)),
    lon = c(0.5, 0.5, 0.5, 0, 0.5),
    lat = c(0.5, 0.5, 0.5, 0, 0.5)
  )

  set.seed(20421)
  pred <- predict(rec, newdata = newdata)

  pred_backend <- stLMM:::build_existing_support_prediction_backend(rec, newdata)
  recover_row <- seq_along(rec$recover_iter)
  draw_index <- rec$recover_iter
  manual <- rec$beta_samples[draw_index, , drop = FALSE] %*% t(pred_backend$X)
  manual <- manual + rec$alpha_samples[draw_index, , drop = FALSE] %*% Matrix::t(pred_backend$Z)
  manual <- as.matrix(manual)

  nngp_map <- pred_backend$process_maps[[1]]
  term_mu <- matrix(0, nrow = length(draw_index), ncol = nrow(newdata))
  existing <- !is.na(nngp_map$map)
  term_mu[, existing] <- rec$w_samples_ordered$nngp_1[, nngp_map$map[existing], drop = FALSE]

  set.seed(20421)
  nngp_new <- stLMM:::simulate_nngp_prediction_nodes(
    object = rec,
    term = rec$backend$process_terms[[1]],
    graph = rec$backend$graphs[[rec$backend$process_terms[[1]]$graph_index]],
    w_samples_ordered = rec$w_samples_ordered,
    recover_row = recover_row,
    draw_index = draw_index,
    new_coords = nngp_map$new_coords,
    neighbor_index = nngp_map$neighbor_index
  )
  term_mu[, !existing] <- nngp_new[, nngp_map$new_id[!existing], drop = FALSE]
  term_w <- term_mu
  term_mu <- sweep(term_w, 2L, nngp_map$scale, `*`)
  manual <- manual + term_mu

  expect_equal(dim(pred$mu_samples), c(length(rec$recover_iter), nrow(newdata)))
  expect_named(pred$w_samples, "nngp_1")
  expect_equal(unname(pred$w_samples$nngp_1), unname(term_w))
  expect_equal(pred$mu_samples[, 1], pred$mu_samples[, 5])
  expect_equal(pred$w_samples$nngp_1[, 1], pred$w_samples$nngp_1[, 3])
  expect_equal(term_mu[, 3], 2 * term_mu[, 1])
  expect_equal(unname(pred$mu_samples), unname(manual))
})

test_that("predict.stLMM_recovery can simulate y samples", {
  set.seed(203)
  dat <- data.frame(y = rep(0, 4), time = seq_len(4))

  fit <- stLMM(
    y ~ 0 + ar1(time),
    data = dat,
    starting = list(tau_sq = 1, ar1_1 = c(sigma_sq = 1, phi = 0.4)),
    tuning = list(tau_sq = 0, ar1_1 = c(sigma_sq = 0, phi = 0)),
    priors = list(
      tau_sq = ig(2, 1),
      ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.9, 0.9))
    ),
    n_samples = 6,
    verbose = FALSE
  )
  rec <- recover(fit, sub_sample = list(start = 2, thin = 2))
  pred <- predict(rec, y_samples = TRUE)

  expect_equal(dim(pred$y_samples), dim(pred$mu_samples))
  expect_equal(pred$draw_index, c(2L, 4L, 6L))
})

test_that("predict.stLMM_recovery simulates new AR1 nodes with repeated-row sharing", {
  set.seed(205)
  dat <- data.frame(y = rep(0, 5), time = seq_len(5))

  fit <- stLMM(
    y ~ 0 + ar1(time),
    data = dat,
    starting = list(tau_sq = 1, ar1_1 = c(sigma_sq = 1, phi = 0.4)),
    tuning = list(tau_sq = 0, ar1_1 = c(sigma_sq = 0, phi = 0)),
    priors = list(
      tau_sq = ig(2, 1),
      ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.9, 0.9))
    ),
    n_samples = 20,
    verbose = FALSE
  )
  rec <- recover(fit, sub_sample = list(start = 2, thin = 2))
  pred <- predict(rec, newdata = data.frame(time = c(2.5, 2.5, 6)))

  expect_equal(dim(pred$mu_samples), c(10L, 3L))
  expect_equal(pred$mu_samples[, 1], pred$mu_samples[, 2])
  expect_false(isTRUE(all.equal(pred$mu_samples[, 1], pred$mu_samples[, 3])))
})

test_that("AR1 new-node prediction uses ordered-support adjacency", {
  set.seed(2051)
  n_draw <- 40000L
  sigma_sq <- 2
  phi <- 0.5
  object <- list(
    sigma_sq_samples = matrix(sigma_sq, n_draw, 1, dimnames = list(NULL, "ar1_1")),
    theta_samples = matrix(phi, n_draw, 1, dimnames = list(NULL, "ar1_1_phi")),
    w_samples = list(ar1_1 = matrix(c(2, 4), n_draw, 2, byrow = TRUE))
  )
  term <- list(name = "ar1_1")
  graph <- list(support = c(1, 10))

  interior <- stLMM:::simulate_ar1_prediction_nodes(
    object, term, graph, seq_len(n_draw), seq_len(n_draw), 5
  )[, 1]
  endpoint <- stLMM:::simulate_ar1_prediction_nodes(
    object, term, graph, seq_len(n_draw), seq_len(n_draw), 100
  )[, 1]

  expect_equal(mean(interior), phi * (2 + 4) / (1 + phi^2), tolerance = 0.02)
  expect_equal(stats::var(interior), sigma_sq * (1 - phi^2) / (1 + phi^2), tolerance = 0.03)
  expect_equal(mean(endpoint), phi * 4, tolerance = 0.02)
  expect_equal(stats::var(endpoint), sigma_sq * (1 - phi^2), tolerance = 0.03)
})

test_that("AR1 new-node prediction handles one fitted support point", {
  set.seed(2052)
  n_draw <- 40000L
  sigma_sq <- 2
  phi <- -0.4
  object <- list(
    sigma_sq_samples = matrix(sigma_sq, n_draw, 1, dimnames = list(NULL, "ar1_1")),
    theta_samples = matrix(phi, n_draw, 1, dimnames = list(NULL, "ar1_1_phi")),
    w_samples = list(ar1_1 = matrix(3, n_draw, 1))
  )
  term <- list(name = "ar1_1")
  graph <- list(support = 10)

  left <- stLMM:::simulate_ar1_prediction_nodes(
    object, term, graph, seq_len(n_draw), seq_len(n_draw), 1
  )[, 1]
  right <- stLMM:::simulate_ar1_prediction_nodes(
    object, term, graph, seq_len(n_draw), seq_len(n_draw), 100
  )[, 1]

  expect_equal(mean(left), phi * 3, tolerance = 0.02)
  expect_equal(stats::var(left), sigma_sq * (1 - phi^2), tolerance = 0.03)
  expect_equal(mean(right), phi * 3, tolerance = 0.02)
  expect_equal(stats::var(right), sigma_sq * (1 - phi^2), tolerance = 0.03)
})

test_that("predict.stLMM_recovery applies SVC scaling to new AR1 nodes", {
  set.seed(2053)
  dat <- data.frame(y = rep(0, 5), x = rep(1, 5), time = seq_len(5))

  fit <- stLMM(
    y ~ 0 + x:ar1(time),
    data = dat,
    starting = list(tau_sq = 1, ar1_1 = c(sigma_sq = 1, phi = 0.4)),
    tuning = list(tau_sq = 0, ar1_1 = c(sigma_sq = 0, phi = 0)),
    priors = list(
      tau_sq = ig(2, 1),
      ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.9, 0.9))
    ),
    n_samples = 20,
    verbose = FALSE
  )
  rec <- recover(fit, sub_sample = list(start = 2, thin = 2))
  pred <- predict(rec, newdata = data.frame(x = c(1, 2), time = c(2.5, 2.5)))

  expect_equal(pred$mu_samples[, 2], 2 * pred$mu_samples[, 1])
})

test_that("predict.stLMM_recovery rebuilds grouped random-effect design for newdata", {
  set.seed(2054)
  dat <- data.frame(
    y = rnorm(12),
    x = rep(c(-1, 0, 1), 4),
    group = factor(rep(letters[1:4], each = 3)),
    time = rep(seq_len(3), 4)
  )

  fit <- suppressWarnings(stLMM(
    y ~ 0 + iid(group) + x:iid(group) + ar1(time),
    data = dat,
     starting = list(tau_sq = 1, iid_1 = list(sigma_sq = 0.5), iid_2 = list(sigma_sq = 0.5), ar1_1 = c(sigma_sq = 1, phi = 0.4)),
    tuning = list(tau_sq = 0, ar1_1 = c(sigma_sq = 0, phi = 0)),
    priors = list(
      tau_sq = ig(2, 1),
      iid_1 = list(sigma_sq = ig(2, 1)),
      iid_2 = list(sigma_sq = ig(2, 1)),
      ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.9, 0.9))
    ),
    n_samples = 8,
    verbose = FALSE
  ))
  rec <- recover(fit, sub_sample = list(start = 2, thin = 2))

  newdata <- data.frame(group = factor(c("b", "d"), levels = levels(dat$group)),
                        x = c(2, -0.5),
                        time = c(1, 3))
  pred <- predict(rec, newdata = newdata)
  pred_backend <- stLMM:::build_existing_support_prediction_backend(rec, newdata)
  draw_index <- rec$recover_iter
  manual <- rec$beta_samples[draw_index, , drop = FALSE] %*% t(pred_backend$X)
  manual <- manual + rec$alpha_samples[draw_index, , drop = FALSE] %*% Matrix::t(pred_backend$Z)
  manual <- as.matrix(manual) +
    rec$w_samples$ar1_1[, pred_backend$process_maps[[1]]$map, drop = FALSE]

  expect_equal(unname(pred$mu_samples), unname(manual))
  bad <- newdata
  bad$group <- as.character(bad$group)
  bad$group[1] <- "new_group"
  expect_error(predict(rec, newdata = bad), "new grouping levels")
})

test_that("predict.stLMM_recovery simulates new dense GP nodes with repeated-row sharing", {
  set.seed(206)
  dat <- data.frame(
    y = rep(0, 6),
    lon = c(0, 1, 2, 0, 1, 2),
    lat = c(0, 0, 0, 1, 1, 1)
  )

  fit <- stLMM(
    y ~ 0 + gp(lon, lat),
    data = dat,
    starting = list(tau_sq = 1, gp_1 = c(sigma_sq = 1, phi = 0.7)),
    tuning = list(tau_sq = 0, gp_1 = c(sigma_sq = 0, phi = 0)),
    priors = list(
      tau_sq = ig(2, 1),
      gp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5))
    ),
    n_samples = 20,
    verbose = FALSE
  )
  rec <- recover(fit, sub_sample = list(start = 2, thin = 2))
  pred <- predict(rec, newdata = data.frame(
    lon = c(0.5, 0.5, 10),
    lat = c(0.5, 0.5, 10)
  ))

  expect_equal(dim(pred$mu_samples), c(10L, 3L))
  expect_equal(pred$mu_samples[, 1], pred$mu_samples[, 2])
  expect_gt(stats::sd(pred$mu_samples[, 3] - pred$mu_samples[, 1]), 0)

  pred_joint <- predict(rec, newdata = data.frame(
    lon = c(0.5, 10),
    lat = c(0.5, 10)
  ), joint = TRUE)

  expect_true(pred_joint$joint)
  expect_equal(pred_joint$joint_method, "full")
  expect_equal(dim(pred_joint$w_samples$gp_1), c(10L, 2L))
})

test_that("predict.stLMM_recovery simulates new NNGP nodes with repeated-row sharing", {
  set.seed(207)
  dat <- data.frame(
    y = rep(0, 6),
    lon = c(0, 1, 2, 0, 1, 2),
    lat = c(0, 0, 0, 1, 1, 1)
  )

  fit <- suppressWarnings(stLMM(
    y ~ 0 + nngp(lon, lat, m = 3),
    data = dat,
    starting = list(tau_sq = 1, nngp_1 = c(sigma_sq = 1, phi = 0.7)),
    tuning = list(tau_sq = 0, nngp_1 = c(sigma_sq = 0, phi = 0)),
    priors = list(
      tau_sq = ig(2, 1),
      nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5))
    ),
    n_samples = 6,
    verbose = FALSE
  ))
  rec <- recover(fit, sub_sample = list(start = 2, thin = 2))
  pred <- predict(rec, newdata = data.frame(
    lon = c(0.5, 0.5, 10),
    lat = c(0.5, 0.5, 10)
  ))

  expect_equal(dim(pred$mu_samples), c(3L, 3L))
  expect_named(pred$w_samples, "nngp_1")
  expect_equal(dim(pred$w_samples$nngp_1), dim(pred$mu_samples))
  expect_equal(pred$w_samples$nngp_1[, 1], pred$w_samples$nngp_1[, 2])
  expect_equal(pred$mu_samples[, 1], pred$mu_samples[, 2])
  expect_gt(stats::sd(pred$mu_samples[, 3] - pred$mu_samples[, 1]), 0)
  expect_match(paste(capture.output(print(pred)), collapse = "\n"), "process samples: nngp_1", fixed = TRUE)
  expect_match(paste(capture.output(summary(pred)), collapse = "\n"), "process samples: nngp_1", fixed = TRUE)

  pred_no_w <- predict(rec, newdata = data.frame(
    lon = c(0.5, 0.5, 10),
    lat = c(0.5, 0.5, 10)
  ), return_w_samples = FALSE)
  expect_null(pred_no_w$w_samples)
  expect_match(paste(capture.output(print(pred_no_w)), collapse = "\n"), "process samples: not retained", fixed = TRUE)
})

test_that("compiled NNGP prediction uses R RNG reproducibly across thread counts", {
  set.seed(2071)
  dat <- data.frame(
    y = rep(0, 8),
    lon = c(0, 1, 2, 3, 0, 1, 2, 3),
    lat = c(0, 0, 0, 0, 1, 1, 1, 1)
  )

  fit <- suppressWarnings(stLMM(
    y ~ 0 + nngp(lon, lat, m = 3),
    data = dat,
    starting = list(tau_sq = 1, nngp_1 = c(sigma_sq = 1, phi = 0.7)),
    tuning = list(tau_sq = 0, nngp_1 = c(sigma_sq = 0, phi = 0)),
    priors = list(
      tau_sq = ig(2, 1),
      nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5))
    ),
    n_samples = 6,
    n_omp_threads = 1,
    verbose = FALSE
  ))
  rec_1 <- recover(fit, sub_sample = list(start = 2, thin = 2))
  rec_2 <- rec_1
  rec_2$backend$n_omp_threads <- 2L
  newdata <- data.frame(lon = c(0.5, 1.5, 2.5), lat = c(0.5, 0.5, 0.5))

  set.seed(99)
  pred_1 <- predict(rec_1, newdata = newdata)
  set.seed(99)
  pred_2 <- predict(rec_2, newdata = newdata)

  expect_equal(pred_1$mu_samples, pred_2$mu_samples)
})

test_that("NNGP joint prediction agrees with non-joint margins and adds cross-node covariance", {
  set.seed(2072)
  dat <- data.frame(
    y = rep(0, 8),
    lon = c(0, 1, 2, 3, 0, 1, 2, 3),
    lat = c(0, 0, 0, 0, 1, 1, 1, 1)
  )

  fit <- suppressWarnings(stLMM(
    y ~ 0 + nngp(lon, lat, m = 4),
    data = dat,
    starting = list(tau_sq = 1, nngp_1 = c(sigma_sq = 1, phi = 0.7)),
    tuning = list(tau_sq = 0, nngp_1 = c(sigma_sq = 0, phi = 0)),
    priors = list(
      tau_sq = ig(2, 1),
      nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5))
    ),
    n_samples = 1200,
    verbose = FALSE
  ))
  rec <- recover(fit, sub_sample = list(start = 101, thin = 2))

  existing <- dat[c(1, 5, 8), , drop = FALSE]
  pred_existing_false <- predict(rec, newdata = existing, joint = FALSE)
  pred_existing_true <- predict(rec, newdata = existing, joint = TRUE)
  pred_existing_vecchia <- predict(rec, newdata = existing, joint = TRUE, joint_method = "vecchia")
  expect_equal(pred_existing_true$mu_samples, pred_existing_false$mu_samples)
  expect_equal(pred_existing_vecchia$mu_samples, pred_existing_false$mu_samples)

  newdata <- data.frame(
    lon = c(0.5, 0.5, 1.5, 2.5),
    lat = c(0.5, 0.5, 0.5, 0.5)
  )

  set.seed(77)
  pred_false <- predict(rec, newdata = newdata, joint = FALSE)
  set.seed(77)
  expect_message(
    pred_true <- predict(rec, newdata = newdata, joint = TRUE),
    "dense covariance"
  )
  set.seed(77)
  pred_vecchia <- predict(rec, newdata = newdata, joint = TRUE, joint_method = "vecchia")

  expect_true(pred_true$joint)
  expect_true(pred_vecchia$joint)
  expect_equal(pred_vecchia$joint_method, "vecchia")
  expect_false(pred_false$joint)
  expect_equal(dim(pred_true$mu_samples), dim(pred_false$mu_samples))
  expect_equal(dim(pred_vecchia$mu_samples), dim(pred_false$mu_samples))
  expect_equal(pred_true$mu_samples[, 1], pred_true$mu_samples[, 2])
  expect_equal(pred_vecchia$mu_samples[, 1], pred_vecchia$mu_samples[, 2])

  false_unique <- pred_false$mu_samples[, c(1, 3, 4), drop = FALSE]
  true_unique <- pred_true$mu_samples[, c(1, 3, 4), drop = FALSE]
  vecchia_unique <- pred_vecchia$mu_samples[, c(1, 3, 4), drop = FALSE]
  expect_lt(max(abs(colMeans(true_unique) - colMeans(false_unique))), 0.18)
  expect_lt(max(abs(apply(true_unique, 2, var) - apply(false_unique, 2, var))), 0.25)
  expect_lt(max(abs(colMeans(vecchia_unique) - colMeans(false_unique))), 0.20)
  expect_lt(max(abs(apply(vecchia_unique, 2, var) - apply(false_unique, 2, var))), 0.30)

  expect_gt(abs(stats::cor(true_unique[, 1], true_unique[, 2])), 0.08)
  expect_gt(abs(stats::cor(vecchia_unique[, 1], vecchia_unique[, 2])), 0.08)
  expect_lt(abs(stats::cor(false_unique[, 1], false_unique[, 2])), 0.25)

  row_order <- c(4, 1, 3, 2)
  set.seed(99)
  pred_ordered <- suppressMessages(predict(rec, newdata = newdata, joint = TRUE))
  set.seed(99)
  pred_reordered <- suppressMessages(predict(rec, newdata = newdata[row_order, ], joint = TRUE))
  expect_equal(pred_ordered$mu_samples[, row_order], pred_reordered$mu_samples)

  set.seed(99)
  pred_vecchia_ordered <- predict(rec, newdata = newdata, joint = TRUE, joint_method = "vecchia")
  set.seed(99)
  pred_vecchia_reordered <- predict(rec, newdata = newdata[row_order, ], joint = TRUE, joint_method = "vecchia")
  expect_equal(pred_vecchia_ordered$mu_samples[, row_order], pred_vecchia_reordered$mu_samples)
})

test_that("Vecchia NNGP prediction graph uses fitted and previous prediction history", {
  support <- matrix(
    c(
      0, 0,
      1, 0,
      0, 1
    ),
    ncol = 2,
    byrow = TRUE
  )
  new_coords <- matrix(
    c(
      0.2, 0.2,
      0.3, 0.2,
      10, 10
    ),
    ncol = 2,
    byrow = TRUE
  )

  graph <- stLMM:::nngp_prediction_vecchia_graph(
    new_coords = new_coords,
    support = support,
    m = 4,
    ordering = "coord",
    cov_model = "exp",
    st_scale = 1,
    term_name = "nngp_1",
    n_omp_threads = 1
  )

  expect_equal(nrow(graph$coords_all), 6L)
  expect_equal(graph$neighbor_count, c(3L, 4L, 4L))
  expect_true(all(graph$neighbor_index[1, seq_len(graph$neighbor_count[1])] <= 3L))
  expect_true(4L %in% graph$neighbor_index[2, seq_len(graph$neighbor_count[2])])
  expect_true(all(graph$neighbor_index[3, seq_len(graph$neighbor_count[3])] < 6L))
  padded <- which(graph$neighbor_count < ncol(graph$neighbor_index))
  expect_true(all(vapply(padded, function(i){
    all(graph$neighbor_index[i, seq.int(graph$neighbor_count[i] + 1L, ncol(graph$neighbor_index))] == 0L)
  }, logical(1))))
})

test_that("Vecchia NNGP simulator reuses previous prediction draws", {
  set.seed(2073)
  n_draw <- 2500L
  support <- matrix(c(0, 10), ncol = 1)
  new_coords <- matrix(c(5, 5.1), ncol = 1)
  coords_all <- rbind(support, new_coords)
  neighbor_index <- matrix(c(1L, 2L, 1L, 3L), nrow = 2, byrow = TRUE)
  neighbor_count <- c(2L, 2L)
  object <- list(
    backend = list(n_omp_threads = 1L),
    sigma_sq_samples = matrix(1, n_draw, 1, dimnames = list(NULL, "nngp_1_sigma_sq")),
    theta_samples = matrix(0.5, n_draw, 1, dimnames = list(NULL, "nngp_1_phi")),
    w_samples_ordered = list(nngp_1 = matrix(0, n_draw, 2))
  )
  term <- list(name = "nngp_1", theta_names = "phi", cov_model = "exp")
  graph <- list(coords_ord = support)

  draws <- stLMM:::simulate_nngp_prediction_nodes_vecchia(
    object = object,
    term = term,
    graph = graph,
    recover_row = seq_len(n_draw),
    draw_index = seq_len(n_draw),
    coords_all = coords_all,
    neighbor_index = neighbor_index,
    neighbor_count = neighbor_count
  )

  expect_equal(dim(draws), c(n_draw, 2L))
  expect_gt(abs(stats::cor(draws[, 1], draws[, 2])), 0.35)
})

test_that("compiled NNGP prediction neighbor search matches scaled R ordering", {
  support <- matrix(
    c(
      0, 0, 0,
      1, 0, 0,
      2, 0, 0,
      0, 1, 1,
      1, 1, 1
    ),
    ncol = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("lon", "lat", "time"))
  )
  new_coords <- matrix(
    c(
      0.4, 0.1, 0.2,
      1.5, 0.4, 0.8
    ),
    ncol = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("lon", "lat", "time"))
  )

  st_scale <- 2
  support_scaled <- support
  new_scaled <- new_coords
  support_scaled[, 3] <- support_scaled[, 3] * st_scale
  new_scaled[, 3] <- new_scaled[, 3] * st_scale

  expected <- matrix(NA_integer_, nrow = nrow(new_coords), ncol = 3)
  for(i in seq_len(nrow(new_coords))){
    d <- rowSums((support_scaled - matrix(new_scaled[i, ], nrow(support_scaled), 3, byrow = TRUE))^2)
    expected[i, ] <- order(d)[seq_len(3)]
  }

  got_1 <- stLMM:::nngp_prediction_neighbors(new_coords, support, 3, "sep_exp", st_scale, "nngp_1", 1)
  got_2 <- stLMM:::nngp_prediction_neighbors(new_coords, support, 3, "sep_exp", st_scale, "nngp_1", 2)

  expect_equal(got_1, expected)
  expect_equal(got_2, expected)

  got_unscaled <- stLMM:::nngp_prediction_neighbors(new_coords, support, 3, "sep_exp", 1, "nngp_1", 1)
  support_time_scaled <- support
  new_time_scaled <- new_coords
  support_time_scaled[, 3] <- support_time_scaled[, 3] * st_scale
  new_time_scaled[, 3] <- new_time_scaled[, 3] * st_scale
  support_space_scaled <- support
  new_space_scaled <- new_coords
  support_space_scaled[, 1:2] <- support_space_scaled[, 1:2] * st_scale
  new_space_scaled[, 1:2] <- new_space_scaled[, 1:2] * st_scale

  expected_time_only <- matrix(NA_integer_, nrow = nrow(new_coords), ncol = 3)
  expected_space_only <- matrix(NA_integer_, nrow = nrow(new_coords), ncol = 3)
  for(i in seq_len(nrow(new_coords))){
    d_time <- rowSums((support_time_scaled - matrix(new_time_scaled[i, ], nrow(support), 3, byrow = TRUE))^2)
    d_space <- rowSums((support_space_scaled - matrix(new_space_scaled[i, ], nrow(support), 3, byrow = TRUE))^2)
    expected_time_only[i, ] <- order(d_time)[seq_len(3)]
    expected_space_only[i, ] <- order(d_space)[seq_len(3)]
  }

  expect_equal(got_1, expected_time_only)
  expect_false(identical(got_1, expected_space_only))
  expect_false(identical(got_1, got_unscaled))
})

test_that("NNGP st_scale affects only space-time covariance neighbor ranking", {
  support <- matrix(
    c(
      0, 0,
      0, 10,
      10, 0
    ),
    ncol = 2,
    byrow = TRUE
  )
  new_coords <- matrix(c(0, 4), nrow = 1)

  exp_unscaled <- stLMM:::nngp_prediction_neighbors(new_coords, support, 2, "exp", 1, "nngp_1", 1)
  exp_scaled <- stLMM:::nngp_prediction_neighbors(new_coords, support, 2, "exp", 100, "nngp_1", 1)

  expect_equal(exp_scaled, exp_unscaled)
})

test_that("predict.stLMM_recovery works for all NNGP ordering options", {
  set.seed(211)
  coords <- data.frame(
    lon = c(0, 1, 2, 0, 1, 2, 0.5, 1.5),
    lat = c(0, 0, 0, 1, 1, 1, 2, 2)
  )
  dat <- data.frame(
    y = rep(0, nrow(coords)),
    lon = coords$lon,
    lat = coords$lat
  )
  newdata <- data.frame(
    lon = c(0.25, 1.25, 10),
    lat = c(0.25, 1.25, 10)
  )
  orderings <- list(
    coord = "coord",
    default = "default",
    maxmin = "maxmin",
    hilbert = "hilbert",
    random = "random",
    user = c(4L, 1L, 7L, 2L, 5L, 8L, 3L, 6L)
  )

  for(ord_name in names(orderings)){
    set.seed(2110 + match(ord_name, names(orderings)))
    form <- if(is.numeric(orderings[[ord_name]])){
      as.formula(paste0(
        "y ~ 0 + nngp(lon, lat, m = 3, ordering = c(",
        paste(orderings[[ord_name]], collapse = ","),
        "))"
      ))
    } else {
      as.formula(paste0(
        "y ~ 0 + nngp(lon, lat, m = 3, ordering = \"",
        orderings[[ord_name]],
        "\")"
      ))
    }

    fit <- suppressWarnings(stLMM(
      form,
      data = dat,
      starting = list(tau_sq = 1, nngp_1 = c(sigma_sq = 1, phi = 0.7)),
      tuning = list(tau_sq = 0, nngp_1 = c(sigma_sq = 0, phi = 0)),
      priors = list(tau_sq = ig(2, 1), nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5))),
      n_samples = 6,
      verbose = FALSE
    ))
    rec <- recover(fit, sub_sample = list(start = 2, thin = 2))
    set.seed(500 + match(ord_name, names(orderings)))
    pred <- predict(rec, newdata = newdata)

    graph <- fit$backend$graphs[[1]]
    expected_ordering <- if(ord_name == "user") "user" else orderings[[ord_name]]
    expect_equal(graph$ordering, expected_ordering)
    expect_equal(dim(pred$mu_samples), c(length(rec$recover_iter), nrow(newdata)))
    expect_equal(attr(rec$w_samples$nngp_1, "node_order"), "support")
    expect_equal(attr(rec$w_samples_ordered$nngp_1, "node_order"), "internal")
    expect_equal(ncol(rec$w_samples$nngp_1), nrow(coords))
    expect_equal(ncol(rec$w_samples_ordered$nngp_1), nrow(coords))
  }
})

test_that("predict.stLMM_recovery applies SVC scaling to new NNGP nodes", {
  set.seed(208)
  dat <- data.frame(
    y = rep(0, 6),
    x_svc = rep(1, 6),
    lon = c(0, 1, 2, 0, 1, 2),
    lat = c(0, 0, 0, 1, 1, 1)
  )

  fit <- suppressWarnings(stLMM(
    y ~ 0 + x_svc:nngp(lon, lat, m = 3),
    data = dat,
    starting = list(tau_sq = 1, nngp_1 = c(sigma_sq = 1, phi = 0.7)),
    tuning = list(tau_sq = 0, nngp_1 = c(sigma_sq = 0, phi = 0)),
    priors = list(
      tau_sq = ig(2, 1),
      nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5))
    ),
    n_samples = 6,
    verbose = FALSE
  ))
  rec <- recover(fit, sub_sample = list(start = 2, thin = 2))
  pred <- predict(rec, newdata = data.frame(
    x_svc = c(1, 2),
    lon = c(0.5, 0.5),
    lat = c(0.5, 0.5)
  ))
  pred_joint <- suppressMessages(predict(rec, newdata = data.frame(
    x_svc = c(1, 2),
    lon = c(0.5, 0.5),
    lat = c(0.5, 0.5)
  ), joint = TRUE))

  expect_equal(pred$mu_samples[, 2], 2 * pred$mu_samples[, 1])
  expect_equal(pred$w_samples$nngp_1[, 2], pred$w_samples$nngp_1[, 1])
  expect_equal(pred_joint$mu_samples[, 2], 2 * pred_joint$mu_samples[, 1])
  expect_equal(pred_joint$w_samples$nngp_1[, 2], pred_joint$w_samples$nngp_1[, 1])
})

test_that("predict.stLMM_recovery maps existing NNGP supports across covariance families", {
  set.seed(204)
  coords <- data.frame(
    lon = c(0, 1, 2, 0, 1, 2, 0.5, 1.5),
    lat = c(0, 0, 0, 1, 1, 1, 2, 2),
    time = c(0, 0, 0, 0.5, 0.5, 0.5, 1, 1)
  )
  idx_node <- rep(seq_len(nrow(coords)), each = 2)
  dat <- data.frame(
    y = rnorm(length(idx_node)),
    lon = coords$lon[idx_node],
    lat = coords$lat[idx_node],
    time = coords$time[idx_node]
  )
  row_idx <- c(12L, 3L, 9L, 2L)

  specs <- list(
    exp = list(
      formula = y ~ 0 + nngp(lon, lat, m = 3, cov_model = "exp", ordering = "maxmin"),
      starting = c(sigma_sq = 1, phi = 0.7),
      priors = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5))
    ),
    matern = list(
      formula = y ~ 0 + nngp(lon, lat, m = 3, cov_model = "matern", ordering = "maxmin"),
      starting = c(sigma_sq = 1, phi = 0.7, nu = 1.1),
      priors = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5), nu = uniform(0.2, 3))
    ),
    sep_exp = list(
      formula = y ~ 0 + nngp(lon, lat, time, m = 3, cov_model = "sep_exp", ordering = "maxmin"),
      starting = c(sigma_sq = 1, phi = 0.7, lambda = 0.4),
      priors = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5), lambda = uniform(0.1, 5))
    ),
    multi_res_sep_exp = list(
      formula = y ~ 0 + nngp(lon, lat, time, m = 3, cov_model = "multi_res_sep_exp", ordering = "maxmin"),
      starting = c(sigma_sq = 1, alpha = 0.5, phi_1 = 0.7, lambda_1 = 0.4, phi_2 = 1.2, lambda_2 = 0.8),
      priors = list(
        sigma_sq = ig(2, 1), alpha = uniform(0.05, 0.95),
        phi_1 = uniform(0.1, 5), lambda_1 = uniform(0.1, 5),
        phi_2 = uniform(0.1, 5), lambda_2 = uniform(0.1, 5)
      )
    ),
    gneiting = list(
      formula = y ~ 0 + nngp(lon, lat, time, m = 3, cov_model = "gneiting", ordering = "maxmin"),
      starting = c(sigma_sq = 1, a = 0.4, c = 0.7, alpha = 1, beta = 0.6, gamma = 0.5, delta = 0),
      priors = list(
        sigma_sq = ig(2, 1), a = uniform(0.1, 5),
        c = uniform(0.1, 5), alpha = uniform(0.05, 1),
        beta = uniform(0.05, 0.95), gamma = uniform(0.05, 1),
        delta = uniform(0, 2)
      )
    )
  )

  for(spec_name in names(specs)){
    spec <- specs[[spec_name]]
    tuning_i <- spec$starting
    tuning_i[] <- 0
    tuning_i["sigma_sq"] <- 0

    fit <- suppressWarnings(stLMM(
      spec$formula,
      data = dat,
      starting = list(tau_sq = 1, nngp_1 = spec$starting),
      tuning = list(tau_sq = 0, nngp_1 = tuning_i),
      priors = list(tau_sq = ig(2, 1), nngp_1 = spec$priors),
      n_samples = 6,
      verbose = FALSE
    ))
    rec <- recover(fit, sub_sample = list(start = 2, thin = 2))
    pred <- predict(rec, newdata = dat[row_idx, setdiff(names(dat), "y"), drop = FALSE])
    full_mu <- fitted(rec, summary = FALSE)
    attr(full_mu, "draw_index") <- NULL

    expect_equal(unname(pred$mu_samples), unname(full_mu[, row_idx, drop = FALSE]))
  }
})

test_that("predict.stLMM_recovery simulates new space-time NNGP nodes across covariance families", {
  set.seed(209)
  coords <- data.frame(
    lon = c(0, 1, 2, 0, 1, 2, 0.5, 1.5),
    lat = c(0, 0, 0, 1, 1, 1, 2, 2),
    time = c(0, 0, 0, 0.5, 0.5, 0.5, 1, 1)
  )
  idx_node <- rep(seq_len(nrow(coords)), each = 2)
  dat <- data.frame(
    y = rep(0, length(idx_node)),
    lon = coords$lon[idx_node],
    lat = coords$lat[idx_node],
    time = coords$time[idx_node]
  )
  newdata <- data.frame(
    lon = c(0.25, 0.25, 1.75),
    lat = c(0.25, 0.25, 1.25),
    time = c(0.25, 0.25, 0.75)
  )

  specs <- list(
    sep_exp = list(
      formula = y ~ 0 + nngp(lon, lat, time, m = 3, cov_model = "sep_exp", ordering = "maxmin"),
      starting = c(sigma_sq = 1, phi = 0.7, lambda = 0.4),
      priors = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5), lambda = uniform(0.1, 5))
    ),
    multi_res_sep_exp = list(
      formula = y ~ 0 + nngp(lon, lat, time, m = 3, cov_model = "multi_res_sep_exp", ordering = "maxmin"),
      starting = c(sigma_sq = 1, alpha = 0.5, phi_1 = 0.7, lambda_1 = 0.4, phi_2 = 1.2, lambda_2 = 0.8),
      priors = list(
        sigma_sq = ig(2, 1), alpha = uniform(0.05, 0.95),
        phi_1 = uniform(0.1, 5), lambda_1 = uniform(0.1, 5),
        phi_2 = uniform(0.1, 5), lambda_2 = uniform(0.1, 5)
      )
    ),
    gneiting = list(
      formula = y ~ 0 + nngp(lon, lat, time, m = 3, cov_model = "gneiting", ordering = "maxmin"),
      starting = c(sigma_sq = 1, a = 0.4, c = 0.7, alpha = 1, beta = 0.6, gamma = 0.5, delta = 0),
      priors = list(
        sigma_sq = ig(2, 1), a = uniform(0.1, 5),
        c = uniform(0.1, 5), alpha = uniform(0.05, 1),
        beta = uniform(0.05, 0.95), gamma = uniform(0.05, 1),
        delta = uniform(0, 2)
      )
    )
  )

  for(spec_name in names(specs)){
    spec <- specs[[spec_name]]
    tuning_i <- spec$starting
    tuning_i[] <- 0
    tuning_i["sigma_sq"] <- 0

    fit <- suppressWarnings(stLMM(
      spec$formula,
      data = dat,
      starting = list(tau_sq = 1, nngp_1 = spec$starting),
      tuning = list(tau_sq = 0, nngp_1 = tuning_i),
      priors = list(tau_sq = ig(2, 1), nngp_1 = spec$priors),
      n_samples = 6,
      verbose = FALSE
    ))
    rec <- recover(fit, sub_sample = list(start = 2, thin = 2))
    pred <- predict(rec, newdata = newdata, st_scale = 2)

    expect_equal(dim(pred$mu_samples), c(3L, 3L))
    expect_equal(pred$mu_samples[, 1], pred$mu_samples[, 2])
  }
})

test_that("predict.stLMM_recovery validates NNGP st_scale", {
  set.seed(210)
  dat <- data.frame(
    y = rep(0, 6),
    lon = c(0, 1, 2, 0, 1, 2),
    lat = c(0, 0, 0, 1, 1, 1),
    time = c(0, 0, 0, 1, 1, 1)
  )

  fit <- suppressWarnings(stLMM(
    y ~ 0 + nngp(lon, lat, time, m = 3, cov_model = "sep_exp"),
    data = dat,
    starting = list(tau_sq = 1, nngp_1 = c(sigma_sq = 1, phi = 0.7, lambda = 0.4)),
    tuning = list(tau_sq = 0, nngp_1 = c(sigma_sq = 0, phi = 0, lambda = 0)),
    priors = list(
      tau_sq = ig(2, 1),
      nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5), lambda = uniform(0.1, 5))
    ),
    n_samples = 4,
    verbose = FALSE
  ))
  rec <- recover(fit, sub_sample = list(start = 2, thin = 2))

  expect_message(
    pred_joint <- predict(rec, newdata = data.frame(lon = 0.5, lat = 0.5, time = 0.5), joint = TRUE),
    "dense covariance"
  )
  expect_equal(dim(pred_joint$mu_samples), c(2L, 1L))

  expect_error(
    predict(rec, newdata = data.frame(lon = 0.5, lat = 0.5, time = 0.5), st_scale = 0),
    "st_scale"
  )
})
