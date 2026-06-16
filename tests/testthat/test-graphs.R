test_that("graph keys preserve coordinate order", {
  spec_a <- list(
    fun = "nngp",
    args = c("lon", "lat", "time"),
    params = list(m = 2, cov_model = "sep_exp", ordering = "coord")
  )
  spec_b <- list(
    fun = "nngp",
    args = c("time", "lon", "lat"),
    params = list(m = 2, cov_model = "sep_exp", ordering = "coord")
  )

  key_a <- stLMM:::build_graph_key(spec_a, stLMM:::process_graph_defaults("nngp"))
  key_b <- stLMM:::build_graph_key(spec_b, stLMM:::process_graph_defaults("nngp"))

  expect_false(identical(key_a, key_b))
})

test_that("NNGP graph reuse requires identical coordinate order", {
  set.seed(4)
  dat <- data.frame(
    y = rnorm(8),
    x = rnorm(8),
    lon = rep(seq_len(4), 2),
    lat = rep(c(1, 2), 4),
    time = rep(c(1, 2), each = 4)
  )

  shared <- suppressWarnings(stLMM(
    y ~ 0 +
      nngp(lon, lat, time, m = 2, cov_model = "sep_exp") +
      x:nngp(lon, lat, time, m = 2, cov_model = "sep_exp"),
    data = dat,
    starting = list(
      tau_sq = 1,
      nngp_1 = c(sigma_sq = 1, phi = 1, lambda = 1),
      nngp_2 = c(sigma_sq = 1, phi = 1, lambda = 1)
    ),
    priors = list(
      tau_sq = ig(2, 1),
      nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5), lambda = uniform(0.1, 5)),
      nngp_2 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5), lambda = uniform(0.1, 5))
    ),
    n_samples = 2,
    verbose = FALSE
  ))

  distinct <- suppressWarnings(stLMM(
    y ~ 0 +
      nngp(lon, lat, time, m = 2, cov_model = "sep_exp") +
      x:nngp(time, lon, lat, m = 2, cov_model = "sep_exp"),
    data = dat,
    starting = list(
      tau_sq = 1,
      nngp_1 = c(sigma_sq = 1, phi = 1, lambda = 1),
      nngp_2 = c(sigma_sq = 1, phi = 1, lambda = 1)
    ),
    priors = list(
      tau_sq = ig(2, 1),
      nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5), lambda = uniform(0.1, 5)),
      nngp_2 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5), lambda = uniform(0.1, 5))
    ),
    n_samples = 2,
    verbose = FALSE
  ))

  expect_equal(shared$term_description$global$n_graphs, 1L)
  expect_equal(distinct$term_description$global$n_graphs, 2L)
})

test_that("NNGP ordering options create expected graph metadata", {
  coords <- cbind(
    lon = c(0, 1, 2, 0, 1, 2, 0.5, 1.5),
    lat = c(0, 0, 0, 1, 1, 1, 2, 2)
  )
  user_order <- c(4L, 1L, 7L, 2L, 5L, 8L, 3L, 6L)
  orderings <- list(
    coord = "coord",
    default = "default",
    maxmin = "maxmin",
    hilbert = "hilbert",
    random = "random",
    user = user_order
  )

  set.seed(301)
  graph_random_1 <- stLMM:::make_nngp_graph(coords, m = 3, ordering = "random")
  set.seed(301)
  graph_random_2 <- stLMM:::make_nngp_graph(coords, m = 3, ordering = "random")
  expect_equal(graph_random_1$ord, graph_random_2$ord)

  for(nm in names(orderings)){
    set.seed(300)
    graph <- stLMM:::make_nngp_graph(coords, m = 3, ordering = orderings[[nm]])

    expect_equal(length(graph$ord), nrow(coords))
    expect_equal(sort(graph$ord), seq_len(nrow(coords)))
    expect_equal(graph$ord_inv[graph$ord], seq_len(nrow(coords)))
    expect_equal(graph$coords_ord, coords[graph$ord, , drop = FALSE])
    expect_equal(length(graph$nnIndxLU), 2L * nrow(coords))
    expect_equal(graph$m, 3L)

    if(nm == "coord")
      expect_equal(graph$ord, do.call(order, as.data.frame(coords)))
    if(nm == "default")
      expect_equal(graph$ord, do.call(order, as.data.frame(coords)))
    if(nm == "user"){
      expect_equal(graph$ordering, "user")
      expect_equal(graph$ord, user_order)
    }
  }

  expect_error(
    stLMM:::make_nngp_graph(coords, m = 3, ordering = c(1, 2, 2, 4, 5, 6, 7, 8)),
    "permutation"
  )
  expect_error(
    stLMM:::make_nngp_graph(coords[, 1, drop = FALSE], m = 3, ordering = "hilbert"),
    "hilbert"
  )
})

test_that("NNGP fast neighbor search matches brute force history-restricted search", {
  set.seed(710)
  coords2 <- cbind(x = runif(60), y = runif(60))
  coords3 <- cbind(x = runif(70), y = runif(70), t = runif(70))
  coords_tie <- cbind(
    x = c(0, 1, -1, 0, 0, 2, -2, 0.5, -0.5, 1.5),
    y = c(0, 0, 0, 1, -1, 0, 0, 0.5, -0.5, 1.5)
  )

  for(coords in list(coords2, coords3, coords_tie)){
    for(m in c(1L, 3L, 8L)){
      fast <- stLMM:::mkNNIndx(coords, m = m, n_omp_threads = 2)
      brute <- stLMM:::mkNNIndxBrute(coords, m = m, n_omp_threads = 2)

      expect_equal(fast$nnIndx, brute$nnIndx)
      expect_equal(fast$nnDist, brute$nnDist)
      expect_equal(fast$nnIndxLU, brute$nnIndxLU)

      for(i in seq_len(nrow(coords))){
        start <- fast$nnIndxLU[i] + 1L
        count <- fast$nnIndxLU[nrow(coords) + i]
        if(count > 0L)
          expect_true(all(fast$nnIndx[start:(start + count - 1L)] < i - 1L))
      }
    }
  }
})

test_that("random NNGP ordering maps observations to retained graph order", {
  set.seed(731)
  n <- 12L
  dat <- data.frame(
    y = rnorm(n),
    lon = runif(n),
    lat = runif(n)
  )

  proc <- stLMM:::build_process_components(
    y ~ 1 + nngp(lon, lat, m = 5, ordering = "random"),
    dat
  )

  graph <- proc$graphs[[1]]
  term <- proc$process_terms[[1]]

  expect_equal(term$map, as.integer(graph$ord_inv))
  expect_equal(term$node_nobs, rep(1L, n))
})

test_that("NNGP graph constructor can select fast or brute neighbor search", {
  set.seed(711)
  coords <- cbind(x = runif(40), y = runif(40))

  graph_fast <- stLMM:::make_nngp_graph(coords, m = 5, nngp_search = "fast")
  graph_brute <- stLMM:::make_nngp_graph(coords, m = 5, nngp_search = "brute")

  expect_equal(graph_fast$nnIndx, graph_brute$nnIndx)
  expect_equal(graph_fast$nnIndxLU, graph_brute$nnIndxLU)
  expect_equal(graph_fast$nngp_search, "fast")
  expect_equal(graph_brute$nngp_search, "brute")
  expect_error(
    stLMM:::make_nngp_graph(coords, m = 5, nngp_search = "approx"),
    "should be one of"
  )
})

test_that("process graphs reject non-finite coordinates", {
  coords <- cbind(
    x = c(0, 1, NA, 3),
    y = c(0, 1, 2, 3)
  )

  expect_error(
    stLMM:::make_nngp_graph(coords, m = 2),
    "coordinates must be finite"
  )
  expect_error(
    stLMM:::make_gp_graph(coords),
    "coordinates must be finite"
  )

  dat <- data.frame(
    y = rnorm(4),
    x = coords[, 1],
    z = coords[, 2]
  )
  expect_error(
    stLMM(
      y ~ nngp(x, z, m = 2),
      data = dat,
      priors = list(
        tau_sq = ig(2, 1),
        nngp_1 = list(sigma_sq = ig(2, 1), phi = uniform(0.1, 5))
      ),
      n_samples = 2,
      verbose = FALSE
    ),
    "process coordinates must be finite"
  )
})

test_that("GP and NNGP coordinate support uses consistent exact row keys", {
  a <- 1 / 3
  b <- 1 - 2 / 3
  coords <- matrix(c(a, b), ncol = 1L)

  expect_false(a == b)
  expect_equal(nrow(stLMM:::unique_coordinate_rows(coords)), 2L)
  expect_equal(stLMM:::match_unique_rows(coords, stLMM:::unique_coordinate_rows(coords)), c(1L, 2L))

  dat <- data.frame(y = c(0, 1), x = c(a, b))
  nngp_comp <- stLMM:::build_nngp_components(
    list(name = "nngp_1", args = "x", params = list(m = 1, ordering = "coord")),
    dat
  )
  gp_comp <- stLMM:::build_gp_components(
    list(name = "gp_1", args = "x", params = list()),
    dat
  )

  expect_equal(nrow(nngp_comp$graph$coords), 2L)
  expect_equal(nngp_comp$term$map, c(1L, 2L))
  expect_equal(nrow(gp_comp$graph$coords), 2L)
  expect_equal(gp_comp$term$map, c(1L, 2L))
})
