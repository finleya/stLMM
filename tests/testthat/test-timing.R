test_that("phase timings distinguish calibration from returned sampling", {
  set.seed(731)
  dat <- data.frame(y = rnorm(12), x = rnorm(12))
  for (warmup in list(FALSE, list(batch_length = 3, min_batches = 2, max_batches = 2))) {
    fit <- stLMM(y ~ x, data = dat, n_samples = 5, warmup = warmup,
                 chain_control = list(seed = 19), verbose = FALSE)
    expect_identical(names(fit$timing), c("backend", "phases"))
    expect_identical(names(fit$timing$phases), c("initialization", "warmup", "sampling"))
    expect_true(all(is.finite(fit$timing$phases)))
    expect_true(all(fit$timing$phases >= 0))
    expect_gt(fit$timing$phases[["sampling"]], 0)
    # proc.time has coarser resolution than the native steady clock.
    expect_lte(sum(fit$timing$phases), fit$timing$backend[["elapsed"]] + 0.02)
    expect_equal(nrow(fit$beta_samples), 5L)
    expected_batches <- if (identical(warmup, FALSE)) 0L else 2L
    expect_equal(fit$adaptive_metropolis$warmup$n_batches, expected_batches)
    expect_output(print(fit), "post-warmup sampling elapsed")
    expect_output(print(fit), "backend elapsed")
  }
})

test_that("multi-chain phase durations preserve chain identity and sum correctly", {
  set.seed(732)
  dat <- data.frame(y = rnorm(16), lon = runif(16), lat = runif(16))
  fit <- stLMM(y ~ nngp(lon, lat, m = 3), data = dat, n_samples = 4,
               priors = list(resid = list(tau_sq = half_t(3, 1)),
                             nngp_1 = list(sigma_sq = half_t(3, 1), phi = uniform(1, 20))),
               chains = 2, chain_control = list(seed = 20),
               warmup = FALSE, verbose = FALSE)
  expect_identical(rownames(fit$timing$phases_by_chain), names(fit$chains))
  expect_identical(colnames(fit$timing$phases_by_chain),
                   c("initialization", "warmup", "sampling"))
  for (i in seq_along(fit$chains))
    expect_equal(fit$timing$phases_by_chain[i, ], fit$chains[[i]]$timing$phases)
  expect_equal(fit$timing$phases_total, colSums(fit$timing$phases_by_chain))
  expect_equal(unname(fit$timing$backend_total[["elapsed"]]),
               sum(fit$timing$backend_by_chain[, "elapsed"]))
  expect_output(print(fit), "post-warmup sampling elapsed")
})
