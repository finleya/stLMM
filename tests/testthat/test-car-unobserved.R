test_that("CAR fits retain graph areas absent from the data", {
  ids <- paste0("a", 1:7)
  adj <- Matrix::sparseMatrix(
    i = c(1:6, 2:7), j = c(2:7, 1:6), x = 1,
    dims = c(7, 7), dimnames = list(ids, ids)
  )
  g <- car_graph(adj)
  dat <- data.frame(
    area_id = ids[c(2, 2, 4, 6, 4)],
    x = c(0.2, 0.8, 1.1, -0.4, 0.6),
    y = c(0.4, 0.7, 1.2, -0.1, 0.9)
  )
  # Include missing areas at the beginning, middle, and end of the graph.
  padded <- rbind(dat, data.frame(
    area_id = ids[c(1, 3, 5, 7)], x = 0, y = NA_real_
  ))
  newdata <- data.frame(area_id = c(ids, ids[7]), x = 0.5)

  for(car_model in c("proper", "leroux")){
    for(svc in c(FALSE, TRUE)){
      formula <- if(svc)
        y ~ x + x:car(area_id, graph = g, car_model = car_model)
      else
        y ~ x + car(area_id, graph = g, car_model = car_model)

      fit_data <- function(data){
        set.seed(731)
        stLMM(
          formula, data = data,
          starting = list(tau_sq = 0.3, car_1 = c(sigma_sq = 0.8, rho = 0.4)),
          tuning = list(tau_sq = 0.1, car_1 = c(sigma_sq = 0.1, rho = 0.1)),
          priors = list(
            tau_sq = ig(2, 0.3),
            car_1 = list(sigma_sq = ig(2, 0.8), rho = uniform(0.05, 0.95))
          ),
          n_samples = 50, warmup = FALSE, verbose = FALSE
        )
      }
      fit <- fit_data(dat)
      fit_padded <- fit_data(padded)

      expect_equal(fit$backend$n_obs, 5L)
      expect_equal(fit$backend$n_missing_response, 0L)
      expect_equal(fit$backend$process_terms[[1]]$n_node, 7L)
      expect_equal(fit$backend$process_terms[[1]]$node_nobs,
                   c(0L, 2L, 0L, 2L, 0L, 1L, 0L))
      expect_equal(as_samples(fit, metadata = FALSE),
                   as_samples(fit_padded, metadata = FALSE), tolerance = 1e-12)

      set.seed(732)
      rec <- recover(fit)
      set.seed(732)
      rec_padded <- recover(fit_padded)
      expect_equal(rec$w_samples, rec_padded$w_samples, tolerance = 1e-12)
      expect_equal(ncol(rec$w_samples$car_1), 7L)

      set.seed(733)
      pred <- predict(rec, newdata = newdata, y_samples = TRUE)
      set.seed(733)
      pred_padded <- predict(rec_padded, newdata = newdata, y_samples = TRUE)
      expect_equal(pred$mu_samples, pred_padded$mu_samples, tolerance = 1e-12)
      expect_equal(pred$y_samples, pred_padded$y_samples, tolerance = 1e-12)
      expect_true(all(is.finite(pred$y_samples)))
      expect_equal(pred$mu_samples[, 7], pred$mu_samples[, 8])

      bad <- dat
      bad$area_id[1] <- "outside_graph"
      expect_error(fit_data(bad), "area value.*not found in graph")
      bad <- newdata
      bad$area_id[1] <- "outside_graph"
      expect_error(predict(rec, newdata = bad), "new CAR area value")
    }
  }
})
