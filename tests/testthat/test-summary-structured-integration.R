test_that("iid summaries, recovery, and prediction compose with structured process families", {
  set.seed(46601)
  ids <- paste0("area", 1:6)
  adj <- Matrix::sparseMatrix(i = c(1:5, 2:6), j = c(2:6, 1:5), x = 1,
    dims = c(6, 6), dimnames = list(ids, ids))
  g <- car_graph(adj)
  d <- expand.grid(area = ids, time = 1:4, stringsAsFactors = FALSE)
  area <- match(d$area, ids)
  d$lon <- (area - 1) %% 3
  d$lat <- (area - 1) %/% 3
  d$parent <- ifelse(area <= 3, "north", "south")
  d$child <- ifelse(d$time %% 2 == 0, "plot1", "plot2")
  d$x <- rnorm(nrow(d)); d$x[c(1, 7)] <- 0
  d$y <- rnorm(nrow(d))
  fixed_cov <- function(...) lapply(list(...), fixed)
  specs <- list(
    ar1 = list(term = "ar1(time)", start = list(ar1_1 = fixed_cov(sigma_sq=.3, phi=.5))),
    gp = list(term = 'gp(lon, lat)', start = list(gp_1 = fixed_cov(sigma_sq=.3, phi=.8))),
    gp_matern = list(term = 'gp(lon, lat, cov_model = "matern")',
                     start = list(gp_1 = fixed_cov(sigma_sq=.3, phi=.8, nu=1.5))),
    nngp = list(term = 'nngp(lon, lat, m = 3, ordering = "maxmin")',
                start = list(nngp_1 = fixed_cov(sigma_sq=.3, phi=.8))),
    nngp_matern = list(term = 'nngp(lon, lat, m = 3, ordering = "maxmin", cov_model = "matern")',
                       start = list(nngp_1 = fixed_cov(sigma_sq=.3, phi=.8, nu=1.5))),
    nngp_slope = list(term = 'x:nngp(lon, lat, m = 3, ordering = "maxmin")',
                      start = list(nngp_1 = fixed_cov(sigma_sq=.3, phi=.8)), slope = "nngp_1"),
    car = list(term = 'car(area, graph = g)', start = list(car_1 = fixed_cov(sigma_sq=.3, rho=.5))),
    leroux = list(term = 'car(area, graph = g, car_model = "leroux")',
                  start = list(car_1 = fixed_cov(sigma_sq=.3, rho=.5))),
    dagar = list(term = 'dagar(area, graph = g, ordering = c(6, 4, 2, 5, 3, 1))',
                 start = list(dagar_1 = fixed_cov(sigma_sq=.3, rho=.5))),
    car_time = list(term = 'car_time(area, time, graph = g)',
                    start = list(car_time_1 = fixed_cov(sigma_sq=.3, rho=.5, phi=.4))),
    dagar_time = list(term = 'dagar_time(area, time, graph = g, time_model = "exp")',
                      start = list(dagar_time_1 = fixed_cov(sigma_sq=.3, rho=.5, lambda=.8))),
    multiple = list(term = 'ar1(time) + x:nngp(lon, lat, m = 3, ordering = "maxmin")',
                    start = list(ar1_1 = fixed_cov(sigma_sq=.2, phi=.4),
                                 nngp_1 = fixed_cov(sigma_sq=.3, phi=.8)), slope = "nngp_1")
  )
  # Repeated supports are intentional. Suppress only their documented notices.
  support_notice <- function(expr) withCallingHandlers(expr, warning = function(w){
    if(grepl("duplicated (time values|coordinates) detected", conditionMessage(w)))
      invokeRestart("muffleWarning")
  })
  for(collapse_iid in c(FALSE, TRUE)) for(family in c("gaussian", "binomial", "probit", "negative_binomial")){
    dat <- d
    if(family != "gaussian") dat$y <- as.numeric(dat$y > 0)
    dat$y[1L] <- NA
    for(name in names(specs)){
      spec <- specs[[name]]
      label <- paste(family, name, collapse_iid)
      formula <- as.formula(paste("y ~ x + iid(parent) + iid(child, ~ 1 + x, by = parent) +", spec$term))
      fit <- support_notice(stLMM(formula, dat, family = family,
        trials = if(family == "binomial") rep(1, nrow(dat)) else NULL,
        size = if(family == "negative_binomial") 2 else NULL,
        starting = c(if(family == "gaussian") list(resid = list(tau_sq = fixed(.4))), spec$start),
        priors = list(iid_1 = list(sigma_sq = ig(2, 1)),
                      iid_2 = list(covariance = iw(4, diag(2)))),
        n_samples = 60, chains = 2, warmup = FALSE, verbose = FALSE, collapse_iid = collapse_iid))
      rec <- recover(fit)
      rng <- .Random.seed
      s <- summary(rec, type = "combined", terms = c("iid_1", "iid_2"), burn = 5,
                    diagnostics = FALSE)
      expect_identical(.Random.seed, rng, info = label)
      expect_equal(nrow(s$parameters), 8L, info = label)
      nd <- dat[, setdiff(names(dat), "y"), drop = FALSE]
      nd$x <- 0
      p0 <- predict(rec, nd, scale = "link", sub_sample = list(start = 6), return_w_samples = TRUE)
      nd$x <- 1
      p1 <- predict(rec, nd, scale = "link", sub_sample = list(start = 6), return_w_samples = TRUE)
      remove_process <- function(pred, x){
        lapply(pred$chains, function(ch){
          ans <- ch$mu_samples
          for(term in names(ch$w_samples))
            ans <- ans - ch$w_samples[[term]] * if(term %in% spec$slope) x else 1
          ans
        })
      }
      a <- remove_process(p0, 0); b <- remove_process(p1, 1)
      for(i in seq_len(nrow(s$parameters))){
        row <- s$quantities$data_row[i]
        draws <- unlist(lapply(1:2, function(ch){
          if(s$quantities$coefficient[i] == "(Intercept)") a[[ch]][, row] else
            b[[ch]][, row] - a[[ch]][, row]
        }))
        expected <- c(mean(draws), sd(draws), quantile(draws, c(.025,.5,.975), names = FALSE))
        expect_equal(unname(s$parameters[i, ]), expected, tolerance = 1e-11, info = label)
      }
      thin <- recover(fit, sub_sample = list(start = 4, thin = 3))
      thin_summary <- summary(thin, type = "combined", terms = c("iid_1", "iid_2"),
                              burn = 5, diagnostics = FALSE)
      if(collapse_iid && family == "gaussian"){
        expect_identical(thin_summary$iteration, seq.int(7L, 60L, 3L), info = label)
        expect_equal(thin_summary$n_used, 18L, info = label)
      } else expect_identical(s, thin_summary, info = label)
      pred <- predict(thin, scale = "link")
      full <- fitted(thin, summary = FALSE, scale = "link")
      for(ch in 1:2){
        expect_equal(pred$chains[[ch]]$draw_index, attr(full[[ch]], "draw_index"), info = label)
        expect_equal(as.numeric(pred$chains[[ch]]$mu_samples), as.numeric(full[[ch]]), info = label)
      }
      # Shared new-child draws at an existing parent and process support.
      nd <- dat[rep(2L, 2L), setdiff(names(dat), "y"), drop = FALSE]
      nd$child <- "new_child"
      expect_error(predict(thin, nd, scale = "link"), "allow_new_levels", info = label)
      new <- predict(thin, nd, allow_new_levels = TRUE, scale = "link")
      for(ch in new$chains) expect_identical(ch$mu_samples[,1], ch$mu_samples[,2], info = label)
      expect_error(predict(thin, transform(nd, parent = "new_parent"), allow_new_levels = TRUE),
                    "new parent", info = label)
      expect_error(predict(thin, transform(nd, child = NA_character_), allow_new_levels = TRUE),
                    "missing", info = label)
      expect_error(predict(thin, transform(nd, parent = NA_character_), allow_new_levels = TRUE),
                    "missing", info = label)
      if(name %in% c("gp", "gp_matern", "nngp", "nngp_matern", "nngp_slope", "multiple")){
        nd$lon <- nd$lon + .17; nd$lat <- nd$lat + .13
        nd <- nd[c(1,2,1), , drop = FALSE]
        nd$lon[3] <- nd$lon[3] + .23
        new <- predict(thin, nd, allow_new_levels = TRUE, scale = "link", joint = TRUE)
        for(ch in new$chains){
          expect_true(all(is.finite(ch$mu_samples)), info = label)
          expect_gt(sd(ch$mu_samples[,1] - ch$mu_samples[,3]), 0, label = label)
          expect_identical(ch$mu_samples[,1], ch$mu_samples[,2], info = label)
        }
      }
    }
  }
})
