process_sigma_sq_draws <- function(object, term, draw_index){
  if(is.null(object$sigma_sq_samples))
    stop("error: process variance samples are missing from fitted object")

  sigma_names <- colnames(object$sigma_sq_samples)
  sigma_col <- paste0(term$name, "_sigma_sq")

  if(!is.null(sigma_names) && sigma_col %in% sigma_names)
    return(object$sigma_sq_samples[draw_index, sigma_col])

  ## Backward compatibility for older fitted objects saved before process
  ## variance columns were named <term>_sigma_sq.
  if(!is.null(sigma_names) && term$name %in% sigma_names)
    return(object$sigma_sq_samples[draw_index, term$name])

  term_names <- vapply(object$backend$process_terms, `[[`, character(1), "name")
  term_index <- match(term$name, term_names)
  if(!is.na(term_index) && term_index <= ncol(object$sigma_sq_samples))
    return(object$sigma_sq_samples[draw_index, term_index])

  stop("error: process variance samples are missing for ", term$name)
}

resolve_n_omp_threads <- function(object, n_omp_threads = NULL){
  if(is.null(n_omp_threads))
    n_omp_threads <- object$backend$n_omp_threads %||% 1L

  if(!is.numeric(n_omp_threads) || length(n_omp_threads) != 1L ||
     is.na(n_omp_threads) || !is.finite(n_omp_threads) ||
     n_omp_threads < 1L || abs(n_omp_threads - round(n_omp_threads)) > 0)
    stop("error: n_omp_threads must be a positive integer")

  as.integer(n_omp_threads)
}

stlmm_progress <- function(verbose, ...){
  if(isTRUE(verbose))
    message(...)
}

prediction_residual_sd <- function(object, newdata, n0, where){
  residual_model <- object$backend$residual_model
  if(is.null(residual_model) || identical(residual_model$type, "global_tau"))
    return(NULL)
  if(!identical(residual_model$type, "fixed_variance"))
    stop("error: unsupported residual model type in fitted object")

  if(is.null(newdata)){
    variance <- residual_model$variance
    if(length(variance) != n0)
      stop("error: fitted fixed residual variance has incompatible length")
  } else {
    dat <- if(is.data.frame(newdata)) newdata else as.data.frame(newdata)
    residual <- structure(
      list(
        type = "fixed_variance",
        variance_expr = residual_model$variance_expr,
        variance_label = residual_model$variance_label,
        env = parent.frame()
      ),
      class = "stLMM_residual"
    )
    variance <- eval_residual_variance(
      residual = residual,
      data = dat,
      n_expected = n0,
      where = where
    )
  }

  if(anyNA(variance) || any(!is.finite(variance)) || any(variance <= 0))
    stop("error: fixed residual variance must be finite and positive for ", where)

  sqrt(as.double(variance))
}

prediction_residual_sd_samples <- function(object, newdata, n0, draw_index, where){
  residual_model <- object$backend$residual_model
  if(is.null(residual_model) || identical(residual_model$type, "global_tau"))
    return(NULL)
  if(identical(residual_model$type, "fixed_variance")){
    sd <- prediction_residual_sd(object, newdata, n0, where)
    return(matrix(rep(sd, each = length(draw_index)), nrow = length(draw_index)))
  }
  if(identical(residual_model$type, "scaled_variance")){
    if(is.null(object$residual_variance_samples))
      stop("error: residual variance samples are missing from fitted object")
    if(is.null(newdata)){
      vhat <- residual_model$vhat
      weight <- residual_model$weight
      if(length(vhat) != n0 || length(weight) != n0)
        stop("error: fitted scaled residual metadata has incompatible length")
    } else {
      dat <- if(is.data.frame(newdata)) newdata else as.data.frame(newdata)
      vhat <- eval_residual_expr(residual_model$vhat_expr, dat, parent.frame())
      if(!is.numeric(vhat) || length(vhat) != n0)
        stop("error: prediction residual vhat must be numeric with length ", n0)
      if(anyNA(vhat) || any(!is.finite(vhat)) || any(vhat <= 0))
        stop("error: prediction residual vhat must be finite and positive for ", where)
      if(isTRUE(residual_model$has_n)){
        n_eff <- eval_residual_expr(residual_model$n_expr, dat, parent.frame())
        if(!is.numeric(n_eff) || length(n_eff) != n0)
          stop("error: prediction residual n must be numeric with length ", n0)
        if(anyNA(n_eff) || any(!is.finite(n_eff)) || any(n_eff <= 1))
          stop("error: prediction residual n must be finite and greater than 1 for ", where)
        weight <- n_eff / (n_eff + residual_model$shrinkage)
      } else {
        weight <- rep(1, n0)
      }
    }

    params <- object$residual_variance_samples[draw_index, , drop = FALSE]
    if(ncol(params) == 1L){
      variance <- params[, 1L, drop = FALSE] %*% t(as.double(vhat))
    } else {
      kappa <- params[, 1L]
      tau0_sq <- params[, 2L]
      base <- exp(matrix(rep(weight * log(vhat), each = length(draw_index)), nrow = length(draw_index)) +
                  matrix(rep((1 - weight), each = length(draw_index)), nrow = length(draw_index)) *
                    log(tau0_sq))
      variance <- base * kappa
    }
    if(anyNA(variance) || any(!is.finite(variance)) || any(variance <= 0))
      stop("error: scaled residual variance must be finite and positive for ", where)
    return(sqrt(variance))
  }
  if(!identical(residual_model$type, "group_ig_variance"))
    stop("error: unsupported residual model type in fitted object")
  if(is.null(object$residual_variance_samples))
    stop("error: residual variance samples are missing from fitted object")

  if(is.null(newdata)){
    group_index <- residual_model$group_index_full
    if(length(group_index) != n0)
      stop("error: fitted residual group index has incompatible length")
  } else {
    dat <- if(is.data.frame(newdata)) newdata else as.data.frame(newdata)
    group <- eval_residual_expr(residual_model$group_expr, dat, parent.frame())
    group_index <- match(as.character(group), residual_model$groups)
    if(length(group_index) != n0)
      stop("error: prediction residual group has incompatible length")
  }
  if(anyNA(group_index))
    stop("error: prediction rows include residual groups not present in the fitted residual model")

  variance <- object$residual_variance_samples[draw_index, group_index, drop = FALSE]
  if(anyNA(variance) || any(!is.finite(variance)) || any(variance <= 0))
    stop("error: sampled residual variance must be finite and positive for ", where)
  sqrt(variance)
}

prediction_trials <- function(object, newdata, n0){
  if(!identical(object$backend$family, "binomial"))
    return(NULL)

  if(is.null(newdata)){
    trials <- object$backend$trials
    if(is.null(trials))
      trials <- object$backend$likelihood$trials
    if(is.null(trials))
      trials <- rep.int(1L, n0)
  } else {
    dat <- if(is.data.frame(newdata)) newdata else as.data.frame(newdata)
    trials <- if("trials" %in% names(dat)) dat$trials else rep.int(1L, n0)
  }

  if(length(trials) != n0)
    stop("error: prediction trials must have length matching prediction rows")
  if(!is.numeric(trials) && !is.integer(trials))
    stop("error: prediction trials must be numeric or integer")
  trials_num <- as.numeric(trials)
  if(anyNA(trials_num) || any(!is.finite(trials_num)) || any(trials_num <= 0))
    stop("error: prediction trials must be finite positive counts")
  if(any(abs(trials_num - round(trials_num)) > sqrt(.Machine$double.eps)))
    stop("error: prediction trials must be integer-valued counts")

  as.integer(round(trials_num))
}

prediction_nb_size <- function(object){
  size <- object$backend$nb_size
  if(is.null(size))
    size <- object$backend$likelihood$size
  if(is.null(size) || length(size) != 1L || is.na(size) ||
     !is.finite(size) || size <= 0)
    stop("error: fitted negative-binomial size is missing or invalid")
  as.double(size)
}

finish_prediction_samples <- function(object,
                                      eta,
                                      y_samples,
                                      scale,
                                      trials,
                                      draw_index,
                                      newdata,
                                      joint,
                                      joint_method = NULL,
                                      w_samples = NULL,
                                      gaussian_y_sampler){
  n_draw <- nrow(eta)
  n0 <- ncol(eta)

  if(identical(object$backend$family, "binomial")){
    prob <- stats::plogis(eta)
    mu <- if(identical(scale, "response")) prob else eta
    out <- list(
      mu_samples = mu,
      y_samples = NULL,
      draw_index = as.integer(draw_index),
      newdata = !is.null(newdata),
      joint = isTRUE(joint),
      joint_method = joint_method,
      scale = scale,
      w_samples = w_samples
    )
    if(isTRUE(y_samples)){
      y <- matrix(0L, nrow = n_draw, ncol = n0)
      colnames(y) <- colnames(eta)
      for(i in seq_len(n_draw))
        y[i, ] <- stats::rbinom(n0, size = trials, prob = prob[i, ])
      out$y_samples <- y
    }
    class(out) <- "stLMM_prediction"
    return(out)
  }

  if(identical(object$backend$family, "negative_binomial")){
    mean <- exp(eta)
    mu <- if(identical(scale, "response")) mean else eta
    out <- list(
      mu_samples = mu,
      y_samples = NULL,
      draw_index = as.integer(draw_index),
      newdata = !is.null(newdata),
      joint = isTRUE(joint),
      joint_method = joint_method,
      scale = scale,
      w_samples = w_samples
    )
    if(isTRUE(y_samples)){
      size <- prediction_nb_size(object)
      y <- matrix(0L, nrow = n_draw, ncol = n0)
      colnames(y) <- colnames(eta)
      for(i in seq_len(n_draw))
        y[i, ] <- stats::rnbinom(n0, size = size, mu = mean[i, ])
      out$y_samples <- y
    }
    class(out) <- "stLMM_prediction"
    return(out)
  }

  out <- list(
    mu_samples = eta,
    y_samples = NULL,
    draw_index = as.integer(draw_index),
    newdata = !is.null(newdata),
    joint = isTRUE(joint),
    joint_method = joint_method,
    w_samples = w_samples
  )

  if(isTRUE(y_samples))
    out$y_samples <- gaussian_y_sampler(eta)

  class(out) <- "stLMM_prediction"
  out
}

predict.stLMM <- function(object,
                         newdata = NULL,
                         y_samples = FALSE,
                         n_omp_threads = NULL,
                         verbose = FALSE,
                         sub_sample = list(start = 1L, thin = 1L),
                         scale = c("response", "link"),
                         ...){

  if(!is.list(object) || is.null(object$backend))
    stop("error: object must be an stLMM fit")
  n_omp_threads <- resolve_n_omp_threads(object, n_omp_threads)
  if(length(object$backend$process_terms) > 0L){
    if(!is.null(object$recover_iter) && length(object$recover_iter) > 0L &&
       !is.null(object$w_samples))
      return(predict.stLMM_recovery(
        object,
        newdata = newdata,
        y_samples = y_samples,
        n_omp_threads = n_omp_threads,
        verbose = verbose,
        sub_sample = sub_sample,
        scale = scale,
        ...
      ))
    stop("error: prediction with process terms requires saved or recovered latent process samples; call recover() on the fitted object first")
  }
  if(!is.list(sub_sample))
    stop("error: sub_sample must be a list with optional entries 'start' and 'thin'")

  scale <- match.arg(scale)
  start <- sub_sample$start
  thin <- sub_sample$thin
  if(is.null(start)) start <- 1L
  if(is.null(thin)) thin <- 1L

  if(!is.numeric(start) || length(start) != 1L || is.na(start) || start < 1)
    stop("error: sub_sample$start must be a positive integer")
  if(!is.numeric(thin) || length(thin) != 1L || is.na(thin) || thin < 1)
    stop("error: sub_sample$thin must be a positive integer")

  start <- as.integer(start)
  thin <- as.integer(thin)
  n_samples <- if(!is.null(object$beta_samples)) nrow(object$beta_samples) else length(object$tau_sq_samples)
  draw_index <- seq.int(start, n_samples, by = thin)
  if(length(draw_index) == 0L)
    stop("error: sub_sample selects no posterior draws")

  if(is.null(newdata)){
    X0 <- object$backend$X
    Z0 <- object$backend$Z
    offset0 <- object$backend$offset
    if(is.null(offset0))
      offset0 <- rep(0, nrow(X0))
    n0 <- as.integer(object$backend$n)
    row_names <- rownames(X0)
  } else {
    stlmm_progress(verbose, "predict: building prediction design matrices")
    pred_backend <- build_existing_support_prediction_backend(
      object,
      newdata,
      st_scale = NULL,
      n_omp_threads = n_omp_threads,
      verbose = verbose
    )
    X0 <- pred_backend$X
    Z0 <- pred_backend$Z
    offset0 <- pred_backend$offset
    n0 <- nrow(X0)
    row_names <- rownames(X0)
  }

  n_draw <- length(draw_index)
  mu <- matrix(0.0, nrow = n_draw, ncol = n0)
  colnames(mu) <- row_names
  stlmm_progress(
    verbose,
    "predict: assembling ", n_draw, " posterior draw(s) for ", n0,
    " prediction row(s) using ", n_omp_threads, " OpenMP thread(s)"
  )

  if(!is.null(object$beta_samples) && ncol(object$beta_samples) > 0L){
    stlmm_progress(verbose, "predict: adding fixed effects")
    beta_draws <- object$beta_samples[draw_index, , drop = FALSE]
    mu <- mu + beta_draws %*% t(X0)
  }

  if(!is.null(object$alpha_samples) && ncol(object$alpha_samples) > 0L){
    stlmm_progress(verbose, "predict: adding grouped random effects")
    alpha_draws <- object$alpha_samples[draw_index, , drop = FALSE]
    mu <- mu + as.matrix(alpha_draws %*% Matrix::t(Z0))
  }

  if(!is.null(offset0)){
    if(length(offset0) != n0)
      stop("error: prediction offset length does not match prediction rows")
    mu <- sweep(mu, 2L, as.numeric(offset0), `+`)
  }

  trials <- prediction_trials(object, newdata, n0)
  gaussian_y_sampler <- function(eta){
    stlmm_progress(verbose, "predict: simulating observation-level predictive samples")
    residual_sd <- prediction_residual_sd_samples(
      object = object,
      newdata = newdata,
      n0 = n0,
      draw_index = draw_index,
      where = "prediction rows"
    )
    y <- eta
    if(is.null(residual_sd)){
      tau <- sqrt(object$tau_sq_samples[draw_index])
      for(i in seq_len(n_draw))
        y[i, ] <- stats::rnorm(n0, mean = eta[i, ], sd = tau[i])
    } else {
      for(i in seq_len(n_draw))
        y[i, ] <- stats::rnorm(n0, mean = eta[i, ], sd = residual_sd[i, ])
    }
    y
  }

  out <- finish_prediction_samples(
    object = object,
    eta = mu,
    y_samples = y_samples,
    scale = scale,
    trials = trials,
    draw_index = draw_index,
    newdata = newdata,
    joint = FALSE,
    gaussian_y_sampler = gaussian_y_sampler
  )
  stlmm_progress(verbose, "predict: done")
  out
}

predict.stLMM_chains <- function(object,
                                 newdata = NULL,
                                 y_samples = FALSE,
                                 n_omp_threads = NULL,
                                 verbose = FALSE,
                                 sub_sample = list(start = 1L, thin = 1L),
                                 scale = c("response", "link"),
                                 ...){

  scale <- match.arg(scale)
  pred <- vector("list", length(object$chains))
  for(i in seq_along(object$chains)){
    stlmm_progress(verbose, "predict: chain ", i, " of ", length(object$chains))
    pred[[i]] <- predict(
      object$chains[[i]],
      newdata = newdata,
      y_samples = y_samples,
      n_omp_threads = n_omp_threads,
      verbose = verbose,
      sub_sample = sub_sample,
      scale = scale,
      ...
    )
  }
  out <- list(
    chains = pred,
    n_chains = object$n_chains,
    call = match.call(),
    fit_call = object$call
  )
  class(out) <- "stLMM_prediction_chains"
  out
}

predict.stLMM_recovery <- function(object,
                                  newdata = NULL,
                                  y_samples = FALSE,
                                  joint = FALSE,
                                  joint_method = c("full", "vecchia"),
                                  pred_m = NULL,
                                  pred_ordering = "maxmin",
                                  st_scale = NULL,
                                  return_w_samples = TRUE,
                                  n_omp_threads = NULL,
                                  verbose = FALSE,
                                  sub_sample = list(start = 1L, thin = 1L),
                                  scale = c("response", "link"),
                                  ...){

  predict_timing <- identical(Sys.getenv("STLMM_PREDICT_TIMING"), "1")

  if(!is.list(object) || is.null(object$backend))
    stop("error: object must be a recovered stLMM fit")
  n_omp_threads <- resolve_n_omp_threads(object, n_omp_threads)

  w_samples_ordered <- object$w_samples_ordered
  if(is.null(w_samples_ordered))
    w_samples_ordered <- object$w_samples

  if(is.null(object$recover_iter) || length(object$recover_iter) == 0L ||
     is.null(w_samples_ordered) || !is.list(w_samples_ordered))
    stop("error: prediction with process terms requires saved or recovered latent process samples; call recover() on the fitted object first")
  if(!identical(joint, FALSE) && !identical(joint, TRUE))
    stop("error: joint must be TRUE or FALSE")
  if(!identical(return_w_samples, FALSE) && !identical(return_w_samples, TRUE))
    stop("error: return_w_samples must be TRUE or FALSE")
  joint_method <- match.arg(joint_method)
  if(!isTRUE(joint))
    joint_method <- "none"

  if(!is.null(pred_m)){
    if(!is.numeric(pred_m) || length(pred_m) != 1L || is.na(pred_m) ||
       pred_m < 1L || abs(pred_m - round(pred_m)) > 0)
      stop("error: pred_m must be a positive integer")
    pred_m <- as.integer(pred_m)
  }

  if(!is.list(sub_sample))
    stop("error: sub_sample must be a list with optional entries 'start' and 'thin'")

  scale <- match.arg(scale)
  start <- sub_sample$start
  thin <- sub_sample$thin
  if(is.null(start)) start <- 1L
  if(is.null(thin)) thin <- 1L

  if(!is.numeric(start) || length(start) != 1L || is.na(start) || start < 1)
    stop("error: sub_sample$start must be a positive integer")
  if(!is.numeric(thin) || length(thin) != 1L || is.na(thin) || thin < 1)
    stop("error: sub_sample$thin must be a positive integer")

  start <- as.integer(start)
  thin <- as.integer(thin)
  recover_row <- which(object$recover_iter >= start & ((object$recover_iter - start) %% thin) == 0L)
  if(length(recover_row) == 0L)
    stop("error: sub_sample selects no recovered latent process draws")

  draw_index <- as.integer(object$recover_iter[recover_row])

  if(is.null(newdata)){
    X0 <- object$backend$X
    Z0 <- object$backend$Z
    offset0 <- object$backend$offset
    if(is.null(offset0))
      offset0 <- rep(0, nrow(X0))
    process_maps <- lapply(object$backend$process_terms, function(term){
      list(map = as.integer(term$map), scale = if(is.null(term$x)) NULL else as.numeric(term$x))
    })
    n0 <- as.integer(object$backend$n)
    row_names <- rownames(X0)
  } else {
    stlmm_progress(
      verbose,
      "predict: building prediction backend using ", n_omp_threads,
      " OpenMP thread(s)"
    )
    if(predict_timing)
      t0 <- proc.time()[["elapsed"]]
    pred_backend <- build_existing_support_prediction_backend(
      object,
      newdata,
      st_scale = st_scale,
      joint_method = joint_method,
      pred_m = pred_m,
      pred_ordering = pred_ordering,
      n_omp_threads = n_omp_threads,
      verbose = verbose
    )
    if(predict_timing)
      message("predict timing: prediction backend = ", round(proc.time()[["elapsed"]] - t0, 3), " sec")
    X0 <- pred_backend$X
    Z0 <- pred_backend$Z
    offset0 <- pred_backend$offset
    process_maps <- pred_backend$process_maps
    n0 <- nrow(X0)
    row_names <- rownames(X0)
  }

  n_draw <- length(draw_index)
  mu <- matrix(0.0, nrow = n_draw, ncol = n0)
  colnames(mu) <- row_names
  stlmm_progress(
    verbose,
    "predict: assembling ", n_draw, " recovered draw(s) for ", n0,
    " prediction row(s) using ", n_omp_threads, " OpenMP thread(s)"
  )

  if(!is.null(object$beta_samples) && ncol(object$beta_samples) > 0L){
    stlmm_progress(verbose, "predict: adding fixed effects")
    if(predict_timing)
      t0 <- proc.time()[["elapsed"]]
    beta_draws <- object$beta_samples[draw_index, , drop = FALSE]
    mu <- mu + beta_draws %*% t(X0)
    if(predict_timing)
      message("predict timing: fixed-effect assembly = ", round(proc.time()[["elapsed"]] - t0, 3), " sec")
  }

  if(!is.null(object$alpha_samples) && ncol(object$alpha_samples) > 0L){
    stlmm_progress(verbose, "predict: adding grouped random effects")
    if(predict_timing)
      t0 <- proc.time()[["elapsed"]]
    alpha_draws <- object$alpha_samples[draw_index, , drop = FALSE]
    mu <- mu + as.matrix(alpha_draws %*% Matrix::t(Z0))
    if(predict_timing)
      message("predict timing: random-effect assembly = ", round(proc.time()[["elapsed"]] - t0, 3), " sec")
  }

  process_names <- vapply(object$backend$process_terms, `[[`, character(1), "name")
  prediction_w_samples <- if(isTRUE(return_w_samples)) {
    out <- vector("list", length(process_names))
    names(out) <- process_names
    out
  } else {
    NULL
  }

  for(i in seq_along(object$backend$process_terms)){
    term <- object$backend$process_terms[[i]]
    term_name <- process_names[i]
    stlmm_progress(
      verbose,
      "predict: processing ", term_name, " (", i, " of ",
      length(object$backend$process_terms), ")"
    )
    w_i <- w_samples_ordered[[term_name]]
    if(is.null(w_i))
      stop("error: saved or recovered latent process samples missing for ", term_name)

    map <- process_maps[[i]]$map
    process_scale <- process_maps[[i]]$scale
    if(predict_timing)
      t0_term <- proc.time()[["elapsed"]]

    term_w_samples <- if(isTRUE(return_w_samples)) {
      out <- matrix(0.0, nrow = n_draw, ncol = n0)
      colnames(out) <- row_names
      out
    } else {
      NULL
    }

    existing_col <- which(!is.na(map))
    if(length(existing_col)){
      existing_w <- w_i[recover_row, map[existing_col], drop = FALSE]
      if(isTRUE(return_w_samples))
        term_w_samples[, existing_col] <- existing_w
      existing_mu <- existing_w
      if(!is.null(process_scale))
        existing_mu <- sweep(existing_mu, 2L, process_scale[existing_col], `*`)
      mu[, existing_col] <- mu[, existing_col, drop = FALSE] + existing_mu
    }

    new_col <- which(is.na(map))
    if(length(new_col)){
      if(term$term_type == "ar1"){
        new_w <- simulate_ar1_prediction_nodes(
          object = object,
          term = term,
          graph = object$backend$graphs[[term$graph_index]],
          w_samples_ordered = w_samples_ordered,
          recover_row = recover_row,
          draw_index = draw_index,
          new_values = process_maps[[i]]$new_values
        )
      } else if(term$term_type == "gp"){
        new_w <- simulate_gp_prediction_nodes(
          object = object,
          term = term,
          graph = object$backend$graphs[[term$graph_index]],
          w_samples_ordered = w_samples_ordered,
          recover_row = recover_row,
          draw_index = draw_index,
          new_coords = process_maps[[i]]$new_coords,
          joint = isTRUE(joint)
        )
      } else if(term$term_type == "nngp"){
        if(predict_timing)
          t0 <- proc.time()[["elapsed"]]
        if(identical(joint_method, "full")){
          message(
            "joint = TRUE for NNGP prediction builds a dense covariance over unique new prediction nodes; ",
            "memory use can be large, and joint = FALSE is recommended for large prediction sets."
          )
          new_w <- simulate_nngp_prediction_nodes_joint(
            object = object,
            term = term,
            graph = object$backend$graphs[[term$graph_index]],
            w_samples_ordered = w_samples_ordered,
            recover_row = recover_row,
            draw_index = draw_index,
            new_coords = process_maps[[i]]$new_coords,
            neighbor_index = process_maps[[i]]$neighbor_index,
            n_omp_threads = n_omp_threads
          )
        } else if(identical(joint_method, "vecchia")){
          stlmm_progress(verbose, "predict: simulating ", term_name, " Vecchia joint NNGP new nodes")
          new_w_ordered <- simulate_nngp_prediction_nodes_vecchia(
            object = object,
            term = term,
            graph = object$backend$graphs[[term$graph_index]],
            w_samples_ordered = w_samples_ordered,
            recover_row = recover_row,
            draw_index = draw_index,
            coords_all = process_maps[[i]]$vecchia_coords_all,
            neighbor_index = process_maps[[i]]$vecchia_neighbor_index,
            neighbor_count = process_maps[[i]]$vecchia_neighbor_count,
            n_omp_threads = n_omp_threads
          )
          new_w <- matrix(0.0, nrow = nrow(new_w_ordered), ncol = ncol(new_w_ordered))
          new_w[, process_maps[[i]]$vecchia_ord] <- new_w_ordered
        } else {
          stlmm_progress(verbose, "predict: simulating ", term_name, " independent NNGP new nodes")
          new_w <- simulate_nngp_prediction_nodes(
            object = object,
            term = term,
            graph = object$backend$graphs[[term$graph_index]],
            w_samples_ordered = w_samples_ordered,
            recover_row = recover_row,
            draw_index = draw_index,
            new_coords = process_maps[[i]]$new_coords,
            neighbor_index = process_maps[[i]]$neighbor_index,
            n_omp_threads = n_omp_threads
          )
        }
        if(predict_timing)
          message("predict timing: ", term_name, " NNGP new-node simulation = ",
                  round(proc.time()[["elapsed"]] - t0, 3), " sec")
      } else {
        stop("error: new-node prediction is not implemented for ", term$term_type)
      }

      new_w_rows <- new_w[, process_maps[[i]]$new_id[new_col], drop = FALSE]
      if(isTRUE(return_w_samples))
        term_w_samples[, new_col] <- new_w_rows
      new_mu <- new_w_rows
      if(!is.null(process_scale))
        new_mu <- sweep(new_mu, 2L, process_scale[new_col], `*`)
      mu[, new_col] <- mu[, new_col, drop = FALSE] + new_mu
    }

    if(isTRUE(return_w_samples))
      prediction_w_samples[[term_name]] <- term_w_samples

    if(predict_timing)
      message("predict timing: ", term_name, " process assembly total = ",
              round(proc.time()[["elapsed"]] - t0_term, 3), " sec")
  }

  if(!is.null(offset0)){
    if(length(offset0) != n0)
      stop("error: prediction offset length does not match prediction rows")
    mu <- sweep(mu, 2L, as.numeric(offset0), `+`)
  }

  trials <- prediction_trials(object, newdata, n0)
  gaussian_y_sampler <- function(eta){
    stlmm_progress(verbose, "predict: simulating observation-level predictive samples")
    if(predict_timing)
      t0 <- proc.time()[["elapsed"]]
    residual_sd <- prediction_residual_sd_samples(
      object = object,
      newdata = newdata,
      n0 = n0,
      draw_index = draw_index,
      where = "prediction rows"
    )
    y <- eta
    if(is.null(residual_sd)){
      tau <- sqrt(object$tau_sq_samples[draw_index])
      for(i in seq_len(n_draw))
        y[i, ] <- stats::rnorm(n0, mean = eta[i, ], sd = tau[i])
    } else {
      for(i in seq_len(n_draw))
        y[i, ] <- stats::rnorm(n0, mean = eta[i, ], sd = residual_sd[i, ])
    }
    if(predict_timing)
      message("predict timing: y sample simulation = ", round(proc.time()[["elapsed"]] - t0, 3), " sec")
    y
  }

  out <- finish_prediction_samples(
    object = object,
    eta = mu,
    y_samples = y_samples,
    scale = scale,
    trials = trials,
    draw_index = draw_index,
    newdata = newdata,
    joint = isTRUE(joint),
    joint_method = if(isTRUE(joint)) joint_method else NULL,
    w_samples = prediction_w_samples,
    gaussian_y_sampler = gaussian_y_sampler
  )
  stlmm_progress(verbose, "predict: done")
  out
}

predict.stLMM_recovery_chains <- function(object,
                                         newdata = NULL,
                                         y_samples = FALSE,
                                         joint = FALSE,
                                         joint_method = c("full", "vecchia"),
                                         pred_m = NULL,
                                         pred_ordering = "maxmin",
                                         st_scale = NULL,
                                         return_w_samples = TRUE,
                                         n_omp_threads = NULL,
                                         verbose = FALSE,
                                         sub_sample = list(start = 1L, thin = 1L),
                                         scale = c("response", "link"),
                                         ...){

  scale <- match.arg(scale)
  joint_method <- match.arg(joint_method)
  pred <- vector("list", length(object$chains))
  for(i in seq_along(object$chains)){
    stlmm_progress(verbose, "predict: chain ", i, " of ", length(object$chains))
    pred[[i]] <- predict(
      object$chains[[i]],
      newdata = newdata,
      y_samples = y_samples,
      joint = joint,
      joint_method = joint_method,
      pred_m = pred_m,
      pred_ordering = pred_ordering,
      st_scale = st_scale,
      return_w_samples = return_w_samples,
      n_omp_threads = n_omp_threads,
      verbose = verbose,
      sub_sample = sub_sample,
      scale = scale,
      ...
    )
  }
  out <- list(
    chains = pred,
    n_chains = object$n_chains,
    call = match.call(),
    recover_call = object$call,
    fit_call = object$fit_call
  )
  class(out) <- "stLMM_prediction_chains"
  out
}

build_existing_support_prediction_backend <- function(object,
                                                      newdata,
                                                      st_scale = NULL,
                                                      joint_method = c("none", "full", "vecchia"),
                                                      pred_m = NULL,
                                                      pred_ordering = "maxmin",
                                                      n_omp_threads = NULL,
                                                      verbose = FALSE){

  predict_timing <- identical(Sys.getenv("STLMM_PREDICT_TIMING"), "1")
  joint_method <- match.arg(joint_method)
  n_omp_threads <- resolve_n_omp_threads(object, n_omp_threads)
  if(predict_timing)
    t0 <- proc.time()[["elapsed"]]

  if(!is.data.frame(newdata))
    newdata <- as.data.frame(newdata)

  dat <- newdata
  response_name <- object$backend$response_name
  if(!is.null(response_name) && nzchar(response_name) && !response_name %in% names(dat))
    dat[[response_name]] <- 0

  re_terms <- object$term_description$random_effects$terms
  if(length(re_terms)){
    for(term in re_terms){
      g <- term$grouping_factor
      lev <- term$levels
      if(!g %in% names(dat))
        stop("error: grouping factor ", g, " missing from newdata")
      group_values <- as.character(dat[[g]])
      if(anyNA(group_values) || any(!group_values %in% lev))
        stop("error: new grouping levels are not currently supported for grouping factor ", g)
      dat[[g]] <- factor(group_values, levels = lev)
    }
  }
  if(predict_timing)
    message("predict timing: backend data prep = ", round(proc.time()[["elapsed"]] - t0, 3), " sec")

  fixed_formula <- object$backend$reduced_formula
  if(predict_timing)
    t0 <- proc.time()[["elapsed"]]
  mf <- stats::model.frame(fixed_formula, dat, na.action = stats::na.fail)
  X0 <- stats::model.matrix(fixed_formula, mf)
  offset0 <- stats::model.offset(mf)
  if(is.null(offset0))
    offset0 <- rep(0, nrow(X0))
  offset0 <- as.numeric(offset0)
  if(length(offset0) != nrow(X0) || anyNA(offset0) || any(!is.finite(offset0)))
    stop("error: prediction offset values must be finite numeric values with length matching newdata")
  if(predict_timing)
    message("predict timing: fixed-effect prediction design = ", round(proc.time()[["elapsed"]] - t0, 3), " sec")

  if(predict_timing)
    t0 <- proc.time()[["elapsed"]]
  Z0 <- build_prediction_Z(dat, re_terms, object$backend$Z)
  if(predict_timing)
    message("predict timing: random-effect prediction design = ", round(proc.time()[["elapsed"]] - t0, 3), " sec")

  if(!identical(colnames(X0), colnames(object$backend$X)))
    stop("error: newdata fixed-effect columns do not match fitted model")
  if(ncol(Z0) != ncol(object$backend$Z))
    stop("error: newdata random-effect columns do not match fitted model")

  process_maps <- vector("list", length(object$backend$process_terms))

  for(i in seq_along(object$backend$process_terms)){
    term <- object$backend$process_terms[[i]]
    graph <- object$backend$graphs[[term$graph_index]]
    if(predict_timing)
      t0 <- proc.time()[["elapsed"]]

    if(term$term_type == "ar1"){
      spec <- parse_process_call(term$label)
      time_name <- spec$args[1]
      if(!time_name %in% names(dat))
        stop("error: AR1 time variable ", time_name, " missing from newdata")
      map <- match(dat[[time_name]], graph$support)
      new_values <- sort(unique(dat[[time_name]][is.na(map)]))
      new_id <- match(dat[[time_name]], new_values)
    } else if(term$term_type %in% c("car", "dagar")){
      spec <- parse_process_call(term$label)
      area_name <- spec$args[1]
      if(!area_name %in% names(dat))
        stop("error: ", toupper(term$term_type), " area variable ", area_name, " missing from newdata")
      ids <- if(term$term_type == "dagar") graph$ids_support %||% graph$ids else graph$ids
      support_map <- match(as.character(dat[[area_name]]), as.character(ids))
      map <- if(term$term_type == "dagar") graph$ord_inv[support_map] else support_map
      if(anyNA(map)){
        bad <- unique(as.character(dat[[area_name]])[is.na(map)])
        stop("error: new ", toupper(term$term_type), " area value(s) not found in fitted graph: ", paste(bad, collapse = ", "))
      }
      new_values <- character(0)
      new_id <- rep(NA_integer_, nrow(dat))
    } else if(term$term_type %in% c("car_time", "dagar_time")){
      spec <- parse_process_call(term$label)
      area_name <- spec$args[1]
      time_name <- spec$args[2]
      if(!area_name %in% names(dat))
        stop("error: ", if(term$term_type == "dagar_time") "DAGAR-time" else "CAR-time", " area variable ", area_name, " missing from newdata")
      if(!time_name %in% names(dat))
        stop("error: ", if(term$term_type == "dagar_time") "DAGAR-time" else "CAR-time", " time variable ", time_name, " missing from newdata")
      ids <- if(term$term_type == "dagar_time") graph$ids_support %||% graph$ids else graph$ids
      area_map_support <- match(as.character(dat[[area_name]]), as.character(ids))
      area_map <- if(term$term_type == "dagar_time") graph$ord_inv[area_map_support] else area_map_support
      if(anyNA(area_map)){
        bad <- unique(as.character(dat[[area_name]])[is.na(area_map)])
        stop("error: new ", if(term$term_type == "dagar_time") "DAGAR-time" else "CAR-time", " area value(s) not found in fitted graph: ", paste(bad, collapse = ", "))
      }
      time_map <- match(dat[[time_name]], graph$time_support)
      if(anyNA(time_map)){
        bad <- unique(as.character(dat[[time_name]])[is.na(time_map)])
        stop("error: new ", if(term$term_type == "dagar_time") "DAGAR-time" else "CAR-time", " value(s) not found in fitted time support: ", paste(bad, collapse = ", "))
      }
      map <- as.integer((area_map - 1L) * graph$n_time + time_map)
      new_values <- character(0)
      new_id <- rep(NA_integer_, nrow(dat))
    } else if(term$term_type %in% c("nngp", "gp")){
      coord_names <- graph$coord_names
      miss <- setdiff(coord_names, names(dat))
      if(length(miss))
        stop("error: missing process coordinate column(s): ", paste(miss, collapse = ", "))
      coords0 <- as.matrix(dat[, coord_names, drop = FALSE])
      storage.mode(coords0) <- "double"
      support <- graph$coords_ord
      map <- match_unique_rows(coords0, support)
      if(anyNA(map)){
        new_coords <- sort_unique_matrix_rows(coords0[is.na(map), , drop = FALSE])
        new_id <- match_unique_rows(coords0, new_coords)
      } else {
        new_coords <- coords0[FALSE, , drop = FALSE]
        new_id <- rep(NA_integer_, nrow(coords0))
      }

      neighbor_index <- NULL
      st_scale_i <- NULL
      if(term$term_type == "nngp" && nrow(new_coords) > 0L){
        st_scale_i <- resolve_prediction_st_scale(
          st_scale = st_scale,
          term_name = term$name,
          graph = graph,
          cov_model = term$cov_model
        )
        stlmm_progress(verbose, "predict: finding ", term$name, " fitted-support neighbors")
        if(predict_timing)
          t0_neighbor <- proc.time()[["elapsed"]]
        neighbor_index <- nngp_prediction_neighbors(
          new_coords = new_coords,
          support = graph$coords_ord,
          m = graph$m,
          cov_model = term$cov_model,
          st_scale = st_scale_i,
          term_name = term$name,
          n_omp_threads = n_omp_threads
        )
        if(predict_timing)
          message("predict timing: ", term$name, " NNGP neighbor search = ",
                  round(proc.time()[["elapsed"]] - t0_neighbor, 3), " sec")
      }
      vecchia_graph <- NULL
      if(term$term_type == "nngp" && nrow(new_coords) > 0L && identical(joint_method, "vecchia")){
        stlmm_progress(verbose, "predict: building ", term$name, " Vecchia prediction graph")
        if(predict_timing)
          t0_vecchia <- proc.time()[["elapsed"]]
        vecchia_graph <- nngp_prediction_vecchia_graph(
          new_coords = new_coords,
          support = graph$coords_ord,
          m = pred_m %||% graph$m,
          ordering = pred_ordering,
          cov_model = term$cov_model,
          st_scale = st_scale_i,
          term_name = term$name,
          n_omp_threads = n_omp_threads
        )
        if(predict_timing)
          message("predict timing: ", term$name, " Vecchia NNGP prediction graph = ",
                  round(proc.time()[["elapsed"]] - t0_vecchia, 3), " sec")
      }
    } else {
      stop("error: unsupported process term type ", term$term_type)
    }

    scale <- NULL
    if(!is.null(term$coef_name)){
      if(!term$coef_name %in% names(dat))
        stop("error: SVC covariate ", term$coef_name, " missing from newdata")
      scale <- as.numeric(dat[[term$coef_name]])
      if(length(scale) != nrow(dat) || anyNA(scale))
        stop("error: SVC covariate ", term$coef_name, " in newdata must be numeric with no missing values")
    }

    if(term$term_type == "ar1" || term$term_type == "car" || term$term_type == "dagar" || term$term_type == "car_time" || term$term_type == "dagar_time"){
      process_maps[[i]] <- list(
        map = as.integer(map),
        scale = scale,
        new_values = new_values,
        new_id = as.integer(new_id)
      )
    } else {
      process_maps[[i]] <- list(
        map = as.integer(map),
        scale = scale,
        new_coords = new_coords,
        new_id = as.integer(new_id),
        neighbor_index = neighbor_index,
        vecchia_coords_all = if(is.null(vecchia_graph)) NULL else vecchia_graph$coords_all,
        vecchia_neighbor_index = if(is.null(vecchia_graph)) NULL else vecchia_graph$neighbor_index,
        vecchia_neighbor_count = if(is.null(vecchia_graph)) NULL else vecchia_graph$neighbor_count,
        vecchia_ord = if(is.null(vecchia_graph)) NULL else vecchia_graph$ord
      )
    }

    if(predict_timing)
      message("predict timing: ", term$name, " process prediction map = ",
              round(proc.time()[["elapsed"]] - t0, 3), " sec")
  }

  list(X = X0, Z = Z0, offset = offset0, process_maps = process_maps)
}

nngp_prediction_neighbors <- function(new_coords,
                                      support,
                                      m,
                                      cov_model,
                                      st_scale,
                                      term_name,
                                      n_omp_threads = 1L){

  st_scale_i <- resolve_st_scale(st_scale, term_name)
  if(!is_space_time_cov_model(cov_model) && !isTRUE(all.equal(st_scale_i, 1)))
    stop("error: st_scale is only valid for space-time NNGP covariance models")

  k <- min(as.integer(m), nrow(support))

  .Call(
    "stLMM_nngp_prediction_neighbors",
    support,
    new_coords,
    as.integer(k),
    cov_model,
    st_scale_i,
    as.integer(n_omp_threads),
    PACKAGE = "stLMM"
  )
}

nngp_prediction_vecchia_graph <- function(new_coords,
                                          support,
                                          m,
                                          ordering,
                                          cov_model,
                                          st_scale,
                                          term_name,
                                          n_omp_threads = 1L){

  new_coords <- as.matrix(new_coords)
  support <- as.matrix(support)
  storage.mode(new_coords) <- "double"
  storage.mode(support) <- "double"

  if(nrow(new_coords) < 1L)
    stop("error: Vecchia NNGP prediction graph requires at least one new node")
  if(nrow(support) < 1L)
    stop("error: Vecchia NNGP prediction graph requires fitted support nodes")
  if(ncol(new_coords) != ncol(support))
    stop("error: prediction and fitted support coordinates have different dimensions")

  space_time <- is_space_time_cov_model(cov_model)
  st_scale_i <- resolve_st_scale(st_scale, term_name)
  if(!space_time && !isTRUE(all.equal(st_scale_i, 1)))
    stop("error: st_scale is only valid for space-time NNGP covariance models")
  if(space_time && identical(ordering, "hilbert"))
    stop("error: hilbert prediction ordering is not supported for space-time NNGP covariance models")

  new_search <- scale_nngp_search_coords(
    new_coords,
    st_scale = st_scale_i,
    space_time = space_time
  )
  support_search <- scale_nngp_search_coords(
    support,
    st_scale = st_scale_i,
    space_time = space_time
  )

  ord_info <- compute_nngp_order(new_search, ordering)
  new_ord <- new_coords[ord_info$ord, , drop = FALSE]
  if(space_time){
    new_search_ord <- new_search[ord_info$ord, , drop = FALSE]
    search_coords <- rbind(support_search, new_search_ord)
  } else {
    search_coords <- rbind(support, new_ord)
  }
  coords_all <- rbind(support, new_ord)

  n_all <- nrow(search_coords)
  m <- min(as.integer(m), n_all - 1L)
  if(m < 1L)
    stop("error: pred_m must leave at least one available neighbor")

  nn <- mkNNIndx(search_coords, m = m, n_omp_threads = n_omp_threads)
  nn_indx <- nn$nnIndx
  nn_lu <- nn$nnIndxLU
  n_fit <- nrow(support)
  n_pred <- nrow(new_ord)
  neighbor_index <- matrix(0L, nrow = n_pred, ncol = m)
  neighbor_count <- integer(n_pred)

  for(j in seq_len(n_pred)){
    row <- n_fit + j
    start <- nn_lu[row] + 1L
    count <- nn_lu[n_all + row]
    neighbor_count[j] <- count
    if(count > 0L)
      neighbor_index[j, seq_len(count)] <- nn_indx[start + seq_len(count) - 1L] + 1L
  }

  list(
    coords_all = coords_all,
    neighbor_index = neighbor_index,
    neighbor_count = as.integer(neighbor_count),
    ord = as.integer(ord_info$ord),
    ordering = ord_info$ordering_type,
    m = m
  )
}

sort_unique_matrix_rows <- function(x){

  x <- as.matrix(x)
  if(nrow(x) == 0L)
    return(x)

  x <- unique(x)
  ord <- do.call(order, as.data.frame(x))
  x[ord, , drop = FALSE]
}

is_space_time_cov_model <- function(cov_model){
  registry <- build_cor_model_registry()
  model_info <- registry[[cov_model]]
  !is.null(model_info) && identical(as.integer(model_info$distance_mode), 2L)
}

resolve_prediction_st_scale <- function(st_scale, term_name, graph, cov_model){

  space_time <- is_space_time_cov_model(cov_model)
  fitted_scale <- graph$st_scale %||% 1

  if(is.null(st_scale)){
    out <- fitted_scale
  } else if(is.list(st_scale)){
    if(!term_name %in% names(st_scale)){
      if(!space_time)
        return(1)
      stop("error: st_scale list is missing an entry for ", term_name)
    }
    out <- st_scale[[term_name]]
  } else if(length(st_scale) > 1L){
    if(is.null(names(st_scale)) || !term_name %in% names(st_scale)){
      if(!space_time)
        return(1)
      stop("error: named st_scale vector must include ", term_name)
    }
    out <- st_scale[[term_name]]
  } else {
    out <- st_scale
  }

  out <- resolve_st_scale(out, term_name)
  if(!space_time && !isTRUE(all.equal(out, 1)))
    stop("error: st_scale is only valid for space-time NNGP covariance models")

  out
}

resolve_st_scale <- function(st_scale, term_name){

  if(is.list(st_scale)){
    if(!term_name %in% names(st_scale))
      stop("error: st_scale list is missing an entry for ", term_name)
    st_scale <- st_scale[[term_name]]
  } else if(length(st_scale) > 1L){
    if(is.null(names(st_scale)) || !term_name %in% names(st_scale))
      stop("error: named st_scale vector must include ", term_name)
    st_scale <- st_scale[[term_name]]
  }

  if(!is.numeric(st_scale) || length(st_scale) != 1L || is.na(st_scale) || st_scale <= 0)
    stop("error: st_scale must be a positive scalar")

  as.numeric(st_scale)
}

build_prediction_Z <- function(dat, re_terms, Z_fit){

  if(!length(re_terms))
    return(Matrix::Matrix(0, nrow(dat), 0, sparse = TRUE))

  n <- nrow(dat)
  q <- ncol(Z_fit)
  n_value <- sum(vapply(re_terms, function(term) n * length(term$coefficients), integer(1)))
  ii <- integer(n_value)
  jj <- integer(n_value)
  xx <- numeric(n_value)
  pos <- 0L
  col_offset <- 0L

  for(term in re_terms){
    group_values <- as.character(dat[[term$grouping_factor]])
    lev <- term$levels
    coef_names <- term$coefficients
    n_coef <- length(coef_names)

    for(i in seq_len(n)){
      level_index <- match(group_values[i], lev)

      for(j in seq_len(n_coef)){
        pos <- pos + 1L
        ii[pos] <- i
        jj[pos] <- col_offset + (level_index - 1L) * n_coef + j
        if(coef_names[j] == "(Intercept)"){
          xx[pos] <- 1.0
        } else {
          if(!coef_names[j] %in% names(dat))
            stop("error: iid slope covariate ", coef_names[j], " missing from newdata")
          if(!is.numeric(dat[[coef_names[j]]]))
            stop("error: iid slope covariate ", coef_names[j], " in newdata must be numeric")
          value <- as.numeric(dat[[coef_names[j]]][i])
          if(anyNA(value))
            stop("error: iid slope covariate ", coef_names[j], " in newdata contains missing values")
          xx[pos] <- value
        }
      }
    }

    col_offset <- col_offset + length(lev) * n_coef
  }

  if(col_offset != q)
    stop("error: newdata random-effect design does not match fitted model")

  Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(n, q))
}

simulate_ar1_prediction_nodes <- function(object,
                                          term,
                                          graph,
                                          recover_row,
                                          draw_index,
                                          new_values,
                                          w_samples_ordered = NULL){

  n_new <- length(new_values)
  n_draw <- length(draw_index)
  out <- matrix(0.0, nrow = n_draw, ncol = n_new)

  if(n_new == 0L)
    return(out)

  sigma_sq <- process_sigma_sq_draws(object, term, draw_index)
  phi <- object$theta_samples[draw_index, paste0(term$name, "_phi")]
  support <- graph$support
  if(is.null(w_samples_ordered))
    w_samples_ordered <- object$w_samples_ordered %||% object$w_samples
  w_fit <- w_samples_ordered[[term$name]][recover_row, , drop = FALSE]

  for(j in seq_len(n_new)){
    value <- new_values[j]
    right <- which(support > value)[1]
    left <- right - 1L

    for(i in seq_len(n_draw)){
      if(is.na(right)){
        mean_j <- phi[i] * w_fit[i, length(support)]
        var_j <- sigma_sq[i] * (1 - phi[i]^2)
      } else if(left < 1L){
        mean_j <- phi[i] * w_fit[i, 1L]
        var_j <- sigma_sq[i] * (1 - phi[i]^2)
      } else {
        denom <- 1 + phi[i]^2
        mean_j <- phi[i] * (w_fit[i, left] + w_fit[i, right]) / denom
        var_j <- sigma_sq[i] * (1 - phi[i]^2) / denom
      }

      out[i, j] <- stats::rnorm(1L, mean = mean_j, sd = sqrt(max(var_j, 0)))
    }
  }

  out
}

simulate_gp_prediction_nodes <- function(object,
                                         term,
                                         graph,
                                         recover_row,
                                         draw_index,
                                         new_coords,
                                         w_samples_ordered = NULL,
                                         joint = FALSE){

  n_new <- nrow(new_coords)
  n_draw <- length(draw_index)
  out <- matrix(0.0, nrow = n_draw, ncol = n_new)

  if(n_new == 0L)
    return(out)

  support <- graph$coords_ord
  if(is.null(w_samples_ordered))
    w_samples_ordered <- object$w_samples_ordered %||% object$w_samples
  w_fit <- w_samples_ordered[[term$name]][recover_row, , drop = FALSE]
  sigma_sq <- process_sigma_sq_draws(object, term, draw_index)
  theta_cols <- paste0(term$name, "_", term$theta_names)

  for(i in seq_len(n_draw)){
    theta <- as.numeric(object$theta_samples[draw_index[i], theta_cols, drop = FALSE])
    names(theta) <- term$theta_names
    C <- gp_covariance_matrix(support, support, term$cov_model, sigma_sq[i], theta)
    C <- C + diag(1e-10, nrow(C))
    C_chol <- chol(C)

    C_no <- gp_covariance_matrix(new_coords, support, term$cov_model, sigma_sq[i], theta)

    if(isTRUE(joint)){
      C_nn <- gp_covariance_matrix(new_coords, new_coords, term$cov_model, sigma_sq[i], theta)
      Coo_inv_Con <- backsolve(C_chol, forwardsolve(t(C_chol), t(C_no)))
      mean_i <- as.numeric(C_no %*% backsolve(C_chol, forwardsolve(t(C_chol), w_fit[i, ])))
      var_i <- C_nn - C_no %*% Coo_inv_Con
      var_i <- 0.5 * (var_i + t(var_i))
      diag(var_i) <- pmax(diag(var_i), 0)
      chol_i <- tryCatch(chol(var_i + diag(1e-10, n_new)), error = identity)
      if(inherits(chol_i, "error"))
        stop("error: joint dense GP prediction covariance is not positive definite")
      out[i, ] <- mean_i + as.numeric(stats::rnorm(n_new) %*% chol_i)
    } else {
      for(j in seq_len(n_new)){
        c0 <- C_no[j, , drop = FALSE]
        v <- backsolve(C_chol, forwardsolve(t(C_chol), as.numeric(c0)))
        mean_j <- sum(v * w_fit[i, ])
        var_j <- sigma_sq[i] - sum(as.numeric(c0) * v)
        out[i, j] <- stats::rnorm(1L, mean = mean_j, sd = sqrt(max(var_j, 0)))
      }
    }
  }

  out
}

simulate_nngp_prediction_nodes <- function(object,
                                           term,
                                           graph,
                                           recover_row,
                                           draw_index,
                                           new_coords,
                                           neighbor_index,
                                           w_samples_ordered = NULL,
                                           n_omp_threads = NULL){

  n_new <- nrow(new_coords)
  n_draw <- length(draw_index)

  if(n_new == 0L)
    return(matrix(0.0, nrow = n_draw, ncol = n_new))
  n_omp_threads <- resolve_n_omp_threads(object, n_omp_threads)

  support <- graph$coords_ord
  if(is.null(w_samples_ordered))
    w_samples_ordered <- object$w_samples_ordered %||% object$w_samples
  w_fit <- w_samples_ordered[[term$name]][recover_row, , drop = FALSE]
  sigma_sq <- process_sigma_sq_draws(object, term, draw_index)
  theta_cols <- paste0(term$name, "_", term$theta_names)
  theta <- object$theta_samples[draw_index, theta_cols, drop = FALSE]

  storage.mode(support) <- "double"
  storage.mode(new_coords) <- "double"
  storage.mode(w_fit) <- "double"
  storage.mode(sigma_sq) <- "double"
  storage.mode(theta) <- "double"
  storage.mode(neighbor_index) <- "integer"

  .Call(
    "stLMM_predict_nngp_joint_false",
    support,
    new_coords,
    neighbor_index,
    w_fit,
    as.numeric(sigma_sq),
    theta,
    term$cov_model,
    as.integer(n_omp_threads),
    PACKAGE = "stLMM"
  )
}

simulate_nngp_prediction_nodes_joint <- function(object,
                                                 term,
                                                 graph,
                                                 recover_row,
                                                 draw_index,
                                                 new_coords,
                                                 neighbor_index,
                                                 w_samples_ordered = NULL,
                                                 n_omp_threads = NULL){

  n_new <- nrow(new_coords)
  n_draw <- length(draw_index)

  if(n_new == 0L)
    return(matrix(0.0, nrow = n_draw, ncol = n_new))
  invisible(resolve_n_omp_threads(object, n_omp_threads))

  support <- graph$coords_ord
  if(is.null(w_samples_ordered))
    w_samples_ordered <- object$w_samples_ordered %||% object$w_samples
  w_fit <- w_samples_ordered[[term$name]][recover_row, , drop = FALSE]
  sigma_sq <- process_sigma_sq_draws(object, term, draw_index)
  theta_cols <- paste0(term$name, "_", term$theta_names)
  theta_draws <- object$theta_samples[draw_index, theta_cols, drop = FALSE]

  support <- as.matrix(support)
  new_coords <- as.matrix(new_coords)
  storage.mode(support) <- "double"
  storage.mode(new_coords) <- "double"
  storage.mode(neighbor_index) <- "integer"

  m <- ncol(neighbor_index)
  out <- matrix(0.0, nrow = n_draw, ncol = n_new)

  for(draw in seq_len(n_draw)){
    theta <- as.numeric(theta_draws[draw, , drop = TRUE])
    names(theta) <- term$theta_names

    B <- matrix(0.0, nrow = n_new, ncol = m)
    mean_draw <- numeric(n_new)

    for(node in seq_len(n_new)){
      nbr <- neighbor_index[node, ]
      C_nn <- gp_correlation_matrix(support[nbr, , drop = FALSE],
                                    support[nbr, , drop = FALSE],
                                    term$cov_model,
                                    theta)
      C_nx <- gp_correlation_matrix(support[nbr, , drop = FALSE],
                                    new_coords[node, , drop = FALSE],
                                    term$cov_model,
                                    theta)
      C_nn <- C_nn + diag(1e-10, nrow(C_nn))
      B[node, ] <- as.numeric(solve(C_nn, C_nx))
      mean_draw[node] <- sum(B[node, ] * w_fit[draw, nbr])
    }

    S <- matrix(0.0, nrow = n_new, ncol = n_new)
    for(a in seq_len(n_new)){
      nbr_a <- neighbor_index[a, ]
      for(b in seq_len(a)){
        nbr_b <- neighbor_index[b, ]
        c_ab <- gp_correlation_matrix(new_coords[a, , drop = FALSE],
                                      new_coords[b, , drop = FALSE],
                                      term$cov_model,
                                      theta)[1, 1]
        c_Na_b <- gp_correlation_matrix(support[nbr_a, , drop = FALSE],
                                        new_coords[b, , drop = FALSE],
                                        term$cov_model,
                                        theta)
        c_a_Nb <- gp_correlation_matrix(new_coords[a, , drop = FALSE],
                                        support[nbr_b, , drop = FALSE],
                                        term$cov_model,
                                        theta)
        c_Na_Nb <- gp_correlation_matrix(support[nbr_a, , drop = FALSE],
                                         support[nbr_b, , drop = FALSE],
                                         term$cov_model,
                                         theta)

        s_ab <- c_ab -
          sum(B[a, ] * as.numeric(c_Na_b)) -
          sum(B[b, ] * as.numeric(c_a_Nb)) +
          as.numeric(B[a, , drop = FALSE] %*% c_Na_Nb %*% B[b, ])

        S[a, b] <- sigma_sq[draw] * s_ab
        S[b, a] <- S[a, b]
      }
    }

    S <- (S + t(S)) / 2
    diag(S) <- pmax(diag(S), 0)
    chol_S <- try(chol(S + diag(1e-10, n_new)), silent = TRUE)
    if(inherits(chol_S, "try-error"))
      stop("error: joint NNGP prediction covariance is not positive definite")

    out[draw, ] <- mean_draw + as.numeric(t(chol_S) %*% stats::rnorm(n_new))
  }

  out
}

simulate_nngp_prediction_nodes_vecchia <- function(object,
                                                   term,
                                                   graph,
                                                   recover_row,
                                                   draw_index,
                                                   coords_all,
                                                   neighbor_index,
                                                   neighbor_count,
                                                   w_samples_ordered = NULL,
                                                   n_omp_threads = NULL){

  n_pred <- length(neighbor_count)
  n_draw <- length(draw_index)

  if(n_pred == 0L)
    return(matrix(0.0, nrow = n_draw, ncol = n_pred))
  n_omp_threads <- resolve_n_omp_threads(object, n_omp_threads)

  if(is.null(w_samples_ordered))
    w_samples_ordered <- object$w_samples_ordered %||% object$w_samples
  w_fit <- w_samples_ordered[[term$name]][recover_row, , drop = FALSE]
  sigma_sq <- process_sigma_sq_draws(object, term, draw_index)
  theta_cols <- paste0(term$name, "_", term$theta_names)
  theta <- object$theta_samples[draw_index, theta_cols, drop = FALSE]

  storage.mode(coords_all) <- "double"
  storage.mode(w_fit) <- "double"
  storage.mode(sigma_sq) <- "double"
  storage.mode(theta) <- "double"
  storage.mode(neighbor_index) <- "integer"
  storage.mode(neighbor_count) <- "integer"

  .Call(
    "stLMM_predict_nngp_vecchia_joint",
    coords_all,
    as.integer(nrow(graph$coords_ord)),
    neighbor_index,
    neighbor_count,
    w_fit,
    as.numeric(sigma_sq),
    theta,
    term$cov_model,
    as.integer(n_omp_threads),
    PACKAGE = "stLMM"
  )
}

gp_correlation_matrix <- function(coords_a, coords_b, cov_model, theta){

  coords_a <- as.matrix(coords_a)
  coords_b <- as.matrix(coords_b)
  n_a <- nrow(coords_a)
  n_b <- nrow(coords_b)
  out <- matrix(0.0, nrow = n_a, ncol = n_b)
  has_time <- cov_model %in% c("sep_exp", "multi_res_sep_exp", "gneiting")

  for(i in seq_len(n_a)){
    for(j in seq_len(n_b)){
      if(has_time){
        d <- ncol(coords_a)
        h <- sqrt(sum((coords_a[i, seq_len(d - 1L)] - coords_b[j, seq_len(d - 1L)])^2))
        u <- abs(coords_a[i, d] - coords_b[j, d])
      } else {
        h <- sqrt(sum((coords_a[i, ] - coords_b[j, ])^2))
        u <- 0
      }

      spatial_dim <- if(has_time) d - 1L else ncol(coords_a)
      out[i, j] <- gp_correlation(cov_model, theta, h, u, spatial_dim)
    }
  }

  out
}

gp_covariance_matrix <- function(coords_a, coords_b, cov_model, sigma_sq, theta){

  coords_a <- as.matrix(coords_a)
  coords_b <- as.matrix(coords_b)
  n_a <- nrow(coords_a)
  n_b <- nrow(coords_b)
  out <- matrix(0.0, nrow = n_a, ncol = n_b)
  has_time <- cov_model %in% c("sep_exp", "multi_res_sep_exp", "gneiting")

  for(i in seq_len(n_a)){
    for(j in seq_len(n_b)){
      if(has_time){
        d <- ncol(coords_a)
        h <- sqrt(sum((coords_a[i, seq_len(d - 1L)] - coords_b[j, seq_len(d - 1L)])^2))
        u <- abs(coords_a[i, d] - coords_b[j, d])
      } else {
        h <- sqrt(sum((coords_a[i, ] - coords_b[j, ])^2))
        u <- 0
      }

      spatial_dim <- if(has_time) d - 1L else ncol(coords_a)
      out[i, j] <- sigma_sq * gp_correlation(cov_model, theta, h, u, spatial_dim)
    }
  }

  out
}

gp_correlation <- function(cov_model, theta, h, u, spatial_dim = 1L){
  if(cov_model == "exp")
    return(exp(-theta["phi"] * h))
  if(cov_model == "sep_exp")
    return(exp(-theta["phi"] * h) * exp(-theta["lambda"] * u))
  if(cov_model == "multi_res_sep_exp"){
    return(
      theta["alpha"] * exp(-theta["phi_1"] * h) * exp(-theta["lambda_1"] * u) +
        (1 - theta["alpha"]) * exp(-theta["phi_2"] * h) * exp(-theta["lambda_2"] * u)
    )
  }
  if(cov_model == "gneiting"){
    t <- 1 + theta["a"] * abs(u)^(2 * theta["alpha"])
    return(t^(-(theta["delta"] + spatial_dim / 2)) *
             exp(-theta["c"] * h^(2 * theta["gamma"]) /
                   t^(theta["beta"] * theta["gamma"])))
  }
  stop("error: unknown covariance model ", cov_model)
}
