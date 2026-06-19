posterior_summary_table <- function(x, probs = c(0.025, 0.5, 0.975)){

  if(is.null(x))
    return(NULL)

  if(is.null(dim(x))){
    x <- matrix(as.numeric(x), ncol = 1L)
    colnames(x) <- "value"
  } else {
    x <- as.matrix(x)
  }

  if(ncol(x) == 0L)
    return(NULL)

  if(is.null(colnames(x)))
    colnames(x) <- paste0("param_", seq_len(ncol(x)))

  out <- matrix(NA_real_, nrow = ncol(x), ncol = 2L + length(probs))
  rownames(out) <- colnames(x)
  colnames(out) <- c("mean", "sd", paste0("q", format(100 * probs, trim = TRUE, scientific = FALSE)))

  out[, "mean"] <- colMeans(x)
  out[, "sd"] <- apply(x, 2L, stats::sd)
  qs <- apply(x, 2L, stats::quantile, probs = probs, names = FALSE)
  if(length(probs) == 1L)
    qs <- matrix(qs, nrow = 1L)
  out[, 2L + seq_along(probs)] <- t(qs)

  out
}

sample_parameter_names <- function(samples, include_w = FALSE){

  if(!isTRUE(include_w))
    samples$w <- NULL

  active <- Filter(Negate(is.null), samples)
  out <- character(0)

  for(nm in names(active)){
    z <- active[[nm]]
    if(is.null(dim(z))){
      out <- c(out, nm)
    } else {
      z <- as.matrix(z)
      cn <- colnames(z)
      if(is.null(cn))
        cn <- paste0("param_", seq_len(ncol(z)))
      out <- c(out, cn)
    }
  }

  out
}

validate_summary_parameters <- function(parameters, available){

  if(is.null(parameters))
    return(NULL)

  if(!is.character(parameters) || anyNA(parameters) || any(!nzchar(parameters)))
    stop("error: parameters must be a non-empty character vector")

  miss <- setdiff(parameters, available)
  if(length(miss))
    stop("error: unknown parameter(s): ", paste(miss, collapse = ", "))

  parameters
}

filter_sample_blocks <- function(samples, parameters = NULL, include_w = FALSE){

  if(is.null(parameters))
    return(samples)

  parameters <- validate_summary_parameters(
    parameters,
    sample_parameter_names(samples, include_w = include_w)
  )

  out <- samples
  for(nm in names(out)){
    if(!isTRUE(include_w) && identical(nm, "w"))
      next

    z <- out[[nm]]
    if(is.null(z))
      next

    if(is.null(dim(z))){
      if(nm %in% parameters)
        out[[nm]] <- matrix(as.numeric(z), ncol = 1L,
                            dimnames = list(NULL, nm))
      else
        out[[nm]] <- NULL
      next
    }

    z <- as.matrix(z)
    cn <- colnames(z)
    if(is.null(cn)){
      cn <- paste0("param_", seq_len(ncol(z)))
      colnames(z) <- cn
    }

    keep <- intersect(parameters, cn)
    if(length(keep))
      out[[nm]] <- z[, keep, drop = FALSE]
    else
      out[[nm]] <- NULL
  }

  out
}

collect_parameter_summaries <- function(samples, probs = c(0.025, 0.5, 0.975)){

  out <- list()

  for(nm in names(samples)){
    if(nm == "w")
      next
    tab <- posterior_summary_table(samples[[nm]], probs = probs)
    if(!is.null(tab))
      out[[nm]] <- tab
  }

  out
}

format_model_formula <- function(x){
  if(is.null(x))
    return(NULL)

  if(is.character(x))
    txt <- paste(x, collapse = " ")
  else
    txt <- paste(deparse(x, width.cutoff = 500L), collapse = " ")

  gsub("[[:space:]]+", " ", txt)
}

print.stLMM <- function(x, ...){

  td <- x$term_description
  cat("stLMM fit\n")
  if(!is.null(td$global$formula))
    cat("  formula: ", format_model_formula(td$global$formula), "\n", sep = "")
  n_missing_response <- x$backend$n_missing_response %||% 0L
  if(n_missing_response > 0L)
    cat("  observations: ", x$backend$n, " (", x$backend$n_obs, " observed, ",
        n_missing_response, " missing response)\n", sep = "")
  else
    cat("  observations: ", x$backend$n, "\n", sep = "")
  family_name <- x$backend$family %||% "gaussian"
  cat("  family: ", family_name, "\n", sep = "")
  cat("  posterior draws: ", x$backend$n_samples, "\n", sep = "")
  cat("  fixed effects: ", x$backend$p, "\n", sep = "")
  cat("  grouped random-effect coefficients: ", x$backend$q, "\n", sep = "")
  cat("  process terms: ", length(x$backend$process_terms), "\n", sep = "")
  residual_type <- x$backend$residual_model$type %||% "global_tau"
  if(identical(family_name, "binomial") ||
     identical(family_name, "negative_binomial"))
    cat("  residual variance: not used for Polya-Gamma likelihood\n", sep = "")
  else if(identical(residual_type, "fixed_variance"))
    cat("  residual variance: fixed row-specific values\n", sep = "")
  else if(identical(residual_type, "group_ig_variance"))
    cat("  residual variance: sampled group-specific values\n", sep = "")
  else if(identical(residual_type, "scaled_variance"))
    cat("  residual variance: scaled direct-estimate values\n", sep = "")
  else
    cat("  residual variance: global tau_sq\n", sep = "")

  active <- names(Filter(Negate(is.null), x$samples))
  if(length(active))
    cat("  samples: ", paste(active, collapse = ", "), "\n", sep = "")

  if(length(x$backend$process_terms) && is.null(x$w_samples)){
    if(identical(family_name, "binomial") ||
       identical(family_name, "negative_binomial"))
      cat("  latent process samples: not retained; refit with save_process enabled\n")
    else
      cat("  latent process samples: not retained; call recover()\n")
  }
  else if(length(x$backend$process_terms) && !is.null(x$w_samples))
    cat("  latent process samples: retained (", length(x$recover_iter), " draws)\n", sep = "")

  sampler_time <- x$timing$sampler
  if(!is.null(sampler_time) && "elapsed" %in% names(sampler_time))
    cat("  sampler elapsed: ", round(unname(sampler_time[["elapsed"]]), 3), " seconds\n", sep = "")

  invisible(x)
}

print.stLMM_chains <- function(x, ...){

  cat("stLMM multi-chain fit\n")
  cat("  chains: ", x$n_chains, "\n", sep = "")
  if(length(x$chains)){
    td <- x$chains[[1L]]$term_description
    if(!is.null(td$global$formula))
      cat("  formula: ", format_model_formula(td$global$formula), "\n", sep = "")
    n_missing_response <- x$chains[[1L]]$backend$n_missing_response %||% 0L
    if(n_missing_response > 0L)
      cat("  observations: ", x$chains[[1L]]$backend$n, " (",
          x$chains[[1L]]$backend$n_obs, " observed, ",
          n_missing_response, " missing response)\n", sep = "")
    else
      cat("  observations: ", x$chains[[1L]]$backend$n, "\n", sep = "")
    family_name <- x$chains[[1L]]$backend$family %||% "gaussian"
    cat("  family: ", family_name, "\n", sep = "")
    cat("  posterior draws per chain: ", x$chains[[1L]]$backend$n_samples, "\n", sep = "")
    cat("  process terms: ", length(x$chains[[1L]]$backend$process_terms), "\n", sep = "")
  }

  accept <- vapply(x$chains, function(z) z$covariance_acceptance %||% NA_real_, numeric(1))
  if(length(accept))
    cat("  covariance acceptance: ", paste(round(accept, 3), collapse = ", "), "\n", sep = "")

  sampler_time <- x$timing$sampler_total
  if(!is.null(sampler_time) && "elapsed" %in% names(sampler_time))
    cat("  sampler elapsed: ", round(unname(sampler_time[["elapsed"]]), 3), " seconds total\n", sep = "")

  invisible(x)
}

summary.stLMM <- function(object,
                         probs = c(0.025, 0.5, 0.975),
                         burn = 0L,
                         parameters = NULL,
                         ...){

  if(!is.numeric(probs) || anyNA(probs) || any(probs <= 0 | probs >= 1))
    stop("error: probs must contain values in (0,1)")
  burn <- validate_burn(burn, object$backend$n_samples)
  samples <- drop_burn_samples(object$samples, burn = burn)
  samples <- filter_sample_blocks(samples, parameters = parameters)

  out <- list(
    call = object$backend$formula,
    n = object$backend$n,
    n_obs = object$backend$n_obs %||% object$backend$n,
    n_missing_response = object$backend$n_missing_response %||% 0L,
    family = object$backend$family %||% "gaussian",
    n_samples = object$backend$n_samples,
    n_fixed = object$backend$p,
    n_random_coef = object$backend$q,
    n_process = length(object$backend$process_terms),
    residual_model = object$backend$residual_model$type %||% "global_tau",
    burn = burn,
    n_used = object$backend$n_samples - burn,
    parameters = collect_parameter_summaries(samples, probs = probs),
    term_description = object$term_description
  )

  class(out) <- "summary_stLMM"
  out
}

print.summary_stLMM <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  cat("stLMM summary\n")
  if(!is.null(x$call))
    cat("  formula: ", format_model_formula(x$call), "\n", sep = "")
  if(x$n_missing_response > 0L)
    cat("  observations: ", x$n, " (", x$n_obs, " observed, ",
        x$n_missing_response, " missing response)\n", sep = "")
  else
    cat("  observations: ", x$n, "\n", sep = "")
  if(!is.null(x$burn) && x$burn > 0L)
    cat("  posterior draws: ", x$n_samples, " (", x$n_used,
        " used after burn = ", x$burn, ")\n", sep = "")
  else
    cat("  posterior draws: ", x$n_samples, "\n", sep = "")
  cat("  family: ", x$family, "\n", sep = "")
  cat("  fixed effects: ", x$n_fixed, "\n", sep = "")
  cat("  grouped random-effect coefficients: ", x$n_random_coef, "\n", sep = "")
  cat("  process terms: ", x$n_process, "\n", sep = "")
  if(identical(x$family, "binomial") ||
     identical(x$family, "negative_binomial"))
    cat("  residual variance: not used for Polya-Gamma likelihood\n", sep = "")
  else if(identical(x$residual_model, "fixed_variance"))
    cat("  residual variance: fixed row-specific values\n", sep = "")
  else if(identical(x$residual_model, "group_ig_variance"))
    cat("  residual variance: sampled group-specific values\n", sep = "")
  else if(identical(x$residual_model, "scaled_variance"))
    cat("  residual variance: scaled direct-estimate values\n", sep = "")
  else
    cat("  residual variance: global tau_sq\n", sep = "")

  for(nm in names(x$parameters)){
    cat("\n", nm, ":\n", sep = "")
    print(round(x$parameters[[nm]], digits))
  }

  invisible(x)
}

summary.stLMM_chains <- function(object,
                                 probs = c(0.025, 0.5, 0.975),
                                 diagnostics = TRUE,
                                 include_w = FALSE,
                                 burn = 0L,
                                 parameters = NULL,
                                 ...){

  if(!is.numeric(probs) || anyNA(probs) || any(probs <= 0 | probs >= 1))
    stop("error: probs must contain values in (0,1)")
  burn <- validate_burn(burn, object$chains[[1L]]$backend$n_samples)

  chains <- as_mcmc(object, include_w = include_w, burn = burn)
  m <- as.matrix(chains)
  if(!is.null(parameters)){
    parameters <- validate_summary_parameters(parameters, colnames(m))
    m <- m[, parameters, drop = FALSE]
  }

  diagnostics_out <- if(isTRUE(diagnostics)){
    d <- chain_diagnostics(object, include_w = include_w, burn = burn)
    if(!is.null(parameters))
      d <- d[parameters, , drop = FALSE]
    d
  } else {
    NULL
  }

  out <- list(
    call = object$chains[[1L]]$backend$formula,
    n_chains = object$n_chains,
    n = object$chains[[1L]]$backend$n,
    n_obs = object$chains[[1L]]$backend$n_obs %||% object$chains[[1L]]$backend$n,
    n_missing_response = object$chains[[1L]]$backend$n_missing_response %||% 0L,
    family = object$chains[[1L]]$backend$family %||% "gaussian",
    n_samples = object$chains[[1L]]$backend$n_samples,
    n_fixed = object$chains[[1L]]$backend$p,
    n_random_coef = object$chains[[1L]]$backend$q,
    n_process = length(object$chains[[1L]]$backend$process_terms),
    residual_model = object$chains[[1L]]$backend$residual_model$type %||% "global_tau",
    burn = burn,
    n_used = object$chains[[1L]]$backend$n_samples - burn,
    parameters = posterior_summary_table(m, probs = probs),
    diagnostics = diagnostics_out,
    term_description = object$term_description
  )

  class(out) <- "summary_stLMM_chains"
  out
}

print.summary_stLMM_chains <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  cat("stLMM multi-chain summary\n")
  if(!is.null(x$call))
    cat("  formula: ", format_model_formula(x$call), "\n", sep = "")
  cat("  chains: ", x$n_chains, "\n", sep = "")
  cat("  family: ", x$family, "\n", sep = "")
  if(x$n_missing_response > 0L)
    cat("  observations: ", x$n, " (", x$n_obs, " observed, ",
        x$n_missing_response, " missing response)\n", sep = "")
  else
    cat("  observations: ", x$n, "\n", sep = "")
  if(!is.null(x$burn) && x$burn > 0L)
    cat("  posterior draws per chain: ", x$n_samples, " (", x$n_used,
        " used after burn = ", x$burn, ")\n", sep = "")
  else
    cat("  posterior draws per chain: ", x$n_samples, "\n", sep = "")
  cat("  process terms: ", x$n_process, "\n", sep = "")

  cat("\nParameters:\n")
  print(round(x$parameters, digits))

  if(!is.null(x$diagnostics)){
    cat("\nChain diagnostics:\n")
    diag <- x$diagnostics
    num_cols <- vapply(diag, is.numeric, logical(1))
    diag[num_cols] <- lapply(diag[num_cols], round, digits = digits)
    print(diag)
  }

  invisible(x)
}

plot_draw_index <- function(n, burnin = 0L, thin = 1L){

  burnin <- as.integer(burnin[1L])
  thin <- as.integer(thin[1L])

  if(is.na(burnin) || burnin < 0L)
    stop("error: burnin must be a nonnegative integer")
  if(is.na(thin) || thin < 1L)
    stop("error: thin must be a positive integer")
  if(burnin >= n)
    stop("error: burnin removes all posterior draws")

  seq.int(burnin + 1L, n, by = thin)
}

plot.stLMM <- function(x,
                      type = c("trace", "density"),
                      parameters = NULL,
                      max_parameters = 12L,
                      burnin = 0L,
                      thin = 1L,
                      ...){

  type <- match.arg(type)

  samples <- x$samples
  samples$w <- NULL
  active <- Filter(Negate(is.null), samples)
  if(!length(active))
    stop("error: no posterior samples available to plot")

  mats <- vector("list", length(active))
  names(mats) <- names(active)
  for(i in seq_along(active)){
    z <- active[[i]]
    if(is.null(dim(z))){
      mats[[i]] <- matrix(z, ncol = 1L, dimnames = list(NULL, names(active)[i]))
    } else {
      mats[[i]] <- as.matrix(z)
    }
  }
  mat <- do.call(cbind, mats)

  if(is.null(colnames(mat)))
    colnames(mat) <- paste0("param_", seq_len(ncol(mat)))

  if(!is.null(parameters)){
    miss <- setdiff(parameters, colnames(mat))
    if(length(miss))
      stop("error: unknown parameter(s): ", paste(miss, collapse = ", "))
    mat <- mat[, parameters, drop = FALSE]
  }

  if(ncol(mat) > max_parameters)
    mat <- mat[, seq_len(max_parameters), drop = FALSE]

  keep <- plot_draw_index(nrow(mat), burnin = burnin, thin = thin)
  mat <- mat[keep, , drop = FALSE]

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = grDevices::n2mfrow(ncol(mat)))

  for(j in seq_len(ncol(mat))){
    if(type == "trace"){
      graphics::plot(mat[, j], type = "l", xlab = "draw", ylab = colnames(mat)[j], main = colnames(mat)[j], ...)
    } else {
      graphics::plot(stats::density(mat[, j]), xlab = colnames(mat)[j], main = colnames(mat)[j], ...)
    }
  }

  invisible(x)
}

plot.stLMM_chains <- function(x,
                             type = c("trace", "density"),
                             parameters = NULL,
                             max_parameters = 12L,
                             include_w = FALSE,
                             n_col = NULL,
                             chain_colors = NULL,
                             lwd = 1,
                             burnin = 0L,
                             thin = 1L,
  ...){

  type <- match.arg(type)
  plot_draw_index(x$chains[[1L]]$backend$n_samples, burnin = burnin, thin = thin)
  m <- as_mcmc(x, include_w = include_w, burn = burnin, thin = thin)

  if(!is.null(parameters)){
    miss <- setdiff(parameters, colnames(as.matrix(m[[1L]])))
    if(length(miss))
      stop("error: unknown parameter(s): ", paste(miss, collapse = ", "))
    m <- m[, parameters, drop = FALSE]
  }

  if(ncol(as.matrix(m[[1L]])) > max_parameters)
    m <- m[, seq_len(max_parameters), drop = FALSE]

  chain_mats <- lapply(m, as.matrix)
  n_chain <- length(chain_mats)
  param_names <- colnames(chain_mats[[1L]])
  n_param <- length(param_names)

  if(is.null(chain_colors)){
    base_cols <- c("#031131", "#0D388A", "#2A55A0", "#F7A13C",
                   "#C15831", "#A03310", "#6C8FB6", "#4B4B4B")
    chain_colors <- rep(base_cols, length.out = n_chain)
  } else {
    chain_colors <- rep(chain_colors, length.out = n_chain)
  }

  if(is.null(n_col)){
    layout <- grDevices::n2mfrow(n_param)
  } else {
    n_col <- as.integer(n_col[1L])
    if(is.na(n_col) || n_col < 1L)
      stop("error: n_col must be a positive integer")
    layout <- c(ceiling(n_param / n_col), n_col)
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = layout)

  for(j in seq_len(n_param)){
    param_j <- param_names[j]

    if(type == "trace"){
      y_range <- range(vapply(chain_mats, function(z) range(z[, j], finite = TRUE), numeric(2L)),
                       finite = TRUE)
      x_range <- c(1, max(vapply(chain_mats, nrow, integer(1))))

      graphics::plot(
        x_range,
        y_range,
        type = "n",
        xlab = "draw",
        ylab = param_j,
        main = param_j,
        ...
      )

      for(k in seq_len(n_chain)){
        graphics::lines(
          seq_len(nrow(chain_mats[[k]])),
          chain_mats[[k]][, j],
          col = chain_colors[k],
          lwd = lwd
        )
      }
    } else {
      dens <- lapply(chain_mats, function(z) stats::density(z[, j]))
      x_range <- range(vapply(dens, function(z) range(z$x, finite = TRUE), numeric(2L)),
                       finite = TRUE)
      y_range <- range(vapply(dens, function(z) range(z$y, finite = TRUE), numeric(2L)),
                       finite = TRUE)

      graphics::plot(
        x_range,
        y_range,
        type = "n",
        xlab = param_j,
        ylab = "density",
        main = param_j,
        ...
      )

      for(k in seq_len(n_chain))
        graphics::lines(dens[[k]], col = chain_colors[k], lwd = lwd)
    }
  }

  invisible(x)
}

print.stLMM_recovery <- function(x, ...){

  cat("stLMM recovery\n")
  td <- x$term_description
  if(!is.null(td$global$formula))
    cat("  formula: ", format_model_formula(td$global$formula), "\n", sep = "")
  n_missing_response <- x$backend$n_missing_response %||% 0L
  if(n_missing_response > 0L)
    cat("  observations: ", x$backend$n, " (", x$backend$n_obs, " observed, ",
        n_missing_response, " missing response)\n", sep = "")
  else
    cat("  observations: ", x$backend$n, "\n", sep = "")
  cat("  recovered draws: ", length(x$recover_iter), "\n", sep = "")
  cat("  recovered process terms: ", paste(names(x$w_samples), collapse = ", "), "\n", sep = "")
  invisible(x)
}

print.stLMM_recovery_chains <- function(x, ...){

  cat("stLMM multi-chain recovery\n")
  cat("  chains: ", x$n_chains, "\n", sep = "")
  draws <- vapply(x$chains, function(z) length(z$recover_iter), integer(1))
  cat("  recovered draws per chain: ", paste(draws, collapse = ", "), "\n", sep = "")
  if(length(x$chains))
    cat("  recovered process terms: ", paste(names(x$chains[[1L]]$w_samples), collapse = ", "), "\n", sep = "")

  invisible(x)
}

summary.stLMM_recovery <- function(object,
                                  probs = c(0.025, 0.5, 0.975),
                                  burn = 0L,
                                  parameters = NULL,
                                  include_w = FALSE,
                                  max_w = 20L,
                                  ...){

  burn <- validate_burn(burn, object$backend$n_samples)
  out <- summary.stLMM(object, probs = probs, burn = burn,
                       parameters = parameters)
  out$recover_iter <- object$recover_iter
  out$recovered_terms <- names(object$w_samples)
  out$w_node_order <- vapply(object$w_samples, function(z) attr(z, "node_order") %||% NA_character_, character(1))

  if(isTRUE(include_w)){
    out$w <- list()
    keep_w <- object$recover_iter > burn
    if(!any(keep_w))
      stop("error: burn removes all recovered process draws")
    for(nm in names(object$w_samples)){
      w_i <- object$w_samples[[nm]]
      w_i <- w_i[keep_w, , drop = FALSE]
      if(ncol(w_i) > max_w)
        w_i <- w_i[, seq_len(max_w), drop = FALSE]
      out$w[[nm]] <- posterior_summary_table(w_i, probs = probs)
    }
  }

  class(out) <- "summary_stLMM_recovery"
  out
}

print.summary_stLMM_recovery <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  class(x) <- "summary_stLMM"
  print(x, digits = digits, ...)
  cat("\nRecovery:\n")
  cat("  recovered draws: ", length(x$recover_iter), "\n", sep = "")
  cat("  process terms: ", paste(x$recovered_terms, collapse = ", "), "\n", sep = "")
  if(length(x$w_node_order))
    cat("  w_samples order: ", paste(names(x$w_node_order), x$w_node_order, sep = "=", collapse = ", "), "\n", sep = "")

  if(!is.null(x$w)){
    for(nm in names(x$w)){
      cat("\nw:", nm, "\n")
      print(round(x$w[[nm]], digits))
    }
  }

  invisible(x)
}

summary.stLMM_recovery_chains <- function(object,
                                         probs = c(0.025, 0.5, 0.975),
                                         diagnostics = TRUE,
                                         include_w = FALSE,
                                         burn = 0L,
                                         parameters = NULL,
                                         ...){

  out <- summary.stLMM_chains(
    object = object,
    probs = probs,
    diagnostics = diagnostics,
    include_w = include_w,
    burn = burn,
    parameters = parameters,
    ...
  )
  out$recover_iter <- lapply(object$chains, `[[`, "recover_iter")
  out$recovered_terms <- if(length(object$chains)) names(object$chains[[1L]]$w_samples) else character(0)
  class(out) <- "summary_stLMM_recovery_chains"
  out
}

print.summary_stLMM_recovery_chains <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  class(x) <- "summary_stLMM_chains"
  print(x, digits = digits, ...)
  cat("\nRecovery:\n")
  cat("  process terms: ", paste(x$recovered_terms, collapse = ", "), "\n", sep = "")

  invisible(x)
}

plot.stLMM_recovery <- function(x,
                               term = NULL,
                               nodes = NULL,
                               type = c("trace", "density", "fitted"),
                               observed = NULL,
                               max_nodes = 12L,
                               burnin = 0L,
                               thin = 1L,
                               ...){

  type <- match.arg(type)

  if(is.null(x$w_samples) || !length(x$w_samples))
    stop("error: no saved or recovered latent process samples available")

  if(type == "fitted"){
    fit <- stats::fitted(x, summary = FALSE)
    keep <- plot_draw_index(nrow(fit), burnin = burnin, thin = thin)
    fit <- fit[keep, , drop = FALSE]
    fit_mean <- colMeans(fit)
    if(is.null(observed))
      observed <- x$backend$y
    observed <- as.numeric(observed)
    if(length(observed) != length(fit_mean) || anyNA(observed))
      stop("error: observed values must match the fitted data rows with no missing values")

    lim <- range(observed, fit_mean, finite = TRUE)
    graphics::plot(observed, fit_mean, xlab = "observed", ylab = "posterior mean fitted",
                   xlim = lim, ylim = lim, pch = 16, ...)
    graphics::abline(0, 1, col = "gray60")
    return(invisible(x))
  }

  if(is.null(term))
    term <- names(x$w_samples)[1L]
  if(!term %in% names(x$w_samples))
    stop("error: unknown recovered process term ", term)

  w <- x$w_samples[[term]]
  if(is.null(nodes))
    nodes <- seq_len(min(ncol(w), max_nodes))

  if(any(nodes < 1L | nodes > ncol(w)))
    stop("error: nodes are out of bounds for ", term)

  w <- w[, nodes, drop = FALSE]
  keep <- plot_draw_index(nrow(w), burnin = burnin, thin = thin)
  w <- w[keep, , drop = FALSE]

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = grDevices::n2mfrow(ncol(w)))

  for(j in seq_len(ncol(w))){
    if(type == "trace"){
      graphics::plot(w[, j], type = "l", xlab = "recovered draw", ylab = colnames(w)[j], main = colnames(w)[j], ...)
    } else {
      graphics::plot(stats::density(w[, j]), xlab = colnames(w)[j], main = colnames(w)[j], ...)
    }
  }

  invisible(x)
}

plot.stLMM_recovery_chains <- function(x,
                                      type = c("trace", "density"),
                                      parameters = NULL,
                                      max_parameters = 12L,
                                      include_w = FALSE,
                                      ...){

  plot.stLMM_chains(
    x = x,
    type = type,
    parameters = parameters,
    max_parameters = max_parameters,
    include_w = include_w,
    ...
  )
}

print.stLMM_prediction <- function(x, ...){

  cat("stLMM prediction\n")
  cat("  mean samples: ", nrow(x$mu_samples), " draws x ", ncol(x$mu_samples), " rows\n", sep = "")
  cat("  newdata: ", x$newdata, "\n", sep = "")
  cat("  joint: ", x$joint, "\n", sep = "")
  cat("  y samples: ", if(is.null(x$y_samples)) "not simulated" else "simulated", "\n", sep = "")
  cat("  process samples: ",
      if(is.null(x$w_samples) || !length(x$w_samples)) "not retained" else paste(names(x$w_samples), collapse = ", "),
      "\n", sep = "")

  invisible(x)
}

print.stLMM_prediction_chains <- function(x, ...){

  cat("stLMM multi-chain prediction\n")
  cat("  chains: ", x$n_chains, "\n", sep = "")
  draws <- vapply(x$chains, function(z) nrow(z$mu_samples), integer(1))
  rows <- if(length(x$chains)) ncol(x$chains[[1L]]$mu_samples) else 0L
  cat("  mean samples per chain: ", paste(draws, collapse = ", "), " draws x ", rows, " rows\n", sep = "")
  if(length(x$chains)){
    cat("  newdata: ", x$chains[[1L]]$newdata, "\n", sep = "")
    cat("  joint: ", x$chains[[1L]]$joint, "\n", sep = "")
    cat("  y samples: ", if(is.null(x$chains[[1L]]$y_samples)) "not simulated" else "simulated", "\n", sep = "")
    w_terms <- if(is.null(x$chains[[1L]]$w_samples) || !length(x$chains[[1L]]$w_samples)) character(0) else names(x$chains[[1L]]$w_samples)
    cat("  process samples: ", if(!length(w_terms)) "not retained" else paste(w_terms, collapse = ", "), "\n", sep = "")
  }

  invisible(x)
}

summary.stLMM_prediction <- function(object,
                                    probs = c(0.025, 0.5, 0.975),
                                    include_y = !is.null(object$y_samples),
                                    ...){

  if(!is.numeric(probs) || anyNA(probs) || any(probs <= 0 | probs >= 1))
    stop("error: probs must contain values in (0,1)")

  out <- list(
    n_draw = nrow(object$mu_samples),
    n_row = ncol(object$mu_samples),
    draw_index = object$draw_index,
    newdata = object$newdata,
    joint = object$joint,
    w_terms = if(is.null(object$w_samples) || !length(object$w_samples)) character(0) else names(object$w_samples),
    mu = posterior_summary_table(object$mu_samples, probs = probs),
    y = if(isTRUE(include_y) && !is.null(object$y_samples)) posterior_summary_table(object$y_samples, probs = probs) else NULL
  )

  class(out) <- "summary_stLMM_prediction"
  out
}

print.summary_stLMM_prediction <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  cat("stLMM prediction summary\n")
  cat("  draws: ", x$n_draw, "\n", sep = "")
  cat("  rows: ", x$n_row, "\n", sep = "")
  cat("  newdata: ", x$newdata, "\n", sep = "")
  cat("  process samples: ", if(!length(x$w_terms)) "not retained" else paste(x$w_terms, collapse = ", "), "\n", sep = "")

  cat("\nmu:\n")
  print(round(x$mu, digits))

  if(!is.null(x$y)){
    cat("\ny:\n")
    print(round(x$y, digits))
  }

  invisible(x)
}

summary.stLMM_prediction_chains <- function(object,
                                           probs = c(0.025, 0.5, 0.975),
                                           include_y = !is.null(object$chains[[1L]]$y_samples),
                                           ...){

  if(!is.numeric(probs) || anyNA(probs) || any(probs <= 0 | probs >= 1))
    stop("error: probs must contain values in (0,1)")

  mu <- do.call(rbind, lapply(object$chains, `[[`, "mu_samples"))
  y <- if(isTRUE(include_y) && !is.null(object$chains[[1L]]$y_samples))
    do.call(rbind, lapply(object$chains, `[[`, "y_samples"))
  else
    NULL

  out <- list(
    n_chains = object$n_chains,
    n_draw = nrow(mu),
    n_row = ncol(mu),
    newdata = object$chains[[1L]]$newdata,
    joint = object$chains[[1L]]$joint,
    mu = posterior_summary_table(mu, probs = probs),
    y = if(!is.null(y)) posterior_summary_table(y, probs = probs) else NULL
  )

  class(out) <- "summary_stLMM_prediction_chains"
  out
}

print.summary_stLMM_prediction_chains <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  cat("stLMM multi-chain prediction summary\n")
  cat("  chains: ", x$n_chains, "\n", sep = "")
  cat("  pooled draws: ", x$n_draw, "\n", sep = "")
  cat("  rows: ", x$n_row, "\n", sep = "")
  cat("  newdata: ", x$newdata, "\n", sep = "")

  cat("\nmu:\n")
  print(round(x$mu, digits))

  if(!is.null(x$y)){
    cat("\ny:\n")
    print(round(x$y, digits))
  }

  invisible(x)
}

plot.stLMM_prediction <- function(x,
                                 sample = c("mu", "y"),
                                 type = c("interval", "density", "scatter"),
                                 rows = NULL,
                                 observed = NULL,
                                 probs = c(0.025, 0.5, 0.975),
                                 max_rows = 200L,
                                 burnin = 0L,
                                 thin = 1L,
                                 ...){

  sample <- match.arg(sample)
  type <- match.arg(type)

  mat <- if(sample == "mu") x$mu_samples else x$y_samples
  if(is.null(mat))
    stop("error: requested prediction sample is not available")

  if(is.null(rows))
    rows <- seq_len(min(ncol(mat), max_rows))
  if(any(rows < 1L | rows > ncol(mat)))
    stop("error: rows are out of bounds")

  mat <- mat[, rows, drop = FALSE]
  keep <- plot_draw_index(nrow(mat), burnin = burnin, thin = thin)
  mat <- mat[keep, , drop = FALSE]

  if(type == "scatter"){
    if(is.null(observed))
      stop("error: observed values are required for scatter plots")
    observed <- as.numeric(observed)
    if(length(observed) != ncol(mat) || anyNA(observed))
      stop("error: observed values must match the selected prediction rows with no missing values")
    pred_mean <- colMeans(mat)
    lim <- range(observed, pred_mean, finite = TRUE)
    graphics::plot(observed, pred_mean, xlab = "observed", ylab = paste("posterior mean", sample),
                   xlim = lim, ylim = lim, pch = 16, ...)
    graphics::abline(0, 1, col = "gray60")
  } else if(type == "interval"){
    center <- apply(mat, 2L, stats::quantile, probs = probs[2L], names = FALSE)
    lower <- apply(mat, 2L, stats::quantile, probs = probs[1L], names = FALSE)
    upper <- apply(mat, 2L, stats::quantile, probs = probs[3L], names = FALSE)
    xx <- seq_along(center)
    ylim <- range(lower, upper, finite = TRUE)
    graphics::plot(xx, center, ylim = ylim, type = "n", xlab = "prediction row", ylab = sample, ...)
    graphics::segments(xx, lower, xx, upper, col = "gray60")
    graphics::points(xx, center, pch = 16)
  } else {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mfrow = grDevices::n2mfrow(ncol(mat)))
    for(j in seq_len(ncol(mat)))
      graphics::plot(stats::density(mat[, j]), xlab = colnames(mat)[j] %||% paste0("row_", rows[j]), main = paste0(sample, " ", rows[j]), ...)
  }

  invisible(x)
}

plot.stLMM_prediction_chains <- function(x,
                                        sample = c("mu", "y"),
                                        type = c("interval", "density"),
                                        rows = NULL,
                                        probs = c(0.025, 0.5, 0.975),
                                        max_rows = 200L,
                                        burnin = 0L,
                                        thin = 1L,
                                        ...){

  sample <- match.arg(sample)
  type <- match.arg(type)
  mat <- if(sample == "mu")
    do.call(rbind, lapply(x$chains, `[[`, "mu_samples"))
  else if(!is.null(x$chains[[1L]]$y_samples))
    do.call(rbind, lapply(x$chains, `[[`, "y_samples"))
  else
    stop("error: y samples are not available")

  obj <- x$chains[[1L]]
  obj[[paste0(sample, "_samples")]] <- mat
  plot.stLMM_prediction(
    x = obj,
    sample = sample,
    type = type,
    rows = rows,
    probs = probs,
    max_rows = max_rows,
    burnin = burnin,
    thin = thin,
    ...
  )
}
