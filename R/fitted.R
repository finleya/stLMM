fitted.stLMM <- function(object,
                        summary = TRUE,
                        sub_sample = list(start = 1L, thin = 1L),
                        scale = c("response", "link"),
                        ...){

  if(!is.list(object) || is.null(object$backend))
    stop("error: object must be an stLMM fit")

  scale <- match.arg(scale)
  n <- as.integer(object$backend$n)
  process_terms <- object$backend$process_terms
  has_process <- length(process_terms) > 0L

  if(!is.list(sub_sample))
    stop("error: sub_sample must be a list with optional entries 'start' and 'thin'")

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

  if(has_process){
    w_samples_ordered <- object$w_samples_ordered
    if(is.null(w_samples_ordered))
      w_samples_ordered <- object$w_samples

    if(is.null(object$recover_iter) || length(object$recover_iter) == 0L ||
       is.null(w_samples_ordered) || !is.list(w_samples_ordered))
      stop("error: fitted values with process terms require saved or recovered latent process samples; call recover() on the fitted object first")

    keep <- which(object$recover_iter >= start & ((object$recover_iter - start) %% thin) == 0L)
    if(length(keep) == 0L)
      stop("error: sub_sample selects no saved or recovered latent process draws")

    draw_index <- as.integer(object$recover_iter[keep])
    recover_row <- keep
  } else {
    n_samples <- nrow(object$beta_samples)
    if(is.null(n_samples))
      n_samples <- length(object$tau_sq_samples)

    draw_index <- seq.int(start, n_samples, by = thin)
    if(length(draw_index) == 0L)
      stop("error: sub_sample selects no posterior draws")

    recover_row <- integer(0)
  }

  n_draw <- length(draw_index)
  mu <- matrix(0.0, nrow = n_draw, ncol = n)
  colnames(mu) <- rownames(object$backend$X)

  if(!is.null(object$beta_samples) && ncol(object$beta_samples) > 0L){
    beta_draws <- object$beta_samples[draw_index, , drop = FALSE]
    mu <- mu + beta_draws %*% t(object$backend$X)
  }

  if(!is.null(object$alpha_samples) && ncol(object$alpha_samples) > 0L){
    alpha_draws <- object$alpha_samples[draw_index, , drop = FALSE]
    mu <- mu + as.matrix(alpha_draws %*% Matrix::t(object$backend$Z))
  }

  if(has_process){
    process_names <- vapply(process_terms, `[[`, character(1), "name")

    for(i in seq_along(process_terms)){
      term <- process_terms[[i]]
      term_name <- process_names[i]
      w_i <- w_samples_ordered[[term_name]]

      if(is.null(w_i))
        stop("error: saved or recovered latent process samples missing for ", term_name)
      if(nrow(w_i) < max(recover_row))
        stop("error: recovered latent process sample rows do not align with recover_iter")

      map <- as.integer(term$map)
      if(anyNA(map) || length(map) != n)
        stop("error: malformed process map for ", term_name)

      term_mu <- w_i[recover_row, map, drop = FALSE]

      if(!is.null(term$x)){
        term_scale <- as.numeric(term$x)
        if(length(term_scale) != n || anyNA(term_scale))
          stop("error: malformed SVC scale vector for ", term_name)
        term_mu <- sweep(term_mu, 2L, term_scale, `*`)
      }

      mu <- mu + term_mu
    }
  }

  offset <- object$backend$offset
  if(!is.null(offset)){
    if(length(offset) != n)
      stop("error: fitted offset length does not match fitted rows")
    mu <- sweep(mu, 2L, as.numeric(offset), `+`)
  }

  if(identical(scale, "response")){
    if(identical(object$backend$family, "binomial"))
      mu <- stats::plogis(mu)
    else if(identical(object$backend$family, "negative_binomial"))
      mu <- exp(mu)
  }

  if(isTRUE(summary))
    return(as.numeric(colMeans(mu)))

  attr(mu, "draw_index") <- draw_index
  if(identical(object$backend$family, "binomial") ||
     identical(object$backend$family, "negative_binomial"))
    attr(mu, "scale") <- scale
  mu
}

fitted.stLMM_recovery_chains <- function(object,
                                        summary = TRUE,
                                        sub_sample = list(start = 1L, thin = 1L),
                                        scale = c("response", "link"),
                                        ...){

  scale <- match.arg(scale)
  fits <- lapply(object$chains, stats::fitted, summary = FALSE,
                 sub_sample = sub_sample, scale = scale, ...)

  if(isTRUE(summary))
    return(as.numeric(colMeans(do.call(rbind, fits))))

  out <- fits
  class(out) <- "stLMM_fitted_chains"
  out
}
