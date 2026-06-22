recover <- function(object, ...)
  UseMethod("recover")

recover.stLMM <- function(object,
                         n_omp_threads = NULL,
                         verbose = FALSE,
                         sub_sample = list(start = 1L, thin = 1L),
                         ...){

  if(!is.list(object) || is.null(object$backend))
    stop("error: object must be an stLMM fit")
  n_omp_threads <- resolve_n_omp_threads(object, n_omp_threads)

  if(length(object$backend$process_terms) == 0L)
    stop("error: recover() requires at least one structured process term")

  if(!is.list(sub_sample))
    stop("error: sub_sample must be a list with optional entries 'start' and 'thin'")
  if("pg_iter" %in% names(sub_sample))
    stop("error: sub_sample$pg_iter is no longer supported; Polya-Gamma process models must save process draws during stLMM() via save_process")

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

  n_samples <- as.integer(object$backend$n_samples)
  draw_index <- seq.int(start, n_samples, by = thin)
  if(length(draw_index) == 0L)
    stop("error: sub_sample selects no posterior draws")
  stlmm_progress(
    verbose,
    "recover: selected ", length(draw_index), " posterior draw(s) using ",
    n_omp_threads, " OpenMP thread(s)"
  )

  object$backend$n_omp_threads <- n_omp_threads

  is_pg_likelihood <- identical(object$backend$family, "binomial") ||
    identical(object$backend$family, "negative_binomial")

  if(is_pg_likelihood){
    stlmm_progress(verbose, "recover: selecting saved in-chain process draws")
    if(is.null(object$recover_iter) || !length(object$recover_iter) ||
       is.null(object$w_samples_stacked) || !is.matrix(object$w_samples_stacked))
      stop(
        "error: Polya-Gamma process recovery requires process draws saved during stLMM(); ",
        "refit with save_process = list(start = ..., thin = ...)"
      )

    keep <- which(object$recover_iter %in% draw_index)
    if(length(keep) == 0L)
      stop("error: sub_sample selects no saved process draws")
    recover_iter <- as.integer(object$recover_iter[keep])

    object$w_samples_stacked <- object$w_samples_stacked[keep, , drop = FALSE]
    object$w_samples_ordered <- unstack_w_samples(
      w_samples_stacked = object$w_samples_stacked,
      process_terms = object$backend$process_terms,
      term_description = object$term_description,
      graphs = object$backend$graphs,
      user_order = FALSE
    )
    object$w_samples <- unstack_w_samples(
      w_samples_stacked = object$w_samples_stacked,
      process_terms = object$backend$process_terms,
      term_description = object$term_description,
      graphs = object$backend$graphs,
      user_order = TRUE
    )
    object$recover_iter <- recover_iter
    object$samples$w <- object$w_samples

    class(object) <- unique(c("stLMM_recovery", class(object)))
    stlmm_progress(verbose, "recover: done")
    return(object)
  }

  alpha_samples <- object$alpha_samples
  if(is.null(alpha_samples))
    alpha_samples <- matrix(0.0, nrow = n_samples, ncol = 0L)

  tau_sq_samples <- object$tau_sq_samples
  if(is.null(tau_sq_samples))
    tau_sq_samples <- rep(1.0, n_samples)

  residual_variance_samples <- object$residual_variance_samples
  if(is.null(residual_variance_samples))
    residual_variance_samples <- matrix(0.0, nrow = n_samples, ncol = 0L)

  stlmm_progress(verbose, "recover: sampling latent process draws")
  out <- .Call(
    "stLMM_recover_w",
    object$backend,
    object$beta_samples,
    alpha_samples,
    tau_sq_samples,
    residual_variance_samples,
    object$sigma_sq_samples,
    object$theta_samples,
    as.integer(draw_index),
    1L,
    PACKAGE = "stLMM"
  )

  object$w_samples_stacked <- out$w_samples
  object$w_samples_ordered <- unstack_w_samples(
    w_samples_stacked = out$w_samples,
    process_terms = object$backend$process_terms,
    term_description = object$term_description,
    graphs = object$backend$graphs,
    user_order = FALSE
  )
  object$w_samples <- unstack_w_samples(
    w_samples_stacked = out$w_samples,
    process_terms = object$backend$process_terms,
    term_description = object$term_description,
    graphs = object$backend$graphs,
    user_order = TRUE
  )
  object$recover_iter <- out$recover_iter
  object$samples$w <- object$w_samples

  class(object) <- unique(c("stLMM_recovery", class(object)))
  stlmm_progress(verbose, "recover: done")
  object
}

recover.stLMM_chains <- function(object,
                                 n_omp_threads = NULL,
                                 verbose = FALSE,
                                 sub_sample = list(start = 1L, thin = 1L),
                                 ...){

  recovered <- vector("list", length(object$chains))
  for(i in seq_along(object$chains)){
    stlmm_progress(verbose, "recover: chain ", i, " of ", length(object$chains))
    recovered[[i]] <- recover(
      object$chains[[i]],
      n_omp_threads = n_omp_threads,
      verbose = verbose,
      sub_sample = sub_sample,
      ...
    )
  }
  out <- list(
    chains = recovered,
    n_chains = object$n_chains,
    chain_control = object$chain_control,
    call = match.call(),
    fit_call = object$call,
    term_description = object$term_description
  )
  class(out) <- "stLMM_recovery_chains"
  out
}

unstack_w_samples <- function(w_samples_stacked,
                              process_terms,
                              term_description,
                              graphs = NULL,
                              user_order = FALSE){

  process_names <- vapply(process_terms, `[[`, character(1), "name")
  w_samples_by_term <- vector("list", length(process_terms))
  names(w_samples_by_term) <- process_names

  td_process <- term_description$process_terms

  for(i in seq_along(process_terms)){
    q_lat_i <- as.integer(td_process[[i]]$q_lat)
    w_offset_i <- as.integer(td_process[[i]]$w_offset)

    if(is.na(q_lat_i) || is.na(w_offset_i) || q_lat_i < 1L || w_offset_i < 1L)
      stop("internal error: missing recovered process indexing information for ", process_names[i])

    idx <- w_offset_i + seq_len(q_lat_i) - 1L

    term_samples <- w_samples_stacked[, idx, drop = FALSE]

    if(isTRUE(user_order) &&
       process_terms[[i]]$term_type %in% c("nngp", "dagar", "dagar_time") &&
       !is.null(graphs)){
      graph_i <- graphs[[process_terms[[i]]$graph_index]]
      if(!is.null(graph_i$ord_inv)){
        if(process_terms[[i]]$term_type == "dagar_time"){
          n_time <- as.integer(graph_i$n_time %||% NA_integer_)
          if(is.na(n_time) || n_time < 1L || length(graph_i$ord_inv) * n_time != ncol(term_samples))
            stop("internal error: malformed DAGAR-time ordering metadata for ", process_names[i])
          reorder <- as.integer(unlist(lapply(
            graph_i$ord_inv,
            function(space_index) (space_index - 1L) * n_time + seq_len(n_time)
          )))
          term_samples <- term_samples[, reorder, drop = FALSE]
        } else {
          if(length(graph_i$ord_inv) != ncol(term_samples))
            stop("internal error: malformed process ordering metadata for ", process_names[i])
          term_samples <- term_samples[, graph_i$ord_inv, drop = FALSE]
        }
      }
    }

    colnames(term_samples) <- paste0(process_names[i], "_", seq_len(ncol(term_samples)))
    attr(term_samples, "node_order") <- if(isTRUE(user_order)) "support" else "internal"
    w_samples_by_term[[i]] <- term_samples
  }

  w_samples_by_term
}
