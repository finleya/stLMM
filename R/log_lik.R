log_lik <- function(object, ...)
  UseMethod("log_lik")

waic <- function(object, ...)
  UseMethod("waic")

waic.default <- function(object, ...){
  require_loo()
  loo::waic(object, ...)
}

require_loo <- function(){
  if(!requireNamespace("loo", quietly = TRUE))
    stop("error: package 'loo' is required for waic(); install it with install.packages(\"loo\")")
  invisible(TRUE)
}

log_lik.stLMM <- function(object,
                          sub_sample = list(start = 1L, thin = 1L),
                          ...){

  if(!is.list(object) || is.null(object$backend))
    stop("error: object must be an stLMM fit")

  has_process <- length(object$backend$process_terms) > 0L
  if(has_process){
    w_samples_ordered <- object$w_samples_ordered
    if(is.null(w_samples_ordered))
      w_samples_ordered <- object$w_samples
    if(is.null(object$recover_iter) || length(object$recover_iter) == 0L ||
       is.null(w_samples_ordered) || !is.list(w_samples_ordered))
      stop("error: log_lik() with process terms requires saved or recovered latent process samples; call recover() on the fitted object first")
  }

  eta <- fitted(object, summary = FALSE, sub_sample = sub_sample, scale = "link")
  draw_index <- attr(eta, "draw_index")
  if(is.null(draw_index))
    draw_index <- seq_len(nrow(eta))

  out <- pointwise_log_lik_from_eta(object, eta, draw_index)
  attr(out, "draw_index") <- as.integer(draw_index)
  attr(out, "observed_index") <- as.integer(object$backend$observed_index)
  out
}

log_lik.stLMM_recovery <- function(object,
                                   sub_sample = list(start = 1L, thin = 1L),
                                   ...)
  log_lik.stLMM(object, sub_sample = sub_sample, ...)

log_lik.stLMM_chains <- function(object,
                                 sub_sample = list(start = 1L, thin = 1L),
                                 ...){
  if(!is.list(object) || is.null(object$chains))
    stop("error: object must be an stLMM_chains fit")

  ll <- lapply(object$chains, log_lik, sub_sample = sub_sample, ...)
  out <- do.call(rbind, ll)
  attr(out, "chain") <- rep(seq_along(ll), vapply(ll, nrow, integer(1)))
  attr(out, "draw_index") <- as.integer(unlist(lapply(ll, attr, "draw_index"),
                                               use.names = FALSE))
  attr(out, "observed_index") <- attr(ll[[1L]], "observed_index")
  out
}

log_lik.stLMM_recovery_chains <- function(object,
                                          sub_sample = list(start = 1L, thin = 1L),
                                          ...)
  log_lik.stLMM_chains(object, sub_sample = sub_sample, ...)

waic.stLMM <- function(object,
                       sub_sample = list(start = 1L, thin = 1L),
                       ...){
  require_loo()
  loo::waic(log_lik(object, sub_sample = sub_sample), ...)
}

waic.stLMM_recovery <- function(object,
                                sub_sample = list(start = 1L, thin = 1L),
                                ...)
  waic.stLMM(object, sub_sample = sub_sample, ...)

waic.stLMM_chains <- function(object,
                              sub_sample = list(start = 1L, thin = 1L),
                              ...){
  require_loo()
  loo::waic(log_lik(object, sub_sample = sub_sample), ...)
}

waic.stLMM_recovery_chains <- function(object,
                                       sub_sample = list(start = 1L, thin = 1L),
                                       ...)
  waic.stLMM_chains(object, sub_sample = sub_sample, ...)

pointwise_log_lik_from_eta <- function(object, eta, draw_index){

  backend <- object$backend
  obs <- as.integer(backend$observed_index)
  if(length(obs) == 0L)
    stop("error: WAIC requires at least one observed response")

  eta_obs <- eta[, obs, drop = FALSE]
  y_obs <- as.numeric(backend$y[obs])
  n_draw <- nrow(eta_obs)
  n_obs <- length(obs)
  family <- backend$family %||% "gaussian"

  y_mat <- matrix(rep(y_obs, each = n_draw), nrow = n_draw, ncol = n_obs)

  if(identical(family, "binomial")){
    trials <- backend$trials_obs
    if(is.null(trials))
      trials <- backend$likelihood$trials_obs
    if(is.null(trials))
      trials <- backend$trials[obs]
    if(is.null(trials))
      trials <- rep.int(1L, n_obs)

    trials_mat <- matrix(rep(as.integer(trials), each = n_draw),
                         nrow = n_draw, ncol = n_obs)
    out <- matrix(
      stats::dbinom(as.vector(y_mat),
                    size = as.vector(trials_mat),
                    prob = stats::plogis(as.vector(eta_obs)),
                    log = TRUE),
      nrow = n_draw,
      ncol = n_obs
    )
  } else if(identical(family, "negative_binomial")) {
    size <- backend$nb_size
    if(is.null(size))
      size <- backend$likelihood$size
    if(is.null(size) || length(size) != 1L || is.na(size) ||
       !is.finite(size) || size <= 0)
      stop("error: fitted negative-binomial size is missing or invalid")

    out <- matrix(
      stats::dnbinom(as.vector(y_mat),
                     size = as.double(size),
                     mu = exp(as.vector(eta_obs)),
                     log = TRUE),
      nrow = n_draw,
      ncol = n_obs
    )
  } else if(identical(family, "gaussian")) {
    sd <- fitted_residual_sd_samples(object, draw_index)
    sd_obs <- sd[, obs, drop = FALSE]
    out <- matrix(
      stats::dnorm(as.vector(y_mat),
                   mean = as.vector(eta_obs),
                   sd = as.vector(sd_obs),
                   log = TRUE),
      nrow = n_draw,
      ncol = n_obs
    )
  } else {
    stop("error: unsupported family for log_lik(): ", family)
  }

  colnames(out) <- rownames(backend$X)[obs]
  out
}

fitted_residual_sd_samples <- function(object, draw_index){

  n <- as.integer(object$backend$n)
  sd <- prediction_residual_sd_samples(
    object = object,
    newdata = NULL,
    n0 = n,
    draw_index = draw_index,
    where = "WAIC"
  )

  if(!is.null(sd))
    return(sd)

  if(is.null(object$tau_sq_samples))
    stop("error: tau_sq samples are missing from fitted object")
  tau_sq <- object$tau_sq_samples[draw_index]
  if(anyNA(tau_sq) || any(!is.finite(tau_sq)) || any(tau_sq <= 0))
    stop("error: tau_sq samples must be finite and positive for WAIC")

  matrix(rep(sqrt(tau_sq), n), nrow = length(draw_index), ncol = n)
}
