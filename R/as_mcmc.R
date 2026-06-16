as_mcmc <- function(object, ...){
  if(inherits(object, "stLMM_recovery_chains"))
    return(as_mcmc.stLMM_recovery_chains(object, ...))
  if(inherits(object, "stLMM_chains"))
    return(as_mcmc.stLMM_chains(object, ...))
  if(inherits(object, "stLMM_recovery"))
    return(as_mcmc.stLMM_recovery(object, ...))
  if(inherits(object, "stLMM"))
    return(as_mcmc.stLMM(object, ...))
  stop("error: no as_mcmc method for object of class ", paste(class(object), collapse = ", "))
}

sample_matrix <- function(object, include_w = FALSE){

  samples <- object$samples
  if(is.null(samples))
    stop("error: object has no posterior samples")

  if(!isTRUE(include_w))
    samples$w <- NULL

  active <- Filter(function(z){
    if(is.null(z))
      return(FALSE)
    if(!is.null(dim(z)) && ncol(as.matrix(z)) == 0L)
      return(FALSE)
    TRUE
  }, samples)
  if(!length(active))
    stop("error: no posterior samples available")

  mats <- vector("list", length(active))
  names(mats) <- names(active)

  for(i in seq_along(active)){
    block <- active[[i]]

    if(is.list(block) && !is.data.frame(block)){
      block <- do.call(cbind, lapply(names(block), function(nm){
        mat <- as.matrix(block[[nm]])
        if(is.null(colnames(mat)))
          colnames(mat) <- paste0(nm, "_", seq_len(ncol(mat)))
        mat
      }))
    } else if(is.null(dim(block))) {
      block <- matrix(as.numeric(block), ncol = 1L)
      colnames(block) <- names(active)[i]
    } else {
      block <- as.matrix(block)
    }

    if(ncol(block) == 0L)
      next

    if(is.null(colnames(block)))
      colnames(block) <- paste0(names(active)[i], "_", seq_len(ncol(block)))

    mats[[i]] <- block
  }

  n_row <- vapply(mats, nrow, integer(1))
  if(length(unique(n_row)) != 1L)
    stop("error: sample blocks do not have a common number of rows")

  do.call(cbind, mats)
}

as_mcmc.stLMM <- function(object,
                         include_w = FALSE,
                         burn = 0L,
                         thin = 1L,
                         ...){

  mat <- sample_matrix(object, include_w = include_w)
  keep <- sample_draw_index(nrow(mat), burn = burn, thin = thin)
  coda::mcmc(mat[keep, , drop = FALSE])
}

as_mcmc.stLMM_recovery <- function(object,
                                  include_w = FALSE,
                                  burn = 0L,
                                  thin = 1L,
                                  ...){
  if(isTRUE(include_w)){
    return(coda::mcmc(as.matrix(as_samples(
      object,
      burn = burn,
      thin = thin,
      metadata = FALSE,
      include_w = TRUE,
      ...
    ))))
  }

  as_mcmc.stLMM(object, include_w = include_w, burn = burn, thin = thin, ...)
}

as_mcmc.stLMM_chains <- function(object,
                                include_w = FALSE,
                                burn = 0L,
                                thin = 1L,
                                ...){

  do.call(coda::mcmc.list, lapply(
    object$chains,
    as_mcmc,
    include_w = include_w,
    burn = burn,
    thin = thin,
    ...
  ))
}

as_mcmc.stLMM_recovery_chains <- function(object,
                                         include_w = FALSE,
                                         burn = 0L,
                                         thin = 1L,
                                         ...){

  do.call(coda::mcmc.list, lapply(
    object$chains,
    as_mcmc,
    include_w = include_w,
    burn = burn,
    thin = thin,
    ...
  ))
}

validate_burn <- function(burn, n = NULL){

  if(!is.numeric(burn) || length(burn) != 1L || is.na(burn) ||
     !is.finite(burn) || burn < 0)
    stop("error: burn must be a nonnegative integer")
  if(abs(burn - round(burn)) > sqrt(.Machine$double.eps))
    stop("error: burn must be a nonnegative integer")

  burn <- as.integer(round(burn))
  if(!is.null(n) && burn >= n)
    stop("error: burn removes all posterior draws")

  burn
}

drop_burn_matrix <- function(x, burn){

  x <- as.matrix(x)
  burn <- validate_burn(burn, nrow(x))
  if(burn == 0L)
    return(x)

  x[seq.int(burn + 1L, nrow(x)), , drop = FALSE]
}

drop_burn_block <- function(x, burn){

  if(is.null(x))
    return(NULL)

  if(is.list(x) && !is.data.frame(x))
    return(lapply(x, drop_burn_block, burn = burn))

  if(is.null(dim(x))){
    burn <- validate_burn(burn, length(x))
    if(burn == 0L)
      return(x)
    return(x[seq.int(burn + 1L, length(x))])
  }

  drop_burn_matrix(x, burn = burn)
}

drop_burn_samples <- function(samples, burn){
  lapply(samples, drop_burn_block, burn = burn)
}

drop_burn_mcmc <- function(x, burn){

  if(inherits(x, "mcmc.list")){
    chains <- lapply(x, drop_burn_mcmc, burn = burn)
    return(do.call(coda::mcmc.list, chains))
  }

  mat <- drop_burn_matrix(x, burn = burn)
  coda::mcmc(mat)
}

chain_diagnostics <- function(object, include_w = FALSE, burn = 0L){

  chains <- as_mcmc(object, include_w = include_w, burn = burn)
  if(!inherits(chains, "mcmc.list"))
    stop("error: chain diagnostics require a multi-chain object")

  rhat <- tryCatch(
    coda::gelman.diag(chains, autoburnin = FALSE, multivariate = FALSE)$psrf[, "Point est."],
    error = function(e) NULL
  )
  ess <- tryCatch(coda::effectiveSize(chains), error = function(e) NULL)

  params <- colnames(as.matrix(chains[[1L]]))
  out <- data.frame(
    parameter = params,
    rhat = NA_real_,
    effective_size = NA_real_,
    row.names = params,
    check.names = FALSE
  )

  if(!is.null(rhat))
    out[names(rhat), "rhat"] <- as.numeric(rhat)
  if(!is.null(ess))
    out[names(ess), "effective_size"] <- as.numeric(ess)

  out
}
