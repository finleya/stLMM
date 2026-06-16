as_samples <- function(object, ...)
  UseMethod("as_samples")

sample_draw_index <- function(n, burn = 0L, thin = 1L){

  burn <- validate_burn(burn, n)
  if(!is.numeric(thin) || length(thin) != 1L || is.na(thin) ||
     !is.finite(thin) || thin < 1)
    stop("error: thin must be a positive integer")
  if(abs(thin - round(thin)) > sqrt(.Machine$double.eps))
    stop("error: thin must be a positive integer")

  thin <- as.integer(round(thin))
  seq.int(burn + 1L, n, by = thin)
}

validate_thin <- function(thin){

  if(!is.numeric(thin) || length(thin) != 1L || is.na(thin) ||
     !is.finite(thin) || thin < 1)
    stop("error: thin must be a positive integer")
  if(abs(thin - round(thin)) > sqrt(.Machine$double.eps))
    stop("error: thin must be a positive integer")

  as.integer(round(thin))
}

sample_available_index <- function(iteration,
                                   burn = 0L,
                                   thin = 1L,
                                   what = "posterior draws"){

  burn <- validate_burn(burn)
  thin <- validate_thin(thin)

  keep <- which(iteration > burn)
  if(length(keep) == 0L)
    stop("error: burn removes all ", what)

  keep[seq.int(1L, length(keep), by = thin)]
}

samples_make_names <- function(prefix, names, n){

  if(is.null(names))
    names <- character(n)
  names <- as.character(names)
  missing <- is.na(names) | !nzchar(names)
  names[missing] <- paste0(prefix, "_", which(missing))
  names
}

samples_matrix_block <- function(block, prefix, index = NULL){

  if(is.null(block))
    return(NULL)

  if(is.null(dim(block))){
    block <- matrix(as.numeric(block), ncol = 1L)
    colnames(block) <- prefix
  } else {
    block <- as.matrix(block)
  }

  if(ncol(block) == 0L)
    return(NULL)

  if(!is.null(index))
    block <- block[index, , drop = FALSE]

  colnames(block) <- samples_make_names(prefix, colnames(block), ncol(block))
  block
}

samples_w_blocks <- function(w_samples, index = NULL){

  if(is.null(w_samples) || !length(w_samples))
    return(NULL)

  out <- list()
  for(nm in names(w_samples)){
    block <- samples_matrix_block(w_samples[[nm]], paste0("w_", nm), index = index)
    if(!is.null(block)){
      colnames(block) <- paste0("w_", nm, "_", seq_len(ncol(block)))
      out[[nm]] <- block
    }
  }

  if(!length(out))
    return(NULL)

  do.call(cbind, out)
}

stLMM_parameter_samples <- function(object, index = NULL){

  samples <- object$samples
  if(is.null(samples))
    stop("error: object has no posterior samples")
  samples$w <- NULL

  blocks <- list()
  for(nm in names(samples)){
    block <- samples_matrix_block(samples[[nm]], nm, index = index)
    if(!is.null(block))
      blocks[[nm]] <- block
  }

  if(!length(blocks))
    stop("error: no posterior samples available")

  n_row <- vapply(blocks, nrow, integer(1))
  if(length(unique(n_row)) != 1L)
    stop("error: sample blocks do not have a common number of rows")

  do.call(cbind, blocks)
}

samples_data_frame <- function(mat,
                               iteration,
                               chain = 1L,
                               metadata = TRUE,
                               check.names = FALSE){

  if(!is.matrix(mat))
    mat <- as.matrix(mat)
  if(nrow(mat) != length(iteration))
    stop("error: sample metadata does not match sample rows")

  out <- as.data.frame(mat, check.names = check.names, optional = TRUE)

  if(isTRUE(metadata)){
    meta <- data.frame(
      .chain = rep.int(as.integer(chain), nrow(mat)),
      .iteration = as.integer(iteration),
      check.names = FALSE
    )
    out <- cbind(meta, out)
  }

  out
}

as_samples.stLMM <- function(object,
                             burn = 0L,
                             thin = 1L,
                             metadata = TRUE,
                             include_w = FALSE,
                             ...){

  if(isTRUE(include_w) && !is.null(object$w_samples))
    return(as_samples.stLMM_recovery(
      object,
      burn = burn,
      thin = thin,
      metadata = metadata,
      include_w = include_w,
      ...
    ))

  n <- as.integer(object$backend$n_samples)
  iteration <- sample_draw_index(n, burn = burn, thin = thin)
  mat <- stLMM_parameter_samples(object, index = iteration)

  samples_data_frame(
    mat = mat,
    iteration = iteration,
    chain = object$chain %||% 1L,
    metadata = metadata
  )
}

as_samples.stLMM_recovery <- function(object,
                                      burn = 0L,
                                      thin = 1L,
                                      metadata = TRUE,
                                      include_w = FALSE,
                                      ...){

  if(is.null(object$recover_iter) || !length(object$recover_iter))
    stop("error: recovered object has no recover_iter metadata")

  burn <- validate_burn(burn, object$backend$n_samples)
  keep <- sample_available_index(object$recover_iter, burn = burn, thin = thin,
                                 what = "recovered draws")

  iteration <- as.integer(object$recover_iter[keep])
  mat <- stLMM_parameter_samples(object, index = iteration)

  if(isTRUE(include_w)){
    w <- samples_w_blocks(object$w_samples, index = keep)
    if(!is.null(w))
      mat <- cbind(mat, w)
  }

  samples_data_frame(
    mat = mat,
    iteration = iteration,
    chain = object$chain %||% 1L,
    metadata = metadata
  )
}

as_samples.stLMM_chains <- function(object,
                                    burn = 0L,
                                    thin = 1L,
                                    metadata = TRUE,
                                    include_w = FALSE,
                                    combine_chains = TRUE,
                                    ...){

  out <- lapply(object$chains, as_samples,
                burn = burn, thin = thin, metadata = metadata,
                include_w = include_w, ...)

  if(!isTRUE(combine_chains))
    return(out)

  do.call(rbind, unname(out))
}

as_samples.stLMM_recovery_chains <- function(object,
                                             burn = 0L,
                                             thin = 1L,
                                             metadata = TRUE,
                                             include_w = FALSE,
                                             combine_chains = TRUE,
                                             ...){

  out <- lapply(object$chains, as_samples,
                burn = burn, thin = thin, metadata = metadata,
                include_w = include_w, ...)

  if(!isTRUE(combine_chains))
    return(out)

  do.call(rbind, unname(out))
}

prediction_sample_matrix <- function(object,
                                     sample = c("mu", "y", "all"),
                                     include_w = FALSE){

  sample <- match.arg(sample)
  blocks <- list()

  if(sample %in% c("mu", "all"))
    blocks$mu <- samples_matrix_block(object$mu_samples, "mu")

  if(sample %in% c("y", "all")){
    y <- samples_matrix_block(object$y_samples, "y")
    if(is.null(y) && sample == "y")
      stop("error: y samples are not available")
    if(!is.null(y))
      blocks$y <- y
  }

  if(isTRUE(include_w)){
    w <- samples_w_blocks(object$w_samples)
    if(!is.null(w))
      blocks$w <- w
  }

  blocks <- Filter(Negate(is.null), blocks)
  if(!length(blocks))
    stop("error: no prediction samples available")

  if(!is.null(blocks$mu))
    colnames(blocks$mu) <- paste0("mu_", seq_len(ncol(blocks$mu)))
  if(!is.null(blocks$y))
    colnames(blocks$y) <- paste0("y_", seq_len(ncol(blocks$y)))

  n_row <- vapply(blocks, nrow, integer(1))
  if(length(unique(n_row)) != 1L)
    stop("error: prediction sample blocks do not have a common number of rows")

  do.call(cbind, blocks)
}

as_samples.stLMM_prediction <- function(object,
                                        sample = c("mu", "y", "all"),
                                        burn = 0L,
                                        thin = 1L,
                                        metadata = TRUE,
                                        include_w = FALSE,
                                        ...){

  sample <- match.arg(sample)
  mat <- prediction_sample_matrix(object, sample = sample, include_w = include_w)
  iteration_all <- if(!is.null(object$draw_index)) object$draw_index else seq_len(nrow(mat))
  keep <- sample_available_index(iteration_all, burn = burn, thin = thin,
                                 what = "prediction draws")
  iteration <- iteration_all[keep]

  samples_data_frame(
    mat = mat[keep, , drop = FALSE],
    iteration = iteration,
    chain = 1L,
    metadata = metadata
  )
}

as_samples.stLMM_prediction_chains <- function(object,
                                               sample = c("mu", "y", "all"),
                                               burn = 0L,
                                               thin = 1L,
                                               metadata = TRUE,
                                               include_w = FALSE,
                                               combine_chains = TRUE,
                                               ...){

  sample <- match.arg(sample)
  out <- vector("list", length(object$chains))
  names(out) <- names(object$chains)

  for(i in seq_along(object$chains)){
    out[[i]] <- as_samples(
      object$chains[[i]],
      sample = sample,
      burn = burn,
      thin = thin,
      metadata = metadata,
      include_w = include_w,
      ...
    )
    if(isTRUE(metadata))
      out[[i]]$.chain <- as.integer(i)
  }

  if(!isTRUE(combine_chains))
    return(out)

  do.call(rbind, unname(out))
}

as_samples.matrix <- function(object,
                              burn = 0L,
                              thin = 1L,
                              metadata = TRUE,
                              chain = 1L,
                              ...){

  iteration <- attr(object, "draw_index") %||% seq_len(nrow(object))
  keep <- sample_available_index(iteration, burn = burn, thin = thin,
                                 what = "sample draws")
  iteration <- iteration[keep]

  samples_data_frame(
    mat = object[keep, , drop = FALSE],
    iteration = iteration,
    chain = chain,
    metadata = metadata
  )
}

as_samples.stLMM_fitted_chains <- function(object,
                                           burn = 0L,
                                           thin = 1L,
                                           metadata = TRUE,
                                           combine_chains = TRUE,
                                           ...){

  out <- vector("list", length(object))

  for(i in seq_along(object))
    out[[i]] <- as_samples.matrix(object[[i]], burn = burn, thin = thin,
                                  metadata = metadata, chain = i, ...)

  if(!isTRUE(combine_chains))
    return(out)

  do.call(rbind, out)
}
