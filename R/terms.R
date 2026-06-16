############################################################
# Term constructors
############################################################

build_obs_index <- function(map){

  map <- as.integer(map)

  if(length(map) == 0L)
    stop("error: map must have positive length")

  if(anyNA(map))
    stop("error: map contains NA values")

  n_obs <- length(map)
  n_node <- max(map)

  if(!setequal(unique(map), seq_len(n_node)))
    stop("error: map must contain contiguous positive node ids starting at 1")

  ord <- order(map)
  map_ord <- map[ord]

  obsIndx <- as.integer(ord)
  obsIndxLU <- matrix(0L, nrow = n_node, ncol = 2L)
  node_nobs <- integer(n_node)

  pos <- 1L

  for(i in seq_len(n_node)){

    start <- pos

    while(pos <= n_obs && map_ord[pos] == i)
      pos <- pos + 1L

    end <- pos - 1L

    obsIndxLU[i, ] <- c(start, end)
    node_nobs[i] <- end - start + 1L
  }

  if(sum(node_nobs) != n_obs)
    stop("error: compressed observation index does not sum to n_obs")

  list(
    obsIndx = obsIndx,
    obsIndxLU = obsIndxLU,
    node_nobs = as.integer(node_nobs),
    n_node = as.integer(n_node)
  )
}

build_obs_index_with_n_node <- function(map, n_node){

  map <- as.integer(map)
  n_node <- as.integer(n_node)

  if(length(n_node) != 1L || is.na(n_node) || n_node < 1L)
    stop("error: n_node must be a positive integer")

  if(anyNA(map))
    stop("error: map contains NA values")

  n_obs <- length(map)

  if(n_obs && (any(map < 1L) || any(map > n_node)))
    stop("error: map contains node ids outside 1:n_node")

  obsIndx <- integer(n_obs)
  obsIndxLU <- matrix(0L, nrow = n_node, ncol = 2L)
  node_nobs <- integer(n_node)

  if(n_obs){
    ord <- order(map)
    map_ord <- map[ord]
    obsIndx <- as.integer(ord)

    pos <- 1L
    for(i in seq_len(n_node)){
      start <- pos
      while(pos <= n_obs && map_ord[pos] == i)
        pos <- pos + 1L
      end <- pos - 1L

      if(end >= start){
        obsIndxLU[i, ] <- c(start, end)
        node_nobs[i] <- end - start + 1L
      }
    }
  }

  if(sum(node_nobs) != n_obs)
    stop("error: compressed observation index does not sum to n_obs")

  list(
    obsIndx = as.integer(obsIndx),
    obsIndxLU = obsIndxLU,
    node_nobs = as.integer(node_nobs),
    n_node = as.integer(n_node)
  )
}

subset_process_term_observed <- function(term, observed_index){

  out <- term
  observed_index <- as.integer(observed_index)

  out$map <- as.integer(term$map[observed_index])
  if(!is.null(term$x))
    out$x <- as.numeric(term$x[observed_index])

  obs_info <- build_obs_index_with_n_node(out$map, term$n_node)
  out$n_obs <- as.integer(length(out$map))
  out$obsIndx <- obs_info$obsIndx
  out$obsIndxLU <- obs_info$obsIndxLU
  out$node_nobs <- obs_info$node_nobs
  out$n_node <- obs_info$n_node

  out
}

make_process_term <- function(term_type,
                              graph_id,
                              label,
                              name = NULL,
                              map,
                              coef_name = NULL,
                              cov_model = NULL,
                              params = list(),
                              data = NULL,
                              n_node = NULL){

  if(anyNA(map))
    stop("error: NA values detected in term map for ", label)

  if(!is.list(params))
    stop("error: params must be a list for ", label)

  ## params stores structure / correlation settings for the process term.
  ## Variance scaling is handled separately from the graph/operator layer.
  map <- as.integer(map)
  n_obs <- length(map)

  if(is.null(n_node)){
    obs_info <- build_obs_index(map)
  } else {
    obs_info <- build_obs_index_with_n_node(map, n_node)
  }

  # ---- attach covariate for SVC terms (applies to ar1, nngp, etc.)
  x <- NULL

  if(!is.null(coef_name)){
    if(is.null(data))
      stop("internal error: data must be supplied when coef_name is used")

    if(!coef_name %in% names(data))
      stop("variable '", coef_name, "' not found in data")

    x <- as.numeric(data[[coef_name]])

    if(length(x) != n_obs)
      stop("length mismatch for covariate '", coef_name, "' in term ", label)
  }

  out <- list(
    term_type = term_type,
    graph_id = graph_id,
    label = label,
    name = if(!is.null(name)) name else label,
    map = map,
    coef_name = coef_name,
    x = x,
    is.svc = !is.null(coef_name),
    cov_model = cov_model,
    params = params,
    n_obs = as.integer(n_obs),
    n_node = obs_info$n_node,
    obsIndx = obs_info$obsIndx,
    obsIndxLU = obs_info$obsIndxLU,
    node_nobs = obs_info$node_nobs
  )

  class(out) <- "stLMM_process_term"
  out
}
