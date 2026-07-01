# Utility functions for the stLMM vignettes.
#
# These helpers are intentionally small and explicit. They build covariance and
# precision matrices using the same forms as the stLMM model terms. The
# vignettes keep the actual random-effect and response construction visible.

stlmm_palette <- function(n = 200, palette = "Spectral", reverse = TRUE){
  colors <- grDevices::hcl.colors(n, palette = palette)

  if(reverse)
    colors <- rev(colors)

  colors
}

stlmm_continuous_palette <- function(n = 200, palette = c("navia", "batlowW"), reverse = FALSE){
  palette <- match.arg(palette)

  anchors <- switch(
    palette,
    navia = c(
      "#021326", "#053059", "#105185", "#236B90", "#327B88", "#408A7F",
      "#539A74", "#70B369", "#A6D278", "#DCE6AD", "#FCF3D8"
    ),
    batlowW = c(
      "#001959", "#0F3C5F", "#175361", "#2D685D", "#52784C", "#7C8637",
      "#B19939", "#DAA66B", "#F2B39E", "#FEDBDA", "#FFFEFE"
    )
  )

  colors <- grDevices::colorRampPalette(anchors, space = "Lab")(n)

  if(reverse)
    colors <- rev(colors)

  colors
}

stlmm_discrete_colors <- function(n, palette = "Spectral", reverse = TRUE){
  if(!is.numeric(n) || length(n) != 1 || n < 1)
    stop("n must be a positive scalar.")

  colors <- stlmm_palette(200, palette = palette, reverse = reverse)
  idx <- switch(
    as.character(as.integer(n)),
    "1" = 25,
    "2" = c(25, 175),
    "3" = c(25, 65, 175),
    "4" = c(25, 65, 145, 180),
    "5" = c(20, 55, 90, 145, 180),
    round(seq(20, 180, length.out = n))
  )

  colors[idx]
}

stlmm_color <- function(name = c("primary", "secondary", "accent", "contrast", "muted")){
  name <- match.arg(name)
  colors <- stlmm_palette(200)
  named <- c(
    primary = colors[25],
    secondary = colors[175],
    accent = colors[65],
    contrast = colors[145],
    muted = colors[90]
  )

  unname(named[name])
}

rmvnorm <- function(mean, Sigma, n = 1, jitter = 0, drop = TRUE){
  mean <- as.numeric(mean)
  Sigma <- as.matrix(Sigma)

  if(nrow(Sigma) != ncol(Sigma))
    stop("Sigma must be a square matrix.")
  if(length(mean) != nrow(Sigma))
    stop("mean and Sigma have incompatible dimensions.")
  if(!is.numeric(n) || length(n) != 1 || n < 1)
    stop("n must be a positive scalar.")

  if(jitter > 0)
    diag(Sigma) <- diag(Sigma) + jitter

  z <- matrix(rnorm(n * length(mean)), nrow = n)
  draws <- sweep(z %*% chol(Sigma), 2, mean, "+")

  if(drop && n == 1)
    as.numeric(draws[1, ])
  else
    draws
}

dist_mat <- function(coords){
  coords <- as.matrix(coords)

  as.matrix(stats::dist(coords))
}

exp_cor <- function(coords, phi){
  if(phi <= 0)
    stop("phi must be positive.")

  exp(-phi * dist_mat(coords))
}

exp_cov <- function(coords, sigma_sq, phi){
  if(sigma_sq <= 0)
    stop("sigma_sq must be positive.")

  C <- sigma_sq * exp_cor(coords, phi = phi)

  C
}

sep_exp_cov <- function(coords, time, sigma_sq, phi, lambda){
  if(sigma_sq <= 0)
    stop("sigma_sq must be positive.")
  if(phi <= 0)
    stop("phi must be positive.")
  if(lambda <= 0)
    stop("lambda must be positive.")

  coords <- as.matrix(coords)
  time <- as.numeric(time)

  if(nrow(coords) != length(time))
    stop("coords and time have incompatible dimensions.")

  space_dist <- dist_mat(coords)
  time_dist <- abs(outer(time, time, "-"))

  sigma_sq * exp(-phi * space_dist - lambda * time_dist)
}

gneiting_cov <- function(coords, time, sigma_sq, a, c, alpha, beta, gamma, delta){
  if(sigma_sq <= 0)
    stop("sigma_sq must be positive.")
  if(a <= 0)
    stop("a must be positive.")
  if(c <= 0)
    stop("c must be positive.")
  if(alpha <= 0 || alpha > 1)
    stop("alpha must be in (0, 1].")
  if(beta < 0 || beta > 1)
    stop("beta must be between 0 and 1.")
  if(gamma <= 0 || gamma > 1)
    stop("gamma must be in (0, 1].")
  if(delta < 0)
    stop("delta must be nonnegative.")

  coords <- as.matrix(coords)
  time <- as.numeric(time)

  if(nrow(coords) != length(time))
    stop("coords and time have incompatible dimensions.")

  space_dist <- dist_mat(coords)
  time_dist <- abs(outer(time, time, "-"))
  scale <- 1 + a * time_dist^(2 * alpha)
  spatial_dim <- ncol(coords)

  sigma_sq * scale^(-(delta + spatial_dim / 2)) *
    exp(-c * space_dist^(2 * gamma) / scale^(beta * gamma))
}

matern_cor <- function(coords, phi, nu){
  if(phi <= 0)
    stop("phi must be positive.")
  if(nu <= 0)
    stop("nu must be positive.")

  h <- dist_mat(coords)
  x <- phi * h
  R <- matrix(1, nrow(h), ncol(h))
  positive <- x > 0

  R[positive] <- exp(
    nu * log(x[positive]) -
      (nu - 1) * log(2) -
      lgamma(nu) -
      x[positive] +
      log(besselK(x[positive], nu, expon.scaled = TRUE))
  )

  R[!is.finite(R)] <- 0
  R
}

matern_cov <- function(coords, sigma_sq, phi, nu){
  if(sigma_sq <= 0)
    stop("sigma_sq must be positive.")

  C <- sigma_sq * matern_cor(coords, phi = phi, nu = nu)

  C
}

ar1_cor <- function(index, phi){
  if(abs(phi) >= 1)
    stop("phi must be between -1 and 1.")

  support <- sort(unique(index))
  lag <- abs(outer(seq_along(support), seq_along(support), "-"))
  phi^lag
}

ar1_prec <- function(index, sigma_sq, phi){
  if(sigma_sq <= 0)
    stop("sigma_sq must be positive.")
  if(abs(phi) >= 1)
    stop("phi must be between -1 and 1.")

  support <- sort(unique(index))
  q <- length(support)
  Q <- matrix(0, q, q, dimnames = list(support, support))

  if(q == 1){
    Q[1, 1] <- 1 / sigma_sq
    return(Q)
  }

  den <- 1 - phi^2
  diag(Q) <- c(1, rep(1 + phi^2, q - 2), 1) / den
  off <- -phi / den
  Q[cbind(seq_len(q - 1), 2:q)] <- off
  Q[cbind(2:q, seq_len(q - 1))] <- off

  Q / sigma_sq
}

exp_time_cor <- function(time, lambda){
  if(lambda <= 0)
    stop("lambda must be positive.")

  support <- sort(unique(time))
  exp(-lambda * abs(outer(support, support, "-")))
}

exp_time_prec <- function(time, sigma_sq, lambda){
  if(sigma_sq <= 0)
    stop("sigma_sq must be positive.")
  if(lambda <= 0)
    stop("lambda must be positive.")

  support <- sort(unique(time))
  q <- length(support)
  Q <- matrix(0, q, q, dimnames = list(support, support))

  if(q == 1){
    Q[1, 1] <- 1 / sigma_sq
    return(Q)
  }

  gap_phi <- exp(-lambda * diff(support))
  den <- 1 - gap_phi^2
  diag(Q)[1] <- 1 / den[1]
  diag(Q)[q] <- 1 / den[q - 1]

  if(q > 2){
    for(i in 2:(q - 1))
      diag(Q)[i] <- 1 / den[i - 1] + gap_phi[i]^2 / den[i]
  }

  off <- -gap_phi / den
  Q[cbind(seq_len(q - 1), 2:q)] <- off
  Q[cbind(2:q, seq_len(q - 1))] <- off

  Q / sigma_sq
}

car_prec <- function(adjacency, sigma_sq, rho){
  if(inherits(adjacency, "stLMM_car_graph"))
    adjacency <- adjacency$adjacency

  adjacency <- as.matrix(adjacency)

  if(nrow(adjacency) != ncol(adjacency))
    stop("adjacency must be a square matrix.")
  if(sigma_sq <= 0)
    stop("sigma_sq must be positive.")
  if(rho <= 0 || rho >= 1)
    stop("rho must be between 0 and 1.")

  degree <- rowSums(adjacency)
  (diag(degree, nrow(adjacency)) - rho * adjacency) / sigma_sq
}

dagar_prec <- function(graph, sigma_sq, rho, ordering = "coord"){
  if(!inherits(graph, "stLMM_car_graph"))
    stop("graph must be a stLMM_car_graph object.")
  if(sigma_sq <= 0)
    stop("sigma_sq must be positive.")
  if(rho < 0 || rho >= 1)
    stop("rho must be in [0, 1).")

  n <- graph$n
  if(is.numeric(ordering)){
    ord <- as.integer(ordering)
    if(length(ord) != n || any(sort(ord) != seq_len(n)))
      stop("numeric ordering must be a permutation of 1:n.")
    ordering_type <- "user"
  } else {
    ordering <- match.arg(ordering, c("coord", "default"))
    ordering_type <- ordering
    if(is.null(graph$geometry)){
      ord <- seq_len(n)
    } else {
      pts <- suppressWarnings(sf::st_point_on_surface(graph$geometry))
      coords <- sf::st_coordinates(pts)
      ord <- do.call(order, c(as.data.frame(coords[, 1:2, drop = FALSE]), list(decreasing = FALSE)))
    }
  }

  ord_inv <- integer(n)
  ord_inv[ord] <- seq_len(n)
  adj <- as.matrix(graph$adjacency)
  parent_lists <- vector("list", n)
  for(i in seq_len(n))
    parent_lists[[i]] <- integer(0)

  edges <- which(adj != 0, arr.ind = TRUE)
  edges <- edges[edges[, "row"] < edges[, "col"], , drop = FALSE]
  for(e in seq_len(nrow(edges))){
    i_support <- edges[e, "row"]
    j_support <- edges[e, "col"]
    i_ord <- ord_inv[i_support]
    j_ord <- ord_inv[j_support]
    if(i_ord < j_ord)
      parent_lists[[j_ord]] <- c(parent_lists[[j_ord]], i_ord)
    else
      parent_lists[[i_ord]] <- c(parent_lists[[i_ord]], j_ord)
  }

  parent_count <- vapply(parent_lists, length, integer(1))
  parent_start <- integer(n)
  if(n > 1L)
    parent_start[-1L] <- cumsum(parent_count)[-n]
  total_parent <- sum(parent_count)
  parent_index <- integer(total_parent)
  pos <- 1L
  for(i in seq_len(n)){
    if(parent_count[i] > 0L){
      parents <- sort(unique(parent_lists[[i]]))
      parent_index[pos:(pos + parent_count[i] - 1L)] <- parents
      pos <- pos + parent_count[i]
    }
  }

  dg <- list(
    graph_type = "dagar",
    n = as.integer(n),
    q = as.integer(n),
    ids = graph$ids[ord],
    ids_support = graph$ids,
    geometry = if(!is.null(graph$geometry)) graph$geometry[ord] else NULL,
    geometry_support = graph$geometry,
    degree = as.double(graph$degree[ord]),
    n_edge = as.integer(nrow(edges)),
    ordering = ordering_type,
    ordering_type = ordering_type,
    ord = as.integer(ord),
    ord_inv = as.integer(ord_inv),
    parent_index = as.integer(parent_index - 1L),
    parent_start = as.integer(parent_start),
    parent_count = as.integer(parent_count),
    zero_parent_nodes = as.integer(which(parent_count == 0L))
  )
  class(dg) <- "stLMM_graph"

  Q <- matrix(0, n, n, dimnames = list(dg$ids, dg$ids))

  for(i in seq_len(n)){
    m <- dg$parent_count[i]
    parents <- if(m > 0L) {
      dg$parent_index[seq.int(dg$parent_start[i] + 1L, length.out = m)] + 1L
    } else {
      integer(0)
    }

    denom <- 1 + (m - 1) * rho^2
    b <- if(m > 0L) rho / denom else 0
    f <- denom / (1 - rho^2)
    ell <- numeric(n)
    ell[i] <- 1
    if(m > 0L)
      ell[parents] <- -b

    Q <- Q + f * tcrossprod(ell)
  }

  Q <- Q / sigma_sq
  attr(Q, "dagar_graph") <- dg
  Q
}

fitted_draw_matrix <- function(object,
                               sub_sample = list(start = 1L, thin = 1L),
                               scale = "response"){
  if(is.list(object) && !is.null(object$chains)){
    x <- lapply(object$chains, stats::fitted, summary = FALSE,
                sub_sample = sub_sample, scale = scale)
    return(do.call(rbind, x))
  }

  x <- stats::fitted(object, summary = FALSE,
                     sub_sample = sub_sample, scale = scale)
  if(inherits(x, "stLMM_fitted_chains"))
    x <- do.call(rbind, unclass(x))

  x
}
