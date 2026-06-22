
############################################################
# Utilities
############################################################

`%||%` <- function(x, y) if(is.null(x)) y else x

is_sf_geometry_column <- function(x){
  inherits(x, "sfc") || inherits(x, "sfc_POINT")
}

extract_process_coords <- function(data, coord_names){

  if(length(coord_names) == 1L){

    nm <- coord_names[1]

    if(!nm %in% names(data))
      stop("error: ", nm, " not found in data")

    col <- data[[nm]]

    ## Case 1: sf geometry column
    if(is_sf_geometry_column(col)){

      if(!requireNamespace("sf", quietly = TRUE))
        stop("error: package 'sf' is required for geometry-based process terms")

      coords <- sf::st_coordinates(col)

      if(ncol(coords) < 2L)
        stop("error: geometry column must contain point coordinates")

      coords <- coords[,1:2,drop=FALSE]
      colnames(coords) <- c("x","y")

    } else {

      ## Case 2: numeric 1-D coordinate
      if(!is.numeric(col))
        stop("error: coordinate column '", nm, "' must be numeric or an sf geometry")

      coords <- matrix(col, ncol = 1)
      colnames(coords) <- nm
    }

  } else {

    ## Case 3: multi-dimensional numeric coordinates
    miss <- setdiff(coord_names, names(data))
    if(length(miss))
      stop("error: missing coordinate column(s): ", paste(miss, collapse = ", "))

    coords <- as.matrix(data[, coord_names, drop = FALSE])

    if(!all(vapply(coords, is.numeric, logical(1))))
      stop("error: coordinate columns must be numeric")
  }

  storage.mode(coords) <- "double"
  if(any(!is.finite(coords))){
    bad <- which(rowSums(!is.finite(coords)) > 0L)
    stop(
      "error: process coordinates must be finite; found ",
      length(bad),
      " row(s) with NA, NaN, or Inf in ",
      paste(colnames(coords) %||% paste0("coord", seq_len(ncol(coords))), collapse = ", "),
      ". Remove or impute those rows before fitting."
    )
  }
  coords
}


car_graph_signature <- function(graph){

  if(inherits(graph, "stLMM_car_graph")){
    adj <- graph$adjacency
  } else if(inherits(graph, "Matrix") || is.matrix(graph)) {
    adj <- Matrix::Matrix(graph, sparse = TRUE)
  } else {
    return(NULL)
  }

  ids <- rownames(adj)
  if(is.null(ids))
    ids <- colnames(adj)
  if(is.null(ids))
    ids <- as.character(seq_len(nrow(adj)))

  dimnames(adj) <- list(as.character(ids), as.character(ids))
  adj <- Matrix::drop0(adj)
  trip <- Matrix::summary(adj)
  keep <- trip$i < trip$j
  trip <- trip[keep, , drop = FALSE]

  if(nrow(trip)){
    ord <- order(trip$i, trip$j)
    trip <- trip[ord, , drop = FALSE]
    edge_tag <- paste(
      paste(ids[trip$i], ids[trip$j], signif(trip$x, 16), sep = ":"),
      collapse = ";"
    )
  } else {
    edge_tag <- ""
  }

  paste0(
    "ids=", paste(ids, collapse = ";"),
    "|edges=", edge_tag
  )
}

car_graph_key_tag <- function(graph_obj, graph_id_col = NULL, queen = TRUE){

  sig <- car_graph_signature(graph_obj)
  if(!is.null(sig))
    return(paste0("car_graph:", sig))

  if(inherits(graph_obj, "sf")){
    if(is.null(graph_id_col) || !is.character(graph_id_col) ||
       length(graph_id_col) != 1L || !nzchar(graph_id_col)){
      return(paste(deparse(graph_obj, nlines = 1L), collapse = ""))
    }

    graph <- car_graph(graph_obj, id = graph_id_col, queen = queen)
    sig <- car_graph_signature(graph)
    return(paste0("car_graph:", sig))
  }

  if(is.character(graph_obj) && length(graph_obj) == 1L)
    return(graph_obj)

  paste(deparse(graph_obj, nlines = 1L), collapse = "")
}

normalize_car_time_model <- function(x, env = NULL){

  if(is.character(x) && length(x) == 1L){
    val <- tolower(gsub('^["\']|["\']$', "", x))
    if(val %in% c("ar1", "exp"))
      return(val)
  }

  if(!is.null(env))
    x <- eval_process_param(x, env)

  if(identical(x, base::exp))
    return("exp")

  if(is.character(x) && length(x) == 1L){
    val <- tolower(gsub('^["\']|["\']$', "", x))
    if(val %in% c("ar1", "exp"))
      return(val)
  }

  x
}

normalize_car_model <- function(x, env = NULL){

  if(is.character(x) && length(x) == 1L){
    val <- tolower(gsub('^["\']|["\']$', "", x))
    if(val %in% c("proper", "leroux"))
      return(val)
  }

  if(!is.null(env))
    x <- eval_process_param(x, env)

  if(is.character(x) && length(x) == 1L){
    val <- tolower(gsub('^["\']|["\']$', "", x))
    if(val %in% c("proper", "leroux"))
      return(val)
  }

  x
}

build_graph_key <- function(spec, defaults = list(), formula_env = NULL){

  fun <- spec$fun

  ## Coordinate order is part of the graph definition. In particular,
  ## space-time covariance models treat the final coordinate as time.
  args <- trimws(spec$args)
  args_str <- paste(args, collapse = ",")

  if(identical(fun, "car")){
    graph_obj <- spec$params$graph %||% "graph"
    graph_id <- spec$params$graph_id %||% spec$params$id %||% ""
    queen <- spec$params$queen %||% TRUE
    if(!is.null(formula_env)){
      graph_obj <- eval_process_param(graph_obj, formula_env)
      graph_id <- eval_process_param(graph_id, formula_env)
      queen <- eval_process_param(queen, formula_env)
    }
    graph_tag <- car_graph_key_tag(graph_obj, graph_id_col = graph_id, queen = queen)
    if(length(graph_id) != 1L)
      graph_id <- paste(deparse(graph_id, nlines = 1L), collapse = "")
    car_model <- normalize_car_model(spec$params$car_model %||% "proper", formula_env)
    return(paste(fun, args_str, graph_tag, graph_id, car_model, sep = "::"))
  }

  if(identical(fun, "car_time")){
    graph_obj <- spec$params$graph %||% "graph"
    graph_id <- spec$params$graph_id %||% spec$params$id %||% ""
    queen <- spec$params$queen %||% TRUE
    if(!is.null(formula_env)){
      graph_obj <- eval_process_param(graph_obj, formula_env)
      graph_id <- eval_process_param(graph_id, formula_env)
      queen <- eval_process_param(queen, formula_env)
    }
    graph_tag <- car_graph_key_tag(graph_obj, graph_id_col = graph_id, queen = queen)
    if(length(graph_id) != 1L)
      graph_id <- paste(deparse(graph_id, nlines = 1L), collapse = "")
    time_model <- normalize_car_time_model(spec$params$time_model %||% "ar1", formula_env)
    car_model <- normalize_car_model(spec$params$car_model %||% "proper", formula_env)
    return(paste(fun, args_str, graph_tag, graph_id, time_model, car_model, sep = "::"))
  }

  if(identical(fun, "dagar")){
    graph_obj <- spec$params$graph %||% "graph"
    graph_id <- spec$params$graph_id %||% spec$params$id %||% ""
    queen <- spec$params$queen %||% TRUE
    ordering <- spec$params$ordering %||% "coord"
    if(!is.null(formula_env)){
      graph_obj <- eval_process_param(graph_obj, formula_env)
      graph_id <- eval_process_param(graph_id, formula_env)
      queen <- eval_process_param(queen, formula_env)
      ordering <- eval_process_param(ordering, formula_env)
    }
    graph_tag <- car_graph_key_tag(graph_obj, graph_id_col = graph_id, queen = queen)
    if(length(graph_id) != 1L)
      graph_id <- paste(deparse(graph_id, nlines = 1L), collapse = "")
    if(is.numeric(ordering)){
      ordering_tag <- paste0("ordering=user[", length(ordering), "]=", paste(as.integer(ordering), collapse = ";"))
    } else {
      ordering_tag <- paste0("ordering=", paste(deparse(ordering, nlines = 1L), collapse = ""))
    }
    return(paste(fun, args_str, graph_tag, graph_id, ordering_tag, sep = "::"))
  }

  if(identical(fun, "dagar_time")){
    graph_obj <- spec$params$graph %||% "graph"
    graph_id <- spec$params$graph_id %||% spec$params$id %||% ""
    queen <- spec$params$queen %||% TRUE
    ordering <- spec$params$ordering %||% "coord"
    if(!is.null(formula_env)){
      graph_obj <- eval_process_param(graph_obj, formula_env)
      graph_id <- eval_process_param(graph_id, formula_env)
      queen <- eval_process_param(queen, formula_env)
      ordering <- eval_process_param(ordering, formula_env)
    }
    graph_tag <- car_graph_key_tag(graph_obj, graph_id_col = graph_id, queen = queen)
    if(length(graph_id) != 1L)
      graph_id <- paste(deparse(graph_id, nlines = 1L), collapse = "")
    if(is.numeric(ordering)){
      ordering_tag <- paste0("ordering=user[", length(ordering), "]=", paste(as.integer(ordering), collapse = ";"))
    } else {
      ordering_tag <- paste0("ordering=", paste(deparse(ordering, nlines = 1L), collapse = ""))
    }
    time_model <- normalize_car_time_model(spec$params$time_model %||% "ar1", formula_env)
    return(paste(fun, args_str, graph_tag, graph_id, ordering_tag, time_model, sep = "::"))
  }

  ## grouping does not affect graph structure for currently implemented
  ## process graphs; it only affects term-level mapping / replication.
  group_str <- "_obs"

  ## merge defaults with supplied params
  params <- modifyList(defaults, spec$params)

  if(identical(fun, "nngp")){
    cov_model <- params$cov_model %||% "exp"
    registry <- build_cor_model_registry()
    if(!is.null(formula_env) &&
       !(is.character(cov_model) && length(cov_model) == 1L && cov_model %in% names(registry)))
      cov_model <- eval_process_param(cov_model, formula_env)
    if(!is.character(cov_model) || length(cov_model) != 1L || !cov_model %in% names(registry))
      stop("error: nngp(): unknown cov_model '", paste(deparse(cov_model, nlines = 1L), collapse = ""), "'")
    params$space_time <- identical(as.integer(registry[[cov_model]]$distance_mode), 2L)
  }

  ## cov_model does NOT otherwise affect graph structure
  params$cov_model <- NULL

  if(length(params) == 0L){

    param_str <- ""

  } else {

    nm <- sort(names(params))

    vals <- vapply(
      nm,
      function(name){

        val <- params[[name]]

        ## canonicalize ordering vectors
        if(identical(name, "ordering") && is.numeric(val)){

          val <- as.integer(val)

          paste0(
            "ordering=user[",
            length(val),
            "]=",
            paste(val, collapse=";")
          )

        } else if(identical(name, "st_scale") && is.numeric(val) && length(val) == 1L) {

          paste0(name, "=", formatC(as.numeric(val), digits = 17L, format = "fg", flag = "#"))

        } else {

          paste0(name,"=",paste(deparse(val,nlines=1),collapse=""))

        }

      },
      character(1)
    )

    param_str <- paste(vals, collapse = ",")
  }

  paste(fun, args_str, group_str, param_str, sep = "::")
}

coordinate_row_keys <- function(x){
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  x[x == 0] <- 0
  x_chr <- matrix(NA_character_, nrow = nrow(x), ncol = ncol(x))
  for(j in seq_len(ncol(x))){
    x_chr[, j] <- formatC(x[, j], digits = 17L, format = "fg", flag = "#")
  }
  do.call(paste, c(as.data.frame(x_chr, stringsAsFactors = FALSE), sep = "\r"))
}

unique_coordinate_rows <- function(x){
  keys <- coordinate_row_keys(x)
  x[!duplicated(keys), , drop = FALSE]
}

match_unique_rows <- function(x, u){
  match(coordinate_row_keys(x), coordinate_row_keys(u))
}

############################################################
# Ordering helpers
############################################################

validate_ordering_vector <- function(ord, q, where = "ordering"){
  if(!is.numeric(ord) || length(ord) != q)
    stop("error: ", where, " must be a numeric permutation of 1:q with q = ", q)

  ord <- as.integer(ord)

  if(anyNA(ord))
    stop("error: ", where, " contains NA values")

  if(!setequal(ord, seq_len(q)))
    stop("error: ", where, " must be a permutation of 1:q")

  ord
}

maxmin_order <- function(coords){

  q <- nrow(coords)

  if(q == 1L)
    return(1L)

  ord <- integer(q)

  ord[1] <- do.call(order, c(as.data.frame(coords), list(decreasing = FALSE)))[1]
  chosen <- logical(q)
  chosen[ord[1]] <- TRUE

  dmin <- rep(Inf, q)

  for(i in 2:q){
    last <- ord[i - 1L]
    d2 <- rowSums((coords - matrix(coords[last, ], q, ncol(coords), byrow = TRUE))^2)
    dmin <- pmin(dmin, d2)
    dmin[chosen] <- -Inf
    ord[i] <- which.max(dmin)
    chosen[ord[i]] <- TRUE
  }

  ord
}

.hilbert_rot <- function(n, x, y, rx, ry){
  if(ry == 0L){
    if(rx == 1L){
      x <- n - 1L - x
      y <- n - 1L - y
    }
    tmp <- x
    x <- y
    y <- tmp
  }
  list(x = x, y = y)
}

.hilbert_xy2d <- function(bits, x, y){
  n <- bitwShiftL(1L, bits)
  d <- 0
  s <- bitwShiftL(1L, bits - 1L)

  while(s > 0L){
    rx <- if(bitwAnd(x, s) > 0L) 1L else 0L
    ry <- if(bitwAnd(y, s) > 0L) 1L else 0L
    d <- d + s * s * bitwXor(3L * rx, ry)
    rot <- .hilbert_rot(s, x, y, rx, ry)
    x <- rot$x
    y <- rot$y
    s <- bitwShiftR(s, 1L)
  }

  d
}

hilbert_order <- function(coords, bits = 10L){

  if(ncol(coords) < 2L)
    stop("error: hilbert ordering currently requires at least two coordinates")

  x <- coords[, 1L]
  y <- coords[, 2L]

  scale01 <- function(z){
    rng <- range(z)
    if(diff(rng) == 0) return(rep(0, length(z)))
    (z - rng[1L]) / diff(rng)
  }

  maxv <- bitwShiftL(1L, bits) - 1L
  xi <- as.integer(round(scale01(x) * maxv))
  yi <- as.integer(round(scale01(y) * maxv))

  h <- vapply(seq_along(xi), function(i) .hilbert_xy2d(bits, xi[i], yi[i]), numeric(1))
  order(h, seq_along(h))
}

compute_nngp_order <- function(coords, ordering){

  q <- nrow(coords)

  if(q < 1L)
    stop("error: empty coordinate matrix")

  if(is.numeric(ordering)){
    ord <- validate_ordering_vector(ordering, q)
    ordering_type <- "user"
  } else {

    if(!is.character(ordering) || length(ordering) != 1L || !nzchar(ordering))
      stop("error: ordering must be one of 'coord', 'default', 'maxmin', 'hilbert', 'random' or a numeric permutation")

    ordering_type <- ordering

    ord <- switch(
      ordering,
      coord = do.call(order, c(as.data.frame(coords), list(decreasing = FALSE))),
      default = do.call(order, c(as.data.frame(coords), list(decreasing = FALSE))),
      maxmin = maxmin_order(coords),
      hilbert = hilbert_order(coords),
      random = sample.int(q),
      stop("error: unknown ordering '", ordering, "'")
    )
  }

  ord <- as.integer(ord)
  ord_inv <- integer(q)
  ord_inv[ord] <- seq_len(q)

  list(
    ord = ord,
    ord_inv = ord_inv,
    ordering_type = ordering_type
  )
}

compute_graph_order <- function(graph, ordering){

  q <- as.integer(graph$n)
  if(q < 1L)
    stop("error: empty graph")

  if(is.numeric(ordering)){
    ord <- validate_ordering_vector(ordering, q)
    ordering_type <- "user"
  } else {
    if(!is.character(ordering) || length(ordering) != 1L || !nzchar(ordering))
      stop("error: ordering must be one of 'coord', 'maxmin', 'hilbert', 'random' or a numeric permutation")

    ordering <- tolower(ordering)
    ok <- c("coord", "default", "maxmin", "hilbert", "random")
    if(!ordering %in% ok)
      stop("error: ordering must be one of ", paste(shQuote(ok), collapse = ", "), " or a numeric permutation")

    ordering_type <- ordering
    if(ordering %in% c("coord", "default", "maxmin", "hilbert")){
      if(is.null(graph$geometry)){
        if(ordering %in% c("coord", "default")){
          ord <- seq_len(q)
        } else {
          stop("error: ", ordering, " ordering for areal graph terms requires an sf-based car_graph() object")
        }
      } else {
        pts <- suppressWarnings(sf::st_point_on_surface(graph$geometry))
        coords <- sf::st_coordinates(pts)
        coords <- as.matrix(coords[, 1:2, drop = FALSE])
        storage.mode(coords) <- "double"
        ord <- switch(
          ordering,
          coord = do.call(order, c(as.data.frame(coords), list(decreasing = FALSE))),
          default = do.call(order, c(as.data.frame(coords), list(decreasing = FALSE))),
          maxmin = maxmin_order(coords),
          hilbert = hilbert_order(coords)
        )
      }
    } else {
      ord <- sample.int(q)
    }
  }

  ord <- as.integer(ord)
  ord_inv <- integer(q)
  ord_inv[ord] <- seq_len(q)

  list(
    ord = ord,
    ord_inv = ord_inv,
    ordering_type = ordering_type
  )
}

scale_nngp_search_coords <- function(coords, st_scale = 1, space_time = FALSE){

  coords <- as.matrix(coords)
  storage.mode(coords) <- "double"

  if(!is.numeric(st_scale) || length(st_scale) != 1L ||
     is.na(st_scale) || st_scale <= 0)
    stop("error: st_scale must be a positive scalar")

  st_scale <- as.numeric(st_scale)
  if(!isTRUE(space_time) && !isTRUE(all.equal(st_scale, 1)))
    stop("error: st_scale is only valid for space-time covariance models")

  if(isTRUE(space_time))
    coords[, ncol(coords)] <- coords[, ncol(coords)] * st_scale

  coords
}

############################################################
# Graph constructors
############################################################

make_nngp_graph <- function(coords,
                            m = 15,
                            ordering = "coord",
                            st_scale = 1,
                            space_time = FALSE,
                            n_omp_threads = 1,
                            nngp_search = c("fast", "brute")){

  if(!is.matrix(coords))
    stop("error: coords must be a matrix")

  storage.mode(coords) <- "double"
  if(any(!is.finite(coords)))
    stop("error: NNGP graph coordinates must be finite")

  q <- nrow(coords)
  d <- ncol(coords)

  if(q == 0L)
    stop("error: empty coordinate matrix")

  if(!is.numeric(m) || length(m) != 1L || is.na(m) || m < 1L || abs(m - round(m)) > 0)
    stop("error: m must be a positive integer")

  m <- as.integer(m)
  space_time <- isTRUE(space_time)
  nngp_search <- match.arg(nngp_search)

  if(m >= q)
    stop(
      "error: m (", m,
      ") must be less than the number of spatial nodes (", q,
      ") for an NNGP graph"
    )

  if(space_time && identical(ordering, "hilbert"))
    stop("error: hilbert ordering is not supported for space-time NNGP graphs")

  search_coords <- scale_nngp_search_coords(
    coords,
    st_scale = st_scale,
    space_time = space_time
  )
  st_scale <- as.numeric(st_scale)

  ord_info <- compute_nngp_order(search_coords, ordering)
  coords_ord <- coords[ord_info$ord, , drop = FALSE]
  search_coords_ord <- search_coords[ord_info$ord, , drop = FALSE]

  indx <- switch(
    nngp_search,
    fast = mkNNIndx(search_coords_ord, m, n_omp_threads),
    brute = mkNNIndxBrute(search_coords_ord, m, n_omp_threads)
  )

  nn_indx <- indx$nnIndx
  nn_indx_lu <- indx$nnIndxLU

  storage.mode(nn_indx) <- "integer"
  storage.mode(nn_indx_lu) <- "integer"

  out <- list(
    graph_type = "nngp",
    n = q,
    q = q,
    dim = d,
    coords = coords,
    coords_ord = coords_ord,
    coord_names = colnames(coords) %||% paste0("coord", seq_len(d)),
    m = m,
    st_scale = st_scale,
    space_time = space_time,
    search_coord_scale = if(space_time) c(rep(1, d - 1L), st_scale) else rep(1, d),
    nngp_search = nngp_search,
    ordering = ord_info$ordering_type,
    ordering_type = ord_info$ordering_type,
    ord = ord_info$ord,
    ord_inv = ord_info$ord_inv,
    params = list(
      m = m,
      st_scale = st_scale,
      space_time = space_time,
      nngp_search = nngp_search,
      ordering = ordering
    ),
    nnIndx = nn_indx,
    nnIndxLU = nn_indx_lu
  )

  class(out) <- "stLMM_graph"
  out
}

make_gp_graph <- function(coords){

  if(!is.matrix(coords))
    stop("error: coords must be a matrix")

  storage.mode(coords) <- "double"
  if(any(!is.finite(coords)))
    stop("error: GP graph coordinates must be finite")

  q <- nrow(coords)
  d <- ncol(coords)

  if(q == 0L)
    stop("error: empty coordinate matrix")

  out <- list(
    graph_type = "gp",
    n = q,
    q = q,
    dim = d,
    coords = coords,
    coords_ord = coords,
    coord_names = colnames(coords) %||% paste0("coord", seq_len(d)),
    params = list()
  )

  class(out) <- "stLMM_graph"
  out
}

make_ar1_graph <- function(time_values){

  if(length(time_values) < 1L)
    stop("error: empty time domain")

  u <- sort(unique(time_values), na.last = NA)

  if(length(u) < 2L)
    warning("ar1 graph has fewer than 2 unique time points")

  if(anyDuplicated(u))
    stop("error: AR1 support points must be unique")

  out <- list(
    graph_type = "ar1",
    n = length(u),
    support = u,
    coord_dim = 1L,
    params = list()
  )

  class(out) <- "stLMM_graph"
  out
}

car_connected_components <- function(adj){
  n <- nrow(adj)
  if(n == 0L)
    return(integer(0))

  adj <- Matrix::drop0(adj)
  trip <- Matrix::summary(adj)
  nb <- vector("list", n)
  for(i in seq_len(n))
    nb[[i]] <- integer(0)
  if(nrow(trip)){
    for(k in seq_len(nrow(trip))){
      nb[[trip$i[k]]] <- c(nb[[trip$i[k]]], trip$j[k])
    }
  }

  comp <- integer(n)
  comp_id <- 0L
  for(start in seq_len(n)){
    if(comp[start] != 0L)
      next
    comp_id <- comp_id + 1L
    queue <- start
    comp[start] <- comp_id
    head <- 1L
    while(head <= length(queue)){
      node <- queue[head]
      head <- head + 1L
      for(next_node in nb[[node]]){
        if(comp[next_node] == 0L){
          comp[next_node] <- comp_id
          queue <- c(queue, next_node)
        }
      }
    }
  }
  comp
}

connect_car_components_nearest <- function(adj, geometry, ids, island_k){
  if(island_k < 1L)
    stop("error: island_k must be a positive integer")

  comp <- car_connected_components(adj)
  n_comp <- max(comp)
  added <- data.frame(from = character(0), to = character(0), distance = numeric(0))
  if(n_comp <= 1L)
    return(list(adjacency = adj, added_edges = added, n_components_initial = n_comp))

  pts <- suppressWarnings(sf::st_point_on_surface(geometry))
  dist_mat <- as.matrix(sf::st_distance(pts))
  diag(dist_mat) <- Inf

  while(max(comp) > 1L){
    tab <- tabulate(comp, nbins = max(comp))
    comp_pick <- which.min(tab)
    in_comp <- which(comp == comp_pick)
    out_comp <- which(comp != comp_pick)

    pair_grid <- expand.grid(i = in_comp, j = out_comp)
    pair_grid$d <- dist_mat[cbind(pair_grid$i, pair_grid$j)]
    pair_grid <- pair_grid[order(pair_grid$d), , drop = FALSE]

    used <- 0L
    for(r in seq_len(nrow(pair_grid))){
      i <- pair_grid$i[r]
      j <- pair_grid$j[r]
      if(adj[i, j] != 0)
        next
      adj[i, j] <- 1
      adj[j, i] <- 1
      added <- rbind(
        added,
        data.frame(
          from = ids[i],
          to = ids[j],
          distance = as.numeric(pair_grid$d[r]),
          stringsAsFactors = FALSE
        )
      )
      used <- used + 1L
      if(used >= island_k)
        break
    }

    if(used == 0L)
      stop("error: unable to connect CAR graph islands with nearest-neighbor edges")

    adj <- Matrix::drop0(adj)
    comp <- car_connected_components(adj)
  }

  list(adjacency = adj, added_edges = added, n_components_initial = n_comp)
}

car_graph <- function(graph, id = NULL, queen = TRUE,
                      island = c("error", "nearest"), island_k = 1L){

  if(inherits(graph, "stLMM_car_graph"))
    return(graph)

  island <- match.arg(island)
  island_k <- as.integer(island_k[1])
  if(is.na(island_k) || island_k < 1L)
    stop("error: island_k must be a positive integer")

  island_added_edges <- data.frame(from = character(0), to = character(0), distance = numeric(0))
  island_components_initial <- 1L

  if(inherits(graph, "sf")){
    if(!requireNamespace("sf", quietly = TRUE))
      stop("error: package 'sf' is required to build a CAR graph from sf polygons")
    if(is.null(id) || !is.character(id) || length(id) != 1L || !nzchar(id))
      stop("error: car_graph() requires id = <column name> for sf polygon inputs")
    if(!id %in% names(graph))
      stop("error: id column ", id, " not found in sf graph")

    ids <- as.character(graph[[id]])
    if(anyNA(ids) || any(!nzchar(ids)))
      stop("error: sf id column contains missing or empty values")
    if(anyDuplicated(ids))
      stop("error: sf id column contains duplicate values")

    nb <- sf::st_touches(graph, sparse = TRUE)
    if(isFALSE(queen)){
      nb <- sf::st_relate(graph, pattern = "F***1****", sparse = TRUE)
    }

    n <- length(ids)
    ii <- rep(seq_len(n), lengths(nb))
    jj <- as.integer(unlist(nb, use.names = FALSE))

    if(length(ii)){
      keep <- ii != jj
      ii <- ii[keep]
      jj <- jj[keep]
    }

    adj <- Matrix::sparseMatrix(
      i = c(ii, jj),
      j = c(jj, ii),
      x = rep(1, length(c(ii, jj))),
      dims = c(n, n),
      dimnames = list(ids, ids)
    )
    adj <- Matrix::drop0(adj)
    adj@x[] <- 1

    comp <- car_connected_components(adj)
    island_components_initial <- max(comp)
    if(island_components_initial > 1L && island == "nearest"){
      connected <- connect_car_components_nearest(
        adj = adj,
        geometry = sf::st_geometry(graph),
        ids = ids,
        island_k = island_k
      )
      adj <- connected$adjacency
      island_added_edges <- connected$added_edges
      island_components_initial <- connected$n_components_initial
    }

  } else if(inherits(graph, "Matrix") || is.matrix(graph)) {
    if(island == "nearest")
      stop("error: island = 'nearest' requires sf polygon input")
    adj <- Matrix::Matrix(graph, sparse = TRUE)
    if(nrow(adj) != ncol(adj))
      stop("error: CAR adjacency matrix must be square")
    ids <- rownames(adj)
    if(is.null(ids))
      ids <- colnames(adj)
    if(is.null(ids))
      ids <- as.character(seq_len(nrow(adj)))
    ids <- as.character(ids)
    if(anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
      stop("error: CAR adjacency matrix row names must be non-missing unique ids")
    dimnames(adj) <- list(ids, ids)
    adj <- Matrix::drop0(adj)
    adj@x[] <- 1
  } else {
    stop("error: car_graph() expects an sf polygon object or adjacency matrix")
  }

  if(!Matrix::isSymmetric(adj))
    stop("error: CAR adjacency graph must be symmetric")

  diag(adj) <- 0
  adj <- Matrix::drop0(adj)
  degree <- as.numeric(Matrix::rowSums(adj != 0))

  if(any(degree <= 0)){
    bad <- ids[degree <= 0]
    stop("error: CAR graph contains isolated area(s): ", paste(bad, collapse = ", "))
  }

  comp <- car_connected_components(adj)
  if(max(comp) > 1L)
    stop(
      "error: CAR graph contains disconnected components; ",
      "use island = 'nearest' with sf polygon input to add nearest-neighbor bridge edges"
    )

  out <- list(
    graph_type = "car",
    ids = ids,
    n = length(ids),
    q = length(ids),
    adjacency = adj,
    geometry = if(inherits(graph, "sf")) sf::st_geometry(graph) else NULL,
    degree = degree,
    style = "binary",
    island_policy = island,
    island_k = island_k,
    island_components_initial = as.integer(island_components_initial),
    island_added_edges = island_added_edges,
    params = list(
      island_policy = island,
      island_k = island_k,
      island_components_initial = as.integer(island_components_initial),
      island_added_edges = island_added_edges
    )
  )

  class(out) <- "stLMM_car_graph"
  out
}

make_dagar_graph <- function(graph, ordering = "coord"){

  if(!inherits(graph, "stLMM_car_graph"))
    stop("error: make_dagar_graph() expects a stLMM_car_graph object")

  ord_info <- compute_graph_order(graph, ordering)
  adj <- Matrix::drop0(graph$adjacency)
  trip <- Matrix::summary(adj)
  keep <- trip$i < trip$j
  trip <- trip[keep, , drop = FALSE]

  n <- as.integer(graph$n)
  parent_lists <- vector("list", n)
  for(i in seq_len(n))
    parent_lists[[i]] <- integer(0)

  for(e in seq_len(nrow(trip))){
    i_support <- trip$i[e]
    j_support <- trip$j[e]
    i_ord <- ord_info$ord_inv[i_support]
    j_ord <- ord_info$ord_inv[j_support]

    if(i_ord == j_ord)
      stop("internal error: DAGAR ordering maps two nodes to the same position")

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
      if(length(parents) != parent_count[i])
        stop("internal error: duplicate DAGAR parent edge detected")
      if(any(parents >= i))
        stop("internal error: DAGAR parent set violates history ordering")
      parent_lists[[i]] <- parents
      parent_index[pos:(pos + parent_count[i] - 1L)] <- parents
      pos <- pos + parent_count[i]
    }
  }

  if(length(parent_index) && any(parent_index < 1L | parent_index > n))
    stop("internal error: DAGAR parent index outside 1:n")

  out <- list(
    graph_type = "dagar",
    n = n,
    q = n,
    ids = graph$ids[ord_info$ord],
    ids_support = graph$ids,
    geometry = if(!is.null(graph$geometry)) graph$geometry[ord_info$ord] else NULL,
    geometry_support = graph$geometry,
    degree = as.double(graph$degree[ord_info$ord]),
    n_edge = as.integer(nrow(trip)),
    ordering = ord_info$ordering_type,
    ordering_type = ord_info$ordering_type,
    ord = ord_info$ord,
    ord_inv = ord_info$ord_inv,
    parent_index = as.integer(parent_index - 1L),
    parent_start = as.integer(parent_start),
    parent_count = as.integer(parent_count),
    zero_parent_nodes = as.integer(which(parent_count == 0L)),
    params = c(graph$params %||% list(), list(ordering = ordering))
  )

  class(out) <- "stLMM_graph"
  out
}

make_car_graph <- function(graph, car_model = "proper"){

  if(!inherits(graph, "stLMM_car_graph"))
    stop("error: make_car_graph() expects a stLMM_car_graph object")

  if(!is.character(car_model) || length(car_model) != 1L || !nzchar(car_model))
    stop("error: make_car_graph() car_model must be a non-empty character string")
  car_model <- tolower(car_model)
  if(!car_model %in% c("proper", "leroux"))
    stop("error: make_car_graph() supports car_model = 'proper' or 'leroux'")

  adj <- Matrix::drop0(graph$adjacency)
  trip <- Matrix::summary(adj)
  keep <- trip$i < trip$j
  trip <- trip[keep, , drop = FALSE]

  out <- list(
    graph_type = "car",
    n = as.integer(graph$n),
    q = as.integer(graph$q),
    ids = graph$ids,
    geometry = graph$geometry,
    degree = as.double(graph$degree),
    edge_i = as.integer(trip$i),
    edge_j = as.integer(trip$j),
    edge_w = as.double(trip$x),
    n_edge = as.integer(nrow(trip)),
    car_model = car_model,
    params = c(graph$params %||% list(), list(car_model = car_model))
  )

  class(out) <- "stLMM_graph"
  out
}

make_car_time_graph <- function(graph, time_support, time_model = "ar1", car_model = "proper"){

  if(!inherits(graph, "stLMM_car_graph"))
    stop("error: make_car_time_graph() expects a stLMM_car_graph object")

  if(!is.character(time_model) || length(time_model) != 1L || !nzchar(time_model))
    stop("error: make_car_time_graph() time_model must be a non-empty character string")
  time_model <- tolower(time_model)
  if(!time_model %in% c("ar1", "exp"))
    stop("error: make_car_time_graph() supports time_model = 'ar1' or 'exp'")

  time_support <- sort(unique(time_support), na.last = NA)
  if(length(time_support) < 1L)
    stop("error: empty time domain for car_time graph")
  if(time_model == "exp" && (!is.numeric(time_support) || any(!is.finite(time_support))))
    stop("error: car_time() time_model = 'exp' requires numeric finite time values")
  time_delta <- if(time_model == "exp") diff(as.numeric(time_support)) else numeric(0)
  if(length(time_delta) && any(!is.finite(time_delta) | time_delta <= 0))
    stop("error: car_time() time_model = 'exp' requires strictly increasing unique time values")

  car_part <- make_car_graph(graph, car_model = car_model)
  n_space <- as.integer(graph$n)
  n_time <- as.integer(length(time_support))

  out <- car_part
  out$graph_type <- "car_time"
  out$n_space <- n_space
  out$n_time <- n_time
  out$n <- as.integer(n_space * n_time)
  out$q <- out$n
  out$time_support <- time_support
  out$time_model <- time_model
  out$time_delta <- as.double(time_delta)
  out$params <- c(out$params %||% list(), list(time_model = time_model))

  class(out) <- "stLMM_graph"
  out
}

make_dagar_time_graph <- function(graph, time_support, ordering = "coord", time_model = "ar1"){

  if(!inherits(graph, "stLMM_car_graph"))
    stop("error: make_dagar_time_graph() expects a stLMM_car_graph object")

  if(!is.character(time_model) || length(time_model) != 1L || !nzchar(time_model))
    stop("error: make_dagar_time_graph() time_model must be a non-empty character string")
  time_model <- tolower(time_model)
  if(!time_model %in% c("ar1", "exp"))
    stop("error: make_dagar_time_graph() supports time_model = 'ar1' or 'exp'")

  time_support <- sort(unique(time_support), na.last = NA)
  if(length(time_support) < 1L)
    stop("error: empty time domain for dagar_time graph")
  if(time_model == "exp" && (!is.numeric(time_support) || any(!is.finite(time_support))))
    stop("error: dagar_time() time_model = 'exp' requires numeric finite time values")
  time_delta <- if(time_model == "exp") diff(as.numeric(time_support)) else numeric(0)
  if(length(time_delta) && any(!is.finite(time_delta) | time_delta <= 0))
    stop("error: dagar_time() time_model = 'exp' requires strictly increasing unique time values")

  dagar_part <- make_dagar_graph(graph, ordering = ordering)
  n_space <- as.integer(graph$n)
  n_time <- as.integer(length(time_support))

  out <- dagar_part
  out$graph_type <- "dagar_time"
  out$n_space <- n_space
  out$n_time <- n_time
  out$n <- as.integer(n_space * n_time)
  out$q <- out$n
  out$time_support <- time_support
  out$time_model <- time_model
  out$time_delta <- as.double(time_delta)
  out$params <- c(out$params %||% list(), list(time_model = time_model))

  class(out) <- "stLMM_graph"
  out
}
