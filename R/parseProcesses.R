
`%||%` <- function(x, y) if(is.null(x)) y else x

is_process_term <- function(lbl){
  grepl("^(nngp|gp|ar1|car|car_time|dagar|dagar_time)\\(", lbl) ||
    grepl(":(nngp|gp|ar1|car|car_time|dagar|dagar_time)\\(", lbl)
}

split_top_level <- function(x, split = ","){

  chars <- strsplit(x, "", fixed = TRUE)[[1]]
  depth <- 0L
  out <- character(0)
  start <- 1L

  for(i in seq_along(chars)){
    ch <- chars[i]
    if(ch == "(") depth <- depth + 1L
    if(ch == ")") depth <- depth - 1L

    if(ch == split && depth == 0L){
      out <- c(out, trimws(substr(x, start, i - 1L)))
      start <- i + 1L
    }
  }

  out <- c(out, trimws(substr(x, start, nchar(x))))
  out[nzchar(out)]
}

split_top_level_colon <- function(x){

  chars <- strsplit(x, "", fixed = TRUE)[[1]]
  depth <- 0L

  for(i in seq_along(chars)){
    ch <- chars[i]
    if(ch == "(") depth <- depth + 1L
    if(ch == ")") depth <- depth - 1L

    if(ch == ":" && depth == 0L){
      return(c(
        trimws(substr(x, 1L, i - 1L)),
        trimws(substr(x, i + 1L, nchar(x)))
      ))
    }
  }

  NULL
}

parse_process_call <- function(lbl){

  colon_parts <- split_top_level_colon(lbl)

  if(is.null(colon_parts)){
    proc_part <- trimws(lbl)
    coef_name <- NULL
  } else {

    left_is.proc  <- grepl("^[A-Za-z0-9_.]+\\(", colon_parts[1])
    right_is.proc <- grepl("^[A-Za-z0-9_.]+\\(", colon_parts[2])

    if(left_is.proc && !right_is.proc){
      proc_part <- colon_parts[1]
      coef_name <- trimws(colon_parts[2])
    } else if(!left_is.proc && right_is.proc){
      proc_part <- colon_parts[2]
      coef_name <- trimws(colon_parts[1])
    } else {
      stop("error: unable to parse process interaction term ", lbl)
    }
  }

  m <- regexec("^([A-Za-z0-9_.]+)\\((.*)\\)$", proc_part)
  hit <- regmatches(proc_part, m)[[1]]

  if(length(hit) != 3L)
    stop("error: malformed process term ", lbl)

  fun <- tolower(hit[2])
  inside <- trimws(hit[3])

  pipe_parts <- split_top_level(inside, split = "|")

  if(length(pipe_parts) > 1L)
    stop(
      "error: use ",
      fun,
      "(...) without | in process terms: ",
      lbl,
      "; use iid(group) for grouped iid random effects"
    )

  param_string <- trimws(pipe_parts[1])

  pieces <- split_top_level(param_string, split = ",")

  args <- character(0)
  params <- list()

  for(p in pieces){

    if(grepl("=", p, fixed = TRUE)){

      kv <- strsplit(p, "=", fixed = TRUE)[[1]]

      if(length(kv) != 2L)
        stop("error: malformed parameter in ", lbl)

      key <- trimws(kv[1])
      val <- trimws(kv[2])

      if(identical(key, "time_model")){
        val <- gsub('^["\']|["\']$', "", val)
      } else {
        val <- tryCatch(eval(parse(text = val)), error = function(e) val)
      }
      params[[key]] <- val

    } else {
      args <- c(args, trimws(p))
    }
  }

  list(
    label = lbl,
    fun = fun,
    args = args,
    params = params,
    group = NULL,
    coef_name = coef_name
  )
}

nngp_graph_defaults <- function(){
  list(
    m = 15L,
    ordering = "coord",
    st_scale = 1
  )
}

nngp_term_defaults <- function(){
  list(
    cov_model = "exp"
  )
}

gp_graph_defaults <- function(){
  list()
}

gp_term_defaults <- function(){
  list(
    cov_model = "exp"
  )
}

ar1_term_defaults <- function(){
  list()
}

car_term_defaults <- function(){
  list(car_model = "proper")
}

car_time_term_defaults <- function(){
  list(time_model = "ar1", car_model = "proper")
}

dagar_graph_defaults <- function(){
  list(ordering = "coord")
}

resolve_nngp_spec <- function(spec,
                              registry = build_cor_model_registry()){

  ####################################################
  ## Coordinate canonicalization
  ####################################################

  coords <- trimws(spec$args)
  coord_count <- length(coords)

  ####################################################
  ## Graph + operator parameter split
  ####################################################

  graph_params <- modifyList(nngp_graph_defaults(), spec$params)
  operator_params <- modifyList(nngp_term_defaults(), spec$params)

  operator_params$m <- NULL
  operator_params$ordering <- NULL
  operator_params$st_scale <- NULL

  ####################################################
  ## Correlation model
  ####################################################

  cov_model <- operator_params$cov_model

  if(!is.character(cov_model) || length(cov_model) != 1L || !nzchar(cov_model))
    stop("error: nngp(): cov_model must be a non-empty character string")

  if(!cov_model %in% names(registry))
    stop("error: nngp(): unknown cov_model '", cov_model, "'")

  operator_params$cov_model <- cov_model

  ####################################################
  ## Coordinate validation based on distance mode
  ####################################################

  dist_mode <- registry[[cov_model]]$distance_mode

  if(dist_mode == 1L){  # COR_SINGLE

    if(coord_count < 1L)
      stop(
        "error: cov_model='", cov_model,
        "' requires at least one coordinate"
      )

  } else if(dist_mode == 2L){  # COR_SPACE_TIME

    if(coord_count < 2L)
      stop(
        "error: cov_model='", cov_model,
        "' requires at least two coordinates (space + time)"
      )

  }

  ####################################################
  ## Graph parameters
  ####################################################

  m <- graph_params$m

  if(!is.numeric(m) || length(m) != 1L || is.na(m) || m < 1 ||
     abs(m - round(m)) > 0)
    stop("error: nngp(): m must be a positive integer")

  graph_params$m <- as.integer(m)

  ordering <- graph_params$ordering

  if(is.character(ordering)){

    if(length(ordering) != 1L || !nzchar(ordering))
      stop("error: nngp(): ordering string must be non-empty")

    ok <- c("coord", "default", "maxmin", "hilbert", "random")

    if(!ordering %in% ok)
      stop(
        "error: nngp(): ordering must be one of ",
        paste(shQuote(ok), collapse = ", "),
        " or a numeric permutation"
      )

  } else if(is.numeric(ordering)){

    if(length(ordering) < 1L)
      stop("error: nngp(): numeric ordering must have positive length")

    ordering <- as.integer(ordering)

  } else {

    stop("error: nngp(): ordering must be a character option or numeric permutation")

  }

  graph_params$ordering <- ordering

  st_scale <- graph_params$st_scale
  if(!is.numeric(st_scale) || length(st_scale) != 1L ||
     is.na(st_scale) || st_scale <= 0)
    stop("error: nngp(): st_scale must be a positive scalar")

  st_scale <- as.numeric(st_scale)
  if(dist_mode != 2L && !isTRUE(all.equal(st_scale, 1)))
    stop("error: nngp(): st_scale is only valid for space-time covariance models")

  if(dist_mode == 2L && identical(ordering, "hilbert"))
    stop("error: nngp(): hilbert ordering is not supported for space-time covariance models; use coord, default, maxmin, random, or a numeric ordering")

  graph_params$st_scale <- st_scale

  ####################################################

  list(
    coords = coords,
    graph_params = graph_params,
    operator_params = operator_params
  )
}

build_nngp_components <- function(spec,
                                  data,
                                  n_omp_threads = 1,
                                  nngp_search = "fast",
                                  warn = TRUE,
                                  registry = build_cor_model_registry()){

  resolved <- resolve_nngp_spec(
    spec = spec,
    registry = registry
  )

  graph_params <- resolved$graph_params
  operator_params <- resolved$operator_params

  coords <- extract_process_coords(data, resolved$coords)

  coords_unique <- unique_coordinate_rows(coords)

  n_obs <- nrow(coords)
  n_node <- nrow(coords_unique)

  if(warn && n_node < n_obs){
    warning(
      sprintf(
        "duplicated coordinates detected in nngp(%s); fitting latent process on %d unique locations (from %d observations)",
        paste(resolved$coords, collapse = ","),
        n_node,
        n_obs
      ),
      call. = FALSE
    )
  }

  map <- as.integer(match_unique_rows(coords, coords_unique))

  if(anyNA(map))
    stop("error: spatial node mapping produced NA values")

  if(is.numeric(graph_params$ordering) && length(graph_params$ordering) != n_node)
    stop("error: nngp(): custom ordering vector must have length q = ", n_node, " unique spatial nodes")

  space_time <- identical(as.integer(registry[[operator_params$cov_model]]$distance_mode), 2L)

  graph <- make_nngp_graph(
    coords = coords_unique,
    m = graph_params$m,
    ordering = graph_params$ordering,
    st_scale = graph_params$st_scale,
    space_time = space_time,
    n_omp_threads = n_omp_threads,
    nngp_search = nngp_search
  )

  ## map is in original unique-coordinate order; remap into graph
  ## (reordered) node order
  map <- as.integer(graph$ord_inv[map])

  if(anyNA(map))
    stop("internal error: reordered spatial node mapping produced NA values")

  # sanity check: node mapping matches graph
  if(max(map) > graph$n)
    stop("internal error: map and graph node count mismatch")

  list(
    graph = graph,
    term = make_process_term(
      term_type = "nngp",
      graph_id = NA_character_,
      label = spec$label,
      name = NULL,
      map = map,
      coef_name = spec$coef_name,
      cov_model = operator_params$cov_model,
      params = operator_params,
      data = data
    )
  )
}

resolve_gp_spec <- function(spec,
                            registry = build_cor_model_registry()){

  coords <- trimws(spec$args)
  coord_count <- length(coords)

  graph_params <- modifyList(gp_graph_defaults(), spec$params)
  operator_params <- modifyList(gp_term_defaults(), spec$params)

  cov_model <- operator_params$cov_model

  if(!is.character(cov_model) || length(cov_model) != 1L || !nzchar(cov_model))
    stop("error: gp(): cov_model must be a non-empty character string")

  if(!cov_model %in% names(registry))
    stop("error: gp(): unknown cov_model '", cov_model, "'")

  operator_params$cov_model <- cov_model

  dist_mode <- registry[[cov_model]]$distance_mode

  if(dist_mode == 1L){
    if(coord_count < 1L)
      stop("error: cov_model='", cov_model, "' requires at least one coordinate")
  } else if(dist_mode == 2L){
    if(coord_count < 2L)
      stop("error: cov_model='", cov_model, "' requires at least two coordinates (space + time)")
  }

  list(
    coords = coords,
    graph_params = graph_params,
    operator_params = operator_params
  )
}

build_gp_components <- function(spec,
                                data,
                                warn = TRUE,
                                registry = build_cor_model_registry()){

  resolved <- resolve_gp_spec(spec = spec, registry = registry)
  graph_params <- resolved$graph_params
  operator_params <- resolved$operator_params

  coords <- extract_process_coords(data, resolved$coords)

  coords_unique <- unique_coordinate_rows(coords)

  n_obs <- nrow(coords)
  n_node <- nrow(coords_unique)

  if(warn && n_node < n_obs){
    warning(
      sprintf(
        "duplicated coordinates detected in gp(%s); fitting latent process on %d unique locations (from %d observations)",
        paste(resolved$coords, collapse = ","),
        n_node,
        n_obs
      ),
      call. = FALSE
    )
  }

  map <- as.integer(match_unique_rows(coords, coords_unique))

  if(anyNA(map))
    stop("error: spatial node mapping produced NA values")

  graph <- make_gp_graph(coords_unique)

  list(
    graph = graph,
    term = make_process_term(
      term_type = "gp",
      graph_id = NA_character_,
      label = spec$label,
      name = NULL,
      map = map,
      coef_name = spec$coef_name,
      cov_model = operator_params$cov_model,
      params = operator_params,
      data = data
    )
  )
}

build_ar1_components <- function(spec, data, warn = TRUE){

  if(length(spec$args) != 1L)
    stop("error: ar1() currently expects exactly one time variable")

  time_name <- spec$args[1]

  if(!time_name %in% names(data))
    stop("error: time variable ", time_name, " not found in data")

  tt <- data[[time_name]]

  if(length(spec$params)){
    bad <- names(spec$params)
    stop(
      "error: ar1() does not currently accept formula-level parameters: ",
      paste(bad, collapse = ", "),
      ". Correlation parameter(s) are handled later, and variance scaling is handled separately."
    )
  }

  if(is.factor(tt))
    tt <- as.character(tt)

  support <- sort(unique(tt), na.last = NA)

  n_obs <- length(tt)
  n_node <- length(support)

  if(warn && n_node < n_obs){
    warning(
      sprintf(
        "duplicated time values detected in ar1(%s); fitting latent process on %d unique time points (from %d observations)",
        time_name,
        n_node,
        n_obs
      ),
      call. = FALSE
    )
  }

  map <- match(tt, support)

  if(anyNA(map))
    stop("error: AR1 time mapping produced NA values")

  graph <- make_ar1_graph(support)

  list(
    graph = graph,
    term = make_process_term(
      term_type = "ar1",
      graph_id = NA_character_,
      label = spec$label,
      name = NULL,
      map = map,
      coef_name = spec$coef_name,
      cov_model = NULL,
      params = ar1_term_defaults(),
      data = data
    )
  )
}

eval_process_param <- function(x, env){
  if(!is.character(x) || length(x) != 1L)
    return(x)
  tryCatch(eval(parse(text = x), envir = env), error = function(e) x)
}

build_car_components <- function(spec, data, formula_env){

  if(length(spec$args) != 1L)
    stop("error: car() currently expects exactly one area id variable")

  area_name <- spec$args[1]

  if(!area_name %in% names(data))
    stop("error: area variable ", area_name, " not found in data")

  params <- spec$params
  if(is.null(params$graph))
    stop("error: car() requires graph = <car_graph object>")

  graph <- eval_process_param(params$graph, formula_env)
  graph_id_col <- eval_process_param(params$graph_id %||% params$id, formula_env)
  queen <- eval_process_param(params$queen %||% TRUE, formula_env)
  car_model <- normalize_car_model(params$car_model %||% "proper", formula_env)

  if(!is.character(car_model) || length(car_model) != 1L || !nzchar(car_model))
    stop("error: car() car_model must be a non-empty character string")
  car_model <- tolower(car_model)
  if(!car_model %in% c("proper", "leroux"))
    stop("error: car() car_model must be 'proper' or 'leroux'")

  if(!inherits(graph, "stLMM_car_graph"))
    graph <- car_graph(graph, id = graph_id_col, queen = queen)

  ids <- as.character(graph$ids)
  area <- as.character(data[[area_name]])

  if(anyNA(area))
    stop("error: area variable ", area_name, " contains missing values")

  map <- match(area, ids)
  if(anyNA(map)){
    bad <- unique(area[is.na(map)])
    stop("error: car() area value(s) not found in graph: ", paste(bad, collapse = ", "))
  }

  bad_params <- setdiff(names(params), c("graph", "graph_id", "id", "queen", "car_model"))
  if(length(bad_params))
    stop("error: car() does not accept formula-level parameter(s): ", paste(bad_params, collapse = ", "))

  list(
    graph = make_car_graph(graph, car_model = car_model),
    term = make_process_term(
      term_type = "car",
      graph_id = NA_character_,
      label = spec$label,
      name = NULL,
      map = as.integer(map),
      coef_name = spec$coef_name,
      cov_model = NULL,
      params = list(car_model = car_model),
      data = data
    )
  )
}

build_dagar_components <- function(spec, data, formula_env){

  if(length(spec$args) != 1L)
    stop("error: dagar() currently expects exactly one area id variable")

  area_name <- spec$args[1]

  if(!area_name %in% names(data))
    stop("error: area variable ", area_name, " not found in data")

  params <- modifyList(dagar_graph_defaults(), spec$params)
  if(is.null(params$graph))
    stop("error: dagar() requires graph = <car_graph object>")

  graph <- eval_process_param(params$graph, formula_env)
  graph_id_col <- eval_process_param(params$graph_id %||% params$id, formula_env)
  queen <- eval_process_param(params$queen %||% TRUE, formula_env)
  ordering <- eval_process_param(params$ordering, formula_env)

  if(!inherits(graph, "stLMM_car_graph"))
    graph <- car_graph(graph, id = graph_id_col, queen = queen)

  ids <- as.character(graph$ids)
  area <- as.character(data[[area_name]])

  if(anyNA(area))
    stop("error: area variable ", area_name, " contains missing values")

  support_map <- match(area, ids)
  if(anyNA(support_map)){
    bad <- unique(area[is.na(support_map)])
    stop("error: dagar() area value(s) not found in graph: ", paste(bad, collapse = ", "))
  }

  bad_params <- setdiff(names(params), c("graph", "graph_id", "id", "queen", "ordering"))
  if(length(bad_params))
    stop("error: dagar() does not accept formula-level parameter(s): ", paste(bad_params, collapse = ", "))

  dagar_graph <- make_dagar_graph(graph, ordering = ordering)
  map <- as.integer(dagar_graph$ord_inv[support_map])

  list(
    graph = dagar_graph,
    term = make_process_term(
      term_type = "dagar",
      graph_id = NA_character_,
      label = spec$label,
      name = NULL,
      map = map,
      coef_name = spec$coef_name,
      cov_model = NULL,
      params = list(ordering = dagar_graph$ordering),
      data = data,
      n_node = length(ids)
    )
  )
}

build_car_time_components <- function(spec, data, formula_env, warn = TRUE){

  if(length(spec$args) != 2L)
    stop("error: car_time() expects exactly two variables: area id and time")

  area_name <- spec$args[1]
  time_name <- spec$args[2]

  if(!area_name %in% names(data))
    stop("error: area variable ", area_name, " not found in data")
  if(!time_name %in% names(data))
    stop("error: time variable ", time_name, " not found in data")

  params <- spec$params
  if(is.null(params$graph))
    stop("error: car_time() requires graph = <car_graph object>")

  graph <- eval_process_param(params$graph, formula_env)
  graph_id_col <- eval_process_param(params$graph_id %||% params$id, formula_env)
  queen <- eval_process_param(params$queen %||% TRUE, formula_env)
  time_model <- normalize_car_time_model(params$time_model %||% "ar1", formula_env)
  car_model <- normalize_car_model(params$car_model %||% "proper", formula_env)

  if(!is.character(time_model) || length(time_model) != 1L || !nzchar(time_model))
    stop("error: car_time() time_model must be a non-empty character string")
  time_model <- tolower(time_model)
  if(!time_model %in% c("ar1", "exp"))
    stop("error: car_time() time_model must be 'ar1' or 'exp'")
  if(!is.character(car_model) || length(car_model) != 1L || !nzchar(car_model))
    stop("error: car_time() car_model must be a non-empty character string")
  car_model <- tolower(car_model)
  if(!car_model %in% c("proper", "leroux"))
    stop("error: car_time() car_model must be 'proper' or 'leroux'")

  if(!inherits(graph, "stLMM_car_graph"))
    graph <- car_graph(graph, id = graph_id_col, queen = queen)

  ids <- as.character(graph$ids)
  area <- as.character(data[[area_name]])
  tt <- data[[time_name]]

  if(anyNA(area))
    stop("error: area variable ", area_name, " contains missing values")
  if(anyNA(tt))
    stop("error: time variable ", time_name, " contains missing values")

  if(is.factor(tt))
    tt <- as.character(tt)
  if(time_model == "exp" && (!is.numeric(tt) || any(!is.finite(tt))))
    stop("error: car_time() time_model = 'exp' requires numeric finite time values")
  time_support <- sort(unique(tt), na.last = NA)

  space_index <- match(area, ids)
  if(anyNA(space_index)){
    bad <- unique(area[is.na(space_index)])
    stop("error: car_time() area value(s) not found in graph: ", paste(bad, collapse = ", "))
  }

  time_index <- match(tt, time_support)
  if(anyNA(time_index))
    stop("error: car_time() time mapping produced NA values")

  bad_params <- setdiff(names(params), c("graph", "graph_id", "id", "queen", "time_model", "car_model"))
  if(length(bad_params))
    stop("error: car_time() does not accept formula-level parameter(s): ", paste(bad_params, collapse = ", "))

  n_space <- length(ids)
  n_time <- length(time_support)
  map <- as.integer((space_index - 1L) * n_time + time_index)

  if(warn && length(unique(map)) < length(map)){
    warning(
      sprintf(
        "duplicated area-time values detected in car_time(%s, %s); fitting latent process on %d area-time nodes (from %d observations)",
        area_name,
        time_name,
        n_space * n_time,
        length(map)
      ),
      call. = FALSE
    )
  }

  list(
    graph = make_car_time_graph(
      graph,
      time_support = time_support,
      time_model = time_model,
      car_model = car_model
    ),
    term = make_process_term(
      term_type = "car_time",
      graph_id = NA_character_,
      label = spec$label,
      name = NULL,
      map = map,
      coef_name = spec$coef_name,
      cov_model = NULL,
      params = list(time_model = time_model, car_model = car_model),
      data = data,
      n_node = n_space * n_time
    )
  )
}

build_dagar_time_components <- function(spec, data, formula_env, warn = TRUE){

  if(length(spec$args) != 2L)
    stop("error: dagar_time() expects exactly two variables: area id and time")

  area_name <- spec$args[1]
  time_name <- spec$args[2]

  if(!area_name %in% names(data))
    stop("error: area variable ", area_name, " not found in data")
  if(!time_name %in% names(data))
    stop("error: time variable ", time_name, " not found in data")

  params <- modifyList(dagar_graph_defaults(), spec$params)
  if(is.null(params$graph))
    stop("error: dagar_time() requires graph = <car_graph object>")

  graph <- eval_process_param(params$graph, formula_env)
  graph_id_col <- eval_process_param(params$graph_id %||% params$id, formula_env)
  queen <- eval_process_param(params$queen %||% TRUE, formula_env)
  ordering <- eval_process_param(params$ordering, formula_env)
  time_model <- normalize_car_time_model(params$time_model %||% "ar1", formula_env)

  if(!is.character(time_model) || length(time_model) != 1L || !nzchar(time_model))
    stop("error: dagar_time() time_model must be a non-empty character string")
  time_model <- tolower(time_model)
  if(!time_model %in% c("ar1", "exp"))
    stop("error: dagar_time() time_model must be 'ar1' or 'exp'")

  if(!inherits(graph, "stLMM_car_graph"))
    graph <- car_graph(graph, id = graph_id_col, queen = queen)

  ids <- as.character(graph$ids)
  area <- as.character(data[[area_name]])
  tt <- data[[time_name]]

  if(anyNA(area))
    stop("error: area variable ", area_name, " contains missing values")
  if(anyNA(tt))
    stop("error: time variable ", time_name, " contains missing values")

  if(is.factor(tt))
    tt <- as.character(tt)
  if(time_model == "exp" && (!is.numeric(tt) || any(!is.finite(tt))))
    stop("error: dagar_time() time_model = 'exp' requires numeric finite time values")
  time_support <- sort(unique(tt), na.last = NA)

  support_map <- match(area, ids)
  if(anyNA(support_map)){
    bad <- unique(area[is.na(support_map)])
    stop("error: dagar_time() area value(s) not found in graph: ", paste(bad, collapse = ", "))
  }

  time_index <- match(tt, time_support)
  if(anyNA(time_index))
    stop("error: dagar_time() time mapping produced NA values")

  bad_params <- setdiff(names(params), c("graph", "graph_id", "id", "queen", "ordering", "time_model"))
  if(length(bad_params))
    stop("error: dagar_time() does not accept formula-level parameter(s): ", paste(bad_params, collapse = ", "))

  dagar_graph <- make_dagar_time_graph(
    graph,
    time_support = time_support,
    ordering = ordering,
    time_model = time_model
  )

  n_space <- length(ids)
  n_time <- length(time_support)
  space_index <- dagar_graph$ord_inv[support_map]
  map <- as.integer((space_index - 1L) * n_time + time_index)

  if(warn && length(unique(map)) < length(map)){
    warning(
      sprintf(
        "duplicated area-time values detected in dagar_time(%s, %s); fitting latent process on %d area-time nodes (from %d observations)",
        area_name,
        time_name,
        n_space * n_time,
        length(map)
      ),
      call. = FALSE
    )
  }

  list(
    graph = dagar_graph,
    term = make_process_term(
      term_type = "dagar_time",
      graph_id = NA_character_,
      label = spec$label,
      name = NULL,
      map = map,
      coef_name = spec$coef_name,
      cov_model = NULL,
      params = list(ordering = dagar_graph$ordering, time_model = time_model),
      data = data,
      n_node = n_space * n_time
    )
  )
}

process_graph_defaults <- function(fun){
  switch(
    fun,
    nngp = nngp_graph_defaults(),
    gp = gp_graph_defaults(),
    ar1 = list(),
    car = list(),
    car_time = list(),
    dagar = dagar_graph_defaults(),
    dagar_time = dagar_graph_defaults(),
    stop("error: graph defaults not implemented for ", fun)
  )
}

default_process_names <- function(specs){
  counts <- setNames(integer(0), character(0))
  out <- character(length(specs))

  for(i in seq_along(specs)){
    fun <- tolower(specs[[i]]$fun)
    if(!fun %in% names(counts))
      counts[fun] <- 0L
    counts[fun] <- counts[fun] + 1L
    out[i] <- paste0(fun, "_", counts[fun])
  }

  out
}

build_process_components <- function(formula,
                                     data,
                                     n_omp_threads = 1,
                                     nngp_search = "fast",
                                     registry = build_cor_model_registry()){

  tt <- terms(formula, data = data, keep.order = TRUE)
  term_labels <- attr(tt, "term.labels")

  proc_idx <- which(vapply(term_labels, is_process_term, logical(1)))
  proc_labels <- term_labels[proc_idx]

  if(!length(proc_labels)){
    return(list(
      process_labels = character(0),
      process_idx = integer(0),
      graphs = list(),
      process_terms = list(),
      reduced_formula = formula,
      specs = list()
    ))
  }

  specs <- lapply(proc_labels, parse_process_call)
  proc_names <- default_process_names(specs)

  keys <- vapply(
    specs,
    function(spec) build_graph_key(
      spec,
      defaults = process_graph_defaults(spec$fun),
      formula_env = environment(formula)
    ),
    character(1)
  )

  ukeys <- unique(keys)
  graph_ids <- paste0("g", seq_along(ukeys))
  key_to_id <- setNames(graph_ids, ukeys)

  graphs <- vector("list", length(ukeys))
  names(graphs) <- graph_ids

  for(i in seq_along(ukeys)){

    spec <- specs[[match(ukeys[i], keys)]]

    comp <- switch(
      spec$fun,
      nngp = build_nngp_components(
        spec,
        data,
        n_omp_threads = n_omp_threads,
        nngp_search = nngp_search,
        warn = TRUE,
        registry = registry
      ),
      gp = build_gp_components(spec, data, warn = TRUE, registry = registry),
      ar1 = build_ar1_components(spec, data, warn = TRUE),
      car = build_car_components(spec, data, environment(formula)),
      dagar = build_dagar_components(spec, data, environment(formula)),
      car_time = build_car_time_components(spec, data, environment(formula), warn = TRUE),
      dagar_time = build_dagar_time_components(spec, data, environment(formula), warn = TRUE),
      stop("error: process builder not implemented for ", spec$fun)
    )

    graphs[[i]] <- comp$graph
    graphs[[i]]$graph_id <- graph_ids[i]
    graphs[[i]]$graph_key <- ukeys[i]
  }

  terms_out <- vector("list", length(specs))
  names(terms_out) <- paste0("s", seq_along(specs))

  for(i in seq_along(specs)){

    spec <- specs[[i]]
    gid <- key_to_id[[keys[i]]]

    comp <- switch(
      spec$fun,
      nngp = build_nngp_components(
        spec,
        data,
        n_omp_threads = n_omp_threads,
        nngp_search = nngp_search,
        warn = FALSE,
        registry = registry
      ),
      gp = build_gp_components(spec, data, warn = FALSE, registry = registry),
      ar1 = build_ar1_components(spec, data, warn = FALSE),
      car = build_car_components(spec, data, environment(formula)),
      dagar = build_dagar_components(spec, data, environment(formula)),
      car_time = build_car_time_components(spec, data, environment(formula), warn = FALSE),
      dagar_time = build_dagar_time_components(spec, data, environment(formula), warn = FALSE),
      stop("error: process builder not implemented for ", spec$fun)
    )

    graph_index <- match(gid, graph_ids)
    if(spec$fun == "nngp"){
      ## Random NNGP graph construction is intentionally stochastic. Graphs are
      ## built once for backend reuse above; if this term build drew a different
      ## random order, rebase the observation map onto the retained graph.
      support_map <- comp$graph$ord[comp$term$map]
      comp$term$map <- as.integer(graphs[[graph_index]]$ord_inv[support_map])
      obs_info <- build_obs_index_with_n_node(comp$term$map, graphs[[graph_index]]$q)
      comp$term$obsIndx <- obs_info$obsIndx
      comp$term$obsIndxLU <- obs_info$obsIndxLU
      comp$term$node_nobs <- obs_info$node_nobs
      comp$term$n_node <- obs_info$n_node
    }

    terms_out[[i]] <- comp$term
    terms_out[[i]]$graph_id <- gid
    terms_out[[i]]$term_id <- names(terms_out)[i]
    terms_out[[i]]$name <- proc_names[i]
  }

  reduced_formula <- formula

  if(length(proc_labels)){
    for(lbl in proc_labels){
      reduced_formula <- update(reduced_formula, paste(". ~ . -", lbl))
    }
  }

  list(
    process_labels = proc_labels,
    process_names = proc_names,
    process_idx = proc_idx,
    graphs = graphs,
    process_terms = terms_out,
    reduced_formula = reduced_formula,
    specs = specs
  )
}
