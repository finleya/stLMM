############################################################
# IID random-effect term support
############################################################

is_iid_term <- function(lbl){
  grepl("^iid\\(", lbl) || grepl(":iid\\(", lbl)
}

validate_iid_variable_name <- function(x, where){
  if(!is.character(x) || length(x) != 1L || !nzchar(x))
    stop("error: ", where, " must be a single variable name")

  if(!grepl("^[A-Za-z.][A-Za-z0-9_.]*$", x))
    stop("error: ", where, " currently supports simple variable names only")

  x
}

parse_iid_call <- function(lbl){

  colon_parts <- split_top_level_colon(lbl)

  if(is.null(colon_parts)){
    iid_part <- trimws(lbl)
    coef_name <- NULL
  } else {
    left_is_iid <- grepl("^iid\\(", colon_parts[1])
    right_is_iid <- grepl("^iid\\(", colon_parts[2])

    if(left_is_iid && !right_is_iid){
      iid_part <- colon_parts[1]
      coef_name <- trimws(colon_parts[2])
    } else if(!left_is_iid && right_is_iid){
      iid_part <- colon_parts[2]
      coef_name <- trimws(colon_parts[1])
    } else {
      stop("error: unable to parse iid random-effect term ", lbl)
    }
  }

  m <- regexec("^iid\\((.*)\\)$", iid_part)
  hit <- regmatches(iid_part, m)[[1]]

  if(length(hit) != 2L)
    stop("error: malformed iid random-effect term ", lbl)

  inside <- trimws(hit[2])
  pieces <- split_top_level(inside, split = ",")

  if(length(pieces) != 1L)
    stop("error: iid() currently accepts exactly one grouping variable")

  group_name <- validate_iid_variable_name(trimws(pieces[1]), "iid() grouping variable")

  if(!is.null(coef_name))
    coef_name <- validate_iid_variable_name(coef_name, "iid slope covariate")

  list(
    label = lbl,
    group_name = group_name,
    coef_name = coef_name,
    coefficient = if(is.null(coef_name)) "(Intercept)" else coef_name
  )
}

default_iid_names <- function(specs){
  paste0("iid_", seq_along(specs))
}

build_iid_components <- function(formula, data){

  tt <- terms(formula, data = data, keep.order = TRUE)
  term_labels <- attr(tt, "term.labels")

  iid_idx <- which(vapply(term_labels, is_iid_term, logical(1)))
  iid_labels <- term_labels[iid_idx]

  if(!length(iid_labels)){
    return(list(
      iid_labels = character(0),
      iid_idx = integer(0),
      terms = list(),
      reduced_formula = formula
    ))
  }

  specs <- lapply(iid_labels, parse_iid_call)
  iid_names <- default_iid_names(specs)

  for(i in seq_along(specs)){
    specs[[i]]$name <- iid_names[i]
    specs[[i]]$term_id <- paste0("iid", i)
  }

  reduced_formula <- formula
  for(lbl in iid_labels)
    reduced_formula <- update(reduced_formula, paste(". ~ . -", lbl))

  list(
    iid_labels = iid_labels,
    iid_names = iid_names,
    iid_idx = iid_idx,
    terms = specs,
    reduced_formula = reduced_formula
  )
}

build_iid_design <- function(iid_info, data, n_obs){

  terms <- iid_info$terms

  if(!length(terms)){
    return(list(
      Z = Matrix::Matrix(0, n_obs, 0, sparse = TRUE),
      re_names = character(0),
      re_nlevels = integer(0),
      re_ncoef = integer(0),
      re.q = integer(0),
      re_block_id = integer(0),
      terms = list()
    ))
  }

  iid_vars <- unique(unlist(lapply(terms, function(term){
    c(term$group_name, term$coef_name)
  }), use.names = FALSE))
  iid_vars <- iid_vars[nzchar(iid_vars)]

  iid_formula <- stats::reformulate(iid_vars)
  environment(iid_formula) <- environment(iid_info$reduced_formula)
  iid_mf <- stats::model.frame(iid_formula, data, na.action = stats::na.fail)

  if(nrow(iid_mf) != n_obs)
    stop("error: iid random-effect data rows do not match fixed-effect model rows")

  for(i in seq_along(terms)){
    g <- terms[[i]]$group_name
    if(!g %in% names(iid_mf))
      stop("error: grouping factor ", g, " not found for ", terms[[i]]$label)

    group <- iid_mf[[g]]
    if(is.factor(group)){
      levels_i <- levels(group)
      group <- factor(group, levels = levels_i)
    } else {
      group <- factor(group)
      levels_i <- levels(group)
    }

    if(anyNA(group))
      stop("error: grouping factor ", g, " contains missing values")
    if(length(levels_i) < 1L)
      stop("error: grouping factor ", g, " has no levels")

    terms[[i]]$levels <- levels_i
    terms[[i]]$n_levels <- length(levels_i)

    if(!is.null(terms[[i]]$coef_name)){
      x <- iid_mf[[terms[[i]]$coef_name]]
      if(!is.numeric(x))
        stop("error: iid slope covariate ", terms[[i]]$coef_name, " must be numeric")
      x <- as.numeric(x)
      if(length(x) != n_obs || anyNA(x))
        stop("error: iid slope covariate ", terms[[i]]$coef_name, " contains missing values")
      terms[[i]]$scale <- x
    } else {
      terms[[i]]$scale <- rep(1.0, n_obs)
    }

    terms[[i]]$group_index <- as.integer(group)
  }

  re.q <- vapply(terms, `[[`, integer(1), "n_levels")
  q <- sum(re.q)
  n_value <- n_obs * length(terms)
  ii <- integer(n_value)
  jj <- integer(n_value)
  xx <- numeric(n_value)
  pos <- 0L
  col_offset <- 0L

  for(i in seq_along(terms)){
    term <- terms[[i]]
    idx <- term$group_index
    scale <- term$scale

    for(row in seq_len(n_obs)){
      pos <- pos + 1L
      ii[pos] <- row
      jj[pos] <- col_offset + idx[row]
      xx[pos] <- scale[row]
    }

    col_offset <- col_offset + term$n_levels
  }

  Z <- Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(n_obs, q))

  z_names <- character(q)
  col_offset <- 0L
  for(i in seq_along(terms)){
    idx <- col_offset + seq_len(terms[[i]]$n_levels)
    z_names[idx] <- paste0(terms[[i]]$name, "_", make.names(terms[[i]]$levels, unique = TRUE))
    col_offset <- col_offset + terms[[i]]$n_levels
  }
  colnames(Z) <- z_names

  list(
    Z = Z,
    re_names = vapply(terms, `[[`, character(1), "name"),
    re_nlevels = as.integer(re.q),
    re_ncoef = rep.int(1L, length(terms)),
    re.q = as.integer(re.q),
    re_block_id = rep.int(seq_along(terms), re.q),
    terms = terms
  )
}
