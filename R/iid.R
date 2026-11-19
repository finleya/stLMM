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

# A one-column matrix (for example, scale(x)) represents one numeric covariate.
# Normalize its storage without changing the supplied values.
iid_numeric_covariate <- function(x, n, where){
  if(!is.numeric(x) || length(x) != n ||
     (!is.null(dim(x)) && (!is.matrix(x) || ncol(x) != 1L || nrow(x) != n)))
    stop("error: ", where, " must be a finite numeric vector or one-column matrix with ", n, " rows")
  if(anyNA(x))
    stop("error: ", where, " contains missing values")
  if(any(!is.finite(x)))
    stop("error: ", where, " must be finite numeric values")
  as.numeric(x)
}

# Covariance draws store the lower triangle by column, including the diagonal.
# Both output labeling and prediction use this order.
iid_covariance_sample_names <- function(term){
  coefficients <- term$coefficients
  k <- length(coefficients)
  ij <- which(lower.tri(matrix(0, k, k), diag = TRUE), arr.ind = TRUE)
  paste0(term$name, "_cov[", coefficients[ij[, 1]], ",", coefficients[ij[, 2]], "]")
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

  if(length(pieces) != 1L){
    spec <- parse_iid_coefficients(iid_part, lbl)
    if(!is.null(coef_name)){
      if(!is.null(spec$cov_model))
        stop("error: an iid coefficient formula cannot be combined with an external slope")
      spec$coef_name <- validate_iid_variable_name(coef_name, "iid slope covariate")
      spec$coefficient <- spec$coef_name
    }
    return(spec)
  }

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

parse_iid_coefficients <- function(iid_part, lbl){
  args <- as.list(str2lang(iid_part))[-1L]
  nms <- names(args)
  if(is.null(nms)) nms <- rep("", length(args))
  allowed <- c("group", "coefficients", "cov_model", "by")
  if(length(args) > length(allowed))
    stop("error: iid() accepts group, coefficients, cov_model, and by")
  named <- nms[nzchar(nms)]
  if(anyDuplicated(named) || any(!named %in% allowed))
    stop("error: iid() accepts group, coefficients, cov_model, and by without duplicate arguments")
  positional <- which(!nzchar(nms))
  nms[positional] <- setdiff(allowed, named)[seq_along(positional)]
  names(args) <- nms
  if(!is.symbol(args$group))
    stop("error: iid() grouping variable must be a simple variable name")
  group <- validate_iid_variable_name(as.character(args$group), "iid() grouping variable")
  by_name <- NULL
  if("by" %in% names(args)){
    if(!is.symbol(args$by))
      stop("error: iid() by must be a simple variable name")
    by_name <- validate_iid_variable_name(as.character(args$by), "iid() by variable")
  }
  f <- args$coefficients
  if(is.null(f) && !is.null(by_name)){
    if(!is.null(args$cov_model))
      stop("error: iid() cov_model requires a coefficient formula")
    return(list(label = lbl, group_name = group, coef_name = NULL,
                coefficient = "(Intercept)", by_name = by_name))
  }
  if(!is.call(f) || !identical(f[[1L]], as.name("~")) || length(f) != 2L)
    stop("error: iid() coefficients must be a one-sided formula, e.g. ~ 1 + x")
  additive <- function(x){
    if(is.symbol(x)) return(as.character(x) != ".")
    if(is.numeric(x)) return(length(x) == 1L && x %in% c(0, 1))
    is.call(x) && as.character(x[[1L]]) %in% c("+", "-") &&
      all(vapply(as.list(x)[-1L], additive, logical(1)))
  }
  if(!additive(f[[2L]]))
    stop("error: iid() coefficients currently supports additive numeric variable names and an optional intercept")
  tt <- stats::terms(stats::as.formula(f))
  vars <- attr(tt, "term.labels")
  for(v in vars) validate_iid_variable_name(v, "iid slope covariate")
  coefficients <- c(if(attr(tt, "intercept") == 1L) "(Intercept)", vars)
  if(!length(coefficients))
    stop("error: iid() coefficient formula must include at least one coefficient")
  cov_model <- args$cov_model %||% "unstructured"
  if(!identical(cov_model, "unstructured") && !identical(cov_model, "independent"))
    stop('error: iid() cov_model must be "unstructured" or "independent"')
  out <- list(label = lbl, group_name = group, coef_name = vars,
              coefficient = coefficients, cov_model = cov_model)
  if(!is.null(by_name)) out$by_name <- by_name
  out
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
    c(term$group_name, term$coef_name, term$by_name)
  }), use.names = FALSE))
  iid_vars <- iid_vars[nzchar(iid_vars)]

  iid_formula <- stats::reformulate(iid_vars)
  environment(iid_formula) <- environment(iid_info$reduced_formula)
  iid_mf <- stats::model.frame(iid_formula, data, na.action = stats::na.fail)

  if(nrow(iid_mf) != n_obs)
    stop("error: iid random-effect data rows do not match fixed-effect model rows")

  slope_names <- unique(unlist(lapply(terms, `[[`, "coef_name"), use.names = FALSE))
  slope_values <- stats::setNames(lapply(slope_names, function(nm)
    iid_numeric_covariate(iid_mf[[nm]], n_obs, paste0("iid slope covariate ", nm))), slope_names)

  # Block contract: name identifies a sampler block; formula_name, when present,
  # identifies its original control-list term. coefficient lists the block's
  # coefficients, while formula_coefficients retains an independent formula's
  # full coefficient list for control validation. rows/group_index/levels map
  # observations to contiguous group-major columns. Public term descriptions
  # expose formula_name as formula_term and coefficient as coefficients.
  # Independent coefficient formulas use the same scalar blocks as separate
  # iid intercept/slope terms, in coefficient order before parent expansion.
  scalar_terms <- list()
  for(term in terms){
    if(!identical(term$cov_model, "independent")){
      scalar_terms[[length(scalar_terms) + 1L]] <- term
      next
    }
    for(coef in term$coefficient){
      block <- term
      block$formula_name <- term$name
      block$formula_coefficients <- term$coefficient
      block$coefficient <- coef
      block$coef_name <- if(coef == "(Intercept)") NULL else coef
      block$name <- paste0(term$name, "[coefficient=", encodeString(coef, quote = '"'), "]")
      scalar_terms[[length(scalar_terms) + 1L]] <- block
    }
  }
  terms <- scalar_terms

  # Each observed parent is an existing sampler block. Its children occupy
  # contiguous columns; inactive observations never enter that block's Z.
  # Keeping parent and child labels separate avoids ambiguous pasted keys.
  blocks <- list()
  for(term in terms){
    if(is.null(term$by_name)){
      blocks[[length(blocks) + 1L]] <- term
      next
    }
    parent <- iid_mf[[term$by_name]]
    if(!is.atomic(parent) || !is.null(dim(parent)) || length(parent) != n_obs || anyNA(parent))
      stop("error: iid() by variable must be a nonmissing vector")
    parent <- droplevels(factor(parent))
    parent_levels <- levels(parent)
    if(!length(parent_levels)) stop("error: iid() by variable has no observed levels")
    for(level in parent_levels){
      block <- term
      block$formula_name <- term$formula_name %||% term$name
      block$by_level <- level
      block$by_levels <- parent_levels
      block$rows <- which(parent == level)
      block$name <- paste0(term$name, "[", term$by_name, "=", encodeString(level, quote = '"'), "]")
      blocks[[length(blocks) + 1L]] <- block
    }
  }
  terms <- blocks

  for(i in seq_along(terms)){
    g <- terms[[i]]$group_name
    if(!g %in% names(iid_mf))
      stop("error: grouping factor ", g, " not found for ", terms[[i]]$label)

    rows <- terms[[i]]$rows
    if(is.null(rows)) rows <- seq_len(n_obs)
    group <- iid_mf[[g]][rows]
    if(!is.null(terms[[i]]$by_name) && is.factor(group)) group <- droplevels(group)
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

    if(identical(terms[[i]]$cov_model, "unstructured")){
      values <- lapply(terms[[i]]$coefficient, function(nm){
        if(nm == "(Intercept)") return(rep(1, length(rows)))
        slope_values[[nm]][rows]
      })
      terms[[i]]$scale <- do.call(cbind, values)
    } else if(!is.null(terms[[i]]$coef_name)){
      terms[[i]]$scale <- slope_values[[terms[[i]]$coef_name]][rows]
    } else {
      terms[[i]]$scale <- rep(1.0, length(rows))
    }

    terms[[i]]$group_index <- as.integer(group)
  }

  re_nlevels <- vapply(terms, `[[`, integer(1), "n_levels")
  re_ncoef <- vapply(terms, function(term) length(term$coefficient), integer(1))
  re.q <- re_nlevels * re_ncoef
  q <- sum(re.q)
  n_value <- sum(vapply(terms, function(term)
    length(term$group_index) * length(term$coefficient), integer(1)))
  ii <- integer(n_value)
  jj <- integer(n_value)
  xx <- numeric(n_value)
  pos <- 0L
  col_offset <- 0L

  for(i in seq_along(terms)){
    term <- terms[[i]]
    idx <- term$group_index
    rows <- term$rows
    if(is.null(rows)) rows <- seq_len(n_obs)
    scale <- matrix(term$scale, nrow = length(rows), ncol = re_ncoef[i])

    for(row in seq_along(rows)){
      for(j in seq_len(re_ncoef[i])){
        pos <- pos + 1L
        ii[pos] <- rows[row]
        jj[pos] <- col_offset + (idx[row] - 1L) * re_ncoef[i] + j
        xx[pos] <- scale[row, j]
      }
    }

    col_offset <- col_offset + re.q[i]
  }

  Z <- Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(n_obs, q))

  z_names <- character(q)
  col_offset <- 0L
  for(i in seq_along(terms)){
    idx <- col_offset + seq_len(re.q[i])
    level_names <- paste0(terms[[i]]$name, "_", make.names(terms[[i]]$levels, unique = TRUE))
    z_names[idx] <- if(identical(terms[[i]]$cov_model, "unstructured"))
      paste0(rep(level_names, each = re_ncoef[i]), "_",
             rep(make.names(terms[[i]]$coefficient, unique = TRUE), times = re_nlevels[i])) else level_names
    z_names[idx] <- make.unique(z_names[idx])
    col_offset <- col_offset + re.q[i]
  }
  colnames(Z) <- z_names

  list(
    Z = Z,
    re_names = vapply(terms, `[[`, character(1), "name"),
    re_nlevels = as.integer(re_nlevels),
    re_ncoef = re_ncoef,
    re.q = as.integer(re.q),
    re_block_id = rep.int(seq_along(terms), re.q),
    terms = terms
  )
}
