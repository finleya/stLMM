resid <- function(model = c("tau_sq", "constant", "fixed", "group", "scaled"),
                  group,
                  variance,
                  vhat,
                  n = NULL,
                  prior = c("default", "ig", "shannon"),
                  method,
                  shape = 4,
                  center = c("mean", "mode"),
                  shrinkage = 10,
                  kappa_log_prior = c(mean = 0, sd = 1),
                  tau0_sq_log_prior = NULL,
                  starting = NULL,
                  tuning = 0.1){
  model <- match.arg(model)
  if(model == "constant")
    model <- "tau_sq"

  if(model == "tau_sq"){
    return(structure(
      list(type = "global_tau", label = "resid()"),
      class = "stLMM_residual"
    ))
  }

  if(model == "fixed"){
    if(missing(variance) && missing(vhat))
      stop("error: resid(model = 'fixed') requires variance")
    variance_expr <- if(missing(variance)) substitute(vhat) else substitute(variance)
    label <- paste(deparse(variance_expr, width.cutoff = 500L), collapse = " ")
    return(structure(
      list(
        type = "fixed_variance",
        variance_expr = variance_expr,
        variance_label = label,
        env = parent.frame()
      ),
      class = "stLMM_residual"
    ))
  }

  if(model == "group"){
    if(missing(group))
      stop("error: resid(model = 'group') requires group")
    prior <- match.arg(prior)
    has_variance <- !missing(variance) || !missing(vhat)
    if(prior == "default" && !has_variance){
      group_expr <- substitute(group)
      return(structure(
        list(
          type = "group_variance",
          group_expr = group_expr,
          label = paste0(
            "resid(model = 'group', group = ",
            paste(deparse(group_expr, width.cutoff = 500L), collapse = " "),
            ")"
          ),
          env = parent.frame()
        ),
        class = "stLMM_residual"
      ))
    }
    if(prior == "default")
      prior <- "ig"
    if(!has_variance)
      stop("error: resid(model = 'group', prior = '", prior, "') requires variance")
    group_expr <- substitute(group)
    variance_expr <- if(missing(variance)) substitute(vhat) else substitute(variance)
    n_expr <- substitute(n)
    method_value <- match.arg(
      if(!missing(method)) method else if(prior == "shannon") "shannon" else "constant",
      c("constant", "shannon")
    )
    center <- match.arg(center)
    if(method_value == "constant"){
      if(!is.numeric(shape) || length(shape) != 1L || is.na(shape) || shape <= 1)
        stop("error: shape must be a numeric scalar greater than 1")
    }else if(!missing(shape)){
      stop("error: shape is only used with prior = 'ig'")
    }
    if(!is.numeric(tuning) || length(tuning) != 1L || is.na(tuning) || tuning < 0)
      stop("error: tuning must be a nonnegative numeric scalar")
    return(structure(
      list(
        type = "group_ig_variance",
        group_expr = group_expr,
        vhat_expr = variance_expr,
        n_expr = n_expr,
        method = method_value,
        shape = if(method_value == "constant") as.double(shape) else NA_real_,
        center = center,
        starting = starting,
        tuning = as.double(tuning),
        label = paste0(
          "resid(model = 'group', group = ",
          paste(deparse(group_expr, width.cutoff = 500L), collapse = " "),
          ", variance = ",
          paste(deparse(variance_expr, width.cutoff = 500L), collapse = " "),
          ")"
        ),
        env = parent.frame()
      ),
      class = "stLMM_residual"
    ))
  }

  if(model == "scaled"){
    if(missing(variance) && missing(vhat))
      stop("error: resid(model = 'scaled') requires variance")
    variance_expr <- if(missing(variance)) substitute(vhat) else substitute(variance)
    n_expr <- substitute(n)
    has_n <- !missing(n) && !identical(n_expr, quote(NULL))

    if(!is.numeric(kappa_log_prior) || length(kappa_log_prior) != 2L ||
       anyNA(kappa_log_prior) || !is.finite(kappa_log_prior[1]) ||
       !is.finite(kappa_log_prior[2]) || kappa_log_prior[2] <= 0)
      stop("error: kappa_log_prior must be c(mean, sd) with finite mean and positive sd")
    if(!is.null(tau0_sq_log_prior) &&
       (!is.numeric(tau0_sq_log_prior) || length(tau0_sq_log_prior) != 2L ||
        anyNA(tau0_sq_log_prior) || !is.finite(tau0_sq_log_prior[1]) ||
        !is.finite(tau0_sq_log_prior[2]) || tau0_sq_log_prior[2] <= 0))
      stop("error: tau0_sq_log_prior must be c(mean, sd) with finite mean and positive sd")
    if(!is.numeric(shrinkage) || length(shrinkage) != 1L || is.na(shrinkage) ||
       !is.finite(shrinkage) || shrinkage <= 0)
      stop("error: shrinkage must be a finite positive numeric scalar")
    if(!is.numeric(tuning) || anyNA(tuning) || any(!is.finite(tuning)) || any(tuning < 0) ||
       !(length(tuning) %in% c(1L, 2L)))
      stop("error: tuning must be a nonnegative numeric scalar or length-two vector")

    return(structure(
      list(
        type = "scaled_variance",
        vhat_expr = variance_expr,
        n_expr = n_expr,
        has_n = has_n,
        shrinkage = as.double(shrinkage),
        kappa_log_prior = as.double(kappa_log_prior),
        tau0_sq_log_prior = if(is.null(tau0_sq_log_prior)) NULL else as.double(tau0_sq_log_prior),
        starting = starting,
        tuning = as.double(tuning),
        label = paste0(
          "resid(model = 'scaled', variance = ",
          paste(deparse(variance_expr, width.cutoff = 500L), collapse = " "),
          if(has_n) paste0(", n = ", paste(deparse(n_expr, width.cutoff = 500L), collapse = " ")) else "",
          ")"
        ),
        env = parent.frame()
      ),
      class = "stLMM_residual"
    ))
  }
}

is_resid_term <- function(lbl){
  grepl("^resid\\(", lbl) || grepl(":resid\\(", lbl) || grepl("^resid\\([^)]*\\):", lbl)
}

build_resid_components <- function(formula, data){
  tt <- terms(formula, data = data, keep.order = TRUE)
  term_labels <- attr(tt, "term.labels")
  resid_idx <- which(vapply(term_labels, is_resid_term, logical(1)))
  resid_labels <- term_labels[resid_idx]

  if(!length(resid_labels))
    return(list(residual = NULL, residual_label = NULL, reduced_formula = formula))

  if(length(resid_labels) > 1L)
    stop("error: only one resid() term is allowed")
  if(grepl(":", resid_labels))
    stop("error: resid() defines the Gaussian residual variance model and cannot be used in interactions")

  eval_env <- list2env(list(resid = resid), parent = environment(formula))
  residual <- eval(parse(text = resid_labels), envir = eval_env)
  if(!inherits(residual, "stLMM_residual") || is.null(residual$type))
    stop("error: resid() did not produce a valid stLMM residual model")

  reduced_formula <- update(formula, paste(". ~ . -", resid_labels))
  list(residual = residual, residual_label = resid_labels, reduced_formula = reduced_formula)
}

eval_residual_expr <- function(expr, data, env){
  value <- eval(expr, data, env)
  if(is.character(value) && length(value) == 1L){
    if(is.data.frame(data) && value %in% names(data))
      value <- data[[value]]
    else if(is.environment(data) && exists(value, envir = data, inherits = FALSE))
      value <- get(value, envir = data, inherits = FALSE)
  }
  value
}

is_fixed_parameter <- function(x){
  inherits(x, "stLMM_fixed_parameter")
}

fixed_parameter_numeric <- function(x, where){
  if(is_fixed_parameter(x))
    x <- x$value
  if(!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x))
    stop("error: fixed value in ", where, " must be a finite numeric scalar")
  as.double(x)
}

expand_residual_group_values <- function(x, group_levels, default, where){
  n_group <- length(group_levels)
  fixed <- rep(FALSE, n_group)
  names(fixed) <- group_levels

  if(is.null(x))
    return(list(value = rep(default, n_group), fixed = fixed))

  if(is.list(x) && !is_fixed_parameter(x) && !is.null(x$tau_sq))
    x <- x$tau_sq

  if(is_fixed_parameter(x)){
    value <- rep(fixed_parameter_numeric(x, where), n_group)
    fixed[] <- TRUE
    return(list(value = value, fixed = fixed))
  }

  if(is.list(x) && !is_fixed_parameter(x)){
    if(is.null(names(x)) || any(!nzchar(names(x))))
      stop("error: ", where, " must be a scalar or named by residual group")
    missing_value <- setdiff(group_levels, names(x))
    extra_value <- setdiff(names(x), group_levels)
    if(length(missing_value))
      stop("error: missing residual group value(s): ", paste(missing_value, collapse = ", "))
    if(length(extra_value))
      stop("error: unknown residual group value(s): ", paste(extra_value, collapse = ", "))
    value <- numeric(n_group)
    for(i in seq_along(group_levels)){
      nm <- group_levels[i]
      fixed[i] <- is_fixed_parameter(x[[nm]])
      value[i] <- fixed_parameter_numeric(x[[nm]], paste0(where, "$", nm))
    }
    names(value) <- group_levels
    return(list(value = value, fixed = fixed))
  }

  if(is.atomic(x)){
    value_names <- names(x)
    value <- as.double(x)
    names(value) <- value_names
    if(length(value) == 1L)
      return(list(value = rep(value, n_group), fixed = fixed))
    if(is.null(value_names))
      stop("error: ", where, " must be a scalar or named by residual group")
    missing_value <- setdiff(group_levels, value_names)
    extra_value <- setdiff(value_names, group_levels)
    if(length(missing_value))
      stop("error: missing residual group value(s): ", paste(missing_value, collapse = ", "))
    if(length(extra_value))
      stop("error: unknown residual group value(s): ", paste(extra_value, collapse = ", "))
    value <- value[group_levels]
    return(list(value = value, fixed = fixed))
  }

  stop("error: ", where, " must be a scalar or named by residual group")
}

get_residual_control_block <- function(x){
  if(is.null(x))
    return(NULL)
  if("resid" %in% names(x))
    return(x[["resid"]])
  NULL
}

eval_residual_variance <- function(residual, data, n_expected, where){
  if(is.null(residual))
    return(NULL)
  if(!inherits(residual, "stLMM_residual") || is.null(residual$type))
    stop("error: residual model must be created by resid()")
  if(residual$type != "fixed_variance")
    stop("error: unsupported residual model type: ", residual$type)

  value <- eval_residual_expr(residual$variance_expr, data, residual$env)

  if(!is.numeric(value) || length(value) != n_expected)
    stop("error: fixed residual variance in ", where, " must be numeric with length ", n_expected)

  as.double(value)
}

validate_fixed_residual_variance <- function(variance, observed_index, where){
  variance_obs <- variance[observed_index]
  if(anyNA(variance_obs))
    stop("error: fixed residual variance contains missing values for observed responses in ", where)
  if(any(!is.finite(variance_obs)) || any(variance_obs <= 0))
    stop("error: fixed residual variance must be finite and positive for observed responses in ", where)
  invisible(TRUE)
}

build_scaled_residual_model <- function(residual, data, observed_index, n_expected, where){
  vhat <- eval_residual_expr(residual$vhat_expr, data, residual$env)
  if(!is.numeric(vhat) || length(vhat) != n_expected)
    stop("error: residual vhat in ", where, " must be numeric with length ", n_expected)

  vhat_obs <- as.double(vhat[observed_index])
  if(anyNA(vhat_obs) || any(!is.finite(vhat_obs)) || any(vhat_obs <= 0))
    stop("error: residual vhat must be finite and positive for observed responses in ", where)

  weight <- rep(1, n_expected)
  n_eff <- rep(NA_real_, n_expected)
  param_names <- "kappa"
  starting_default <- c(kappa = 1)
  prior_meanlog <- c(kappa = residual$kappa_log_prior[1])
  prior_sdlog <- c(kappa = residual$kappa_log_prior[2])

  if(residual$has_n){
    n_eff_value <- eval_residual_expr(residual$n_expr, data, residual$env)
    if(!is.numeric(n_eff_value) || length(n_eff_value) != n_expected)
      stop("error: residual n in ", where, " must be numeric with length ", n_expected)
    n_obs <- as.double(n_eff_value[observed_index])
    if(anyNA(n_obs) || any(!is.finite(n_obs)) || any(n_obs <= 1))
      stop("error: residual n must be finite and greater than 1 for observed responses in ", where)
    n_eff <- as.double(n_eff_value)
    weight <- n_eff / (n_eff + residual$shrinkage)

    tau0_start <- stats::median(vhat_obs)
    tau0_prior <- residual$tau0_sq_log_prior
    if(is.null(tau0_prior))
      tau0_prior <- c(log(tau0_start), 1)
    param_names <- c("kappa", "tau0_sq")
    starting_default <- c(kappa = 1, tau0_sq = tau0_start)
    prior_meanlog <- c(kappa = residual$kappa_log_prior[1], tau0_sq = tau0_prior[1])
    prior_sdlog <- c(kappa = residual$kappa_log_prior[2], tau0_sq = tau0_prior[2])
  }

  starting <- residual$starting
  if(is.null(starting)){
    starting <- starting_default
  } else {
    starting <- as.double(starting)
    if(is.null(names(starting))){
      if(length(starting) == 1L && length(param_names) == 2L)
        starting <- c(starting, starting_default["tau0_sq"])
      names(starting) <- param_names[seq_along(starting)]
    }
    if(!all(param_names %in% names(starting)) || any(!names(starting) %in% param_names))
      stop("error: residual scaled starting values must be named ", paste(param_names, collapse = ", "))
    starting <- starting[param_names]
  }
  if(length(starting) != length(param_names) || any(!is.finite(starting)) || any(starting <= 0))
    stop("error: residual scaled starting values must be finite and positive")

  tuning <- residual$tuning
  if(length(tuning) == 1L)
    tuning <- rep(tuning, length(param_names))
  if(length(tuning) != length(param_names))
    stop("error: residual scaled tuning must have length 1 or ", length(param_names))
  names(tuning) <- param_names

  list(
    type = "scaled_variance",
    label = residual$label,
    vhat_expr = residual$vhat_expr,
    n_expr = residual$n_expr,
    has_n = residual$has_n,
    shrinkage = residual$shrinkage,
    parameter_names = param_names,
    starting = as.double(starting),
    tuning = as.double(tuning),
    prior_meanlog = as.double(prior_meanlog),
    prior_sdlog = as.double(prior_sdlog),
    vhat = as.double(vhat),
    vhat_obs = as.double(vhat_obs),
    n = as.double(n_eff),
    weight = as.double(weight),
    weight_obs = as.double(weight[observed_index])
  )
}

build_group_ig_residual_model <- function(residual, data, observed_index, n_expected, where){
  group <- eval_residual_expr(residual$group_expr, data, residual$env)
  vhat <- eval_residual_expr(residual$vhat_expr, data, residual$env)
  if(!is.null(residual$n_expr) && !identical(residual$n_expr, quote(NULL)))
    n_eff <- eval_residual_expr(residual$n_expr, data, residual$env)
  else
    n_eff <- NULL

  if(length(group) != n_expected)
    stop("error: residual group in ", where, " must have length ", n_expected)
  if(!is.numeric(vhat) || length(vhat) != n_expected)
    stop("error: residual vhat in ", where, " must be numeric with length ", n_expected)
  if(!is.null(n_eff) && (!is.numeric(n_eff) || length(n_eff) != n_expected))
    stop("error: residual n in ", where, " must be numeric with length ", n_expected)

  group_chr <- as.character(group)
  if(anyNA(group_chr[observed_index]))
    stop("error: residual group contains missing values for observed responses in ", where)

  group_levels <- unique(group_chr[!is.na(group_chr)])
  group_index_full <- match(group_chr, group_levels)
  group_index_obs <- as.integer(group_index_full[observed_index])
  n_group <- length(group_levels)

  vhat_group <- rep(NA_real_, n_group)
  n_group_eff <- rep(NA_real_, n_group)
  shape <- rep(NA_real_, n_group)
  scale <- rep(NA_real_, n_group)

  for(g in seq_len(n_group)){
    idx_obs <- observed_index[group_index_obs == g]
    idx_full <- which(group_index_full == g)
    v_obs <- vhat[idx_obs]
    v_all <- vhat[idx_full]
    v_good <- unique(v_obs[is.finite(v_obs) & v_obs > 0])

    if(length(idx_obs)){
      if(length(v_good) == 0L)
        stop("error: residual vhat is missing or nonpositive for observed group ", group_levels[g])
      if(length(v_good) > 1L)
        stop("error: residual vhat must be constant within observed group ", group_levels[g])
      vhat_group[g] <- v_good
    } else {
      v_full_good <- unique(v_all[is.finite(v_all) & v_all > 0])
      if(length(v_full_good) > 1L)
        stop("error: residual vhat must be constant within group ", group_levels[g])
      if(length(v_full_good) == 1L)
        vhat_group[g] <- v_full_good
    }

    if(!is.null(n_eff)){
      n_vals <- unique(n_eff[idx_full][is.finite(n_eff[idx_full]) & n_eff[idx_full] > 0])
      if(length(n_vals) > 1L)
        stop("error: residual n must be constant within group ", group_levels[g])
      if(length(n_vals) == 1L)
        n_group_eff[g] <- n_vals
    }

    if(is.na(vhat_group[g])){
      shape[g] <- NA_real_
      scale[g] <- NA_real_
    } else if(residual$method == "constant"){
      shape[g] <- residual$shape
      scale[g] <- if(residual$center == "mean")
        (shape[g] - 1) * vhat_group[g]
      else
        (shape[g] + 1) * vhat_group[g]
    } else {
      if(is.na(n_group_eff[g]) || n_group_eff[g] <= 1)
        stop("error: Shannon residual prior requires n > 1 for group ", group_levels[g])
      shape[g] <- n_group_eff[g] / 2
      scale[g] <- (n_group_eff[g] - 1) * vhat_group[g] / 2
      if(!(shape[g] > 0 && scale[g] > 0))
        stop("error: Shannon residual prior produced invalid parameters for group ", group_levels[g])
    }
  }

  active_group <- which(!is.na(shape) & !is.na(scale))
  if(!all(unique(group_index_obs) %in% active_group))
    stop("error: every observed residual group must have a finite positive prior")

  starting <- residual$starting
  if(is.null(starting)){
    starting <- vhat_group
  } else if(length(starting) == 1L){
    starting <- rep(as.double(starting), n_group)
  } else {
    starting <- as.double(starting)
  }
  if(length(starting) != n_group)
    stop("error: residual starting values must have length 1 or the number of residual groups")
  starting[is.na(starting)] <- vhat_group[is.na(starting)]
  starting[is.na(starting)] <- 1

  tuning <- rep(residual$tuning, n_group)
  tuning[is.na(shape)] <- 0

  if(any(!is.finite(starting[active_group]) | starting[active_group] <= 0))
    stop("error: residual starting variances must be finite and positive")

  list(
    type = "group_ig_variance",
    label = residual$label,
    group_expr = residual$group_expr,
    vhat_expr = residual$vhat_expr,
    n_expr = residual$n_expr,
    method = residual$method,
    center = residual$center,
    groups = group_levels,
    group_index_full = as.integer(group_index_full),
    group_index_obs = as.integer(group_index_obs),
    starting = as.double(starting),
    tuning = as.double(tuning),
    shape = as.double(shape),
    scale = as.double(scale),
    vhat = as.double(vhat_group),
    n = as.double(n_group_eff)
  )
}

build_group_residual_model <- function(residual,
                                       data,
                                       observed_index,
                                       n_expected,
                                       where,
                                       starting,
                                       tuning,
                                       priors,
                                       validate_variance_prior,
                                       prior_family_code){
  group <- eval_residual_expr(residual$group_expr, data, residual$env)

  if(length(group) != n_expected)
    stop("error: residual group in ", where, " must have length ", n_expected)

  group_chr <- as.character(group)
  if(anyNA(group_chr[observed_index]))
    stop("error: residual group contains missing values for observed responses in ", where)

  group_levels <- unique(group_chr[!is.na(group_chr)])
  group_index_full <- match(group_chr, group_levels)
  group_index_obs <- as.integer(group_index_full[observed_index])
  n_group <- length(group_levels)

  starting_entry <- get_residual_control_block(starting)
  starting_expanded <- expand_residual_group_values(
    starting_entry, group_levels, default = 1,
    where = "starting$resid"
  )
  starting_value <- starting_expanded$value
  fixed_value <- starting_expanded$fixed

  tuning_entry <- get_residual_control_block(tuning)
  if(is.null(tuning_entry)){
    tuning_value <- rep(0.1, n_group)
    names(tuning_value) <- group_levels
  } else {
    tuning_expanded <- expand_residual_group_values(
      tuning_entry, group_levels, default = 0.1,
      where = "tuning$resid"
    )
    tuning_value <- tuning_expanded$value
  }
  tuning_value[fixed_value] <- 0

  if(!is.null(tuning_entry)){
    tuning_check <- expand_residual_group_values(
      tuning_entry, group_levels, default = 0.1,
      where = "tuning$resid"
    )$value
    if(any(fixed_value & tuning_check > 0))
      stop("error: fixed residual group variance cannot have positive tuning")
  }

  if(length(starting_value) != n_group || any(!is.finite(starting_value)) || any(starting_value <= 0))
    stop("error: residual group starting variances must be finite and positive")
  if(length(tuning_value) != n_group || any(!is.finite(tuning_value)) || any(tuning_value < 0))
    stop("error: residual group tuning values must be finite and nonnegative")

  free_group <- group_levels[tuning_value > 0]
  residual_prior <- get_residual_control_block(priors)
  if(is.null(residual_prior) && length(free_group))
    stop("error: missing residual variance prior. Use priors = list(resid = list(tau_sq = ig(shape, scale)))")

  prior_entries <- vector("list", n_group)
  names(prior_entries) <- group_levels

  if(is.null(residual_prior)){
    prior_entries[] <- list(NULL)
  } else if(inherits(residual_prior, "stLMM_prior")){
    prior_entries[] <- list(residual_prior)
  } else if(is.list(residual_prior) && inherits(residual_prior$tau_sq, "stLMM_prior")){
    prior_entries[] <- list(residual_prior$tau_sq)
  } else if(is.list(residual_prior)){
    missing_prior <- setdiff(free_group, names(residual_prior))
    extra_prior <- setdiff(names(residual_prior), group_levels)
    if(length(missing_prior))
      stop("error: missing residual group prior(s): ", paste(missing_prior, collapse = ", "))
    if(length(extra_prior))
      stop("error: unknown residual group prior(s): ", paste(extra_prior, collapse = ", "))
    prior_entries[names(residual_prior)] <- residual_prior
  } else {
    stop("error: priors$resid must be a prior, list(tau_sq = <prior>), or a named list by residual group")
  }

  prior <- matrix(NA_real_, n_group, 6,
                  dimnames = list(group_levels, c("family", "p1", "p2", "lower", "upper", "scale")))
  shape <- rep(NA_real_, n_group)
  scale <- rep(NA_real_, n_group)
  for(g in seq_len(n_group)){
    if(is.null(prior_entries[[g]])){
      prior[g, ] <- c(prior_family_code[["ig"]], 1, 1, 0, 0, 0)
    } else {
      prior[g, ] <- validate_variance_prior(
        prior_entries[[g]],
        paste0("priors$resid$", group_levels[g])
      )
    }
    if(prior[g, 1] == prior_family_code[["ig"]]){
      shape[g] <- prior[g, 2]
      scale[g] <- prior[g, 3]
    }
  }

  list(
    type = "group_ig_variance",
    label = residual$label,
    group_expr = residual$group_expr,
    method = "prior",
    center = NA_character_,
    groups = group_levels,
    group_index_full = as.integer(group_index_full),
    group_index_obs = as.integer(group_index_obs),
    starting = as.double(starting_value),
    tuning = as.double(tuning_value),
    shape = as.double(shape),
    scale = as.double(scale),
    prior = unname(prior),
    vhat = rep(NA_real_, n_group),
    n = rep(NA_real_, n_group)
  )
}
