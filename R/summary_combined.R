# Diagnostics accept iterations by chains for one scalar quantity. Derived
# quantities are constructed before this function is called; chains stay separate.
summary_draw_diagnostics <- function(draws, method){
  if(method == "rank"){
    if(!requireNamespace("posterior", quietly = TRUE))
      stop('error: diagnostic_method = "rank" requires the posterior package')
    if(nrow(draws) < 4L || any(!is.finite(draws)) ||
       any(vapply(seq_len(ncol(draws)), function(i)
         all(draws[, i] == draws[1L, i]), logical(1))))
      return(c(rhat = NA_real_, ess_bulk = NA_real_, ess_tail = NA_real_))
    return(c(rhat = if(ncol(draws) > 1L) posterior::rhat(draws) else NA_real_,
             ess_bulk = posterior::ess_bulk(draws),
             ess_tail = posterior::ess_tail(draws)))
  }
  chains <- do.call(coda::mcmc.list, lapply(seq_len(ncol(draws)), function(i)
    coda::mcmc(matrix(draws[, i], ncol = 1L))))
  rhat <- if(ncol(draws) > 1L && nrow(draws) > 1L)
    tryCatch(coda::gelman.diag(chains, autoburnin = FALSE,
      multivariate = FALSE)$psrf[1L, "Point est."], error = function(e) NA_real_) else NA_real_
  ess <- tryCatch(unname(coda::effectiveSize(chains)[1L]), error = function(e) NA_real_)
  c(rhat = rhat, effective_size = ess)
}

summary_matrix_diagnostics <- function(chains, method, parameters = NULL){
  if(inherits(chains, "mcmc")) chains <- list(chains)
  available <- colnames(chains[[1L]])
  parameters <- validate_summary_parameters(parameters, available) %||% available
  ans <- lapply(parameters, function(nm){
    draws <- vapply(chains, function(ch) as.numeric(ch[, nm]), numeric(nrow(chains[[1L]])))
    # vapply simplifies to a vector when only one iteration is retained.
    draws <- matrix(draws, nrow = nrow(chains[[1L]]), ncol = length(chains))
    summary_draw_diagnostics(draws, method)
  })
  values <- do.call(rbind, ans)
  data.frame(parameter = parameters, values, row.names = parameters, check.names = FALSE)
}

validate_combined_only_args <- function(coefficients, terms){
  if(!is.null(coefficients) || !is.null(terms))
    stop('error: coefficients and terms require type = "combined"')
}

combined_names <- function(x, what){
  if(!is.character(x) || !length(x) || anyNA(x) || any(!nzchar(x)) || anyDuplicated(x))
    stop("error: ", what, " must be a nonempty character vector without duplicates")
  x
}

# Read membership from the stored sparse design, including explicit zero entries
# for zero-valued slope covariates. Never divide Z by covariates or densify Z.
combined_iid_membership <- function(object, block, offset){
  Z <- object$backend$Z
  if(!inherits(Z, "dgCMatrix"))
    stop("error: combined summaries require the stored sparse iid design")
  membership <- integer(nrow(Z))
  ncoef <- length(block$coefficients)
  for(g in seq_along(block$levels)){
    col <- offset + (g - 1L) * ncoef + 1L
    first <- Z@p[col] + 1L
    last <- Z@p[col + 1L]
    if(first <= last){
      rows <- Z@i[seq.int(first, last)] + 1L
      if(any(membership[rows] != 0L))
        stop("error: stored iid design has ambiguous group membership")
      membership[rows] <- g
    }
  }
  membership
}

combined_iid_plan <- function(object, coefficients, terms){
  blocks <- object$term_description$random_effects$terms
  ids <- vapply(blocks, function(b) b$formula_term %||% b$name, character(1))
  available <- unique(ids)
  process <- vapply(object$backend$process_terms, `[[`, character(1), "name")
  if(is.null(terms))
    stop('error: type = "combined" requires explicit terms, or terms = "all"')
  terms <- combined_names(terms, "terms")
  if(identical(terms, "all")){
    if(length(process))
      stop('error: terms = "all" includes structured effects; combined summaries require explicit iid term selection in this model')
    terms <- available
  }
  if(any(terms %in% process))
    stop("error: combined summaries accept iid terms, not structured effects; select iid terms")
  unknown <- setdiff(terms, available)
  if(length(unknown))
    stop("error: unknown iid formula term(s): ", paste(unknown, collapse = ", "),
         "; available: ", paste(available, collapse = ", "))
  if(!length(terms)) stop("error: no iid terms available to combine")
  selected <- which(ids %in% terms)
  choices <- unique(unlist(lapply(blocks[selected], `[[`, "coefficients"), use.names = FALSE))
  if(is.null(coefficients)) coefficients <- choices
  coefficients <- combined_names(coefficients, "coefficients")
  unknown <- setdiff(coefficients, choices)
  if(length(unknown))
    stop("error: no selected iid term contains coefficient(s): ", paste(unknown, collapse = ", "))

  sizes <- vapply(blocks, function(b) length(b$levels) * length(b$coefficients), integer(1))
  offsets <- cumsum(c(0L, sizes))[seq_along(sizes)]
  if(sum(sizes) != ncol(object$backend$Z) || sum(sizes) != ncol(object$alpha_samples))
    stop("error: iid term metadata do not match stored coefficient draws")
  membership <- lapply(selected, function(i) combined_iid_membership(object, blocks[[i]], offsets[i]))
  names(membership) <- as.character(selected)
  n <- nrow(object$backend$Z)
  group_vars <- unique(unlist(lapply(blocks[selected], function(b) c(b$by, b$grouping_factor))))
  beta_names <- colnames(object$beta_samples)
  quantities <- list()
  alpha_columns <- list()

  for(coef in coefficients){
    contributors <- terms[vapply(terms, function(id)
      any(vapply(blocks[ids == id], function(b) coef %in% b$coefficients, logical(1))), logical(1))]
    columns <- matrix(0L, nrow = n, ncol = length(contributors))
    groups <- as.data.frame(setNames(rep(list(rep(NA_character_, n)), length(group_vars)),
                                     paste0("group_", group_vars)), stringsAsFactors = FALSE)
    for(k in seq_along(contributors)){
      for(i in selected[ids[selected] == contributors[k]]){
        b <- blocks[[i]]
        j <- match(coef, b$coefficients)
        if(is.na(j)) next
        member <- membership[[as.character(i)]]
        rows <- which(member > 0L)
        if(any(columns[rows, k] != 0L))
          stop("error: iid parent blocks overlap in the stored design")
        columns[rows, k] <- offsets[i] + (member[rows] - 1L) * length(b$coefficients) + j
        labels <- setNames(list(b$levels[member[rows]]), b$grouping_factor)
        if(!is.null(b$by)) labels[[b$by]] <- rep(b$by_level, length(rows))
        for(v in names(labels)){
          old <- groups[[paste0("group_", v)]][rows]
          if(any(!is.na(old) & old != labels[[v]]))
            stop("error: selected iid terms have inconsistent grouping labels")
          groups[[paste0("group_", v)]][rows] <- labels[[v]]
        }
      }
    }
    if(any(columns == 0L))
      stop("error: stored iid design is missing group membership; use the original fit with its explicit sparse entries")
    # Integer column tuples define represented combinations, with no pasted-label
    # collisions and no Cartesian expansion of crossed or nested groups.
    rows <- which(!duplicated(as.data.frame(columns)))
    for(row in rows){
      idx <- columns[row, ]
      fixed <- coef %in% beta_names
      expression <- c(if(fixed) paste0("beta[", encodeString(coef, quote = '"'), "]") else "0",
                      paste0("alpha[", encodeString(colnames(object$alpha_samples)[idx], quote = '"'), "]"))
      quantities[[length(quantities) + 1L]] <- cbind(
        data.frame(coefficient = coef, fixed_effect = fixed,
                   terms = paste(contributors, collapse = " + "), data_row = row,
                   expression = paste(expression, collapse = " + "), stringsAsFactors = FALSE),
        groups[row, , drop = FALSE])
      alpha_columns[[length(alpha_columns) + 1L]] <- as.integer(idx)
    }
  }
  info <- do.call(rbind, quantities)
  rownames(info) <- paste0("combined_", seq_len(nrow(info)))
  list(info = info, alpha_columns = alpha_columns, terms = terms, coefficients = coefficients)
}

summary_iid_combined <- function(object, probs, burn, coefficients, terms,
                                 diagnostics, diagnostic_method, parameters, include_w = FALSE){
  if(!is.null(parameters))
    stop('error: use coefficients and terms, not parameters, with type = "combined"')
  if(isTRUE(include_w))
    stop('error: include_w does not select contributions to combined coefficients; use terms')
  chains <- if(inherits(object, "stLMM_chains") || inherits(object, "stLMM_recovery_chains"))
    object$chains else list(object)
  reference <- chains[[1L]]
  burn <- validate_burn(burn, reference$backend$n_samples)
  plan <- combined_iid_plan(reference, coefficients, terms)
  for(k in seq_along(chains)){
    ch <- chains[[k]]
    # Starting values may differ; compare the actual coefficient/group maps.
    if(k > 1L && (!identical(combined_iid_plan(ch, coefficients, terms), plan) ||
                  ch$backend$n_samples != reference$backend$n_samples))
      stop("error: combined summaries require matching iid groups, coefficients, and iterations across chains")
    if(!identical(colnames(ch$beta_samples), colnames(reference$beta_samples)) ||
       !identical(colnames(ch$alpha_samples), colnames(reference$alpha_samples)) ||
       nrow(ch$beta_samples) != reference$backend$n_samples ||
       nrow(ch$alpha_samples) != reference$backend$n_samples)
      stop("error: combined coefficient samples are not aligned across chains")
  }
  iteration <- seq.int(burn + 1L, reference$backend$n_samples)
  summaries <- vector("list", nrow(plan$info))
  diags <- if(isTRUE(diagnostics)) vector("list", nrow(plan$info)) else NULL
  for(i in seq_len(nrow(plan$info))){
    draws <- vapply(chains, function(ch){
      value <- if(plan$info$fixed_effect[i]) ch$beta_samples[iteration, plan$info$coefficient[i]] else rep(0, length(iteration))
      for(col in plan$alpha_columns[[i]]) value <- value + ch$alpha_samples[iteration, col]
      value
    }, numeric(length(iteration)))
    draws <- matrix(draws, nrow = length(iteration), ncol = length(chains))
    summaries[[i]] <- posterior_summary_table(as.numeric(draws), probs = probs)
    if(isTRUE(diagnostics)) diags[[i]] <- summary_draw_diagnostics(draws, diagnostic_method)
  }
  tab <- do.call(rbind, summaries)
  rownames(tab) <- rownames(plan$info)
  if(isTRUE(diagnostics)){
    diags <- data.frame(parameter = rownames(plan$info), do.call(rbind, diags),
                        row.names = rownames(plan$info), check.names = FALSE)
  }
  out <- list(call = reference$backend$formula, type = "combined", scale = "link coefficient",
              terms = plan$terms, coefficients = plan$coefficients,
              quantities = plan$info, parameters = tab, diagnostics = diags,
              diagnostic_method = if(isTRUE(diagnostics)) diagnostic_method else NULL,
              n_chains = length(chains), n_samples = reference$backend$n_samples,
              burn = burn, n_used = length(iteration), iteration = iteration,
              chain_ids = vapply(seq_along(chains), function(i) as.integer(chains[[i]]$chain %||% i), integer(1)))
  class(out) <- "summary_stLMM_combined"
  out
}

print.summary_stLMM_combined <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  cat("stLMM combined iid coefficient summary\n")
  cat("  scale: ", x$scale, "; covariate values are not multiplied in\n", sep = "")
  cat("  terms: ", paste(x$terms, collapse = ", "), "\n", sep = "")
  cat("  chains: ", x$n_chains, "; draws used per chain: ", x$n_used,
      " (burn = ", x$burn, ")\n", sep = "")
  cat("  fixed_effect indicates whether a matching fixed coefficient is included\n")
  if(!is.null(x$diagnostics))
    cat("  diagnostics: ", if(x$diagnostic_method == "rank")
      "rank-normalized split/folded R-hat, bulk ESS, tail ESS (posterior)" else
      "classical R-hat and effective size (coda)", "\n", sep = "")
  cat("  rows describe represented group combinations; no averaging over observations\n\n")
  info <- x$quantities[setdiff(names(x$quantities), c("data_row", "expression"))]
  tab <- cbind(info, as.data.frame(x$parameters))
  if(!is.null(x$diagnostics)) tab <- cbind(tab, x$diagnostics[setdiff(names(x$diagnostics), "parameter")])
  numeric <- vapply(tab, is.numeric, logical(1))
  tab[numeric] <- lapply(tab[numeric], round, digits = digits)
  print(tab)
  invisible(x)
}
