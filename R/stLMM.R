stLMM <- function(formula, data = parent.frame(),
                 starting = list(), tuning = NULL, priors = NULL,
                 family = "gaussian",
                 trials = NULL,
                 size = NULL,
                 n_samples,
                 n_omp_threads = 1,
                 nngp_search = c("fast", "brute"),
                 verbose = TRUE,
                 n_report = 100,
                 warmup = TRUE,
                 metropolis = list(),
                 cholmod_control = list(),
                 save_process = NULL,
                 chains = 1L,
                 chain_control = list(),
                 describe_terms = FALSE,
                 ...){
    
    ## stLMM() is intentionally the R-side orchestration point. The flow below
    ## is linear: parse formula terms, build model matrices, resolve controls,
    ## construct the backend contract consumed by C++, run the sampler, then
    ## label and package returned samples.

    ####################################################
    ## Helpers
    ####################################################

    `%||%` <- function(x, y){
        if(is.null(x)) y else x
    }

    require_positive_pair <- function(x, name){
        if(!is.numeric(x) || length(x) != 2)
            stop("error: ", name, " must be numeric length 2")
        if(any(is.na(x)) || any(x <= 0))
            stop("error: ", name, " must be positive")
        x
    }
    
    require_positive_scalar <- function(x, name, where){
        val <- x[1]
        if(!is.numeric(val) || length(val) != 1 || is.na(val) || val <= 0)
            stop("error: ", name, " must be positive in ", where)
        val
    }

    require_nonnegative_scalar <- function(x, name, where){
        val <- x[1]
        if(!is.numeric(val) || length(val) != 1 || is.na(val) || val < 0)
            stop("error: ", name, " must be non-negative in ", where)
        val
    }

    normalize_save_process <- function(x, is_pg_likelihood, has_process){
        default_enabled <- isTRUE(is_pg_likelihood) && isTRUE(has_process)
        if(is.null(x)){
            return(list(
                enabled = default_enabled,
                start = 1L,
                thin = 1L,
                defaulted = TRUE
            ))
        }

        if(is.logical(x) && length(x) == 1L && !is.na(x)){
            enabled <- isTRUE(x)
            if(enabled && (!isTRUE(is_pg_likelihood) || !isTRUE(has_process)))
                stop("error: save_process can only be enabled for Polya-Gamma models with structured process terms")
            return(list(
                enabled = enabled,
                start = 1L,
                thin = 1L,
                defaulted = FALSE
            ))
        }

        if(!is.list(x))
            stop("error: save_process must be NULL, TRUE/FALSE, or a list with optional entries 'start' and 'thin'")

        enabled <- x$enabled %||% x$process %||% TRUE
        if(!is.logical(enabled) || length(enabled) != 1L || is.na(enabled))
            stop("error: save_process$enabled must be TRUE or FALSE")
        if(isTRUE(enabled) && (!isTRUE(is_pg_likelihood) || !isTRUE(has_process)))
            stop("error: save_process can only be enabled for Polya-Gamma models with structured process terms")

        start <- x$start %||% 1L
        thin <- x$thin %||% 1L
        if(!is.numeric(start) || length(start) != 1L || is.na(start) || start < 1)
            stop("error: save_process$start must be a positive integer")
        if(!is.numeric(thin) || length(thin) != 1L || is.na(thin) || thin < 1)
            stop("error: save_process$thin must be a positive integer")

        bad <- setdiff(names(x), c("enabled", "process", "start", "thin"))
        if(length(bad))
            stop("error: unsupported save_process field(s): ", paste(bad, collapse = ", "))

        list(
            enabled = isTRUE(enabled),
            start = as.integer(start),
            thin = as.integer(thin),
            defaulted = FALSE
        )
    }

    is_fixed_parameter <- function(x){
        inherits(x, "stLMM_fixed_parameter")
    }

    fixed_parameter_value <- function(x, name, where){
        if(!is_fixed_parameter(x))
            return(x)
        require_positive <- FALSE
        if(!is.null(name))
            require_positive <- name %in% c("tau_sq", "sigma_sq")
        val <- x$value
        if(!is.numeric(val) || length(val) != 1L || is.na(val) || !is.finite(val))
            stop("error: fixed value for ", name, " in ", where, " must be finite")
        if(require_positive && val <= 0)
            stop("error: fixed value for ", name, " in ", where, " must be positive")
        val
    }

    is_fixed_entry <- function(block, name){
        !is.null(block) && name %in% names(block) && is_fixed_parameter(block[[name]])
    }

    no_op_variance_prior <- function(){
        encode_prior(ig(1, 1))
    }

    no_op_theta_prior <- function(){
        encode_prior(uniform(0, 1))
    }

    default_positive_scale <- function(x){
        val <- stats::var(as.numeric(x))
        if(!is.finite(val) || val <= 0)
            val <- 1
        val
    }

    normalize_n_report <- function(x){
        if(is.null(x) || length(x) == 0)
            return(0L)
        val <- x[1]
        if(!is.numeric(val) || length(val) != 1 || is.na(val) || val < 0)
            stop("error: n_report must be a non-negative integer")
        val <- as.integer(val)
        if(val < 0L)
            stop("error: n_report must be a non-negative integer")
        val
    }

    normalize_warmup <- function(x, metropolis){
        default_target <- if(identical(metropolis$blocking, "scalar")) c(0.30, 0.60) else c(0.15, 0.45)
        defaults <- list(
            enabled = TRUE,
            batch_length = 25L,
            min_batches = 0L,
            max_batches = 20L,
            target = default_target,
            near_zero = 0.02
        )
        if(is.null(x))
            return(modifyList(defaults, list(enabled = FALSE)))
        if(is.logical(x) && length(x) == 1L)
            return(modifyList(defaults, list(enabled = isTRUE(x))))
        if(!is.list(x))
            stop("error: warmup must be TRUE, FALSE, NULL, or a list")

        out <- modifyList(defaults, x)
        if(!is.logical(out$enabled) || length(out$enabled) != 1L || is.na(out$enabled))
            stop("error: warmup$enabled must be TRUE or FALSE")
        out$batch_length <- require_positive_scalar(out$batch_length, "batch_length", "warmup")
        out$min_batches <- require_nonnegative_scalar(out$min_batches, "min_batches", "warmup")
        out$max_batches <- require_nonnegative_scalar(out$max_batches, "max_batches", "warmup")
        if(out$min_batches > out$max_batches)
            stop("error: warmup$min_batches must be less than or equal to warmup$max_batches")
        if(!is.numeric(out$target) || length(out$target) != 2L || any(is.na(out$target)) ||
           out$target[1] <= 0 || out$target[2] >= 1 || out$target[1] >= out$target[2])
            stop("error: warmup$target must be numeric length 2 with 0 < lower < upper < 1")
        out$near_zero <- require_nonnegative_scalar(out$near_zero, "near_zero", "warmup")
        if(out$near_zero >= out$target[1])
            stop("error: warmup$near_zero must be less than warmup$target[1]")
        out$batch_length <- as.integer(out$batch_length)
        out$min_batches <- as.integer(out$min_batches)
        out$max_batches <- as.integer(out$max_batches)
        out
    }

    normalize_metropolis <- function(x){
        if(is.null(x))
            x <- list()
        if(is.character(x) && length(x) == 1L)
            x <- list(blocking = x)
        if(!is.list(x))
            stop("error: metropolis must be a list or a blocking string")
        defaults <- list(blocking = "joint", target_accept = NULL, batch_length = 25L)
        out <- modifyList(defaults, x)
        if(!is.character(out$blocking) || length(out$blocking) != 1L || is.na(out$blocking))
            stop("error: metropolis$blocking must be a single string")
        choices <- c(
            joint = 0L,
            by_term = 1L,
            residual_process = 2L,
            variance_theta = 3L,
            process_variance = 4L,
            scalar = 5L
        )
        if(!out$blocking %in% names(choices))
            stop(
                "error: unknown metropolis$blocking; use one of ",
                paste(names(choices), collapse = ", ")
            )
        if(is.null(out$target_accept))
            out$target_accept <- if(identical(out$blocking, "scalar")) 0.44 else 0.234
        if(!is.numeric(out$target_accept) || length(out$target_accept) != 1L || is.na(out$target_accept))
            stop("error: metropolis$target_accept must be a single number")
        out$target_accept <- require_positive_scalar(out$target_accept, "target_accept", "metropolis")
        if(out$target_accept >= 1)
            stop("error: metropolis$target_accept must be less than 1")
        if(!is.numeric(out$batch_length) || length(out$batch_length) != 1L || is.na(out$batch_length))
            stop("error: metropolis$batch_length must be a single number")
        out$batch_length <- require_positive_scalar(out$batch_length, "batch_length", "metropolis")
        out$batch_length <- as.integer(out$batch_length)
        out$blocking_code <- unname(choices[[out$blocking]])
        out
    }

    normalize_chains <- function(x){
        if(!is.numeric(x) || length(x) != 1L || is.na(x) || x < 1 || abs(x - round(x)) > 0)
            stop("error: chains must be a positive integer")
        as.integer(x)
    }

    normalize_chain_control <- function(x){
        if(is.null(x))
            x <- list()
        if(!is.list(x))
            stop("error: chain_control must be a list")
        defaults <- list(
            seed = NULL,
            dispersion = 1
        )
        out <- modifyList(defaults, x)
        if(!is.null(out$seed)){
            if(!is.numeric(out$seed) || length(out$seed) != 1L || is.na(out$seed))
                stop("error: chain_control$seed must be a numeric scalar")
            out$seed <- as.integer(out$seed)
        }
        if(!is.numeric(out$dispersion) || length(out$dispersion) != 1L ||
           is.na(out$dispersion) || !is.finite(out$dispersion) || out$dispersion < 0)
            stop("error: chain_control$dispersion must be a nonnegative numeric scalar")
        out
    }

    normalize_nngp_search <- function(x){
        if(is.null(x))
            x <- "fast"
        match.arg(x, c("fast", "brute"))
    }

    normalize_cholmod_control <- function(x){
        if(is.null(x))
            x <- list()
        if(is.character(x) && length(x) == 1L)
            x <- list(ordering = x)
        if(!is.list(x))
            stop("error: cholmod_control must be a list or an ordering string")
        defaults <- list(
            ordering = "auto",
            postorder = TRUE
        )
        out <- modifyList(defaults, x)
        if(!is.character(out$ordering) || length(out$ordering) != 1L || is.na(out$ordering))
            stop("error: cholmod_control$ordering must be a single string")
        out$ordering <- match.arg(
            out$ordering,
            c("auto", "best", "natural", "amd", "metis", "nesdis", "colamd")
        )
        if(!is.logical(out$postorder) || length(out$postorder) != 1L || is.na(out$postorder))
            stop("error: cholmod_control$postorder must be TRUE or FALSE")
        out
    }

    normalize_family <- function(x){
        if(is.function(x))
            x <- x()
        if(is.character(x) && length(x) == 1L){
            fam <- tolower(x)
        } else if(is.list(x) && !is.null(x$family)) {
            fam <- tolower(as.character(x$family[1L]))
        } else {
            stop("error: family must be \"gaussian\", \"binomial\", or a corresponding stats family")
        }
        if(fam %in% c("gaussian", "normal"))
            return("gaussian")
        if(fam == "binomial")
            return("binomial")
        if(fam %in% c("negative_binomial", "negbin", "nb"))
            return("negative_binomial")
        stop("error: unsupported family '", fam, "'; currently use \"gaussian\", \"binomial\", or \"negative_binomial\"")
    }

    normalize_trials <- function(x, data, n_expected){
        if(is.null(x))
            return(rep.int(1L, n_expected))

        if(is.character(x) && length(x) == 1L){
            if(is.data.frame(data) && x %in% names(data)){
                val <- data[[x]]
            } else {
                val <- eval(as.name(x), envir = data, enclos = parent.frame())
            }
        } else {
            val <- x
        }

        if(length(val) != n_expected)
            stop("error: trials must have length matching the model frame")
        if(!is.numeric(val) && !is.integer(val))
            stop("error: trials must be numeric or integer")
        if(anyNA(val) || any(!is.finite(as.numeric(val))) || any(as.numeric(val) <= 0))
            stop("error: trials must be finite positive counts")
        if(any(abs(as.numeric(val) - round(as.numeric(val))) > sqrt(.Machine$double.eps)))
            stop("error: trials must be integer-valued counts")
        as.integer(round(as.numeric(val)))
    }

    normalize_nb_size <- function(x){
        if(is.null(x))
            stop("error: size must be supplied for family = \"negative_binomial\"")
        if(is_fixed_parameter(x))
            x <- fixed_parameter_value(x, "size", "size")
        val <- x[1]
        if(!is.numeric(val) || length(val) != 1L || is.na(val) ||
           !is.finite(val) || val <= 0)
            stop("error: size must be a finite positive scalar for family = \"negative_binomial\"")
        as.double(val)
    }

    chains <- normalize_chains(chains)
    chain_control <- normalize_chain_control(chain_control)
    nngp_search <- normalize_nngp_search(nngp_search)
    cholmod_control <- normalize_cholmod_control(cholmod_control)
    if(!is.null(chain_control$seed))
        set.seed(chain_control$seed)

    check_missing_extra <- function(x_names, required_names, where){
        
        missing_names <- setdiff(required_names, x_names)
        
        if(length(missing_names))
            stop(
                "error: missing parameter(s) ",
                paste(missing_names, collapse = ", "),
                " in ", where
            )
        
        extra_names <- setdiff(x_names, required_names)
        
        if(length(extra_names))
            stop(
                "error: unknown parameter(s) ",
                paste(extra_names, collapse = ", "),
                " in ", where
            )
    }
    
    normalize_iid_sigma_block <- function(x, where){
        if(inherits(x, "stLMM_prior") || is.atomic(x))
            stop("error: ", where, " must be a list with entry sigma_sq")
        block <- normalize_param_block(x, where)
        if(is.null(block) || !"sigma_sq" %in% names(block))
            stop("error: ", where, " must include sigma_sq")
        extra <- setdiff(names(block), "sigma_sq")
        if(length(extra))
            stop("error: unknown parameter(s) ", paste(extra, collapse = ", "), " in ", where)
        block$sigma_sq
    }

    collect_re_pairs <- function(x_list, re_names, where){
        
        out <- matrix(NA_real_, length(re_names), 2,
                      dimnames = list(re_names, c("shape", "scale")))
        
        for(i in seq_along(re_names)){
            xi <- x_list[[re_names[i]]]
            
            if(is.null(xi))
                stop("error: missing entry for ", re_names[i], " in ", where)
            
            xi <- normalize_iid_sigma_block(xi, paste0(where, "$", re_names[i]))
            xi <- validate_variance_prior(xi, paste0(where, "$", re_names[i], "$sigma_sq"), allow_only_ig = TRUE)
            out[i, ] <- xi[2:3]
        }
        
        out
    }
    
    collect_re_scalars <- function(x_list, re_names, where){
        
        out <- numeric(length(re_names))
        names(out) <- re_names
        
        for(i in seq_along(re_names)){
            xi <- x_list[[re_names[i]]]
            
            if(is.null(xi))
                stop("error: missing entry for ", re_names[i], " in ", where)
            
            xi <- normalize_iid_sigma_block(xi, paste0(where, "$", re_names[i]))
            out[i] <- require_positive_scalar(xi, "sigma_sq", paste0(where, "$", re_names[i]))
        }
        
        out
    }
    

    normalize_named_list <- function(x, where){
        if(is.null(x))
            return(list())
        if(!is.list(x))
            stop("error: ", where, " must be a list")
        if(length(x) == 0L)
            return(x)
        nms <- names(x)
        if(is.null(nms))
            return(x)
        names(x) <- tolower(nms)
        x
    }

    get_resid_control <- function(x, field = NULL){
        block <- NULL
        if("resid" %in% names(x))
            block <- x[["resid"]]
        if(is.null(field))
            return(block)
        if(is.list(block) && !is.null(block[[field]]))
            return(block[[field]])
        NULL
    }

    status_message <- function(...){
        if(isTRUE(verbose))
            message(...)
    }

    ## Convert a user-supplied parameter block into a named list whose
    ## parameter names are case-insensitive. This supports either a named
    ## numeric vector, e.g. c(sigma_sq = 1, phi = 0.4), or a named list,
    ## e.g. list(sigma_sq = 1, phi = 0.4).
    normalize_param_block <- function(x, where){
        if(is.null(x))
            return(NULL)

        if(is.atomic(x) && !is.list(x)){
            nms <- names(x)
            if(is.null(nms) || any(!nzchar(nms)))
                stop("error: ", where, " must be a named vector or named list")
            out <- as.list(x)
            names(out) <- tolower(nms)
            return(out)
        }

        if(is.list(x)){
            nms <- names(x)
            if(is.null(nms) || any(!nzchar(nms)))
                stop("error: ", where, " must have named entries")
            names(x) <- tolower(nms)
            return(x)
        }

        stop("error: ", where, " must be a named vector or named list")
    }

    prior_family_code <- c(
        ig = 0,
        uniform = 1,
        log_normal = 2,
        gamma_dist = 3,
        half_normal = 4,
        half_t = 5,
        beta_dist = 6,
        flat = 7,
        normal = 8
    )

    encode_prior <- function(pr){
        p <- rep(0, 6)
        p[1] <- prior_family_code[[pr$family]]
        if(length(pr$parameters) >= 1L)
            p[2] <- as.double(pr$parameters[1])
        if(length(pr$parameters) >= 2L)
            p[3] <- as.double(pr$parameters[2])
        if(!is.null(pr$support)){
            p[4] <- as.double(pr$support[1])
            p[5] <- as.double(pr$support[2])
        }
        p[6] <- if(identical(pr$scale, "sd")) 1 else 0
        p
    }

    legacy_ig_pair <- function(prior_desc){
        if(prior_desc[1] == prior_family_code[["ig"]])
            as.double(prior_desc[2:3])
        else
            c(2, 1)
    }

    require_prior_object <- function(x, where){
        if(!inherits(x, "stLMM_prior"))
            stop("error: ", where, " must use a prior constructor such as ig(), uniform(), log_normal(), gamma_dist(), half_normal(), half_t(), beta_dist(), flat(), or normal()")
        if(!is.character(x$family) || length(x$family) != 1L || is.na(x$family) ||
           !x$family %in% names(prior_family_code))
            stop("error: invalid prior family in ", where)
        if(is.null(x$parameters))
            x$parameters <- numeric(0)
        if(!identical(x$family, "normal")){
            x$parameters <- as.double(x$parameters)
            if(any(!is.finite(x$parameters)))
                stop("error: prior hyperparameters must be finite in ", where)
        }
        if(!is.null(x$support)){
            if(!is.numeric(x$support) || length(x$support) != 2L || any(!is.finite(x$support)) ||
               x$support[1] >= x$support[2])
                stop("error: prior support must be finite numeric length 2 with lower < upper in ", where)
            x$support <- c(lower = as.double(x$support[1]), upper = as.double(x$support[2]))
        }
        if(is.null(x$scale))
            x$scale <- "value"
        if(!x$scale %in% c("value", "sd"))
            stop("error: invalid prior scale in ", where)
        x
    }

    validate_prior_hyperparameters <- function(pr, where){
        fam <- pr$family
        par <- pr$parameters
        if(fam == "ig"){
            if(length(par) != 2L || any(par <= 0))
                stop("error: ig() requires positive shape and scale in ", where)
        } else if(fam == "uniform"){
            if(length(par) != 0L || is.null(pr$support))
                stop("error: uniform() prior is malformed in ", where)
        } else if(fam == "log_normal"){
            if(length(par) != 2L || par[2] <= 0)
                stop("error: log_normal() requires finite meanlog and positive sdlog in ", where)
        } else if(fam == "gamma_dist"){
            if(length(par) != 2L || any(par <= 0))
                stop("error: gamma_dist() requires positive shape and rate in ", where)
        } else if(fam == "half_normal"){
            if(length(par) != 1L || par[1] <= 0)
                stop("error: half_normal() requires positive scale in ", where)
            if(!identical(pr$scale, "sd"))
                stop("error: half_normal() is defined on the standard-deviation scale in ", where)
        } else if(fam == "half_t"){
            if(length(par) != 2L || any(par <= 0))
                stop("error: half_t() requires positive df and scale in ", where)
            if(!identical(pr$scale, "sd"))
                stop("error: half_t() is defined on the standard-deviation scale in ", where)
        } else if(fam == "beta_dist"){
            if(length(par) != 2L || any(par <= 0))
                stop("error: beta_dist() requires positive shape1 and shape2 in ", where)
        } else if(fam == "flat"){
            if(length(par) != 0L)
                stop("error: flat() prior is malformed in ", where)
        } else if(fam == "normal"){
            if(!is.list(pr$parameters) || !all(c("mean", "sd") %in% names(pr$parameters)))
                stop("error: normal() prior is malformed in ", where)
        }
        pr
    }

    validate_beta_prior <- function(x, p, where, default_family){
        if(p <= 0L)
            return(list(type = 0L, mean = numeric(0), precision = numeric(0)))

        if(is.null(x))
            x <- if(identical(default_family, "normal")) normal(mean = 0, sd = 10) else flat()

        if(!inherits(x, "stLMM_prior"))
            stop("error: ", where, " must use flat() or normal(mean, sd)")
        if(!is.character(x$family) || length(x$family) != 1L || is.na(x$family))
            stop("error: invalid beta prior family in ", where)

        if(identical(x$family, "flat")){
            if(length(x$parameters) != 0L)
                stop("error: flat() prior is malformed in ", where)
            return(list(type = 0L, mean = numeric(p), precision = numeric(p)))
        }

        if(!identical(x$family, "normal"))
            stop("error: ", where, " must use flat() or normal(mean, sd)")
        if(!is.list(x$parameters) || !all(c("mean", "sd") %in% names(x$parameters)))
            stop("error: normal() prior is malformed in ", where)

        mean <- as.numeric(x$parameters$mean)
        sd <- as.numeric(x$parameters$sd)
        if(length(mean) == 1L)
            mean <- rep(mean, p)
        if(length(sd) == 1L)
            sd <- rep(sd, p)
        if(length(mean) != p)
            stop("error: normal() beta prior mean must have length 1 or length matching the fixed effects in ", where)
        if(length(sd) != p)
            stop("error: normal() beta prior sd must have length 1 or length matching the fixed effects in ", where)
        if(any(!is.finite(mean)))
            stop("error: normal() beta prior mean must be finite in ", where)
        if(any(!is.finite(sd)) || any(sd <= 0))
            stop("error: normal() beta prior sd must be finite and positive in ", where)

        list(type = 1L, mean = as.double(mean), precision = as.double(1 / sd^2))
    }

    validate_variance_prior <- function(x, where, allow_only_ig = FALSE){
        pr <- validate_prior_hyperparameters(require_prior_object(x, where), where)
        allowed <- if(allow_only_ig) "ig" else c("ig", "log_normal", "gamma_dist", "half_normal", "half_t")
        if(!pr$family %in% allowed){
            if(allow_only_ig)
                stop("error: ", where, " must use ig(); grouped random-effect variances are currently Gibbs-updated with conjugate inverse-gamma priors")
            stop("error: ", where, " uses unsupported variance prior ", pr$family)
        }
        if(!is.null(pr$support))
            stop("error: variance priors do not currently accept finite support in ", where)
        encode_prior(pr)
    }

    theta_domain <- function(theta_name, theta_type, term){
        if(term$term_type %in% c("nngp", "gp") &&
           identical(term$cov_model, "gneiting") &&
           theta_name == "delta")
            return("nonnegative")
        if(term$term_type == "ar1")
            return("bounded")
        if(term$term_type %in% c("car", "dagar"))
            return("unit")
        if(term$term_type %in% c("car_time", "dagar_time")){
            if(theta_name == "rho")
                return("unit")
            if(theta_name == "phi")
                return("bounded")
            if(theta_name == "lambda")
                return("positive")
        }
        if(theta_type == 1L)
            return("positive")
        if(theta_type == 2L)
            return("unit")
        if(theta_type == 3L)
            return("bounded")
        "bounded"
    }

    validate_theta_prior <- function(x, theta_name, theta_type, term, where){
        pr <- validate_prior_hyperparameters(require_prior_object(x, where), where)
        domain <- theta_domain(theta_name, theta_type, term)

        allowed <- switch(
            domain,
            positive = c("uniform", "log_normal", "gamma_dist"),
            nonnegative = "uniform",
            unit = c("uniform", "beta_dist"),
            bounded = "uniform",
            "uniform"
        )
        if(!pr$family %in% allowed)
            stop("error: ", where, " uses ", pr$family, ", which is not supported for ", domain, " theta parameter ", theta_name)

        if(pr$family == "beta_dist"){
            pr$support <- c(lower = 0, upper = 1)
        } else if(is.null(pr$support)){
            stop("error: ", where, " must declare finite theta support; use uniform(lower, upper) or support = c(lower, upper)")
        }

        lower <- pr$support[1]
        upper <- pr$support[2]

        if(domain == "positive" && lower <= 0)
            stop("error: ", where, " support lower bound must be > 0 for positive theta parameter ", theta_name)
        if(domain == "nonnegative" && lower < 0)
            stop("error: ", where, " support lower bound must be >= 0 for nonnegative theta parameter ", theta_name)
        if(domain == "unit" && (lower < 0 || upper > 1))
            stop("error: ", where, " support must lie within [0,1] for unit theta parameter ", theta_name)
        if(domain == "bounded" && (lower <= -1 || upper >= 1))
            stop("error: ", where, " support must lie strictly inside (-1,1) for bounded theta parameter ", theta_name)
        list(prior = encode_prior(pr), bounds = c(lower = lower, upper = upper))
    }

    default_theta_names_types <- function(term){
        if(term$term_type %in% c("nngp", "gp")){
            registry <- build_cor_model_registry()
            info <- registry[[term$cov_model]]
            return(list(names = info$theta_names, types = as.integer(info$theta_types)))
        }
        if(term$term_type == "ar1")
            return(list(names = "phi", types = as.integer(3L)))
        if(term$term_type %in% c("car", "dagar"))
            return(list(names = "rho", types = as.integer(2L)))
        if(term$term_type %in% c("car_time", "dagar_time")){
            time_model <- tolower(term$params$time_model %||% "ar1")
            if(time_model == "ar1")
                return(list(names = c("rho", "phi"), types = as.integer(c(2L, 3L))))
            if(time_model == "exp")
                return(list(names = c("rho", "lambda"), types = as.integer(c(2L, 1L))))
            stop("error: unsupported ", term$term_type, " time_model '", time_model, "'")
        }
        stop("error: theta defaults not implemented for term_type '", term$term_type, "'")
    }

    validate_theta_vector <- function(val, theta_names, theta_types, term, where = NULL,
                                      fixed = NULL){
        context <- where %||% term$name %||% term$label

        val <- as.numeric(val)
        if(length(val) != length(theta_names))
            stop("error: theta vector for ", context, " must have length ", length(theta_names))
        if(is.null(fixed))
            fixed <- rep(FALSE, length(theta_names))
        if(length(fixed) != length(theta_names))
            stop("error: internal fixed theta flag length mismatch for ", context)

        for(ii in seq_along(val)){
            if(is.na(val[ii]))
                stop("error: theta contains NA for ", context)

            if(term$term_type == "ar1"){
                if(val[ii] <= -1 || val[ii] >= 1)
                    stop("error: ar1 phi must lie in (-1,1) for ", context)
            } else if(term$term_type %in% c("car", "dagar")){
                if(val[ii] <= 0 || val[ii] >= 1)
                    stop("error: ", term$term_type, " rho must lie in (0,1) for ", context)
            } else if(term$term_type %in% c("car_time", "dagar_time")){
                if(theta_names[ii] == "rho" && (val[ii] <= 0 || val[ii] >= 1))
                    stop("error: ", term$term_type, " rho must lie in (0,1) for ", context)
                if(theta_names[ii] == "phi" && (val[ii] <= -1 || val[ii] >= 1))
                    stop("error: ", term$term_type, " phi must lie in (-1,1) for ", context)
                if(theta_names[ii] == "lambda" && val[ii] <= 0)
                    stop("error: ", term$term_type, " lambda must be positive for ", context)
            } else if(term$term_type %in% c("nngp", "gp") &&
                      identical(term$cov_model, "gneiting") &&
                      theta_names[ii] == "delta"){
                if(isTRUE(fixed[ii])){
                    if(val[ii] < 0)
                        stop("error: fixed nonnegative theta required for ", context, ": ", theta_names[ii])
                } else if(val[ii] <= 0) {
                    stop("error: nonnegative theta must be positive when free for ", context, ": ", theta_names[ii])
                }
            } else if(theta_types[ii] == 1L && val[ii] <= 0){
                stop("error: positive theta required for ", context, ": ", theta_names[ii])
            } else if(theta_types[ii] == 2L){
                if(isTRUE(fixed[ii])){
                    if(val[ii] < 0 || val[ii] > 1)
                        stop("error: fixed unit-interval theta required in [0,1] for ", context, ": ", theta_names[ii])
                } else if(val[ii] <= 0 || val[ii] >= 1) {
                    stop("error: unit-interval theta required for ", context, ": ", theta_names[ii])
                }
            } else if(theta_types[ii] == 3L && (val[ii] <= -1 || val[ii] >= 1)){
                stop("error: bounded theta required in (-1,1) for ", context, ": ", theta_names[ii])
            }
        }

        stats::setNames(as.double(val), theta_names)
    }

    validate_theta_priors <- function(pr, theta_names, theta_types, term, where = NULL){
        context <- where %||% paste0("priors$", term$name %||% term$label)
        pr <- normalize_param_block(pr, context)

        check_missing_extra(names(pr), theta_names, context)

        out <- matrix(NA_real_, length(theta_names), 2,
                      dimnames = list(theta_names, c("lower", "upper")))
        prior_desc <- matrix(NA_real_, length(theta_names), 6,
                             dimnames = list(theta_names, c("family", "p1", "p2", "lower", "upper", "scale")))

        for(j in seq_along(theta_names)){
            nm <- theta_names[j]
            val <- pr[[nm]]

            if(is.null(val))
                stop("error: missing theta prior for ", nm, " in ", context)

            checked <- validate_theta_prior(
                val, nm, theta_types[j], term,
                where = paste0(context, "$", nm)
            )
            out[j, ] <- checked$bounds
            prior_desc[j, ] <- checked$prior
        }

        list(bounds = out, prior = prior_desc)
    }

    get_term_block <- function(cfg, term){
        if(length(cfg) == 0L)
            return(NULL)

        term_name <- tolower(term$name %||% term$label)
        term_label <- tolower(term$label %||% term_name)

        if(term_name %in% names(cfg))
            return(cfg[[term_name]])
        if(term_label %in% names(cfg))
            return(cfg[[term_label]])

        NULL
    }

    validate_term_block_names <- function(block, required_names, where){
        check_missing_extra(names(block), required_names, where)
        invisible(NULL)
    }

    ## Parse user controls for one process term using only the
    ## current term-keyed interface, e.g.
    ##   starting = list(nngp_1 = c(sigma_sq = 2, phi = 6))
    resolve_process_controls <- function(term, starting, tuning, priors){

        term_name <- term$name %||% term$label
        theta_info <- if(term$term_type %in% c("nngp", "gp", "ar1", "car", "car_time", "dagar", "dagar_time")) default_theta_names_types(term) else list(names = character(0), types = integer(0))
        theta_names <- theta_info$names
        theta_types <- theta_info$types
        required_names <- c("sigma_sq", theta_names)

        start_where <- paste0("starting$", term_name)
        tune_where  <- paste0("tuning$", term_name)
        prior_where <- paste0("priors$", term_name)

        start_block <- normalize_param_block(get_term_block(starting, term), start_where)
        tune_block  <- normalize_param_block(get_term_block(tuning, term), tune_where)
        prior_block <- normalize_param_block(get_term_block(priors, term), prior_where)

        if(!is.null(start_block)){
            extra_start_names <- setdiff(names(start_block), required_names)
            if(length(extra_start_names))
                stop(
                    "error: unknown parameter(s) ",
                    paste(extra_start_names, collapse = ", "),
                    " in ", start_where
                )
        }
        if(!is.null(tune_block))
            validate_term_block_names(tune_block, required_names, tune_where)
        if(!is.null(prior_block)){
            extra_prior_names <- setdiff(names(prior_block), required_names)
            if(length(extra_prior_names))
                stop(
                    "error: unknown parameter(s) ",
                    paste(extra_prior_names, collapse = ", "),
                    " in ", prior_where
                )
        }

        sigma_start_entry <- if(!is.null(start_block) && "sigma_sq" %in% names(start_block)) {
            require_positive_scalar(
                fixed_parameter_value(start_block[["sigma_sq"]], "sigma_sq", start_where),
                "sigma_sq", start_where
            )
        } else {
            1
        }

        sigma_fixed <- is_fixed_entry(start_block, "sigma_sq")
        if(sigma_fixed && !is.null(tune_block) && "sigma_sq" %in% names(tune_block) &&
           require_nonnegative_scalar(tune_block[["sigma_sq"]], "sigma_sq", tune_where) > 0)
            stop("error: fixed process parameter sigma_sq in ", start_where,
                 " cannot have positive tuning in ", tune_where)

        if(sigma_fixed){
            sigma_tune_entry <- 0
        } else if(!is.null(tune_block)){
            sigma_tune_entry <- require_nonnegative_scalar(tune_block[["sigma_sq"]], "sigma_sq", tune_where)
        } else {
            sigma_tune_entry <- 0.1
        }

        sigma_prior_entry <- if(!is.null(prior_block) && "sigma_sq" %in% names(prior_block)) {
            validate_variance_prior(prior_block[["sigma_sq"]], paste0(prior_where, "$sigma_sq"))
        } else if(sigma_tune_entry == 0) {
            no_op_variance_prior()
        } else {
            stop("error: missing prior for free process parameter sigma_sq in ", prior_where)
        }

        out <- list(
            sigma_sq_starting = as.double(sigma_start_entry),
            sigma_sq_tuning   = as.double(sigma_tune_entry),
            sigma_sq_IG       = legacy_ig_pair(sigma_prior_entry),
            sigma_sq_prior    = as.double(sigma_prior_entry)
        )

        if(length(theta_names) == 0L)
            return(out)

        theta_bounds <- matrix(NA_real_, length(theta_names), 2,
                               dimnames = list(theta_names, c("lower", "upper")))
        theta_prior <- matrix(NA_real_, length(theta_names), 6,
                              dimnames = list(theta_names, c("family", "p1", "p2", "lower", "upper", "scale")))
        theta_tune_entry <- stats::setNames(rep(NA_real_, length(theta_names)), theta_names)
        theta_start_entry <- stats::setNames(rep(NA_real_, length(theta_names)), theta_names)

        for(j in seq_along(theta_names)){
            nm <- theta_names[j]
            fixed_j <- is_fixed_entry(start_block, nm)

            if(fixed_j && !is.null(tune_block) && nm %in% names(tune_block) &&
               require_nonnegative_scalar(tune_block[[nm]], nm, tune_where) > 0)
                stop("error: fixed process parameter ", nm, " in ", start_where,
                     " cannot have positive tuning in ", tune_where)

            theta_tune_entry[nm] <- if(fixed_j) {
                0
            } else if(!is.null(tune_block)) {
                require_nonnegative_scalar(tune_block[[nm]], nm, tune_where)
            } else {
                0.1
            }

            has_prior <- !is.null(prior_block) && nm %in% names(prior_block)
            if(has_prior){
                checked <- validate_theta_prior(
                    prior_block[[nm]], nm, theta_types[j], term,
                    where = paste0(prior_where, "$", nm)
                )
                theta_bounds[nm, ] <- checked$bounds
                theta_prior[nm, ] <- checked$prior
            } else if(theta_tune_entry[nm] == 0) {
                theta_prior[nm, ] <- no_op_theta_prior()
                val <- if(!is.null(start_block) && nm %in% names(start_block))
                    fixed_parameter_value(start_block[[nm]], nm, start_where)
                else
                    NA_real_
                if(is.na(val))
                    stop("error: fixed process parameter ", nm,
                         " without a prior must have a starting value in ", start_where)
                pad <- max(abs(val) * 1e-8, 1e-8)
                theta_bounds[nm, ] <- c(val - pad, val + pad)
            } else {
                stop("error: missing prior for free process parameter ", nm, " in ", prior_where)
            }

            theta_start_entry[nm] <- if(!is.null(start_block) && nm %in% names(start_block)) {
                fixed_parameter_value(start_block[[nm]], nm, start_where)
            } else {
                mean(theta_bounds[nm, ])
            }
        }

        theta_start_entry <- validate_theta_vector(
            theta_start_entry,
            theta_names, theta_types, term, where = start_where,
            fixed = theta_tune_entry == 0
        )

        for(j in seq_along(theta_names)){
            nm <- theta_names[j]
            if(theta_tune_entry[nm] > 0 &&
               (theta_start_entry[nm] <= theta_bounds[nm, 1] ||
                theta_start_entry[nm] >= theta_bounds[nm, 2]))
                stop("error: theta starting value for ", term_name, "$", nm,
                     " must lie strictly inside the corresponding prior bounds")
        }

        out$theta_starting <- as.double(theta_start_entry)
        out$theta_tuning   <- as.double(theta_tune_entry)
        out$theta_names    <- theta_names
        out$theta_bounds   <- unname(theta_bounds)
        out$theta_prior    <- unname(theta_prior)
        out
    }

    ####################################################
    ## Check unused args
    ####################################################
    
    formal_args <- names(formals())
    elip_args <- list(...)
    
    for(i in names(elip_args))
        if(!i %in% formal_args)
            warning("'", i, "' is not an argument")
    
    ####################################################
    ## Formula checks
    ####################################################
    
    if(missing(formula))
        stop("error: formula must be specified")
    
    if(!inherits(formula, "formula"))
        stop("error: formula is misspecified")

    likelihood_family <- normalize_family(family)
    is_binomial <- identical(likelihood_family, "binomial")
    is_negbin <- identical(likelihood_family, "negative_binomial")
    is_pg_likelihood <- is_binomial || is_negbin

    resid_info <- build_resid_components(formula = formula, data = data)
    residual <- resid_info$residual
    formula_no_resid <- resid_info$reduced_formula
    
    ####################################################
    ## Structured process terms
    ####################################################

    formula_text <- paste(deparse(formula_no_resid), collapse = " ")
    has_nngp_term <- grepl("\\bnngp\\s*\\(", formula_text)
    has_gp_term <- grepl("\\bgp\\s*\\(", formula_text)
    has_ar1_term <- grepl("\\bar1\\s*\\(", formula_text)
    has_car_term <- grepl("\\bcar\\s*\\(", formula_text)
    has_car_time_term <- grepl("\\bcar_time\\s*\\(", formula_text)
    has_dagar_term <- grepl("\\bdagar\\s*\\(", formula_text)
    has_dagar_time_term <- grepl("\\bdagar_time\\s*\\(", formula_text)

    if(has_nngp_term)
        status_message("Building NNGP graph(s).")
    if(has_gp_term || has_ar1_term || has_car_term || has_car_time_term || has_dagar_term || has_dagar_time_term)
        status_message("Building structured process support.")
    
    proc_info <- build_process_components(
        formula = formula_no_resid,
        data = data,
        n_omp_threads = n_omp_threads,
        nngp_search = nngp_search
    )
    
    if(length(proc_info$process_terms))
        status_message("Built ", length(proc_info$process_terms), " structured process term(s).")
    
    iid_info <- build_iid_components(
        formula = proc_info$reduced_formula,
        data = data
    )

    if(length(iid_info$terms))
        status_message("Built ", length(iid_info$terms), " iid random-effect term(s).")

    work_formula <- iid_info$reduced_formula

    tt <- terms(work_formula, keep.order = TRUE)
    term_labels <- attr(tt, "term.labels")

    if(any(grepl("\\|", term_labels)))
        stop(
            "error: use iid(group) for random intercepts and x:iid(group) for random slopes"
        )
    
    ####################################################
    ## Fixed-effect model frame and response partition
    ####################################################

    ## Predictors must be complete. Missing responses are allowed and are
    ## carried in the fitted object, but only observed rows enter the collapsed
    ## likelihood passed to C++.
    mf <- model.frame(work_formula, data, na.action = stats::na.pass)

    response_col <- attr(terms(work_formula), "response")
    predictor_cols <- if(ncol(mf)) setdiff(seq_len(ncol(mf)), response_col) else integer(0)
    if(length(predictor_cols) && anyNA(mf[, predictor_cols, drop = FALSE]))
        stop("error: missing data in model predictors")

    y <- model.response(mf)
    if(is_pg_likelihood && is.matrix(y))
        stop("error: matrix responses are not yet supported for Polya-Gamma likelihoods; use a numeric response vector")
    if(anyNA(y) && !is.numeric(y))
        stop("error: missing response values are only supported for numeric responses")

    X <- model.matrix(work_formula, mf)
    x_names <- colnames(X)
    offset <- model.offset(mf)
    if(is.null(offset))
        offset <- rep(0, nrow(X))
    offset <- as.numeric(offset)
    if(length(offset) != nrow(X) || anyNA(offset) || any(!is.finite(offset)))
        stop("error: offset() values must be finite numeric values with length matching the data")
    has_offset <- any(offset != 0)
    
    ####################################################
    ## Observed likelihood rows and iid random-effect design
    ####################################################
    
    n <- nrow(X)
    observed_index <- which(!is.na(y))
    missing_index <- which(is.na(y))
    n_full <- n
    n_obs <- length(observed_index)
    n_missing_response <- length(missing_index)
    if(n_obs < 1L)
        stop("error: at least one observed response is required")

    trials_vec <- if(is_binomial) normalize_trials(trials, data, n_full) else NULL
    nb_size <- if(is_negbin) normalize_nb_size(size) else NULL
    if(!is_binomial && !is.null(trials))
        stop("error: trials is used only with family = \"binomial\"")
    if(!is_negbin && !is.null(size))
        stop("error: size is used only with family = \"negative_binomial\"")
    if(is_binomial){
        y_num <- as.numeric(y)
        y_obs_check <- y_num[observed_index]
        trials_obs_check <- trials_vec[observed_index]
        if(any(!is.finite(y_obs_check)) || any(y_obs_check < 0))
            stop("error: binomial responses must be finite non-negative counts")
        if(any(abs(y_obs_check - round(y_obs_check)) > sqrt(.Machine$double.eps)))
            stop("error: binomial responses must be integer-valued counts")
        if(any(y_obs_check > trials_obs_check))
            stop("error: binomial responses must be less than or equal to trials")
    }
    if(is_negbin){
        y_num <- as.numeric(y)
        y_obs_check <- y_num[observed_index]
        if(any(!is.finite(y_obs_check)) || any(y_obs_check < 0))
            stop("error: negative-binomial responses must be finite non-negative counts")
        if(any(abs(y_obs_check - round(y_obs_check)) > sqrt(.Machine$double.eps)))
            stop("error: negative-binomial responses must be integer-valued counts")
    }

    ####################################################
    ## Normalize user controls
    ####################################################

    starting <- normalize_named_list(starting, "starting")
    tuning   <- normalize_named_list(tuning,   "tuning")
    priors   <- normalize_named_list(priors,   "priors")

    residual_model <- list(type = "global_tau", label = "tau_sq")
    if(is_pg_likelihood && !is.null(residual))
        stop("error: resid() is not used with family = \"", likelihood_family, "\"")
    if(!is.null(residual) && !identical(residual$type, "global_tau")){
        if("tau_sq" %in% names(starting))
            stop("error: starting$tau_sq is not used with this residual variance model; use starting$resid for residual variance controls")
        if("tau_sq" %in% names(tuning))
            stop("error: tuning$tau_sq is not used with this residual variance model; use tuning$resid for residual variance controls")
        if("tau_sq" %in% names(priors))
            stop("error: priors$tau_sq is not used with this residual variance model; use priors$resid for residual variance controls")
    }
    if(!is.null(residual)){
        if(!inherits(residual, "stLMM_residual") || is.null(residual$type))
            stop("error: residual model must be created by resid()")
        if(residual$type == "global_tau"){
            residual_model <- list(type = "global_tau", label = residual$label %||% "tau_sq")
        } else if(residual$type == "fixed_variance"){
            residual_variance <- eval_residual_variance(
                residual = residual,
                data = data,
                n_expected = n,
                where = "fitted data"
            )
            validate_fixed_residual_variance(
                variance = residual_variance,
                observed_index = observed_index,
                where = "fitted data"
            )
            residual_model <- list(
                type = "fixed_variance",
                variance_label = residual$variance_label,
                variance_expr = residual$variance_expr,
                variance = as.double(residual_variance),
                variance_obs = as.double(residual_variance[observed_index]),
                obs_precision_obs = as.double(1 / residual_variance[observed_index])
            )
        } else if(residual$type == "group_ig_variance"){
            residual_model <- build_group_ig_residual_model(
                residual = residual,
                data = data,
                observed_index = observed_index,
                n_expected = n,
                where = "fitted data"
            )
        } else if(residual$type == "group_variance"){
            residual_model <- build_group_residual_model(
                residual = residual,
                data = data,
                observed_index = observed_index,
                n_expected = n,
                where = "fitted data",
                starting = starting,
                tuning = tuning,
                priors = priors,
                validate_variance_prior = validate_variance_prior,
                prior_family_code = prior_family_code
            )
        } else if(residual$type == "scaled_variance"){
            residual_model <- build_scaled_residual_model(
                residual = residual,
                data = data,
                observed_index = observed_index,
                n_expected = n,
                where = "fitted data"
            )
        } else {
            stop("error: unsupported residual model type: ", residual$type)
        }
    }
    fixed_residual_variance <- identical(residual_model$type, "fixed_variance")
    sampled_residual_variance <- residual_model$type %in% c("group_ig_variance", "scaled_variance")

    p <- ncol(X)
    iid_design <- build_iid_design(iid_info, data, n)
    Z <- iid_design$Z
    q <- ncol(Z)
    
    ####################################################
    ## Explicit iid random-effect structure
    ####################################################
    
    re_names   <- iid_design$re_names
    re_nlevels <- iid_design$re_nlevels
    re_ncoef   <- iid_design$re_ncoef
    re.q       <- iid_design$re.q
    re_block_id <- iid_design$re_block_id
    
    n_re <- length(re_names)
    
    if(length(re_block_id) != q)
        stop("error: internal re_block_id length mismatch")

    ####################################################
    ## Structured process term naming and validation
    ####################################################
    
    graphs <- proc_info$graphs
    terms  <- proc_info$process_terms

    if(length(terms)){
        type_counter <- integer(0)
        for(i in seq_along(terms)){
            type_i <- terms[[i]]$term_type
            type_counter[type_i] <- if(type_i %in% names(type_counter)) type_counter[type_i] + 1L else 1L
            if(is.null(terms[[i]]$name) || !nzchar(terms[[i]]$name))
                terms[[i]]$name <- paste0(type_i, "_", type_counter[type_i])
            terms[[i]]$name <- tolower(terms[[i]]$name)
        }
    }

    for(i in seq_along(terms)){
        if(!is.null(terms[[i]]$x) && anyNA(terms[[i]]$x))
            stop("error: SVC covariate ", terms[[i]]$coef_name, " contains missing values")
    }

    beta_prior <- validate_beta_prior(
        priors$beta,
        p = p,
        where = "priors$beta",
        default_family = if(is_pg_likelihood) "normal" else "flat"
    )

    if(is_pg_likelihood){
        if("tau_sq" %in% names(starting))
            stop("error: starting$tau_sq is not used with family = \"", likelihood_family, "\"")
        if("resid" %in% names(starting))
            stop("error: starting$resid is not used with family = \"", likelihood_family, "\"")
        if("tau_sq" %in% names(tuning))
            stop("error: tuning$tau_sq is not used with family = \"", likelihood_family, "\"")
        if("resid" %in% names(tuning))
            stop("error: tuning$resid is not used with family = \"", likelihood_family, "\"")
        if("tau_sq" %in% names(priors))
            stop("error: priors$tau_sq is not used with family = \"", likelihood_family, "\"")
        if("resid" %in% names(priors))
            stop("error: priors$resid is not used with family = \"", likelihood_family, "\"")
    }

    tau_sq_start_control <- get_resid_control(starting, "tau_sq")
    tau_sq_tune_control <- get_resid_control(tuning, "tau_sq")
    tau_sq_prior_control <- get_resid_control(priors, "tau_sq")
    if(is.null(tau_sq_start_control) && "tau_sq" %in% names(starting))
        tau_sq_start_control <- starting$tau_sq
    if(is.null(tau_sq_tune_control) && "tau_sq" %in% names(tuning))
        tau_sq_tune_control <- tuning$tau_sq
    if(is.null(tau_sq_prior_control) && "tau_sq" %in% names(priors))
        tau_sq_prior_control <- priors$tau_sq

    if((fixed_residual_variance || sampled_residual_variance) && "tau_sq" %in% names(starting))
        stop("error: starting$tau_sq is not used with this residual variance model")
    if((fixed_residual_variance || sampled_residual_variance) && "tau_sq" %in% names(tuning))
        stop("error: tuning$tau_sq is not used with this residual variance model")
    if((fixed_residual_variance || sampled_residual_variance) && "tau_sq" %in% names(priors))
        stop("error: priors$tau_sq is not used with this residual variance model")

    if(fixed_residual_variance || sampled_residual_variance){
        tau_sq_start_control <- NULL
        tau_sq_tune_control <- NULL
        tau_sq_prior_control <- NULL
    }

    tau_sq_fixed <- !is.null(tau_sq_start_control) && is_fixed_parameter(tau_sq_start_control)

    tau_sq_starting <- if(fixed_residual_variance || sampled_residual_variance) {
        1
    } else if(is_pg_likelihood) {
        1
    } else if(!is.null(tau_sq_start_control)) {
        require_positive_scalar(
            fixed_parameter_value(tau_sq_start_control, "tau_sq", "starting$resid"),
            "tau_sq_starting", "starting"
        )
    } else {
        default_positive_scale(as.numeric(y[observed_index]) - offset[observed_index])
    }

    if(tau_sq_fixed && !is.null(tau_sq_tune_control) &&
       require_nonnegative_scalar(tau_sq_tune_control, "tau_sq", "tuning$resid") > 0)
        stop("error: fixed tau_sq in starting cannot have positive tuning")

    tau_sq_tuning <- if(fixed_residual_variance || sampled_residual_variance) {
        0
    } else if(is_pg_likelihood) {
        0
    } else if(tau_sq_fixed) {
        0
    } else if(!is.null(tau_sq_tune_control)) {
        require_nonnegative_scalar(tau_sq_tune_control, "tau_sq", "tuning$resid")
    } else 0.2

    tau_sq_IG <- if(fixed_residual_variance || sampled_residual_variance) {
        c(1, 1)
    } else if(is_pg_likelihood) {
        c(1, 1)
    } else if(!is.null(tau_sq_prior_control)) {
        tau_prior_desc <- validate_variance_prior(tau_sq_prior_control, "priors$resid$tau_sq")
        legacy_ig_pair(tau_prior_desc)
    } else c(2, 1)

    tau_sq_prior <- if(fixed_residual_variance || sampled_residual_variance) {
        encode_prior(ig(1, 1))
    } else if(is_pg_likelihood) {
        encode_prior(ig(1, 1))
    } else if(!is.null(tau_sq_prior_control)) {
        validate_variance_prior(tau_sq_prior_control, "priors$resid$tau_sq")
    } else {
        encode_prior(ig(2, 1))
    }

    if(n_re){
        if("sigma_sq_re" %in% names(starting))
            stop("error: starting$sigma_sq_re is no longer supported; use starting = list(iid_1 = list(sigma_sq = value))")
        if("sigma_sq_re" %in% names(priors))
            stop("error: priors$sigma_sq_re is no longer supported; use priors = list(iid_1 = list(sigma_sq = ig(shape, scale)))")

        sigma_sq_re_start_list <- stats::setNames(lapply(re_names, function(nm) starting[[nm]]), re_names)
        start_where_re <- "starting"

        missing_re_start <- re_names[vapply(sigma_sq_re_start_list[re_names], is.null, logical(1))]
        if(length(missing_re_start)){
            for(nm in missing_re_start)
                sigma_sq_re_start_list[[nm]] <- list(sigma_sq = 1)
        }

        sigma_sq_re_starting <- collect_re_scalars(sigma_sq_re_start_list, re_names, start_where_re)

        sigma_sq_re_prior_list <- stats::setNames(vector("list", n_re), re_names)
        for(i in seq_along(re_names))
            sigma_sq_re_prior_list[[re_names[i]]] <- priors[[re_names[i]]]
        prior_where_re <- "priors"

        missing_re_prior <- re_names[vapply(sigma_sq_re_prior_list[re_names], is.null, logical(1))]
        if(length(missing_re_prior))
            stop("error: missing prior(s) for grouped random-effect variance ",
                 paste(missing_re_prior, collapse = ", "), " in ", prior_where_re)

        sigma_sq_re_IG <- collect_re_pairs(
            sigma_sq_re_prior_list,
            re_names,
            prior_where_re
        )
    } else {
        sigma_sq_re_starting <- numeric(0)
        sigma_sq_re_IG <- matrix(0, 0, 2)
    }
    
    ## assign graph indices and term-level parameter metadata
    graph_names <- names(graphs)
    
    for(i in seq_along(terms)){
        gid <- terms[[i]]$graph_id
        gi  <- match(gid, graph_names)
        
        if(is.na(gi))
            stop("error: graph_id '", gid, "' not found")
        
        terms[[i]]$graph_index <- as.integer(gi)

        ctrl <- resolve_process_controls(
            term = terms[[i]],
            starting = starting,
            tuning = tuning,
            priors = priors
        )

        terms[[i]]$sigma_sq_starting <- ctrl$sigma_sq_starting
        terms[[i]]$sigma_sq_tuning   <- ctrl$sigma_sq_tuning
        terms[[i]]$sigma_sq_IG       <- ctrl$sigma_sq_IG
        terms[[i]]$sigma_sq_prior    <- ctrl$sigma_sq_prior

        if(!is.null(ctrl$theta_starting)){
            terms[[i]]$theta_starting <- ctrl$theta_starting
            terms[[i]]$theta_tuning   <- ctrl$theta_tuning
            terms[[i]]$theta_names    <- ctrl$theta_names
            terms[[i]]$theta_bounds   <- ctrl$theta_bounds
            terms[[i]]$theta_prior    <- ctrl$theta_prior
        }
        
        ## enforce integer storage where needed
        ## NOTE:
        ##   - graph$nnIndx and graph$nnIndxLU are already 0-based for direct C use
        ##   - term-level observation/node indexing is intentionally kept 1-based
        ##     because the same objects are also consumed on the R side
        terms[[i]]$map <- as.integer(terms[[i]]$map)
        terms[[i]]$obsIndx <- as.integer(terms[[i]]$obsIndx)
        terms[[i]]$obsIndxLU <- matrix(as.integer(terms[[i]]$obsIndxLU), ncol = 2L)
        terms[[i]]$node_nobs <- as.integer(terms[[i]]$node_nobs)
        terms[[i]]$n_obs <- as.integer(terms[[i]]$n_obs)
        terms[[i]]$n_node <- as.integer(terms[[i]]$n_node)
    }

    storage.mode(y) <- "double"
    storage.mode(X) <- "double"
    #Z = Z dgCMatrix
    n <- as.integer(n)
    n_full <- as.integer(n_full)
    n_obs <- as.integer(n_obs)
    observed_index <- as.integer(observed_index)
    missing_index <- as.integer(missing_index)
    p <- as.integer(p)
    q <- as.integer(q)

    ## Keep full-data objects for fitted values, recovery, summaries, and
    ## prediction. The sampler sees the observed-row contract below.
    y_obs <- y[observed_index]
    if(is_pg_likelihood)
        y_obs <- as.double(round(as.numeric(y_obs)))
    offset_obs <- offset[observed_index]
    y_model_obs <- y_obs
    if(!is_pg_likelihood)
        y_model_obs <- as.double(as.numeric(y_obs) - offset_obs)
    X_obs <- X[observed_index, , drop = FALSE]
    storage.mode(y_obs) <- "double"
    storage.mode(y_model_obs) <- "double"
    storage.mode(offset) <- "double"
    storage.mode(offset_obs) <- "double"
    storage.mode(X_obs) <- "double"
    Z_obs <- Z[observed_index, , drop = FALSE]
    terms_obs <- lapply(terms, subset_process_term_observed, observed_index = observed_index)
    for(i in seq_along(terms_obs)){
        terms_obs[[i]]$map <- as.integer(terms_obs[[i]]$map)
        terms_obs[[i]]$obsIndx <- as.integer(terms_obs[[i]]$obsIndx)
        terms_obs[[i]]$obsIndxLU <- matrix(as.integer(terms_obs[[i]]$obsIndxLU), ncol = 2L)
        terms_obs[[i]]$node_nobs <- as.integer(terms_obs[[i]]$node_nobs)
        terms_obs[[i]]$n_obs <- as.integer(terms_obs[[i]]$n_obs)
        terms_obs[[i]]$n_node <- as.integer(terms_obs[[i]]$n_node)
        if(!is.null(terms_obs[[i]]$x))
            storage.mode(terms_obs[[i]]$x) <- "double"
    }

    ## metadata used by the C side to print the pre-MCMC term report and
    ## returned afterward as fit$term_description.
    formula_txt <- paste(deparse(formula, width.cutoff = 500L), collapse = " ")

    re_term_description <- vector("list", n_re)
    if(n_re){
        for(i in seq_len(n_re)){
            iid_term <- iid_design$terms[[i]]
            coef_names <- iid_term$coefficient
            coef_tag <- if(length(coef_names)) paste(coef_names, collapse = ", ") else "(none)"
            re_term_description[[i]] <- list(
                name = re_names[i],
                label = iid_term$label,
                grouping_factor = iid_term$group_name,
                coefficients = coef_names,
                levels = iid_term$levels,
                n_levels = as.integer(re_nlevels[i]),
                n_coefficients = as.integer(re_ncoef[i]),
                q_contribution = as.integer(re.q[i]),
                sigma_sq = list(
                    current = as.double(sigma_sq_re_starting[i]),
                    tuning = NA_real_,
                    prior_family = "inverse-gamma",
                    prior_hyperparameters = c(shape = sigma_sq_re_IG[i, 1], scale = sigma_sq_re_IG[i, 2])
                ),
                sampler = list(
                    block_acceptance = NA_real_,
                    status = "gibbs"
                )
            )
        }
        names(re_term_description) <- re_names
    }

    process_term_description <- vector("list", length(terms))
    if(length(terms)){
        for(i in seq_along(terms)){
            type_i <- terms[[i]]$term_type

            sigma_start <- if(!is.null(terms[[i]]$sigma_sq_starting)) as.double(terms[[i]]$sigma_sq_starting) else 1.0
            sigma_tune  <- if(!is.null(terms[[i]]$sigma_sq_tuning)) as.double(terms[[i]]$sigma_sq_tuning) else 0.1
            sigma_ig    <- if(!is.null(terms[[i]]$sigma_sq_IG)) as.double(terms[[i]]$sigma_sq_IG) else c(2, 1)

            theta_info <- vector("list", length(terms[[i]]$theta_names))
            if(length(theta_info)){
                for(j in seq_along(theta_info)){
                    bounds_j <- terms[[i]]$theta_bounds[j, ]
                    theta_info[[j]] <- list(
                        name = terms[[i]]$theta_names[j],
                        current = as.double(terms[[i]]$theta_starting[j]),
                        tuning = as.double(terms[[i]]$theta_tuning[j]),
                        lower = as.double(bounds_j[1]),
                        upper = as.double(bounds_j[2]),
                        prior_type = "uniform",
                        transform = "bounded-logit"
                    )
                }
                names(theta_info) <- terms[[i]]$theta_names
            }

            graph_i <- graphs[[terms[[i]]$graph_index]]
            repeated_flag <- NA
            unique_count <- as.integer(NA)
            coord_cols <- NULL
            coord_dim <- as.integer(NA)
            time_col <- NULL
            neighbor_m <- as.integer(NA)
            ordering_method <- NULL
            repeated_collapsed <- NA

            if(type_i == "nngp"){
                repeated_flag <- if(!is.null(terms[[i]]$n_obs) && !is.null(terms[[i]]$n_node)) (terms[[i]]$n_obs > terms[[i]]$n_node) else NA
                unique_count <- if(!is.null(terms[[i]]$n_node)) as.integer(terms[[i]]$n_node) else as.integer(NA)
                coord_cols <- graph_i$coord_cols %||% terms[[i]]$coord_cols
                coord_dim <- if(!is.null(graph_i$dim)) as.integer(graph_i$dim) else as.integer(NA)
                neighbor_m <- if(!is.null(graph_i$m)) as.integer(graph_i$m) else as.integer(NA)
                ordering_method <- graph_i$order %||% graph_i$ordering %||% terms[[i]]$order %||% terms[[i]]$ordering
                repeated_collapsed <- repeated_flag
            } else if(type_i == "gp"){
                repeated_flag <- if(!is.null(terms[[i]]$n_obs) && !is.null(terms[[i]]$n_node)) (terms[[i]]$n_obs > terms[[i]]$n_node) else NA
                unique_count <- if(!is.null(terms[[i]]$n_node)) as.integer(terms[[i]]$n_node) else as.integer(NA)
                coord_cols <- graph_i$coord_cols %||% terms[[i]]$coord_cols
                coord_dim <- if(!is.null(graph_i$dim)) as.integer(graph_i$dim) else as.integer(NA)
                repeated_collapsed <- repeated_flag
            } else if(type_i == "ar1"){
                repeated_flag <- if(!is.null(terms[[i]]$n_obs) && !is.null(terms[[i]]$n_node)) (terms[[i]]$n_obs > terms[[i]]$n_node) else NA
                unique_count <- if(!is.null(terms[[i]]$n_node)) as.integer(terms[[i]]$n_node) else as.integer(NA)
                time_col <- graph_i$index_col %||% graph_i$time_col %||% terms[[i]]$index_col %||% terms[[i]]$time_col
                repeated_collapsed <- repeated_flag
            } else if(type_i == "car_time"){
                repeated_flag <- if(!is.null(terms[[i]]$n_obs) && !is.null(terms[[i]]$n_node)) (terms[[i]]$n_obs > terms[[i]]$n_node) else NA
                unique_count <- if(!is.null(terms[[i]]$n_node)) as.integer(terms[[i]]$n_node) else as.integer(NA)
                time_col <- graph_i$time_col %||% terms[[i]]$time_col
                repeated_collapsed <- repeated_flag
            } else if(type_i == "dagar"){
                repeated_flag <- if(!is.null(terms[[i]]$n_obs) && !is.null(terms[[i]]$n_node)) (terms[[i]]$n_obs > terms[[i]]$n_node) else NA
                unique_count <- if(!is.null(terms[[i]]$n_node)) as.integer(terms[[i]]$n_node) else as.integer(NA)
                ordering_method <- graph_i$ordering %||% terms[[i]]$ordering
                repeated_collapsed <- repeated_flag
            } else if(type_i == "dagar_time"){
                repeated_flag <- if(!is.null(terms[[i]]$n_obs) && !is.null(terms[[i]]$n_node)) (terms[[i]]$n_obs > terms[[i]]$n_node) else NA
                unique_count <- if(!is.null(terms[[i]]$n_node)) as.integer(terms[[i]]$n_node) else as.integer(NA)
                time_col <- graph_i$time_col %||% terms[[i]]$time_col
                ordering_method <- graph_i$ordering %||% terms[[i]]$ordering
                repeated_collapsed <- repeated_flag
            }

            process_term_description[[i]] <- list(
                name = terms[[i]]$name,
                label = terms[[i]]$label,
                type = type_i,
                graph_index = as.integer(terms[[i]]$graph_index),
                cov_model = if(!is.null(terms[[i]]$cov_model)) terms[[i]]$cov_model else NA_character_,
                is_svc = !is.null(terms[[i]]$x),
                svc_covariate = if(!is.null(terms[[i]]$covariate)) terms[[i]]$covariate else NA_character_,
                n_obs = as.integer(terms[[i]]$n_obs),
                n_node = as.integer(terms[[i]]$n_node),
                q_lat = as.integer(NA),
                w_offset = as.integer(NA),
                repeated_index_collapsed = repeated_flag,
                unique_count = unique_count,
                sigma_sq = list(
                    current = sigma_start,
                    tuning = sigma_tune,
                    prior_family = "inverse-gamma",
                    prior_hyperparameters = c(shape = sigma_ig[1], scale = sigma_ig[2])
                ),
                theta = theta_info,
                diagnostics = list(
                    coord_columns = coord_cols,
                    coord_dim = coord_dim,
                    neighbor_count_m = neighbor_m,
                    ordering = ordering_method,
                    time_column = time_col,
                    repeated_collapsed = repeated_collapsed,
                    zero_parent_nodes = as.integer(NA),
                    prior_precision_nnz = as.integer(NA),
                    neighbor_summary = c(min = NA_real_, mean = NA_real_, max = NA_real_)
                ),
                sampler = list(
                    block_acceptance = NA_real_,
                    status = "updating"
                )
            )
        }
        names(process_term_description) <- vapply(process_term_description, `[[`, character(1), "name")
    }
    
    ####################################################
    ## Backend contract for C++
    ####################################################

    status_message("Preparing sampler backend.")

    n_report <- normalize_n_report(n_report)
    metropolis <- normalize_metropolis(metropolis)
    warmup <- normalize_warmup(warmup, metropolis)
    save_process_control <- normalize_save_process(
        save_process,
        is_pg_likelihood = is_pg_likelihood,
        has_process = length(terms) > 0L
    )

    backend <- list(
        
        y = y,
        y_obs = y_obs,
        y_model_obs = y_model_obs,
        offset = offset,
        offset_obs = offset_obs,
        has_offset = isTRUE(has_offset),
        X = X,
        X_obs = X_obs,
        Z = Z,
        Z_obs = Z_obs,
        n = n,
        n_full = n_full,
        n_obs = n_obs,
        observed_index = observed_index,
        missing_index = missing_index,
        n_missing_response = as.integer(n_missing_response),
        p = p,
        q = q,
        n_samples = as.integer(n_samples),
        n_report = as.integer(n_report),
        warmup = list(
            enabled = isTRUE(warmup$enabled),
            batch_length = as.integer(warmup$batch_length),
            min_batches = as.integer(warmup$min_batches),
            max_batches = as.integer(warmup$max_batches),
            target = as.double(warmup$target),
            near_zero = as.double(warmup$near_zero)
        ),
        metropolis_blocking = as.integer(metropolis$blocking_code),
        metropolis_batch_length = as.integer(metropolis$batch_length),
        metropolis_target_accept = as.double(metropolis$target_accept),
        metropolis = metropolis,
        cholmod_control = cholmod_control,
        n_omp_threads = as.integer(n_omp_threads),
        nngp_search = nngp_search,
        verbose = verbose,
        describe_terms = isTRUE(describe_terms),
        family = likelihood_family,
        likelihood = list(
            family = likelihood_family,
            trials = if(is_binomial) as.integer(trials_vec) else NULL,
            trials_obs = if(is_binomial) as.integer(trials_vec[observed_index]) else NULL,
            size = if(is_negbin) as.double(nb_size) else NULL
        ),
        recover_process = save_process_control$enabled,
        recover_start = save_process_control$start,
        recover_thin = save_process_control$thin,
        save_process = save_process_control,
        formula = formula,
        reduced_formula = work_formula,
        response_name = all.vars(formula)[1],
        trials = if(is_binomial) as.integer(trials_vec) else NULL,
        trials_obs = if(is_binomial) as.integer(trials_vec[observed_index]) else NULL,
        nb_size = if(is_negbin) as.double(nb_size) else NULL,
        beta_starting = if(is_pg_likelihood) {
            beta_start <- if("beta" %in% names(starting)) {
                val <- as.numeric(starting$beta)
                if(length(val) != p || any(!is.finite(val)))
                    stop("error: starting$beta must be a finite numeric vector with length matching the fixed effects")
                val
            } else {
                rep(0, p)
            }
            as.double(beta_start)
        } else NULL,
        beta_prior_type = as.integer(beta_prior$type),
        beta_prior_mean = as.double(beta_prior$mean),
        beta_prior_precision = as.double(beta_prior$precision),

        tau_sq_starting = as.double(tau_sq_starting),
        tau_sq_tuning = as.double(tau_sq_tuning),
        tau_sq_IG = as.double(tau_sq_IG),
        tau_sq_prior = as.double(tau_sq_prior),
        residual_model = residual_model,
        
        ## random effects
        sigma_sq_re_starting = as.double(sigma_sq_re_starting),
        sigma_sq_re_IG = sigma_sq_re_IG,
        re = list(
            q = as.integer(re.q),
            block_id = as.integer(re_block_id)
        ),
        
        ## process info
        graphs = graphs,
        process_terms = terms,
        process_terms_obs = terms_obs,
        term_description_meta = list(
            global = list(
                formula = formula_txt,
                family = likelihood_family,
                n = as.integer(n),
                n_obs = as.integer(n_obs),
                n_missing_response = as.integer(n_missing_response),
                residual = list(
                    type = residual_model$type,
                    label = residual_model$variance_label %||% residual_model$label %||% NA_character_
                ),
                p = as.integer(p),
                q = as.integer(q),
                qLatTotal = as.integer(NA),
                n_process_terms = as.integer(length(terms)),
                n_graphs = as.integer(length(graphs)),
                M_dim = as.integer(c(NA, NA)),
                M_nnz = as.integer(NA),
                cholmod_requested_ordering = cholmod_control$ordering,
                cholmod_postorder = isTRUE(cholmod_control$postorder),
                cholmod_ordering = NA_character_,
                cholmod_fill_ratio = NA_real_,
                cholmod_lnz = NA_real_,
                cholmod_flops = NA_real_,
                factorization_status = NA_character_
            ),
            fixed_effects = list(
                p = as.integer(p),
                names = x_names,
                has_intercept = "(Intercept)" %in% x_names
            ),
            random_effects = list(
                q = as.integer(q),
                n_terms = as.integer(n_re),
                terms = re_term_description
            ),
            process_terms = process_term_description
        )
        
    )
    
    class(backend) <- "stLMM_backend"

    process_names <- vapply(backend$process_terms, `[[`, character(1), "name")

    chain_value <- function(x, chain, chains, name, positive = FALSE, lower = -Inf, upper = Inf){
        if(is_fixed_parameter(x))
            return(fixed_parameter_value(x, name, "starting"))
        if(length(x) == 0L)
            stop("error: missing starting value for ", name)
        if(!(length(x) %in% c(1L, chains)))
            stop("error: starting value for ", name, " must have length 1 or chains")
        val <- as.numeric(x[if(length(x) == 1L) 1L else chain])
        if(length(val) != 1L || is.na(val) || !is.finite(val))
            stop("error: starting value for ", name, " must be finite")
        if(positive && val <= 0)
            stop("error: starting value for ", name, " must be positive")
        if(val <= lower || val >= upper)
            stop("error: starting value for ", name, " must lie in (", lower, ", ", upper, ")")
        val
    }

    log_jitter <- function(center, dispersion){
        center * exp(stats::rnorm(1L, 0, dispersion))
    }

    bounded_jitter <- function(lower, upper){
        eps <- 0.05
        stats::runif(1L, lower + eps * (upper - lower), upper - eps * (upper - lower))
    }

    process_start_block <- function(starting, term){
        block <- get_term_block(starting, term)
        if(is.null(block))
            return(NULL)
        normalize_param_block(block, paste0("starting$", term$name))
    }

    build_chain_backend <- function(chain){
        b <- backend
        if(!(fixed_residual_variance || sampled_residual_variance)){
            if(!is.null(tau_sq_start_control))
                b$tau_sq_starting <- chain_value(tau_sq_start_control, chain, chains, "tau_sq", positive = TRUE)
            else if(chains > 1L)
                b$tau_sq_starting <- log_jitter(tau_sq_starting, chain_control$dispersion)
        }

        if(n_re){
            for(i in seq_along(re_names)){
                start_entry <- NULL
                if(re_names[i] %in% names(starting))
                    start_entry <- normalize_iid_sigma_block(starting[[re_names[i]]], paste0("starting$", re_names[i]))
                if(!is.null(start_entry))
                    b$sigma_sq_re_starting[i] <- chain_value(start_entry, chain, chains, re_names[i], positive = TRUE)
                else if(chains > 1L)
                    b$sigma_sq_re_starting[i] <- log_jitter(sigma_sq_re_starting[i], chain_control$dispersion)
            }
        }

        for(i in seq_along(b$process_terms)){
            term <- b$process_terms[[i]]
            block <- process_start_block(starting, term)
            if(!is.null(block) && "sigma_sq" %in% names(block))
                b$process_terms[[i]]$sigma_sq_starting <- chain_value(block[["sigma_sq"]], chain, chains, paste0(term$name, "$sigma_sq"), positive = TRUE)
            else if(chains > 1L && term$sigma_sq_tuning > 0)
                b$process_terms[[i]]$sigma_sq_starting <- log_jitter(term$sigma_sq_starting, chain_control$dispersion)

            theta_names_i <- term$theta_names %||% character(0)
            if(length(theta_names_i)){
                for(j in seq_along(theta_names_i)){
                    nm <- theta_names_i[j]
                    bounds <- term$theta_bounds[j, ]
                    if(!is.null(block) && nm %in% names(block)){
                        b$process_terms[[i]]$theta_starting[j] <- chain_value(
                            block[[nm]], chain, chains, paste0(term$name, "$", nm),
                            lower = bounds[1], upper = bounds[2]
                        )
                    } else if(chains > 1L && term$theta_tuning[j] > 0) {
                        b$process_terms[[i]]$theta_starting[j] <- bounded_jitter(bounds[1], bounds[2])
                    }
                }
            }
        }
        b$process_terms_obs <- lapply(b$process_terms, subset_process_term_observed, observed_index = observed_index)
        for(i in seq_along(b$process_terms_obs)){
            b$process_terms_obs[[i]]$map <- as.integer(b$process_terms_obs[[i]]$map)
            b$process_terms_obs[[i]]$obsIndx <- as.integer(b$process_terms_obs[[i]]$obsIndx)
            b$process_terms_obs[[i]]$obsIndxLU <- matrix(as.integer(b$process_terms_obs[[i]]$obsIndxLU), ncol = 2L)
            b$process_terms_obs[[i]]$node_nobs <- as.integer(b$process_terms_obs[[i]]$node_nobs)
            b$process_terms_obs[[i]]$n_obs <- as.integer(b$process_terms_obs[[i]]$n_obs)
            b$process_terms_obs[[i]]$n_node <- as.integer(b$process_terms_obs[[i]]$n_node)
            if(!is.null(b$process_terms_obs[[i]]$x))
                storage.mode(b$process_terms_obs[[i]]$x) <- "double"
        }

        if(identical(b$residual_model$type, "group_ig_variance") && is.null(residual$starting) && chains > 1L){
            active <- b$residual_model$tuning > 0
            b$residual_model$starting[active] <- vapply(
                b$residual_model$starting[active], log_jitter, numeric(1),
                dispersion = chain_control$dispersion
            )
        }
        if(identical(b$residual_model$type, "scaled_variance") && is.null(residual$starting) && chains > 1L){
            active <- b$residual_model$tuning > 0
            b$residual_model$starting[active] <- vapply(
                b$residual_model$starting[active], log_jitter, numeric(1),
                dispersion = chain_control$dispersion
            )
        }

        class(b) <- "stLMM_backend"
        b
    }

    run_collapsed_sampler <- function(backend_i){
        sampler_start <- proc.time()
        out <- .Call("stLMM_collapsed_sampler", backend_i, PACKAGE = "stLMM")
        sampler_time <- proc.time() - sampler_start
        out$timing <- list(sampler = sampler_time)
        out
    }

    finalize_stLMM_output <- function(out, backend_i, chain = NULL){
        process_names_i <- vapply(backend_i$process_terms, `[[`, character(1), "name")

        if(is.null(out$term_description))
            out$term_description <- backend_i$term_description_meta

        if(!is.null(out$beta_samples) && is.matrix(out$beta_samples))
            colnames(out$beta_samples) <- x_names

        if(!is.null(out$alpha_samples) && is.matrix(out$alpha_samples) && ncol(out$alpha_samples) > 0L){
            z_names <- colnames(Z)
            if(is.null(z_names) || length(z_names) != ncol(out$alpha_samples))
                z_names <- paste0("alpha_", seq_len(ncol(out$alpha_samples)))
            colnames(out$alpha_samples) <- z_names
        }

        if(!is.null(out$sigma_sq_re_samples) && is.matrix(out$sigma_sq_re_samples) && ncol(out$sigma_sq_re_samples) > 0L)
            colnames(out$sigma_sq_re_samples) <- paste0(re_names, "_sigma_sq")
        out$iid_sigma_sq_samples <- out$sigma_sq_re_samples

        if(!is.null(out$sigma_sq_samples) && is.matrix(out$sigma_sq_samples) && ncol(out$sigma_sq_samples) > 0L)
            colnames(out$sigma_sq_samples) <- paste0(process_names_i, "_sigma_sq")

        if(!is.null(out$term_param_accept))
            names(out$term_param_accept) <- process_names_i

        if(!is.null(out$theta_samples) && is.matrix(out$theta_samples) && ncol(out$theta_samples) > 0L){
            theta_sample_names <- character(0)
            for(i in seq_along(backend_i$process_terms)){
                theta_names_i <- backend_i$process_terms[[i]]$theta_names %||% character(0)
                if(length(theta_names_i))
                    theta_sample_names <- c(theta_sample_names, paste0(process_names_i[i], "_", theta_names_i))
            }
            if(length(theta_sample_names) == ncol(out$theta_samples))
                colnames(out$theta_samples) <- theta_sample_names
        }

        if(is_pg_likelihood || fixed_residual_variance || sampled_residual_variance)
            out$tau_sq_samples <- NULL

        if(sampled_residual_variance && !is.null(out$residual_variance_samples) &&
           is.matrix(out$residual_variance_samples)){
            if(identical(backend_i$residual_model$type, "group_ig_variance"))
                colnames(out$residual_variance_samples) <- backend_i$residual_model$groups
            else if(identical(backend_i$residual_model$type, "scaled_variance"))
                colnames(out$residual_variance_samples) <- backend_i$residual_model$parameter_names
        }

        out$samples <- list(
            beta = out$beta_samples,
            alpha = if(is.matrix(out$alpha_samples) && ncol(out$alpha_samples) > 0L) out$alpha_samples else NULL,
            tau_sq = out$tau_sq_samples,
            residual_variance = if(is.matrix(out$residual_variance_samples) && ncol(out$residual_variance_samples) > 0L)
                out$residual_variance_samples else NULL,
            iid_sigma_sq = if(is.matrix(out$iid_sigma_sq_samples) && ncol(out$iid_sigma_sq_samples) > 0L) out$iid_sigma_sq_samples else NULL,
            sigma_sq = if(is.matrix(out$sigma_sq_samples) && ncol(out$sigma_sq_samples) > 0L) out$sigma_sq_samples else NULL,
            theta = if(is.matrix(out$theta_samples) && ncol(out$theta_samples) > 0L) out$theta_samples else NULL
        )

        out$sigma_sq_re_samples <- NULL

        if(is.matrix(out$alpha_samples) && ncol(out$alpha_samples) == 0L)
            out$alpha_samples <- NULL
        if(is.matrix(out$iid_sigma_sq_samples) && ncol(out$iid_sigma_sq_samples) == 0L)
            out$iid_sigma_sq_samples <- NULL
        if(is.matrix(out$sigma_sq_samples) && ncol(out$sigma_sq_samples) == 0L)
            out$sigma_sq_samples <- NULL
        if(is.matrix(out$theta_samples) && ncol(out$theta_samples) == 0L)
            out$theta_samples <- NULL
        if(is.matrix(out$residual_variance_samples) && ncol(out$residual_variance_samples) == 0L)
            out$residual_variance_samples <- NULL

        if(!is.null(out$w_samples) && is.matrix(out$w_samples) &&
           nrow(out$w_samples) > 0L && length(backend_i$process_terms)){
            w_samples_stacked <- out$w_samples
            out$w_samples_stacked <- w_samples_stacked
            out$w_samples_ordered <- unstack_w_samples(
                w_samples_stacked = w_samples_stacked,
                process_terms = backend_i$process_terms,
                term_description = out$term_description,
                graphs = backend_i$graphs,
                user_order = FALSE
            )
            out$w_samples <- unstack_w_samples(
                w_samples_stacked = w_samples_stacked,
                process_terms = backend_i$process_terms,
                term_description = out$term_description,
                graphs = backend_i$graphs,
                user_order = TRUE
            )
            out$samples$w <- out$w_samples
        } else {
            out$w_samples <- NULL
            out$w_samples_ordered <- NULL
            out$w_samples_stacked <- NULL
            out$recover_iter <- NULL
            out$samples$w <- NULL
        }

        out$backend <- backend_i
        if(!is.null(chain))
            out$chain <- as.integer(chain)
        class(out) <- c("stLMM", class(out))
        out
    }

    if(chains == 1L){
        status_message("Calling C++ sampler.")
        out <- run_collapsed_sampler(backend)
        return(finalize_stLMM_output(out, backend))
    }

    chain_seeds <- sample.int(.Machine$integer.max, chains)
    chain_backends <- vector("list", chains)
    chain_rng_states <- vector("list", chains)
    for(chain in seq_len(chains)){
        set.seed(chain_seeds[chain])
        chain_backends[[chain]] <- build_chain_backend(chain)
        chain_rng_states[[chain]] <- get(".Random.seed", envir = .GlobalEnv)
    }

    chain_fits <- vector("list", chains)
    for(chain in seq_len(chains)){
        backend_i <- chain_backends[[chain]]
        assign(".Random.seed", chain_rng_states[[chain]], envir = .GlobalEnv)
        if(isTRUE(verbose))
            message("chain ", chain, " of ", chains)
        out_i <- run_collapsed_sampler(backend_i)
        chain_fits[[chain]] <- finalize_stLMM_output(out_i, backend_i, chain = chain)
    }
    names(chain_fits) <- paste0("chain_", seq_len(chains))

    sampler_timing <- do.call(
        rbind,
        lapply(chain_fits, function(z) as.numeric(z$timing$sampler))
    )
    colnames(sampler_timing) <- names(chain_fits[[1L]]$timing$sampler)

    out <- list(
        chains = chain_fits,
        n_chains = as.integer(chains),
        chain_control = chain_control,
        call = match.call(),
        backend = backend,
        term_description = chain_fits[[1]]$term_description,
        timing = list(
            sampler_by_chain = sampler_timing,
            sampler_total = structure(colSums(sampler_timing), class = "proc_time")
        )
    )
    class(out) <- "stLMM_chains"
    out

}
