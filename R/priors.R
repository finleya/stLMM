make_stLMM_prior <- function(family, parameters = numeric(0), support = NULL,
                             scale = "value"){
  structure(
    list(
      family = family,
      parameters = parameters,
      support = support,
      scale = scale
    ),
    class = "stLMM_prior"
  )
}

fixed <- function(value){
  if(!is.numeric(value) || length(value) != 1L || is.na(value) || !is.finite(value))
    stop("error: fixed() value must be a finite numeric scalar")
  structure(
    list(value = as.double(value)),
    class = "stLMM_fixed_parameter"
  )
}

flat <- function(){
  make_stLMM_prior(
    "flat",
    numeric(0)
  )
}

normal <- function(mean = 0, sd = 1){
  make_stLMM_prior(
    "normal",
    list(mean = mean, sd = sd)
  )
}

ig <- function(shape, scale){
  make_stLMM_prior(
    "ig",
    c(shape = shape, scale = scale)
  )
}

uniform <- function(lower, upper){
  make_stLMM_prior(
    "uniform",
    numeric(0),
    support = c(lower = lower, upper = upper)
  )
}

log_normal <- function(meanlog, sdlog, support = NULL){
  make_stLMM_prior(
    "log_normal",
    c(meanlog = meanlog, sdlog = sdlog),
    support = support
  )
}

gamma_dist <- function(shape, rate, support = NULL){
  make_stLMM_prior(
    "gamma_dist",
    c(shape = shape, rate = rate),
    support = support
  )
}

half_normal <- function(scale){
  make_stLMM_prior(
    "half_normal",
    c(scale = scale),
    scale = "sd"
  )
}

half_t <- function(df, scale){
  make_stLMM_prior(
    "half_t",
    c(df = df, scale = scale),
    scale = "sd"
  )
}

beta_dist <- function(shape1, shape2){
  make_stLMM_prior(
    "beta_dist",
    c(shape1 = shape1, shape2 = shape2),
    support = c(lower = 0, upper = 1)
  )
}
