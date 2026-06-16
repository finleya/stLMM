get_cor_models <- function() {
  .Call("get_cor_models", PACKAGE = "stLMM")
}

build_cor_model_registry <- function() {

  raw <- get_cor_models()

  out <- list()

  for(name in names(raw)) {

    entry <- raw[[name]]

    out[[name]] <- list(
      name = name,
      theta_names = entry$names,
      theta_types = entry$types,
      theta_domains = entry$domains,
      distance_mode = entry$distance_mode,
      n_theta = length(entry$names)
    )
  }

  class(out) <- "cor_model_registry"

  out
}
