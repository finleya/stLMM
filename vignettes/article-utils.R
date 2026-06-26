article_helper_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)

source(file.path(article_helper_dir, "utils.R"))
source(file.path(article_helper_dir, "article-cache.R"))

show_table <- function(x, caption = NULL, ...){
  knitr::kable(
    x,
    format = "html",
    escape = FALSE,
    caption = caption,
    table.attr = "class=\"stlmm-table\"",
    ...
  ) |>
    kableExtra::kable_styling(
      bootstrap_options = c("striped", "hover", "condensed", "responsive"),
      full_width = FALSE,
      position = "left",
      font_size = 13
    )
}

posterior_samples <- function(prediction, sample_name = "mu_samples") {
  if (inherits(prediction, "stLMM_prediction_chains")) {
    return(do.call(rbind, lapply(prediction$chains, `[[`, sample_name)))
  }

  prediction[[sample_name]]
}

summarize_draw_matrix <- function(draw_matrix, prefix = "") {
  tibble::tibble(
    prediction_row = seq_len(ncol(draw_matrix)),
    "{prefix}mean" := colMeans(draw_matrix, na.rm = TRUE),
    "{prefix}median" := apply(draw_matrix, 2, median, na.rm = TRUE),
    "{prefix}lower" := apply(draw_matrix, 2, quantile, probs = 0.025, na.rm = TRUE),
    "{prefix}upper" := apply(draw_matrix, 2, quantile, probs = 0.975, na.rm = TRUE),
    "{prefix}sd" := apply(draw_matrix, 2, sd, na.rm = TRUE)
  )
}

fmt_int <- function(x) {
  format(x, big.mark = ",", scientific = FALSE)
}

fmt_num <- function(x, digits = 1) {
  format(round(x, digits), big.mark = ",", scientific = FALSE)
}

manifest_value <- function(manifest, key = NULL) {
  if (is.null(key)) {
    key <- manifest
    manifest <- get("manifest", envir = parent.frame())
  }

  manifest$value[match(key, manifest$key)]
}
