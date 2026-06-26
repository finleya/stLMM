#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1 || length(args) > 2) {
  stop(
    "Usage: Rscript scripts/extract_article_r.R ARTICLE.qmd [OUTPUT.R]",
    call. = FALSE
  )
}

input <- args[[1]]
output <- if (length(args) == 2) {
  args[[2]]
} else {
  sub("\\.qmd$", ".R", input)
}

if (!file.exists(input)) {
  stop("Input file does not exist: ", input, call. = FALSE)
}

tmp <- tempfile(fileext = ".R")
invisible(knitr::purl(input, output = tmp, quiet = TRUE, documentation = 0))

code <- readLines(tmp, warn = FALSE)

# Quarto chunk options are valid comments, but they distract from the
# downloadable article script.
code <- code[!grepl("^#\\|", code)]

writeLines(code, output)
message("Wrote ", normalizePath(output, mustWork = FALSE))
