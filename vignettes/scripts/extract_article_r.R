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

old_extract <- Sys.getenv("STLMM_EXTRACT_R", unset = NA_character_)
on.exit({
  if (is.na(old_extract)) {
    Sys.unsetenv("STLMM_EXTRACT_R")
  } else {
    Sys.setenv(STLMM_EXTRACT_R = old_extract)
  }
}, add = TRUE)

Sys.setenv(STLMM_EXTRACT_R = "true")

tmp <- tempfile(fileext = ".R")
knitr::purl(input, output = tmp, quiet = TRUE, documentation = 0)

code <- readLines(tmp, warn = FALSE)

# Quarto chunk options are valid comments, but they distract from the
# downloadable teaching script.
code <- code[!grepl("^#\\|", code)]

writeLines(code, output)
message("Wrote ", normalizePath(output, mustWork = FALSE))
