article_cache_root <- function() {
  Sys.getenv("STLMM_ARTICLE_CACHE_DIR", unset = file.path("cache", "articles"))
}

article_cache_file <- function(article_id, cache_name) {
  file.path(article_cache_root(), article_id, paste0(cache_name, ".rds"))
}

run_long_articles <- function() {
  identical(tolower(Sys.getenv("STLMM_RUN_LONG_ARTICLES", unset = "false")), "true")
}

use_article_cache <- function(article_id, cache_name, code) {
  cache_file <- article_cache_file(article_id, cache_name)

  if (!run_long_articles() && file.exists(cache_file)) {
    return(readRDS(cache_file))
  }

  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
  out <- force(code)
  saveRDS(out, cache_file)
  out
}
