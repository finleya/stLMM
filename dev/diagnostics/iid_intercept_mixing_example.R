## Teaching example for iid random-intercept diagnostics.
##
## This reproduces the small example from GitHub issue #1 and compares MCMC
## diagnostics for:
##   1. the raw intercept decomposition, beta_0 and alpha_j;
##   2. the group-specific intercepts, beta_0 + alpha_j.
##
## Run from the package root with:
##   Rscript dev/diagnostics/iid_intercept_mixing_example.R

suppressPackageStartupMessages({
  library(coda)
})

if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", quiet = TRUE)
} else {
  suppressPackageStartupMessages(library(stLMM))
}

simulate_issue_data <- function(seed = 1) {
  set.seed(seed)

  n_group <- 5
  n_per_group <- 10

  dat <- data.frame(
    group = factor(rep(letters[1:n_group], each = n_per_group)),
    x = rep(seq(-1, 1, length.out = n_per_group), n_group)
  )

  sigma_sq_alpha <- 0.5^2
  alpha_true <- rnorm(n_group, sd = sqrt(sigma_sq_alpha))
  names(alpha_true) <- levels(dat$group)

  beta <- c("(Intercept)" = 0.5, x = -0.5)
  tau_sq <- 0.25^2

  mu <- beta[["(Intercept)"]] +
    beta[["x"]] * dat$x +
    alpha_true[as.character(dat$group)]

  dat$y <- mu + rnorm(nrow(dat), sd = sqrt(tau_sq))
  dat
}

iid_alpha_columns <- function(x) {
  alpha_cols <- grep("^iid_1_", colnames(x), value = TRUE)
  setdiff(alpha_cols, "iid_1_sigma_sq")
}

add_group_intercepts <- function(chain) {
  chain <- as.matrix(chain)
  alpha_cols <- iid_alpha_columns(chain)
  group_intercepts <- sweep(
    chain[, alpha_cols, drop = FALSE],
    1L,
    chain[, "(Intercept)"],
    "+"
  )
  colnames(group_intercepts) <- sub("^iid_1_", "beta0_plus_iid_1_", alpha_cols)

  coda::mcmc(
    cbind(
      chain[, c("(Intercept)", "x", alpha_cols, "iid_1_sigma_sq", "tau_sq"), drop = FALSE],
      group_intercepts
    )
  )
}

diagnostics_table <- function(chains, parameters) {
  chains <- lapply(chains, function(chain) {
    chain <- as.matrix(chain)
    coda::mcmc(chain[, parameters, drop = FALSE])
  })
  chains <- do.call(coda::mcmc.list, chains)

  rhat <- coda::gelman.diag(
    chains,
    autoburnin = FALSE,
    multivariate = FALSE
  )$psrf[, "Point est."]

  ess <- coda::effectiveSize(chains)

  out <- data.frame(
    parameter = parameters,
    rhat = unname(rhat[parameters]),
    effective_size = unname(ess[parameters]),
    row.names = NULL,
    check.names = FALSE
  )

  out[order(out$effective_size), ]
}

plot_trace_set <- function(chains, parameters, file, title) {
  chain_values <- lapply(chains, function(chain) {
    chain <- as.matrix(chain)
    chain[, parameters, drop = FALSE]
  })

  grDevices::png(file, width = 1800, height = 1100, res = 160)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  })

  n_col <- 2L
  n_row <- ceiling(length(parameters) / n_col)
  chain_col <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")

  graphics::par(
    mfrow = c(n_row, n_col),
    mar = c(3.1, 4.0, 2.2, 1.0),
    oma = c(0, 0, 3.0, 0)
  )

  for (parameter in parameters) {
    y <- lapply(chain_values, function(x) x[, parameter])
    y_range <- range(unlist(y, use.names = FALSE), finite = TRUE)

    graphics::plot(
      y[[1L]],
      type = "l",
      col = chain_col[1L],
      ylim = y_range,
      xlab = "post-burn draw",
      ylab = parameter,
      main = parameter
    )

    if (length(y) > 1L) {
      for (chain in 2:length(y)) {
        graphics::lines(
          y[[chain]],
          col = chain_col[(chain - 1L) %% length(chain_col) + 1L]
        )
      }
    }
  }

  if (length(parameters) %% n_col != 0L)
    graphics::plot.new()

  graphics::mtext(title, outer = TRUE, cex = 1.1, font = 2)
  invisible(file)
}

dat <- simulate_issue_data(seed = 1)

fit <- stLMM(
  y ~ x + iid(group),
  data = dat,
  priors = list(
    resid = list(tau_sq = half_t(df = 3, scale = 0.5)),
    iid_1 = list(sigma_sq = ig(shape = 3, scale = 0.25))
  ),
  n_samples = 2500,
  chains = 4,
  chain_control = list(seed = 0),
  verbose = FALSE
)

burn <- 500
chains <- as_mcmc(fit, burn = burn)
chains_with_derived <- lapply(chains, add_group_intercepts)

first_chain <- as.matrix(chains_with_derived[[1L]])
alpha_cols <- iid_alpha_columns(first_chain)
group_intercept_cols <- grep("^beta0_plus_iid_1_", colnames(first_chain), value = TRUE)

raw_intercept_parameters <- c("(Intercept)", alpha_cols)
comparison_parameters <- c(
  "(Intercept)",
  alpha_cols,
  group_intercept_cols,
  "x",
  "iid_1_sigma_sq",
  "tau_sq"
)

raw_trace_file <- file.path(
  "dev",
  "diagnostics",
  "iid_intercept_mixing_raw_trace.png"
)
derived_trace_file <- file.path(
  "dev",
  "diagnostics",
  "iid_intercept_mixing_group_intercepts_trace.png"
)

plot_trace_set(
  chains_with_derived,
  raw_intercept_parameters,
  raw_trace_file,
  "Raw intercept decomposition: beta_0 and alpha_j"
)
plot_trace_set(
  chains_with_derived,
  group_intercept_cols,
  derived_trace_file,
  "Group-specific intercepts: beta_0 + alpha_j"
)

cat("\nRaw intercept decomposition: beta_0 and alpha_j\n")
print(
  diagnostics_table(chains_with_derived, raw_intercept_parameters),
  row.names = FALSE,
  digits = 3
)

cat("\nSame fit, but diagnosing beta_0 + alpha_j\n")
print(
  diagnostics_table(chains_with_derived, group_intercept_cols),
  row.names = FALSE,
  digits = 3
)

cat("\nBroader comparison\n")
print(
  diagnostics_table(chains_with_derived, comparison_parameters),
  row.names = FALSE,
  digits = 3
)

cat("\nTrace plots written to:\n")
cat("  ", raw_trace_file, "\n", sep = "")
cat("  ", derived_trace_file, "\n", sep = "")

cat(
  "\nInterpretation:\n",
  "  beta_0 and alpha_j can have poor ESS/Rhat because the sampler updates\n",
  "  them separately and the posterior strongly couples the global intercept\n",
  "  with the average random-intercept level. The group-specific intercepts\n",
  "  beta_0 + alpha_j are usually the better diagnostics for fitted means and\n",
  "  convergence of the identifiable intercept quantities.\n",
  sep = ""
)
