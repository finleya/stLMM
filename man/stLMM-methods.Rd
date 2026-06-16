\name{stLMM-methods}
\alias{print.stLMM}
\alias{summary.stLMM}
\alias{plot.stLMM}
\alias{print.stLMM_chains}
\alias{summary.stLMM_chains}
\alias{plot.stLMM_chains}
\alias{print.stLMM_recovery}
\alias{summary.stLMM_recovery}
\alias{summary.stLMM_recovery_chains}
\alias{plot.stLMM_recovery}
\alias{plot.stLMM_recovery_chains}
\alias{print.stLMM_prediction}
\alias{summary.stLMM_prediction}
\alias{plot.stLMM_prediction}
\alias{plot.stLMM_prediction_chains}
\title{Print, summarize, and plot stLMM objects}
\usage{
\method{print}{stLMM}(x, ...)
\method{summary}{stLMM}(object, probs = c(0.025, 0.5, 0.975),
  burn = 0L, ...)
\method{plot}{stLMM}(x, type = c("trace", "density"), parameters = NULL,
  max_parameters = 12L, burnin = 0L, thin = 1L, ...)
\method{print}{stLMM_chains}(x, ...)
\method{summary}{stLMM_chains}(object, probs = c(0.025, 0.5, 0.975),
  diagnostics = TRUE, include_w = FALSE, burn = 0L, ...)
\method{plot}{stLMM_chains}(x, type = c("trace", "density"),
  parameters = NULL, max_parameters = 12L, include_w = FALSE,
  n_col = NULL, chain_colors = NULL, lwd = 1, burnin = 0L,
  thin = 1L, ...)
\method{print}{stLMM_recovery}(x, ...)
\method{summary}{stLMM_recovery}(object, probs = c(0.025, 0.5, 0.975),
  burn = 0L, include_w = FALSE, max_w = 20L, ...)
\method{summary}{stLMM_recovery_chains}(object,
  probs = c(0.025, 0.5, 0.975), diagnostics = TRUE,
  include_w = FALSE, burn = 0L, ...)
\method{plot}{stLMM_recovery}(x, term = NULL, nodes = NULL,
  type = c("trace", "density", "fitted"), observed = NULL,
  max_nodes = 12L, burnin = 0L, thin = 1L, ...)
\method{plot}{stLMM_recovery_chains}(x, type = c("trace", "density"),
  parameters = NULL, max_parameters = 12L, include_w = FALSE, ...)
\method{print}{stLMM_prediction}(x, ...)
\method{summary}{stLMM_prediction}(object, probs = c(0.025, 0.5, 0.975),
  include_y = !is.null(object$y_samples), ...)
\method{plot}{stLMM_prediction}(x, sample = c("mu", "y"),
  type = c("interval", "density", "scatter"), rows = NULL,
  observed = NULL, probs = c(0.025, 0.5, 0.975), max_rows = 200L,
  burnin = 0L, thin = 1L, ...)
\method{plot}{stLMM_prediction_chains}(x, sample = c("mu", "y"),
  type = c("interval", "density"), rows = NULL,
  probs = c(0.025, 0.5, 0.975), max_rows = 200L,
  burnin = 0L, thin = 1L, ...)
}
\arguments{
\item{x, object}{An \code{stLMM}, \code{stLMM_recovery}, or
\code{stLMM_prediction} object, or the corresponding multi-chain object.}
\item{probs}{Posterior quantiles to report. Prediction interval plots expect
three probabilities: lower, center, and upper.}
\item{type}{Plot type. For fit objects, \code{"trace"} or \code{"density"}.
For recovery objects, \code{"trace"}, \code{"density"}, or \code{"fitted"}.
For prediction objects, \code{"interval"}, \code{"density"}, or
\code{"scatter"}.}
\item{parameters}{Optional parameter names to plot for an \code{stLMM} fit.}
\item{max_parameters}{Maximum number of parameters to plot by default.}
\item{burn}{Number of initial posterior draws to discard before computing
fit summaries and multi-chain diagnostics. The default \code{0L} uses all
retained draws. This affects only the summary call; it does not modify the
stored posterior samples.}
\item{diagnostics}{Logical; for multi-chain summaries, compute \pkg{coda}
diagnostics such as potential scale reduction and effective sample size.}
\item{burnin}{Number of initial plotted draws to discard. This affects only
the plotting call; it does not modify the stored posterior samples.}
\item{thin}{Keep every \code{thin}-th plotted draw after \code{burnin}. This
affects only the plotting call.}
\item{include_w}{Logical; include summaries of saved or recovered latent process
samples. Defaults to \code{FALSE} because these can be large.}
\item{n_col}{Optional number of columns for multi-chain fit trace or density
plot panels. Defaults to a compact base graphics layout.}
\item{chain_colors}{Optional vector of colors for multi-chain fit plots. Values
are recycled to the number of chains.}
\item{lwd}{Line width for multi-chain fit trace or density plots.}
\item{max_w}{Maximum number of latent process nodes summarized when
\code{include_w = TRUE}.}
\item{term}{Recovered process term to plot. Defaults to the first recovered
term.}
\item{nodes}{Latent process nodes to plot.}
\item{max_nodes}{Maximum number of recovered nodes to plot by default.}
\item{include_y}{Logical; include posterior predictive observation summaries
when \code{y_samples} were simulated.}
\item{sample}{Prediction sample matrix to plot: \code{"mu"} or \code{"y"}.}
\item{rows}{Prediction rows to plot.}
\item{observed}{Observed values for recovery fitted-value scatter plots or
prediction scatter plots. Recovery plots use the fitted response by default.}
\item{max_rows}{Maximum number of prediction rows plotted by default.}
\item{...}{Additional arguments passed to base plotting functions.}
}
\description{
S3 methods for compact inspection of fitted, recovered, and prediction objects.
Summaries report posterior means, standard deviations, and quantiles for active
sample blocks only. Plot methods use base graphics and choose conservative
defaults to avoid drawing very large posterior or prediction objects by
accident.
}
\value{
Print methods return \code{x} invisibly. Summary methods return summary objects.
Plot methods return \code{x} invisibly.
}
