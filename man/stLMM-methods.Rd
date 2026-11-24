\name{stLMM-methods}
\alias{print.stLMM}
\alias{summary.stLMM}
\alias{plot.stLMM}
\alias{print.stLMM_chains}
\alias{summary.stLMM_chains}
\alias{print.summary_stLMM_combined}
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
  burn = 0L, parameters = NULL, type = c("parameters", "combined"),
  coefficients = NULL, terms = NULL, diagnostics = FALSE,
  diagnostic_method = c("coda", "rank"), ...)
\method{plot}{stLMM}(x, type = c("trace", "density"), parameters = NULL,
  max_parameters = 12L, burnin = 0L, thin = 1L, ...)
\method{print}{stLMM_chains}(x, ...)
\method{summary}{stLMM_chains}(object, probs = c(0.025, 0.5, 0.975),
  diagnostics = TRUE, include_w = FALSE, burn = 0L,
  parameters = NULL, type = c("parameters", "combined"),
  coefficients = NULL, terms = NULL,
  diagnostic_method = c("coda", "rank"), ...)
\method{print}{summary_stLMM_combined}(x,
  digits = max(3L, getOption("digits") - 3L), ...)
\method{plot}{stLMM_chains}(x, type = c("trace", "density"),
  parameters = NULL, max_parameters = 12L, include_w = FALSE,
  n_col = NULL, chain_colors = NULL, lwd = 1, burnin = 0L,
  thin = 1L, ...)
\method{print}{stLMM_recovery}(x, ...)
\method{summary}{stLMM_recovery}(object, probs = c(0.025, 0.5, 0.975),
  burn = 0L, parameters = NULL, include_w = FALSE, max_w = 20L,
  type = c("parameters", "combined"), coefficients = NULL, terms = NULL,
  diagnostics = FALSE, diagnostic_method = c("coda", "rank"), ...)
\method{summary}{stLMM_recovery_chains}(object,
  probs = c(0.025, 0.5, 0.975), diagnostics = TRUE,
  include_w = FALSE, burn = 0L, parameters = NULL,
  type = c("parameters", "combined"), coefficients = NULL, terms = NULL,
  diagnostic_method = c("coda", "rank"), ...)
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
\item{type}{For fit and recovery summaries, \code{"parameters"} reports the
ordinary parameter summaries; \code{"combined"} summarizes matching fixed
and selected iid coefficients. For plots of fit objects, \code{"trace"} or \code{"density"}.
For recovery objects, \code{"trace"}, \code{"density"}, or \code{"fitted"}.
For prediction objects, \code{"interval"}, \code{"density"}, or
\code{"scatter"}.}
\item{parameters}{Optional exact parameter names to summarize or plot for an
\code{stLMM} fit. For summaries, unlisted parameters are omitted from the
returned parameter table or tables.}
\item{coefficients}{For \code{type = "combined"}, exact coefficient names such
as \code{"(Intercept)"} or \code{"x"}. The default selects all coefficients
appearing in the selected iid terms. Each requested coefficient must appear
in at least one selected iid term.}
\item{terms}{For \code{type = "combined"}, required iid formula term identifiers
such as \code{c("iid_1", "iid_2")}. A parent-specific \code{by} term is selected
by its formula identifier, e.g., \code{"iid_2"}, which includes its parent blocks.
For \code{cov_model = "independent"}, the formula identifier includes all of
its coefficient blocks; \code{coefficients} selects which coefficients to summarize.
\code{"all"} explicitly selects all terms in an iid-only model. Combined
summaries select iid terms. In a model with structured effects, select the iid
terms explicitly.}
\item{max_parameters}{Maximum number of parameters to plot by default.}
\item{burn}{Number of initial posterior draws to discard before computing
fit summaries and multi-chain diagnostics. The default \code{0L} uses all
retained draws. This affects only the summary call; it does not modify the
stored posterior samples.}
\item{diagnostics}{Logical; compute diagnostics using \code{diagnostic_method}.
Defaults to \code{TRUE} for multiple chains and \code{FALSE} for one chain.
Single-chain diagnostic output reports \code{NA} for R-hat.}
\item{diagnostic_method}{\code{"coda"} reports classical R-hat
and \code{effective_size}. \code{"rank"} requires the optional
\pkg{posterior} package and reports \code{rhat} (rank-normalized split/folded
R-hat), \code{ess_bulk}, and \code{ess_tail}. This selection is available for
both ordinary and combined coefficient summaries.}
\item{digits}{Number of digits used when printing a combined summary.}
\item{burnin}{Number of initial plotted draws to discard. This affects only
the plotting call; it does not modify the stored posterior samples. For
single-chain recovery plots this counts recovered rows, rather than original
iteration numbers. Multi-chain fit and recovery plots use an original
posterior iteration cutoff, as in \code{\link{as_mcmc}}.}
\item{thin}{Keep every \code{thin}-th plotted draw after \code{burnin}. This
affects only the plotting call.}
\item{include_w}{Logical; include summaries of saved or recovered latent process
samples in ordinary parameter summaries. Defaults to \code{FALSE} because these
can be large. This argument does not select contributions to combined coefficients.}
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
\item{...}{Additional arguments, including arguments passed to base plotting functions.}
}
\description{
S3 methods for compact inspection of fitted, recovered, and prediction objects.
For IID-only Gaussian recovery objects, trace and density plots show parameter
and recovered IID draws on the selected recovery grid; \code{parameters} can
select columns as for an ordinary fit. The \code{max_parameters} argument can
be passed through \code{...}, defaulting to \code{max_nodes} for these plots.
Fitted-value plots use the recovered
IID coefficients.
Summaries report posterior means, standard deviations, and quantiles for active
sample blocks only. Plot methods use base graphics and choose conservative
defaults to avoid drawing very large posterior or prediction objects by
accident.
}
\details{
Combined summaries construct each coefficient separately within every retained
iteration and chain. For example, selecting a stand intercept and a parent-specific
plot intercept gives \eqn{\beta_0 + u_{stand} + v_{stand,plot}}. A slope summary
gives the matching fixed slope plus selected random slopes, without multiplication
by covariate values. If no matching fixed coefficient exists, its contribution is
zero. The scale is the linear-predictor coefficient scale for every family;
offsets, other coefficients, and response transformations are not included.
When a covariate also has a structured coefficient process, the combined iid
slope excludes that process contribution. Its intervals and diagnostics describe
the selected iid sum. Use \code{fitted()} or \code{predict()} for the full
linear predictor or response quantity.

Each row describes a represented combination of the contributing grouping
variables. Repeated observations are not averaged, unused factor levels are
omitted, and unrepresented crossed or nested combinations are not created.
For a coefficient absent from some selected terms, only the matching terms
contribute; the contributing identifiers are recorded per row. Group labels
are retained separately, including parent labels for \code{by} terms.

For collapsed Gaussian fits, call \code{recover()} before requesting IID sums.
These sums pair recovered IID draws with fixed-effect draws at the original
\code{recover_iter} indices after \code{burn}. Parameter summaries and diagnostics
on these recovery objects also use that selected grid; summarize the original
fit to use every retained parameter draw. With separate IID sampling or augmented
likelihoods, IID sums use all original fixed and IID draws after \code{burn}, even
when structured draws were saved at fewer iterations. No
draws are simulated and no latent recovery is performed by a summary. Group
membership is read from the original sparse design, including explicit zero
entries for zero-valued slope covariates; do not remove those entries from the
stored fit. Only one combined quantity's draws are materialized at a time.

Diagnostics apply to the actual combined draws, retaining chain identity.
Good diagnostics for a sum do not establish good mixing of its separate
components. The rank method follows Vehtari et al. (2021); \code{ess_tail}
is the smaller ESS for the 5th and 95th percentiles, not an ESS specific to the
requested \code{probs}. Rank diagnostics are \code{NA} for fewer than four
retained iterations, nonfinite draws, or a constant chain. R-hat is also
\code{NA} when only one chain is available. The default summary type is
\code{"parameters"}, with \code{"coda"} diagnostics for multiple chains.
}
\value{
Print methods return \code{x} invisibly. Summary methods return summary objects.
For \code{type = "combined"}, the returned \code{summary_stLMM_combined} list
contains \code{quantities} (coefficient, contributing terms, fixed-effect flag,
\code{group_} columns, a representative original \code{data_row}, and an exact
draw expression), \code{parameters} (posterior summaries), \code{diagnostics},
\code{diagnostic_method}, retained \code{iteration} indices, \code{chain_ids},
and draw counts. Corresponding tables share row names. For ordinary summaries,
opting into the rank method records \code{diagnostic_method = "rank"}.
Plot methods return \code{x} invisibly.
}
\references{
Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., and Bürkner, P.-C. (2021).
Rank-normalization, folding, and localization: An improved R-hat for assessing
convergence of MCMC (with discussion). \emph{Bayesian Analysis}, 16(2), 667--718.
\doi{10.1214/20-BA1221}.

Bürkner, P.-C., Gabry, J., Kay, M., and Vehtari, A. \emph{posterior: Tools for
Working with Posterior Distributions}. R package.
\url{https://mc-stan.org/posterior/}.
}
\examples{
\dontrun{
# fit has iid(stand) + iid(plot, ~ 1 + x, by = stand)
summary(fit, type = "combined", coefficients = "(Intercept)",
        terms = c("iid_1", "iid_2"), burn = 1000,
        diagnostic_method = "rank")
summary(fit, type = "combined", coefficients = "x",
        terms = "iid_2", burn = 1000, diagnostic_method = "rank")
summary(fit, parameters = c("(Intercept)", "x"), burn = 1000,
        diagnostic_method = "rank")
}
}
