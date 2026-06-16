\name{fitted.stLMM}
\alias{fitted.stLMM}
\title{Posterior fitted values for an stLMM fit}
\usage{
\method{fitted}{stLMM}(object, summary = TRUE,
  sub_sample = list(start = 1L, thin = 1L), scale = c("response", "link"),
  ...)
}
\arguments{
\item{object}{An object returned by \code{stLMM} or \code{recover}.}
\item{summary}{Logical; if \code{TRUE}, return posterior mean fitted values.
If \code{FALSE}, return the posterior fitted-value sample matrix.}
\item{sub_sample}{List with optional \code{start} and \code{thin} entries used
to subset posterior draws. For models with structured process terms, this
subsets the saved or recovered latent-process draw indices in
\code{object$recover_iter}.}
\item{scale}{For \code{family = "binomial"}, return fitted probabilities on
the \code{"response"} scale or the linear predictor on the \code{"link"} scale.
For \code{family = "negative_binomial"}, return fitted mean counts on the
\code{"response"} scale or the log mean on the \code{"link"} scale. For
Gaussian models these scales are identical.}
\item{...}{Currently unused.}
}
\value{
If \code{summary = TRUE}, a numeric vector of posterior mean fitted values for
the original fitted data rows, including rows whose response was missing at fit
time. If \code{summary = FALSE}, a matrix with posterior draws in rows and
original data rows in columns. The matrix has a \code{draw_index}
attribute giving the original MCMC iteration indices used and, for
\code{summary = FALSE} Polya-Gamma fits, a \code{scale} attribute.
}
\description{
Computes posterior fitted values for the original data used to fit an
\code{stLMM} model. These fitted values include all model components that were
fit: formula offsets, fixed effects, explicit grouped random effects, and
saved or recovered structured latent process effects. Binomial fitted values default to
posterior fitted probabilities, and negative-binomial fitted values default to
posterior fitted mean counts.
}
\details{
This method does not predict at new locations. For models with structured
process terms, fitted values require saved or recovered latent process samples.
For Polya-Gamma process models, \code{stLMM()} saves in-chain process draws by
default as \code{save_process = list(start = 1, thin = 1)} and
\code{recover()} subsets them. For Gaussian process models, call
\code{recover()} on the fitted object to make those samples available.
}
\examples{
\dontrun{
fit <- stLMM(...)
rec <- recover(fit)
mu_hat <- fitted(rec)
mu_draws <- fitted(rec, summary = FALSE)
}
}
