\name{log_lik}
\alias{log_lik}
\alias{log_lik.stLMM}
\alias{log_lik.stLMM_chains}
\alias{log_lik.stLMM_recovery}
\alias{log_lik.stLMM_recovery_chains}
\alias{waic}
\alias{waic.default}
\alias{waic.stLMM}
\alias{waic.stLMM_chains}
\alias{waic.stLMM_recovery}
\alias{waic.stLMM_recovery_chains}
\title{Pointwise log likelihood and WAIC}
\usage{
log_lik(object, ...)

\method{log_lik}{stLMM}(object,
  sub_sample = list(start = 1L, thin = 1L), ...)

\method{log_lik}{stLMM_chains}(object,
  sub_sample = list(start = 1L, thin = 1L), ...)

\method{log_lik}{stLMM_recovery}(object,
  sub_sample = list(start = 1L, thin = 1L), ...)

\method{log_lik}{stLMM_recovery_chains}(object,
  sub_sample = list(start = 1L, thin = 1L), ...)

waic(object, ...)

\method{waic}{stLMM}(object,
  sub_sample = list(start = 1L, thin = 1L), ...)

\method{waic}{stLMM_chains}(object,
  sub_sample = list(start = 1L, thin = 1L), ...)

\method{waic}{stLMM_recovery}(object,
  sub_sample = list(start = 1L, thin = 1L), ...)

\method{waic}{stLMM_recovery_chains}(object,
  sub_sample = list(start = 1L, thin = 1L), ...)
}
\arguments{
\item{object}{An \code{stLMM} fit, recovered \code{stLMM_recovery} object,
or corresponding multi-chain object. Models with structured process terms
require saved or recovered latent process samples.}
\item{sub_sample}{A list with optional integer entries \code{start} and
\code{thin}. For ordinary fits these select posterior draw indices. For
recovered objects and Polya-Gamma process fits with saved process draws, these
select saved or recovered draws through \code{recover_iter}.}
\item{...}{Additional arguments. For \code{waic()}, these are passed to
\code{loo::waic()}.}
}
\description{
\code{log_lik()} returns the pointwise log-likelihood matrix used for
model-fit diagnostics. \code{waic()} computes Watanabe-Akaike information
criterion using \code{loo::waic()}.
}
\details{
The returned log-likelihood matrix has posterior draws in rows and observed
response rows in columns. Missing response rows are excluded because they do
not contribute to the fitted likelihood.

The calculation is conditional on the latent model expression. For Gaussian
models, each entry is evaluated with the fitted linear predictor and the
draw-specific residual standard deviation. For binomial models, entries are
computed with the inverse-logit probability and fitted trial count. For
fixed-size negative-binomial models, entries are computed with the log-mean
linear predictor and the fitted size parameter.

For models with structured process terms, \code{log_lik()} requires saved or
recovered latent process samples so that the latent linear predictor can be
evaluated draw by draw. Gaussian process models should therefore be passed to
\code{recover()} before calling \code{log_lik()} or \code{waic()}. Polya-Gamma
process fits can be used directly when process draws were saved during
\code{stLMM()}.

\code{waic()} requires the optional \pkg{loo} package. \pkg{loo} is suggested,
not imported, so ordinary model fitting does not require it.
}
\value{
\code{log_lik()} returns a numeric matrix. \code{waic()} returns the object
produced by \code{loo::waic()}.
}
\examples{
\dontrun{
fit <- stLMM(y ~ x, data = dat, n_samples = 1000)
ll <- log_lik(fit, sub_sample = list(start = 501, thin = 2))
waic(fit, sub_sample = list(start = 501, thin = 2))

fit_p <- stLMM(y ~ x + ar1(time), data = dat, n_samples = 1000)
rec_p <- recover(fit_p, sub_sample = list(start = 501, thin = 2))
waic(rec_p)
}
}
\seealso{\code{\link{stLMM}}, \code{\link{recover}},
\code{\link{fitted.stLMM}}, \code{\link{predict.stLMM_recovery}}}
