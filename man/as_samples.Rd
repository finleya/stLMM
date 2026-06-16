\name{as_samples}
\alias{as_samples}
\alias{as_samples.stLMM}
\alias{as_samples.stLMM_chains}
\alias{as_samples.stLMM_recovery}
\alias{as_samples.stLMM_recovery_chains}
\alias{as_samples.stLMM_prediction}
\alias{as_samples.stLMM_prediction_chains}
\alias{as_samples.stLMM_fitted_chains}
\alias{as_samples.matrix}
\title{Collect stLMM posterior samples in a data frame}
\usage{
as_samples(object, ...)
}
\arguments{
\item{object}{An \code{stLMM}, \code{stLMM_chains}, \code{stLMM_recovery},
\code{stLMM_recovery_chains}, \code{stLMM_prediction},
\code{stLMM_prediction_chains}, \code{stLMM_fitted_chains}, or posterior
sample matrix.}
\item{...}{Additional arguments passed to methods. Supported arguments include
\code{burn}, \code{thin}, \code{metadata}, \code{include_w},
\code{combine_chains}, and, for prediction objects, \code{sample}.}
}
\description{
Collect posterior sample matrices from fitted, recovered, prediction, or
fitted-value objects into a plain data frame with samples in rows and named
sample columns across the top.
}
\details{
For fitted \code{stLMM} objects, columns contain active parameter sample blocks
such as fixed effects, grouped random effects, residual variance parameters,
process variances, and process correlation parameters. Structured latent
process draws are available when they were saved during fitting or added by
\code{recover()}, and are included only when \code{include_w = TRUE}.

For recovery objects, parameter samples are aligned to
\code{object$recover_iter} before saved or recovered latent process draws are appended.
This ensures parameter columns and \code{w_} columns describe the same retained
posterior iterations.

For multi-chain objects, samples are combined by row by default and the
\code{.chain} metadata column preserves chain membership. Set
\code{combine_chains = FALSE} to return one data frame per chain.

Prediction methods accept \code{sample = "mu"}, \code{"y"}, or \code{"all"}.
Prediction columns are named \code{mu_1}, \code{mu_2}, and so on, or
\code{y_1}, \code{y_2}, and so on.
}
\value{
A data frame, or a list of data frames when \code{combine_chains = FALSE}.
When \code{metadata = TRUE}, the first columns are \code{.chain} and
\code{.iteration}. \code{.iteration} is the original MCMC iteration within
chain when that metadata is available.
}
\examples{
\dontrun{
fit <- stLMM(y ~ x, data = dat, n_samples = 1000, chains = 4)
draws <- as_samples(fit, burn = 500, thin = 2)

rec <- recover(fit, sub_sample = list(start = 501, thin = 2))
w_draws <- as_samples(rec, include_w = TRUE)

pred <- predict(rec, newdata = newdat, y_samples = TRUE)
pred_draws <- as_samples(pred, sample = "all")
}
}
\seealso{\code{\link{as_mcmc}}, \code{\link{recover}},
\code{\link{predict.stLMM_recovery}}, \code{\link{fitted.stLMM}}}
