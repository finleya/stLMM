\name{as_mcmc}
\alias{as_mcmc}
\alias{as_mcmc.stLMM}
\alias{as_mcmc.stLMM_chains}
\alias{as_mcmc.stLMM_recovery}
\alias{as_mcmc.stLMM_recovery_chains}
\title{Convert stLMM samples to coda objects}
\usage{
as_mcmc(object, ...)
}
\arguments{
\item{object}{An \code{stLMM}, \code{stLMM_chains}, \code{stLMM_recovery}, or
\code{stLMM_recovery_chains} object.}
\item{...}{Additional arguments passed to methods. Supported arguments include
\code{include_w}, logical, \code{burn}, a nonnegative integer, and
\code{thin}, a positive integer.}
}
\description{
Convert retained posterior samples to \pkg{coda} objects for chain diagnostics
and plotting.
}
\details{
Single-chain fits are returned as \code{coda::mcmc}. Multi-chain fits are
returned as \code{coda::mcmc.list}. By default, saved or recovered latent
process samples \code{w} are omitted because they can be high-dimensional and
are often not needed for covariance-parameter diagnostics.

For ordinary fitted objects, \code{burn} and \code{thin} select retained MCMC
rows by row number. For recovered objects, \code{include_w = TRUE} aligns
parameter samples to \code{recover_iter} before appending saved or recovered latent
process draws. In that case, \code{burn} is interpreted on the original
posterior iteration scale and \code{thin} is applied to the recovered draws
that remain after burn-in. Multi-chain recovery objects apply the same rule
within each chain and return a \code{coda::mcmc.list}.
}
\value{
A \code{coda::mcmc} or \code{coda::mcmc.list} object containing active posterior
sample blocks such as fixed effects, grouped random effects, residual
variance parameters, process variances, and process correlation parameters.
}
\examples{
\dontrun{
fit <- stLMM(y ~ x, data = dat, n_samples = 1000, chains = 4)
m <- as_mcmc(fit, burn = 500)
coda::gelman.diag(m)
coda::effectiveSize(m)

rec <- recover(fit, sub_sample = list(start = 501, thin = 5))
rec_m <- as_mcmc(rec, include_w = TRUE)
}
}
\seealso{\code{\link{as_samples}}, \code{\link{recover}},
\code{\link{stLMM}}, \code{\link{summary.stLMM}}}
