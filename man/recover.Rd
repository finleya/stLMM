\name{recover}
\alias{recover}
\alias{recover.stLMM}
\alias{recover.stLMM_chains}
\title{Recover latent process draws}
\usage{
recover(object, ...)

\method{recover}{stLMM}(object, n_omp_threads = NULL, verbose = FALSE,
  sub_sample = list(start = 1L, thin = 1L), ...)

\method{recover}{stLMM_chains}(object, n_omp_threads = NULL, verbose = FALSE,
  sub_sample = list(start = 1L, thin = 1L), ...)
}
\arguments{
\item{object}{An object returned by \code{stLMM()}. Multi-chain fits are
handled chain by chain.}
\item{n_omp_threads}{Optional positive integer number of OpenMP threads for
threaded recovery setup kernels. The default \code{NULL} uses the thread count
stored in the fitted object. This mainly affects NNGP covariance-coefficient
updates needed while reconstructing latent process draws; sparse CHOLMOD
factorization and random-number generation are not OpenMP-threaded by
\code{recover()}.}
\item{verbose}{Logical; if \code{TRUE}, print phase-level recovery progress
messages. Progress is reported from R before and after major recovery phases
rather than from individual OpenMP worker threads.}
\item{sub_sample}{A list with optional integer entries \code{start} and
\code{thin} selecting posterior draws to recover or select from saved process
draws.}
\item{...}{Currently unused.}
}
\description{
Recover posterior draws of latent process terms from a fitted collapsed
\code{stLMM} model.
}
\details{
For Gaussian process models, the collapsed sampler integrates latent process
values out during fitting. This function reconstructs the posterior
conditional distribution of the latent process for selected posterior draws and
samples from it using the stored backend, posterior parameter samples, and
CHOLMOD sparse factorization.

By default, recovery uses the OpenMP thread count stored in the fitted object.
Use \code{n_omp_threads} to override that value for a recovery call without
refitting. The override is most relevant for NNGP terms because their
covariance-coefficient updates are threaded. Other major recovery work,
including sparse factorization and latent Gaussian random draws, is not
parallelized by this argument.

For Polya-Gamma process models, process draws are saved during \code{stLMM()}
by default as \code{save_process = list(start = 1, thin = 1)}. In that case,
\code{recover()} only selects and labels the saved
in-chain process draws. If a Polya-Gamma process model was fit with
\code{save_process = FALSE}, \code{recover()} errors and the model must be
refit with \code{save_process = TRUE} or
\code{save_process = list(start = ..., thin = ...)}.

Recovery is only needed for models with structured process terms such as
\code{ar1()}, \code{gp()}, \code{nngp()}, \code{car()}, \code{dagar()}, or
\code{car_time()}, or \code{dagar_time()}.
Explicit iid grouped random effects are sampled during fitting and do not need
this recovery step. When the fitted data include rows with missing numeric
responses, recovered process samples are on the full latent support, including
latent nodes associated only with those missing-response rows.

For Polya-Gamma process models, the sampler uses warm-running internal
Polya-Gamma latent-process draws to update the augmentation variables. Those
same in-chain process draws are the process samples used by \code{recover()}.
Post-fit Polya-Gamma reconstruction is intentionally not supported.
}
\value{
The fitted object with recovered latent process draws added. User-facing
\code{w_samples} are returned by process term; NNGP and DAGAR terms are in
original support order. The internal ordered copy used by fitted values and
prediction is stored in \code{w_samples_ordered}. The object also includes
\code{w_samples_stacked} and \code{recover_iter}. The returned object has class
\code{stLMM_recovery} and can be passed to \code{fitted()} or \code{predict()}.
For multi-chain input, the returned object has class
\code{stLMM_recovery_chains} and contains \code{chains}, a list of ordinary
\code{stLMM_recovery} objects. The \code{sub_sample} rule is applied within
each chain.
}
\examples{
\dontrun{
fit <- stLMM(...)
rec <- recover(fit, sub_sample = list(start = 501, thin = 5))
fitted(rec)
}
}
\seealso{\code{\link{stLMM}}, \code{\link{fitted.stLMM}},
\code{\link{predict.stLMM_recovery}}, \code{\link{summary.stLMM_recovery}}}
