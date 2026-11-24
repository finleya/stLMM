\name{recover}
\alias{recover}
\alias{recover.stLMM}
\alias{recover.stLMM_chains}
\alias{recover.stLMM_recovery_chains}
\title{Recover latent effect draws}
\usage{
recover(object, ...)

\method{recover}{stLMM}(object, n_omp_threads = NULL, verbose = FALSE,
  sub_sample = list(start = 1L, thin = 1L), ...)

\method{recover}{stLMM_chains}(object, n_omp_threads = NULL, verbose = FALSE,
  sub_sample = list(start = 1L, thin = 1L), ...)

\method{recover}{stLMM_recovery_chains}(object, n_omp_threads = NULL, verbose = FALSE,
  sub_sample = list(start = 1L, thin = 1L), ...)
}
\arguments{
\item{object}{An object returned by \code{stLMM()} or \code{recover()}. Multi-chain objects are
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
Recover posterior draws of latent effects from a fitted collapsed
\code{stLMM} model.
}
\details{
For Gaussian process models, the collapsed sampler integrates latent process
values out during parameter updates. This function reconstructs the posterior
conditional distribution of the latent process for selected posterior draws and
samples from it using the stored backend, posterior parameter samples, and
CHOLMOD sparse factorization.

By default, recovery uses the OpenMP thread count stored in the fitted object.
Use \code{n_omp_threads} to override that value for a recovery call without
refitting. The override is most relevant for NNGP terms because their
covariance-coefficient updates are threaded. Other major recovery work,
including sparse factorization and latent Gaussian random draws, is not
parallelized by this argument.

For augmented non-Gaussian process models, process draws are saved during \code{stLMM()}
by default as \code{save_process = list(start = 1, thin = 1)}. In that case,
\code{recover()} only selects and labels the saved
in-chain process draws. If an augmented non-Gaussian process model was fit with
\code{save_process = FALSE}, \code{recover()} errors and the model must be
refit with \code{save_process = TRUE} or
\code{save_process = list(start = ..., thin = ...)}.

For Gaussian models with \code{collapse_iid = TRUE}, recovery draws IID and
structured coefficients jointly, including IID-only models. Fitting integrates
out both coefficient blocks and retains their covariance parameters. Recovery
uses the combined design and sparse precision to preserve posterior dependence
between IID and structured effects. Recovered \code{alpha_samples} and structured
samples have one row per selected iteration, identified by \code{recover_iter}.
The parameter matrices retain all original iterations. Prediction, fitted values,
and combined coefficient summaries require recovery first.
Calling \code{recover()} again on a Gaussian recovery object, including a
multi-chain object, draws fresh latent effects for the requested original
iterations. For augmented models, repeated recovery selects from the process
draws still present in the object.

With separate IID sampling, coefficients are sampled during fitting and recovery
draws structured effects conditional on the stored \code{alpha_samples}.
Previously saved fits containing full in-chain IID draws retain this behavior.
For augmented models, recovery selects the paired in-chain draws and does not
resample IID coefficients. When the fitted
data include rows with missing numeric
responses, recovered process samples are on the full latent support, including
latent nodes associated only with those missing-response rows.

For augmented non-Gaussian process models, the sampler uses warm-running internal
latent-process draws to update the augmentation variables. Those
same in-chain process draws are the process samples used by \code{recover()}.
Post-fit augmented non-Gaussian process reconstruction is intentionally not supported.
}
\value{
The fitted object with recovered latent effect draws added. Joint Gaussian
recovery adds \code{alpha_samples} for the selected iterations. User-facing
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
set.seed(1)
dat <- data.frame(
  y = rnorm(8),
  x = rnorm(8),
  time = seq_len(8)
)

fit <- stLMM(
  y ~ x + ar1(time),
  data = dat,
  n_samples = 8,
  priors = list(
    resid = list(tau_sq = ig(2, 1)),
    ar1_1 = list(sigma_sq = ig(2, 1), phi = uniform(-0.8, 0.8))
  ),
  warmup = FALSE,
  verbose = FALSE
)
rec <- recover(fit, sub_sample = list(start = 5, thin = 1))
stats::fitted(rec)
}
\seealso{\code{\link{stLMM}}, \code{\link[stats]{fitted}},
\code{\link{predict.stLMM_recovery}}, \code{\link{summary.stLMM_recovery}}}
