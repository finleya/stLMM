\name{stLMM}
\alias{stLMM}
\title{Fit a Bayesian mixed model with structured latent processes}
\usage{
stLMM(formula, data = parent.frame(), starting = list(), tuning = NULL,
  priors = NULL, family = "gaussian", trials = NULL,
  size = NULL, n_samples, n_omp_threads = 1, nngp_search = c("fast", "brute"),
  verbose = TRUE, n_report = 100, warmup = TRUE, metropolis = list(),
  cholmod_control = list(), save_process = NULL, chains = 1L,
  chain_control = list(), describe_terms = FALSE, collapse_iid = FALSE, ...)
}
\arguments{
\item{formula}{A model formula with fixed effects, optional \code{offset()}
terms, optional iid grouped
random-effect terms such as \code{iid(group)}, \code{x:iid(group)}, or
\code{iid(group, ~ 1 + x)} or \code{iid(plot, ~ 1 + x, by = stand)}, and
optional structured process terms such as \code{ar1()}, \code{gp()},
\code{nngp()}, \code{car()}, \code{dagar()}, \code{car_time()}, or
\code{dagar_time()}.}
\item{data}{A data frame or environment containing variables in \code{formula}.}
\item{starting}{Optional named list of starting values. The public convention
uses underscores. Default residual variance is supplied as
\code{starting = list(resid = list(tau_sq = value))}. Process starts use term-specific
blocks, for example \code{starting = list(nngp_1 = c(sigma_sq = 1, phi = 1))}.
IID random-effect variance starts use the same term-block convention, for
example \code{starting = list(iid_1 = list(sigma_sq = 1))}. Unstructured coefficient-formula
iid terms use \code{starting = list(iid_1 = list(covariance = Sigma))}, with
a symmetric positive definite matrix in coefficient order. If supplied,
dimnames must exactly match that order; a single matrix is used for every
chain. Omitted covariance starts use the inverse-Wishart prior mode, with
positive diagonal congruence perturbations for multiple chains. Process covariance and
residual variance parameters can be fixed with \code{\link{fixed}}, for
example \code{starting = list(resid = list(tau_sq = fixed(0.5)))}.}
\item{tuning}{Optional named list of initial Metropolis proposal scales for
structured-process covariance, residual variance, and estimated negative-binomial
size parameters. Scalar IID variances with non-IG priors also use Metropolis
updates, as do all IID variances and covariances with \code{collapse_iid = TRUE};
see the IID sampling strategies below. A zero Metropolis tuning value fixes the
corresponding parameter. Parameters supplied with \code{\link{fixed}} are
assigned zero tuning automatically.}
\item{priors}{Named list of prior objects created with constructors such as
\code{\link{flat}}, \code{\link{normal}}, \code{\link{ig}}, \code{\link{iw}},
\code{\link{uniform}}, \code{\link{log_normal}},
\code{\link{gamma_dist}}, \code{\link{half_normal}}, \code{\link{half_t}}, and
\code{\link{beta_dist}}. Fixed effects use \code{priors = list(beta = flat())}
or \code{priors = list(beta = normal(mean, sd))}. Default residual variance
uses \code{resid$tau_sq}; for example,
\code{priors = list(resid = list(tau_sq = ig(shape, scale)))} puts an inverse-gamma prior
directly on \code{tau_sq}, while
\code{priors = list(resid = list(tau_sq = half_t(df, scale)))} defines a half-\eqn{t} prior
on the residual standard deviation \eqn{\sqrt{\tau^2}} and transforms it
internally to the variance scale. Process terms require
term-specific blocks, for example
\code{priors = list(nngp_1 = list(sigma_sq = ig(shape, scale), phi = uniform(lower, upper)))}.
IID random-effect variance priors use the same term-block convention, for
example \code{priors = list(iid_1 = list(sigma_sq = ig(shape, scale)))} or
\code{priors = list(iid_1 = list(sigma_sq = half_t(df, scale)))}. Scalar IID
variances accept \code{ig()}, \code{log_normal()}, \code{gamma_dist()},
\code{half_normal()}, and \code{half_t()}; the last two specify a prior on
the standard deviation. There is no default IID variance prior.
Unstructured coefficient-formula iid terms require
\code{priors = list(iid_1 = list(covariance = iw(df, scale)))}. These
covariances use Gibbs updates when IID effects are separate and Metropolis
updates when collapsed. Use zero tuning to hold an IID covariance fixed under
collapsing; \code{fixed()} is not accepted for IID covariance matrices.
Independent coefficient-formula terms use \code{cov_model = "independent"}
and require a scalar variance prior or a coefficient-named list
of variance priors under their formula term name. A shared specification sets
the same prior independently on each variance. Starting variances may be a
positive scalar or a coefficient-named numeric vector or list under
\code{starting$iid_1$sigma_sq}. Independent variances use the selected IID sampling strategy.
Priors
are not required for covariance or residual variance parameters supplied with
\code{\link{fixed}}.}
\item{family}{Model likelihood family. Use \code{"gaussian"} for Gaussian
linear mixed models or \code{"binomial"} for binomial logistic models using
Polya-Gamma Gibbs updates. Use \code{"probit"} for binary probit models using
Albert-Chib latent-normal augmentation. Use \code{"negative_binomial"} for
negative-binomial log-mean count models using Polya-Gamma Gibbs updates. The
default is \code{"gaussian"}. For compatibility with familiar R model syntax,
\code{gaussian()}, \code{binomial()}, \code{binomial(link = "logit")}, and
\code{binomial(link = "probit")} are also accepted. Other binomial links are
not currently supported.}
\item{trials}{Optional binomial trial counts for \code{family = "binomial"}.
The default \code{NULL} uses one trial per row, i.e. Bernoulli data. A single
character string names a column in \code{data}; otherwise supply a positive
integer-valued vector with length matching the model frame. The response should
contain successes and must be between zero and \code{trials}.}
\item{size}{Positive fixed negative-binomial size parameter for
\code{family = "negative_binomial"}; \code{fixed(value)} is also accepted.
To estimate size, leave this argument \code{NULL} and supply
\code{priors = list(size = gamma_dist(shape, rate))}. A fixed size and a size
prior cannot both be supplied. The fitted linear predictor is
\eqn{\eta_i = \log(\mu_i)}, where \eqn{\mu_i} is the mean count.}
\item{n_samples}{Number of MCMC samples.}
\item{n_omp_threads}{Number of OpenMP threads for supported kernels.}
\item{nngp_search}{Nearest-neighbor search used when constructing NNGP graphs.
\code{"fast"} uses the exact history-restricted k-d tree implementation.
\code{"brute"} checks every eligible historical node and is mainly useful for
testing and debugging.}
\item{verbose}{Logical; print setup, sampler-start, and progress reports.}
\item{n_report}{Progress report interval for retained sampling and warmup.
For warmup, reports are rounded to completed warmup batches. Use \code{0} to
suppress progress reports while keeping other \code{verbose} setup messages.}
\item{warmup}{Logical or list controlling discarded Metropolis proposal-scale
calibration before retained sampling. \code{TRUE} uses default calibration;
\code{FALSE} disables it. A list may contain \code{enabled},
\code{batch_length}, \code{min_batches}, \code{max_batches}, \code{target},
and \code{near_zero}. \code{min_batches} forces at least that many warmup
batches before early stopping is allowed; \code{max_batches} remains the
ceiling.
Defaults are \code{enabled = TRUE}, \code{batch_length = 25},
\code{min_batches = 0}, \code{max_batches = 20}, and
\code{near_zero = 0.02}. When \code{target} is omitted, the default
acceptance interval is \code{c(0.15, 0.45)} for block updates and
\code{c(0.30, 0.60)} for scalar updates.}
\item{metropolis}{Optional sampler-control list for Metropolis updates of
covariance, residual variance, and estimated negative-binomial size parameters.
The default \code{list(blocking = "joint")} keeps active covariance and
residual variance parameters in one adaptive covariance block. Estimated size
always has its own scalar block. Built-in alternatives are \code{"by_term"},
\code{"residual_process"}, \code{"variance_theta"},
\code{"process_variance"}, and \code{"scalar"}. A single string such as
\code{metropolis = "by_term"} is also accepted. A list may also contain
\code{target_accept} for retained-sampling adaptation and \code{batch_length}
for retained adaptation batches. Defaults are \code{target_accept = 0.234} for
block updates, \code{0.44} for scalar updates, and
\code{batch_length = 25}.}
\item{cholmod_control}{Optional sparse factorization control list. The
\code{ordering} element controls the CHOLMOD fill-reducing ordering used for
latent sparse matrix factorizations. Supported values are \code{"auto"} (the
default CHOLMOD strategy), \code{"best"} (try CHOLMOD's full built-in method
set), \code{"natural"}, \code{"amd"}, \code{"metis"}, \code{"nesdis"}, and
\code{"colamd"}. A single string such as \code{cholmod_control = "amd"} is
also accepted. The \code{postorder} element is logical and defaults to
\code{TRUE}. METIS-based choices require the CHOLMOD library exposed by
\pkg{Matrix} to have been built with METIS/partitioning support.}
\item{save_process}{Controls whether structured latent process draws are saved
during MCMC. Use \code{NULL} for the default: process draws are saved
automatically for augmented non-Gaussian process models as
\code{list(start = 1, thin = 1)} and not saved otherwise. Use
\code{TRUE} as shorthand for \code{list(start = 1, thin = 1)}, use
\code{FALSE} to turn saving off, or use a list with \code{start} and
\code{thin}, for example
\code{save_process = list(start = 501, thin = 5)}. For augmented non-Gaussian process
models, saved in-chain process draws are required for \code{recover()},
\code{fitted()}, and process prediction. Enabling \code{save_process} is only
supported for augmented non-Gaussian models with structured process terms.}
\item{chains}{Positive integer number of MCMC chains to run. The default
\code{1L} returns an ordinary \code{stLMM} object. Values greater than one run
sequential chains and return an \code{stLMM_chains} object.}
\item{chain_control}{Optional list controlling multi-chain initialization.
\code{seed} initializes reproducible chain-specific random number streams when
supplied. \code{dispersion} controls log-scale jitter for omitted positive
starting values.}
\item{describe_terms}{Logical; print model term diagnostics before sampling.}
\item{collapse_iid}{IID sampling strategy, applied to all IID terms. The default
\code{FALSE} keeps IID coefficients separate from the marginalized structured
effects. IG variances and IW covariances use conjugate Gibbs updates; scalar
variances with other priors use conditional Metropolis updates. \code{TRUE} integrates IID
and structured coefficients together for fixed-effect and covariance updates.
Gaussian coefficients are jointly drawn afterward by \code{recover()}; augmented
likelihoods require joint coefficient draws within each iteration.
This experimental strategy uses Metropolis
updates for IID covariances and requires at least one IID term. \code{NULL}
selects collapsing for IID-only models and separate IID updates when structured
process terms are present. Grouped residual variances do not count as structured
process terms. See Details for tuning and computational trade-offs.}
\item{...}{Currently unused.}
}
\value{
With \code{chains = 1}, an object of class \code{stLMM}. With
\code{chains > 1}, an object of class \code{stLMM_chains} with component
\code{chains}, a list of ordinary \code{stLMM} objects, one per chain.

For each ordinary fit, posterior samples are available both as top-level
entries such as \code{beta_samples}, \code{tau_sq_samples},
\code{sigma_sq_samples}, \code{theta_samples}, \code{alpha_samples},
\code{iid_sigma_sq_samples} (scalar and independent coefficient-formula iid variances), and
\code{iid_covariance_samples} (unstructured coefficient-formula terms), and
through the active-component list \code{samples}. Empty sample blocks for absent
model components are \code{NULL}. Covariance Metropolis diagnostics are
returned as \code{covariance_acceptance} and \code{adaptive_metropolis}. The
default is one global covariance block; optional blocking and scalar-update
diagnostics are stored in \code{adaptive_metropolis$blocks}. The object also
contains \code{term_description} and the stored \code{backend} used by recovery
and prediction. \code{backend$collapse_iid_requested} records the supplied
strategy (including \code{NULL}), and \code{backend$collapse_iid} records the
resolved logical choice. Both are retained within every chain.
For unstructured coefficient-formula iid terms, \code{iid_covariance_samples} and
\code{samples$iid_covariance} contain the lower triangle of each covariance,
including the diagonal, column by column, in iid term order. Columns are named
\code{iid_1_cov[x,(Intercept)]}, for example. Coefficients in
\code{alpha_samples} are ordered by sampler block, then group, then coefficient.
The old scalar iid output names and ordering are preserved.
Independent coefficient formulas expand into scalar blocks in coefficient order;
their output names include \code{[coefficient="x"]}, for example.
Each scalar or unstructured block with \code{by} expands by represented parent, ordered
by the parent's factor levels (or the levels of \code{factor(parent)}).
Within each block, child levels and coefficients retain their usual order.
Output names identify the parent, for example
\code{iid_1[stand="north"]_sigma_sq} or
\code{iid_2[stand="north"]_cov[x,(Intercept)]}.
Each block has an entry in \code{term_description$random_effects$terms},
with \code{formula_term}, \code{by}, \code{by_level}, and \code{by_levels}
identifying its original term and parent support. Input priors and starting
values remain keyed by the original formula term name.

\code{timing$backend} contains \code{proc.time()} CPU and elapsed wall-time
differences around the entire compiled call. \code{timing$phases} is a named
numeric vector of monotonic elapsed seconds for compiled \code{initialization},
proposal-calibration \code{warmup}, and post-warmup \code{sampling}. The sampling
timer surrounds exactly the \code{n_samples} MCMC iterations, including repeated
covariance updates, factorizations, sample storage, and any in-chain latent-effect
recovery. It excludes R-side preparation, graph construction, compiled
initialization, calibration warmup, and final result assembly. Later user-side
burn-in or thinning does not change the recorded duration. A skipped warmup may
record a very small duration for the warmup function's bookkeeping. These are
wall times, not CPU times or effective-sample-size measurements.

Multi-chain fits provide \code{timing$backend_by_chain},
\code{timing$backend_total}, \code{timing$phases_by_chain}, and
\code{timing$phases_total}. Totals sum the durations across chains. The backend
timer includes finalization, which is outside the three reported phases.
In earlier builds without phase timing, the backend fields were named \code{sampler},
\code{sampler_by_chain}, and \code{sampler_total}; those names did not denote
sampling-only time and are no longer used. Earlier saved fits must be refitted
to obtain phase timings. The DESCRIPTION field
\code{Config/stLMM/Timing: phases-v1} identifies this capability; an earlier
0.0.5 installation may lack it.

Latent process draws are saved during fitting only for augmented non-Gaussian likelihoods
with structured process terms. The default in that case is
\code{save_process = list(start = 1, thin = 1)}, because those in-chain process
draws are already needed by the likelihood augmentation and post-fit
augmented process reconstruction is not used. For Gaussian process models,
structured draws are not retained during fitting and
\code{recover()} is the intended post-fit step for drawing latent process
samples.
Multi-chain fits can be converted to a
\code{coda::mcmc.list} with \code{\link{as_mcmc}} and summarized with
\code{summary()}, which includes chain diagnostics when possible.
}
\description{
Fits Bayesian linear mixed models with fixed effects, optional explicit iid
grouped random effects, and optional structured latent process terms. Structured
process terms are integrated out during parameter updates. For Gaussian process
models, their latent draws are recovered afterward with \code{\link{recover}}.
For augmented non-Gaussian process models, the in-chain process draws already used by the
augmentation step are saved and reused by fitted-value, recovery, and
prediction methods.

For \code{family = "binomial"}, the sampler uses Polya-Gamma latent variables.
Conditional on the Polya-Gamma weights, the model is fit as a Gaussian working
model with diagonal observation precision equal to the current weights.
Polya-Gamma random variates are generated through BayesLogit's hybrid sampler.
Binomial trials default to one per row and can be supplied through
\code{trials}. The Gaussian \code{resid()} residual variance term and
\code{tau_sq} controls are not
used for binomial fits. For structured binomial process models, the sampler
saves in-chain process draws by default; \code{\link{recover}} selects from
those saved draws.

For \code{family = "probit"}, the response must be binary 0/1. The sampler
uses Albert-Chib latent-normal variables with unit latent residual variance.
Conditional on these latent responses, the model is fit as a Gaussian working
model with unit observation precision. The Gaussian \code{resid()} residual
variance term and \code{tau_sq} controls are not used for probit fits. For
structured probit process models, the sampler saves in-chain process draws by
default; \code{\link{recover}} selects from those saved draws.

For \code{family = "negative_binomial"}, the response must be non-negative
integer counts. Supply a positive fixed \code{size}, or estimate one shared
size parameter using an explicit \code{priors$size = gamma_dist(shape, rate)}.
There is no automatic size prior. The Gamma prior is on size itself and uses
shape and rate; truncated supports are not supported for this parameter.
For estimated size, \code{starting$size} accepts a positive scalar or one value
per chain. The default center is the prior mean (shape/rate), with the usual
log-scale dispersion for omitted multi-chain starts. \code{tuning$size} is the
positive initial log-size proposal standard deviation, default 0.25. These
controls are rejected for fixed size; zero tuning does not select fixed size.

Estimated size has its own scalar Metropolis block, regardless of covariance
blocking. It shares the scalar adaptation implementation and the \code{warmup}
and \code{metropolis} controls. Its default acceptance target is 0.44 and its
default warmup target interval is \code{c(0.30, 0.60)}; explicitly supplied
\code{metropolis$target_accept} or \code{warmup$target} overrides these defaults
for size as well as the covariance blocks. Adaptation continues during retained
sampling. Size-only models also run calibration warmup.

Size draws are returned as \code{size_samples} and \code{samples$size};
\code{adaptive_metropolis$size} records its acceptance, proposal histories, and
warmup diagnostics separately from covariance-block diagnostics. Summaries and
sample extraction include the parameter \code{size}.

Size estimation is allowed with IID and structured random effects. The fit
issues one warning because size and random-effect variances can be weakly
identified. Check chain diagnostics and sensitivity to their priors. Within-group
replication and dependence patterns help distinguish these sources of variation.
For example, \code{gamma_dist(1, 0.01)} has mean 100 and a central 95 percent
interval approximately 2.5--369; it is not a noninformative prior.

Extremely small sizes can cause numerical underflow in the underlying
Polya-Gamma generator, particularly for zero counts. The sampler stops with
the failing shape and tilt rather than replacing the draw or bounding size.

The sampler uses the log-mean parameterization \eqn{\eta_i=\log(\mu_i)} and an
internal logit tilt \eqn{\psi_i=\eta_i-\log(\code{size})}. Conditional on the
Polya-Gamma weights, the model again becomes a Gaussian working model with
diagonal observation precision equal to the current weights. Residual
\code{resid()} terms and \code{tau_sq} controls are not used for this
likelihood. Larger \code{size} values make the negative-binomial distribution
closer to Poisson, but can slow the sampler because the Polya-Gamma shape is
the observed count plus \code{size}.
}
\details{
Predictor variables must be complete. Missing numeric responses are allowed:
only observed response rows enter the likelihood, while the full row support is
retained for fitted values, recovery, and prediction metadata. The number of
missing responses is reported by \code{print()} and \code{summary()}.
Formula offsets are supported through standard R \code{offset()} terms. The
offset is a known additive component of the linear predictor. For
negative-binomial count models, \code{offset(log_exposure)} gives
\eqn{\log(\mu_i) = \log(\mathrm{exposure}_i) + \cdots}; for binomial models it
is a known log-odds contribution. Fitted values and predictions include the
offset on both link and response scales. When predicting with \code{newdata},
variables used inside \code{offset()} must be present in \code{newdata}.

The intended workflow is:
\preformatted{
fit <- stLMM(...)
summary(fit)
rec <- recover(fit)
stats::fitted(rec)
predict(rec, newdata = ...)
}
For models without structured process terms, \code{fitted()} and
\code{predict()} can be called directly on the fit object, except for Gaussian
models with collapsed IID effects, which require \code{recover()} first.

Prior lists use explicit prior constructors. Raw numeric vectors such as
\code{c(2, 1)} are not accepted as priors. Common patterns are:
\preformatted{
## residual variance
priors = list(resid = list(tau_sq = ig(shape, scale)))
priors = list(resid = list(tau_sq = half_t(df = 3, scale = 0.5)))

## fixed effects
priors = list(beta = flat())
priors = list(beta = normal(mean = 0, sd = 10))

## known additive offset
y ~ x + offset(log_exposure) + car(area, graph = g)

## fixed direct-estimate sampling variances
y ~ x + resid(model = "fixed", variance = vhat)

## sampled group-specific residual variances
y ~ x + resid(model = "group", group = area)
priors = list(resid = list(tau_sq = half_t(df = 3, scale = 1)))

## direct-estimate-informed group-specific residual variances
y ~ x + resid(model = "group", group = area, variance = vhat,
              prior = "ig", shape = 6)

## scaled direct-estimate sampling variances
y ~ x + resid(model = "scaled", variance = vhat, n = n_eff)

## iid random-effect variance for iid(group)
priors = list(iid_1 = list(sigma_sq = ig(shape, scale)))

## fixed covariance parameters
starting = list(
  resid = list(tau_sq = fixed(0.5)),
  nngp_1 = list(sigma_sq = fixed(1), phi = fixed(3))
)

## process covariance controls
priors = list(nngp_1 = list(sigma_sq = ig(shape, scale),
                            phi = uniform(lower, upper)))

## CAR controls
priors = list(car_1 = list(sigma_sq = ig(shape, scale),
                           rho = uniform(lower, upper)))

## DAGAR controls
priors = list(dagar_1 = list(sigma_sq = ig(shape, scale),
                             rho = uniform(lower, upper)))

## separable CAR-time controls
priors = list(car_time_1 = list(sigma_sq = ig(shape, scale),
                                rho = uniform(lower, upper),
                                phi = uniform(lower, upper)))

## separable CAR-time controls with time_model = "exp"
priors = list(car_time_1 = list(sigma_sq = ig(shape, scale),
                                rho = uniform(lower, upper),
                                lambda = uniform(lower, upper)))

## separable DAGAR-time controls
priors = list(dagar_time_1 = list(sigma_sq = ig(shape, scale),
                                  rho = uniform(lower, upper),
                                  phi = uniform(lower, upper)))
}
\itemize{
\item \strong{Starting values, tuning, and fixed parameters.}
Omitting \code{starting} uses data-based or conservative defaults. Omitting
\code{tuning} uses default initial Metropolis proposal scales. Setting a tuning
value to zero fixes that structured-process covariance or residual variance parameter. Parameters
marked with \code{fixed()} use their supplied values and are omitted from
Metropolis proposal blocks.

IID parameters updated by Metropolis accept zero tuning to hold the supplied
starting variance or covariance fixed. If its starting
value is omitted, the default is held fixed across chains. IID priors remain
required, including when tuning is zero; \code{fixed()} is not an IID control.
With separate IID effects, IG variances and IW covariances use Gibbs updates
and cannot be fixed through tuning. In an independent coefficient formula
mixing IG and other priors, tuning applies only to the non-IG variances.

\item \strong{IID sampling strategies (experimental joint collapsing).}
The default \code{collapse_iid = FALSE} keeps IID coefficients separate and
the structured sparse factorization smaller. IG variances and IW covariances
retain conjugate updates. Other scalar variance priors use Metropolis updates
conditional on the IID coefficients. With
\code{collapse_iid = TRUE}, a single sparse system includes both IID and
structured coefficients. Integrating them together can reduce dependence in the
fixed-effect updates, but IID covariance updates then require Metropolis steps.
Groups connecting many structured nodes, or crossed IID groups, can increase
Cholesky fill and factorization cost. Neither strategy is uniformly faster or
better mixing. Compare posterior agreement, R-hat, effective sample sizes per
second, and memory for the quantities relevant to the analysis. Joint collapsing
is experimental; it changes the sampling strategy, not the model or prior.

For scalar or independent variances updated by Metropolis, use
\code{tuning = list(iid_1 = list(sigma_sq = 0.2))}. Independent coefficient
formulas also accept coefficient-named numeric vectors or lists, just as their
starting values do. For an unstructured covariance, use
\code{tuning = list(iid_1 = list(covariance = 0.1))}; this scalar supplies the
initial scale for every log-Cholesky coordinate. These are also the default
scales. A \code{by} term reuses its tuning specification for each parent's
separate matrix or variance. Tuning entries must be finite and nonnegative.

IID variance Metropolis updates use log coordinates. Collapsed unstructured covariances use a
lower Cholesky factor with log diagonal and unrestricted off-diagonal entries.
The prior densities include the corresponding Jacobians in the Metropolis
target. All existing \code{metropolis} blocking options apply. \code{"by_term"}
uses a separate block for each independent IID variance or covariance, including
each parent-specific block. \code{"variance_theta"} and
\code{"process_variance"} place all IID coordinates in their variance block;
\code{"residual_process"} places them with process covariance parameters.

Gaussian fits with \code{collapse_iid = TRUE} defer IID and structured coefficient
draws to \code{recover()}. This includes IID-only fits: call \code{recover()}
before prediction, fitted values, or combined coefficient summaries. Joint
recovery returns \code{alpha_samples} and structured effects on the selected
\code{recover_iter} grid; the parameter samples retain all fitting iterations.
Separate IID sampling retains \code{alpha_samples} during fitting, and Gaussian
structured recovery conditions on those coefficients.
Augmented non-Gaussian models use paired in-chain IID and process draws on the
\code{save_process} grid. Formula syntax, combined summaries, and explicit
\code{allow_new_levels} prediction rules are the same for both strategies.

\item \strong{Priors and parameter transforms.}
Process terms require priors for parameters that are sampled. The Gaussian
likelihood defaults to a flat fixed effect prior, while augmented non-Gaussian likelihoods
default to \code{normal(mean = 0, sd = 10)} for fixed effects. The
\code{flat()} and \code{normal()} constructors are accepted only for
\code{priors$beta}. Variance parameters are positive. The
\code{half_normal()} and \code{half_t()} variance-prior constructors are
specified on the corresponding standard deviation scale and transformed
internally to the variance scale. Covariance theta parameters use a
bounded-logit proposal transform and therefore require finite support. For
positive theta parameters, \code{log_normal()} and \code{gamma_dist()} must
include \code{support = c(lower, upper)}. See
\code{\link{stLMM_priors}} and \code{\link{stLMM_terms}} for details.

\item \strong{Multiple chains.}
For \code{chains > 1}, omitted starting values are dispersed automatically
from the defaults. With \code{chain_control$seed}, chain-specific random
number streams are generated from the supplied seed, so dispersed starts are
reproducible across different Metropolis blocking choices. If a starting value
is supplied by the user, it controls that parameter: scalar values are recycled
to all chains, and vectors with length \code{chains} provide chain-specific
starts. Other lengths are rejected. For process terms, chain-specific starts
should be supplied in list blocks, for example
\code{starting = list(nngp_1 = list(sigma_sq = c(1, 2, 3),
phi = c(0.2, 0.4, 0.6)))}.

For an independent iid coefficient formula, names under \code{sigma_sq} identify
coefficients, not chains. A coefficient-named list may contain a scalar or a
chain-specific vector for each coefficient. A shared scalar applies to every
coefficient and chain. Unstructured iid terms use one starting covariance matrix
for every chain.

\item \strong{Warmup and retained-sampling adaptation.}
By default, the Metropolis sampler runs a short discarded warmup phase that
calibrates proposal scales before retained samples are collected. Warmup draws
are not returned in posterior sample matrices; warmup diagnostics are stored in
\code{fit$adaptive_metropolis$warmup}. Retained-sampling adaptation is
controlled separately through \code{metropolis$target_accept} and
\code{metropolis$batch_length}.

\item \strong{Metropolis blocking.}
Metropolis blocking changes only the proposal partition used for covariance and
residual variance parameters; it does not change the model. Blocking can
improve mixing for multi-term or residual-variance models, but each block may
require a separate collapsed likelihood evaluation and sparse factorization.
The default \code{"joint"} block is usually the fastest choice. Treat the
other blocking schemes as diagnostic or model-specific alternatives: compare
trace plots, R-hat, effective sample sizes, and runtime before changing the
default for routine analyses. \code{metropolis = "scalar"} updates one active
covariance or residual variance parameter at a time using scalar random-walk
adaptation. It can be useful when block proposals mix poorly, but it may be
substantially slower for models with expensive collapsed likelihood
evaluations.

\item \strong{Sparse factorization controls.}
The \code{cholmod_control} argument is intended for advanced diagnostics and
performance benchmarking. The default \code{ordering = "auto"} uses CHOLMOD's
default ordering strategy. Alternative orderings can change sparse Cholesky
fill-in, floating-point work, and runtime without changing the statistical
model. Use \code{describe_terms = TRUE} or inspect \code{fit$term_description}
to compare the selected CHOLMOD ordering, fill ratio, \code{lnz}, and flops.
The \code{"best"} option can make symbolic analysis slower because CHOLMOD
tries more orderings before factorization.

\item \strong{Fixed residual variances.}
When \code{resid(model = "fixed", variance = vhat)} is used, \code{tau_sq} is
not sampled and should not be supplied in \code{starting}, \code{tuning}, or
\code{priors}.
}
See \code{\link{stLMM_terms}} for the package-native formula terms recognized
by \code{stLMM()}.
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
\seealso{\code{\link{stLMM_terms}}, \code{\link{stLMM_priors}},
\code{\link{recover}}, \code{\link[stats]{fitted}},
\code{\link{predict.stLMM_recovery}}, \code{\link{summary.stLMM}}}
