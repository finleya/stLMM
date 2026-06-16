\name{stLMM-priors}
\alias{stLMM-priors}
\alias{stLMM_priors}
\alias{flat}
\alias{normal}
\alias{ig}
\alias{uniform}
\alias{log_normal}
\alias{gamma_dist}
\alias{half_normal}
\alias{half_t}
\alias{beta_dist}
\alias{fixed}
\title{Prior constructors for stLMM}
\usage{
flat()

normal(mean = 0, sd = 1)

ig(shape, scale)

uniform(lower, upper)

log_normal(meanlog, sdlog, support = NULL)

gamma_dist(shape, rate, support = NULL)

half_normal(scale)

half_t(df, scale)

beta_dist(shape1, shape2)

fixed(value)
}
\arguments{
\item{mean, sd}{Mean and standard deviation for independent normal
fixed-effect priors created by \code{normal()}. Scalars are recycled across
fixed effects; vectors must match the number of fixed effects.}
\item{shape, scale}{For \code{ig()}, positive inverse-gamma shape and scale
parameters. For \code{half_normal()}, positive standard-deviation scale.}
\item{lower, upper}{Finite support limits for \code{uniform()}.}
\item{meanlog, sdlog}{Mean and standard deviation on the log scale for
\code{log_normal()}; \code{sdlog} must be positive.}
\item{support}{Optional finite support \code{c(lower, upper)}. This is required
when \code{log_normal()} or \code{gamma_dist()} is used for a covariance theta
parameter.}
\item{rate}{Positive gamma rate parameter.}
\item{df}{Positive degrees of freedom for \code{half_t()}.}
\item{shape1, shape2}{Positive beta shape parameters.}
\item{value}{Finite numeric scalar value at which a covariance or variance
parameter should be fixed.}
}
\description{
Constructors for prior and fixed-parameter objects used by \code{\link{stLMM}}.
}
\details{
Raw numeric vectors such as \code{c(2, 1)} are not accepted as priors. Use an
explicit constructor so the prior family and parameter support are clear.

\code{fixed(value)} marks a covariance or variance parameter as fixed at
\code{value}. It is used in \code{starting}, not \code{priors}. Internally this
sets the starting value and sets the corresponding Metropolis tuning value to
zero. Fixed parameters are omitted from Metropolis proposal blocks and do not
require priors. Their values still define the likelihood and latent-process
precision. If a fixed parameter also has an explicitly positive tuning value,
\code{stLMM()} errors.

\code{flat()} and \code{normal(mean, sd)} are fixed-effect priors and are
currently accepted only as \code{priors = list(beta = ...)}. A flat beta prior
is the improper prior used by earlier versions of the sampler. A normal beta
prior is an independent Gaussian prior on the fixed-effect coefficients. The
Gaussian likelihood defaults to \code{flat()}; Polya-Gamma likelihoods default
to \code{normal(mean = 0, sd = 10)}.

The inverse-gamma constructor uses the shape/scale parameterization
\deqn{
p(x \mid a, b) =
\frac{b^a}{\Gamma(a)} x^{-(a+1)} \exp\{-b/x\}, \quad x > 0,
}
where \code{shape = a} and \code{scale = b}. Thus the second argument to
\code{ig()} is a scale parameter for the inverse-gamma distribution, not a rate.

Other positive variance priors use the following parameterizations. The
\code{log_normal(meanlog, sdlog)} prior is the usual log-normal density with
\code{meanlog} and \code{sdlog} on \eqn{\log(x)}. The
\code{gamma_dist(shape, rate)} prior uses the gamma shape/rate density
\eqn{p(x) \propto x^{a-1}\exp\{-r x\}}. The \code{half_normal(scale)} and
\code{half_t(df, scale)} priors are specified on the standard deviation
\eqn{s = \sqrt{x}} and transformed internally to the variance scale by the
Jacobian \eqn{1/(2\sqrt{x})}. Up to normalizing constants, this gives
\eqn{p(x) \propto x^{-1/2}\exp\{-x/(2A^2)\}} for the half-normal and
\eqn{p(x) \propto x^{-1/2}\{1 + x/(\nu A^2)\}^{-(\nu+1)/2}} for the half-t.

\code{beta_dist(shape1, shape2)} uses the standard beta density
\eqn{p(x) \propto x^{a-1}(1-x)^{b-1}} on \eqn{0 < x < 1}.

Variance parameters, including \code{tau_sq} and process \code{sigma_sq}, are
strictly positive. Sampled residual and process variances may use
\code{ig()}, \code{log_normal()}, \code{gamma_dist()}, \code{half_normal()}, or
\code{half_t()}. The half-normal and half-t priors are defined on the standard
deviation scale and transformed internally to the variance scale. IID grouped
random-effect variances currently use conjugate Gibbs updates and must use
\code{ig()}.

When \code{resid(model = "group", group = group)} is used, the grouped residual
variance priors are supplied through \code{priors$resid}. A single prior can
be recycled over groups, for example
\code{priors = list(resid = list(tau_sq = half_t(df = 3, scale = 1)))}.
Alternatively, \code{priors$resid} may be a named list with one prior for
each observed residual group.

Covariance theta parameters use a bounded-logit proposal transform. Therefore,
all theta priors must have finite support. Positive theta parameters such as
\code{phi}, \code{lambda}, \code{nu}, \code{a}, \code{c}, and Gneiting
\code{delta} when free may
use \code{uniform()}, \code{log_normal(..., support = c(lower, upper))}, or
\code{gamma_dist(..., support = c(lower, upper))}. Unit-interval theta
parameters such as mixture \code{alpha} and Gneiting \code{alpha},
\code{beta}, and \code{gamma} may use \code{uniform()} or \code{beta_dist()}.
Bounded AR1 parameters currently use \code{uniform()}.

For exponential covariance \eqn{\exp(-\phi h)}, the distance where correlation
equals 0.05 is approximately \eqn{-\log(0.05) / \phi}. This is often a useful
way to choose finite support for \code{phi} after coordinates have been scaled
or projected into meaningful units. For binomial and other weak-information
settings, spatial decay parameters can be weakly identified and may be
sensitive to the chosen support, neighbor size, and coordinate scale; inspect
trace plots and posterior mass near prior boundaries.
}
\examples{
\dontrun{
fit <- stLMM(
  y ~ x + nngp(lon, lat, cov_model = "matern"),
  data = dat,
  priors = list(
    beta = normal(mean = 0, sd = 10),
    resid = list(tau_sq = ig(2, 0.1)),
    nngp_1 = list(
      sigma_sq = half_t(df = 4, scale = 1),
      phi = log_normal(log(5), 0.75, support = c(0.5, 20)),
      nu = uniform(0.1, 2)
    )
  ),
  n_samples = 1000
)

fit_fixed <- stLMM(
  y ~ x + nngp(lon, lat),
  data = dat,
  starting = list(
    resid = list(tau_sq = fixed(0.5)),
    nngp_1 = list(sigma_sq = fixed(1), phi = fixed(3))
  ),
  n_samples = 1000
)
}
}
\seealso{\code{\link{stLMM}}, \code{\link{stLMM_terms}},
\code{\link{get_cor_models}}}
