\name{resid}
\alias{resid}
\title{Residual variance formula term}
\usage{
resid(model = c("tau_sq", "constant", "fixed", "group", "scaled"), group,
  variance, vhat, n = NULL, prior = c("default", "ig", "shannon"),
  method, shape = 4, center = c("mean", "mode"), shrinkage = 10,
  kappa_log_prior = c(mean = 0, sd = 1), tau0_sq_log_prior = NULL,
  starting = NULL, tuning = 0.1)
}
\arguments{
\item{model}{Residual variance model. \code{"tau_sq"} and \code{"constant"}
estimate one global Gaussian residual variance. \code{"fixed"} uses fixed
row-specific variances. \code{"group"} estimates one residual variance per
group. \code{"scaled"} estimates a low-dimensional scaling model for supplied
row-specific variances.}
\item{group}{Grouping variable for sampled group-specific residual variances.}
\item{variance}{Numeric vector or expression evaluated in the model data giving
one residual variance per row. This is used by \code{model = "fixed"},
\code{model = "group"} with direct-estimate-informed priors, and
\code{model = "scaled"}. Variances for observed response rows must be finite
and positive. Missing response rows may have missing variances because they do
not enter the likelihood.}
\item{vhat}{Alias for \code{variance}, useful when the supplied variances are
direct-estimate variances. For grouped direct-estimate-informed residual
variances, this is supplied per row and must be constant within group for
observed rows.}
\item{n}{Optional effective sample size. For \code{prior = "shannon"}, this is
group-level. For \code{model = "scaled"}, this is row-level and controls
shrinkage toward a common residual scale.}
\item{prior}{Prior construction for \code{model = "group"} when
\code{variance} is supplied. \code{"ig"} uses a constant-shape inverse-gamma
prior centered on the supplied variance. \code{"shannon"} uses the
effective-sample-size construction described below.}
\item{method}{Alias for \code{prior}; retained only to make the argument name
read naturally when comparing constant-shape and Shannon-style constructions.}
\item{shape}{Inverse-gamma shape for \code{prior = "ig"}.}
\item{center}{Whether the constant-shape inverse-gamma prior mean or mode is
centered on \code{variance}.}
\item{shrinkage}{Positive constant controlling how quickly the sample-size
aware scaled model trusts \code{variance}. The weight is
\code{n / (n + shrinkage)}.}
\item{kappa_log_prior}{Mean and standard deviation of the normal prior on
\code{log(kappa)}.}
\item{tau0_sq_log_prior}{Optional mean and standard deviation of the normal
prior on \code{log(tau0_sq)} for the sample-size aware scaled model. If
omitted, the prior is centered on the median observed \code{variance}.}
\item{starting}{Optional starting values for residual variance parameters used
by the grouped direct-estimate-informed and scaled models. For
\code{model = "group"} without \code{variance}, supply starting values through
\code{starting = list(resid = list(tau_sq = value))} in \code{\link{stLMM}}.}
\item{tuning}{Initial log-scale Metropolis proposal scale for sampled residual
variance parameters.}
}
\description{
Defines the Gaussian residual variance model as part of an \code{\link{stLMM}}
formula. Omitting the term is equivalent to writing \code{resid()} or
\code{resid(model = "tau_sq")}: one global residual variance \code{tau_sq} is
estimated.
}
\details{
\code{resid()} is an \code{stLMM} formula term, not the usual fitted-model
\code{resid()} or \code{residuals()} extractor. Only one \code{resid()} term
is allowed, and it cannot be used in interactions. The formula term selects
the residual variance model; corresponding priors, starting values, and tuning
values are supplied through the \code{resid} block in \code{priors},
\code{starting}, and \code{tuning}. Gaussian residual variance models are not
used with Polya-Gamma likelihoods such as \code{family = "binomial"} or
\code{family = "negative_binomial"}.

\code{resid(model = "tau_sq")} estimates one scalar residual variance. Use
\code{priors = list(resid = list(tau_sq = ig(shape, scale)))} or
\code{priors = list(resid = list(tau_sq = half_t(df, scale)))} to place a
prior on this variance model.

\code{resid(model = "fixed", variance = vhat)} uses observed-row precision
\code{obs_precision_i = 1 / vhat_i}. When this model is used, \code{tau_sq} is
not sampled. Prediction of new response draws with \code{y_samples = TRUE}
requires residual variance information for the prediction rows.

\code{resid(model = "group", group = g)} estimates one residual variance per
group. Priors are supplied through \code{priors$resid}; use
\code{priors = list(resid = list(tau_sq = ig(shape, scale)))} to apply one
prior to all groups, or supply a named list with one prior per group. Starting
values and tuning values are supplied through
\code{starting$resid$tau_sq} and \code{tuning$resid$tau_sq}; scalars are
recycled over groups and named vectors or lists may be used when
group-specific values are needed. Values marked with \code{\link{fixed}} are
held fixed, assigned zero tuning, and do not require residual group priors.

\code{resid(model = "group", group = g, variance = vhat, prior = "ig")}
also estimates one residual variance per group, but constructs inverse-gamma
priors from externally supplied group-level variance information. For
\code{prior = "ig"}, \code{shape} is the inverse-gamma shape and the scale is
chosen so either the prior mean or prior mode equals \code{vhat}. With
\code{prior = "shannon"}, the prior is
\code{IG(n / 2, (n - 1) * vhat / 2)} and requires \code{n > 1}; larger
\code{n} gives a more concentrated prior around the supplied direct-estimate
variance.

\code{resid(model = "scaled", variance = vhat)} estimates a low-dimensional
multiplier for direct-estimate variances. Without \code{n}, the model is
\code{tau_i^2 = kappa * vhat_i}. With \code{n}, the model uses
\code{log(tau_i^2) = log(kappa) + w_i log(vhat_i) + (1 - w_i) log(tau0_sq)}
where \code{w_i = n_i / (n_i + shrinkage)}. The positive parameters
\code{kappa} and \code{tau0_sq} are sampled on the log scale, so zero or
negative residual variances are not possible.
}
\examples{
\dontrun{
fit <- stLMM(
  yhat ~ x + car(area, graph = g) + resid(model = "fixed", variance = vhat),
  data = dat,
  n_samples = 1000,
  priors = list(car_1 = list(sigma_sq = ig(2, 1), rho = uniform(0.05, 0.95)))
)

fit2 <- stLMM(
  yhat ~ x + car(area, graph = g) + resid(model = "group", group = area),
  data = dat,
  n_samples = 1000,
  priors = list(
    resid = list(tau_sq = half_t(df = 3, scale = 1)),
    car_1 = list(sigma_sq = ig(2, 1), rho = uniform(0.05, 0.95))
  )
)

fit3 <- stLMM(
  yhat ~ x + car(area, graph = g) +
    resid(model = "group", group = area, variance = vhat, prior = "ig", shape = 6),
  data = dat,
  n_samples = 1000,
  priors = list(car_1 = list(sigma_sq = ig(2, 1), rho = uniform(0.05, 0.95)))
)

fit4 <- stLMM(
  yhat ~ x + car(area, graph = g) + resid(model = "scaled", variance = vhat, n = n_eff),
  data = dat,
  n_samples = 1000,
  priors = list(car_1 = list(sigma_sq = ig(2, 1), rho = uniform(0.05, 0.95)))
)
}
}
\seealso{\code{\link{stLMM}}, \code{\link{predict.stLMM_recovery}}}
