\name{get_cor_models}
\alias{get_cor_models}
\title{List supported covariance models}
\usage{
get_cor_models()
}
\value{
A named list describing available covariance models. Each entry contains
the model's parameter names, parameter type codes, theta domains, and distance
mode code.
}
\description{
Returns metadata for covariance models supported by the current sampler.
This is a lightweight discovery helper for valid \code{cov_model} values and
the theta parameters expected by each model.
}
\details{
The \code{domains} entry is used to validate theta priors. Positive theta
parameters require finite positive support. Unit theta parameters have support
within \eqn{[0,1]}.

The returned metadata applies to covariance-model labels used by \code{gp()}
and \code{nngp()} terms. CAR-family terms have fixed parameter names documented
in \code{\link{stLMM_terms}}.
}
\examples{
names(get_cor_models())
get_cor_models()[["exp"]]
get_cor_models()[["matern"]]
}
