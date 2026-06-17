\name{stLMM-terms}
\alias{stLMM-terms}
\alias{stLMM_terms}
\alias{iid}
\alias{nngp}
\alias{gp}
\alias{ar1}
\alias{car}
\alias{car_time}
\alias{dagar}
\alias{dagar_time}
\title{stLMM formula terms}
\description{
Formula-term constructors recognized by \code{\link{stLMM}}.
}
\details{
The symbols documented here are package-native model terms used inside an
\code{stLMM()} formula. They are not exported as ordinary R functions to call
outside a model formula.
}

\section{IID grouped random effects}{
\code{iid(group)} adds an iid Gaussian random intercept for the observed
levels of \code{group}. It contributes one column per observed group level to
the explicit random-effect design matrix \code{Z}.

\code{x:iid(group)} adds an iid Gaussian group-varying slope for the numeric
covariate \code{x}. It contributes one column per observed group level, with
row entries scaled by \code{x}.

\code{iid(group) + x:iid(group)} creates two separate scalar iid random-effect
terms with separate variance parameters, typically named \code{iid_1} and
\code{iid_2}. This is not a correlated intercept/slope block.

IID random effects are sampled through the explicit \code{Z alpha} random-effect
path. They do not require \code{\link{recover}}. Prediction from a fit object is
available for existing grouping levels. New grouping levels in \code{newdata}
currently error.
}

\section{Structured latent process terms}{
\code{ar1(index)} adds a one-dimensional AR1 latent process over the sorted
unique values of \code{index}. Repeated observations at the same index share a
latent process node. The implemented AR1 support uses adjacency in the sorted
unique fitted support.

\code{gp(coord1, coord2, ..., cov_model = "exp")} adds a dense Gaussian
process over unique coordinate rows. Spatial covariance models include
\code{cov_model = "exp"} with parameter \code{phi} and
\code{cov_model = "matern"} with parameters \code{phi} and \code{nu}.

\code{nngp(coord1, coord2, ..., m = 15, ordering = "coord",
cov_model = "exp")} adds a nearest-neighbor Gaussian process over unique
coordinate rows. Supported \code{ordering} values are \code{"coord"},
\code{"default"}, \code{"maxmin"}, \code{"hilbert"}, \code{"random"}, or a
numeric permutation. The \code{"maxmin"} ordering can be very slow for large
numbers of unique coordinate rows because it uses an exact iterative
max-min construction. The \code{"hilbert"} ordering uses only the first two
coordinate columns to rank nodes; covariance distances still use all coordinate
columns according to the selected \code{cov_model}.
Space-time covariance models include \code{cov_model = "sep_exp"},
\code{"multi_res_sep_exp"}, and \code{"gneiting"}.

\code{car(area_id, graph = g, car_model = "proper")} adds a CAR latent
process over areal nodes defined by \code{g}, where \code{g} is produced by
\code{\link{car_graph}}.
Rows with the same \code{area_id} share a latent process node. Areas in the
graph that do not have observed responses remain in the latent support and can
be recovered when they are represented by rows with missing response values.
Prediction is implemented for existing graph areas. New CAR areas currently
error. For polygon inputs, \code{car_graph()} can optionally add auditable
nearest-neighbor bridge edges for disconnected island components.
\code{car_model = "proper"} uses the package's original degree-scaled proper
CAR precision. \code{car_model = "leroux"} uses the Leroux CAR precision, which
bridges iid area effects and intrinsic-CAR smoothing.

\code{dagar(area_id, graph = g, ordering = "coord")} adds an ordered directed
acyclic graph autoregressive (DAGAR) latent process over the areal nodes in
\code{g}. The same \code{\link{car_graph}} object used by \code{car()} supplies
the undirected adjacency graph; \code{ordering} orients each adjacency edge
from the earlier area to the later area. Supported orderings are
\code{"coord"}, \code{"default"}, \code{"maxmin"}, \code{"hilbert"},
\code{"random"}, or a numeric permutation. Coordinate-based orderings require
an \code{sf}-based graph; for adjacency-matrix graphs, \code{"coord"} and
\code{"default"} use the supplied row order. Prediction is implemented for
existing graph areas. New DAGAR areas currently error.

\code{car_time(area_id, time, graph = g, time_model = "ar1",
car_model = "proper")} adds a separable CAR-by-time latent process over the
full Cartesian product of graph areas and sorted unique fitted time values.
The latent support is stacked
location-major, so all time points for the first area come first. The current
precision is \code{sigma_sq^-1 * (Q_space \%x\% Q_time)}, with one variance
term for the whole space-time process. The default temporal model is
\code{time_model = "ar1"}, which uses adjacency in the sorted fitted time
support. \code{time_model = "exp"} uses
\code{Corr(t_i, t_j) = exp(-lambda * abs(t_i - t_j))}; fitted time values must
be numeric and finite. New prediction rows must use existing areas and existing
fitted time values.

\code{dagar_time(area_id, time, graph = g, ordering = "coord",
time_model = "ar1")} adds a separable DAGAR-by-time latent process over the
full product of ordered DAGAR graph areas and sorted unique fitted time values.
The undirected graph supplied by \code{g} is first oriented by \code{ordering}
as in \code{dagar()}, and the resulting ordered DAGAR spatial precision is
crossed with the temporal precision. The latent support is stacked
ordered-area-major, so all time points for the first ordered area come first.
The user-facing recovered samples are returned in the original graph-area order
crossed with fitted time order. The default temporal model is
\code{time_model = "ar1"}. \code{time_model = "exp"} uses
\code{Corr(t_i, t_j) = exp(-lambda * abs(t_i - t_j))}; fitted time values must
be numeric and finite. New prediction rows must use existing graph areas and
existing fitted time values.

For space-time covariance models, coordinate order matters and the final
coordinate is treated as time. For example,
\code{nngp(lon, lat, time, cov_model = "sep_exp")} is not equivalent to
\code{nngp(time, lon, lat, cov_model = "sep_exp")}.

For new-node space-time NNGP prediction, \code{\link{predict.stLMM_recovery}}
has an \code{st_scale} argument that affects nearest-neighbor ranking only.
It multiplies the final coordinate before ranking fitted-support neighbors; it
does not standardize spatial and temporal coordinates, does not change the
fitted NNGP graph, and does not change covariance calculations after the
neighbor set has been selected. Users should choose \code{st_scale} so the
spatial and temporal extents are balanced in meaningful units for neighbor
selection.

Process terms can be covariate-scaled with the same colon convention used for
spatially varying coefficients, for example \code{x:nngp(lon, lat)} or
\code{x:ar1(day)}.
\code{x:car(area_id, graph = g)} is the CAR spatially varying coefficient
form.
\code{x:car_time(area_id, time, graph = g)} is the CAR-time spatial-temporal
varying coefficient form.
\code{x:dagar_time(area_id, time, graph = g)} is the DAGAR-time
spatial-temporal varying coefficient form.

Structured process terms are collapsed during parameter updates. Gaussian
process fits use \code{\link{recover}} after fitting to draw latent process
samples before using process contributions in \code{\link{fitted.stLMM}} or
\code{\link{predict.stLMM_recovery}}. Polya-Gamma process fits save in-chain
process draws by default as \code{save_process = list(start = 1, thin = 1)} and
\code{\link{recover}} selects from those draws.
}

\section{Covariance functions}{
Dense GP and NNGP terms use the process variance \code{sigma_sq} and the
theta parameters named by \code{\link{get_cor_models}}.

For \code{cov_model = "exp"}, with spatial distance \eqn{h},
\deqn{
C(h) = \sigma^2 \exp(-\phi h).
}

For \code{cov_model = "matern"},
\deqn{
C(h) =
\sigma^2
\frac{(\phi h)^\nu}{2^{\nu - 1}\Gamma(\nu)}
K_\nu(\phi h),
}
where \eqn{K_\nu(\cdot)} is the modified Bessel function of the second kind.

For \code{cov_model = "sep_exp"}, with spatial distance \eqn{h} and temporal
distance \eqn{u},
\deqn{
C(h, u) = \sigma^2 \exp(-\phi h)\exp(-\lambda u).
}

For \code{cov_model = "multi_res_sep_exp"},
\deqn{
C(h, u) =
\sigma^2\left[
\alpha \exp(-\phi_1 h)\exp(-\lambda_1 u) +
(1 - \alpha)\exp(-\phi_2 h)\exp(-\lambda_2 u)
\right].
}
This two-scale separable covariance follows the spatial-temporal forest
inventory model used by May and Finley (2025).

For \code{cov_model = "gneiting"},
\deqn{
C(h, u) =
\sigma^2
(1 + a |u|^{2\alpha})^{-(\delta + d/2)}
\exp\left\{
-\frac{c h^{2\gamma}}{(1 + a |u|^{2\alpha})^{\beta\gamma}}
\right\}.
}
Here \eqn{h} is Euclidean spatial distance, \eqn{u} is temporal distance, and
\eqn{d} is the number of spatial coordinate columns before the final time
coordinate. This is Gneiting (2002, equation 12). Parameter supports are
\eqn{a > 0}, \eqn{c > 0}, \eqn{0 < \alpha \le 1},
\eqn{0 \le \beta \le 1}, \eqn{0 < \gamma \le 1}, and
\eqn{\delta \ge 0}. Free parameters use the strict interior of these supports
where needed by the Metropolis transformation; boundary values such as
\code{beta = fixed(0)}, \code{alpha = fixed(1)}, \code{gamma = fixed(1)}, or
\code{delta = fixed(0)} can be fixed explicitly.

The covariance used by Datta et al. (2016, equation 6.1) is recovered in two
spatial dimensions by setting \code{alpha = fixed(1)},
\code{gamma = fixed(0.5)}, \code{delta = fixed(0)}, and using \code{beta} as
Datta's interaction parameter. A conservative starting point is to estimate
\code{a}, \code{c}, and optionally \code{beta}, while fixing \code{alpha},
\code{gamma}, and \code{delta}.

For \code{ar1(index)}, the fitted support is the sorted unique values of
\code{index}. The correlation between support positions is
\deqn{
\mathrm{Cor}(w_i, w_j) = \phi^{|i-j|}.
}
This is an ordered-support AR1, not a numeric time-gap model.

For \code{car(area, graph = g)}, the proper CAR precision is
\deqn{
Q = \sigma^{-2}(D - \rho G),
}
where \eqn{G} is the binary adjacency matrix stored by \code{car_graph()} and
\eqn{D} is the diagonal degree matrix. The graph is not stored as a
row-standardized matrix; rather, this degree-scaled precision is equivalent to
a CAR conditional specification whose mean uses the neighbor average,
\deqn{
E(w_i \mid w_{-i}) = \frac{\rho}{d_i}\sum_{j \sim i} w_j,
}
with conditional variance \eqn{\sigma^2/d_i}. The current implementation
restricts \eqn{\rho} to \eqn{(0,1)}.

For \code{car(area, graph = g, car_model = "leroux")}, the Leroux CAR
precision is
\deqn{
Q = \sigma^{-2}\{(1-\rho)I + \rho(D - G)\}.
}
This model has the same graph input and parameter names as the proper CAR
model. As \eqn{\rho} approaches zero, the prior approaches iid area effects
with variance \eqn{\sigma^2}; as \eqn{\rho} approaches one, it approaches an
intrinsic-CAR smoothing structure. The implementation keeps \eqn{\rho} strictly
inside \eqn{(0,1)}.

For \code{dagar(area, graph = g)}, let \eqn{\pi} be the ordering of the graph
areas and let
\deqn{
N_\pi(i) = \{j : i \sim j,\ \pi^{-1}(j) < \pi^{-1}(i)\}
}
be the directed parent set for area \eqn{i}. Let \eqn{n_\pi(i)=|N_\pi(i)|}.
The ordered DAGAR conditional specification is
\deqn{
w_i \mid w_{<i,\pi}
\sim
N\left(
\frac{\rho}{1 + \{n_\pi(i)-1\}\rho^2}
\sum_{j \in N_\pi(i)} w_j,\,
\frac{1 + \{n_\pi(i)-1\}\rho^2}{1-\rho^2}
\right),
}
where the second normal argument is precision. Nodes with no directed parents
have conditional mean zero. In stLMM scale convention,
\deqn{
Q = \sigma^{-2} L_\pi' F_\pi L_\pi,
}
where \eqn{L_\pi} is unit lower triangular and \eqn{F_\pi} contains the
conditional precisions. The ordering is part of the model specification, not a
CHOLMOD or display-only ordering. This is the ordered spatial DAGAR model of
Datta, Banerjee, Hodges, and Gao (2019).

For \code{car_time(area, time, graph = g)}, the separable precision is
\deqn{
Q = \sigma^{-2}\{Q_{\mathrm{space}}(\rho) \otimes Q_{\mathrm{time}}\}.
}
\code{car_model} selects whether \code{Q_space} is the proper CAR precision
\eqn{D - \rho G} or the Leroux precision \eqn{(1-\rho)I + \rho(D-G)}.
With \code{time_model = "ar1"}, \code{Q_time} is the ordered-support AR1
precision with parameter \code{phi}. With \code{time_model = "exp"}, the
temporal correlation is
\deqn{
\mathrm{Cor}(w_t, w_{t'}) = \exp(-\lambda |t - t'|),
}
represented internally through its Markov precision on sorted numeric fitted
time values.

For \code{dagar_time(area, time, graph = g)}, the separable precision is
\deqn{
Q = \sigma^{-2}\{Q_{\mathrm{dagar}}(\rho) \otimes Q_{\mathrm{time}}\}.
}
\code{Q_dagar} is the scale-free ordered DAGAR precision defined above, after
orienting graph edges by \code{ordering}. \code{Q_time} is the same AR1 or
exponential-time Markov precision used by \code{car_time()}. The log
determinant is evaluated by the Kronecker identity,
\deqn{
\log|Q| = -ST\log(\sigma^2) +
T\log|Q_{\mathrm{dagar}}(\rho)| +
S\log|Q_{\mathrm{time}}|,
}
where \eqn{S} is the number of graph areas and \eqn{T} is the number of fitted
time values.
}

\section{Controls and priors}{
Term-specific covariance controls use generated names such as \code{iid_1},
\code{ar1_1}, \code{gp_1}, \code{nngp_1}, \code{car_1}, \code{dagar_1},
\code{car_time_1}, and \code{dagar_time_1}. The names are assigned by term
type in formula order.

IID variance priors are inverse-gamma prior objects, for example:
\preformatted{
priors = list(iid_1 = list(sigma_sq = ig(shape, scale)))
}

Process terms require a process variance prior plus finite-support priors for
their theta parameters, for example:
\preformatted{
priors = list(
  nngp_1 = list(
    sigma_sq = ig(shape, scale),
    phi = uniform(lower, upper)
  )
)
}

CAR terms use a process variance prior and CAR \code{rho} bounds:
\preformatted{
priors = list(
  car_1 = list(
    sigma_sq = ig(shape, scale),
    rho = uniform(lower, upper)
  )
)
}

DAGAR terms use the same process variance and \code{rho} prior structure:
\preformatted{
priors = list(
  dagar_1 = list(
    sigma_sq = ig(shape, scale),
    rho = uniform(lower, upper)
  )
)
}

CAR-time terms use one process variance, spatial CAR \code{rho}, and a temporal
parameter. For \code{time_model = "ar1"}, the temporal parameter is AR1
\code{phi}:
\preformatted{
priors = list(
  car_time_1 = list(
    sigma_sq = ig(shape, scale),
    rho = uniform(lower, upper),
    phi = uniform(lower, upper)
  )
)
}
For \code{time_model = "exp"}, replace \code{phi} with positive temporal decay
\code{lambda}.

DAGAR-time terms use the same process variance, ordered-DAGAR \code{rho}, and
temporal parameter structure:
\preformatted{
priors = list(
  dagar_time_1 = list(
    sigma_sq = ig(shape, scale),
    rho = uniform(lower, upper),
    phi = uniform(lower, upper)
  )
)
}
For \code{time_model = "exp"}, replace \code{phi} with positive temporal decay
\code{lambda}.

Use \code{\link{get_cor_models}} to inspect supported \code{cov_model} values
and their correlation parameter names.
}
\references{
Gneiting, T. (2002). Nonseparable, stationary covariance functions for
space-time data. \emph{Journal of the American Statistical Association},
97(458), 590--600. doi:10.1198/016214502760047113.

May, P. B. and Finley, A. O. (2025). Spatial-temporal prediction of forest
attributes using latent Gaussian models and inventory data.
\emph{Spatial Statistics}, 69, 100917. doi:10.1016/j.spasta.2025.100917.
}
\examples{
adj <- matrix(
  c(
    0, 1, 0,
    1, 0, 1,
    0, 1, 0
  ),
  nrow = 3,
  byrow = TRUE,
  dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
)
g <- car_graph(adj)

\dontrun{
fit <- stLMM(
  y ~ x + car(area, graph = g),
  data = dat,
  priors = list(
    resid = list(tau_sq = ig(2, 1)),
    car_1 = list(sigma_sq = ig(2, 1), rho = uniform(0.05, 0.95))
  ),
  n_samples = 1000
)
rec <- recover(fit)
fitted(rec)
}

\dontrun{
fit_st <- stLMM(
  y ~ x + car_time(area, time, graph = g),
  data = dat_st,
  priors = list(
    resid = list(tau_sq = ig(2, 1)),
    car_time_1 = list(
      sigma_sq = ig(2, 1),
      rho = uniform(0.05, 0.95),
      phi = uniform(-0.95, 0.95)
    )
  ),
  n_samples = 1000
)
}
}
\seealso{
\code{\link{stLMM}}, \code{\link{recover}}, \code{\link{fitted.stLMM}},
\code{\link{predict.stLMM_recovery}}, \code{\link{get_cor_models}},
\code{\link{stLMM_priors}}
}
