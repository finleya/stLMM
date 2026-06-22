\name{predict.stLMM_recovery}
\alias{predict.stLMM}
\alias{predict.stLMM_chains}
\alias{predict.stLMM_recovery}
\alias{predict.stLMM_recovery_chains}
\title{Posterior prediction from recovered stLMM fits}
\usage{
\method{predict}{stLMM}(object, newdata = NULL, y_samples = FALSE,
  n_omp_threads = NULL, verbose = FALSE, sub_sample = list(start = 1L,
  thin = 1L), scale = c("response", "link"), ...)

\method{predict}{stLMM_chains}(object, newdata = NULL, y_samples = FALSE,
  n_omp_threads = NULL, verbose = FALSE, sub_sample = list(start = 1L,
  thin = 1L), scale = c("response", "link"), ...)

\method{predict}{stLMM_recovery}(object, newdata = NULL, y_samples = FALSE,
  joint = FALSE, joint_method = c("full", "vecchia"), pred_m = NULL,
  pred_ordering = "maxmin", st_scale = 1,
  return_w_samples = TRUE, n_omp_threads = NULL, verbose = FALSE,
  sub_sample = list(start = 1L, thin = 1L), scale = c("response", "link"),
  ...)

\method{predict}{stLMM_recovery_chains}(object, newdata = NULL,
  y_samples = FALSE, joint = FALSE, joint_method = c("full", "vecchia"),
  pred_m = NULL, pred_ordering = "maxmin", st_scale = 1,
  return_w_samples = TRUE, n_omp_threads = NULL, verbose = FALSE,
  sub_sample = list(start = 1L, thin = 1L), scale = c("response", "link"),
  ...)
}
\arguments{
\item{object}{An \code{stLMM} fit without structured process terms, a
Polya-Gamma process fit with saved in-chain process draws, a recovered object
returned by \code{recover()}, or the corresponding multi-chain object.}
\item{newdata}{Optional data frame. Grouped random-effect levels must already
exist in the fitted model support. For recovered process models, process
coordinates may be existing support locations or new AR1, GP, or NNGP nodes.
CAR and CAR-time prediction rows must use existing graph areas; CAR-time rows
must also use fitted time values.}
\item{y_samples}{Logical; if \code{TRUE}, also simulate observation-level
posterior predictive samples. Gaussian fits use the fitted residual variance.
Binomial fits draw counts from the fitted probabilities and prediction trial
counts.}
\item{joint}{Logical; for new GP or NNGP nodes, whether to draw unique new
nodes with cross-location dependence. The default \code{FALSE} simulates
unique new nodes independently conditional on the fitted support or
fitted-support neighbors and is recommended when only marginal prediction
summaries are needed. For dense \code{gp()} terms, \code{joint = TRUE} uses
the exact dense GP conditional covariance among unique new nodes.}
\item{joint_method}{Method used for NNGP terms when \code{joint = TRUE}.
\code{"full"} keeps the historical behavior: fitted-support neighbors are
used for each new node and the new-node residual covariance is drawn from a
dense matrix. \code{"vecchia"} orders the unique new prediction nodes and
simulates them sequentially from a Vecchia/NNGP-style history containing all
fitted support nodes plus previously simulated prediction nodes. Dense
\code{gp()} terms ignore \code{joint_method} and use exact joint Gaussian
conditioning when \code{joint = TRUE}.}
\item{pred_m}{Optional positive integer neighbor count for
\code{joint_method = "vecchia"}. The default \code{NULL} uses the fitted
NNGP term's \code{m}. The first prediction nodes always have the fitted support
available; if fewer than \code{pred_m} history nodes are available, all
available history nodes are used.}
\item{pred_ordering}{Ordering for unique new prediction nodes under
\code{joint_method = "vecchia"}. May be one of the NNGP ordering options
\code{"coord"}, \code{"default"}, \code{"maxmin"}, \code{"hilbert"},
\code{"random"}, or a numeric permutation of the unique new prediction nodes.}
\item{st_scale}{Positive scalar, or named positive vector/list by process term.
For space-time NNGP new-node prediction, the final coordinate is multiplied by
\code{st_scale} before nearest-neighbor ranking. This is not automatic
standardization: spatial coordinates are left on their original scale and the
time coordinate is left on its original scale except for this single
multiplier during neighbor search. The covariance calculation still uses the
original coordinates.}
\item{return_w_samples}{Logical; for recovered process prediction, whether to
return raw latent process samples for each process term in \code{w_samples}.
The default \code{TRUE} is useful for diagnostics and process maps. Set to
\code{FALSE} for large prediction problems when memory use matters. For
covariate-scaled process terms, returned values are the raw latent process
draws, not the covariate-multiplied row-level contribution.}
\item{n_omp_threads}{Optional positive integer number of OpenMP threads for
threaded prediction kernels. The default \code{NULL} uses the thread count
stored in the fitted or recovered object. This affects threaded NNGP
prediction neighbor searches and compiled NNGP prediction simulation kernels;
R-level dense joint GP and dense joint NNGP prediction are not OpenMP-threaded.}
\item{verbose}{Logical; if \code{TRUE}, print phase-level prediction progress
messages. Progress is reported from R before and after major prediction phases
rather than from individual OpenMP worker threads.}
\item{sub_sample}{A list with optional integer entries \code{start} and
\code{thin}. For \code{predict.stLMM()} without process terms, these select
posterior draw indices. For Polya-Gamma process fits with saved process draws,
and for \code{predict.stLMM_recovery()}, these select saved or recovered draw
indices through \code{object$recover_iter}.}
\item{scale}{For \code{family = "binomial"}, return fitted probabilities on
the \code{"response"} scale or the linear predictor on the \code{"link"} scale.
For \code{family = "negative_binomial"}, return fitted mean counts on the
\code{"response"} scale or the log mean on the \code{"link"} scale. For
Gaussian models these scales are identical.}
\item{...}{Currently unused.}
}
\description{
Returns posterior samples for the mean response at fitted rows or at
\code{newdata} rows. For binomial models, the default mean response is the
posterior fitted probability. For negative-binomial models, the default mean
response is the posterior fitted mean count.
}
\details{
Prediction always represents the model as fit. Prediction-time dropping of
model components is not supported. New grouped random-effect levels are not
supported.

By default, prediction uses the OpenMP thread count stored in the fitted or
recovered object. Use \code{n_omp_threads} to override that value for a
prediction call without refitting or re-running \code{recover()}. Because the
long-running prediction kernels can be OpenMP threaded, \code{verbose = TRUE}
reports progress at the R phase level rather than printing from worker threads.

With \code{newdata = NULL}, prediction returns fitted-row posterior mean
samples for the original data rows. With \code{newdata}, fixed-effect and iid
random-effect design matrices are rebuilt from the fitted formula metadata.
If the fitted formula contains \code{offset()}, the same offset expression is
evaluated in \code{newdata} and added to the prediction linear predictor.
New AR1, GP, and NNGP process nodes are supported. Repeated prediction rows
sharing a new process node reuse the same simulated latent draw within a
posterior sample.

When \code{return_w_samples = TRUE}, recovered-process prediction objects also
contain \code{w_samples}, a named list with one matrix per process term. Each
matrix has posterior draws in rows and prediction rows in columns, aligned with
\code{mu_samples}. These are raw latent process values on the prediction rows.
For spatially varying coefficient terms, \code{w_samples} does not include the
covariate scaling used when adding the process contribution to
\code{mu_samples}.

For \code{family = "binomial"}, prediction first constructs posterior samples
of the linear predictor and then applies the inverse-logit transform when
\code{scale = "response"}. If \code{y_samples = TRUE}, prediction simulates
binomial counts. With \code{newdata = NULL}, the fitted trial counts are used.
With \code{newdata}, a column named \code{trials} is used when present;
otherwise prediction defaults to one trial per row.

For \code{family = "negative_binomial"}, prediction first constructs posterior
samples of the log mean and then exponentiates them when
\code{scale = "response"}. If \code{y_samples = TRUE}, prediction simulates
negative-binomial counts using the fitted fixed \code{size}.

For dense \code{gp()} terms, \code{joint = FALSE} draws each unique new node
from its marginal conditional distribution given the fitted support. With
\code{joint = TRUE}, unique new GP nodes are drawn jointly from the exact
conditional Gaussian distribution using the dense covariance among prediction
nodes and fitted nodes. This is appropriate for small prediction sets and is
dense in the number of unique new GP nodes.

For NNGP terms, \code{joint = FALSE} uses the independent new-node prediction
approach based on nearest fitted-support neighbors, following the conditional
NNGP prediction machinery described by Finley et al. (2019). Each unique new
node is simulated independently conditional on its fitted-support neighbors.
With \code{joint = TRUE, joint_method = "full"}, nearest neighbors are still
chosen only from the fitted support, but unique new NNGP nodes are drawn
jointly from a dense residual covariance. This can require substantial memory
for large prediction sets.

With \code{joint = TRUE, joint_method = "vecchia"}, unique new NNGP nodes are
ordered by \code{pred_ordering}. Neighbors for each ordered prediction node are
chosen from the combined history made of all fitted support nodes and earlier
prediction nodes in that ordering. The prediction draw is then generated
sequentially, so simulated earlier prediction nodes can enter the conditional
mean of later prediction nodes. This avoids a dense new-node covariance matrix
and is usually the more scalable joint prediction method. This option follows
the observed-prediction ordering idea and latent-first full-conditioning
construction described by Katzfuss et al. (2020, Sects. 3.1, 4.2, and 4.3.4),
adapted to the package workflow in which fitted latent process draws are
saved or recovered first and prediction nodes are then simulated sequentially
conditional on those draws. The method is ordering dependent, as in other
Vecchia/NNGP approximations.

The independent, dense-joint, and Vecchia-joint approaches should have similar
posterior means and medians for each prediction location when the neighbor
sets are adequate. They differ mainly in the cross-location dependence of the
simulated prediction draws.

For space-time NNGP prediction, \code{st_scale} controls the relative influence
of the final coordinate when choosing nearest neighbors. The
neighbor-search distance is equivalent to Euclidean distance after replacing
\code{time} by \code{st_scale * time}. Larger values make temporal separation
matter more and tend to choose neighbors closer in time; smaller values make
spatial separation dominate. Choose \code{st_scale} based on the units of the
coordinates and the desired space-time neighbor extent. For example, if spatial
coordinates are kilometers and time is years, \code{st_scale = 50} treats one
year of time separation like about 50 kilometers for neighbor ranking. This
choice does not change the covariance parameters or covariance evaluations
once the neighbors have been selected.
}
\references{
Finley, A. O., Datta, A., Cook, B. D., Morton, D. C., Andersen, H. E.,
and Banerjee, S. (2019). Efficient Algorithms for Bayesian Nearest Neighbor
Gaussian Processes. \emph{Journal of Computational and Graphical Statistics},
28(2), 401--414. \doi{10.1080/10618600.2018.1537924}.

Katzfuss, M., Guinness, J., Gong, W., and Zilber, D. (2020). Vecchia
Approximations of Gaussian-Process Predictions. \emph{Journal of Agricultural,
Biological and Environmental Statistics}, 25(3), 383--414.
\doi{10.1007/s13253-020-00401-7}.
}
\value{
An object of class \code{stLMM_prediction}, a list with \code{mu_samples},
optional \code{y_samples}, optional \code{w_samples}, \code{draw_index}, and
metadata. Binomial prediction objects also include \code{scale}. For
multi-chain input, the returned object has class
\code{stLMM_prediction_chains} and contains \code{chains}, a list of ordinary
\code{stLMM_prediction} objects. The \code{sub_sample} rule is applied within
each chain.
}
\examples{
\dontrun{
rec <- recover(fit)
pred_fit <- predict(rec)
pred_new <- predict(rec, newdata = newdat, y_samples = TRUE)
summary(pred_new)
}
}
\seealso{\code{\link{recover}}, \code{\link{fitted.stLMM}},
\code{\link{summary.stLMM_prediction}}}
