#include "stLMM_internal.h"

#include <algorithm>
#include <chrono>

#include <BayesLogit.h>

/* Main collapsed MCMC driver and sampler reporting.
 *
 * R/stLMM.R assembles the backend contract; stLMM_state.cpp validates and
 * parses it; stLMM_precision.cpp maintains sparse M. This file runs the Gibbs
 * and adaptive Metropolis updates, records posterior samples, and returns the
 * fit object pieces consumed by the R methods.
 */

static inline int sampler_space_time_node(int space, int time, int nTime)
{
  return space * nTime + time;
}

void scatter_Z_col_to_dense(SamplerState *s, int col, double *out)
{
  int k;

  for(k = 0; k < s->n; k++)
    out[k] = 0.0;

  if(s->q <= 0 || s->Zp == NULL)
    return;

  for(k = s->Zp[col]; k < s->Zp[col + 1]; k++)
    out[s->Zi[k]] = s->Zx[k];
}

void sparse_Z_mult(SamplerState *s, const double *alpha, double *out)
{
  int i, j, k;

  for(i = 0; i < s->n; i++)
    out[i] = 0.0;

  if(s->q <= 0 || s->Zp == NULL)
    return;

  for(j = 0; j < s->q; j++){
    if(alpha[j] == 0.0)
      continue;
    for(k = s->Zp[j]; k < s->Zp[j + 1]; k++)
      out[s->Zi[k]] += s->Zx[k] * alpha[j];
  }
}

void sparse_Zt_mult(SamplerState *s, const double *v, double *out)
{
  int j, k;

  for(j = 0; j < s->q; j++)
    out[j] = 0.0;

  if(s->q <= 0 || s->Zp == NULL)
    return;

  for(j = 0; j < s->q; j++){
    for(k = s->Zp[j]; k < s->Zp[j + 1]; k++)
      out[j] += s->Zx[k] * v[s->Zi[k]];
  }
}

double sparse_Z_col_dot_dense(SamplerState *s, int col, const double *v)
{
  int k;
  double out;

  out = 0.0;
  for(k = s->Zp[col]; k < s->Zp[col + 1]; k++)
    out += s->Zx[k] * v[s->Zi[k]];
  return out;
}

void update_mean_without_alpha(SamplerState *s, const double *beta, double *out)
{
  int inc;
  double one, minusOne;

  inc      = 1;
  one      = 1.0;
  minusOne = -1.0;

  F77_CALL(dcopy)(&s->n, s->y, &inc, out, &inc);
  if(s->p > 0){
    F77_CALL(dgemv)("N", &s->n, &s->p,
                    &minusOne, s->X, &s->n, beta, &inc,
                    &one, out, &inc FCONE);
  }
}

/* Cache the small coefficient precision for repeated group blocks. */
void refresh_iid_covariance_precision(SamplerState *s, int b)
{
  int k, j, col, info;
  char lower;
  double *L, *P;

  k = s->reCovDim[b];
  lower = 'L';
  L = s->reCovChol[b];
  P = s->reCovPrecision[b];

  for(j = 0; j < k * k; j++){
    L[j] = s->reCovariance[b][j];
    P[j] = 0.0;
  }
  for(j = 0; j < k; j++)
    P[j + k * j] = 1.0;

  F77_CALL(dpotrf)(&lower, &k, L, &k, &info FCONE);
  if(info != 0)
    Rf_error("iid covariance Cholesky failed");
  s->reCovLogDet[b] = 0.0;
  for(j = 0; j < k; j++)
    s->reCovLogDet[b] += 2.0 * std::log(L[j + k * j]);
  F77_CALL(dpotrs)(&lower, &k, &k, L, &k, P, &k, &info FCONE);
  if(info != 0)
    Rf_error("iid covariance precision solve failed");
  for(col = 0; col < k; col++)
    for(j = 0; j < col; j++)
      L[j + k * col] = 0.0;
}

static void update_iid_covariance_gibbs(SamplerState *s, int b)
{
  int k, ng, j, g, col, row, info;
  char left, lower, noTrans, trans, nonUnit;
  double one, zero;
  double *L, *A, *T;
  const double *alpha;

  k = s->reCovDim[b];
  ng = s->reBlockSize[b] / k;
  left = 'L';
  lower = 'L';
  noTrans = 'N';
  trans = 'T';
  nonUnit = 'N';
  one = 1.0;
  zero = 0.0;

  L = (double*)R_alloc(k * k, sizeof(double));
  A = (double*)R_alloc(k * k, sizeof(double));
  T = (double*)R_alloc(k * k, sizeof(double));

  for(j = 0; j < k * k; j++){
    L[j] = s->reCovScale[b][j];
    A[j] = 0.0;
  }
  for(g = 0; g < ng; g++){
    alpha = s->alpha + s->reOffset[b] + g * k;
    for(col = 0; col < k; col++)
      for(row = col; row < k; row++)
        L[row + k * col] += alpha[row] * alpha[col];
  }
  F77_CALL(dpotrf)(&lower, &k, L, &k, &info FCONE);
  if(info != 0)
    Rf_error("iid posterior scale Cholesky failed");

  /* A A' ~ Wishart(df + ng, I); solve A T = L', then Sigma = T' T. */
  for(col = 0; col < k; col++){
    A[col + k * col] = std::sqrt(Rf_rchisq(s->reCovDF[b] + ng - col));
    for(row = col + 1; row < k; row++)
      A[row + k * col] = norm_rand();
    for(row = 0; row < k; row++)
      T[row + k * col] = row <= col ? L[col + k * row] : 0.0;
  }
  F77_CALL(dtrsm)(&left, &lower, &noTrans, &nonUnit,
                  &k, &k, &one, A, &k, T, &k FCONE FCONE FCONE FCONE);
  F77_CALL(dsyrk)(&lower, &trans, &k, &k,
                  &one, T, &k, &zero, s->reCovariance[b], &k FCONE FCONE);
  for(col = 0; col < k; col++)
    for(row = 0; row < col; row++)
      s->reCovariance[b][row + k * col] = s->reCovariance[b][col + k * row];
  refresh_iid_covariance_precision(s, b);
}

void update_sigmaSqRE_gibbs(SamplerState *s)
{
  int b, j, k, inc, mb;
  double ss, shape, scale;

  if(s->nRE <= 0 || s->q <= 0)
    return;

  inc = 1;

  for(b = 0; b < s->nRE; b++){
    if(s->reCovDim[b] > 0){
      update_iid_covariance_gibbs(s, b);
      continue;
    }
    mb = s->reBlockSize[b];
    k = 0;
    for(j = 0; j < s->q; j++){
      if(s->reBlockID[j] == b)
        s->alphaBlockBuf[k++] = s->alpha[j];
    }

    if(k != mb)
      Rf_error("update_sigmaSqRE_gibbs: block gather size mismatch");

    ss = 0.0;
    if(mb > 0)
      ss = F77_CALL(ddot)(&mb, s->alphaBlockBuf, &inc, s->alphaBlockBuf, &inc);

    shape = s->sigmaSqREShape[b] + 0.5 * (double)mb;
    scale = s->sigmaSqREScale[b]  + 0.5 * ss;

    s->sigmaSqRE[b] = 1.0 / Rf_rgamma(shape, 1.0 / scale);
  }
}

void refactor_current_M(SamplerState *s, const char *where)
{
  assemble_M_lat_numeric(s);
  s->cm.status = CHOLMOD_OK;
  if(!M_cholmod_factorize(s->M_lat, s->M_lat_fac, &s->cm) || s->cm.status != CHOLMOD_OK)
    Rf_error("%s: cholmod factorize failed for latent M", where);
}

void update_linear_predictor(SamplerState *s,
                             const double *beta,
                             const double *alpha,
                             const double *w,
                             double *eta)
{
  int i, inc;
  double one, zero;

  inc = 1;
  one = 1.0;
  zero = 0.0;

  for(i = 0; i < s->n; i++)
    eta[i] = 0.0;

  if(s->p > 0)
    F77_CALL(dgemv)("N", &s->n, &s->p,
                    &one, s->X, &s->n, beta, &inc,
                    &zero, eta, &inc FCONE);

  if(s->q > 0){
    sparse_Z_mult(s, alpha, s->Za);
    for(i = 0; i < s->n; i++)
      eta[i] += s->Za[i];
  }

  if(s->qLatTotal > 0){
    apply_A(s, w, s->nWork1);
    for(i = 0; i < s->n; i++)
      eta[i] += s->nWork1[i];
  }
}

int is_pg_likelihood(const SamplerState *s)
{
  return s->likelihoodFamily == LIKELIHOOD_BINOMIAL ||
         s->likelihoodFamily == LIKELIHOOD_NEGATIVE_BINOMIAL;
}

int is_probit_likelihood(const SamplerState *s)
{
  return s->likelihoodFamily == LIKELIHOOD_PROBIT;
}

int is_augmented_likelihood(const SamplerState *s)
{
  return is_pg_likelihood(s) || is_probit_likelihood(s);
}

static double rstdnorm_lower_trunc(double lower)
{
  double x, alpha, u, logu;

  if(!R_FINITE(lower))
    Rf_error("non-finite probit truncation point");

  if(lower <= 0.0){
    do {
      x = rnorm(0.0, 1.0);
    } while(x < lower);
    return x;
  }

  alpha = 0.5 * (lower + std::sqrt(lower * lower + 4.0));
  for(;;){
    u = runif(0.0, 1.0);
    if(u <= 0.0)
      continue;
    x = lower - std::log(u) / alpha;
    logu = std::log(runif(0.0, 1.0));
    if(logu <= -0.5 * (x - alpha) * (x - alpha))
      return x;
  }
}

static double rnorm_probit_latent(double mean, int y)
{
  double lower;

  if(!R_FINITE(mean))
    Rf_error("non-finite probit linear predictor");

  lower = -mean;
  if(y)
    return mean + rstdnorm_lower_trunc(lower);

  return mean - rstdnorm_lower_trunc(-lower);
}

void update_probit_working_model(SamplerState *s,
                                 const double *beta,
                                 const double *alpha,
                                 const double *w)
{
  int i, y;
  double obsOffset, latentFull;

  if(!is_probit_likelihood(s))
    return;

  /*
    Albert-Chib augmentation:

      y_i = 1(z_i > 0),  z_i = offset_i + eta_i + e_i, e_i ~ N(0, 1).

    The collapsed Gaussian sampler works on eta_i without the offset, so the
    mutable response below stores z_i - offset_i and unit observation precision.
  */
  update_linear_predictor(s, beta, alpha, w, s->nWork2);

  for(i = 0; i < s->n; i++){
    obsOffset = (s->offset != NULL) ? s->offset[i] : 0.0;
    y = s->yObserved[i] > 0.0 ? 1 : 0;
    latentFull = rnorm_probit_latent(s->nWork2[i] + obsOffset, y);
    s->y[i] = latentFull - obsOffset;
    s->obsPrecision[i] = 1.0;
  }

  refactor_current_M(s, "update_probit_working_model");
}

void update_pg_working_model(SamplerState *s,
                             BayesLogit_rpg_hybrid_t pg,
                             const double *beta,
                             const double *alpha,
                             const double *w)
{
  int i;
  double omega, kappa, shape, tilt, nbOffset, obsOffset;

  if(!is_pg_likelihood(s))
    return;

  /*
    Polya-Gamma augmentation turns supported logit-type likelihoods into a
    Gaussian working model conditional on omega.

    Binomial:

      kappa_i = y_i - n_i / 2
      z_i     = kappa_i / omega_i
      z_i     = eta_i + e_i,  e_i ~ N(0, 1 / omega_i)

    Negative binomial with fixed size r and eta_i = log(mu_i):

      psi_i   = eta_i - log(r)
      kappa_i = (y_i - r) / 2
      z_i     = kappa_i / omega_i + log(r)
      z_i     = eta_i + e_i,  e_i ~ N(0, 1 / omega_i)

    If the formula has a known offset o_i, the PG draw uses the full
    eta_i = o_i + sampled linear predictor. The Gaussian working model below
    remains on the sampled linear predictor, so o_i is subtracted from z_i.

    The existing collapsed Gaussian code already works with a diagonal
    observation precision W. For PG models we store omega_i in that same
    obsPrecision slot. Thus omega_i replaces tau_sq^{-1}; it is not tau_sq_i and
    it is not a residual variance.

    After the marginal updates, the iteration schedule draws compatible latent
    effects so that the next augmentation uses eta = X beta + Z alpha + A w.
  */
  update_linear_predictor(s, beta, alpha, w, s->nWork2);

  for(i = 0; i < s->n; i++){
    obsOffset = (s->offset != NULL) ? s->offset[i] : 0.0;
    if(s->likelihoodFamily == LIKELIHOOD_BINOMIAL){
      shape = (double)s->nTrial[i];
      tilt = s->nWork2[i] + obsOffset;
      kappa = s->yObserved[i] - 0.5 * shape;
      nbOffset = 0.0;
    } else if(s->likelihoodFamily == LIKELIHOOD_NEGATIVE_BINOMIAL){
      shape = s->yObserved[i] + s->nbSize;
      tilt = s->nWork2[i] + obsOffset - std::log(s->nbSize);
      kappa = 0.5 * (s->yObserved[i] - s->nbSize);
      nbOffset = std::log(s->nbSize);
    } else {
      Rf_error("unsupported Polya-Gamma likelihood");
    }

    omega = pg(shape, tilt);
    if(!R_FINITE(omega) || !(omega > 0.0))
      Rf_error("Polya-Gamma draw produced a non-positive observation precision (shape=%.17g, tilt=%.17g)", shape, tilt);

    s->obsPrecision[i] = omega;
    s->y[i] = kappa / omega + nbOffset - obsOffset;
  }

  refactor_current_M(s, "update_pg_working_model");
}

double log_ig_kernel_collapsed(double x, double shape, double scale)
{
  return -(shape + 1.0) * std::log(x) - scale / x;
}

double log_prior_kernel_original(double x, PriorSpec prior)
{
  double z;

  if(!(x > 0.0) && prior.family != PRIOR_UNIFORM)
    return R_NegInf;

  switch(prior.family){
  case PRIOR_IG:
    return log_ig_kernel_collapsed(x, prior.p1, prior.p2);
  case PRIOR_UNIFORM:
    return 0.0;
  case PRIOR_LOG_NORMAL:
    z = (std::log(x) - prior.p1) / prior.p2;
    return -std::log(x) - 0.5 * z * z;
  case PRIOR_GAMMA:
    return (prior.p1 - 1.0) * std::log(x) - prior.p2 * x;
  case PRIOR_HALF_NORMAL:
    return -0.5 * x / (prior.p1 * prior.p1) - 0.5 * std::log(x);
  case PRIOR_HALF_T:
    return -0.5 * std::log(x) -
      0.5 * (prior.p1 + 1.0) * std::log1p(x / (prior.p1 * prior.p2 * prior.p2));
  case PRIOR_BETA:
    if(!(x > 0.0 && x < 1.0))
      return R_NegInf;
    return (prior.p1 - 1.0) * std::log(x) +
      (prior.p2 - 1.0) * std::log(1.0 - x);
  default:
    return R_NegInf;
  }
}

double log_cov_params_target_collapsed(SamplerState *s, const double *resid)
{
  int i, j;
  double out;

  out = compute_logDetV(s);
  out = -0.5 * (out + quadform_Vinv(s, resid));

  if(is_augmented_likelihood(s)){
    /* No tau_sq term: augmented likelihoods own their observation precision. */
  } else if(s->residualModel == 0 && s->tauTune > 0.0){
    out += log_prior_kernel_original(s->tauSq, s->tauSqPrior)
         + std::log(s->tauSq);
  } else if(s->residualModel == 2){
    for(i = 0; i < s->nResidualGroup; i++)
      if(s->residualVarianceTune[i] > 0.0)
        out += log_prior_kernel_original(s->residualVariance[i],
                                         s->residualVariancePrior[i])
             + std::log(s->residualVariance[i]);
  } else if(s->residualModel == 3){
    for(i = 0; i < s->nResidualGroup; i++)
      if(s->residualVarianceTune[i] > 0.0){
        double z = (std::log(s->residualVariance[i]) - s->residualVarianceMeanLog[i]) /
                   s->residualVarianceSdLog[i];
        out += -0.5 * z * z;
      }
  }

  for(i = 0; i < s->nTerms; i++){
    TermState *term = s->terms + i;
    if(term->sigmaSqTune > 0.0)
      out += log_prior_kernel_original(term->sigmaSq, term->sigmaSqPrior)
           + std::log(term->sigmaSq);
    for(j = 0; j < term->thetaDim; j++)
      if(term->thetaTune[j] > 0.0)
        out += log_prior_kernel_original(term->theta[j], term->thetaPrior[j])
             + theta_log_jacobian(term->theta[j], term->thetaLower[j], term->thetaUpper[j]);
  }

  if(s->collapseIid)
    out += log_iid_covariance_target(s);

  return out;
}

int log_cov_params_target_collapsed_try(SamplerState *s, const double *resid, double *out)
{
  int i, j;
  double logdetQw, logdetR, logdetV;

  if(!logdet_Qw_blocks_try(s, &logdetQw))
    return 0;
  if(s->collapseIid)
    logdetQw += logdet_Qiid_blocks(s);

  logdetR = 0.0;
  for(i = 0; i < s->n; i++)
    logdetR -= std::log(s->obsPrecision[i]);
  logdetV = logdetR + logDetFactor(s->M_lat_fac) - logdetQw;

  *out = -0.5 * (logdetV + quadform_Vinv(s, resid));

  if(is_augmented_likelihood(s)){
    /* No tau_sq term: augmented likelihoods own their observation precision. */
  } else if(s->residualModel == 0 && s->tauTune > 0.0){
    *out += log_prior_kernel_original(s->tauSq, s->tauSqPrior)
          + std::log(s->tauSq);
  } else if(s->residualModel == 2){
    for(i = 0; i < s->nResidualGroup; i++)
      if(s->residualVarianceTune[i] > 0.0)
        *out += log_prior_kernel_original(s->residualVariance[i],
                                          s->residualVariancePrior[i])
              + std::log(s->residualVariance[i]);
  } else if(s->residualModel == 3){
    for(i = 0; i < s->nResidualGroup; i++)
      if(s->residualVarianceTune[i] > 0.0){
        double z = (std::log(s->residualVariance[i]) - s->residualVarianceMeanLog[i]) /
                   s->residualVarianceSdLog[i];
        *out += -0.5 * z * z;
      }
  }

  for(i = 0; i < s->nTerms; i++){
    TermState *term = s->terms + i;
    if(term->sigmaSqTune > 0.0)
      *out += log_prior_kernel_original(term->sigmaSq, term->sigmaSqPrior)
            + std::log(term->sigmaSq);
    for(j = 0; j < term->thetaDim; j++)
      if(term->thetaTune[j] > 0.0)
        *out += log_prior_kernel_original(term->theta[j], term->thetaPrior[j])
              + theta_log_jacobian(term->theta[j], term->thetaLower[j], term->thetaUpper[j]);
  }

  if(s->collapseIid)
    *out += log_iid_covariance_target(s);

  return R_FINITE(*out);
}

typedef struct {
  double tauSq;
  double *residualVariance;
  cholmod_sparse *M_lat;
  cholmod_factor *M_lat_fac;
  double *sigmaSq;
  double **theta;
  double **B;
  double **F;
  double **Q;
  double *logDetQ;
  double *sigmaSqRE;
  double **reCovariance, **reCovPrecision, **reCovChol;
  double *reCovLogDet;
} CovBindingSave;

CovBindingSave bind_cov_proposal_state(SamplerState *s, CovProposalState *p)
{
  int t;
  CovBindingSave save;

  save.tauSq = s->tauSq;
  save.residualVariance = (double*)R_alloc(s->nResidualGroup > 0 ? s->nResidualGroup : 1, sizeof(double));
  for(t = 0; t < s->nResidualGroup; t++)
    save.residualVariance[t] = s->residualVariance[t];
  save.sigmaSqRE = s->sigmaSqRE;
  save.reCovariance = s->reCovariance;
  save.reCovPrecision = s->reCovPrecision;
  save.reCovChol = s->reCovChol;
  save.reCovLogDet = s->reCovLogDet;
  save.M_lat = s->M_lat;
  save.M_lat_fac = s->M_lat_fac;
  save.sigmaSq = (double*)R_alloc(s->nTerms > 0 ? s->nTerms : 1, sizeof(double));
  save.theta = (double**)R_alloc(s->nTerms > 0 ? s->nTerms : 1, sizeof(double*));
  save.B = (double**)R_alloc(s->nTerms > 0 ? s->nTerms : 1, sizeof(double*));
  save.F = (double**)R_alloc(s->nTerms > 0 ? s->nTerms : 1, sizeof(double*));
  save.Q = (double**)R_alloc(s->nTerms > 0 ? s->nTerms : 1, sizeof(double*));
  save.logDetQ = (double*)R_alloc(s->nTerms > 0 ? s->nTerms : 1, sizeof(double));

  if(s->collapseIid){
    s->sigmaSqRE = p->sigmaSqRE;
    s->reCovariance = p->reCovariance;
    s->reCovPrecision = p->reCovPrecision;
    s->reCovChol = p->reCovChol;
    s->reCovLogDet = p->reCovLogDet;
  }
  s->tauSq = p->tauSq;
  for(t = 0; t < s->nResidualGroup; t++)
    s->residualVariance[t] = p->residualVariance[t];
  s->M_lat = p->M_lat;
  s->M_lat_fac = p->M_lat_fac;
  refresh_observation_precision(s);

  for(t = 0; t < s->nTerms; t++){
    TermState *term = s->terms + t;

    save.sigmaSq[t] = term->sigmaSq;
    save.theta[t] = term->theta;
    save.B[t] = term->B;
    save.F[t] = term->F;
    save.Q[t] = term->Q;
    save.logDetQ[t] = term->logDetQ;

    term->sigmaSq = p->sigmaSq[t];
    term->theta = p->theta[t];
    if(p->B[t] != NULL) term->B = p->B[t];
    if(p->F[t] != NULL) term->F = p->F[t];
    if(p->Q[t] != NULL) term->Q = p->Q[t];
    term->logDetQ = p->logDetQ[t];
  }

  return save;
}

void restore_cov_binding(SamplerState *s, CovBindingSave *save)
{
  int t;

  s->sigmaSqRE = save->sigmaSqRE;
  s->reCovariance = save->reCovariance;
  s->reCovPrecision = save->reCovPrecision;
  s->reCovChol = save->reCovChol;
  s->reCovLogDet = save->reCovLogDet;
  s->tauSq = save->tauSq;
  for(t = 0; t < s->nResidualGroup; t++)
    s->residualVariance[t] = save->residualVariance[t];
  s->M_lat = save->M_lat;
  s->M_lat_fac = save->M_lat_fac;
  refresh_observation_precision(s);

  for(t = 0; t < s->nTerms; t++){
    TermState *term = s->terms + t;
    term->sigmaSq = save->sigmaSq[t];
    term->theta = save->theta[t];
    term->B = save->B[t];
    term->F = save->F[t];
    term->Q = save->Q[t];
    term->logDetQ = save->logDetQ[t];
  }
}

void accept_cov_proposal_state(SamplerState *s, CovProposalState *p)
{
  int t;
  double tmpDouble;
  double *tmpPtr;
  double **tmpMatrices;
  cholmod_sparse *tmpSparse;
  cholmod_factor *tmpFactor;

  if(s->collapseIid){
    tmpPtr = s->sigmaSqRE;
    s->sigmaSqRE = p->sigmaSqRE;
    p->sigmaSqRE = tmpPtr;
    tmpMatrices = s->reCovariance;
    s->reCovariance = p->reCovariance;
    p->reCovariance = tmpMatrices;
    tmpMatrices = s->reCovPrecision;
    s->reCovPrecision = p->reCovPrecision;
    p->reCovPrecision = tmpMatrices;
    tmpMatrices = s->reCovChol;
    s->reCovChol = p->reCovChol;
    p->reCovChol = tmpMatrices;
    tmpPtr = s->reCovLogDet;
    s->reCovLogDet = p->reCovLogDet;
    p->reCovLogDet = tmpPtr;
  }

  tmpDouble = s->tauSq;
  s->tauSq = p->tauSq;
  p->tauSq = tmpDouble;
  for(t = 0; t < s->nResidualGroup; t++){
    tmpDouble = s->residualVariance[t];
    s->residualVariance[t] = p->residualVariance[t];
    p->residualVariance[t] = tmpDouble;
  }
  refresh_observation_precision(s);

  tmpSparse = s->M_lat;
  s->M_lat = p->M_lat;
  p->M_lat = tmpSparse;

  tmpFactor = s->M_lat_fac;
  s->M_lat_fac = p->M_lat_fac;
  p->M_lat_fac = tmpFactor;

  for(t = 0; t < s->nTerms; t++){
    TermState *term = s->terms + t;

    tmpDouble = term->sigmaSq;
    term->sigmaSq = p->sigmaSq[t];
    p->sigmaSq[t] = tmpDouble;

    tmpPtr = term->theta;
    term->theta = p->theta[t];
    p->theta[t] = tmpPtr;

    tmpPtr = term->B;
    term->B = p->B[t];
    p->B[t] = tmpPtr;

    tmpPtr = term->F;
    term->F = p->F[t];
    p->F[t] = tmpPtr;

    tmpPtr = term->Q;
    term->Q = p->Q[t];
    p->Q[t] = tmpPtr;

    tmpDouble = term->logDetQ;
    term->logDetQ = p->logDetQ[t];
    p->logDetQ[t] = tmpDouble;
  }
}

void reset_proposal_factor(SamplerState *s, CovProposalState *p)
{
  if(p->M_lat_fac != NULL)
    M_cholmod_free_factor(&p->M_lat_fac, &s->cm);

  s->cm.status = CHOLMOD_OK;
  p->M_lat_fac = M_cholmod_copy_factor(s->M_lat_sym, &s->cm);
  if(p->M_lat_fac == NULL || s->cm.status != CHOLMOD_OK)
    Rf_error("reset_proposal_factor: failed to copy symbolic M factor");
}

void pack_cov_eta(SamplerState *s, double *eta)
{
  int k, t, j;
  CovProposalState *p;

  p = &s->covProp;
  for(k = 0; k < p->dim; k++){
    if(p->paramType[k] == 0){
      eta[k] = std::log(s->tauSq);
    } else if(p->paramType[k] == 3){
      t = p->paramTerm[k];
      eta[k] = std::log(s->residualVariance[t]);
    } else if(p->paramType[k] == 1){
      t = p->paramTerm[k];
      eta[k] = std::log(s->terms[t].sigmaSq);
    } else if(p->paramType[k] == 2){
      t = p->paramTerm[k];
      j = p->paramTheta[k];
      eta[k] = theta_forward(s->terms[t].theta[j],
                             s->terms[t].thetaLower[j],
                             s->terms[t].thetaUpper[j]);
    } else if(p->paramType[k] == 4){
      eta[k] = std::log(s->sigmaSqRE[p->paramTerm[k]]);
    } else if(p->paramType[k] == 5){
      t = p->paramTerm[k];
      j = p->paramTheta[k];
      eta[k] = s->reCovChol[t][j];
      if(j % (s->reCovDim[t] + 1) == 0)
        eta[k] = std::log(eta[k]);
    } else {
      Rf_error("pack_cov_eta: unknown adaptive parameter type");
    }
  }
}

void pack_cov_block_eta(SamplerState *s, CovProposalBlock *b, double *eta)
{
  int i;
  double *allEta;

  allEta = (double*)R_alloc(s->covProp.dim > 0 ? s->covProp.dim : 1, sizeof(double));
  pack_cov_eta(s, allEta);
  for(i = 0; i < b->dim; i++)
    eta[i] = allEta[b->paramIndex[i]];
}

void unpack_cov_eta_to_proposal(SamplerState *s, const double *eta)
{
  int i, j, k, t;
  int *changed = NULL;
  CovProposalState *p;

  p = &s->covProp;
  if(s->collapseIid){
    copy_iid_covariance_proposal(s, p);
    changed = (int*)R_alloc(s->nRE, sizeof(int));
    for(i = 0; i < s->nRE; i++)
      changed[i] = 0;
  }
  p->tauSq = s->tauSq;
  for(i = 0; i < s->nResidualGroup; i++)
    p->residualVariance[i] = s->residualVariance[i];
  for(i = 0; i < s->nTerms; i++){
    TermState *term = s->terms + i;
    p->sigmaSq[i] = term->sigmaSq;
    for(j = 0; j < term->thetaDim; j++)
      p->theta[i][j] = term->theta[j];
  }

  for(k = 0; k < p->dim; k++){
    if(p->paramType[k] == 0){
      p->tauSq = std::exp(eta[k]);
    } else if(p->paramType[k] == 3){
      t = p->paramTerm[k];
      p->residualVariance[t] = std::exp(eta[k]);
    } else if(p->paramType[k] == 1){
      t = p->paramTerm[k];
      p->sigmaSq[t] = std::exp(eta[k]);
    } else if(p->paramType[k] == 2){
      t = p->paramTerm[k];
      j = p->paramTheta[k];
      p->theta[t][j] = theta_inverse(eta[k],
                                     s->terms[t].thetaLower[j],
                                     s->terms[t].thetaUpper[j]);
      if(!(s->terms[t].thetaLower[j] < p->theta[t][j] &&
           p->theta[t][j] < s->terms[t].thetaUpper[j]))
        p->theta[t][j] = s->terms[t].theta[j];
    } else if(p->paramType[k] == 4){
      p->sigmaSqRE[p->paramTerm[k]] = std::exp(eta[k]);
    } else if(p->paramType[k] == 5){
      t = p->paramTerm[k];
      j = p->paramTheta[k];
      if(j % (s->reCovDim[t] + 1) == 0){
        if(eta[k] != std::log(s->reCovChol[t][j])){
          p->reCovChol[t][j] = std::exp(eta[k]);
          changed[t] = 1;
        }
      } else if(eta[k] != s->reCovChol[t][j]){
        p->reCovChol[t][j] = eta[k];
        changed[t] = 1;
      }
    } else {
      Rf_error("unpack_cov_eta_to_proposal: unknown adaptive parameter type");
    }
  }
  if(s->collapseIid){
    for(t = 0; t < s->nRE; t++)
      if(changed[t] && !refresh_iid_covariance_proposal(s, p, t))
        p->iidValid = 0;
  }
}

void unpack_cov_block_eta_to_proposal(SamplerState *s, CovProposalBlock *b, const double *eta)
{
  int i;
  double *allEta;

  allEta = (double*)R_alloc(s->covProp.dim > 0 ? s->covProp.dim : 1, sizeof(double));
  pack_cov_eta(s, allEta);
  for(i = 0; i < b->dim; i++)
    allEta[b->paramIndex[i]] = eta[i];
  unpack_cov_eta_to_proposal(s, allEta);
}

static void draw_adaptive_block_proposal(CovProposalBlock *b)
{
  int i, j;

  for(i = 0; i < b->dim; i++)
    b->z[i] = rnorm(0.0, 1.0);

  for(i = 0; i < b->dim; i++){
    b->etaProp[i] = b->eta[i];
    for(j = 0; j <= i; j++)
      b->etaProp[i] += b->proposalChol[i + b->dim * j] * b->z[j];
  }
}

void draw_adaptive_cov_block_proposal(SamplerState *s, CovProposalBlock *b)
{
  pack_cov_block_eta(s, b, b->eta);
  draw_adaptive_block_proposal(b);
  unpack_cov_block_eta_to_proposal(s, b, b->etaProp);
}

void factor_adaptive_block_proposal_covariance(CovProposalBlock *b)
{
  int i, j, info;
  char lower = 'L';

  for(j = 0; j < b->dim * b->dim; j++){
    if(b->scalarMode)
      b->proposalCov[j] = b->sigmaSqM * b->sigmaSqM * b->Sigma0[j];
    else
      b->proposalCov[j] = b->sigmaSqM * b->Sigma0[j];
    b->proposalChol[j] = b->proposalCov[j];
  }

  F77_CALL(dpotrf)(&lower, &b->dim, b->proposalChol, &b->dim, &info FCONE);
  if(info != 0)
    Rf_error("adaptive covariance proposal is not positive definite");

  for(j = 0; j < b->dim; j++)
    for(i = 0; i < j; i++)
      b->proposalChol[i + b->dim * j] = 0.0;
}

void update_adaptive_covariance_block_state(SamplerState *s, CovProposalBlock *b, int accepted)
{
  int i, j, n;
  double rhat, gamma1, gamma2, mean_i, mean_j, cov_ij;
  int bb;

  pack_cov_block_eta(s, b, b->eta);

  for(i = 0; i < b->dim; i++)
    b->batchSamples[b->batchPos + b->batchLength * i] = b->eta[i];
  b->batchAccept += accepted;
  b->batchPos++;
  b->acceptCount += accepted;

  (void)s;

  if(b->batchPos < b->batchLength)
    return;

  b->batchIndex++;
  n = b->batchLength;
  rhat = (double)b->batchAccept / (double)n;
  b->lastBatchAccept = rhat;
  if(b->batchIndex <= b->maxBatchHistory)
    b->batchAcceptHistory[b->batchIndex - 1] = rhat;

  gamma1 = 1.0 / std::pow((double)(b->batchIndex + 1), b->c1);
  gamma2 = b->c0 * gamma1;
  b->sigmaSqM = std::exp(std::log(b->sigmaSqM) + gamma2 * (rhat - b->targetAccept));
  if(b->batchIndex <= b->maxBatchHistory)
    b->proposalScaleHistory[b->batchIndex - 1] = b->sigmaSqM;

  for(j = 0; j < b->dim; j++){
    mean_j = 0.0;
    for(bb = 0; bb < n; bb++)
      mean_j += b->batchSamples[bb + b->batchLength * j];
    mean_j /= (double)n;

    for(i = j; i < b->dim; i++){
      mean_i = 0.0;
      for(bb = 0; bb < n; bb++)
        mean_i += b->batchSamples[bb + b->batchLength * i];
      mean_i /= (double)n;

      cov_ij = 0.0;
      for(bb = 0; bb < n; bb++)
        cov_ij += (b->batchSamples[bb + b->batchLength * i] - mean_i) *
                  (b->batchSamples[bb + b->batchLength * j] - mean_j);
      cov_ij /= (double)(n - 1);

      b->Sigma0[i + b->dim * j] += gamma1 * (cov_ij - b->Sigma0[i + b->dim * j]);
      b->Sigma0[j + b->dim * i] = b->Sigma0[i + b->dim * j];
    }
  }

  factor_adaptive_block_proposal_covariance(b);
  b->batchPos = 0;
  b->batchAccept = 0;
}

void update_adaptive_scalar_block_state(CovProposalBlock *b, int accepted)
{
  int n;
  double rhat, step;

  if(b->dim != 1)
    Rf_error("update_adaptive_scalar_block_state: scalar block has wrong dimension");

  b->batchSamples[b->batchPos] = b->eta[0];
  b->batchAccept += accepted;
  b->batchPos++;
  b->acceptCount += accepted;

  if(b->batchPos < b->batchLength)
    return;

  b->batchIndex++;
  n = b->batchLength;
  rhat = (double)b->batchAccept / (double)n;
  b->lastBatchAccept = rhat;
  if(b->batchIndex <= b->maxBatchHistory)
    b->batchAcceptHistory[b->batchIndex - 1] = rhat;

  step = std::min(0.01, 1.0 / std::sqrt((double)b->batchIndex));
  if(rhat > b->targetAccept)
    b->sigmaSqM = std::exp(std::log(b->sigmaSqM) + step);
  else
    b->sigmaSqM = std::exp(std::log(b->sigmaSqM) - step);

  if(b->batchIndex <= b->maxBatchHistory)
    b->proposalScaleHistory[b->batchIndex - 1] = b->sigmaSqM;

  factor_adaptive_block_proposal_covariance(b);
  b->batchPos = 0;
  b->batchAccept = 0;
}

static double softplus(double x)
{
  return x > 0.0 ? x + std::log1p(std::exp(-x)) : std::log1p(std::exp(x));
}

static double log_nb_size_target(SamplerState *s, double logSize, const double *eta)
{
  int i;
  double size, out, logMu, mu, value, shift;

  size = std::exp(logSize);
  if(!R_FINITE(size) || !(size > 0.0))
    return R_NegInf;
  out = log_prior_kernel_original(size, s->nbSizePrior) + logSize;
  for(i = 0; i < s->n; i++){
    logMu = eta[i] + (s->offset != NULL ? s->offset[i] : 0.0);
    if(!R_FINITE(logMu))
      return R_NegInf;
    mu = std::exp(logMu);
    if(R_FINITE(mu) && mu > 0.0){
      value = Rf_dnbinom_mu(s->yObserved[i], size, mu, 1);
    } else {
      /* Shift a stable Rmath density at mean one when exp(logMu) is not
         representable. This avoids differences of large log-gamma values. */
      shift = softplus(-logSize) - softplus(logMu - logSize);
      value = Rf_dnbinom_mu(s->yObserved[i], size, 1.0, 1) +
        size * shift + s->yObserved[i] * (logMu + shift);
    }
    out += value;
  }
  return R_FINITE(out) ? out : R_NegInf;
}

static int update_nb_size(SamplerState *s, const double *beta,
                          const double *alpha, const double *w)
{
  CovProposalBlock *b;
  double current, proposal, logAccept;
  int accepted;

  if(!s->nbSizeEstimate)
    return 0;
  b = &s->nbSizeProposal;
  update_linear_predictor(s, beta, alpha, w, s->nWork2);
  b->eta[0] = std::log(s->nbSize);
  draw_adaptive_block_proposal(b);
  current = log_nb_size_target(s, b->eta[0], s->nWork2);
  if(!R_FINITE(current))
    Rf_error("non-finite negative-binomial size target at the current state");
  proposal = log_nb_size_target(s, b->etaProp[0], s->nWork2);
  logAccept = proposal - current;
  accepted = std::log(runif(0.0, 1.0)) < logAccept;
  if(accepted)
    s->nbSize = std::exp(b->etaProp[0]);
  b->eta[0] = std::log(s->nbSize);
  return accepted;
}

int cov_params_block_mh_collapsed_step(SamplerState *s, CovProposalBlock *block, const double *resid)
{
  int i, factorized, factorStatus, tryCatch;
  double logPostOld, logPostProp, logu;
  CovProposalState *p;
  CovBindingSave save;
  int accepted;

  if(!s->covProp.initialized)
    Rf_error("cov_params_block_mh_collapsed_step: proposal state is not initialized");

  p = &s->covProp;
  if(p->dim == 0 || block->dim == 0)
    return 0;

  accepted = 0;
  logPostOld = log_cov_params_target_collapsed(s, resid);
  if(s->collapseIid && !R_FINITE(logPostOld))
    Rf_error("non-finite collapsed IID covariance target at the current state");

  draw_adaptive_cov_block_proposal(s, block);
  if(s->collapseIid && !p->iidValid)
    return 0;
  if(s->collapseIid)
    for(i = 0; i < s->nRE; i++)
      if(!R_FINITE(p->sigmaSqRE[i]) || p->sigmaSqRE[i] <= 0.0 ||
         !R_FINITE(1.0 / p->sigmaSqRE[i]))
        return 0;

  for(i = 0; i < s->nTerms; i++){
    TermState *term = s->terms + i;
    GraphState *g = s->graphs + term->graphIndex;
    TermState tmpTerm;

    tmpTerm = *term;
    tmpTerm.sigmaSq = p->sigmaSq[i];
    tmpTerm.theta = p->theta[i];
    tmpTerm.B = p->B[i];
    tmpTerm.F = p->F[i];
    tmpTerm.Q = p->Q[i];
    tmpTerm.logDetQ = p->logDetQ[i];

    if(g->type == GRAPH_NNGP)
      update_nngp_BF(s, g, &tmpTerm);
    else if(g->type == GRAPH_GP)
      update_gp_Q(s, g, &tmpTerm);

    p->logDetQ[i] = tmpTerm.logDetQ;
  }

  save = bind_cov_proposal_state(s, p);
  assemble_M_lat_numeric(s);
  s->cm.status = CHOLMOD_OK;
  /* Extreme IID proposals can lose positive definiteness in floating point.
     Handle that status as a rejection; retain errors for other failures. */
  tryCatch = s->cm.try_catch;
  if(s->collapseIid)
    s->cm.try_catch = 1;
  factorized = M_cholmod_factorize(s->M_lat, s->M_lat_fac, &s->cm);
  factorStatus = s->cm.status;
  s->cm.try_catch = tryCatch;
  if(!factorized || factorStatus != CHOLMOD_OK){
    s->cm.status = CHOLMOD_OK;
    restore_cov_binding(s, &save);
    reset_proposal_factor(s, p);
    if(s->collapseIid && factorStatus != CHOLMOD_NOT_POSDEF && factorStatus != CHOLMOD_DSMALL)
      Rf_error("collapsed IID proposal factorization failed (CHOLMOD status %d)", factorStatus);
    return 0;
  }

  if(!log_cov_params_target_collapsed_try(s, resid, &logPostProp)){
    restore_cov_binding(s, &save);
    reset_proposal_factor(s, p);
    return 0;
  }
  restore_cov_binding(s, &save);

  logu = std::log(runif(0.0, 1.0));
  if(logu < (logPostProp - logPostOld)){
    accept_cov_proposal_state(s, p);
    accepted = 1;
  }

  return accepted;
}

int update_cov_params_mh_collapsed(SamplerState *s, const double *resid)
{
  int b, accepted, totalAccepted;
  CovProposalState *p;

  p = &s->covProp;
  if(p->dim <= 0 || p->nBlocks <= 0)
    return 0;

  totalAccepted = 0;
  for(b = 0; b < p->nBlocks; b++){
    accepted = cov_params_block_mh_collapsed_step(s, p->blocks + b, resid);
    if(s->metropolisBlocking == 5){
      pack_cov_block_eta(s, p->blocks + b, p->blocks[b].eta);
      update_adaptive_scalar_block_state(p->blocks + b, accepted);
    } else
      update_adaptive_covariance_block_state(s, p->blocks + b, accepted);
    totalAccepted += accepted;
    p->totalBlockAccept += accepted;
    s->covParamAttempts++;
  }

  return totalAccepted;
}

void update_resid_for_current_fixed_random(SamplerState *s,
                                           const double *betaCurrent,
                                           const double *alphaCurrent,
                                           double *resid)
{
  int i, inc;
  double minusOne, one;

  minusOne = -1.0;
  one = 1.0;
  inc = 1;

  for(i = 0; i < s->n; i++)
    resid[i] = s->y[i];

  if(s->p > 0){
    F77_CALL(dgemv)("N", &s->n, &s->p,
                    &minusOne, s->X, &s->n, betaCurrent, &inc,
                    &one, resid, &inc FCONE);
  }

  if(s->q > 0)
    sparse_Z_mult(s, alphaCurrent, s->Za);
  else
    for(i = 0; i < s->n; i++) s->Za[i] = 0.0;

  for(i = 0; i < s->n; i++)
    resid[i] -= s->Za[i];
}

struct SamplerIterationResult {
  int covAccepted;
  int sizeAccepted;
  int wCurrentFresh;
};

static int update_covariance_iteration(SamplerState *s, const double *resid, int warmup)
{
  CovProposalState *p;
  int block, accepted, acceptedTotal;

  p = &s->covProp;
  acceptedTotal = 0;
  if(warmup){
    for(block = 0; block < p->nBlocks; block++){
      accepted = cov_params_block_mh_collapsed_step(s, p->blocks + block, resid);
      acceptedTotal += accepted;
      p->blocks[block].batchAccept += accepted;
      p->blocks[block].warmupNAccepted += accepted;
      p->blocks[block].warmupNAttempted++;
    }
    p->warmupNAccepted += acceptedTotal;
    p->warmupNAttempted += p->nBlocks;
  } else {
    acceptedTotal = update_cov_params_mh_collapsed(s, resid);
    s->covParamAccept += acceptedTotal;
  }
  return acceptedTotal;
}

static int update_likelihood_iteration(SamplerState *s,
                                       const double *betaCurrent,
                                       const double *alphaCurrent,
                                       const double *wCurrent,
                                       BayesLogit_rpg_hybrid_t pg,
                                       int warmup)
{
  int accepted;

  accepted = 0;
  /* The size step integrates out omega; redraw omega immediately before
     any Gaussian working-model update. */
  if(s->nbSizeEstimate){
    accepted = update_nb_size(s, betaCurrent, alphaCurrent, wCurrent);
    if(warmup){
      s->nbSizeProposal.warmupNAttempted++;
      s->nbSizeProposal.warmupNAccepted += accepted;
    } else
      update_adaptive_scalar_block_state(&s->nbSizeProposal, accepted);
  }
  if(is_pg_likelihood(s))
    update_pg_working_model(s, pg, betaCurrent, alphaCurrent, wCurrent);
  else if(is_probit_likelihood(s))
    update_probit_working_model(s, betaCurrent, alphaCurrent, wCurrent);
  return accepted;
}

static SamplerIterationResult update_separate_iid_iteration(SamplerState *s,
                                                           double *betaCurrent,
                                                           double *alphaCurrent,
                                                           double *resid,
                                                           double *wCurrent,
                                                           BayesLogit_rpg_hybrid_t pg,
                                                           int warmup)
{
  SamplerIterationResult result;
  int i;

  result.covAccepted = 0;
  result.wCurrentFresh = 0;
  result.sizeAccepted = update_likelihood_iteration(s, betaCurrent, alphaCurrent,
                                                     wCurrent, pg, warmup);

  if(s->p > 0)
    update_beta_draw_collapsed(s, alphaCurrent, betaCurrent);
  if(s->q > 0){
    update_alpha_draw_collapsed(s, betaCurrent, alphaCurrent);
    for(i = 0; i < s->q; i++)
      s->alpha[i] = alphaCurrent[i];
  }
  if(s->nRE > 0)
    update_sigmaSqRE_gibbs(s);

  update_resid_for_current_fixed_random(s, betaCurrent, alphaCurrent, resid);
  result.covAccepted = update_covariance_iteration(s, resid, warmup);

  if(is_augmented_likelihood(s) && s->qLatTotal > 0){
    update_resid_for_current_fixed_random(s, betaCurrent, alphaCurrent, resid);
    recover_w_draw_collapsed(s, resid, wCurrent);
    result.wCurrentFresh = 1;
  }

  return result;
}

static SamplerIterationResult update_collapsed_iid_iteration(SamplerState *s,
                                                            double *betaCurrent,
                                                            double *alphaCurrent,
                                                            double *resid,
                                                            double *wCurrent,
                                                            BayesLogit_rpg_hybrid_t pg,
                                                            int warmup)
{
  SamplerIterationResult result;
  double *latent;
  int i;

  result.sizeAccepted = update_likelihood_iteration(s, betaCurrent, alphaCurrent,
                                                     wCurrent, pg, warmup);
  result.wCurrentFresh = 1;
  update_beta_draw_collapsed(s, alphaCurrent, betaCurrent);
  update_mean_without_alpha(s, betaCurrent, resid);
  result.covAccepted = update_covariance_iteration(s, resid, warmup);

  latent = (double*)R_alloc(s->qLatSystem, sizeof(double));
  recover_latent_draw_collapsed(s, resid, latent);
  for(i = 0; i < s->q; i++)
    s->alpha[i] = alphaCurrent[i] = latent[s->qLatTotal + i];
  for(i = 0; i < s->qLatTotal; i++)
    wCurrent[i] = latent[i];
  return result;
}

static int adapt_warmup_block(CovProposalBlock *b, double acceptance, int batch)
{
  b->warmupBatchAcceptHistory[batch] = acceptance;
  b->warmupProposalScaleHistory[batch] = b->sigmaSqM;
  b->warmupBatches = batch + 1;
  if(acceptance < b->warmupNearZero)
    b->sigmaSqM *= 0.10;
  else if(acceptance < 0.10)
    b->sigmaSqM *= 0.25;
  else if(acceptance < b->warmupTargetLower)
    b->sigmaSqM *= 0.50;
  else if(acceptance > 0.65)
    b->sigmaSqM *= 3.00;
  else if(acceptance > b->warmupTargetUpper)
    b->sigmaSqM *= 1.50;
  factor_adaptive_block_proposal_covariance(b);
  return acceptance >= b->warmupTargetLower && acceptance <= b->warmupTargetUpper;
}

static void finish_warmup_block(CovProposalBlock *b)
{
  int k;
  double ratio;

  if(R_FINITE(b->sigmaSqM) && R_FINITE(b->baseSigmaSqM) && b->baseSigmaSqM > 0.0){
    ratio = b->sigmaSqM / b->baseSigmaSqM;
    if(b->scalarMode)
      ratio *= ratio;
    for(k = 0; k < b->dim * b->dim; k++)
      b->Sigma0[k] *= ratio;
    b->sigmaSqM = b->baseSigmaSqM;
    factor_adaptive_block_proposal_covariance(b);
  }
  for(k = 0; k < b->dim * b->dim; k++)
    b->warmupEndingProposalCov[k] = b->proposalCov[k];
}

void run_covariance_warmup(SamplerState *s,
                           double *betaCurrent,
                           double *alphaCurrent,
                           double *resid,
                           double *wCurrent,
                           BayesLogit_rpg_hybrid_t pg,
                           int verbose,
                           int nReport)
{
  CovProposalState *p;
  SamplerIterationResult iteration;
  int b, j, k, block, batchAccepted, allBlocksInTarget, canStop, reportEveryBatches;
  double rhat;
  int sizeAccepted;
  CovProposalBlock *sizeBlock;
  const void *vmax_iter;

  p = &s->covProp;
  sizeBlock = &s->nbSizeProposal;
  if(s->nbSizeEstimate){
    sizeBlock->warmupStartingEta[0] = std::log(s->nbSize);
    sizeBlock->warmupStartingProposalCov[0] = sizeBlock->proposalCov[0];
  }
  if(!p->initialized || (p->nBlocks <= 0 && !s->nbSizeEstimate) || !p->warmupEnabled || p->warmupMaxBatches <= 0){
    if(s->nbSizeEstimate){
      sizeBlock->warmupStoppedReason = 3;
      sizeBlock->warmupEndingEta[0] = std::log(s->nbSize);
      sizeBlock->warmupEndingProposalCov[0] = sizeBlock->proposalCov[0];
    }
    p->warmupStoppedReason = 3;
    for(block = 0; block < p->nBlocks; block++){
      CovProposalBlock *blk = p->blocks + block;
      if(blk->dim > 0){
        pack_cov_block_eta(s, blk, blk->warmupStartingEta);
        pack_cov_block_eta(s, blk, blk->warmupEndingEta);
        for(k = 0; k < blk->dim * blk->dim; k++){
          blk->warmupStartingProposalCov[k] = blk->proposalCov[k];
          blk->warmupEndingProposalCov[k] = blk->proposalCov[k];
        }
      }
    }
    return;
  }

  reportEveryBatches = 0;
  if(nReport > 0)
    reportEveryBatches = (nReport + p->warmupBatchLength - 1) / p->warmupBatchLength;
  if(reportEveryBatches < 1 && nReport > 0)
    reportEveryBatches = 1;

  for(block = 0; block < p->nBlocks; block++){
    CovProposalBlock *blk = p->blocks + block;
    pack_cov_block_eta(s, blk, blk->warmupStartingEta);
    for(k = 0; k < blk->dim * blk->dim; k++)
      blk->warmupStartingProposalCov[k] = blk->proposalCov[k];
  }

  p->warmupStoppedReason = 1;
  for(b = 0; b < p->warmupMaxBatches; b++){
    R_CheckUserInterrupt();
    batchAccepted = 0;
    sizeAccepted = 0;
    for(block = 0; block < p->nBlocks; block++)
      p->blocks[block].batchAccept = 0;

    for(j = 0; j < p->warmupBatchLength; j++){
      R_CheckUserInterrupt();
      vmax_iter = vmaxget();

      if(s->collapseIid)
        iteration = update_collapsed_iid_iteration(s, betaCurrent, alphaCurrent,
                                                   resid, wCurrent, pg, 1);
      else
        iteration = update_separate_iid_iteration(s, betaCurrent, alphaCurrent,
                                                  resid, wCurrent, pg, 1);
      sizeAccepted += iteration.sizeAccepted;
      batchAccepted += iteration.covAccepted;

      vmaxset(vmax_iter);
    }

    rhat = p->nBlocks > 0 ? (double)batchAccepted / (double)(p->warmupBatchLength * p->nBlocks) : NA_REAL;
    p->warmupBatches = b + 1;

    allBlocksInTarget = 1;
    for(block = 0; block < p->nBlocks; block++){
      CovProposalBlock *blk = p->blocks + block;
      double blockRhat = (double)blk->batchAccept / (double)p->warmupBatchLength;
      if(!adapt_warmup_block(blk, blockRhat, b))
        allBlocksInTarget = 0;
    }
    if(s->nbSizeEstimate &&
       !adapt_warmup_block(sizeBlock, (double)sizeAccepted / p->warmupBatchLength, b))
      allBlocksInTarget = 0;

    canStop = allBlocksInTarget && (b + 1 >= p->warmupMinBatches);

    if(verbose && nReport > 0 &&
       (((b + 1) % reportEveryBatches) == 0 ||
        canStop ||
        b == p->warmupMaxBatches - 1)){
      if(p->nBlocks > 0 || !s->nbSizeEstimate)
        Rprintf("warmup batch %d: covariance block acc=%.3f, sigma_sq_m=%.4g\n",
                b + 1, rhat, p->sigmaSqM);
      else
        Rprintf("warmup batch %d\n", b + 1);
      if(s->nbSizeEstimate)
        Rprintf("  size acc=%.3f\n", (double)sizeAccepted / p->warmupBatchLength);
      R_FlushConsole();
    }

    if(canStop){
      p->warmupStoppedReason = 2;
      break;
    }
  }

  for(block = 0; block < p->nBlocks; block++){
    CovProposalBlock *blk = p->blocks + block;
    pack_cov_block_eta(s, blk, blk->warmupEndingEta);
    finish_warmup_block(blk);
  }
  if(s->nbSizeEstimate){
    sizeBlock->warmupStoppedReason = p->warmupStoppedReason;
    sizeBlock->warmupEndingEta[0] = std::log(s->nbSize);
    finish_warmup_block(sizeBlock);
    if(p->nBlocks == 0){
      p->warmupBatches = 0;
      p->warmupStoppedReason = 3;
    }
  }
}

void draw_gaussian_from_precision(double *Prec, double *rhs, int dim, double *draw)
{
  int i, info, inc;
  char lower, trans, diag;
  double *L, *mu, *z;

  if(dim <= 0)
    return;

  lower = 'L';
  trans = 'T';
  diag = 'N';
  inc = 1;

  L = (double*)R_alloc(dim * dim, sizeof(double));
  mu = (double*)R_alloc(dim, sizeof(double));
  z = (double*)R_alloc(dim, sizeof(double));

  for(i = 0; i < dim * dim; i++)
    L[i] = Prec[i];
  for(i = 0; i < dim; i++)
    mu[i] = rhs[i];

  F77_CALL(dpotrf)(&lower, &dim, L, &dim, &info FCONE);
  if(info != 0)
    Rf_error("draw_gaussian_from_precision: dpotrf failed");

  F77_CALL(dpotrs)(&lower, &dim, &inc, L, &dim, mu, &dim, &info FCONE);
  if(info != 0)
    Rf_error("draw_gaussian_from_precision: dpotrs failed");

  for(i = 0; i < dim; i++)
    z[i] = rnorm(0.0, 1.0);

  F77_CALL(dtrsv)(&lower, &trans, &diag, &dim, L, &dim, z, &inc FCONE FCONE FCONE);

  for(i = 0; i < dim; i++)
    draw[i] = mu[i] + z[i];
}

void update_beta_draw_collapsed(SamplerState *s, const double *alpha, double *betaDraw)
{
  int i;
  double *ytilde, *Prec, *rhs;

  if(s->p <= 0)
    return;

  ytilde = (double*)R_alloc(s->n, sizeof(double));
  Prec = (double*)R_alloc(s->p * s->p, sizeof(double));
  rhs = (double*)R_alloc(s->p, sizeof(double));

  if(s->q > 0 && !s->collapseIid)
    sparse_Z_mult(s, alpha, s->Za);
  else
    for(i = 0; i < s->n; i++) s->Za[i] = 0.0;

  for(i = 0; i < s->n; i++)
    ytilde[i] = s->y[i] - s->Za[i];

  form_XtVinvX_and_rhs(s, ytilde, Prec, rhs);
  if(s->betaPriorType == 1){
    for(i = 0; i < s->p; i++){
      Prec[i + s->p * i] += s->betaPriorPrecision[i];
      rhs[i] += s->betaPriorPrecision[i] * s->betaPriorMean[i];
    }
  }
  draw_gaussian_from_precision(Prec, rhs, s->p, betaDraw);
}

void update_alpha_draw_collapsed(SamplerState *s, const double *beta, double *alphaDraw)
{
  int i, j, b, k, g, offset, col, row;
  double *ytilde, *Prec, *rhs;

  if(s->q <= 0)
    return;
  if(s->collapseIid)
    Rf_error("separate IID coefficient draw requires the conditional structured system");

  ytilde = (double*)R_alloc(s->n, sizeof(double));
  Prec = (double*)R_alloc(s->q * s->q, sizeof(double));
  rhs = (double*)R_alloc(s->q, sizeof(double));

  for(i = 0; i < s->n; i++)
    ytilde[i] = s->y[i];
  if(s->p > 0){
    double minusOne = -1.0, one = 1.0;
    int inc = 1;
    F77_CALL(dgemv)("N", &s->n, &s->p,
                    &minusOne, s->X, &s->n, beta, &inc,
                    &one, ytilde, &inc FCONE);
  }

  for(i = 0; i < s->q * s->q; i++)
    Prec[i] = 0.0;
  form_ZtVinvZ_and_rhs(s, ytilde, Prec, rhs);

  for(j = 0; j < s->q; j++)
    if(s->reCovDim[s->reBlockID[j]] == 0)
      Prec[j + s->q * j] += 1.0 / s->sigmaSqRE[s->reBlockID[j]];
  for(b = 0; b < s->nRE; b++){
    k = s->reCovDim[b];
    if(k == 0)
      continue;
    for(g = 0; g < s->reBlockSize[b] / k; g++){
      offset = s->reOffset[b] + g * k;
      for(col = 0; col < k; col++)
        for(row = 0; row < k; row++)
          Prec[offset + row + s->q * (offset + col)] +=
            s->reCovPrecision[b][row + k * col];
    }
  }

  draw_gaussian_from_precision(Prec, rhs, s->q, alphaDraw);
}

SEXP char_vector_from_term_labels(SamplerState *s)
{
  int t;
  SEXP out;

  PROTECT(out = Rf_allocVector(STRSXP, s->nTerms));
  for(t = 0; t < s->nTerms; t++)
    SET_STRING_ELT(out, t, Rf_mkChar(s->terms[t].label));
  UNPROTECT(1);
  return out;
}

SEXP list_B_by_term(SamplerState *s)
{
  int t, j, total_nnbr;
  SEXP out, Bj;

  PROTECT(out = Rf_allocVector(VECSXP, s->nTerms));

  for(t = 0; t < s->nTerms; t++){
    TermState *term = s->terms + t;
    GraphState *g = s->graphs + term->graphIndex;

    if(g->type != GRAPH_NNGP){
      SET_VECTOR_ELT(out, t, R_NilValue);
      continue;
    }

    total_nnbr = g->totalNnbr;
    PROTECT(Bj = Rf_allocVector(REALSXP, total_nnbr));
    for(j = 0; j < total_nnbr; j++)
      REAL(Bj)[j] = term->B[j];
    SET_VECTOR_ELT(out, t, Bj);
    UNPROTECT(1);
  }

  UNPROTECT(1);
  return out;
}

SEXP list_F_by_term(SamplerState *s)
{
  int t, i;
  SEXP out, Fj;

  PROTECT(out = Rf_allocVector(VECSXP, s->nTerms));

  for(t = 0; t < s->nTerms; t++){
    TermState *term = s->terms + t;
    GraphState *g = s->graphs + term->graphIndex;

    if(g->type != GRAPH_NNGP){
      SET_VECTOR_ELT(out, t, R_NilValue);
      continue;
    }

    PROTECT(Fj = Rf_allocVector(REALSXP, term->nNode));
    for(i = 0; i < term->nNode; i++)
      REAL(Fj)[i] = term->F[i];
    SET_VECTOR_ELT(out, t, Fj);
    UNPROTECT(1);
  }

  UNPROTECT(1);
  return out;
}

static void cov_param_label(SamplerState *s, CovProposalState *p, int idx, char *label, size_t labelSize)
{
  int t, j;

  if(p->paramType[idx] == 0){
    std::snprintf(label, labelSize, "log(tau_sq)");
  } else if(p->paramType[idx] == 3){
    t = p->paramTerm[idx];
    if(s->residualModel == 3 && t == 0)
      std::snprintf(label, labelSize, "log(kappa)");
    else if(s->residualModel == 3 && t == 1)
      std::snprintf(label, labelSize, "log(tau0_sq)");
    else
      std::snprintf(label, labelSize, "log(residual_variance_%d)", t + 1);
  } else if(p->paramType[idx] == 1){
    t = p->paramTerm[idx];
    std::snprintf(label, labelSize, "log(%s_sigma_sq)", s->terms[t].name);
  } else if(p->paramType[idx] == 2){
    t = p->paramTerm[idx];
    j = p->paramTheta[idx];
    std::snprintf(label, labelSize, "%s_theta_%d", s->terms[t].name, j + 1);
  } else if(p->paramType[idx] == 4){
    std::snprintf(label, labelSize, "log(%s_sigma_sq)", s->reNames[p->paramTerm[idx]]);
  } else if(p->paramType[idx] == 5){
    t = p->paramTerm[idx];
    j = p->paramTheta[idx];
    if(j % (s->reCovDim[t] + 1) == 0)
      std::snprintf(label, labelSize, "log(%s_chol[%d,%d])", s->reNames[t],
                    j % s->reCovDim[t] + 1, j / s->reCovDim[t] + 1);
    else
      std::snprintf(label, labelSize, "%s_chol[%d,%d]", s->reNames[t],
                    j % s->reCovDim[t] + 1, j / s->reCovDim[t] + 1);
  } else {
    std::snprintf(label, labelSize, "unknown_%d", idx + 1);
  }
}

SEXP build_adaptive_metropolis_blocks_object(SamplerState *s)
{
  CovProposalState *p;
  SEXP out, block_r, names, labels_r, idx_r, eta_r, batch_r, scale_r;
  int b, i, nHistory, gi;
  char label[256], blockName[64];

  p = &s->covProp;
  PROTECT(out = Rf_allocVector(VECSXP, p->nBlocks));

  for(b = 0; b < p->nBlocks; b++){
    CovProposalBlock *blk = p->blocks + b;
    nHistory = blk->batchIndex;
    if(nHistory > blk->maxBatchHistory)
      nHistory = blk->maxBatchHistory;

    PROTECT(block_r = Rf_allocVector(VECSXP, 8));
    PROTECT(names = Rf_allocVector(STRSXP, 8));
    PROTECT(labels_r = Rf_allocVector(STRSXP, blk->dim));
    PROTECT(idx_r = Rf_allocVector(INTSXP, blk->dim));
    PROTECT(eta_r = Rf_allocVector(REALSXP, blk->dim));
    PROTECT(batch_r = Rf_allocVector(REALSXP, nHistory));
    PROTECT(scale_r = Rf_allocVector(REALSXP, nHistory));

    pack_cov_block_eta(s, blk, blk->eta);
    for(i = 0; i < blk->dim; i++){
      gi = blk->paramIndex[i];
      INTEGER(idx_r)[i] = gi + 1;
      REAL(eta_r)[i] = blk->eta[i];
      cov_param_label(s, p, gi, label, sizeof(label));
      SET_STRING_ELT(labels_r, i, Rf_mkChar(label));
    }
    Rf_namesgets(eta_r, labels_r);
    for(i = 0; i < nHistory; i++){
      REAL(batch_r)[i] = blk->batchAcceptHistory[i];
      REAL(scale_r)[i] = blk->proposalScaleHistory[i];
    }

    std::snprintf(blockName, sizeof(blockName), "block_%d", b + 1);
    SET_VECTOR_ELT(block_r, 0, Rf_mkString(blockName));
    SET_VECTOR_ELT(block_r, 1, Rf_ScalarInteger(blk->dim));
    SET_VECTOR_ELT(block_r, 2, idx_r);
    SET_VECTOR_ELT(block_r, 3, labels_r);
    SET_VECTOR_ELT(block_r, 4, eta_r);
    SET_VECTOR_ELT(block_r, 5, Rf_ScalarReal(s->nSamples > 0 ? (double)blk->acceptCount / (double)s->nSamples : NA_REAL));
    SET_VECTOR_ELT(block_r, 6, batch_r);
    SET_VECTOR_ELT(block_r, 7, scale_r);

    SET_STRING_ELT(names, 0, Rf_mkChar("name"));
    SET_STRING_ELT(names, 1, Rf_mkChar("dimension"));
    SET_STRING_ELT(names, 2, Rf_mkChar("parameter_index"));
    SET_STRING_ELT(names, 3, Rf_mkChar("parameter_labels"));
    SET_STRING_ELT(names, 4, Rf_mkChar("current_eta"));
    SET_STRING_ELT(names, 5, Rf_mkChar("acceptance"));
    SET_STRING_ELT(names, 6, Rf_mkChar("batch_acceptance_history"));
    SET_STRING_ELT(names, 7, Rf_mkChar("proposal_scale_history"));
    Rf_namesgets(block_r, names);
    SET_VECTOR_ELT(out, b, block_r);
    UNPROTECT(7);
  }

  UNPROTECT(1);
  return out;
}

SEXP build_adaptive_metropolis_object(SamplerState *s)
{
  CovProposalState *p;
  SEXP out, names, Sigma0_r, proposal_cov_r, eta_r, label_r, batch_accept_r, scale_history_r;
  SEXP warmup_r, warmup_names, warmup_start_eta_r, warmup_end_eta_r;
  SEXP warmup_start_cov_r, warmup_end_cov_r, warmup_batch_accept_r, warmup_scale_history_r, warmup_target_r;
  SEXP blocks_r;
  int i, t, j, nHistory, nWarmupHistory;
  char label[256];
  const char *warmupReason;

  p = &s->covProp;
  nHistory = (p->nBlocks == 1) ? p->blocks[0].batchIndex : 0;
  if(p->nBlocks == 1 && nHistory > p->blocks[0].maxBatchHistory)
    nHistory = p->blocks[0].maxBatchHistory;
  nWarmupHistory = (p->nBlocks == 1) ? p->blocks[0].warmupBatches : 0;
  if(p->nBlocks == 1 && nWarmupHistory > p->blocks[0].warmupMaxBatches)
    nWarmupHistory = p->blocks[0].warmupMaxBatches;

  PROTECT(out = Rf_allocVector(VECSXP, 15));
  PROTECT(names = Rf_allocVector(STRSXP, 15));
  PROTECT(Sigma0_r = Rf_allocMatrix(REALSXP, p->dim, p->dim));
  PROTECT(proposal_cov_r = Rf_allocMatrix(REALSXP, p->dim, p->dim));
  PROTECT(eta_r = Rf_allocVector(REALSXP, p->dim));
  PROTECT(label_r = Rf_allocVector(STRSXP, p->dim));
  PROTECT(batch_accept_r = Rf_allocVector(REALSXP, nHistory));
  PROTECT(scale_history_r = Rf_allocVector(REALSXP, nHistory));
  PROTECT(warmup_r = Rf_allocVector(VECSXP, 16));
  PROTECT(warmup_names = Rf_allocVector(STRSXP, 16));
  PROTECT(warmup_start_eta_r = Rf_allocVector(REALSXP, p->dim));
  PROTECT(warmup_end_eta_r = Rf_allocVector(REALSXP, p->dim));
  PROTECT(warmup_start_cov_r = Rf_allocMatrix(REALSXP, p->dim, p->dim));
  PROTECT(warmup_end_cov_r = Rf_allocMatrix(REALSXP, p->dim, p->dim));
  PROTECT(warmup_batch_accept_r = Rf_allocVector(REALSXP, nWarmupHistory));
  PROTECT(warmup_scale_history_r = Rf_allocVector(REALSXP, nWarmupHistory));
  PROTECT(warmup_target_r = Rf_allocVector(REALSXP, 2));
  PROTECT(blocks_r = build_adaptive_metropolis_blocks_object(s));

  if(p->dim > 0)
    pack_cov_eta(s, p->eta);

  for(i = 0; i < p->dim * p->dim; i++){
    REAL(Sigma0_r)[i] = 0.0;
    REAL(proposal_cov_r)[i] = 0.0;
  }
  for(t = 0; t < p->nBlocks; t++){
    CovProposalBlock *blk = p->blocks + t;
    for(i = 0; i < blk->dim; i++){
      int gi = blk->paramIndex[i];
      for(j = 0; j < blk->dim; j++){
        int gj = blk->paramIndex[j];
        REAL(Sigma0_r)[gi + p->dim * gj] = blk->Sigma0[i + blk->dim * j];
        REAL(proposal_cov_r)[gi + p->dim * gj] = blk->proposalCov[i + blk->dim * j];
      }
    }
  }
  for(i = 0; i < p->dim; i++){
    REAL(eta_r)[i] = p->eta[i];
    REAL(warmup_start_eta_r)[i] = NA_REAL;
    REAL(warmup_end_eta_r)[i] = NA_REAL;
    cov_param_label(s, p, i, label, sizeof(label));
    SET_STRING_ELT(label_r, i, Rf_mkChar(label));
  }
  for(i = 0; i < p->dim * p->dim; i++){
    REAL(warmup_start_cov_r)[i] = 0.0;
    REAL(warmup_end_cov_r)[i] = 0.0;
  }
  for(t = 0; t < p->nBlocks; t++){
    CovProposalBlock *blk = p->blocks + t;
    for(i = 0; i < blk->dim; i++){
      int gi = blk->paramIndex[i];
      REAL(warmup_start_eta_r)[gi] = blk->warmupStartingEta[i];
      REAL(warmup_end_eta_r)[gi] = blk->warmupEndingEta[i];
      for(j = 0; j < blk->dim; j++){
        int gj = blk->paramIndex[j];
        REAL(warmup_start_cov_r)[gi + p->dim * gj] = blk->warmupStartingProposalCov[i + blk->dim * j];
        REAL(warmup_end_cov_r)[gi + p->dim * gj] = blk->warmupEndingProposalCov[i + blk->dim * j];
      }
    }
  }
  for(i = 0; i < nHistory; i++){
    REAL(batch_accept_r)[i] = p->blocks[0].batchAcceptHistory[i];
    REAL(scale_history_r)[i] = p->blocks[0].proposalScaleHistory[i];
  }
  for(i = 0; i < nWarmupHistory; i++){
    REAL(warmup_batch_accept_r)[i] = p->blocks[0].warmupBatchAcceptHistory[i];
    REAL(warmup_scale_history_r)[i] = p->blocks[0].warmupProposalScaleHistory[i];
  }
  REAL(warmup_target_r)[0] = p->warmupTargetLower;
  REAL(warmup_target_r)[1] = p->warmupTargetUpper;
  Rf_namesgets(eta_r, label_r);
  Rf_namesgets(warmup_start_eta_r, label_r);
  Rf_namesgets(warmup_end_eta_r, label_r);

  if(p->warmupStoppedReason == 2)
    warmupReason = "target_reached";
  else if(p->warmupStoppedReason == 3)
    warmupReason = "skipped";
  else
    warmupReason = "max_batches";

  SET_VECTOR_ELT(warmup_r, 0, Rf_ScalarLogical(p->warmupEnabled && p->dim > 0));
  SET_VECTOR_ELT(warmup_r, 1, Rf_ScalarInteger(p->warmupBatchLength));
  SET_VECTOR_ELT(warmup_r, 2, Rf_ScalarInteger(p->warmupMinBatches));
  SET_VECTOR_ELT(warmup_r, 3, Rf_ScalarInteger(p->warmupMaxBatches));
  SET_VECTOR_ELT(warmup_r, 4, Rf_ScalarInteger(p->warmupBatches));
  SET_VECTOR_ELT(warmup_r, 5, Rf_ScalarInteger(p->warmupNAttempted));
  SET_VECTOR_ELT(warmup_r, 6, Rf_ScalarInteger(p->warmupNAccepted));
  SET_VECTOR_ELT(warmup_r, 7, Rf_ScalarReal(p->warmupNAttempted > 0 ? (double)p->warmupNAccepted / (double)p->warmupNAttempted : NA_REAL));
  SET_VECTOR_ELT(warmup_r, 8, warmup_target_r);
  SET_VECTOR_ELT(warmup_r, 9, warmup_batch_accept_r);
  SET_VECTOR_ELT(warmup_r, 10, warmup_scale_history_r);
  SET_VECTOR_ELT(warmup_r, 11, warmup_start_eta_r);
  SET_VECTOR_ELT(warmup_r, 12, warmup_end_eta_r);
  SET_VECTOR_ELT(warmup_r, 13, warmup_start_cov_r);
  SET_VECTOR_ELT(warmup_r, 14, warmup_end_cov_r);
  SET_VECTOR_ELT(warmup_r, 15, Rf_mkString(warmupReason));

  SET_STRING_ELT(warmup_names, 0, Rf_mkChar("enabled"));
  SET_STRING_ELT(warmup_names, 1, Rf_mkChar("batch_length"));
  SET_STRING_ELT(warmup_names, 2, Rf_mkChar("min_batches"));
  SET_STRING_ELT(warmup_names, 3, Rf_mkChar("max_batches"));
  SET_STRING_ELT(warmup_names, 4, Rf_mkChar("n_batches"));
  SET_STRING_ELT(warmup_names, 5, Rf_mkChar("n_attempted"));
  SET_STRING_ELT(warmup_names, 6, Rf_mkChar("n_accepted"));
  SET_STRING_ELT(warmup_names, 7, Rf_mkChar("acceptance"));
  SET_STRING_ELT(warmup_names, 8, Rf_mkChar("target"));
  SET_STRING_ELT(warmup_names, 9, Rf_mkChar("batch_acceptance"));
  SET_STRING_ELT(warmup_names, 10, Rf_mkChar("proposal_scale_history"));
  SET_STRING_ELT(warmup_names, 11, Rf_mkChar("starting_transformed"));
  SET_STRING_ELT(warmup_names, 12, Rf_mkChar("ending_transformed"));
  SET_STRING_ELT(warmup_names, 13, Rf_mkChar("starting_proposal_cov"));
  SET_STRING_ELT(warmup_names, 14, Rf_mkChar("ending_proposal_cov"));
  SET_STRING_ELT(warmup_names, 15, Rf_mkChar("stopped_reason"));
  Rf_namesgets(warmup_r, warmup_names);

  SET_VECTOR_ELT(out, 0, Rf_ScalarInteger(p->dim));
  SET_VECTOR_ELT(out, 1, Rf_ScalarInteger(p->batchLength));
  SET_VECTOR_ELT(out, 2, Rf_ScalarInteger(p->batchIndex));
  SET_VECTOR_ELT(out, 3, Rf_ScalarReal(p->targetAccept));
  SET_VECTOR_ELT(out, 4, Rf_ScalarReal(p->lastBatchAccept));
  SET_VECTOR_ELT(out, 5, Rf_ScalarReal(p->sigmaSqM));
  SET_VECTOR_ELT(out, 6, Sigma0_r);
  SET_VECTOR_ELT(out, 7, proposal_cov_r);
  SET_VECTOR_ELT(out, 8, eta_r);
  SET_VECTOR_ELT(out, 9, label_r);
  SET_VECTOR_ELT(out, 10, Rf_ScalarReal(s->covParamAttempts > 0 ? (double)s->covParamAccept / (double)s->covParamAttempts : NA_REAL));
  SET_VECTOR_ELT(out, 11, batch_accept_r);
  SET_VECTOR_ELT(out, 12, scale_history_r);
  SET_VECTOR_ELT(out, 13, warmup_r);
  SET_VECTOR_ELT(out, 14, blocks_r);

  SET_STRING_ELT(names, 0, Rf_mkChar("dimension"));
  SET_STRING_ELT(names, 1, Rf_mkChar("batch_length"));
  SET_STRING_ELT(names, 2, Rf_mkChar("batch_index"));
  SET_STRING_ELT(names, 3, Rf_mkChar("target_accept"));
  SET_STRING_ELT(names, 4, Rf_mkChar("last_batch_accept"));
  SET_STRING_ELT(names, 5, Rf_mkChar("sigma_sq_m"));
  SET_STRING_ELT(names, 6, Rf_mkChar("Sigma0"));
  SET_STRING_ELT(names, 7, Rf_mkChar("proposal_cov"));
  SET_STRING_ELT(names, 8, Rf_mkChar("current_eta"));
  SET_STRING_ELT(names, 9, Rf_mkChar("parameter_labels"));
  SET_STRING_ELT(names, 10, Rf_mkChar("acceptance"));
  SET_STRING_ELT(names, 11, Rf_mkChar("batch_acceptance_history"));
  SET_STRING_ELT(names, 12, Rf_mkChar("proposal_scale_history"));
  SET_STRING_ELT(names, 13, Rf_mkChar("warmup"));
  SET_STRING_ELT(names, 14, Rf_mkChar("blocks"));
  Rf_namesgets(out, names);

  UNPROTECT(18);
  return out;
}



static SEXP build_size_metropolis_object(SamplerState *s)
{
  CovProposalBlock *b;
  SEXP out, names, batch, scale, warmup, warmupNames, warmupBatch, warmupScale, target;
  int i, nHistory;
  const char *labels[] = {"parameter_labels", "batch_length", "batch_index", "target_accept",
    "current_eta", "proposal_cov", "acceptance", "batch_acceptance_history",
    "proposal_scale_history", "warmup"};
  const char *warmupLabels[] = {"enabled", "n_batches", "n_attempted", "n_accepted",
    "acceptance", "target", "batch_acceptance", "proposal_scale_history",
    "starting_transformed", "ending_transformed", "starting_proposal_cov",
    "ending_proposal_cov", "stopped_reason"};

  b = &s->nbSizeProposal;
  nHistory = b->batchIndex;
  PROTECT(out = Rf_allocVector(VECSXP, 10));
  PROTECT(names = Rf_allocVector(STRSXP, 10));
  PROTECT(batch = Rf_allocVector(REALSXP, nHistory));
  PROTECT(scale = Rf_allocVector(REALSXP, nHistory));
  PROTECT(warmup = Rf_allocVector(VECSXP, 13));
  PROTECT(warmupNames = Rf_allocVector(STRSXP, 13));
  PROTECT(warmupBatch = Rf_allocVector(REALSXP, b->warmupBatches));
  PROTECT(warmupScale = Rf_allocVector(REALSXP, b->warmupBatches));
  PROTECT(target = Rf_allocVector(REALSXP, 2));
  for(i = 0; i < nHistory; i++){
    REAL(batch)[i] = b->batchAcceptHistory[i];
    REAL(scale)[i] = b->proposalScaleHistory[i];
  }
  for(i = 0; i < b->warmupBatches; i++){
    REAL(warmupBatch)[i] = b->warmupBatchAcceptHistory[i];
    REAL(warmupScale)[i] = b->warmupProposalScaleHistory[i];
  }
  REAL(target)[0] = b->warmupTargetLower;
  REAL(target)[1] = b->warmupTargetUpper;
  SET_VECTOR_ELT(warmup, 0, Rf_ScalarLogical(b->warmupEnabled));
  SET_VECTOR_ELT(warmup, 1, Rf_ScalarInteger(b->warmupBatches));
  SET_VECTOR_ELT(warmup, 2, Rf_ScalarInteger(b->warmupNAttempted));
  SET_VECTOR_ELT(warmup, 3, Rf_ScalarInteger(b->warmupNAccepted));
  SET_VECTOR_ELT(warmup, 4, Rf_ScalarReal(b->warmupNAttempted > 0 ?
    (double)b->warmupNAccepted / b->warmupNAttempted : NA_REAL));
  SET_VECTOR_ELT(warmup, 5, target);
  SET_VECTOR_ELT(warmup, 6, warmupBatch);
  SET_VECTOR_ELT(warmup, 7, warmupScale);
  SET_VECTOR_ELT(warmup, 8, Rf_ScalarReal(b->warmupStartingEta[0]));
  SET_VECTOR_ELT(warmup, 9, Rf_ScalarReal(b->warmupEndingEta[0]));
  SET_VECTOR_ELT(warmup, 10, Rf_ScalarReal(b->warmupStartingProposalCov[0]));
  SET_VECTOR_ELT(warmup, 11, Rf_ScalarReal(b->warmupEndingProposalCov[0]));
  SET_VECTOR_ELT(warmup, 12, Rf_mkString(b->warmupStoppedReason == 3 ? "skipped" :
    (b->warmupStoppedReason == 2 ? "target_reached" : "max_batches")));
  for(i = 0; i < 13; i++)
    SET_STRING_ELT(warmupNames, i, Rf_mkChar(warmupLabels[i]));
  Rf_namesgets(warmup, warmupNames);
  SET_VECTOR_ELT(out, 0, Rf_mkString("log(size)"));
  SET_VECTOR_ELT(out, 1, Rf_ScalarInteger(b->batchLength));
  SET_VECTOR_ELT(out, 2, Rf_ScalarInteger(b->batchIndex));
  SET_VECTOR_ELT(out, 3, Rf_ScalarReal(b->targetAccept));
  SET_VECTOR_ELT(out, 4, Rf_ScalarReal(std::log(s->nbSize)));
  SET_VECTOR_ELT(out, 5, Rf_ScalarReal(b->proposalCov[0]));
  SET_VECTOR_ELT(out, 6, Rf_ScalarReal(s->nSamples > 0 ? (double)b->acceptCount / s->nSamples : NA_REAL));
  SET_VECTOR_ELT(out, 7, batch);
  SET_VECTOR_ELT(out, 8, scale);
  SET_VECTOR_ELT(out, 9, warmup);
  for(i = 0; i < 10; i++)
    SET_STRING_ELT(names, i, Rf_mkChar(labels[i]));
  Rf_namesgets(out, names);
  UNPROTECT(9);
  return out;
}

extern "C" {

  SEXP stLMM_collapsed_sampler(SEXP backend_r)
  {
    // Read the monotonic wall clock only at phase boundaries. No timer enters
    // an MCMC update or changes the random-number stream.
    typedef std::chrono::steady_clock SamplerClock;
    const SamplerClock::time_point initialization_start = SamplerClock::now();
    SamplerState s;
    SamplerIterationResult iteration;
    SEXP out, outNames, term_description_r, iid_cov_samples_r;
    SEXP phase_seconds_r, phase_names_r;
    SEXP beta_samples_r, alpha_samples_r, tau_samples_r, residual_variance_samples_r, sigma_re_samples_r, theta_samples_r, process_sigma_samples_r;
    SEXP term_labels_r, theta_accept_r, covariance_acceptance_r, w_samples_r, recover_iter_r, adaptive_metropolis_r;
    double *iidCovSamples, *sizeSamples;
    double *betaSamples, *alphaSamples, *tauSamples, *residualVarianceSamples, *sigmaSqRESamples, *thetaSamples, *processSigmaSamples;
    double *betaCurrent, *alphaCurrent, *resid, *wSamples, *wCurrent;
    int *recoverIter;
    int i, j, iter, iter1, b, thetaOffset, describeTerms, nReport, reportLen, tauAcceptReport, doVerbose;
    int recoverIdx, covColumn, covDim, covCol, covRow;
    size_t z;
    const void *vmax_iter;
    BayesLogit_rpg_hybrid_t pg = BayesLogit_rpg_hybrid();

    SEXP collapse_r = getListElement(backend_r, "collapse_iid");
    init_sampler_state(&s, backend_r, collapse_r == R_NilValue ? 0 :
                        as_flag_scalar(collapse_r, "backend$collapse_iid"));
    
    for(i = 0; i < s.nTerms; i++){
      GraphState *g = s.graphs + s.terms[i].graphIndex;
      TermState *term = s.terms + i;
      
      if(g->type == GRAPH_NNGP)
	update_nngp_BF(&s, g, term);
      else if(g->type == GRAPH_GP)
	update_gp_Q(&s, g, term);
    }
    
    build_M_lat_pattern_sparse(&s);
    s.M_lat_sym = M_cholmod_analyze(s.M_lat, &s.cm);
    if(s.M_lat_sym == NULL || s.cm.status != CHOLMOD_OK)
      Rf_error("cholmod analyze failed for latent M");
    
    s.M_lat_fac = M_cholmod_copy_factor(s.M_lat_sym, &s.cm);
    if(s.M_lat_fac == NULL || s.cm.status != CHOLMOD_OK)
      Rf_error("cholmod copy_factor failed for latent M");
    
    build_M_lat_index_caches(&s);

    assemble_M_lat_numeric(&s);
    if(!M_cholmod_factorize(s.M_lat, s.M_lat_fac, &s.cm) || s.cm.status != CHOLMOD_OK)
      Rf_error("cholmod factorize failed for latent M");

    init_cov_proposal_state(&s);
    
    describeTerms = as_flag_scalar(getListElement(backend_r, "describe_terms"), "backend$describe_terms");
    nReport = as_nonneg_int_scalar(getListElement(backend_r, "n_report"), "backend$n_report");
    doVerbose = as_flag_scalar(getListElement(backend_r, "verbose"), "backend$verbose");
    
    tauAcceptReport = 0;
    reportLen = 0;
    
    PROTECT(term_description_r = build_term_description_object(&s, backend_r, 0));
    if(describeTerms)
      print_term_description(term_description_r);
    
    sizeSamples = s.nbSizeEstimate ? (double*)R_alloc(s.nSamples, sizeof(double)) : NULL;
    iidCovSamples = (double*)R_alloc((size_t)s.nSamples * s.nRECovElements, sizeof(double));
    betaSamples = (double*)R_alloc((s.nSamples > 0 ? s.nSamples : 1) * (s.p > 0 ? s.p : 1), sizeof(double));
    alphaSamples = (double*)R_alloc((s.nSamples > 0 ? s.nSamples : 1) * (s.q > 0 ? s.q : 1), sizeof(double));
    sigmaSqRESamples = (double*)R_alloc((s.nSamples > 0 ? s.nSamples : 1) * (s.nRE > 0 ? s.nRE : 1), sizeof(double));
    thetaSamples = (double*)R_alloc((s.nSamples > 0 ? s.nSamples : 1) * (s.nThetaTotal > 0 ? s.nThetaTotal : 1), sizeof(double));
    processSigmaSamples = (double*)R_alloc((s.nSamples > 0 ? s.nSamples : 1) * (s.nTerms > 0 ? s.nTerms : 1), sizeof(double));
    tauSamples = (double*)R_alloc(s.nSamples > 0 ? s.nSamples : 1, sizeof(double));
    residualVarianceSamples = (double*)R_alloc((s.nSamples > 0 ? s.nSamples : 1) *
                                               (s.nResidualGroup > 0 ? s.nResidualGroup : 1),
                                               sizeof(double));
    betaCurrent = (double*)R_alloc(s.p > 0 ? s.p : 1, sizeof(double));
    alphaCurrent = (double*)R_alloc(s.q > 0 ? s.q : 1, sizeof(double));
    resid = (double*)R_alloc(s.n > 0 ? s.n : 1, sizeof(double));
    wSamples = (double*)R_alloc((s.nRecover > 0 ? s.nRecover : 1) *
                                (s.qLatTotal > 0 ? s.qLatTotal : 1), sizeof(double));
    wCurrent = (double*)R_alloc(s.qLatTotal > 0 ? s.qLatTotal : 1, sizeof(double));
    recoverIter = (int*)R_alloc(s.nRecover > 0 ? s.nRecover : 1, sizeof(int));
    
    for(i = 0; i < (s.nSamples > 0 ? s.nSamples : 1) * (s.p > 0 ? s.p : 1); i++) betaSamples[i] = 0.0;
    for(i = 0; i < (s.nSamples > 0 ? s.nSamples : 1) * (s.q > 0 ? s.q : 1); i++) alphaSamples[i] = 0.0;
    for(i = 0; i < (s.nSamples > 0 ? s.nSamples : 1) * (s.nRE > 0 ? s.nRE : 1); i++) sigmaSqRESamples[i] = 0.0;
    for(i = 0; i < (s.nSamples > 0 ? s.nSamples : 1) * (s.nThetaTotal > 0 ? s.nThetaTotal : 1); i++) thetaSamples[i] = 0.0;
    for(i = 0; i < (s.nSamples > 0 ? s.nSamples : 1) * (s.nTerms > 0 ? s.nTerms : 1); i++) processSigmaSamples[i] = 0.0;
    for(i = 0; i < (s.nSamples > 0 ? s.nSamples : 1); i++) tauSamples[i] = s.tauSq;
    for(i = 0; i < (s.nSamples > 0 ? s.nSamples : 1) * (s.nResidualGroup > 0 ? s.nResidualGroup : 1); i++)
      residualVarianceSamples[i] = NA_REAL;
    for(i = 0; i < (s.nRecover > 0 ? s.nRecover : 1) * (s.qLatTotal > 0 ? s.qLatTotal : 1); i++) wSamples[i] = 0.0;
    for(i = 0; i < (s.qLatTotal > 0 ? s.qLatTotal : 1); i++) wCurrent[i] = 0.0;
    for(i = 0; i < (s.nRecover > 0 ? s.nRecover : 1); i++) recoverIter[i] = 0;
    
    for(i = 0; i < s.p; i++) betaCurrent[i] = 0.0;
    if(is_augmented_likelihood(&s)){
      SEXP beta_start_r = getListElement(backend_r, "beta_starting");
      if(beta_start_r != R_NilValue){
        if(TYPEOF(beta_start_r) != REALSXP || LENGTH(beta_start_r) != s.p)
          Rf_error("backend$beta_starting malformed");
        for(i = 0; i < s.p; i++)
          betaCurrent[i] = REAL(beta_start_r)[i];
      }
    }
    for(i = 0; i < s.q; i++) alphaCurrent[i] = s.alpha[i];
    recoverIdx = 0;

    if(doVerbose){
      Rprintf("Starting MCMC sampling: n_samples=%d, process_terms=%d, covariance_blocks=%d\n",
              s.nSamples, s.nTerms, s.covProp.nBlocks);
      R_FlushConsole();
    }
    
    GetRNGstate();
    const SamplerClock::time_point warmup_start = SamplerClock::now();
    run_covariance_warmup(&s, betaCurrent, alphaCurrent, resid, wCurrent, pg, doVerbose, nReport);

    const SamplerClock::time_point sampling_start = SamplerClock::now();
    for(iter = 0; iter < s.nSamples; iter++){
      R_CheckUserInterrupt();
      /*
	Many sampler subroutines use R_alloc() for temporary work arrays.
	Without resetting the transient allocation stack, those allocations
	accumulate until the .Call returns, which shows up as iteration-by-
	iteration RAM growth. Rewinding to vmax_iter at the end of each loop
	iteration keeps the transient workspace bounded while preserving all
	persistent objects allocated before the loop.
      */
      vmax_iter = vmaxget();
      if(s.collapseIid)
        iteration = update_collapsed_iid_iteration(&s, betaCurrent, alphaCurrent,
                                                   resid, wCurrent, pg, 0);
      else
        iteration = update_separate_iid_iteration(&s, betaCurrent, alphaCurrent,
                                                  resid, wCurrent, pg, 0);
      tauAcceptReport += iteration.covAccepted;
      
      reportLen++;
      if(nReport > 0 && doVerbose && (reportLen >= nReport || iter == s.nSamples - 1)){
	Rprintf("iter %d", iter + 1);
        if(s.covProp.nBlocks > 0 || !s.nbSizeEstimate)
          Rprintf(": covariance block acc=%.3f", s.covProp.nBlocks > 0 ? (double) tauAcceptReport / (double)(reportLen * s.covProp.nBlocks) : NA_REAL);
        if(s.nbSizeEstimate)
          Rprintf(": size=%.4g, size acc=%.3f", s.nbSize,
                  (double)s.nbSizeProposal.acceptCount / (iter + 1));
	Rprintf("\n");
	tauAcceptReport = 0;
	reportLen = 0;
      }
      
      for(i = 0; i < s.p; i++)
	betaSamples[iter + s.nSamples * i] = betaCurrent[i];
      for(i = 0; i < s.q; i++)
	alphaSamples[iter + s.nSamples * i] = alphaCurrent[i];
      for(b = 0; b < s.nRE; b++)
	sigmaSqRESamples[iter + s.nSamples * b] = s.sigmaSqRE[b];
      covColumn = 0;
      for(b = 0; b < s.nRE; b++){
        covDim = s.reCovDim[b];
        for(covCol = 0; covCol < covDim; covCol++)
          for(covRow = covCol; covRow < covDim; covRow++)
            iidCovSamples[iter + s.nSamples * covColumn++] =
              s.reCovariance[b][covRow + covDim * covCol];
      }
      for(i = 0; i < s.nTerms; i++)
	processSigmaSamples[iter + s.nSamples * i] = s.terms[i].sigmaSq;
      thetaOffset = 0;
      for(i = 0; i < s.nTerms; i++){
	for(j = 0; j < s.terms[i].thetaDim; j++)
	  thetaSamples[iter + s.nSamples * (thetaOffset + j)] = s.terms[i].theta[j];
	thetaOffset += s.terms[i].thetaDim;
      }
      tauSamples[iter] = s.tauSq;
      if(s.nbSizeEstimate)
        sizeSamples[iter] = s.nbSize;
      for(i = 0; i < s.nResidualGroup; i++)
        residualVarianceSamples[iter + s.nSamples * i] = s.residualVariance[i];

      iter1 = iter + 1;
      if(should_recover_iteration(&s, iter1)){
        /*
          Augmented non-Gaussian fits need a fresh latent-process draw each
          iteration for the next latent-data update. On recovery iterations, reuse
          that draw instead of drawing another iid conditional sample only
          to store it.
        */
        if(!iteration.wCurrentFresh)
          recover_w_draw_collapsed(&s, resid, wCurrent);

        if(recoverIdx >= s.nRecover)
          Rf_error("internal error: recoverIdx exceeds nRecover");

        for(i = 0; i < s.qLatTotal; i++)
          wSamples[recoverIdx + s.nRecover * i] = wCurrent[i];

        recoverIter[recoverIdx] = iter1;
        recoverIdx++;
      }
      vmaxset(vmax_iter);
    }
    const SamplerClock::time_point sampling_end = SamplerClock::now();
    if(recoverIdx != s.nRecover)
      Rf_error("internal error: recovered sample count mismatch");
    PutRNGstate();
    
    update_resid_for_current_fixed_random(&s, betaCurrent, alphaCurrent, resid);
    
    UNPROTECT(1);
    PROTECT(term_description_r = build_term_description_object(&s, backend_r, s.nSamples));
    
    PROTECT(out = Rf_allocVector(VECSXP, s.nbSizeEstimate ? 18 : 16));
    PROTECT(outNames = Rf_allocVector(STRSXP, s.nbSizeEstimate ? 18 : 16));
    PROTECT(iid_cov_samples_r = Rf_allocMatrix(REALSXP, s.nSamples, s.nRECovElements));
    for(z = 0; z < (size_t)s.nSamples * s.nRECovElements; z++)
      REAL(iid_cov_samples_r)[z] = iidCovSamples[z];
    PROTECT(beta_samples_r = Rf_allocMatrix(REALSXP, s.nSamples, s.p > 0 ? s.p : 0));
    PROTECT(alpha_samples_r = Rf_allocMatrix(REALSXP, s.nSamples, s.q > 0 ? s.q : 0));
    PROTECT(sigma_re_samples_r = Rf_allocMatrix(REALSXP, s.nSamples, s.nRE > 0 ? s.nRE : 0));
    PROTECT(theta_samples_r = Rf_allocMatrix(REALSXP, s.nSamples, s.nThetaTotal > 0 ? s.nThetaTotal : 0));
    PROTECT(process_sigma_samples_r = Rf_allocMatrix(REALSXP, s.nSamples, s.nTerms > 0 ? s.nTerms : 0));
    PROTECT(tau_samples_r = Rf_allocVector(REALSXP, s.nSamples));
    PROTECT(residual_variance_samples_r = Rf_allocMatrix(REALSXP, s.nSamples, s.nResidualGroup > 0 ? s.nResidualGroup : 0));
    PROTECT(term_labels_r = char_vector_from_term_labels(&s));
    PROTECT(theta_accept_r = Rf_allocVector(REALSXP, s.nTerms));
    PROTECT(covariance_acceptance_r = Rf_ScalarReal(s.covParamAttempts > 0 ? (double)s.covParamAccept / (double)s.covParamAttempts : NA_REAL));
    PROTECT(w_samples_r = Rf_allocMatrix(REALSXP, s.nRecover, s.qLatTotal > 0 ? s.qLatTotal : 0));
    PROTECT(recover_iter_r = Rf_allocVector(INTSXP, s.nRecover));
    PROTECT(adaptive_metropolis_r = build_adaptive_metropolis_object(&s));
    PROTECT(phase_seconds_r = Rf_allocVector(REALSXP, 3));
    PROTECT(phase_names_r = Rf_allocVector(STRSXP, 3));
    REAL(phase_seconds_r)[0] = std::chrono::duration<double>(warmup_start - initialization_start).count();
    REAL(phase_seconds_r)[1] = std::chrono::duration<double>(sampling_start - warmup_start).count();
    REAL(phase_seconds_r)[2] = std::chrono::duration<double>(sampling_end - sampling_start).count();
    SET_STRING_ELT(phase_names_r, 0, Rf_mkChar("initialization"));
    SET_STRING_ELT(phase_names_r, 1, Rf_mkChar("warmup"));
    SET_STRING_ELT(phase_names_r, 2, Rf_mkChar("sampling"));
    Rf_namesgets(phase_seconds_r, phase_names_r);
    
    for(i = 0; i < s.nSamples * s.p; i++) REAL(beta_samples_r)[i] = betaSamples[i];
    for(i = 0; i < s.nSamples * s.q; i++) REAL(alpha_samples_r)[i] = alphaSamples[i];
    for(i = 0; i < s.nSamples * s.nRE; i++) REAL(sigma_re_samples_r)[i] = sigmaSqRESamples[i];
    for(i = 0; i < s.nSamples * s.nThetaTotal; i++) REAL(theta_samples_r)[i] = thetaSamples[i];
    for(i = 0; i < s.nSamples * s.nTerms; i++) REAL(process_sigma_samples_r)[i] = processSigmaSamples[i];
    for(i = 0; i < s.nSamples; i++) REAL(tau_samples_r)[i] = tauSamples[i];
    for(i = 0; i < s.nSamples * s.nResidualGroup; i++)
      REAL(residual_variance_samples_r)[i] = residualVarianceSamples[i];
    for(i = 0; i < s.nRecover * s.qLatTotal; i++) REAL(w_samples_r)[i] = wSamples[i];
    for(i = 0; i < s.nRecover; i++) INTEGER(recover_iter_r)[i] = recoverIter[i];
    for(i = 0; i < s.nTerms; i++)
      REAL(theta_accept_r)[i] = NA_REAL;
    
    SET_VECTOR_ELT(out, 0, beta_samples_r);
    SET_VECTOR_ELT(out, 1, alpha_samples_r);
    SET_VECTOR_ELT(out, 2, sigma_re_samples_r);
    SET_VECTOR_ELT(out, 3, theta_samples_r);
    SET_VECTOR_ELT(out, 4, process_sigma_samples_r);
    SET_VECTOR_ELT(out, 5, tau_samples_r);
    SET_VECTOR_ELT(out, 6, residual_variance_samples_r);
    SET_VECTOR_ELT(out, 7, term_labels_r);
    SET_VECTOR_ELT(out, 8, theta_accept_r);
    SET_VECTOR_ELT(out, 9, term_description_r);
    SET_VECTOR_ELT(out, 10, covariance_acceptance_r);
    SET_VECTOR_ELT(out, 11, w_samples_r);
    SET_VECTOR_ELT(out, 12, recover_iter_r);
    SET_VECTOR_ELT(out, 13, adaptive_metropolis_r);
    SET_VECTOR_ELT(out, 14, phase_seconds_r);
    SET_VECTOR_ELT(out, 15, iid_cov_samples_r);
    
    SET_STRING_ELT(outNames, 0, Rf_mkChar("beta_samples"));
    SET_STRING_ELT(outNames, 1, Rf_mkChar("alpha_samples"));
    SET_STRING_ELT(outNames, 2, Rf_mkChar("sigma_sq_re_samples"));
    SET_STRING_ELT(outNames, 3, Rf_mkChar("theta_samples"));
    SET_STRING_ELT(outNames, 4, Rf_mkChar("sigma_sq_samples"));
    SET_STRING_ELT(outNames, 5, Rf_mkChar("tau_sq_samples"));
    SET_STRING_ELT(outNames, 6, Rf_mkChar("residual_variance_samples"));
    SET_STRING_ELT(outNames, 7, Rf_mkChar("term_labels"));
    SET_STRING_ELT(outNames, 8, Rf_mkChar("term_param_accept"));
    SET_STRING_ELT(outNames, 9, Rf_mkChar("term_description"));
    SET_STRING_ELT(outNames, 10, Rf_mkChar("covariance_acceptance"));
    SET_STRING_ELT(outNames, 11, Rf_mkChar("w_samples"));
    SET_STRING_ELT(outNames, 12, Rf_mkChar("recover_iter"));
    SET_STRING_ELT(outNames, 13, Rf_mkChar("adaptive_metropolis"));
    SET_STRING_ELT(outNames, 14, Rf_mkChar("phase_seconds"));
    SET_STRING_ELT(outNames, 15, Rf_mkChar("iid_covariance_samples"));
    if(s.nbSizeEstimate){
      SEXP size_samples_r, size_metropolis_r;
      PROTECT(size_samples_r = Rf_allocVector(REALSXP, s.nSamples));
      PROTECT(size_metropolis_r = build_size_metropolis_object(&s));
      for(i = 0; i < s.nSamples; i++)
        REAL(size_samples_r)[i] = sizeSamples[i];
      SET_VECTOR_ELT(out, 16, size_samples_r);
      SET_VECTOR_ELT(out, 17, size_metropolis_r);
      SET_STRING_ELT(outNames, 16, Rf_mkChar("size_samples"));
      SET_STRING_ELT(outNames, 17, Rf_mkChar("size_metropolis"));
      UNPROTECT(2);
    }
    Rf_namesgets(out, outNames);

    UNPROTECT(19);
    free_sampler_state(&s);
    return out;
  }


  
} // extern "C"

int term_prior_block_nnz(TermState *term, GraphState *g)
{
  PatternCol *cols;
  int q, row, start, m, supportSize, a, b, col, k, nnz, t;
  int *nodes;

  if(g->type == GRAPH_AR1){
    if(term->nNode <= 0)
      return 0;
    return term->nNode + (term->nNode > 1 ? term->nNode - 1 : 0);
  }

  if(g->type == GRAPH_GP){
    if(term->nNode <= 0)
      return 0;
    return term->nNode * (term->nNode + 1) / 2;
  }

  if(g->type == GRAPH_CAR){
    if(term->nNode <= 0)
      return 0;
    return term->nNode + g->nEdge;
  }

  if(g->type == GRAPH_CAR_TIME){
    if(term->nNode <= 0)
      return 0;
    return g->nSpace * (2 * g->nTime - 1) + g->nEdge * (3 * g->nTime - 2);
  }

  if(g->type == GRAPH_DAGAR_TIME){
    if(term->nNode <= 0)
      return 0;
    q = term->nNode;
    cols = R_Calloc(q > 0 ? q : 1, PatternCol);
    for(col = 0; col < q; col++){
      cols[col].rows = NULL;
      cols[col].len = 0;
      cols[col].cap = 0;
    }

    nodes = (int*)R_alloc(1 + (g->totalParent > 0 ? g->totalParent : 1), sizeof(int));
    for(row = 0; row < g->nSpace; row++){
      start = g->parentStart[row];
      m = g->parentCount[row];
      supportSize = m + 1;

      nodes[0] = row;
      for(a = 0; a < m; a++)
        nodes[a + 1] = g->parentIndx[start + a];

      for(a = 0; a < supportSize; a++){
        for(b = 0; b <= a; b++){
          int s1 = nodes[b];
          int s2 = nodes[a];
          for(t = 0; t < g->nTime; t++)
            pattern_append_lower(cols, sampler_space_time_node(s2, t, g->nTime), sampler_space_time_node(s1, t, g->nTime));
          for(t = 0; t < g->nTime - 1; t++){
            if(s1 == s2){
              pattern_append_lower(cols, sampler_space_time_node(s1, t + 1, g->nTime), sampler_space_time_node(s1, t, g->nTime));
            } else {
              pattern_append_lower(cols, sampler_space_time_node(s2, t + 1, g->nTime), sampler_space_time_node(s1, t, g->nTime));
              pattern_append_lower(cols, sampler_space_time_node(s2, t, g->nTime), sampler_space_time_node(s1, t + 1, g->nTime));
            }
          }
        }
      }
    }

    nnz = 0;
    for(col = 0; col < q; col++){
      if(cols[col].len > 1)
        qsort(cols[col].rows, (size_t)cols[col].len, sizeof(int), cmp_int_asc);

      if(cols[col].len > 0){
        int u;
        u = 1;
        for(k = 1; k < cols[col].len; k++){
          if(cols[col].rows[k] != cols[col].rows[u - 1]){
            cols[col].rows[u] = cols[col].rows[k];
            u++;
          }
        }
        cols[col].len = u;
      }
      nnz += cols[col].len;
    }

    for(col = 0; col < q; col++){
      if(cols[col].rows != NULL)
        R_Free(cols[col].rows);
    }
    R_Free(cols);
    return nnz;
  }

  if(g->type == GRAPH_DAGAR){
    if(term->nNode <= 0)
      return 0;
    q = term->nNode;
    cols = R_Calloc(q > 0 ? q : 1, PatternCol);
    for(col = 0; col < q; col++){
      cols[col].rows = NULL;
      cols[col].len = 0;
      cols[col].cap = 0;
    }

    nodes = (int*)R_alloc(1 + (g->totalParent > 0 ? g->totalParent : 1), sizeof(int));
    for(row = 0; row < g->nNode; row++){
      start = g->parentStart[row];
      m = g->parentCount[row];
      supportSize = m + 1;

      nodes[0] = row;
      for(a = 0; a < m; a++)
        nodes[a + 1] = g->parentIndx[start + a];

      for(a = 0; a < supportSize; a++){
        for(b = 0; b <= a; b++){
          if(nodes[a] >= nodes[b])
            pattern_col_append(cols, nodes[b], nodes[a]);
          else
            pattern_col_append(cols, nodes[a], nodes[b]);
        }
      }
    }

    nnz = 0;
    for(col = 0; col < q; col++){
      if(cols[col].len > 1)
        qsort(cols[col].rows, (size_t)cols[col].len, sizeof(int), cmp_int_asc);

      if(cols[col].len > 0){
        int u;
        u = 1;
        for(k = 1; k < cols[col].len; k++){
          if(cols[col].rows[k] != cols[col].rows[u - 1]){
            cols[col].rows[u] = cols[col].rows[k];
            u++;
          }
        }
        cols[col].len = u;
      }
      nnz += cols[col].len;
    }

    for(col = 0; col < q; col++){
      if(cols[col].rows != NULL)
        R_Free(cols[col].rows);
    }
    R_Free(cols);
    return nnz;
  }

  if(g->type != GRAPH_NNGP)
    return NA_INTEGER;

  q = term->nNode;
  cols = R_Calloc(q > 0 ? q : 1, PatternCol);
  for(col = 0; col < q; col++){
    cols[col].rows = NULL;
    cols[col].len = 0;
    cols[col].cap = 0;
  }

  nodes = (int*)R_alloc(1 + (g->totalNnbr > 0 ? g->totalNnbr : 1), sizeof(int));

  for(row = 0; row < g->nNode; row++){
    start = g->nnStart[row];
    m = g->nnCount[row];
    supportSize = m + 1;

    nodes[0] = row;
    for(a = 0; a < m; a++)
      nodes[a + 1] = g->nnIndx[start + a];

    for(a = 0; a < supportSize; a++){
      for(b = 0; b <= a; b++){
        if(nodes[a] >= nodes[b])
          pattern_col_append(cols, nodes[b], nodes[a]);
        else
          pattern_col_append(cols, nodes[a], nodes[b]);
      }
    }
  }

  nnz = 0;
  for(col = 0; col < q; col++){
    if(cols[col].len > 1)
      qsort(cols[col].rows, (size_t)cols[col].len, sizeof(int), cmp_int_asc);

    if(cols[col].len > 0){
      int u;
      u = 1;
      for(k = 1; k < cols[col].len; k++){
        if(cols[col].rows[k] != cols[col].rows[u - 1]){
          cols[col].rows[u] = cols[col].rows[k];
          u++;
        }
      }
      cols[col].len = u;
    }

    nnz += cols[col].len;
  }

  for(col = 0; col < q; col++){
    if(cols[col].rows != NULL)
      R_Free(cols[col].rows);
  }
  R_Free(cols);

  return nnz;
}

void print_rule(char ch, int n)
{
  int k;
  for(k = 0; k < n; k++)
    Rprintf("%c", ch);
  Rprintf("\n");
}

void print_term_description(SEXP td_r)
{
  SEXP global_r, fixed_r, random_r, process_terms_r, re_terms_r;
  int i, j;

  global_r = getListElement(td_r, "global");
  fixed_r = getListElement(td_r, "fixed_effects");
  random_r = getListElement(td_r, "random_effects");
  process_terms_r = getListElement(td_r, "process_terms");
  re_terms_r = getListElement(random_r, "terms");

  print_rule('=', 72);
  Rprintf("stLMM term description\n");
  print_rule('=', 72);

  Rprintf("Model summary\n");
  print_rule('-', 72);
  Rprintf("  formula: %s\n", CHAR(STRING_ELT(getListElement(global_r, "formula"), 0)));
  Rprintf("  n: %d\n", INTEGER(getListElement(global_r, "n"))[0]);
  Rprintf("  p: %d\n", INTEGER(getListElement(global_r, "p"))[0]);
  Rprintf("  q: %d\n", INTEGER(getListElement(global_r, "q"))[0]);
  Rprintf("  qLatTotal: %d\n", INTEGER(getListElement(global_r, "qLatTotal"))[0]);
  Rprintf("  process terms: %d\n", INTEGER(getListElement(global_r, "n_process_terms"))[0]);
  Rprintf("  explicit random-effect terms: %d\n", INTEGER(getListElement(random_r, "n_terms"))[0]);
  Rprintf("  graphs: %d\n", INTEGER(getListElement(global_r, "n_graphs"))[0]);
  Rprintf("  M dimensions: %d x %d\n",
          INTEGER(getListElement(global_r, "M_dim"))[0],
          INTEGER(getListElement(global_r, "M_dim"))[1]);
  Rprintf("  M nnz: %d\n", INTEGER(getListElement(global_r, "M_nnz"))[0]);
  Rprintf("  CHOLMOD ordering: %s\n",
          CHAR(STRING_ELT(getListElement(global_r, "cholmod_ordering"), 0)));
  Rprintf("  CHOLMOD fill ratio: %.3f\n",
          REAL(getListElement(global_r, "cholmod_fill_ratio"))[0]);
  Rprintf("  CHOLMOD lnz: %.6g\n",
          REAL(getListElement(global_r, "cholmod_lnz"))[0]);
  Rprintf("  CHOLMOD flops: %.6g\n",
          REAL(getListElement(global_r, "cholmod_flops"))[0]);
  Rprintf("  factorization status: %s\n",
          CHAR(STRING_ELT(getListElement(global_r, "factorization_status"), 0)));

  Rprintf("\n");
  Rprintf("Fixed effects (X)\n");
  print_rule('-', 72);
  Rprintf("  columns (p): %d\n", INTEGER(getListElement(fixed_r, "p"))[0]);
  {
    SEXP xnames = getListElement(fixed_r, "names");
    Rprintf("  names: ");
    for(i = 0; i < LENGTH(xnames); i++)
      Rprintf("%s%s",
              CHAR(STRING_ELT(xnames, i)),
              (i == LENGTH(xnames) - 1) ? "" : ", ");
    Rprintf("\n");
  }

  Rprintf("\n");
  Rprintf("Explicit random effects (Z)\n");
  print_rule('-', 72);
  if(re_terms_r == R_NilValue || LENGTH(re_terms_r) == 0){
    Rprintf("  none\n");
  } else {
    for(i = 0; i < LENGTH(re_terms_r); i++){
      SEXP re_i = VECTOR_ELT(re_terms_r, i);
      SEXP sigma_r = getListElement(re_i, "sigma_sq");
      SEXP coeffs_r = getListElement(re_i, "coefficients");

      if(i > 0)
        print_rule('.', 60);

      Rprintf("  [%d] %s\n", i + 1, CHAR(STRING_ELT(getListElement(re_i, "name"), 0)));
      Rprintf("    label: %s\n", CHAR(STRING_ELT(getListElement(re_i, "label"), 0)));
      Rprintf("    grouping factor: %s\n",
              CHAR(STRING_ELT(getListElement(re_i, "grouping_factor"), 0)));
      Rprintf("    coefficients: ");
      for(j = 0; j < LENGTH(coeffs_r); j++)
        Rprintf("%s%s",
                CHAR(STRING_ELT(coeffs_r, j)),
                (j == LENGTH(coeffs_r) - 1) ? "" : ", ");
      Rprintf("\n");
      Rprintf("    levels: %d\n", INTEGER(getListElement(re_i, "n_levels"))[0]);
      Rprintf("    q contribution: %d\n", INTEGER(getListElement(re_i, "q_contribution"))[0]);
      if(sigma_r == R_NilValue){
        SEXP cov = getListElement(re_i, "covariance");
        SEXP hp = getListElement(cov, "prior_hyperparameters");
        const char *status = CHAR(STRING_ELT(getListElement(getListElement(re_i, "sampler"), "status"), 0));
        Rprintf("    covariance: unstructured, inverse-Wishart(df=%.6g), %s\n",
                REAL(getListElement(hp, "df"))[0],
                std::strcmp(status, "gibbs") == 0 ? "Gibbs" : status);
        continue;
      }
      Rprintf("    sigma^2 current: %.6g\n", REAL(getListElement(sigma_r, "current"))[0]);
      {
        SEXP hp = getListElement(sigma_r, "prior_hyperparameters");
        Rprintf("    sigma^2 prior: inverse-gamma(shape=%.6g, scale=%.6g)\n",
                REAL(hp)[0], REAL(hp)[1]);
      }
    }
  }

  Rprintf("\n");
  Rprintf("Process terms\n");
  print_rule('-', 72);
  if(process_terms_r == R_NilValue || LENGTH(process_terms_r) == 0){
    Rprintf("  none\n");
    print_rule('=', 72);
    return;
  }

  for(i = 0; i < LENGTH(process_terms_r); i++){
    SEXP term_i = VECTOR_ELT(process_terms_r, i);
    SEXP sigma_r = getListElement(term_i, "sigma_sq");
    SEXP theta_r = getListElement(term_i, "theta");
    SEXP diag_r = getListElement(term_i, "diagnostics");
    const char *type = CHAR(STRING_ELT(getListElement(term_i, "type"), 0));

    if(i > 0)
      print_rule('-', 72);

    Rprintf("  [%d] %s\n", i + 1, CHAR(STRING_ELT(getListElement(term_i, "name"), 0)));
    Rprintf("    label: %s\n", CHAR(STRING_ELT(getListElement(term_i, "label"), 0)));
    Rprintf("    type: %s\n", type);
    Rprintf("    graph index: %d\n", INTEGER(getListElement(term_i, "graph_index"))[0]);

    if(getListElement(term_i, "cov_model") != R_NilValue &&
       TYPEOF(getListElement(term_i, "cov_model")) == STRSXP)
      Rprintf("    covariance model: %s\n",
              CHAR(STRING_ELT(getListElement(term_i, "cov_model"), 0)));

    Rprintf("    SVC term: %s\n",
            LOGICAL(getListElement(term_i, "is_svc"))[0] ? "yes" : "no");
    Rprintf("    n linked observations: %d\n",
            INTEGER(getListElement(term_i, "n_obs"))[0]);
    Rprintf("    latent nodes: %d\n",
            INTEGER(getListElement(term_i, "n_node"))[0]);
    Rprintf("    q_lat: %d\n",
            INTEGER(getListElement(term_i, "q_lat"))[0]);
    Rprintf("    w_offset: %d\n",
            INTEGER(getListElement(term_i, "w_offset"))[0]);

    if(TYPEOF(getListElement(term_i, "repeated_index_collapsed")) == LGLSXP)
      Rprintf("    repeated indices collapsed: %s\n",
              LOGICAL(getListElement(term_i, "repeated_index_collapsed"))[0] ? "yes" : "no");

    if(std::strcmp(type, "nngp") == 0 || std::strcmp(type, "gp") == 0)
      Rprintf("    unique locations: %d\n",
              INTEGER(getListElement(term_i, "unique_count"))[0]);
    else if(std::strcmp(type, "ar1") == 0)
      Rprintf("    unique time points: %d\n",
              INTEGER(getListElement(term_i, "unique_count"))[0]);

    Rprintf("    sigma^2 current: %.6g\n", REAL(getListElement(sigma_r, "current"))[0]);
    Rprintf("    sigma^2 tuning: %.6g\n", REAL(getListElement(sigma_r, "tuning"))[0]);
    {
      SEXP hp = getListElement(sigma_r, "prior_hyperparameters");
      Rprintf("    sigma^2 prior: inverse-gamma(shape=%.6g, scale=%.6g)\n",
              REAL(hp)[0], REAL(hp)[1]);
    }

    if(theta_r != R_NilValue && LENGTH(theta_r) > 0){
      Rprintf("    theta parameters:\n");
      for(j = 0; j < LENGTH(theta_r); j++){
        SEXP theta_j = VECTOR_ELT(theta_r, j);
        Rprintf("      - %s: current=%.6g, tuning=%.6g, bounds=(%.6g, %.6g), prior=%s, transform=%s\n",
                CHAR(STRING_ELT(getListElement(theta_j, "name"), 0)),
                REAL(getListElement(theta_j, "current"))[0],
                REAL(getListElement(theta_j, "tuning"))[0],
                REAL(getListElement(theta_j, "lower"))[0],
                REAL(getListElement(theta_j, "upper"))[0],
                CHAR(STRING_ELT(getListElement(theta_j, "prior_type"), 0)),
                CHAR(STRING_ELT(getListElement(theta_j, "transform"), 0)));
      }
    }

    if(std::strcmp(type, "nngp") == 0){
      SEXP ns = getListElement(diag_r, "neighbor_summary");

      Rprintf("    coord dimension: %d\n",
              INTEGER(getListElement(diag_r, "coord_dim"))[0]);
      Rprintf("    neighbor count m: %d\n",
              INTEGER(getListElement(diag_r, "neighbor_count_m"))[0]);

      if(getListElement(diag_r, "ordering") != R_NilValue &&
         TYPEOF(getListElement(diag_r, "ordering")) == STRSXP)
        Rprintf("    ordering: %s\n",
                CHAR(STRING_ELT(getListElement(diag_r, "ordering"), 0)));

      if(TYPEOF(getListElement(diag_r, "repeated_collapsed")) == LGLSXP)
        Rprintf("    repeated coordinates collapsed: %s\n",
                LOGICAL(getListElement(diag_r, "repeated_collapsed"))[0] ? "yes" : "no");

      Rprintf("    neighbor summary: min=%.3f, mean=%.3f, max=%.3f\n",
              REAL(ns)[0], REAL(ns)[1], REAL(ns)[2]);

    } else if(std::strcmp(type, "dagar") == 0 || std::strcmp(type, "dagar_time") == 0){
      SEXP ns = getListElement(diag_r, "neighbor_summary");

      if(getListElement(diag_r, "ordering") != R_NilValue &&
         TYPEOF(getListElement(diag_r, "ordering")) == STRSXP)
        Rprintf("    ordering: %s\n",
                CHAR(STRING_ELT(getListElement(diag_r, "ordering"), 0)));

      if(std::strcmp(type, "dagar_time") == 0 &&
         getListElement(diag_r, "time_column") != R_NilValue &&
         TYPEOF(getListElement(diag_r, "time_column")) == STRSXP)
        Rprintf("    time/index column: %s\n",
                CHAR(STRING_ELT(getListElement(diag_r, "time_column"), 0)));

      Rprintf("    parent summary: min=%.3f, mean=%.3f, max=%.3f\n",
              REAL(ns)[0], REAL(ns)[1], REAL(ns)[2]);
      Rprintf("    zero-parent nodes: %d\n",
              INTEGER(getListElement(diag_r, "zero_parent_nodes"))[0]);

    } else if(std::strcmp(type, "gp") == 0){
      SEXP cc = getListElement(diag_r, "coord_columns");

      Rprintf("    coordinate columns: ");
      if(cc == R_NilValue || LENGTH(cc) == 0){
        Rprintf("n/a");
      } else {
        for(j = 0; j < LENGTH(cc); j++)
          Rprintf("%s%s",
                  CHAR(STRING_ELT(cc, j)),
                  (j == LENGTH(cc) - 1) ? "" : ", ");
      }
      Rprintf("\n");

      Rprintf("    coord dimension: %d\n",
              INTEGER(getListElement(diag_r, "coord_dim"))[0]);

      if(TYPEOF(getListElement(diag_r, "repeated_collapsed")) == LGLSXP)
        Rprintf("    repeated coordinates collapsed: %s\n",
                LOGICAL(getListElement(diag_r, "repeated_collapsed"))[0] ? "yes" : "no");

    } else if(std::strcmp(type, "ar1") == 0){
      if(getListElement(diag_r, "time_column") != R_NilValue &&
         TYPEOF(getListElement(diag_r, "time_column")) == STRSXP)
        Rprintf("    time/index column: %s\n",
                CHAR(STRING_ELT(getListElement(diag_r, "time_column"), 0)));

      if(TYPEOF(getListElement(diag_r, "repeated_collapsed")) == LGLSXP)
        Rprintf("    duplicated times collapsed: %s\n",
                LOGICAL(getListElement(diag_r, "repeated_collapsed"))[0] ? "yes" : "no");
    }

    Rprintf("    prior precision nnz: %d\n",
            INTEGER(getListElement(diag_r, "prior_precision_nnz"))[0]);
  }

  print_rule('=', 72);
}

SEXP build_term_description_object(SamplerState *s, SEXP backend_r, int nSamplesDone)
{
  SEXP meta_r, td_r, global_r, random_r, re_terms_r, process_terms_r;
  int t, i;
  double fill;

  meta_r = getListElement(backend_r, "term_description_meta");
  if(meta_r == R_NilValue || TYPEOF(meta_r) != VECSXP)
    Rf_error("backend$term_description_meta must be a list");

  PROTECT(td_r = Rf_duplicate(meta_r));

  global_r = getListElement(td_r, "global");
  random_r = getListElement(td_r, "random_effects");
  process_terms_r = getListElement(td_r, "process_terms");
  re_terms_r = getListElement(random_r, "terms");

  fill = NA_REAL;
  if(s->M_lat != NULL && s->M_lat->nzmax > 0)
    fill = s->cm.lnz / (double)s->M_lat->nzmax;

  setListElementByName(global_r, "qLatTotal", Rf_ScalarInteger(s->qLatTotal));
  {
    SEXP dim_r = PROTECT(Rf_allocVector(INTSXP, 2));
    INTEGER(dim_r)[0] = s->qLatSystem;
    INTEGER(dim_r)[1] = s->qLatSystem;
    setListElementByName(global_r, "M_dim", dim_r);
    UNPROTECT(1);
  }
  setListElementByName(global_r, "M_nnz", Rf_ScalarInteger((int)s->M_lat->nzmax));
  setListElementByName(global_r, "cholmod_ordering", Rf_ScalarString(Rf_mkChar(ordering_label_from_common(&s->cm))));
  setListElementByName(global_r, "cholmod_fill_ratio", Rf_ScalarReal(fill));
  setListElementByName(global_r, "cholmod_lnz", Rf_ScalarReal(s->cm.lnz));
  setListElementByName(global_r, "cholmod_flops", Rf_ScalarReal(s->cm.fl));
  setListElementByName(global_r, "factorization_status", Rf_ScalarString(Rf_mkChar("success")));

  if(re_terms_r != R_NilValue && TYPEOF(re_terms_r) == VECSXP){
    for(i = 0; i < LENGTH(re_terms_r) && i < s->nRE; i++){
      SEXP re_i = VECTOR_ELT(re_terms_r, i);
      SEXP sigma_r = getListElement(re_i, "sigma_sq");
      SEXP sampler_r = getListElement(re_i, "sampler");
      if(s->reCovDim[i] > 0){
        int k, j;
        SEXP value;

        k = s->reCovDim[i];
        PROTECT(value = Rf_allocMatrix(REALSXP, k, k));
        for(j = 0; j < k * k; j++)
          REAL(value)[j] = s->reCovariance[i][j];
        setListElementByName(getListElement(re_i, "covariance"), "current", value);
        UNPROTECT(1);
      } else {
        setListElementByName(sigma_r, "current", Rf_ScalarReal(s->sigmaSqRE[i]));
      }
      if(s->collapseIid){
        int j, active;

        active = 0;
        for(j = 0; j < s->covProp.dim; j++)
          if(s->covProp.paramType[j] >= 4 && s->covProp.paramTerm[j] == i)
            active = 1;
        setListElementByName(sampler_r, "status",
                            Rf_mkString(active ? "collapsed Metropolis" : "fixed"));
        if(sigma_r != R_NilValue)
          setListElementByName(sigma_r, "tuning", Rf_ScalarReal(s->sigmaSqRETune[i]));
      }
      if(nSamplesDone > 0)
        setListElementByName(sampler_r, "block_acceptance", Rf_ScalarReal(NA_REAL));
    }
  }

  if(process_terms_r == R_NilValue || TYPEOF(process_terms_r) != VECSXP)
    Rf_error("term_description_meta$process_terms must be a list");

  for(t = 0; t < LENGTH(process_terms_r) && t < s->nTerms; t++){
    SEXP term_r = VECTOR_ELT(process_terms_r, t);
    SEXP sigma_r = getListElement(term_r, "sigma_sq");
    SEXP diag_r = getListElement(term_r, "diagnostics");
    SEXP sampler_r = getListElement(term_r, "sampler");
    SEXP theta_list_r = getListElement(term_r, "theta");
    GraphState *g = s->graphs + s->terms[t].graphIndex;
    double nmin = NA_REAL, nmean = NA_REAL, nmax = NA_REAL;

    setListElementByName(term_r, "q_lat", Rf_ScalarInteger(s->terms[t].qLat));
    setListElementByName(term_r, "w_offset", Rf_ScalarInteger(s->terms[t].wOffset + 1));
    setListElementByName(sigma_r, "current", Rf_ScalarReal(s->terms[t].sigmaSq));
    setListElementByName(diag_r, "prior_precision_nnz", Rf_ScalarInteger(term_prior_block_nnz(s->terms + t, g)));

    if((g->type == GRAPH_NNGP || g->type == GRAPH_DAGAR || g->type == GRAPH_DAGAR_TIME) && g->nNode > 0){
      int *count = (g->type == GRAPH_NNGP) ? g->nnCount : g->parentCount;
      int n_count = (g->type == GRAPH_DAGAR_TIME) ? g->nSpace : g->nNode;
      int minv = count[0], maxv = count[0], total = 0, zero = 0;
      for(i = 0; i < n_count; i++){
        if(count[i] < minv) minv = count[i];
        if(count[i] > maxv) maxv = count[i];
        if(count[i] == 0) zero++;
        total += count[i];
      }
      nmin = (double)minv;
      nmax = (double)maxv;
      nmean = (double)total / (double)n_count;
      if(g->type == GRAPH_DAGAR || g->type == GRAPH_DAGAR_TIME)
        setListElementByName(diag_r, "zero_parent_nodes", Rf_ScalarInteger(zero));
    }
    {
      SEXP neigh_r = PROTECT(Rf_allocVector(REALSXP, 3));
      REAL(neigh_r)[0] = nmin;
      REAL(neigh_r)[1] = nmean;
      REAL(neigh_r)[2] = nmax;
      SEXP neigh_names = PROTECT(Rf_allocVector(STRSXP, 3));
      SET_STRING_ELT(neigh_names, 0, Rf_mkChar("min"));
      SET_STRING_ELT(neigh_names, 1, Rf_mkChar("mean"));
      SET_STRING_ELT(neigh_names, 2, Rf_mkChar("max"));
      Rf_setAttrib(neigh_r, R_NamesSymbol, neigh_names);
      setListElementByName(diag_r, "neighbor_summary", neigh_r);
      UNPROTECT(2);
    }

    setListElementByName(sampler_r, "block_acceptance", Rf_ScalarReal(NA_REAL));
    setListElementByName(sampler_r, "status", Rf_mkString("global covariance block"));

    if(theta_list_r != R_NilValue && TYPEOF(theta_list_r) == VECSXP){
      SEXP theta_names_r = Rf_getAttrib(theta_list_r, R_NamesSymbol);
      if(LENGTH(theta_list_r) != s->terms[t].thetaDim)
        Rf_error("term_description theta list length does not match thetaDim for term %d", t + 1);
      for(i = 0; i < s->terms[t].thetaDim; i++){
        SEXP theta_i = VECTOR_ELT(theta_list_r, i);
        if(theta_names_r != R_NilValue && TYPEOF(theta_names_r) == STRSXP && LENGTH(theta_names_r) == LENGTH(theta_list_r)){
          SEXP theta_name_field = getListElement(theta_i, "name");
          if(theta_name_field != R_NilValue && TYPEOF(theta_name_field) == STRSXP && LENGTH(theta_name_field) == 1){
            const char *theta_name_attr = CHAR(STRING_ELT(theta_names_r, i));
            const char *theta_name_val = CHAR(STRING_ELT(theta_name_field, 0));
            if(std::strcmp(theta_name_attr, theta_name_val) != 0)
              Rf_error("term_description theta metadata mismatch for term %d, theta %d", t + 1, i + 1);
          }
        }
        setListElementByName(theta_i, "current", Rf_ScalarReal(s->terms[t].theta[i]));
      }
    }
  }

  UNPROTECT(1);
  return td_r;
}
