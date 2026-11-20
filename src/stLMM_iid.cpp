#include "stLMM_internal.h"

/* Scalar variance density in log coordinates. Separate IID updates also
   condition on the current coefficients; collapsed updates integrate them out. */
double log_iid_variance_target(SamplerState *s, int b, double variance)
{
  int j;
  double out, ss;

  if(!R_FINITE(variance) || variance <= 0.0 || !R_FINITE(1.0 / variance))
    return R_NegInf;
  out = log_prior_kernel_original(variance, s->sigmaSqREPrior[b]) + std::log(variance);
  if(!s->collapseIid){
    ss = 0.0;
    for(j = 0; j < s->q; j++)
      if(s->reBlockID[j] == b)
        ss += s->alpha[j] * s->alpha[j];
    out -= 0.5 * (s->reBlockSize[b] * std::log(variance) + ss / variance);
  }
  return out;
}

/* A separate-IID variance block does not change the structured system.
   Reuse proposal/adaptation state without assembling or factoring M. */
int update_iid_variance_block(SamplerState *s, CovProposalBlock *block)
{
  int j, b;
  double current, proposal;
  CovProposalState *p;

  p = &s->covProp;
  current = 0.0;
  for(j = 0; j < block->dim; j++){
    b = p->paramTerm[block->paramIndex[j]];
    current += log_iid_variance_target(s, b, s->sigmaSqRE[b]);
  }
  if(!R_FINITE(current))
    Rf_error("non-finite conditional IID variance target at the current state");
  draw_adaptive_cov_block_proposal(s, block);
  proposal = 0.0;
  for(j = 0; j < block->dim; j++){
    b = p->paramTerm[block->paramIndex[j]];
    proposal += log_iid_variance_target(s, b, p->sigmaSqRE[b]);
  }
  if(!R_FINITE(proposal))
    return 0;
  if(std::log(runif(0.0, 1.0)) >= proposal - current)
    return 0;
  for(j = 0; j < block->dim; j++){
    b = p->paramTerm[block->paramIndex[j]];
    s->sigmaSqRE[b] = p->sigmaSqRE[b];
  }
  return 1;
}

/* Densities in proposal coordinates: log variance for scalar blocks and lower
   Cholesky entries with log diagonals for IW blocks. Each covariance contributes
   its prior and Jacobian once, irrespective of the number of groups. */
double log_iid_covariance_target(SamplerState *s)
{
  int b, k, j;
  double out, trace;

  out = 0.0;
  for(b = 0; b < s->nRE; b++){
    k = s->reCovDim[b];
    if(k == 0){
      if(s->collapseIid || s->sigmaSqREPrior[b].family != PRIOR_IG)
        out += log_iid_variance_target(s, b, s->sigmaSqRE[b]);
    } else if(s->collapseIid){
      trace = 0.0;
      for(j = 0; j < k * k; j++)
        trace += s->reCovScale[b][j] * s->reCovPrecision[b][j];
      out -= 0.5 * ((s->reCovDF[b] + k + 1.0) * s->reCovLogDet[b] + trace);
      out += k * std::log(2.0);
      for(j = 0; j < k; j++)
        out += (k - j + 1.0) * std::log(s->reCovChol[b][j + k * j]);
    }
  }
  return out;
}

void copy_iid_covariance_proposal(SamplerState *s, CovProposalState *p)
{
  int b, k;

  p->iidValid = 1;
  for(b = 0; b < s->nRE; b++){
    p->sigmaSqRE[b] = s->sigmaSqRE[b];
    k = s->reCovDim[b];
    if(k > 0){
      std::memcpy(p->reCovariance[b], s->reCovariance[b], k * k * sizeof(double));
      std::memcpy(p->reCovPrecision[b], s->reCovPrecision[b], k * k * sizeof(double));
      std::memcpy(p->reCovChol[b], s->reCovChol[b], k * k * sizeof(double));
      p->reCovLogDet[b] = s->reCovLogDet[b];
    }
  }
}

void init_iid_covariance_proposal(SamplerState *s, CovProposalState *p)
{
  int b, k;

  p->sigmaSqRE = (double*)R_alloc(s->nRE, sizeof(double));
  p->reCovariance = (double**)R_alloc(s->nRE, sizeof(double*));
  p->reCovPrecision = (double**)R_alloc(s->nRE, sizeof(double*));
  p->reCovChol = (double**)R_alloc(s->nRE, sizeof(double*));
  p->reCovLogDet = (double*)R_alloc(s->nRE, sizeof(double));
  for(b = 0; b < s->nRE; b++){
    k = s->reCovDim[b];
    p->reCovariance[b] = p->reCovPrecision[b] = p->reCovChol[b] = NULL;
    p->reCovLogDet[b] = 0.0;
    if(k > 0){
      p->reCovariance[b] = (double*)R_alloc(k * k, sizeof(double));
      p->reCovPrecision[b] = (double*)R_alloc(k * k, sizeof(double));
      p->reCovChol[b] = (double*)R_alloc(k * k, sizeof(double));
    }
  }
  copy_iid_covariance_proposal(s, p);
}

/* Form the small covariance and precision from a proposed Cholesky factor.
   Invalid floating-point proposals are rejected without covariance jitter. */
int refresh_iid_covariance_proposal(SamplerState *s, CovProposalState *p, int b)
{
  int k, row, col, j, info;
  double *L, *C, *P, *check, one, zero;
  char lower, noTrans;

  k = s->reCovDim[b];
  L = p->reCovChol[b];
  C = p->reCovariance[b];
  P = p->reCovPrecision[b];
  p->reCovLogDet[b] = 0.0;
  for(j = 0; j < k * k; j++){
    if(!R_FINITE(L[j]))
      return 0;
    P[j] = 0.0;
  }
  for(j = 0; j < k; j++){
    if(L[j + k * j] <= 0.0)
      return 0;
    p->reCovLogDet[b] += 2.0 * std::log(L[j + k * j]);
    P[j + k * j] = 1.0;
  }
  one = 1.0;
  zero = 0.0;
  lower = 'L';
  noTrans = 'N';
  F77_CALL(dsyrk)(&lower, &noTrans, &k, &k, &one, L, &k, &zero, C, &k FCONE FCONE);
  for(col = 0; col < k; col++)
    for(row = col; row < k; row++){
      if(!R_FINITE(C[row + k * col]))
        return 0;
      C[col + k * row] = C[row + k * col];
    }
  F77_CALL(dpotrs)(&lower, &k, &k, L, &k, P, &k, &info FCONE);
  if(info != 0)
    return 0;
  for(j = 0; j < k * k; j++)
    if(!R_FINITE(P[j]))
      return 0;
  /* L is nonsingular mathematically, but its product can lose positive
     definiteness in floating point for extreme proposals. */
  check = (double*)R_alloc(k * k, sizeof(double));
  std::memcpy(check, C, k * k * sizeof(double));
  F77_CALL(dpotrf)(&lower, &k, check, &k, &info FCONE);
  return info == 0 && R_FINITE(p->reCovLogDet[b]);
}
