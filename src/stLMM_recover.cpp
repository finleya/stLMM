#include "stLMM_internal.h"
#include <BayesLogit.h>

int should_recover_iteration(SamplerState *s, int iter1)
{
  int delta;

  if(!s->recoverProcess)
    return 0;

  if(iter1 < s->recoverStart)
    return 0;

  delta = iter1 - s->recoverStart;
  return (delta % s->recoverThin) == 0;
}

void recover_latent_draw_collapsed(SamplerState *s,
                                    const double *resid,
                                    double *draw)
{
  cholmod_dense *rhs_r, *sol_r, *z_r, *noise_perm_r, *noise_r;
  double *rhsx, *solx, *zx, *noisex;
  int i;

  if(s->qLatSystem <= 0)
    return;

  rhs_r = M_cholmod_allocate_dense(s->qLatSystem, 1, s->qLatSystem, CHOLMOD_REAL, &s->cm);
  if(rhs_r == NULL || s->cm.status != CHOLMOD_OK)
    Rf_error("recover_latent_draw_collapsed: failed to allocate rhs");

  rhsx = (double*)rhs_r->x;
  apply_H_transpose_weighted(s, resid, rhsx);

  sol_r = M_cholmod_solve(CHOLMOD_A, s->M_lat_fac, rhs_r, &s->cm);
  if(sol_r == NULL || s->cm.status != CHOLMOD_OK){
    M_cholmod_free_dense(&rhs_r, &s->cm);
    Rf_error("recover_latent_draw_collapsed: cholmod_solve(A) failed for posterior mean");
  }

  solx = (double*)sol_r->x;
  for(i = 0; i < s->qLatSystem; i++)
    draw[i] = solx[i];

  z_r = M_cholmod_allocate_dense(s->qLatSystem, 1, s->qLatSystem, CHOLMOD_REAL, &s->cm);
  if(z_r == NULL || s->cm.status != CHOLMOD_OK){
    M_cholmod_free_dense(&sol_r, &s->cm);
    M_cholmod_free_dense(&rhs_r, &s->cm);
    Rf_error("recover_latent_draw_collapsed: failed to allocate z");
  }

  zx = (double*)z_r->x;
  for(i = 0; i < s->qLatSystem; i++)
    zx[i] = rnorm(0.0, 1.0);

  noise_perm_r = M_cholmod_solve(CHOLMOD_Lt, s->M_lat_fac, z_r, &s->cm);
  if(noise_perm_r == NULL || s->cm.status != CHOLMOD_OK){
    M_cholmod_free_dense(&z_r, &s->cm);
    M_cholmod_free_dense(&sol_r, &s->cm);
    M_cholmod_free_dense(&rhs_r, &s->cm);
    Rf_error("recover_latent_draw_collapsed: cholmod_solve(Lt) failed for posterior noise draw");
  }

  noise_r = M_cholmod_solve(CHOLMOD_Pt, s->M_lat_fac, noise_perm_r, &s->cm);
  if(noise_r == NULL || s->cm.status != CHOLMOD_OK){
    M_cholmod_free_dense(&noise_perm_r, &s->cm);
    M_cholmod_free_dense(&z_r, &s->cm);
    M_cholmod_free_dense(&sol_r, &s->cm);
    M_cholmod_free_dense(&rhs_r, &s->cm);
    Rf_error("recover_latent_draw_collapsed: cholmod_solve(Pt) failed for posterior noise draw");
  }

  noisex = (double*)noise_r->x;
  for(i = 0; i < s->qLatSystem; i++)
    draw[i] += noisex[i];

  M_cholmod_free_dense(&noise_r, &s->cm);
  M_cholmod_free_dense(&noise_perm_r, &s->cm);
  M_cholmod_free_dense(&z_r, &s->cm);
  M_cholmod_free_dense(&sol_r, &s->cm);
  M_cholmod_free_dense(&rhs_r, &s->cm);
}

void recover_w_draw_collapsed(SamplerState *s, const double *resid, double *wDraw)
{
  if(s->collapseIid)
    Rf_error("structured recovery requires the conditional separate-IID system");
  recover_latent_draw_collapsed(s, resid, wDraw);
}

static void recover_fixed_random_resid(SamplerState *s,
                                       const double *beta,
                                       const double *alpha,
                                       double *resid)
{
  int i, inc;
  double minusOne, one;

  inc = 1;
  minusOne = -1.0;
  one = 1.0;

  for(i = 0; i < s->n; i++)
    resid[i] = s->y[i];

  if(s->p > 0){
    F77_CALL(dgemv)("N", &s->n, &s->p,
                    &minusOne, s->X, &s->n, beta, &inc,
                    &one, resid, &inc FCONE);
  }

  if(s->q > 0)
    sparse_Z_mult(s, alpha, s->Za);
  else
    for(i = 0; i < s->n; i++) s->Za[i] = 0.0;

  for(i = 0; i < s->n; i++)
    resid[i] -= s->Za[i];
}

static void recover_pg_w_draw(SamplerState *s,
                              BayesLogit_rpg_hybrid_t pg,
                              const double *beta,
                              const double *alpha,
                              int nPgIter,
                              double *resid,
                              double *wDraw)
{
  int i, k;

  /*
    Polya-Gamma fits do not retain the in-chain PG weights or auxiliary w draws.
    For each retained posterior parameter draw we regenerate a short
    conditional PG/w chain:

      omega | beta, alpha, w, y
      w     | beta, alpha, omega, y

    The final wDraw is therefore on the same latent-process scale as Gaussian
    recovery, while obsPrecision temporarily stores omega_i for the current
    working model. The top-level fitted sampler remains collapsed; this
    auxiliary w is used only for PG recovery.
  */
  for(i = 0; i < s->qLatTotal; i++)
    wDraw[i] = 0.0;

  for(k = 0; k < nPgIter; k++){
    update_pg_working_model(s, pg, beta, alpha, wDraw);
    recover_fixed_random_resid(s, beta, alpha, resid);
    recover_w_draw_collapsed(s, resid, wDraw);
  }
}


typedef struct {
  SamplerState sampler;
  SEXP betaSamples, alphaSamples, iidSigmaSamples, iidCovSamples, tauSamples;
  SEXP residualVarianceSamples, processSigmaSamples, thetaSamples, drawIndex;
  int pgIter, rngActive;
} RecoveryState;

/* R errors and interrupts must release the same state as an ordinary return. */
static void cleanup_recovery(void *data, Rboolean jump)
{
  RecoveryState *d = (RecoveryState*)data;
  if(d->rngActive)
    PutRNGstate();
  free_sampler_state(&d->sampler);
}

static SEXP recover_w_samples(void *data)
{
  RecoveryState *d = (RecoveryState*)data;
  SamplerState *s = &d->sampler;
  SEXP out, outNames, w_samples_r, recover_iter_r, alpha_draws_r;
  double *wSamples, *wCurrent, *betaCurrent, *alphaCurrent, *resid;
  int *recoverIter;
  int nRecover, r, i, j, draw, thetaOffset;
  int jointIid, nScalar, scalarCol, covCol, b, k, row, col;
  double value;
  const void *vmax_iter;
  BayesLogit_rpg_hybrid_t pg = BayesLogit_rpg_hybrid();

  nRecover = LENGTH(d->drawIndex);
  jointIid = s->collapseIid;
  if(s->qLatSystem <= 0 || (jointIid && is_augmented_likelihood(s))){
    Rf_error("recover requires structured effects or collapsed Gaussian IID effects");
  }

  if(!Rf_isMatrix(d->betaSamples) || Rf_nrows(d->betaSamples) != s->nSamples || Rf_ncols(d->betaSamples) != s->p){
    Rf_error("beta_samples has incompatible dimensions");
  }
  if(!jointIid && (!Rf_isMatrix(d->alphaSamples) || Rf_nrows(d->alphaSamples) != s->nSamples || Rf_ncols(d->alphaSamples) != s->q)){
    Rf_error("alpha_samples has incompatible dimensions");
  }
  nScalar = 0;
  for(b = 0; b < s->nRE; b++)
    if(s->reCovDim[b] == 0) nScalar++;
  if(jointIid && (!Rf_isReal(d->iidSigmaSamples) || !Rf_isMatrix(d->iidSigmaSamples) ||
                 Rf_nrows(d->iidSigmaSamples) != s->nSamples || Rf_ncols(d->iidSigmaSamples) != nScalar ||
                 !Rf_isReal(d->iidCovSamples) || !Rf_isMatrix(d->iidCovSamples) ||
                 Rf_nrows(d->iidCovSamples) != s->nSamples || Rf_ncols(d->iidCovSamples) != s->nRECovElements)){
    Rf_error("IID covariance samples have incompatible dimensions");
  }
  if(!Rf_isReal(d->tauSamples) || LENGTH(d->tauSamples) != s->nSamples){
    Rf_error("tau_sq_samples has incompatible length");
  }
  if(!Rf_isMatrix(d->residualVarianceSamples) ||
     Rf_nrows(d->residualVarianceSamples) != s->nSamples ||
     Rf_ncols(d->residualVarianceSamples) != s->nResidualGroup){
    Rf_error("residual_variance_samples has incompatible dimensions");
  }
  if(!Rf_isMatrix(d->processSigmaSamples) ||
     Rf_nrows(d->processSigmaSamples) != s->nSamples ||
     Rf_ncols(d->processSigmaSamples) != s->nTerms){
    Rf_error("sigma_sq_samples has incompatible dimensions");
  }
  if(!Rf_isMatrix(d->thetaSamples) ||
     Rf_nrows(d->thetaSamples) != s->nSamples ||
     Rf_ncols(d->thetaSamples) != s->nThetaTotal){
    Rf_error("theta_samples has incompatible dimensions");
  }

  for(i = 0; i < s->nTerms; i++){
    GraphState *g = s->graphs + s->terms[i].graphIndex;
    TermState *term = s->terms + i;

    if(g->type == GRAPH_NNGP)
      update_nngp_BF(s, g, term);
    else if(g->type == GRAPH_GP)
      update_gp_Q(s, g, term);
  }

  build_M_lat_pattern_sparse(s);
  s->M_lat_sym = M_cholmod_analyze(s->M_lat, &s->cm);
  if(s->M_lat_sym == NULL || s->cm.status != CHOLMOD_OK){
    Rf_error("cholmod analyze failed for latent M");
  }

  s->M_lat_fac = M_cholmod_copy_factor(s->M_lat_sym, &s->cm);
  if(s->M_lat_fac == NULL || s->cm.status != CHOLMOD_OK){
    Rf_error("cholmod copy_factor failed for latent M");
  }

  build_M_lat_index_caches(s);

  PROTECT(out = Rf_allocVector(VECSXP, 3));
  PROTECT(outNames = Rf_allocVector(STRSXP, 3));
  PROTECT(w_samples_r = Rf_allocMatrix(REALSXP, nRecover, s->qLatTotal));
  PROTECT(recover_iter_r = Rf_allocVector(INTSXP, nRecover));
  PROTECT(alpha_draws_r = Rf_allocMatrix(REALSXP, nRecover, jointIid ? s->q : 0));

  wSamples = REAL(w_samples_r);
  recoverIter = INTEGER(recover_iter_r);
  wCurrent = (double*)R_alloc(s->qLatSystem, sizeof(double));
  betaCurrent = (double*)R_alloc(s->p > 0 ? s->p : 1, sizeof(double));
  alphaCurrent = (double*)R_alloc(s->q > 0 ? s->q : 1, sizeof(double));
  resid = (double*)R_alloc(s->n > 0 ? s->n : 1, sizeof(double));

  GetRNGstate();
  d->rngActive = 1;
  for(r = 0; r < nRecover; r++){
    R_CheckUserInterrupt();
    vmax_iter = vmaxget();
    draw = INTEGER(d->drawIndex)[r];
    if(draw < 1 || draw > s->nSamples)
      Rf_error("draw_index contains an out-of-range posterior draw");
    draw--;

    s->tauSq = REAL(d->tauSamples)[draw];
    if(s->residualModel == 2 || s->residualModel == 3){
      for(i = 0; i < s->nResidualGroup; i++){
        s->residualVariance[i] = REAL(d->residualVarianceSamples)[draw + s->nSamples * i];
        if(!(s->residualVariance[i] > 0.0))
          Rf_error("residual_variance_samples contains a non-positive value");
      }
    } else if(!(s->tauSq > 0.0)){
      Rf_error("tau_sq_samples contains a non-positive value");
    }
    thetaOffset = 0;
    for(i = 0; i < s->nTerms; i++){
      TermState *term = s->terms + i;
      GraphState *g = s->graphs + term->graphIndex;

      term->sigmaSq = REAL(d->processSigmaSamples)[draw + s->nSamples * i];
      if(!(term->sigmaSq > 0.0))
        Rf_error("sigma_sq_samples contains a non-positive value");

      for(j = 0; j < term->thetaDim; j++)
        term->theta[j] = REAL(d->thetaSamples)[draw + s->nSamples * (thetaOffset + j)];
      thetaOffset += term->thetaDim;

      if(g->type == GRAPH_NNGP)
        update_nngp_BF(s, g, term);
      else if(g->type == GRAPH_GP)
        update_gp_Q(s, g, term);
    }

    for(i = 0; i < s->p; i++)
      betaCurrent[i] = REAL(d->betaSamples)[draw + s->nSamples * i];
    for(i = 0; !jointIid && i < s->q; i++)
      alphaCurrent[i] = REAL(d->alphaSamples)[draw + s->nSamples * i];

    if(jointIid){
      scalarCol = covCol = 0;
      for(b = 0; b < s->nRE; b++){
        k = s->reCovDim[b];
        if(k == 0){
          value = REAL(d->iidSigmaSamples)[draw + s->nSamples * scalarCol++];
          if(!R_FINITE(value) || value <= 0.0)
            Rf_error("IID variance samples must be finite and positive");
          s->sigmaSqRE[b] = value;
        } else {
          for(col = 0; col < k; col++)
            for(row = col; row < k; row++){
              value = REAL(d->iidCovSamples)[draw + s->nSamples * covCol++];
              if(!R_FINITE(value))
                Rf_error("IID covariance samples must be finite");
              s->reCovariance[b][row + k * col] = s->reCovariance[b][col + k * row] = value;
            }
          refresh_iid_covariance_precision(s, b);
        }
      }
    }

    if(is_pg_likelihood(s)){
      recover_pg_w_draw(s, pg, betaCurrent, alphaCurrent,
                        d->pgIter, resid, wCurrent);
    } else {
      refresh_observation_precision(s);
      refactor_current_M(s, "stLMM_recover_w");
      if(jointIid){
        update_mean_without_alpha(s, betaCurrent, resid);
        recover_latent_draw_collapsed(s, resid, wCurrent);
      } else {
        recover_fixed_random_resid(s, betaCurrent, alphaCurrent, resid);
        recover_w_draw_collapsed(s, resid, wCurrent);
      }
    }

    for(i = 0; i < s->qLatTotal; i++)
      wSamples[r + nRecover * i] = wCurrent[i];
    for(i = 0; jointIid && i < s->q; i++)
      REAL(alpha_draws_r)[r + nRecover * i] = wCurrent[s->qLatTotal + i];

    recoverIter[r] = draw + 1;
    vmaxset(vmax_iter);
  }

  SET_VECTOR_ELT(out, 0, w_samples_r);
  SET_VECTOR_ELT(out, 1, recover_iter_r);
  SET_VECTOR_ELT(out, 2, alpha_draws_r);
  SET_STRING_ELT(outNames, 0, Rf_mkChar("w_samples"));
  SET_STRING_ELT(outNames, 1, Rf_mkChar("recover_iter"));
  SET_STRING_ELT(outNames, 2, Rf_mkChar("alpha_samples"));
  Rf_namesgets(out, outNames);

  UNPROTECT(5);
  return out;
}

extern "C" {

  SEXP stLMM_recover_w(SEXP backend_r,
                      SEXP beta_samples_r,
                      SEXP alpha_samples_r,
                      SEXP iid_sigma_samples_r,
                      SEXP iid_cov_samples_r,
                      SEXP tau_samples_r,
                      SEXP residual_variance_samples_r,
                      SEXP process_sigma_samples_r,
                      SEXP theta_samples_r,
                      SEXP draw_index_r,
                      SEXP pg_iter_r)
  {
    RecoveryState d;
    SEXP joint_r;
    int jointIid;

    if(!Rf_isInteger(draw_index_r))
      Rf_error("draw_index must be an integer vector");
    if(LENGTH(draw_index_r) < 1)
      Rf_error("draw_index must select at least one posterior draw");
    d.pgIter = as_int_scalar(pg_iter_r, "recover pg_iter");
    if(d.pgIter < 1)
      Rf_error("recover pg_iter must be a positive integer");

    d.betaSamples = beta_samples_r;
    d.alphaSamples = alpha_samples_r;
    d.iidSigmaSamples = iid_sigma_samples_r;
    d.iidCovSamples = iid_cov_samples_r;
    d.tauSamples = tau_samples_r;
    d.residualVarianceSamples = residual_variance_samples_r;
    d.processSigmaSamples = process_sigma_samples_r;
    d.thetaSamples = theta_samples_r;
    d.drawIndex = draw_index_r;
    d.rngActive = 0;
    joint_r = getListElement(backend_r, "defer_iid_recovery");
    jointIid = joint_r == R_NilValue ? 0 : as_flag_scalar(joint_r, "backend$defer_iid_recovery");
    init_sampler_state(&d.sampler, backend_r, jointIid);
    return R_UnwindProtect(recover_w_samples, &d, cleanup_recovery, &d, NULL);
  }

}
