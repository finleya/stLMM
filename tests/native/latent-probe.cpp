/* Validation-only entry point, compiled into an isolated test installation. */
#include "stLMM_internal.h"

extern "C" SEXP stLMM_latent_probe(SEXP backend, SEXP ndraw_r)
{
  SamplerState s;
  int i, j, k, nDraw = Rf_asInteger(ndraw_r), np = 0;
  SEXP out, names, H, M, VInv, moments, draws, Xprec, Xrhs, logs, restored;
  double *id, *rhs, *draw, *before;
  const void *vmax;
  const char *labels[] = {"H", "M", "Vinv", "mean", "draws", "Xprec", "Xrhs", "logs", "restored"};

  init_sampler_state(&s, backend);
  init_joint_iid_mapping(&s);
  for(i = 0; i < s.nTerms; i++){
    GraphState *g = s.graphs + s.terms[i].graphIndex;
    if(g->type == GRAPH_NNGP) update_nngp_BF(&s, g, s.terms + i);
    else if(g->type == GRAPH_GP) update_gp_Q(&s, g, s.terms + i);
  }
  build_M_lat_pattern_sparse(&s);
  s.M_lat_sym = M_cholmod_analyze(s.M_lat, &s.cm);
  s.M_lat_fac = M_cholmod_copy_factor(s.M_lat_sym, &s.cm);
  build_M_lat_index_caches(&s);
  assemble_M_lat_numeric(&s);
  if(!M_cholmod_factorize(s.M_lat, s.M_lat_fac, &s.cm)) Rf_error("probe factor failed");

  PROTECT(out = Rf_allocVector(VECSXP, 9)); np++;
  PROTECT(names = Rf_allocVector(STRSXP, 9)); np++;
  PROTECT(H = Rf_allocMatrix(REALSXP, s.n, s.qLatSystem)); np++;
  PROTECT(M = Rf_allocMatrix(REALSXP, s.qLatSystem, s.qLatSystem)); np++;
  PROTECT(VInv = Rf_allocMatrix(REALSXP, s.n, s.n)); np++;
  PROTECT(moments = Rf_allocVector(REALSXP, s.qLatSystem)); np++;
  PROTECT(draws = Rf_allocMatrix(REALSXP, nDraw, s.qLatSystem)); np++;
  PROTECT(Xprec = Rf_allocMatrix(REALSXP, s.p, s.p)); np++;
  PROTECT(Xrhs = Rf_allocVector(REALSXP, s.p)); np++;
  PROTECT(logs = Rf_allocVector(REALSXP, 7)); np++;
  PROTECT(restored = Rf_allocVector(LGLSXP, 1)); np++;

  id = (double*)R_alloc(s.qLatSystem * s.qLatSystem, sizeof(double));
  for(j = 0; j < s.qLatSystem; j++)
    for(i = 0; i < s.qLatSystem; i++) id[i + s.qLatSystem*j] = i == j;
  apply_H_multiple(&s, id, s.qLatSystem, REAL(H));
  for(i = 0; i < s.qLatSystem * s.qLatSystem; i++) REAL(M)[i] = 0;
  for(j = 0; j < s.qLatSystem; j++){
    for(k = ((int*)s.M_lat->p)[j]; k < ((int*)s.M_lat->p)[j+1]; k++){
      i = ((int*)s.M_lat->i)[k];
      REAL(M)[i + s.qLatSystem*j] = REAL(M)[j + s.qLatSystem*i] = ((double*)s.M_lat->x)[k];
    }
  }
  id = (double*)R_alloc(s.n * s.n, sizeof(double));
  for(j = 0; j < s.n; j++)
    for(i = 0; i < s.n; i++) id[i + s.n*j] = i == j;
  apply_Vinv_multiple(&s, id, s.n, REAL(VInv));
  rhs = (double*)R_alloc(s.qLatSystem, sizeof(double));
  draw = (double*)R_alloc(s.qLatSystem, sizeof(double));
  apply_H_transpose_weighted(&s, s.y, rhs);
  solve_M_lat(&s, rhs, REAL(moments));
  form_XtVinvX_and_rhs(&s, s.y, REAL(Xprec), REAL(Xrhs));
  for(j = 0; j < s.p; j++)
    for(i = j+1; i < s.p; i++) REAL(Xprec)[j + s.p*i] = REAL(Xprec)[i + s.p*j];
  REAL(logs)[0] = compute_logDetV(&s);
  REAL(logs)[1] = quadform_Vinv(&s, s.y);
  REAL(logs)[2] = logdet_Qiid_blocks(&s);
  REAL(logs)[3] = s.M_lat->nzmax;
  REAL(logs)[4] = s.cm.lnz;
  REAL(logs)[5] = log_cov_params_target_collapsed(&s, s.y);
  if(!log_cov_params_target_collapsed_try(&s, s.y, REAL(logs)+6))
    Rf_error("probe target failed");

  // Mutate and restore every IID covariance using the same pattern and caches.
  before = (double*)R_alloc(s.M_lat->nzmax, sizeof(double));
  std::memcpy(before, s.M_lat->x, s.M_lat->nzmax*sizeof(double));
  for(i = 0; i < s.nRE; i++){
    if(s.reCovDim[i] == 0) s.sigmaSqRE[i] *= 2;
    else {
      k = s.reCovDim[i];
      for(j = 0; j < k*k; j++) s.reCovariance[i][j] *= 2;
      refresh_iid_covariance_precision(&s, i);
    }
  }
  assemble_M_lat_numeric(&s);
  if(!M_cholmod_factorize(s.M_lat, s.M_lat_fac, &s.cm)) Rf_error("probe factor failed");
  for(i = 0; i < s.nRE; i++){
    if(s.reCovDim[i] == 0) s.sigmaSqRE[i] /= 2;
    else {
      k = s.reCovDim[i];
      for(j = 0; j < k*k; j++) s.reCovariance[i][j] /= 2;
      refresh_iid_covariance_precision(&s, i);
    }
  }
  assemble_M_lat_numeric(&s);
  LOGICAL(restored)[0] = std::memcmp(before, s.M_lat->x, s.M_lat->nzmax*sizeof(double)) == 0;
  if(!M_cholmod_factorize(s.M_lat, s.M_lat_fac, &s.cm)) Rf_error("probe factor failed");

  GetRNGstate();
  vmax = vmaxget();
  for(j = 0; j < nDraw; j++){
    recover_latent_draw_collapsed(&s, s.y, draw);
    for(i = 0; i < s.qLatSystem; i++) REAL(draws)[j + nDraw*i] = draw[i];
    vmaxset(vmax);
  }
  PutRNGstate();
  SEXP values[] = {H, M, VInv, moments, draws, Xprec, Xrhs, logs, restored};
  for(i = 0; i < 9; i++){
    SET_VECTOR_ELT(out, i, values[i]);
    SET_STRING_ELT(names, i, Rf_mkChar(labels[i]));
  }
  Rf_namesgets(out, names);
  free_sampler_state(&s);
  UNPROTECT(np);
  return out;
}

extern "C" SEXP stLMM_known_covariance_chain(SEXP backend, SEXP ndraw_r, SEXP joint_r)
{
  SamplerState s;
  SEXP out;
  int iter, j, nDraw, ncol, joint;
  double *beta, *alpha, *w, *latent, *resid;
  const void *vmax;

  nDraw = Rf_asInteger(ndraw_r);
  joint = Rf_asLogical(joint_r);
  init_sampler_state(&s, backend);
  if(joint) init_joint_iid_mapping(&s);
  for(j = 0; j < s.nTerms; j++){
    GraphState *g = s.graphs + s.terms[j].graphIndex;
    if(g->type == GRAPH_NNGP) update_nngp_BF(&s, g, s.terms + j);
    else if(g->type == GRAPH_GP) update_gp_Q(&s, g, s.terms + j);
  }
  build_M_lat_pattern_sparse(&s);
  s.M_lat_sym = M_cholmod_analyze(s.M_lat, &s.cm);
  s.M_lat_fac = M_cholmod_copy_factor(s.M_lat_sym, &s.cm);
  build_M_lat_index_caches(&s);
  assemble_M_lat_numeric(&s);
  if(!M_cholmod_factorize(s.M_lat, s.M_lat_fac, &s.cm)) Rf_error("chain factor failed");
  beta = (double*)R_alloc(s.p > 0 ? s.p : 1, sizeof(double));
  alpha = (double*)R_alloc(s.q, sizeof(double));
  w = (double*)R_alloc(s.qLatTotal > 0 ? s.qLatTotal : 1, sizeof(double));
  latent = (double*)R_alloc(s.qLatSystem > 0 ? s.qLatSystem : 1, sizeof(double));
  resid = (double*)R_alloc(s.n, sizeof(double));
  for(j = 0; j < s.p; j++) beta[j] = 0;
  for(j = 0; j < s.q; j++) alpha[j] = 0;
  ncol = s.p + s.q + s.qLatTotal;
  PROTECT(out = Rf_allocMatrix(REALSXP, nDraw, ncol));
  GetRNGstate();
  vmax = vmaxget();
  for(iter = 0; iter < nDraw; iter++){
    update_beta_draw_collapsed(&s, alpha, beta);
    if(joint){
      update_mean_without_alpha(&s, beta, resid);
      recover_latent_draw_collapsed(&s, resid, latent);
      for(j = 0; j < s.q; j++) alpha[j] = latent[s.qLatTotal + j];
      for(j = 0; j < s.qLatTotal; j++) w[j] = latent[j];
    } else {
      update_alpha_draw_collapsed(&s, beta, alpha);
      if(s.qLatTotal > 0){
        update_mean_without_alpha(&s, beta, resid);
        sparse_Z_mult(&s, alpha, s.Za);
        for(j = 0; j < s.n; j++) resid[j] -= s.Za[j];
        recover_w_draw_collapsed(&s, resid, w);
      }
    }
    for(j = 0; j < s.p; j++) REAL(out)[iter + nDraw*j] = beta[j];
    for(j = 0; j < s.q; j++) REAL(out)[iter + nDraw*(s.p+j)] = alpha[j];
    for(j = 0; j < s.qLatTotal; j++) REAL(out)[iter + nDraw*(s.p+s.q+j)] = w[j];
    vmaxset(vmax);
  }
  PutRNGstate();
  free_sampler_state(&s);
  UNPROTECT(1);
  return out;
}
