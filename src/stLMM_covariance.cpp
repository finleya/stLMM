#include "stLMM_internal.h"

double graph_distance(GraphState *g, TermState *term, int i, int j, double *u)
{
  double h;
  int d, ds, a;
  double diff, ss;

  h = 0.0;
  *u = 0.0;

  if(term->distanceMode == COR_SINGLE){
    ss = 0.0;
    for(a = 0; a < g->dim; a++){
      diff = g->coords[i + g->nNode * a] - g->coords[j + g->nNode * a];
      ss += diff * diff;
    }
    return std::sqrt(ss);
  }

  d = g->dim;
  ds = d - 1;

  ss = 0.0;
  for(a = 0; a < ds; a++){
    diff = g->coords[i + g->nNode * a] - g->coords[j + g->nNode * a];
    ss += diff * diff;
  }
  h = std::sqrt(ss);
  *u = std::fabs(g->coords[i + g->nNode * (d - 1)] - g->coords[j + g->nNode * (d - 1)]);

  return h;
}

static int graph_spatial_dim(GraphState *g, TermState *term)
{
  if(term->distanceMode == COR_SPACE_TIME)
    return g->dim > 1 ? g->dim - 1 : 1;
  return g->dim > 0 ? g->dim : 1;
}

static double nngp_matern_cor_eval(TermState *term, double h, double *bk)
{
  double phi, nu, x, val;

  phi = term->theta[0];
  nu = term->theta[1];
  x = phi * h;

  if(x <= 0.0)
    return 1.0;

  val = pow(x, nu) / (pow(2.0, nu - 1.0) * gammafn(nu)) *
    bessel_k_ex(x, nu, 1.0, bk);

  if(!R_FINITE(val) || val < 0.0)
    return 0.0;
  return val;
}

void update_nngp_BF(SamplerState *s, GraphState *g, TermState *term)
{
  int i, j, k, l;
  int m, start, info, inc, ldc;
  char lower;
  double one, zero;
  double h, u;
  double *c, *C, *bk;
  double sigmaRef;
  int status, statusNode, nThreads, isMatern, spatialDim;

  if(g->type != GRAPH_NNGP)
    return;

  if(term->corFun == NULL)
    Rf_error("NNGP term missing correlation function");

  /*
    B and F are computed on the unit-variance correlation scale.
    The process variance sigmaSq is applied later during observation-level
    assembly into Q.
  */
  sigmaRef = 1.0;
  inc = 1;
  lower = 'L';
  one = 1.0;
  zero = 0.0;
  ldc = s->scratch_BF_m;
  status = 0;
  statusNode = -1;
  nThreads = s->nOmpThreads;
  isMatern = std::strcmp(get_cor_model_info(term->covModelIndex)->name, "matern") == 0;
  spatialDim = graph_spatial_dim(g, term);

#ifdef _OPENMP
#pragma omp parallel for num_threads(nThreads) private(j, k, l, m, start, info, h, u, c, C, bk)
#endif
  for(i = 0; i < g->nNode; i++){
    int threadID = 0;
#ifdef _OPENMP
    threadID = omp_get_thread_num();
#endif
    start = g->nnStart[i];
    m = g->nnCount[i];

    if(m <= 0){
      term->F[i] = 1.0 / sigmaRef;
      continue;
    }

    if(m > s->scratch_BF_m)
    {
#ifdef _OPENMP
#pragma omp critical
#endif
      {
        if(status == 0){
          status = -1;
          statusNode = i;
        }
      }
      continue;
    }

    c = s->scratch_BF_c + (size_t)threadID * (size_t)s->scratch_BF_m;
    C = s->scratch_BF_C + (size_t)threadID * (size_t)s->scratch_BF_m * (size_t)s->scratch_BF_m;
    bk = s->scratch_BF_bk + (size_t)threadID * (size_t)s->scratch_BF_bk_n;

    for(k = 0; k < ldc * m; k++)
      C[k] = 0.0;

    for(k = 0; k < m; k++){
      int jl;
      double hl, ul;

      j = g->nnIndx[start + k];
      h = graph_distance(g, term, i, j, &u);
      c[k] = sigmaRef * (isMatern ?
                         nngp_matern_cor_eval(term, h, bk) :
                         term->corFun(term->theta, h, u, spatialDim));

      for(l = 0; l <= k; l++){
        jl = g->nnIndx[start + l];
        hl = graph_distance(g, term, j, jl, &ul);
        C[k + ldc * l] = sigmaRef * (isMatern ?
                                     nngp_matern_cor_eval(term, hl, bk) :
                                     term->corFun(term->theta, hl, ul, spatialDim));
      }
    }

    F77_CALL(dpotrf)(&lower, &m, C, &ldc, &info FCONE);
    if(info != 0){
#ifdef _OPENMP
#pragma omp critical
#endif
      {
        if(status == 0){
          status = info;
          statusNode = i;
        }
      }
      continue;
    }

    F77_CALL(dpotri)(&lower, &m, C, &ldc, &info FCONE);
    if(info != 0){
#ifdef _OPENMP
#pragma omp critical
#endif
      {
        if(status == 0){
          status = info;
          statusNode = i;
        }
      }
      continue;
    }

    F77_CALL(dsymv)(&lower, &m, &one, C, &ldc, c, &inc, &zero, term->B + start, &inc FCONE);
    term->F[i] = 1.0 / (sigmaRef - F77_CALL(ddot)(&m, term->B + start, &inc, c, &inc));
  }

  if(status == -1)
    Rf_error("update_nngp_BF: neighbor count exceeds scratch_BF_m at node %d", statusNode + 1);
  if(status != 0)
    Rf_error("update_nngp_BF: local LAPACK factorization failed at node %d with info %d", statusNode + 1, status);
}

void update_gp_Q(SamplerState *s, GraphState *g, TermState *term)
{
  int i, j, info, nNode, ldq;
  double h, u, corVal, sigmaSqInv, logDetR;
  char uplo = 'L';

  (void)s;

  if(g == NULL || term == NULL)
    Rf_error("update_gp_Q: NULL graph or term pointer");

  if(g->type != GRAPH_GP)
    Rf_error("update_gp_Q requires a GP graph");

  if(term->nNode != g->nNode)
    Rf_error("update_gp_Q: term->nNode does not match graph->nNode");

  if(term->corFun == NULL)
    Rf_error("update_gp_Q: NULL correlation function pointer");

  if(term->sigmaSq <= 0.0)
    Rf_error("update_gp_Q requires sigmaSq > 0");

  nNode = term->nNode;
  ldq = nNode;

  if(nNode <= 0)
    Rf_error("update_gp_Q: nNode must be positive");

  /* Allocate persistent dense GP precision storage once per term.
     I'm not using R_alloc here: Q lives on TermState and is reused across
     covariance updates and sampler iterations, so it must survive until
     free_sampler_state(). */
  if(term->Q == NULL)
    term->Q = R_Calloc((size_t)term->nNode * (size_t)term->nNode, double);

  /* Build dense correlation matrix R in term->Q (full matrix, column major) */
  for(j = 0; j < nNode; j++){
    term->Q[j + ldq * j] = 1.0;

    for(i = j + 1; i < nNode; i++){
      h = graph_distance(g, term, i, j, &u);
      corVal = term->corFun(term->theta, h, u, graph_spatial_dim(g, term));

      if(!R_finite(corVal))
        Rf_error("update_gp_Q: non-finite correlation value");

      term->Q[i + ldq * j] = corVal;
      term->Q[j + ldq * i] = corVal;
    }
  }

  /* Cholesky factorization of R */
  F77_CALL(dpotrf)(&uplo, &nNode, term->Q, &ldq, &info FCONE);
  if(info != 0){
    if(info > 0)
      Rf_error("update_gp_Q: dense GP correlation matrix is not positive definite (leading minor %d)", info);
    else
      Rf_error("update_gp_Q: dpotrf failed with info=%d", info);
  }

  /* log|R| from Cholesky factor */
  logDetR = 0.0;
  for(i = 0; i < nNode; i++){
    if(term->Q[i + ldq * i] <= 0.0)
      Rf_error("update_gp_Q: non-positive Cholesky diagonal");
    logDetR += 2.0 * std::log(term->Q[i + ldq * i]);
  }

  /* Invert R in place; lower triangle becomes valid */
  F77_CALL(dpotri)(&uplo, &nNode, term->Q, &ldq, &info FCONE);
  if(info != 0)
    Rf_error("update_gp_Q: dpotri failed with info=%d", info);

  /* Scale lower triangle to Q = sigmaSq^{-1} R^{-1} */
  sigmaSqInv = 1.0 / term->sigmaSq;
  for(j = 0; j < nNode; j++){
    for(i = j; i < nNode; i++)
      term->Q[i + ldq * j] *= sigmaSqInv;
  }

  term->logDetQ = -((double)nNode) * std::log(term->sigmaSq) - logDetR;
}

/*==========================================================================*/
/* sparse symbolic pattern helpers                                          */

int find_sparse_entry_index(cholmod_sparse *A, int row, int col)
{
  int left, right, mid;
  int *Ap, *Ai;

  if(row < col)
    Rf_error("find_sparse_entry_index called with row < col");

  Ap = (int*)A->p;
  Ai = (int*)A->i;

  left = Ap[col];
  right = Ap[col + 1] - 1;

  while(left <= right){
    mid = left + (right - left) / 2;

    if(Ai[mid] == row)
      return mid;
    else if(Ai[mid] < row)
      left = mid + 1;
    else
      right = mid - 1;
  }

  Rf_error("sparse entry (%d,%d) missing from symbolic pattern", row + 1, col + 1);
  return -1;
}

void car_Q_pattern_append(PatternCol *cols, int col, int row)
{
  int newCap;
  int *tmp;

  if(row < col)
    Rf_error("car_Q_pattern_append requires row >= col");

  if(cols[col].len >= cols[col].cap){
    newCap = (cols[col].cap == 0) ? 4 : 2 * cols[col].cap;
    tmp = R_Realloc(cols[col].rows, newCap, int);
    if(tmp == NULL)
      Rf_error("car_Q_pattern_append: realloc failed");
    cols[col].rows = tmp;
    cols[col].cap = newCap;
  }

  cols[col].rows[cols[col].len] = row;
  cols[col].len++;
}

void init_car_Q_sparse_cache(SamplerState *s, TermState *term, GraphState *g)
{
  PatternCol *cols;
  int n, col, edge, k, nz, a, b;
  int *Qp, *Qi;
  double *Qx;

  if(term->carQ != NULL && term->carQFac != NULL)
    return;

  n = (g->type == GRAPH_CAR_TIME) ? g->nSpace : term->nNode;
  cols = R_Calloc(n > 0 ? n : 1, PatternCol);
  for(col = 0; col < n; col++){
    cols[col].rows = NULL;
    cols[col].len = 0;
    cols[col].cap = 0;
  }

  for(col = 0; col < n; col++)
    car_Q_pattern_append(cols, col, col);

  for(edge = 0; edge < g->nEdge; edge++){
    a = g->edgeI[edge];
    b = g->edgeJ[edge];
    if(a >= b)
      car_Q_pattern_append(cols, b, a);
    else
      car_Q_pattern_append(cols, a, b);
  }

  nz = 0;
  for(col = 0; col < n; col++){
    if(cols[col].len > 1)
      qsort(cols[col].rows, (size_t)cols[col].len, sizeof(int), cmp_int_asc);

    if(cols[col].len > 0){
      int u = 1;
      for(k = 1; k < cols[col].len; k++){
        if(cols[col].rows[k] != cols[col].rows[u - 1]){
          cols[col].rows[u] = cols[col].rows[k];
          u++;
        }
      }
      cols[col].len = u;
    }
    nz += cols[col].len;
  }

  term->carQ = M_cholmod_allocate_sparse(n, n, nz, 1, 1, -1, CHOLMOD_REAL, &s->cm);
  if(term->carQ == NULL || s->cm.status != CHOLMOD_OK)
    Rf_error("failed to allocate sparse CAR precision");

  Qp = (int*)term->carQ->p;
  Qi = (int*)term->carQ->i;
  Qx = (double*)term->carQ->x;

  nz = 0;
  Qp[0] = 0;
  for(col = 0; col < n; col++){
    for(k = 0; k < cols[col].len; k++){
      Qi[nz] = cols[col].rows[k];
      Qx[nz] = 0.0;
      nz++;
    }
    Qp[col + 1] = nz;
  }

  for(col = 0; col < n; col++){
    if(cols[col].rows != NULL)
      R_Free(cols[col].rows);
  }
  R_Free(cols);

  term->carQDiagCacheIdx = R_Calloc(n > 0 ? n : 1, int);
  term->carQOffCacheIdx = R_Calloc(g->nEdge > 0 ? g->nEdge : 1, int);

  for(col = 0; col < n; col++)
    term->carQDiagCacheIdx[col] = find_sparse_entry_index(term->carQ, col, col);

  for(edge = 0; edge < g->nEdge; edge++){
    a = g->edgeI[edge];
    b = g->edgeJ[edge];
    if(a >= b)
      term->carQOffCacheIdx[edge] = find_sparse_entry_index(term->carQ, a, b);
    else
      term->carQOffCacheIdx[edge] = find_sparse_entry_index(term->carQ, b, a);
  }

  term->carQFac = M_cholmod_analyze(term->carQ, &s->cm);
  if(term->carQFac == NULL || s->cm.status != CHOLMOD_OK)
    Rf_error("cholmod analyze failed for sparse CAR precision");
}

int logdet_car_Q_sparse_try(SamplerState *s, TermState *term, GraphState *g, double *out)
{
  int i, edge, nz;
  double rho, sigmaSqInv, diagVal, *Qx;

  init_car_Q_sparse_cache(s, term, g);

  rho = term->theta[0];
  if(rho <= 0.0 || rho >= 1.0)
    return 0;

  sigmaSqInv = 1.0 / term->sigmaSq;
  Qx = (double*)term->carQ->x;
  nz = (int)term->carQ->nzmax;
  for(i = 0; i < nz; i++)
    Qx[i] = 0.0;

  for(i = 0; i < ((g->type == GRAPH_CAR_TIME) ? g->nSpace : term->nNode); i++){
    if(g->carModel == CAR_MODEL_LEROUX)
      diagVal = (1.0 - rho) + rho * g->degree[i];
    else
      diagVal = g->degree[i];
    Qx[term->carQDiagCacheIdx[i]] += sigmaSqInv * diagVal;
  }

  for(edge = 0; edge < g->nEdge; edge++)
    Qx[term->carQOffCacheIdx[edge]] += -sigmaSqInv * rho * g->edgeW[edge];

  if(!M_cholmod_factorize(term->carQ, term->carQFac, &s->cm) || s->cm.status != CHOLMOD_OK){
    s->cm.status = CHOLMOD_OK;
    return 0;
  }

  *out = logDetFactor(term->carQFac);
  return 1;
}

double logdet_car_Q_sparse(SamplerState *s, TermState *term, GraphState *g)
{
  double out;

  if(!logdet_car_Q_sparse_try(s, term, g, &out))
    Rf_error("CAR precision is not positive definite for current rho");
  return out;
}

int logdet_car_time_Q_try(SamplerState *s, TermState *term, GraphState *g, double *out)
{
  double thetaTime, phi, den, logdetSpaceScaled, logdetSpaceUnscaled, logdetTime;
  int nSpace, nTime, i;

  thetaTime = term->theta[1];
  if(g->timeModel == TIME_MODEL_AR1){
    if(std::fabs(thetaTime) >= 1.0)
      return 0;
  } else if(g->timeModel == TIME_MODEL_EXP){
    if(thetaTime <= 0.0)
      return 0;
  }

  nSpace = g->nSpace;
  nTime = g->nTime;

  if(!logdet_car_Q_sparse_try(s, term, g, &logdetSpaceScaled))
    return 0;
  logdetSpaceUnscaled = logdetSpaceScaled + nSpace * std::log(term->sigmaSq);

  if(nTime <= 1){
    logdetTime = 0.0;
  } else if(g->timeModel == TIME_MODEL_AR1){
    phi = thetaTime;
    den = 1.0 - phi * phi;
    if(den <= 0.0 || !std::isfinite(den))
      return 0;
    logdetTime = -(nTime - 1) * std::log(den);
  } else {
    logdetTime = 0.0;
    for(i = 0; i < nTime - 1; i++){
      phi = std::exp(-thetaTime * g->timeDelta[i]);
      den = 1.0 - phi * phi;
      if(den <= 0.0 || !std::isfinite(den))
        return 0;
      logdetTime -= std::log(den);
    }
  }

  *out = -(nSpace * nTime) * std::log(term->sigmaSq) +
    nTime * logdetSpaceUnscaled +
    nSpace * logdetTime;
  return 1;
}

double logdet_car_time_Q(SamplerState *s, TermState *term, GraphState *g)
{
  double out;

  if(!logdet_car_time_Q_try(s, term, g, &out))
    Rf_error("CAR-time precision is not positive definite for current parameters");
  return out;
}

int logdet_dagar_Q_try(TermState *term, GraphState *g, double *out)
{
  int i, m;
  double rho, rho2, denom;

  rho = term->theta[0];
  if(rho <= 0.0 || rho >= 1.0)
    return 0;

  rho2 = rho * rho;
  *out = -((double)term->nNode) * std::log(term->sigmaSq);

  for(i = 0; i < g->nNode; i++){
    m = g->parentCount[i];
    denom = 1.0 + ((double)m - 1.0) * rho2;
    if(denom <= 0.0 || !std::isfinite(denom))
      return 0;
    *out += std::log(denom) - std::log(1.0 - rho2);
  }

  return 1;
}

double logdet_dagar_Q(TermState *term, GraphState *g)
{
  double out;

  if(!logdet_dagar_Q_try(term, g, &out))
    Rf_error("DAGAR precision is not positive definite for current rho");
  return out;
}

int logdet_time_precision_try(GraphState *g, double thetaTime, double *out)
{
  double phi, den;
  int i;

  if(g->timeModel == TIME_MODEL_AR1){
    if(std::fabs(thetaTime) >= 1.0)
      return 0;
    if(g->nTime <= 1){
      *out = 0.0;
      return 1;
    }
    den = 1.0 - thetaTime * thetaTime;
    if(den <= 0.0 || !std::isfinite(den))
      return 0;
    *out = -(g->nTime - 1) * std::log(den);
    return 1;
  }

  if(g->timeModel == TIME_MODEL_EXP){
    if(thetaTime <= 0.0)
      return 0;
    *out = 0.0;
    for(i = 0; i < g->nTime - 1; i++){
      phi = std::exp(-thetaTime * g->timeDelta[i]);
      den = 1.0 - phi * phi;
      if(den <= 0.0 || !std::isfinite(den))
        return 0;
      *out -= std::log(den);
    }
    return 1;
  }

  return 0;
}

int logdet_dagar_time_Q_try(TermState *term, GraphState *g, double *out)
{
  int i, m;
  double rho, rho2, denom, logdetSpaceUnscaled, logdetTime;

  rho = term->theta[0];
  if(rho <= 0.0 || rho >= 1.0)
    return 0;

  if(!logdet_time_precision_try(g, term->theta[1], &logdetTime))
    return 0;

  rho2 = rho * rho;
  logdetSpaceUnscaled = 0.0;
  for(i = 0; i < g->nSpace; i++){
    m = g->parentCount[i];
    denom = 1.0 + ((double)m - 1.0) * rho2;
    if(denom <= 0.0 || !std::isfinite(denom))
      return 0;
    logdetSpaceUnscaled += std::log(denom) - std::log(1.0 - rho2);
  }

  *out = -((double)g->nSpace * (double)g->nTime) * std::log(term->sigmaSq) +
    ((double)g->nTime) * logdetSpaceUnscaled +
    ((double)g->nSpace) * logdetTime;
  return 1;
}

double logdet_dagar_time_Q(TermState *term, GraphState *g)
{
  double out;

  if(!logdet_dagar_time_Q_try(term, g, &out))
    Rf_error("DAGAR-time precision is not positive definite for current parameters");
  return out;
}

double theta_forward(double x, double lower, double upper)
{
  return std::log((x - lower) / (upper - x));
}

double theta_inverse(double z, double lower, double upper)
{
  double ez;
  ez = std::exp(z);
  return lower + (upper - lower) * ez / (1.0 + ez);
}

double theta_log_jacobian(double x, double lower, double upper)
{
  return std::log(x - lower) + std::log(upper - x);
}

/*==========================================================================*/
/* log determinant of structured latent-process precision blocks             */
/*==========================================================================*/

double logdet_Qw_blocks(SamplerState *s)
{
  int t, i;
  double out;
  TermState *term;
  GraphState *g;

  out = 0.0;
  for(t = 0; t < s->nTerms; t++){
    term = s->terms + t;
    g = s->graphs + term->graphIndex;

    if(g->type == GRAPH_NNGP){
      out -= term->nNode * std::log(term->sigmaSq);
      for(i = 0; i < term->nNode; i++)
        out += std::log(term->F[i]);
    } else if(g->type == GRAPH_GP){
      out += term->logDetQ;
    } else if(g->type == GRAPH_AR1){
      double phi, den;
      int m;
      
      phi = term->theta[0];
      den = 1.0 - phi * phi;
      m = term->nNode;
      
      out -= m * std::log(term->sigmaSq);
      out -= (m - 1) * std::log(den);
    } else if(g->type == GRAPH_CAR){
      out += logdet_car_Q_sparse(s, term, g);
    } else if(g->type == GRAPH_CAR_TIME){
      out += logdet_car_time_Q(s, term, g);
    } else if(g->type == GRAPH_DAGAR){
      out += logdet_dagar_Q(term, g);
    } else if(g->type == GRAPH_DAGAR_TIME){
      out += logdet_dagar_time_Q(term, g);
    } else {
      Rf_error("unsupported graph type in logdet_Qw_blocks");
    }
  }
  return out;
}

int logdet_Qw_blocks_try(SamplerState *s, double *out)
{
  int t, i;
  double val;
  TermState *term;
  GraphState *g;

  *out = 0.0;
  for(t = 0; t < s->nTerms; t++){
    term = s->terms + t;
    g = s->graphs + term->graphIndex;

    if(g->type == GRAPH_NNGP){
      *out -= term->nNode * std::log(term->sigmaSq);
      for(i = 0; i < term->nNode; i++)
        *out += std::log(term->F[i]);
    } else if(g->type == GRAPH_GP){
      *out += term->logDetQ;
    } else if(g->type == GRAPH_AR1){
      double phi, den;
      int m;

      phi = term->theta[0];
      if(std::fabs(phi) >= 1.0)
        return 0;
      den = 1.0 - phi * phi;
      m = term->nNode;

      *out -= m * std::log(term->sigmaSq);
      *out -= (m - 1) * std::log(den);
    } else if(g->type == GRAPH_CAR){
      if(!logdet_car_Q_sparse_try(s, term, g, &val))
        return 0;
      *out += val;
    } else if(g->type == GRAPH_CAR_TIME){
      if(!logdet_car_time_Q_try(s, term, g, &val))
        return 0;
      *out += val;
    } else if(g->type == GRAPH_DAGAR){
      if(!logdet_dagar_Q_try(term, g, &val))
        return 0;
      *out += val;
    } else if(g->type == GRAPH_DAGAR_TIME){
      if(!logdet_dagar_time_Q_try(term, g, &val))
        return 0;
      *out += val;
    } else {
      Rf_error("unsupported graph type in logdet_Qw_blocks_try");
    }
  }
  return 1;
}
