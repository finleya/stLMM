#include "stLMM_internal.h"

#ifdef _OPENMP
#include <omp.h>
#endif

static double coord_distance(const double *a,
                             int ia,
                             const double *b,
                             int ib,
                             int na,
                             int nb,
                             int dim,
                             int distanceMode,
                             double *u)
{
  int k, d_space;
  double ss, diff;

  *u = 0.0;
  ss = 0.0;

  if(distanceMode == COR_SINGLE){
    for(k = 0; k < dim; k++){
      diff = a[ia + na * k] - b[ib + nb * k];
      ss += diff * diff;
    }
    return std::sqrt(ss);
  }

  d_space = dim - 1;
  for(k = 0; k < d_space; k++){
    diff = a[ia + na * k] - b[ib + nb * k];
    ss += diff * diff;
  }
  *u = std::fabs(a[ia + na * (dim - 1)] - b[ib + nb * (dim - 1)]);

  return std::sqrt(ss);
}

extern "C" SEXP stLMM_nngp_prediction_neighbors(SEXP support_r,
                                                SEXP new_coords_r,
                                                SEXP m_r,
                                                SEXP cov_model_r,
                                                SEXP st_scale_r,
                                                SEXP n_omp_threads_r)
{
  int n_fit, n_new, dim, m, k, i, j, r, pos, n_threads;
  int cov_model_index, has_time;
  double st_scale, diff, d2;
  const double *support, *new_coords;
  double *scratch_d;
  int *scratch_i;
  int *out;
  SEXP out_r;
  const CorModelInfo *model_info;

  if(!Rf_isMatrix(support_r) || !Rf_isReal(support_r))
    Rf_error("support must be a numeric matrix");
  if(!Rf_isMatrix(new_coords_r) || !Rf_isReal(new_coords_r))
    Rf_error("new_coords must be a numeric matrix");

  n_fit = INTEGER(Rf_getAttrib(support_r, R_DimSymbol))[0];
  dim = INTEGER(Rf_getAttrib(support_r, R_DimSymbol))[1];
  n_new = INTEGER(Rf_getAttrib(new_coords_r, R_DimSymbol))[0];
  if(INTEGER(Rf_getAttrib(new_coords_r, R_DimSymbol))[1] != dim)
    Rf_error("new_coords and support must have the same number of columns");
  if(n_fit < 1)
    Rf_error("support must have at least one row");
  if(n_new < 1)
    Rf_error("new_coords must have at least one row");

  m = as_int_scalar(m_r, "m");
  if(m < 1)
    Rf_error("m must be positive");
  if(m > n_fit)
    m = n_fit;

  cov_model_index = get_cor_model_index(as_char_scalar(cov_model_r, "cov_model"));
  if(cov_model_index < 0)
    Rf_error("unknown covariance model");
  model_info = get_cor_model_info(cov_model_index);
  if(model_info == NULL)
    Rf_error("invalid covariance model metadata");
  has_time = model_info->distanceMode != COR_SINGLE;

  if(!Rf_isReal(st_scale_r) && !Rf_isInteger(st_scale_r))
    Rf_error("st_scale must be numeric");
  if(Rf_length(st_scale_r) != 1)
    Rf_error("st_scale must be a scalar");
  st_scale = Rf_asReal(st_scale_r);
  if(!R_FINITE(st_scale) || st_scale <= 0.0)
    Rf_error("st_scale must be a positive finite scalar");

  n_threads = as_int_scalar(n_omp_threads_r, "n_omp_threads");
  if(n_threads < 1)
    n_threads = 1;

  support = REAL(support_r);
  new_coords = REAL(new_coords_r);

  PROTECT(out_r = Rf_allocMatrix(INTSXP, n_new, m));
  out = INTEGER(out_r);
  scratch_d = (double*)R_alloc((size_t)n_threads * (size_t)m, sizeof(double));
  scratch_i = (int*)R_alloc((size_t)n_threads * (size_t)m, sizeof(int));

#ifdef _OPENMP
#pragma omp parallel for private(i, j, k, r, pos, diff, d2) num_threads(n_threads) schedule(static)
#endif
  for(i = 0; i < n_new; i++){
    int tid;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif

    double *best_d = scratch_d + (size_t)tid * (size_t)m;
    int *best_i = scratch_i + (size_t)tid * (size_t)m;

    for(k = 0; k < m; k++){
      best_d[k] = R_PosInf;
      best_i[k] = -1;
    }

    for(j = 0; j < n_fit; j++){
      d2 = 0.0;
      for(k = 0; k < dim; k++){
        diff = new_coords[i + n_new * k] - support[j + n_fit * k];
        if(has_time && k == dim - 1)
          diff *= st_scale;
        d2 += diff * diff;
      }

      if(d2 > best_d[m - 1] || (d2 == best_d[m - 1] && j > best_i[m - 1]))
        continue;

      pos = m - 1;
      while(pos > 0 && (d2 < best_d[pos - 1] || (d2 == best_d[pos - 1] && j < best_i[pos - 1]))){
        best_d[pos] = best_d[pos - 1];
        best_i[pos] = best_i[pos - 1];
        pos--;
      }
      best_d[pos] = d2;
      best_i[pos] = j;
    }

    for(r = 0; r < m; r++)
      out[i + n_new * r] = best_i[r] + 1;
  }

  UNPROTECT(1);
  return out_r;
}

extern "C" SEXP stLMM_predict_nngp_joint_false(SEXP support_r,
                                               SEXP new_coords_r,
                                               SEXP neighbor_index_r,
                                               SEXP w_fit_r,
                                               SEXP sigma_sq_r,
                                               SEXP theta_r,
                                               SEXP cov_model_r,
                                               SEXP n_omp_threads_r)
{
  int n_fit, n_new, dim, n_draw, theta_dim, m, n_threads, total;
  int i, j, k, l, idx, draw, node, tid, info;
  int cov_model_index, distance_mode, error_code, error_info;
  double one, h, u, mean, rBr, var;
  const double *support, *new_coords, *w_fit, *sigma_sq, *theta;
  const int *neighbor_index;
  double *out, *z, *scratch_C, *scratch_c, *scratch_B, *scratch_theta;
  SEXP out_r, z_r;
  corFunPtr cor_fun;
  const CorModelInfo *model_info;

  if(!Rf_isMatrix(support_r) || !Rf_isReal(support_r))
    Rf_error("support must be a numeric matrix");
  if(!Rf_isMatrix(new_coords_r) || !Rf_isReal(new_coords_r))
    Rf_error("new_coords must be a numeric matrix");
  if(!Rf_isMatrix(neighbor_index_r) || !Rf_isInteger(neighbor_index_r))
    Rf_error("neighbor_index must be an integer matrix");
  if(!Rf_isMatrix(w_fit_r) || !Rf_isReal(w_fit_r))
    Rf_error("w_fit must be a numeric matrix");
  if(!Rf_isReal(sigma_sq_r))
    Rf_error("sigma_sq must be numeric");
  if(!Rf_isMatrix(theta_r) || !Rf_isReal(theta_r))
    Rf_error("theta must be a numeric matrix");

  n_fit = INTEGER(Rf_getAttrib(support_r, R_DimSymbol))[0];
  dim = INTEGER(Rf_getAttrib(support_r, R_DimSymbol))[1];
  n_new = INTEGER(Rf_getAttrib(new_coords_r, R_DimSymbol))[0];
  if(INTEGER(Rf_getAttrib(new_coords_r, R_DimSymbol))[1] != dim)
    Rf_error("new_coords and support must have the same number of columns");

  if(INTEGER(Rf_getAttrib(w_fit_r, R_DimSymbol))[1] != n_fit)
    Rf_error("w_fit column count must match support rows");
  n_draw = INTEGER(Rf_getAttrib(w_fit_r, R_DimSymbol))[0];

  if(Rf_length(sigma_sq_r) != n_draw)
    Rf_error("sigma_sq length must match w_fit rows");
  if(INTEGER(Rf_getAttrib(theta_r, R_DimSymbol))[0] != n_draw)
    Rf_error("theta row count must match w_fit rows");
  theta_dim = INTEGER(Rf_getAttrib(theta_r, R_DimSymbol))[1];

  if(INTEGER(Rf_getAttrib(neighbor_index_r, R_DimSymbol))[0] != n_new)
    Rf_error("neighbor_index row count must match new_coords rows");
  m = INTEGER(Rf_getAttrib(neighbor_index_r, R_DimSymbol))[1];
  if(m <= 0)
    Rf_error("neighbor_index must have at least one column");

  cov_model_index = get_cor_model_index(as_char_scalar(cov_model_r, "cov_model"));
  if(cov_model_index < 0)
    Rf_error("unknown covariance model");
  model_info = get_cor_model_info(cov_model_index);
  if(model_info == NULL || model_info->fun == NULL)
    Rf_error("invalid covariance model metadata");
  if(model_info->nTheta != theta_dim)
    Rf_error("theta column count does not match covariance model");
  cor_fun = model_info->fun;
  distance_mode = model_info->distanceMode;

  n_threads = as_int_scalar(n_omp_threads_r, "n_omp_threads");
  if(n_threads < 1)
    n_threads = 1;

  support = REAL(support_r);
  new_coords = REAL(new_coords_r);
  neighbor_index = INTEGER(neighbor_index_r);
  w_fit = REAL(w_fit_r);
  sigma_sq = REAL(sigma_sq_r);
  theta = REAL(theta_r);

  PROTECT(out_r = Rf_allocMatrix(REALSXP, n_draw, n_new));
  PROTECT(z_r = Rf_allocVector(REALSXP, (R_xlen_t)n_draw * (R_xlen_t)n_new));
  out = REAL(out_r);
  z = REAL(z_r);

  GetRNGstate();
  for(i = 0; i < n_draw * n_new; i++)
    z[i] = norm_rand();
  PutRNGstate();

  scratch_C = (double*)R_alloc((size_t)n_threads * (size_t)m * (size_t)m, sizeof(double));
  scratch_c = (double*)R_alloc((size_t)n_threads * (size_t)m, sizeof(double));
  scratch_B = (double*)R_alloc((size_t)n_threads * (size_t)m, sizeof(double));
  scratch_theta = (double*)R_alloc((size_t)n_threads * (size_t)theta_dim, sizeof(double));

  total = n_draw * n_new;
  error_code = 0;
  error_info = 0;

  for(i = 0; i < n_new * m; i++){
    if(neighbor_index[i] < 1 || neighbor_index[i] > n_fit)
      Rf_error("neighbor_index contains an out-of-range index");
  }

#ifdef _OPENMP
#pragma omp parallel for private(idx, draw, node, tid, k, l, j, h, u, info, mean, rBr, var) num_threads(n_threads) schedule(static)
#endif
  for(idx = 0; idx < total; idx++){
    draw = idx % n_draw;
    node = idx / n_draw;

#ifdef _OPENMP
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif

    double *C = scratch_C + (size_t)tid * (size_t)m * (size_t)m;
    double *c = scratch_c + (size_t)tid * (size_t)m;
    double *B = scratch_B + (size_t)tid * (size_t)m;
    double *theta_i = scratch_theta + (size_t)tid * (size_t)theta_dim;

    for(k = 0; k < theta_dim; k++)
      theta_i[k] = theta[draw + n_draw * k];

    for(k = 0; k < m * m; k++)
      C[k] = 0.0;

    for(k = 0; k < m; k++){
      int nbr_k = neighbor_index[node + n_new * k] - 1;

      h = coord_distance(new_coords, node, support, nbr_k,
                         n_new, n_fit, dim, distance_mode, &u);
      c[k] = cor_fun(theta_i, h, u, distance_mode == COR_SPACE_TIME ? dim - 1 : dim);

      for(l = 0; l <= k; l++){
        int nbr_l = neighbor_index[node + n_new * l] - 1;

        h = coord_distance(support, nbr_k, support, nbr_l,
                           n_fit, n_fit, dim, distance_mode, &u);
        C[k + m * l] = cor_fun(theta_i, h, u, distance_mode == COR_SPACE_TIME ? dim - 1 : dim);
      }
    }

    info = small_chol_lower(C, m);
    if(info != 0){
#ifdef _OPENMP
#pragma omp critical
#endif
      {
        if(error_code == 0){
          error_code = 1;
          error_info = info;
        }
      }
      out[draw + n_draw * node] = NA_REAL;
      continue;
    }

    info = small_chol_solve_lower(C, c, B, m);
    if(info != 0){
#ifdef _OPENMP
#pragma omp critical
#endif
      {
        if(error_code == 0){
          error_code = 2;
          error_info = info;
        }
      }
      out[draw + n_draw * node] = NA_REAL;
      continue;
    }

    mean = 0.0;
    for(k = 0; k < m; k++){
      j = neighbor_index[node + n_new * k] - 1;
      mean += B[k] * w_fit[draw + n_draw * j];
    }

    rBr = 0.0;
    for(k = 0; k < m; k++)
      rBr += c[k] * B[k];
    var = sigma_sq[draw] * (1.0 - rBr);
    if(var < 0.0)
      var = 0.0;

    out[draw + n_draw * node] = mean + std::sqrt(var) * z[draw + n_draw * node];
  }

  if(error_code == 1)
    Rf_error("stLMM_predict_nngp_joint_false: local Cholesky failed with info=%d", error_info);
  if(error_code == 2)
    Rf_error("stLMM_predict_nngp_joint_false: local Cholesky solve failed with info=%d", error_info);

  UNPROTECT(2);
  return out_r;
}

extern "C" SEXP stLMM_predict_nngp_vecchia_joint(SEXP coords_all_r,
                                                 SEXP n_fit_r,
                                                 SEXP neighbor_index_r,
                                                 SEXP neighbor_count_r,
                                                 SEXP w_fit_r,
                                                 SEXP sigma_sq_r,
                                                 SEXP theta_r,
                                                 SEXP cov_model_r,
                                                 SEXP n_omp_threads_r)
{
  int n_all, n_fit, n_pred, dim, n_draw, theta_dim, m_max, n_threads;
  int draw, node, k, l, j, idx, tid, info, total;
  int cov_model_index, distance_mode, error_code, error_info;
  double h, u, mean, rBr, var;
  const double *coords_all, *w_fit, *sigma_sq, *theta;
  const int *neighbor_index, *neighbor_count;
  double *out, *z, *scratch_C, *scratch_c, *scratch_B, *scratch_theta, *scratch_w;
  SEXP out_r, z_r;
  corFunPtr cor_fun;
  const CorModelInfo *model_info;

  if(!Rf_isMatrix(coords_all_r) || !Rf_isReal(coords_all_r))
    Rf_error("coords_all must be a numeric matrix");
  if(!Rf_isMatrix(neighbor_index_r) || !Rf_isInteger(neighbor_index_r))
    Rf_error("neighbor_index must be an integer matrix");
  if(!Rf_isInteger(neighbor_count_r))
    Rf_error("neighbor_count must be an integer vector");
  if(!Rf_isMatrix(w_fit_r) || !Rf_isReal(w_fit_r))
    Rf_error("w_fit must be a numeric matrix");
  if(!Rf_isReal(sigma_sq_r))
    Rf_error("sigma_sq must be numeric");
  if(!Rf_isMatrix(theta_r) || !Rf_isReal(theta_r))
    Rf_error("theta must be a numeric matrix");

  n_all = INTEGER(Rf_getAttrib(coords_all_r, R_DimSymbol))[0];
  dim = INTEGER(Rf_getAttrib(coords_all_r, R_DimSymbol))[1];
  n_fit = as_int_scalar(n_fit_r, "n_fit");
  if(n_fit < 1 || n_fit >= n_all)
    Rf_error("n_fit must be positive and less than coords_all rows");
  n_pred = n_all - n_fit;

  if(INTEGER(Rf_getAttrib(w_fit_r, R_DimSymbol))[1] != n_fit)
    Rf_error("w_fit column count must match n_fit");
  n_draw = INTEGER(Rf_getAttrib(w_fit_r, R_DimSymbol))[0];

  if(Rf_length(sigma_sq_r) != n_draw)
    Rf_error("sigma_sq length must match w_fit rows");
  if(INTEGER(Rf_getAttrib(theta_r, R_DimSymbol))[0] != n_draw)
    Rf_error("theta row count must match w_fit rows");
  theta_dim = INTEGER(Rf_getAttrib(theta_r, R_DimSymbol))[1];

  if(INTEGER(Rf_getAttrib(neighbor_index_r, R_DimSymbol))[0] != n_pred)
    Rf_error("neighbor_index row count must match prediction node count");
  m_max = INTEGER(Rf_getAttrib(neighbor_index_r, R_DimSymbol))[1];
  if(m_max <= 0)
    Rf_error("neighbor_index must have at least one column");
  if(Rf_length(neighbor_count_r) != n_pred)
    Rf_error("neighbor_count length must match prediction node count");

  cov_model_index = get_cor_model_index(as_char_scalar(cov_model_r, "cov_model"));
  if(cov_model_index < 0)
    Rf_error("unknown covariance model");
  model_info = get_cor_model_info(cov_model_index);
  if(model_info == NULL || model_info->fun == NULL)
    Rf_error("invalid covariance model metadata");
  if(model_info->nTheta != theta_dim)
    Rf_error("theta column count does not match covariance model");
  cor_fun = model_info->fun;
  distance_mode = model_info->distanceMode;

  n_threads = as_int_scalar(n_omp_threads_r, "n_omp_threads");
  if(n_threads < 1)
    n_threads = 1;

  coords_all = REAL(coords_all_r);
  neighbor_index = INTEGER(neighbor_index_r);
  neighbor_count = INTEGER(neighbor_count_r);
  w_fit = REAL(w_fit_r);
  sigma_sq = REAL(sigma_sq_r);
  theta = REAL(theta_r);

  for(node = 0; node < n_pred; node++){
    int count = neighbor_count[node];
    if(count < 1 || count > m_max)
      Rf_error("neighbor_count contains an invalid value");
    for(k = 0; k < count; k++){
      int nbr = neighbor_index[node + n_pred * k] - 1;
      if(nbr < 0 || nbr >= n_fit + node)
        Rf_error("neighbor_index violates Vecchia prediction history ordering");
    }
    for(k = count; k < m_max; k++){
      if(neighbor_index[node + n_pred * k] != 0)
        Rf_error("neighbor_index padding must be zero");
    }
  }

  PROTECT(out_r = Rf_allocMatrix(REALSXP, n_draw, n_pred));
  PROTECT(z_r = Rf_allocVector(REALSXP, (R_xlen_t)n_draw * (R_xlen_t)n_pred));
  out = REAL(out_r);
  z = REAL(z_r);

  GetRNGstate();
  for(idx = 0; idx < n_draw * n_pred; idx++)
    z[idx] = norm_rand();
  PutRNGstate();

  scratch_C = (double*)R_alloc((size_t)n_threads * (size_t)m_max * (size_t)m_max, sizeof(double));
  scratch_c = (double*)R_alloc((size_t)n_threads * (size_t)m_max, sizeof(double));
  scratch_B = (double*)R_alloc((size_t)n_threads * (size_t)m_max, sizeof(double));
  scratch_theta = (double*)R_alloc((size_t)n_threads * (size_t)theta_dim, sizeof(double));
  scratch_w = (double*)R_alloc((size_t)n_threads * (size_t)n_all, sizeof(double));

  error_code = 0;
  error_info = 0;
  total = n_draw;

#ifdef _OPENMP
#pragma omp parallel for private(draw, node, tid, k, l, j, h, u, info, mean, rBr, var) num_threads(n_threads) schedule(static)
#endif
  for(draw = 0; draw < total; draw++){
#ifdef _OPENMP
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif

    double *C = scratch_C + (size_t)tid * (size_t)m_max * (size_t)m_max;
    double *c = scratch_c + (size_t)tid * (size_t)m_max;
    double *B = scratch_B + (size_t)tid * (size_t)m_max;
    double *theta_i = scratch_theta + (size_t)tid * (size_t)theta_dim;
    double *w_all = scratch_w + (size_t)tid * (size_t)n_all;

    for(k = 0; k < n_fit; k++)
      w_all[k] = w_fit[draw + n_draw * k];
    for(k = n_fit; k < n_all; k++)
      w_all[k] = NA_REAL;
    for(k = 0; k < theta_dim; k++)
      theta_i[k] = theta[draw + n_draw * k];

    for(node = 0; node < n_pred; node++){
      int count = neighbor_count[node];
      int target = n_fit + node;

      for(k = 0; k < count * count; k++)
        C[k] = 0.0;

      for(k = 0; k < count; k++){
        int nbr_k = neighbor_index[node + n_pred * k] - 1;

        h = coord_distance(coords_all, target, coords_all, nbr_k,
                           n_all, n_all, dim, distance_mode, &u);
        c[k] = cor_fun(theta_i, h, u, distance_mode == COR_SPACE_TIME ? dim - 1 : dim);

        for(l = 0; l <= k; l++){
          int nbr_l = neighbor_index[node + n_pred * l] - 1;

          h = coord_distance(coords_all, nbr_k, coords_all, nbr_l,
                             n_all, n_all, dim, distance_mode, &u);
          C[k + count * l] = cor_fun(theta_i, h, u, distance_mode == COR_SPACE_TIME ? dim - 1 : dim);
        }
      }

      info = small_chol_lower(C, count);
      if(info != 0){
#ifdef _OPENMP
#pragma omp critical
#endif
        {
          if(error_code == 0){
            error_code = 1;
            error_info = info;
          }
        }
        out[draw + n_draw * node] = NA_REAL;
        continue;
      }

      info = small_chol_solve_lower(C, c, B, count);
      if(info != 0){
#ifdef _OPENMP
#pragma omp critical
#endif
        {
          if(error_code == 0){
            error_code = 2;
            error_info = info;
          }
        }
        out[draw + n_draw * node] = NA_REAL;
        continue;
      }

      mean = 0.0;
      for(k = 0; k < count; k++){
        j = neighbor_index[node + n_pred * k] - 1;
        mean += B[k] * w_all[j];
      }

      rBr = 0.0;
      for(k = 0; k < count; k++)
        rBr += c[k] * B[k];
      var = sigma_sq[draw] * (1.0 - rBr);
      if(var < 0.0)
        var = 0.0;

      w_all[target] = mean + std::sqrt(var) * z[draw + n_draw * node];
      out[draw + n_draw * node] = w_all[target];
    }
  }

  if(error_code == 1)
    Rf_error("stLMM_predict_nngp_vecchia_joint: local Cholesky failed with info=%d", error_info);
  if(error_code == 2)
    Rf_error("stLMM_predict_nngp_vecchia_joint: local Cholesky solve failed with info=%d", error_info);

  UNPROTECT(2);
  return out_r;
}
