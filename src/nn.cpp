#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "nn.h"
#include "utils.h"

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>


inline void getNNIndx(int i, int m, int &iNNIndx, int &iNN){

  if(i == 0){
    iNNIndx = 0;
    iNN = 0;
  }
  else if(i < m){
    iNN = i;
    iNNIndx = i*(i-1)/2;
  }
  else{
    iNN = m;
    iNNIndx = m*(m-1)/2 + (i-m)*m;
  }
}

static void consider_neighbor(double *best_d, int *best_i, int m, double d, int j)
{
  int pos;

  if(m <= 0)
    return;

  if(d > best_d[m - 1] || (d == best_d[m - 1] && j >= best_i[m - 1]))
    return;

  pos = m - 1;
  while(pos > 0 &&
        (d < best_d[pos - 1] ||
         (d == best_d[pos - 1] && j < best_i[pos - 1]))){
    best_d[pos] = best_d[pos - 1];
    best_i[pos] = best_i[pos - 1];
    pos--;
  }

  best_d[pos] = d;
  best_i[pos] = j;
}

static SEXP make_nn_output(SEXP nnIndx_r, SEXP nnDist_r, SEXP nnIndxLU_r)
{
  SEXP out, names;

  PROTECT(out = Rf_allocVector(VECSXP, 3));

  SET_VECTOR_ELT(out, 0, nnIndx_r);
  SET_VECTOR_ELT(out, 1, nnDist_r);
  SET_VECTOR_ELT(out, 2, nnIndxLU_r);

  PROTECT(names = Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(names, 0, Rf_mkChar("nnIndx"));
  SET_STRING_ELT(names, 1, Rf_mkChar("nnDist"));
  SET_STRING_ELT(names, 2, Rf_mkChar("nnIndxLU"));

  Rf_setAttrib(out, R_NamesSymbol, names);

  UNPROTECT(2);
  return out;
}

///////////////////////////////////////////////////////////////////
// Brute force ordered NNGP neighbor search
///////////////////////////////////////////////////////////////////

SEXP mkNNIndxBrute(SEXP coords_r,
                   SEXP n_r,
                   SEXP m_r,
                   SEXP r_r,
                   SEXP nThreads_r){

  int i, j, k, iNNIndx, iNN;
  double d, diff;

  int n = INTEGER(n_r)[0];
  int m = INTEGER(m_r)[0];
  int r = INTEGER(r_r)[0];

  if(m > n)
    Rf_error("m must be less than or equal to n");

  if(r < 1)
    Rf_error("invalid coordinate dimension r");

  
  int nThreads = INTEGER(nThreads_r)[0];

  double *coords = REAL(coords_r);

#ifdef _OPENMP
  omp_set_num_threads(nThreads);
#else
  if(nThreads > 1){
    Rf_warning("OpenMP not available, using single thread.");
    nThreads = 1;
  }
#endif

  ////////////////////////////////////////////////////////////////
  // number of neighbor index entries
  ////////////////////////////////////////////////////////////////

  int nIndx = (n <= m) ? (n * (n - 1)) / 2 : (m*(m+1))/2 + (n-m-1)*m;

  ////////////////////////////////////////////////////////////////
  // allocate output
  ////////////////////////////////////////////////////////////////

  SEXP nnIndx_r   = PROTECT(Rf_allocVector(INTSXP, nIndx));
  SEXP nnDist_r   = PROTECT(Rf_allocVector(REALSXP, nIndx));
  SEXP nnIndxLU_r = PROTECT(Rf_allocVector(INTSXP, 2*n));

  int *nnIndx = INTEGER(nnIndx_r);
  double *nnDist = REAL(nnDist_r);
  int *nnIndxLU = INTEGER(nnIndxLU_r);

  ////////////////////////////////////////////////////////////////
  // initialize neighbor distances
  ////////////////////////////////////////////////////////////////

  for(i = 0; i < nIndx; i++){
    nnDist[i] = std::numeric_limits<double>::infinity();
    nnIndx[i] = -1;
  }

  ////////////////////////////////////////////////////////////////
  // neighbor search
  ////////////////////////////////////////////////////////////////

#ifdef _OPENMP
#pragma omp parallel for private(j,k,iNNIndx,iNN,d,diff)
#endif
  for(i = 0; i < n; i++){
    
    getNNIndx(i, m, iNNIndx, iNN);
    
    nnIndxLU[i]   = iNNIndx;
    nnIndxLU[n+i] = iNN;
    
    if(i == 0)
      continue;
    
    for(j = 0; j < i; j++){
      
      d = 0.0;
      
      for(k = 0; k < r; k++){
	diff = coords[k*n + i] - coords[k*n + j];
	d += diff*diff;
      }
      
      d = sqrt(d);

      if(iNN > 0 && d < nnDist[iNNIndx + iNN - 1]){
        nnDist[iNNIndx + iNN - 1] = d;
        nnIndx[iNNIndx + iNN - 1] = j;

        rsort_with_index(
                         &nnDist[iNNIndx],
                         &nnIndx[iNNIndx],
                         iNN
                         );
      }
    }
    
  }

  ////////////////////////////////////////////////////////////////
  // return result to R
  ////////////////////////////////////////////////////////////////

  SEXP out = PROTECT(make_nn_output(nnIndx_r, nnDist_r, nnIndxLU_r));

  UNPROTECT(4);

  return out;
}

///////////////////////////////////////////////////////////////////
// Exact ordered NNGP neighbor search with history-restricted k-d tree
///////////////////////////////////////////////////////////////////

struct KDNode {
  int point;
  int axis;
  int left;
  int right;
  int minOrder;
};

static int build_kd_tree(std::vector<int> &idx,
                         int lo,
                         int hi,
                         const double *coords,
                         int n,
                         int r,
                         std::vector<KDNode> &nodes,
                         std::vector<double> &bboxMin,
                         std::vector<double> &bboxMax)
{
  int i, k, mid, axis, nodeIndex, left, right, point, minOrder;
  double spread, bestSpread;

  if(lo >= hi)
    return -1;

  axis = 0;
  bestSpread = -1.0;
  for(k = 0; k < r; k++){
    double mn = coords[idx[lo] + n * k];
    double mx = mn;
    for(i = lo + 1; i < hi; i++){
      double val = coords[idx[i] + n * k];
      if(val < mn) mn = val;
      if(val > mx) mx = val;
    }
    spread = mx - mn;
    if(spread > bestSpread){
      bestSpread = spread;
      axis = k;
    }
  }

  mid = lo + (hi - lo) / 2;
  std::nth_element(
    idx.begin() + lo,
    idx.begin() + mid,
    idx.begin() + hi,
    [coords, n, axis](int a, int b){
      double da = coords[a + n * axis];
      double db = coords[b + n * axis];
      if(da < db) return true;
      if(da > db) return false;
      return a < b;
    }
  );

  point = idx[mid];
  nodeIndex = (int)nodes.size();
  nodes.push_back({point, axis, -1, -1, point});
  bboxMin.resize((size_t)(nodeIndex + 1) * (size_t)r);
  bboxMax.resize((size_t)(nodeIndex + 1) * (size_t)r);

  left = build_kd_tree(idx, lo, mid, coords, n, r, nodes, bboxMin, bboxMax);
  right = build_kd_tree(idx, mid + 1, hi, coords, n, r, nodes, bboxMin, bboxMax);

  nodes[nodeIndex].left = left;
  nodes[nodeIndex].right = right;

  minOrder = point;
  for(k = 0; k < r; k++){
    double mn = coords[point + n * k];
    double mx = mn;
    if(left >= 0){
      if(nodes[left].minOrder < minOrder)
        minOrder = nodes[left].minOrder;
      if(bboxMin[(size_t)left * (size_t)r + (size_t)k] < mn)
        mn = bboxMin[(size_t)left * (size_t)r + (size_t)k];
      if(bboxMax[(size_t)left * (size_t)r + (size_t)k] > mx)
        mx = bboxMax[(size_t)left * (size_t)r + (size_t)k];
    }
    if(right >= 0){
      if(nodes[right].minOrder < minOrder)
        minOrder = nodes[right].minOrder;
      if(bboxMin[(size_t)right * (size_t)r + (size_t)k] < mn)
        mn = bboxMin[(size_t)right * (size_t)r + (size_t)k];
      if(bboxMax[(size_t)right * (size_t)r + (size_t)k] > mx)
        mx = bboxMax[(size_t)right * (size_t)r + (size_t)k];
    }
    bboxMin[(size_t)nodeIndex * (size_t)r + (size_t)k] = mn;
    bboxMax[(size_t)nodeIndex * (size_t)r + (size_t)k] = mx;
  }
  nodes[nodeIndex].minOrder = minOrder;

  return nodeIndex;
}

static double bbox_distance2(int nodeIndex,
                             int target,
                             const double *coords,
                             int n,
                             int r,
                             const std::vector<double> &bboxMin,
                             const std::vector<double> &bboxMax)
{
  int k;
  double d2, diff, x, mn, mx;

  d2 = 0.0;
  for(k = 0; k < r; k++){
    x = coords[target + n * k];
    mn = bboxMin[(size_t)nodeIndex * (size_t)r + (size_t)k];
    mx = bboxMax[(size_t)nodeIndex * (size_t)r + (size_t)k];
    if(x < mn){
      diff = mn - x;
      d2 += diff * diff;
    } else if(x > mx){
      diff = x - mx;
      d2 += diff * diff;
    }
  }

  return d2;
}

static void query_kd_tree(int nodeIndex,
                          int target,
                          int m,
                          const double *coords,
                          int n,
                          int r,
                          const std::vector<KDNode> &nodes,
                          const std::vector<double> &bboxMin,
                          const std::vector<double> &bboxMax,
                          double *best_d2,
                          int *best_i)
{
  const KDNode *node;
  int k, point, nearChild, farChild;
  double diff, d2, splitDiff, farBound;

  if(nodeIndex < 0)
    return;

  node = &nodes[nodeIndex];

  if(node->minOrder >= target)
    return;

  if(best_i[m - 1] >= 0 && bbox_distance2(nodeIndex, target, coords, n, r, bboxMin, bboxMax) > best_d2[m - 1])
    return;

  point = node->point;
  if(point < target){
    d2 = 0.0;
    for(k = 0; k < r; k++){
      diff = coords[target + n * k] - coords[point + n * k];
      d2 += diff * diff;
    }
    consider_neighbor(best_d2, best_i, m, d2, point);
  }

  splitDiff = coords[target + n * node->axis] - coords[point + n * node->axis];
  if(splitDiff <= 0.0){
    nearChild = node->left;
    farChild = node->right;
  } else {
    nearChild = node->right;
    farChild = node->left;
  }

  query_kd_tree(nearChild, target, m, coords, n, r, nodes, bboxMin, bboxMax, best_d2, best_i);

  if(farChild >= 0){
    farBound = bbox_distance2(farChild, target, coords, n, r, bboxMin, bboxMax);
    if(best_i[m - 1] < 0 || farBound <= best_d2[m - 1])
      query_kd_tree(farChild, target, m, coords, n, r, nodes, bboxMin, bboxMax, best_d2, best_i);
  }
}

SEXP mkNNIndx(SEXP coords_r,
              SEXP n_r,
              SEXP m_r,
              SEXP r_r,
              SEXP nThreads_r){

  int i, k, root, iNNIndx, iNN;
  int n = INTEGER(n_r)[0];
  int m = INTEGER(m_r)[0];
  int r = INTEGER(r_r)[0];
  int nThreads = INTEGER(nThreads_r)[0];
  double *coords = REAL(coords_r);

  if(m > n)
    Rf_error("m must be less than or equal to n");

  if(r < 1)
    Rf_error("invalid coordinate dimension r");

#ifdef _OPENMP
  omp_set_num_threads(nThreads);
#else
  if(nThreads > 1){
    Rf_warning("OpenMP not available, using single thread.");
    nThreads = 1;
  }
#endif

  int nIndx = (n <= m) ? (n * (n - 1)) / 2 : (m*(m+1))/2 + (n-m-1)*m;

  SEXP nnIndx_r   = PROTECT(Rf_allocVector(INTSXP, nIndx));
  SEXP nnDist_r   = PROTECT(Rf_allocVector(REALSXP, nIndx));
  SEXP nnIndxLU_r = PROTECT(Rf_allocVector(INTSXP, 2*n));

  int *nnIndx = INTEGER(nnIndx_r);
  double *nnDist = REAL(nnDist_r);
  int *nnIndxLU = INTEGER(nnIndxLU_r);

  for(i = 0; i < nIndx; i++){
    nnDist[i] = std::numeric_limits<double>::infinity();
    nnIndx[i] = -1;
  }

  std::vector<int> idx((size_t)n);
  for(i = 0; i < n; i++)
    idx[(size_t)i] = i;

  std::vector<KDNode> nodes;
  std::vector<double> bboxMin;
  std::vector<double> bboxMax;
  nodes.reserve((size_t)n);
  bboxMin.reserve((size_t)n * (size_t)r);
  bboxMax.reserve((size_t)n * (size_t)r);
  root = build_kd_tree(idx, 0, n, coords, n, r, nodes, bboxMin, bboxMax);

#ifdef _OPENMP
#pragma omp parallel for private(k,iNNIndx,iNN) schedule(static)
#endif
  for(i = 0; i < n; i++){
    getNNIndx(i, m, iNNIndx, iNN);

    nnIndxLU[i]   = iNNIndx;
    nnIndxLU[n+i] = iNN;

    if(i == 0)
      continue;

    query_kd_tree(root, i, iNN, coords, n, r, nodes, bboxMin, bboxMax,
                  &nnDist[iNNIndx], &nnIndx[iNNIndx]);

    for(k = 0; k < iNN; k++)
      nnDist[iNNIndx + k] = sqrt(nnDist[iNNIndx + k]);
  }

  SEXP out = PROTECT(make_nn_output(nnIndx_r, nnDist_r, nnIndxLU_r));

  UNPROTECT(4);

  return out;
}
