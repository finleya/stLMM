#include "stLMM_internal.h"
#include <climits>

/* The joint latent order is (w, alpha), preserving the structured term offsets.
   The row mapping contains only supplied designs and membership indices. */
void init_joint_iid_mapping(SamplerState *s)
{
  int i, j, t, pos, nnz;
  int *next;
  size_t count;

  if(s->q <= 0 || s->qLatTotal > INT_MAX - s->q)
    Rf_error("invalid joint latent dimension");
  s->collapseIid = 1;
  s->qLatSystem = s->qLatTotal + s->q;
  count = (size_t)s->n * s->nTerms + s->Z_nnz;
  if(count > INT_MAX)
    Rf_error("joint latent mapping exceeds integer sparse index capacity");
  nnz = (int)count;
  s->Hp = (int*)R_alloc((size_t)s->n + 1, sizeof(int));
  s->Hj = (int*)R_alloc(nnz > 0 ? nnz : 1, sizeof(int));
  s->Hx = (double*)R_alloc(nnz > 0 ? nnz : 1, sizeof(double));
  next = (int*)R_alloc(s->n, sizeof(int));

  for(i = 0; i <= s->n; i++)
    s->Hp[i] = 0;
  for(i = 0; i < s->Z_nnz; i++)
    s->Hp[s->Zi[i] + 1]++;
  for(i = 0; i < s->n; i++){
    s->Hp[i + 1] += s->Hp[i] + s->nTerms;
    next[i] = s->Hp[i];
    for(t = 0; t < s->nTerms; t++){
      TermState *term = s->terms + t;
      pos = next[i]++;
      s->Hj[pos] = term->wOffset + term->map[i];
      s->Hx[pos] = term->scale[i];
    }
  }
  for(j = 0; j < s->q; j++){
    for(i = s->Zp[j]; i < s->Zp[j + 1]; i++){
      pos = next[s->Zi[i]]++;
      s->Hj[pos] = s->qLatTotal + j;
      s->Hx[pos] = s->Zx[i];
    }
  }
  s->latentSolveRhs = (double*)R_alloc(s->qLatSystem, sizeof(double));
  s->latentNoiseZ = (double*)R_alloc(s->qLatSystem, sizeof(double));
}

void append_joint_iid_pattern(PatternCol *cols, SamplerState *s)
{
  int i, a, b, block, k, g, row, col, off;

  for(i = 0; i < s->q; i++)
    pattern_append_lower(cols, s->qLatTotal + i, s->qLatTotal + i);
  for(block = 0; block < s->nRE; block++){
    k = s->reCovDim[block];
    if(k == 0)
      continue;
    for(g = 0; g < s->reBlockSize[block] / k; g++){
      off = s->qLatTotal + s->reOffset[block] + g * k;
      for(col = 0; col < k; col++)
        for(row = col; row < k; row++)
          pattern_append_lower(cols, off + row, off + col);
    }
  }
  for(i = 0; i < s->n; i++)
    for(a = s->Hp[i]; a < s->Hp[i + 1]; a++)
      for(b = s->Hp[i]; b <= a; b++)
        pattern_append_lower(cols, s->Hj[a], s->Hj[b]);
}

void build_joint_iid_caches(SamplerState *s)
{
  int i, a, b, block, k, g, row, col, off;
  size_t count, width, pos;

  count = 0;
  for(i = 0; i < s->n; i++){
    width = s->Hp[i + 1] - s->Hp[i];
    count += width * (width + 1) / 2;
  }
  s->HtHIdx = (int*)R_alloc(count > 0 ? count : 1, sizeof(int));
  pos = 0;
  for(i = 0; i < s->n; i++)
    for(a = s->Hp[i]; a < s->Hp[i + 1]; a++)
      for(b = s->Hp[i]; b <= a; b++)
        s->HtHIdx[pos++] = find_M_lat_entry_index(s, s->Hj[a], s->Hj[b]);

  s->iidDiagIdx = (int*)R_alloc(s->q, sizeof(int));
  for(i = 0; i < s->q; i++)
    s->iidDiagIdx[i] = find_M_lat_entry_index(s, s->qLatTotal + i, s->qLatTotal + i);
  s->iidCovIdx = (int**)R_alloc(s->nRE, sizeof(int*));
  for(block = 0; block < s->nRE; block++){
    k = s->reCovDim[block];
    s->iidCovIdx[block] = NULL;
    if(k == 0)
      continue;
    count = (size_t)s->reBlockSize[block] * (k + 1) / 2;
    s->iidCovIdx[block] = (int*)R_alloc(count, sizeof(int));
    pos = 0;
    for(g = 0; g < s->reBlockSize[block] / k; g++){
      off = s->qLatTotal + s->reOffset[block] + g * k;
      for(col = 0; col < k; col++)
        for(row = col; row < k; row++)
          s->iidCovIdx[block][pos++] = find_M_lat_entry_index(s, off + row, off + col);
    }
  }
}

/* Add D_alpha^{-1} and H' W H after the structured prior has been assembled. */
void assemble_joint_iid_into_M(SamplerState *s)
{
  int i, a, b, block, k, g, row, col;
  size_t pos;
  double *Mx, *P;

  Mx = (double*)s->M_lat->x;
  for(i = 0; i < s->q; i++){
    block = s->reBlockID[i];
    if(s->reCovDim[block] == 0)
      Mx[s->iidDiagIdx[i]] += 1.0 / s->sigmaSqRE[block];
  }
  for(block = 0; block < s->nRE; block++){
    k = s->reCovDim[block];
    if(k == 0)
      continue;
    P = s->reCovPrecision[block];
    pos = 0;
    for(g = 0; g < s->reBlockSize[block] / k; g++)
      for(col = 0; col < k; col++)
        for(row = col; row < k; row++)
          Mx[s->iidCovIdx[block][pos++]] += P[row + k * col];
  }
  pos = 0;
  for(i = 0; i < s->n; i++)
    for(a = s->Hp[i]; a < s->Hp[i + 1]; a++)
      for(b = s->Hp[i]; b <= a; b++)
        Mx[s->HtHIdx[pos++]] += s->obsPrecision[i] * s->Hx[a] * s->Hx[b];
}

double logdet_Qiid_blocks(SamplerState *s)
{
  int b, k;
  double out;

  out = 0.0;
  for(b = 0; b < s->nRE; b++){
    k = s->reCovDim[b];
    if(k == 0)
      out -= s->reBlockSize[b] * std::log(s->sigmaSqRE[b]);
    else
      out -= (s->reBlockSize[b] / k) * s->reCovLogDet[b];
  }
  return out;
}

void apply_H_multiple(SamplerState *s, const double *vLat, int ncol, double *outObs)
{
  int i, j, a;

  if(!s->collapseIid){
    apply_A_multiple(s, vLat, ncol, outObs);
    return;
  }
  for(j = 0; j < ncol; j++){
    for(i = 0; i < s->n; i++){
      outObs[i + s->n * j] = 0.0;
      for(a = s->Hp[i]; a < s->Hp[i + 1]; a++)
        outObs[i + s->n * j] += s->Hx[a] * vLat[s->Hj[a] + s->qLatSystem * j];
    }
  }
}

void apply_H(SamplerState *s, const double *vLat, double *outObs)
{
  if(!s->collapseIid)
    apply_A(s, vLat, outObs);
  else
    apply_H_multiple(s, vLat, 1, outObs);
}

void apply_H_transpose_weighted_multiple(SamplerState *s, const double *vObs, int ncol, double *outLat)
{
  int i, j, a;
  double v;

  if(!s->collapseIid){
    apply_A_transpose_weighted_multiple(s, vObs, ncol, outLat);
    return;
  }
  for(j = 0; j < ncol; j++){
    for(i = 0; i < s->qLatSystem; i++)
      outLat[i + s->qLatSystem * j] = 0.0;
    for(i = 0; i < s->n; i++){
      v = s->obsPrecision[i] * vObs[i + s->n * j];
      for(a = s->Hp[i]; a < s->Hp[i + 1]; a++)
        outLat[s->Hj[a] + s->qLatSystem * j] += s->Hx[a] * v;
    }
  }
}

void apply_H_transpose_weighted(SamplerState *s, const double *vObs, double *outLat)
{
  if(!s->collapseIid)
    apply_A_transpose_weighted(s, vObs, outLat);
  else
    apply_H_transpose_weighted_multiple(s, vObs, 1, outLat);
}
