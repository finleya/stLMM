#include "stLMM_internal.h"

/* Sparse precision assembly for the collapsed likelihood.
 *
 * This file owns the CHOLMOD sparsity pattern for
 *   M = Q_w + A' W A
 * plus cache indices that let each MCMC iteration overwrite numeric values
 * without rebuilding the sparse structure. The current public model still
 * uses W = tau_sq^{-1} I; representing it as a vector is the foundation for
 * direct-estimate and other diagonal observation precision models.
 */

int cmp_int_asc(const void *a, const void *b)
{
  int ia, ib;

  ia = *((const int*)a);
  ib = *((const int*)b);

  if(ia < ib)
    return -1;
  if(ia > ib)
    return 1;
  return 0;
}

void pattern_col_append(PatternCol *cols, int col, int row)
{
  int newCap;
  int *tmp;

  if(row < col)
    Rf_error("pattern_col_append requires row >= col");

  if(cols[col].len >= cols[col].cap){
    newCap = (cols[col].cap == 0) ? 8 : 2 * cols[col].cap;
    tmp = R_Realloc(cols[col].rows, newCap, int);
    if(tmp == NULL)
      Rf_error("pattern_col_append: realloc failed");
    cols[col].rows = tmp;
    cols[col].cap = newCap;
  }

  cols[col].rows[cols[col].len] = row;
  cols[col].len++;
}

void pattern_append_lower(PatternCol *cols, int i, int j)
{
  if(i >= j)
    pattern_col_append(cols, j, i);
  else
    pattern_col_append(cols, i, j);
}

/*==========================================================================*/
/* sparse symbolic pattern construction for M_lat                           */
/*==========================================================================*/

void append_nngp_block_pattern_sparse(PatternCol *cols,
                                             TermState *term,
                                             GraphState *graph)
{
  int row, start, m, supportSize, a, b;
  int *nodes;
  int nodeA, nodeB;
  int off;

  off = term->wOffset;
  nodes = (int*)R_alloc(1 + (graph->totalNnbr > 0 ? graph->totalNnbr : 1),
                        sizeof(int));

  for(row = 0; row < graph->nNode; row++){
    start = graph->nnStart[row];
    m = graph->nnCount[row];
    supportSize = m + 1;

    nodes[0] = row;
    for(a = 0; a < m; a++)
      nodes[a + 1] = graph->nnIndx[start + a];

    for(a = 0; a < supportSize; a++){
      nodeA = off + nodes[a];
      for(b = 0; b <= a; b++){
        nodeB = off + nodes[b];
        pattern_append_lower(cols, nodeA, nodeB);
      }
    }
  }
}

void append_dagar_block_pattern_sparse(PatternCol *cols,
                                       TermState *term,
                                       GraphState *graph)
{
  int row, start, m, supportSize, a, b;
  int *nodes;
  int nodeA, nodeB;
  int off;

  off = term->wOffset;
  nodes = (int*)R_alloc(1 + (graph->totalParent > 0 ? graph->totalParent : 1),
                        sizeof(int));

  for(row = 0; row < graph->nNode; row++){
    start = graph->parentStart[row];
    m = graph->parentCount[row];
    supportSize = m + 1;

    nodes[0] = row;
    for(a = 0; a < m; a++)
      nodes[a + 1] = graph->parentIndx[start + a];

    for(a = 0; a < supportSize; a++){
      nodeA = off + nodes[a];
      for(b = 0; b <= a; b++){
        nodeB = off + nodes[b];
        pattern_append_lower(cols, nodeA, nodeB);
      }
    }
  }
}

void append_gp_block_pattern_sparse(PatternCol *cols,
                                           TermState *term,
                                           GraphState *graph)
{
  int i, j, off;
  (void)graph;

  off = term->wOffset;
  for(j = 0; j < term->nNode; j++){
    for(i = j; i < term->nNode; i++)
      pattern_col_append(cols, off + j, off + i);
  }
}

void append_ar1_block_pattern_sparse(PatternCol *cols,
                                            TermState *term,
                                            GraphState *graph)
{
  int node, off;
  (void)graph;

  off = term->wOffset;
  for(node = 0; node < term->nNode; node++)
    pattern_col_append(cols, off + node, off + node);

  for(node = 0; node < term->nNode - 1; node++)
    pattern_col_append(cols, off + node, off + node + 1);
}

void append_car_block_pattern_sparse(PatternCol *cols,
                                      TermState *term,
                                      GraphState *graph)
{
  int node, edge, off;

  off = term->wOffset;
  for(node = 0; node < term->nNode; node++)
    pattern_col_append(cols, off + node, off + node);

  for(edge = 0; edge < graph->nEdge; edge++)
    pattern_append_lower(cols, off + graph->edgeI[edge], off + graph->edgeJ[edge]);
}

static inline int car_time_node(int space, int time, int nTime)
{
  return space * nTime + time;
}

void append_car_time_block_pattern_sparse(PatternCol *cols,
                                          TermState *term,
                                          GraphState *graph)
{
  int s, t, edge, s1, s2, off, nTime;

  off = term->wOffset;
  nTime = graph->nTime;

  for(s = 0; s < graph->nSpace; s++){
    for(t = 0; t < nTime; t++)
      pattern_append_lower(cols, off + car_time_node(s, t, nTime), off + car_time_node(s, t, nTime));
    for(t = 0; t < nTime - 1; t++)
      pattern_append_lower(cols, off + car_time_node(s, t + 1, nTime), off + car_time_node(s, t, nTime));
  }

  for(edge = 0; edge < graph->nEdge; edge++){
    s1 = graph->edgeI[edge];
    s2 = graph->edgeJ[edge];
    for(t = 0; t < nTime; t++)
      pattern_append_lower(cols, off + car_time_node(s2, t, nTime), off + car_time_node(s1, t, nTime));
    for(t = 0; t < nTime - 1; t++){
      pattern_append_lower(cols, off + car_time_node(s2, t + 1, nTime), off + car_time_node(s1, t, nTime));
      pattern_append_lower(cols, off + car_time_node(s2, t, nTime), off + car_time_node(s1, t + 1, nTime));
    }
  }
}

static inline void append_space_time_pair_pattern(PatternCol *cols,
                                                  int off,
                                                  int nTime,
                                                  int s1,
                                                  int s2)
{
  int t;

  for(t = 0; t < nTime; t++)
    pattern_append_lower(cols, off + car_time_node(s2, t, nTime), off + car_time_node(s1, t, nTime));

  for(t = 0; t < nTime - 1; t++){
    if(s1 == s2){
      pattern_append_lower(cols, off + car_time_node(s1, t + 1, nTime), off + car_time_node(s1, t, nTime));
    } else {
      pattern_append_lower(cols, off + car_time_node(s2, t + 1, nTime), off + car_time_node(s1, t, nTime));
      pattern_append_lower(cols, off + car_time_node(s2, t, nTime), off + car_time_node(s1, t + 1, nTime));
    }
  }
}

void append_dagar_time_block_pattern_sparse(PatternCol *cols,
                                            TermState *term,
                                            GraphState *graph)
{
  int row, start, m, supportSize, a, b;
  int *nodes;
  int off, nTime;

  off = term->wOffset;
  nTime = graph->nTime;
  nodes = (int*)R_alloc(1 + (graph->totalParent > 0 ? graph->totalParent : 1),
                        sizeof(int));

  for(row = 0; row < graph->nSpace; row++){
    start = graph->parentStart[row];
    m = graph->parentCount[row];
    supportSize = m + 1;

    nodes[0] = row;
    for(a = 0; a < m; a++)
      nodes[a + 1] = graph->parentIndx[start + a];

    for(a = 0; a < supportSize; a++){
      for(b = 0; b <= a; b++)
        append_space_time_pair_pattern(cols, off, nTime, nodes[b], nodes[a]);
    }
  }
}

void append_AtA_pattern_sparse(PatternCol *cols, SamplerState *s)
{
  int i, t1, t2;
  int idx1, idx2;
  TermState *term1, *term2;

  for(i = 0; i < s->n; i++){
    for(t1 = 0; t1 < s->nTerms; t1++){
      term1 = s->terms + t1;
      idx1 = term1->wOffset + term1->map[i];

      for(t2 = 0; t2 <= t1; t2++){
        term2 = s->terms + t2;
        idx2 = term2->wOffset + term2->map[i];
        pattern_append_lower(cols, idx1, idx2);
      }
    }
  }
}

void build_M_lat_pattern_sparse(SamplerState *s)
{
  PatternCol *cols;
  int q, t, col, k, nnz;
  int *Mp, *Mi;
  double *Mx;
  TermState *term;
  GraphState *g;

  q = s->qLatTotal;

  cols = R_Calloc(q > 0 ? q : 1, PatternCol);
  for(col = 0; col < q; col++){
    cols[col].rows = NULL;
    cols[col].len = 0;
    cols[col].cap = 0;
  }

  for(t = 0; t < s->nTerms; t++){
    term = s->terms + t;
    g = s->graphs + term->graphIndex;

    if(g->type == GRAPH_NNGP)
      append_nngp_block_pattern_sparse(cols, term, g);
    else if(g->type == GRAPH_GP)
      append_gp_block_pattern_sparse(cols, term, g);
    else if(g->type == GRAPH_AR1)
      append_ar1_block_pattern_sparse(cols, term, g);
    else if(g->type == GRAPH_CAR)
      append_car_block_pattern_sparse(cols, term, g);
    else if(g->type == GRAPH_CAR_TIME)
      append_car_time_block_pattern_sparse(cols, term, g);
    else if(g->type == GRAPH_DAGAR)
      append_dagar_block_pattern_sparse(cols, term, g);
    else if(g->type == GRAPH_DAGAR_TIME)
      append_dagar_time_block_pattern_sparse(cols, term, g);
    else
      Rf_error("unsupported graph type in sparse latent pattern build");
  }

  append_AtA_pattern_sparse(cols, s);

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

  s->M_lat = M_cholmod_allocate_sparse(q, q, nnz, 1, 1, -1, CHOLMOD_REAL, &s->cm);
  if(s->M_lat == NULL || s->cm.status != CHOLMOD_OK)
    Rf_error("failed to allocate latent M pattern");

  Mp = (int*)s->M_lat->p;
  Mi = (int*)s->M_lat->i;
  Mx = (double*)s->M_lat->x;

  s->latentLowerMap = NULL;

  nnz = 0;
  Mp[0] = 0;
  for(col = 0; col < q; col++){
    for(k = 0; k < cols[col].len; k++){
      Mi[nnz] = cols[col].rows[k];
      Mx[nnz] = 0.0;
      nnz++;
    }
    Mp[col + 1] = nnz;
  }

  for(col = 0; col < q; col++){
    if(cols[col].rows != NULL)
      R_Free(cols[col].rows);
  }
  R_Free(cols);
}

/*==========================================================================*/
/* backend parsing                                                          */
/*==========================================================================*/


void zero_M_lat_numeric(SamplerState *s)
{
  double *Mx;
  int k, nz;

  Mx = (double*)s->M_lat->x;
  nz = (int)s->M_lat->nzmax;
  for(k = 0; k < nz; k++)
    Mx[k] = 0.0;
}

/*==========================================================================*/
/* one-time lookup helper used only while building caches                   */
/*==========================================================================*/

int find_M_lat_entry_index(SamplerState *s, int row, int col)
{
  int left, right, mid;
  int *Mp, *Mi;

  if(row < col)
    Rf_error("find_M_lat_entry_index called with row < col");

  Mp = (int*)s->M_lat->p;
  Mi = (int*)s->M_lat->i;

  left = Mp[col];
  right = Mp[col + 1] - 1;

  while(left <= right){
    mid = left + (right - left) / 2;

    if(Mi[mid] == row)
      return mid;
    else if(Mi[mid] < row)
      left = mid + 1;
    else
      right = mid - 1;
  }

  Rf_error("global M entry (%d,%d) missing from symbolic pattern", row + 1, col + 1);
  return -1;
}

/*==========================================================================*/
/* build caches from term-local contributions into global M@x               */
/*==========================================================================*/

void build_nngp_term_cache(SamplerState *s, TermState *term, GraphState *graph)
{
  int row, start, m, supportSize, a, b, off;
  int *nodes;
  int cachePos, totalPairs;
  int nodeA, nodeB;

  off = term->wOffset;

  term->nngpCacheStart = (int*)R_alloc(graph->nNode > 0 ? graph->nNode : 1, sizeof(int));

  totalPairs = 0;
  for(row = 0; row < graph->nNode; row++){
    m = graph->nnCount[row];
    supportSize = m + 1;
    totalPairs += supportSize * (supportSize + 1) / 2;
  }

  term->nngpCacheN = totalPairs;
  term->nngpCacheIdx = (int*)R_alloc(totalPairs > 0 ? totalPairs : 1, sizeof(int));

  nodes = (int*)R_alloc(1 + (graph->totalNnbr > 0 ? graph->totalNnbr : 1), sizeof(int));

  cachePos = 0;
  for(row = 0; row < graph->nNode; row++){
    start = graph->nnStart[row];
    m = graph->nnCount[row];
    supportSize = m + 1;

    term->nngpCacheStart[row] = cachePos;

    nodes[0] = row;
    for(a = 0; a < m; a++)
      nodes[a + 1] = graph->nnIndx[start + a];

    for(a = 0; a < supportSize; a++){
      nodeA = off + nodes[a];
      for(b = 0; b <= a; b++){
        nodeB = off + nodes[b];

        if(nodeA >= nodeB)
          term->nngpCacheIdx[cachePos] = find_M_lat_entry_index(s, nodeA, nodeB);
        else
          term->nngpCacheIdx[cachePos] = find_M_lat_entry_index(s, nodeB, nodeA);

        cachePos++;
      }
    }
  }
}

void build_dagar_term_cache(SamplerState *s, TermState *term, GraphState *graph)
{
  int row, start, m, supportSize, a, b, off;
  int *nodes;
  int cachePos, totalPairs;
  int nodeA, nodeB;

  off = term->wOffset;

  term->dagarCacheStart = (int*)R_alloc(graph->nNode > 0 ? graph->nNode : 1, sizeof(int));

  totalPairs = 0;
  for(row = 0; row < graph->nNode; row++){
    m = graph->parentCount[row];
    supportSize = m + 1;
    totalPairs += supportSize * (supportSize + 1) / 2;
  }

  term->dagarCacheN = totalPairs;
  term->dagarCacheIdx = (int*)R_alloc(totalPairs > 0 ? totalPairs : 1, sizeof(int));

  nodes = (int*)R_alloc(1 + (graph->totalParent > 0 ? graph->totalParent : 1), sizeof(int));

  cachePos = 0;
  for(row = 0; row < graph->nNode; row++){
    start = graph->parentStart[row];
    m = graph->parentCount[row];
    supportSize = m + 1;

    term->dagarCacheStart[row] = cachePos;

    nodes[0] = row;
    for(a = 0; a < m; a++)
      nodes[a + 1] = graph->parentIndx[start + a];

    for(a = 0; a < supportSize; a++){
      nodeA = off + nodes[a];
      for(b = 0; b <= a; b++){
        nodeB = off + nodes[b];

        if(nodeA >= nodeB)
          term->dagarCacheIdx[cachePos] = find_M_lat_entry_index(s, nodeA, nodeB);
        else
          term->dagarCacheIdx[cachePos] = find_M_lat_entry_index(s, nodeB, nodeA);

        cachePos++;
      }
    }
  }
}

void build_gp_term_cache(SamplerState *s, TermState *term, GraphState *graph)
{
  int i, j, off, nNode, pos;
  (void)graph;

  off = term->wOffset;
  nNode = term->nNode;

  term->gpCacheN = nNode * (nNode + 1) / 2;
  term->gpCacheIdx = (int*)R_alloc(term->gpCacheN > 0 ? term->gpCacheN : 1, sizeof(int));

  pos = 0;
  for(j = 0; j < nNode; j++){
    for(i = j; i < nNode; i++){
      term->gpCacheIdx[pos] = find_M_lat_entry_index(s, off + i, off + j);
      pos++;
    }
  }
}

void build_ar1_term_cache(SamplerState *s, TermState *term, GraphState *graph)
{
  int node, off;
  (void)graph;

  off = term->wOffset;

  term->ar1DiagCacheIdx = (int*)R_alloc(term->nNode > 0 ? term->nNode : 1, sizeof(int));
  term->ar1OffCacheIdx = (int*)R_alloc(term->nNode > 1 ? (term->nNode - 1) : 1, sizeof(int));

  for(node = 0; node < term->nNode; node++)
    term->ar1DiagCacheIdx[node] = find_M_lat_entry_index(s, off + node, off + node);

  for(node = 0; node < term->nNode - 1; node++)
    term->ar1OffCacheIdx[node] = find_M_lat_entry_index(s, off + node + 1, off + node);
}

void build_car_term_cache(SamplerState *s, TermState *term, GraphState *graph)
{
  int node, edge, off, i, j;

  off = term->wOffset;

  term->carDiagCacheIdx = (int*)R_alloc(term->nNode > 0 ? term->nNode : 1, sizeof(int));
  term->carOffCacheIdx = (int*)R_alloc(graph->nEdge > 0 ? graph->nEdge : 1, sizeof(int));

  for(node = 0; node < term->nNode; node++)
    term->carDiagCacheIdx[node] = find_M_lat_entry_index(s, off + node, off + node);

  for(edge = 0; edge < graph->nEdge; edge++){
    i = off + graph->edgeI[edge];
    j = off + graph->edgeJ[edge];
    if(i >= j)
      term->carOffCacheIdx[edge] = find_M_lat_entry_index(s, i, j);
    else
      term->carOffCacheIdx[edge] = find_M_lat_entry_index(s, j, i);
  }
}

static inline void cache_car_time_entry(SamplerState *s, TermState *term, int *pos, int a, int b)
{
  if(a >= b)
    term->carTimeCacheIdx[*pos] = find_M_lat_entry_index(s, a, b);
  else
    term->carTimeCacheIdx[*pos] = find_M_lat_entry_index(s, b, a);
  (*pos)++;
}

void build_car_time_term_cache(SamplerState *s, TermState *term, GraphState *graph)
{
  int sidx, tidx, edge, s1, s2, off, nTime, pos;

  off = term->wOffset;
  nTime = graph->nTime;
  term->carTimeCacheN =
    graph->nSpace * (2 * nTime - 1) +
    graph->nEdge * (3 * nTime - 2);
  term->carTimeCacheIdx = (int*)R_alloc(term->carTimeCacheN > 0 ? term->carTimeCacheN : 1, sizeof(int));

  pos = 0;
  for(sidx = 0; sidx < graph->nSpace; sidx++){
    for(tidx = 0; tidx < nTime; tidx++)
      cache_car_time_entry(s, term, &pos, off + car_time_node(sidx, tidx, nTime), off + car_time_node(sidx, tidx, nTime));
    for(tidx = 0; tidx < nTime - 1; tidx++)
      cache_car_time_entry(s, term, &pos, off + car_time_node(sidx, tidx + 1, nTime), off + car_time_node(sidx, tidx, nTime));
  }

  for(edge = 0; edge < graph->nEdge; edge++){
    s1 = graph->edgeI[edge];
    s2 = graph->edgeJ[edge];
    for(tidx = 0; tidx < nTime; tidx++)
      cache_car_time_entry(s, term, &pos, off + car_time_node(s2, tidx, nTime), off + car_time_node(s1, tidx, nTime));
    for(tidx = 0; tidx < nTime - 1; tidx++){
      cache_car_time_entry(s, term, &pos, off + car_time_node(s2, tidx + 1, nTime), off + car_time_node(s1, tidx, nTime));
      cache_car_time_entry(s, term, &pos, off + car_time_node(s2, tidx, nTime), off + car_time_node(s1, tidx + 1, nTime));
    }
  }

  if(pos != term->carTimeCacheN)
    Rf_error("internal error: car_time cache length mismatch");
}

void build_dagar_time_term_cache(SamplerState *s, TermState *term, GraphState *graph)
{
  int row, start, m, supportSize, a, b, off, nTime, pos, totalEntries;
  int *nodes;

  off = term->wOffset;
  nTime = graph->nTime;

  totalEntries = 0;
  for(row = 0; row < graph->nSpace; row++){
    m = graph->parentCount[row];
    supportSize = m + 1;
    for(a = 0; a < supportSize; a++){
      for(b = 0; b <= a; b++)
        totalEntries += (a == b) ? (2 * nTime - 1) : (3 * nTime - 2);
    }
  }
  term->carTimeCacheN = totalEntries;
  term->carTimeCacheIdx = (int*)R_alloc(term->carTimeCacheN > 0 ? term->carTimeCacheN : 1, sizeof(int));

  nodes = (int*)R_alloc(1 + (graph->totalParent > 0 ? graph->totalParent : 1), sizeof(int));

  pos = 0;
  for(row = 0; row < graph->nSpace; row++){
    start = graph->parentStart[row];
    m = graph->parentCount[row];
    supportSize = m + 1;

    nodes[0] = row;
    for(a = 0; a < m; a++)
      nodes[a + 1] = graph->parentIndx[start + a];

    for(a = 0; a < supportSize; a++){
      for(b = 0; b <= a; b++){
        int s1 = nodes[b];
        int s2 = nodes[a];
        int tidx;

        for(tidx = 0; tidx < nTime; tidx++)
          cache_car_time_entry(s, term, &pos, off + car_time_node(s2, tidx, nTime), off + car_time_node(s1, tidx, nTime));

        for(tidx = 0; tidx < nTime - 1; tidx++){
          if(s1 == s2){
            cache_car_time_entry(s, term, &pos, off + car_time_node(s1, tidx + 1, nTime), off + car_time_node(s1, tidx, nTime));
          } else {
            cache_car_time_entry(s, term, &pos, off + car_time_node(s2, tidx + 1, nTime), off + car_time_node(s1, tidx, nTime));
            cache_car_time_entry(s, term, &pos, off + car_time_node(s2, tidx, nTime), off + car_time_node(s1, tidx + 1, nTime));
          }
        }
      }
    }
  }

  if(pos != term->carTimeCacheN)
    Rf_error("internal error: dagar_time cache length mismatch");
}

void build_AtA_cache(SamplerState *s)
{
  int i, t1, t2, pair, idx1, idx2;
  TermState *term1, *term2;

  s->nTermPairs = s->nTerms * (s->nTerms + 1) / 2;
  s->termPair1 = (int*)R_alloc(s->nTermPairs > 0 ? s->nTermPairs : 1, sizeof(int));
  s->termPair2 = (int*)R_alloc(s->nTermPairs > 0 ? s->nTermPairs : 1, sizeof(int));
  s->AtAIdx = (int*)R_alloc((s->n > 0 ? s->n : 1) * (s->nTermPairs > 0 ? s->nTermPairs : 1), sizeof(int));

  pair = 0;
  for(t1 = 0; t1 < s->nTerms; t1++){
    for(t2 = 0; t2 <= t1; t2++){
      s->termPair1[pair] = t1;
      s->termPair2[pair] = t2;
      pair++;
    }
  }

  for(i = 0; i < s->n; i++){
    for(pair = 0; pair < s->nTermPairs; pair++){
      t1 = s->termPair1[pair];
      t2 = s->termPair2[pair];

      term1 = s->terms + t1;
      term2 = s->terms + t2;

      idx1 = term1->wOffset + term1->map[i];
      idx2 = term2->wOffset + term2->map[i];

      if(idx1 >= idx2)
        s->AtAIdx[i + s->n * pair] = find_M_lat_entry_index(s, idx1, idx2);
      else
        s->AtAIdx[i + s->n * pair] = find_M_lat_entry_index(s, idx2, idx1);
    }
  }
}

void build_M_lat_index_caches(SamplerState *s)
{
  int t;
  GraphState *g;
  TermState *term;

  for(t = 0; t < s->nTerms; t++){
    term = s->terms + t;
    g = s->graphs + term->graphIndex;

    if(g->type == GRAPH_NNGP)
      build_nngp_term_cache(s, term, g);
    else if(g->type == GRAPH_GP)
      build_gp_term_cache(s, term, g);
    else if(g->type == GRAPH_AR1)
      build_ar1_term_cache(s, term, g);
    else if(g->type == GRAPH_CAR)
      build_car_term_cache(s, term, g);
    else if(g->type == GRAPH_CAR_TIME)
      build_car_time_term_cache(s, term, g);
    else if(g->type == GRAPH_DAGAR)
      build_dagar_term_cache(s, term, g);
    else if(g->type == GRAPH_DAGAR_TIME)
      build_dagar_time_term_cache(s, term, g);
    else
      Rf_error("unsupported graph type while building M caches");
  }

  build_AtA_cache(s);
}

/*==========================================================================*/
/* numeric assembly using cached M@x indices                                */
/*==========================================================================*/

void assemble_nngp_block_into_M(SamplerState *s, TermState *term, GraphState *graph)
{
  int row, a, b, start, m, supportSize;
  int *cols;
  double *w;
  double fj, wa, wb, qab, sigmaSqInv;
  int maxSupport;
  int cachePos;
  double *Mx;

  sigmaSqInv = 1.0 / term->sigmaSq;
  Mx = (double*)s->M_lat->x;

  maxSupport = 1;
  if(graph->totalNnbr > 0)
    maxSupport += graph->totalNnbr;

  cols = (int*)R_alloc(maxSupport, sizeof(int));
  w = (double*)R_alloc(maxSupport, sizeof(double));

  for(row = 0; row < graph->nNode; row++){
    start = graph->nnStart[row];
    m = graph->nnCount[row];
    supportSize = m + 1;

    cols[0] = row;
    w[0] = 1.0;
    for(a = 0; a < m; a++){
      cols[a + 1] = graph->nnIndx[start + a];
      w[a + 1] = -term->B[start + a];
    }

    fj = term->F[row];
    cachePos = term->nngpCacheStart[row];

    for(a = 0; a < supportSize; a++){
      wa = w[a];
      for(b = 0; b <= a; b++){
        wb = w[b];
        qab = sigmaSqInv * fj * wa * wb;
        if(qab != 0.0)
          Mx[term->nngpCacheIdx[cachePos]] += qab;
        cachePos++;
      }
    }
  }
}

void assemble_dagar_block_into_M(SamplerState *s, TermState *term, GraphState *graph)
{
  int row, a, b, start, m, supportSize;
  int *cols;
  double *w;
  double rho, rho2, denom, bi, fi, wa, wb, qab, sigmaSqInv;
  int maxSupport, cachePos;
  double *Mx;

  rho = term->theta[0];
  if(rho <= 0.0 || rho >= 1.0)
    Rf_error("assemble_dagar_block_into_M requires rho in (0,1)");

  rho2 = rho * rho;
  sigmaSqInv = 1.0 / term->sigmaSq;
  Mx = (double*)s->M_lat->x;

  maxSupport = 1;
  if(graph->totalParent > 0)
    maxSupport += graph->totalParent;

  cols = (int*)R_alloc(maxSupport, sizeof(int));
  w = (double*)R_alloc(maxSupport, sizeof(double));

  for(row = 0; row < graph->nNode; row++){
    start = graph->parentStart[row];
    m = graph->parentCount[row];
    supportSize = m + 1;

    denom = 1.0 + ((double)m - 1.0) * rho2;
    if(denom <= 0.0 || !std::isfinite(denom))
      Rf_error("DAGAR conditional precision denominator became invalid");
    bi = (m > 0) ? rho / denom : 0.0;
    fi = denom / (1.0 - rho2);

    cols[0] = row;
    w[0] = 1.0;
    for(a = 0; a < m; a++){
      cols[a + 1] = graph->parentIndx[start + a];
      w[a + 1] = -bi;
    }

    cachePos = term->dagarCacheStart[row];
    for(a = 0; a < supportSize; a++){
      wa = w[a];
      for(b = 0; b <= a; b++){
        wb = w[b];
        qab = sigmaSqInv * fi * wa * wb;
        if(qab != 0.0)
          Mx[term->dagarCacheIdx[cachePos]] += qab;
        cachePos++;
      }
    }
  }
}

void assemble_gp_block_into_M(SamplerState *s, TermState *term, GraphState *graph)
{
  int i, j, pos, nNode;
  double qij;
  double *Mx;
  (void)graph;

  Mx = (double*)s->M_lat->x;
  nNode = term->nNode;

  pos = 0;
  for(j = 0; j < nNode; j++){
    for(i = j; i < nNode; i++){
      qij = term->Q[i + nNode * j];
      if(qij != 0.0)
        Mx[term->gpCacheIdx[pos]] += qij;
      pos++;
    }
  }
}

void assemble_ar1_block_into_M(SamplerState *s, TermState *term, GraphState *graph)
{
  int node;
  double phi, den, cdiagEnd, cdiagMid, coff, sigmaSqInv, diagVal;
  double *Mx;
  (void)graph;

  phi = term->theta[0];
  if(std::fabs(phi) >= 1.0)
    Rf_error("assemble_ar1_block_into_M requires |phi| < 1");

  Mx = (double*)s->M_lat->x;
  sigmaSqInv = 1.0 / term->sigmaSq;

  if(term->nNode == 1){
    Mx[term->ar1DiagCacheIdx[0]] += sigmaSqInv;
    return;
  }

  den = 1.0 - phi * phi;
  cdiagEnd = 1.0 / den;
  cdiagMid = (1.0 + phi * phi) / den;
  coff = -phi / den;

  for(node = 0; node < term->nNode; node++){
    diagVal = (node == 0 || node == term->nNode - 1) ? cdiagEnd : cdiagMid;
    Mx[term->ar1DiagCacheIdx[node]] += sigmaSqInv * diagVal;
  }

  for(node = 0; node < term->nNode - 1; node++)
    Mx[term->ar1OffCacheIdx[node]] += sigmaSqInv * coff;
}

void assemble_car_block_into_M(SamplerState *s, TermState *term, GraphState *graph)
{
  int node, edge;
  double rho, sigmaSqInv, diagVal;
  double *Mx;

  rho = term->theta[0];
  if(rho <= 0.0 || rho >= 1.0)
    Rf_error("assemble_car_block_into_M requires rho in (0,1)");

  Mx = (double*)s->M_lat->x;
  sigmaSqInv = 1.0 / term->sigmaSq;

  for(node = 0; node < term->nNode; node++){
    if(graph->carModel == CAR_MODEL_LEROUX)
      diagVal = (1.0 - rho) + rho * graph->degree[node];
    else
      diagVal = graph->degree[node];
    Mx[term->carDiagCacheIdx[node]] += sigmaSqInv * diagVal;
  }

  for(edge = 0; edge < graph->nEdge; edge++)
    Mx[term->carOffCacheIdx[edge]] += -sigmaSqInv * rho * graph->edgeW[edge];
}

static inline double car_time_ar1_diag(int time, int nTime, double phi)
{
  double den;
  if(nTime == 1)
    return 1.0;
  den = 1.0 - phi * phi;
  if(time == 0 || time == nTime - 1)
    return 1.0 / den;
  return (1.0 + phi * phi) / den;
}

static inline double car_time_ar1_off(double phi)
{
  return -phi / (1.0 - phi * phi);
}

static inline double car_time_exp_gap_phi(GraphState *graph, int gap, double lambda)
{
  return std::exp(-lambda * graph->timeDelta[gap]);
}

static inline double car_time_exp_diag(GraphState *graph, int time, double lambda)
{
  double phi, den, prevPhi, prevDen;
  int nTime;

  nTime = graph->nTime;
  if(nTime == 1)
    return 1.0;
  if(time == 0){
    phi = car_time_exp_gap_phi(graph, 0, lambda);
    den = 1.0 - phi * phi;
    return 1.0 / den;
  }
  if(time == nTime - 1){
    prevPhi = car_time_exp_gap_phi(graph, time - 1, lambda);
    prevDen = 1.0 - prevPhi * prevPhi;
    return 1.0 / prevDen;
  }

  prevPhi = car_time_exp_gap_phi(graph, time - 1, lambda);
  prevDen = 1.0 - prevPhi * prevPhi;
  phi = car_time_exp_gap_phi(graph, time, lambda);
  den = 1.0 - phi * phi;
  return 1.0 / prevDen + phi * phi / den;
}

static inline double car_time_exp_off(GraphState *graph, int gap, double lambda)
{
  double phi;
  phi = car_time_exp_gap_phi(graph, gap, lambda);
  return -phi / (1.0 - phi * phi);
}

static inline double car_time_diag(GraphState *graph, int time, double thetaTime)
{
  if(graph->timeModel == TIME_MODEL_EXP)
    return car_time_exp_diag(graph, time, thetaTime);
  return car_time_ar1_diag(time, graph->nTime, thetaTime);
}

static inline double car_time_off(GraphState *graph, int gap, double thetaTime)
{
  if(graph->timeModel == TIME_MODEL_EXP)
    return car_time_exp_off(graph, gap, thetaTime);
  return car_time_ar1_off(thetaTime);
}

void assemble_car_time_block_into_M(SamplerState *s, TermState *term, GraphState *graph)
{
  int sidx, tidx, edge, pos;
  double rho, thetaTime, sigmaSqInv, tval, spval;
  double *Mx;

  rho = term->theta[0];
  thetaTime = term->theta[1];
  if(rho <= 0.0 || rho >= 1.0)
    Rf_error("assemble_car_time_block_into_M requires rho in (0,1)");
  if(graph->timeModel == TIME_MODEL_AR1){
    if(std::fabs(thetaTime) >= 1.0)
      Rf_error("assemble_car_time_block_into_M requires |phi| < 1");
  } else if(graph->timeModel == TIME_MODEL_EXP){
    if(thetaTime <= 0.0)
      Rf_error("assemble_car_time_block_into_M requires lambda > 0");
  }

  Mx = (double*)s->M_lat->x;
  sigmaSqInv = 1.0 / term->sigmaSq;

  pos = 0;
  for(sidx = 0; sidx < graph->nSpace; sidx++){
    if(graph->carModel == CAR_MODEL_LEROUX)
      spval = (1.0 - rho) + rho * graph->degree[sidx];
    else
      spval = graph->degree[sidx];
    for(tidx = 0; tidx < graph->nTime; tidx++){
      tval = car_time_diag(graph, tidx, thetaTime);
      Mx[term->carTimeCacheIdx[pos++]] += sigmaSqInv * spval * tval;
    }
    for(tidx = 0; tidx < graph->nTime - 1; tidx++){
      tval = car_time_off(graph, tidx, thetaTime);
      Mx[term->carTimeCacheIdx[pos++]] += sigmaSqInv * spval * tval;
    }
  }

  for(edge = 0; edge < graph->nEdge; edge++){
    spval = -rho * graph->edgeW[edge];
    for(tidx = 0; tidx < graph->nTime; tidx++){
      tval = car_time_diag(graph, tidx, thetaTime);
      Mx[term->carTimeCacheIdx[pos++]] += sigmaSqInv * spval * tval;
    }
    for(tidx = 0; tidx < graph->nTime - 1; tidx++){
      tval = car_time_off(graph, tidx, thetaTime);
      Mx[term->carTimeCacheIdx[pos++]] += sigmaSqInv * spval * tval;
      Mx[term->carTimeCacheIdx[pos++]] += sigmaSqInv * spval * tval;
    }
  }

  if(pos != term->carTimeCacheN)
    Rf_error("internal error: car_time assembly length mismatch");
}

void assemble_dagar_time_block_into_M(SamplerState *s, TermState *term, GraphState *graph)
{
  int row, a, b, start, m, supportSize, tidx, pos;
  int *nodes;
  double *w;
  double rho, rho2, denom, bi, fi, wa, wb, spval, tval, sigmaSqInv;
  double thetaTime;
  double *Mx;

  rho = term->theta[0];
  thetaTime = term->theta[1];
  if(rho <= 0.0 || rho >= 1.0)
    Rf_error("assemble_dagar_time_block_into_M requires rho in (0,1)");
  if(graph->timeModel == TIME_MODEL_AR1){
    if(std::fabs(thetaTime) >= 1.0)
      Rf_error("assemble_dagar_time_block_into_M requires |phi| < 1");
  } else if(graph->timeModel == TIME_MODEL_EXP){
    if(thetaTime <= 0.0)
      Rf_error("assemble_dagar_time_block_into_M requires lambda > 0");
  }

  rho2 = rho * rho;
  sigmaSqInv = 1.0 / term->sigmaSq;
  Mx = (double*)s->M_lat->x;

  nodes = (int*)R_alloc(1 + (graph->totalParent > 0 ? graph->totalParent : 1), sizeof(int));
  w = (double*)R_alloc(1 + (graph->totalParent > 0 ? graph->totalParent : 1), sizeof(double));

  pos = 0;
  for(row = 0; row < graph->nSpace; row++){
    start = graph->parentStart[row];
    m = graph->parentCount[row];
    supportSize = m + 1;

    denom = 1.0 + ((double)m - 1.0) * rho2;
    if(denom <= 0.0 || !std::isfinite(denom))
      Rf_error("DAGAR-time conditional precision denominator became invalid");
    bi = (m > 0) ? rho / denom : 0.0;
    fi = denom / (1.0 - rho2);

    nodes[0] = row;
    w[0] = 1.0;
    for(a = 0; a < m; a++){
      nodes[a + 1] = graph->parentIndx[start + a];
      w[a + 1] = -bi;
    }

    for(a = 0; a < supportSize; a++){
      wa = w[a];
      for(b = 0; b <= a; b++){
        wb = w[b];
        spval = sigmaSqInv * fi * wa * wb;

        for(tidx = 0; tidx < graph->nTime; tidx++){
          tval = car_time_diag(graph, tidx, thetaTime);
          Mx[term->carTimeCacheIdx[pos++]] += spval * tval;
        }

        for(tidx = 0; tidx < graph->nTime - 1; tidx++){
          tval = car_time_off(graph, tidx, thetaTime);
          if(a == b){
            Mx[term->carTimeCacheIdx[pos++]] += spval * tval;
          } else {
            Mx[term->carTimeCacheIdx[pos++]] += spval * tval;
            Mx[term->carTimeCacheIdx[pos++]] += spval * tval;
          }
        }
      }
    }
  }

  if(pos != term->carTimeCacheN)
    Rf_error("internal error: dagar_time assembly length mismatch");
}

void refresh_observation_precision(SamplerState *s)
{
  int i;
  double tauInv;

  if(is_pg_likelihood(s)){
    /*
      Polya-Gamma likelihoods set obsPrecision[i] = omega_i directly.
      In the collapsed Gaussian algebra this is the same slot occupied by
      tau_sq^{-1} (or a row-specific residual precision) for Gaussian models:

        M = Q_w + A' W A,  W = diag(obsPrecision).

      It is not a residual variance, and refreshing Gaussian residual
      precision here would overwrite the current PG augmentation.
    */
    return;
  }

  if(s->residualModel == 1){
    for(i = 0; i < s->n; i++)
      s->obsPrecision[i] = s->obsPrecisionFixed[i];
    return;
  }

  if(s->residualModel == 2){
    for(i = 0; i < s->n; i++)
      s->obsPrecision[i] = 1.0 / s->residualVariance[s->residualGroupIndex[i]];
    return;
  }

  if(s->residualModel == 3){
    double kappa, tau0Sq, w, variance;
    kappa = s->residualVariance[0];
    tau0Sq = (s->nResidualGroup > 1) ? s->residualVariance[1] : 1.0;
    for(i = 0; i < s->n; i++){
      w = s->residualScaledWeight[i];
      if(s->nResidualGroup > 1)
        variance = kappa * std::exp(w * std::log(s->residualScaledVhat[i]) +
                                    (1.0 - w) * std::log(tau0Sq));
      else
        variance = kappa * s->residualScaledVhat[i];
      if(!R_FINITE(variance) || !(variance > 0.0))
        Rf_error("scaled residual variance became non-finite or non-positive");
      s->obsPrecision[i] = 1.0 / variance;
    }
    return;
  }

  tauInv = 1.0 / s->tauSq;
  for(i = 0; i < s->n; i++)
    s->obsPrecision[i] = tauInv;
}

void add_observation_precision_AtA_into_M(SamplerState *s)
{
  int i, pair, t1, t2;
  double p_i, v1, v2;
  double *Mx;
  TermState *term1, *term2;

  Mx = (double*)s->M_lat->x;

  for(pair = 0; pair < s->nTermPairs; pair++){
    t1 = s->termPair1[pair];
    t2 = s->termPair2[pair];
    term1 = s->terms + t1;
    term2 = s->terms + t2;

    for(i = 0; i < s->n; i++){
      p_i = s->obsPrecision[i];
      v1 = term1->scale[i];
      v2 = term2->scale[i];
      Mx[s->AtAIdx[i + s->n * pair]] += p_i * v1 * v2;
    }
  }
}

void assemble_M_lat_numeric(SamplerState *s)
{
  int t;
  TermState *term;
  GraphState *g;

  refresh_observation_precision(s);
  zero_M_lat_numeric(s);

  for(t = 0; t < s->nTerms; t++){
    term = s->terms + t;
    g = s->graphs + term->graphIndex;

    if(g->type == GRAPH_NNGP)
      assemble_nngp_block_into_M(s, term, g);
    else if(g->type == GRAPH_GP)
      assemble_gp_block_into_M(s, term, g);
    else if(g->type == GRAPH_AR1)
      assemble_ar1_block_into_M(s, term, g);
    else if(g->type == GRAPH_CAR)
      assemble_car_block_into_M(s, term, g);
    else if(g->type == GRAPH_CAR_TIME)
      assemble_car_time_block_into_M(s, term, g);
    else if(g->type == GRAPH_DAGAR)
      assemble_dagar_block_into_M(s, term, g);
    else if(g->type == GRAPH_DAGAR_TIME)
      assemble_dagar_time_block_into_M(s, term, g);
    else
      Rf_error("unsupported graph type in M_lat numeric assembly");
  }

  add_observation_precision_AtA_into_M(s);
}

void apply_A_transpose(SamplerState *s, const double *vObs, double *outLat)
{
  int i, t, idx;
  TermState *term;

  for(i = 0; i < s->qLatTotal; i++)
    outLat[i] = 0.0;

  for(i = 0; i < s->n; i++){
    for(t = 0; t < s->nTerms; t++){
      term = s->terms + t;
      idx = term->wOffset + term->map[i];
      outLat[idx] += term->scale[i] * vObs[i];
    }
  }
}

void apply_A_transpose_weighted(SamplerState *s, const double *vObs, double *outLat)
{
  int i, t, idx;
  double wi;
  TermState *term;

  for(i = 0; i < s->qLatTotal; i++)
    outLat[i] = 0.0;

  for(i = 0; i < s->n; i++){
    wi = s->obsPrecision[i] * vObs[i];
    for(t = 0; t < s->nTerms; t++){
      term = s->terms + t;
      idx = term->wOffset + term->map[i];
      outLat[idx] += term->scale[i] * wi;
    }
  }
}

void apply_A(SamplerState *s, const double *vLat, double *outObs)
{
  int i, t, idx;
  TermState *term;

  for(i = 0; i < s->n; i++)
    outObs[i] = 0.0;

  for(i = 0; i < s->n; i++){
    for(t = 0; t < s->nTerms; t++){
      term = s->terms + t;
      idx = term->wOffset + term->map[i];
      outObs[i] += term->scale[i] * vLat[idx];
    }
  }
}

void solve_M_lat(SamplerState *s, const double *b, double *x)
{
  cholmod_dense rhs_view, *sol;
  double *solx;
  int i;

  /*
    Reuse sampler-owned dense storage for the RHS. We still rely on
    M_cholmod_solve() to allocate the returned solution object, but avoid
    allocating/copying into a fresh CHOLMOD RHS object on every call.
  */
  for(i = 0; i < s->qLatTotal; i++)
    s->latentSolveRhs[i] = b[i];

  std::memset(&rhs_view, 0, sizeof(cholmod_dense));
  rhs_view.nrow = s->qLatTotal;
  rhs_view.ncol = 1;
  rhs_view.nzmax = s->qLatTotal;
  rhs_view.d = s->qLatTotal;
  rhs_view.xtype = CHOLMOD_REAL;
  rhs_view.dtype = CHOLMOD_DOUBLE;
  rhs_view.x = (void*)s->latentSolveRhs;
  rhs_view.z = NULL;

  sol = M_cholmod_solve(CHOLMOD_A, s->M_lat_fac, &rhs_view, &s->cm);
  if(sol == NULL || s->cm.status != CHOLMOD_OK)
    Rf_error("solve_M_lat: cholmod_solve(A) failed");

  solx = (double*)sol->x;
  for(i = 0; i < s->qLatTotal; i++)
    x[i] = solx[i];

  M_cholmod_free_dense(&sol, &s->cm);
}

void apply_Vinv(SamplerState *s, const double *vObs, double *outObs)
{
  int i;
  double wi;
  double *Atv, *MinvAtv, *A_Minv_Atv;

  Atv = (double*)R_alloc(s->qLatTotal > 0 ? s->qLatTotal : 1, sizeof(double));
  MinvAtv = (double*)R_alloc(s->qLatTotal > 0 ? s->qLatTotal : 1, sizeof(double));
  A_Minv_Atv = (double*)R_alloc(s->n > 0 ? s->n : 1, sizeof(double));

  apply_A_transpose_weighted(s, vObs, Atv);
  solve_M_lat(s, Atv, MinvAtv);
  apply_A(s, MinvAtv, A_Minv_Atv);

  for(i = 0; i < s->n; i++){
    wi = s->obsPrecision[i];
    outObs[i] = wi * vObs[i] - wi * A_Minv_Atv[i];
  }
}

double compute_logDetV(SamplerState *s)
{
  int i;
  double logDetR;

  logDetR = 0.0;
  for(i = 0; i < s->n; i++)
    logDetR -= std::log(s->obsPrecision[i]);

  return logDetR + logDetFactor(s->M_lat_fac) - logdet_Qw_blocks(s);
}

double quadform_Vinv(SamplerState *s, const double *vObs)
{
  int i;
  double out;
  double *tmp;

  tmp = (double*)R_alloc(s->n > 0 ? s->n : 1, sizeof(double));
  apply_Vinv(s, vObs, tmp);

  out = 0.0;
  for(i = 0; i < s->n; i++)
    out += vObs[i] * tmp[i];

  return out;
}

void apply_A_transpose_multiple(SamplerState *s, const double *vObs, int ncol, double *outLat)
{
  int i, t, idx, j;

  for(j = 0; j < s->qLatTotal * ncol; j++)
    outLat[j] = 0.0;

  for(i = 0; i < s->n; i++){
    for(t = 0; t < s->nTerms; t++){
      TermState *term = s->terms + t;
      idx = term->wOffset + term->map[i];
      for(j = 0; j < ncol; j++)
        outLat[idx + s->qLatTotal * j] += term->scale[i] * vObs[i + s->n * j];
    }
  }
}

void apply_A_transpose_weighted_multiple(SamplerState *s, const double *vObs, int ncol, double *outLat)
{
  int i, t, idx, j;
  double wi;

  for(j = 0; j < s->qLatTotal * ncol; j++)
    outLat[j] = 0.0;

  for(i = 0; i < s->n; i++){
    wi = s->obsPrecision[i];
    for(t = 0; t < s->nTerms; t++){
      TermState *term = s->terms + t;
      idx = term->wOffset + term->map[i];
      for(j = 0; j < ncol; j++)
        outLat[idx + s->qLatTotal * j] += term->scale[i] * wi * vObs[i + s->n * j];
    }
  }
}

void apply_A_multiple(SamplerState *s, const double *vLat, int ncol, double *outObs)
{
  int i, t, idx, j;

  for(j = 0; j < s->n * ncol; j++)
    outObs[j] = 0.0;

  for(i = 0; i < s->n; i++){
    for(t = 0; t < s->nTerms; t++){
      TermState *term = s->terms + t;
      idx = term->wOffset + term->map[i];
      for(j = 0; j < ncol; j++)
        outObs[i + s->n * j] += term->scale[i] * vLat[idx + s->qLatTotal * j];
    }
  }
}

void solve_M_lat_multiple(SamplerState *s, const double *B, int ncol, double *X)
{
  cholmod_dense *rhs, *sol;
  double *rhsx, *solx;
  int i, j;

  rhs = M_cholmod_allocate_dense(s->qLatTotal, ncol, s->qLatTotal, CHOLMOD_REAL, &s->cm);
  if(rhs == NULL || s->cm.status != CHOLMOD_OK)
    Rf_error("solve_M_lat_multiple: failed to allocate rhs");

  rhsx = (double*)rhs->x;
  for(j = 0; j < ncol; j++)
    for(i = 0; i < s->qLatTotal; i++)
      rhsx[i + s->qLatTotal * j] = B[i + s->qLatTotal * j];

  sol = M_cholmod_solve(CHOLMOD_A, s->M_lat_fac, rhs, &s->cm);
  if(sol == NULL || s->cm.status != CHOLMOD_OK){
    M_cholmod_free_dense(&rhs, &s->cm);
    Rf_error("solve_M_lat_multiple: cholmod_solve(A) failed");
  }

  solx = (double*)sol->x;
  for(j = 0; j < ncol; j++)
    for(i = 0; i < s->qLatTotal; i++)
      X[i + s->qLatTotal * j] = solx[i + s->qLatTotal * j];

  M_cholmod_free_dense(&sol, &s->cm);
  M_cholmod_free_dense(&rhs, &s->cm);
}

void apply_Vinv_multiple(SamplerState *s, const double *Vobs, int ncol, double *outObs)
{
  int i, j;
  double wi;
  double *AtV, *MinvAtV, *A_Minv_AtV;

  AtV = (double*)R_alloc((s->qLatTotal > 0 ? s->qLatTotal : 1) * (ncol > 0 ? ncol : 1), sizeof(double));
  MinvAtV = (double*)R_alloc((s->qLatTotal > 0 ? s->qLatTotal : 1) * (ncol > 0 ? ncol : 1), sizeof(double));
  A_Minv_AtV = (double*)R_alloc((s->n > 0 ? s->n : 1) * (ncol > 0 ? ncol : 1), sizeof(double));

  apply_A_transpose_weighted_multiple(s, Vobs, ncol, AtV);
  solve_M_lat_multiple(s, AtV, ncol, MinvAtV);
  apply_A_multiple(s, MinvAtV, ncol, A_Minv_AtV);

  for(j = 0; j < ncol; j++)
    for(i = 0; i < s->n; i++){
      wi = s->obsPrecision[i];
      outObs[i + s->n * j] = wi * Vobs[i + s->n * j] - wi * A_Minv_AtV[i + s->n * j];
    }
}

void form_XtVinvX_and_rhs(SamplerState *s, const double *ytilde, double *XtVinvX, double *XtVinvy)
{
  int i, j, inc;
  double one, zero;
  double *VinvX, *Vinvy;

  if(s->p <= 0)
    return;

  one = 1.0;
  zero = 0.0;
  inc = 1;

  VinvX = (double*)R_alloc(s->n * s->p, sizeof(double));
  Vinvy = (double*)R_alloc(s->n, sizeof(double));

  apply_Vinv_multiple(s, s->X, s->p, VinvX);
  apply_Vinv(s, ytilde, Vinvy);

  F77_CALL(dgemm)("T", "N", &s->p, &s->p, &s->n,
                  &one, s->X, &s->n, VinvX, &s->n,
                  &zero, XtVinvX, &s->p FCONE FCONE);

  F77_CALL(dgemv)("T", &s->n, &s->p,
                  &one, s->X, &s->n, Vinvy, &inc,
                  &zero, XtVinvy, &inc FCONE);
}

void form_ZtVinvZ_and_rhs(SamplerState *s, const double *ytilde, double *ZtVinvZ, double *ZtVinvy)
{
  int i, j;
  if(s->q <= 0)
    return;

  for(j = 0; j < s->q; j++){
    scatter_Z_col_to_dense(s, j, s->nWork1);
    apply_Vinv(s, s->nWork1, s->nWork2);

    for(i = j; i < s->q; i++)
      ZtVinvZ[i + s->q * j] = sparse_Z_col_dot_dense(s, i, s->nWork2);

  }

  apply_Vinv(s, ytilde, s->nWork1);
  sparse_Zt_mult(s, s->nWork1, ZtVinvy);
}
