#ifndef NN_H
#define NN_H

#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

inline void getNNIndx(int i, int m, int &iNNIndx, int &iNN);
  
SEXP mkNNIndx(SEXP coords_r,
              SEXP n_r,
              SEXP m_r,
              SEXP r_r,
              SEXP nThreads_r);

SEXP mkNNIndxBrute(SEXP coords_r,
                   SEXP n_r,
                   SEXP m_r,
                   SEXP r_r,
                   SEXP nThreads_r);

#ifdef __cplusplus
}
#endif

#endif
