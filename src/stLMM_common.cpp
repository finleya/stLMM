#include "stLMM_internal.h"

SEXP getListElement(SEXP list, const char *name)
{
  SEXP names;
  int i, n;

  names = Rf_getAttrib(list, R_NamesSymbol);
  if(names == R_NilValue)
    return R_NilValue;

  n = LENGTH(list);
  for(i = 0; i < n; i++){
    if(std::strcmp(CHAR(STRING_ELT(names, i)), name) == 0)
      return VECTOR_ELT(list, i);
  }

  return R_NilValue;
}

int getListIndex(SEXP list, const char *name)
{
  SEXP names;
  int i, n;

  names = Rf_getAttrib(list, R_NamesSymbol);
  if(names == R_NilValue)
    return -1;

  n = LENGTH(list);
  for(i = 0; i < n; i++){
    if(std::strcmp(CHAR(STRING_ELT(names, i)), name) == 0)
      return i;
  }

  return -1;
}

void setListElementByName(SEXP list, const char *name, SEXP value)
{
  int idx;

  PROTECT(value);
  idx = getListIndex(list, name);
  if(idx < 0){
    UNPROTECT(1);
    Rf_error("list element '%s' not found", name);
  }
  SET_VECTOR_ELT(list, idx, value);
  UNPROTECT(1);
}

int as_flag_scalar(SEXP x, const char *where)
{
  if(x == R_NilValue)
    Rf_error("%s is required", where);

  if(LENGTH(x) != 1)
    Rf_error("%s must be a scalar", where);

  if(TYPEOF(x) == LGLSXP)
    return LOGICAL(x)[0] != 0;
  if(TYPEOF(x) == INTSXP)
    return INTEGER(x)[0] != 0;
  if(TYPEOF(x) == REALSXP)
    return REAL(x)[0] != 0.0;

  Rf_error("%s must be logical, integer, or numeric", where);
  return 0;
}

int as_nonneg_int_scalar(SEXP x, const char *where)
{
  int out;

  if(x == R_NilValue)
    Rf_error("%s is required", where);

  if(LENGTH(x) != 1)
    Rf_error("%s must be a scalar", where);

  if(TYPEOF(x) == INTSXP)
    out = INTEGER(x)[0];
  else if(TYPEOF(x) == REALSXP)
    out = (int) REAL(x)[0];
  else if(TYPEOF(x) == LGLSXP)
    out = LOGICAL(x)[0] ? 1 : 0;
  else
    Rf_error("%s must be integer, numeric, or logical", where);

  if(out < 0)
    Rf_error("%s must be nonnegative", where);

  return out;
}

const char *ordering_label_from_common(cholmod_common *cm)
{
  int m, ordtype;

  m = cm->selected;
  ordtype = cm->method[m].ordering;

  switch(ordtype){
    case CHOLMOD_NATURAL: return "natural";
    case CHOLMOD_AMD:     return "amd";
    case CHOLMOD_METIS:   return "metis";
    case CHOLMOD_NESDIS:  return "nesdis";
    case CHOLMOD_COLAMD:  return "colamd";
    default:              return "unknown";
  }
}

void configure_cholmod_control(cholmod_common *cm, SEXP backend_r)
{
  SEXP ctrl_r, ordering_r, postorder_r;
  const char *ordering;
  int ord_code;

  ctrl_r = getListElement(backend_r, "cholmod_control");
  if(ctrl_r == R_NilValue)
    return;
  if(TYPEOF(ctrl_r) != VECSXP)
    Rf_error("backend$cholmod_control must be a list");

  ordering_r = getListElement(ctrl_r, "ordering");
  postorder_r = getListElement(ctrl_r, "postorder");
  ordering = as_char_scalar(ordering_r, "backend$cholmod_control$ordering");
  cm->postorder = as_flag_scalar(postorder_r, "backend$cholmod_control$postorder");

  if(std::strcmp(ordering, "auto") == 0){
    cm->nmethods = 0;
    return;
  }

  if(std::strcmp(ordering, "best") == 0){
    cm->nmethods = CHOLMOD_MAXMETHODS;
    return;
  }

  if(std::strcmp(ordering, "natural") == 0)
    ord_code = CHOLMOD_NATURAL;
  else if(std::strcmp(ordering, "amd") == 0)
    ord_code = CHOLMOD_AMD;
  else if(std::strcmp(ordering, "metis") == 0)
    ord_code = CHOLMOD_METIS;
  else if(std::strcmp(ordering, "nesdis") == 0)
    ord_code = CHOLMOD_NESDIS;
  else if(std::strcmp(ordering, "colamd") == 0)
    ord_code = CHOLMOD_COLAMD;
  else
    Rf_error("unknown CHOLMOD ordering '%s'", ordering);

  cm->nmethods = 1;
  cm->method[0].ordering = ord_code;
}

int as_int_scalar(SEXP x, const char *where)
{
  if(x == R_NilValue || LENGTH(x) != 1)
    Rf_error("missing or malformed scalar at %s", where);

  if(TYPEOF(x) == INTSXP)
    return INTEGER(x)[0];

  if(TYPEOF(x) == REALSXP)
    return (int)REAL(x)[0];

  Rf_error("expected integer/numeric scalar at %s", where);
  return 0;
}

double as_real_scalar_strict(SEXP x, const char *where)
{
  if(x == R_NilValue || TYPEOF(x) != REALSXP || LENGTH(x) != 1)
    Rf_error("expected numeric scalar at %s", where);
  return REAL(x)[0];
}

const char *as_char_scalar(SEXP x, const char *where)
{
  if(x == R_NilValue || TYPEOF(x) != STRSXP || LENGTH(x) != 1)
    Rf_error("expected character scalar at %s", where);

  return CHAR(STRING_ELT(x, 0));
}

int compute_n_recover(int nSamples, int recoverProcess, int recoverStart, int recoverThin)
{
  int first, last, span;

  if(!recoverProcess)
    return 0;

  if(recoverStart > nSamples)
    return 0;

  first = recoverStart;
  last = nSamples;
  span = last - first;

  return 1 + span / recoverThin;
}

void require_matrix_real(SEXP x, int nr, int nc, const char *where)
{
  SEXP dim;

  if(TYPEOF(x) != REALSXP)
    Rf_error("%s must be numeric", where);

  dim = Rf_getAttrib(x, R_DimSymbol);
  if(dim == R_NilValue || LENGTH(dim) != 2)
    Rf_error("%s must be a matrix", where);

  if(INTEGER(dim)[0] != nr || INTEGER(dim)[1] != nc)
    Rf_error("%s has wrong dimensions", where);
}

void require_matrix_int(SEXP x, int nr, int nc, const char *where)
{
  SEXP dim;

  if(TYPEOF(x) != INTSXP)
    Rf_error("%s must be integer", where);

  dim = Rf_getAttrib(x, R_DimSymbol);
  if(dim == R_NilValue || LENGTH(dim) != 2)
    Rf_error("%s must be a matrix", where);

  if(INTEGER(dim)[0] != nr || INTEGER(dim)[1] != nc)
    Rf_error("%s has wrong dimensions", where);
}

double logDetFactor(cholmod_factor *L)
{
  double det;
  int i, j;
  int *lpi, *lsup, *lpx;
  double *x;
  int nrp1, nc;

  if(!L->is_super)
    Rf_error("factor must be super in logDetFactor");

  det = 0.0;
  lpi = (int*)L->pi;
  lsup = (int*)L->super;
  lpx = (int*)L->px;

  for(i = 0; i < L->nsuper; i++){
    nrp1 = 1 + lpi[i + 1] - lpi[i];
    nc = lsup[i + 1] - lsup[i];
    x = (double*)L->x + lpx[i];

    for(j = 0; j < nc; j++)
      det += 2.0 * std::log(std::fabs(x[j * nrp1]));
  }

  return det;
}

/*==========================================================================*/
/* sampler state structs                                                    */
/*==========================================================================*/
