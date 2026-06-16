#include <math.h>
#include <string.h>

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include "corModels.h"

////////////////////////////////////////////////////////
//// Correlation functions
////////////////////////////////////////////////////////

double cor_exp(const double *theta, double h, double u, int spatialDim){
  (void)u;
  (void)spatialDim;
  double phi = theta[0];
  return exp(-phi*h);
}

double cor_matern(const double *theta, double h, double u, int spatialDim){
  (void)u;
  (void)spatialDim;

  double phi = theta[0];
  double nu = theta[1];
  double x, k_unscaled, log_cor;

  if(h <= 0.0)
    return 1.0;

  x = phi * h;
  if(x <= 0.0)
    return 1.0;

  k_unscaled = bessel_k(x, nu, 1.0);
  if(k_unscaled <= 0.0 || !R_FINITE(k_unscaled))
    return 0.0;

  log_cor = nu * log(x) -
    (nu - 1.0) * log(2.0) -
    lgammafn(nu) +
    log(k_unscaled);

  if(!R_FINITE(log_cor))
    return 0.0;
  if(log_cor < -745.0)
    return 0.0;

  return exp(log_cor);
}

double cor_sep_exp(const double *theta, double h, double u, int spatialDim){
  (void)spatialDim;

  double phi    = theta[0];
  double lambda = theta[1];

  return exp(-phi*h) * exp(-lambda*u);
}

double cor_multires_sep_exp(const double *theta, double h, double u, int spatialDim){
  (void)spatialDim;

  double alpha   = theta[0];
  double phi1    = theta[1];
  double lambda1 = theta[2];
  double phi2    = theta[3];
  double lambda2 = theta[4];

  return alpha*exp(-phi1*h)*exp(-lambda1*u) +
         (1.0-alpha)*exp(-phi2*h)*exp(-lambda2*u);
}

double cor_gneiting(const double *theta, double h, double u, int spatialDim){

  double a     = theta[0];
  double c     = theta[1];
  double alpha = theta[2];
  double beta  = theta[3];
  double gamma = theta[4];
  double delta = theta[5];

  double t = 1.0 + a * pow(fabs(u), 2.0 * alpha);
  double log_cor;

  if(spatialDim < 1)
    spatialDim = 1;

  log_cor = -(delta + 0.5 * (double)spatialDim) * log(t) -
    c * pow(h, 2.0 * gamma) / pow(t, beta * gamma);

  if(!R_FINITE(log_cor))
    return 0.0;
  if(log_cor < -745.0)
    return 0.0;

  return exp(log_cor);
}

////////////////////////////////////////////////////////
//// Parameter metadata
////////////////////////////////////////////////////////

static const char *theta_exp[] = {
  "phi"
};

static const ThetaType theta_exp_types[] = {
  THETA_POS
};

static const char *theta_matern[] = {
  "phi",
  "nu"
};

static const ThetaType theta_matern_types[] = {
  THETA_POS,
  THETA_POS
};

static const char *theta_sep_exp[] = {
  "phi",
  "lambda"
};

static const ThetaType theta_sep_exp_types[] = {
  THETA_POS,
  THETA_POS
};

static const char *theta_multires_sep_exp[] = {
  "alpha",
  "phi_1",
  "lambda_1",
  "phi_2",
  "lambda_2"
};

static const ThetaType theta_multires_sep_exp_types[] = {
  THETA_UNIT,
  THETA_POS,
  THETA_POS,
  THETA_POS,
  THETA_POS
};

static const char *theta_gneiting[] = {
  "a",
  "c",
  "alpha",
  "beta",
  "gamma",
  "delta"
};

static const ThetaType theta_gneiting_types[] = {
  THETA_POS,
  THETA_POS,
  THETA_UNIT,
  THETA_UNIT,
  THETA_UNIT,
  THETA_POS
};

////////////////////////////////////////////////////////
//// Correlation model registry
////////////////////////////////////////////////////////

const CorModelInfo corModelsInfo[] = {

  {
    "exp",
    1,
    theta_exp,
    theta_exp_types,
    COR_SINGLE,
    cor_exp
  },

  {
    "matern",
    2,
    theta_matern,
    theta_matern_types,
    COR_SINGLE,
    cor_matern
  },

  {
    "sep_exp",
    2,
    theta_sep_exp,
    theta_sep_exp_types,
    COR_SPACE_TIME,
    cor_sep_exp
  },

  {
    "multi_res_sep_exp",
    5,
    theta_multires_sep_exp,
    theta_multires_sep_exp_types,
    COR_SPACE_TIME,
    cor_multires_sep_exp
  },

  {
    "gneiting",
    6,
    theta_gneiting,
    theta_gneiting_types,
    COR_SPACE_TIME,
    cor_gneiting
  }

};

const int nCorModels =
  sizeof(corModelsInfo) / sizeof(CorModelInfo);

////////////////////////////////////////////////////////
//// Lookup helpers
////////////////////////////////////////////////////////

int get_cor_model_index(const char *name){

  for(int i=0; i<nCorModels; i++){
    if(strcmp(name, corModelsInfo[i].name) == 0)
      return i;
  }

  return -1;
}

int get_cor_model_nTheta(int indx){

  if(indx < 0 || indx >= nCorModels)
    return -1;

  return corModelsInfo[indx].nTheta;
}

corFunPtr get_cor_fun(int indx){

  if(indx < 0 || indx >= nCorModels)
    return NULL;

  return corModelsInfo[indx].fun;
}

const CorModelInfo* get_cor_model_info(int indx){

  if(indx < 0 || indx >= nCorModels)
    return NULL;

  return &corModelsInfo[indx];
}

const char **get_cor_model_thetaNames(int indx){

  if(indx < 0 || indx >= nCorModels)
    return NULL;

  return corModelsInfo[indx].thetaNames;
}

const ThetaType *get_cor_model_thetaTypes(int indx){

  if(indx < 0 || indx >= nCorModels)
    return NULL;

  return corModelsInfo[indx].thetaTypes;
}

////////////////////////////////////////////////////////
//// R interface
////////////////////////////////////////////////////////

extern "C"
SEXP get_cor_models(){

  SEXP result, names;

  PROTECT(result = Rf_allocVector(VECSXP, nCorModels));
  PROTECT(names  = Rf_allocVector(STRSXP, nCorModels));

  for(int i = 0; i < nCorModels; i++){

    int nTheta = corModelsInfo[i].nTheta;

    SEXP params, types, domains, mode;

    PROTECT(params = Rf_allocVector(STRSXP, nTheta));
    PROTECT(types  = Rf_allocVector(INTSXP, nTheta));
    PROTECT(domains = Rf_allocVector(STRSXP, nTheta));

    for(int j = 0; j < nTheta; j++){

      SET_STRING_ELT(params, j,
        Rf_mkChar(corModelsInfo[i].thetaNames[j]));

      INTEGER(types)[j] =
        corModelsInfo[i].thetaTypes[j];

      if(corModelsInfo[i].thetaTypes[j] == THETA_UNIT)
        SET_STRING_ELT(domains, j, Rf_mkChar("unit"));
      else
        SET_STRING_ELT(domains, j, Rf_mkChar("positive"));
    }

    PROTECT(mode =
      Rf_ScalarInteger(corModelsInfo[i].distanceMode));

    SEXP entry;
    PROTECT(entry = Rf_allocVector(VECSXP,4));

    SET_VECTOR_ELT(entry,0,params);
    SET_VECTOR_ELT(entry,1,types);
    SET_VECTOR_ELT(entry,2,domains);
    SET_VECTOR_ELT(entry,3,mode);

    SEXP entryNames;
    PROTECT(entryNames = Rf_allocVector(STRSXP,4));

    SET_STRING_ELT(entryNames,0,Rf_mkChar("names"));
    SET_STRING_ELT(entryNames,1,Rf_mkChar("types"));
    SET_STRING_ELT(entryNames,2,Rf_mkChar("domains"));
    SET_STRING_ELT(entryNames,3,Rf_mkChar("distance_mode"));

    Rf_namesgets(entry,entryNames);

    SET_VECTOR_ELT(result,i,entry);

    SET_STRING_ELT(names,i,
      Rf_mkChar(corModelsInfo[i].name));

    UNPROTECT(6);
  }

  Rf_namesgets(result,names);

  UNPROTECT(2);

  return result;
}
