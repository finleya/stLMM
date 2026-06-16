#include "stLMM_internal.h"

/* R backend -> C++ sampler state.
 *
 * This file validates the list contract assembled in R/stLMM.R and translates
 * it into SamplerState, GraphState, and TermState structures. It should not
 * build precision matrices or run sampler updates; those belong in the
 * precision/covariance/sampler files.
 */

PriorSpec parse_prior_spec(SEXP x, const char *name, int defaultFamily, double p1, double p2)
{
  PriorSpec out;

  out.family = defaultFamily;
  out.p1 = p1;
  out.p2 = p2;

  if(x == R_NilValue)
    return out;
  if(TYPEOF(x) != REALSXP || LENGTH(x) < 3)
    Rf_error("%s malformed", name);

  out.family = (int)REAL(x)[0];
  out.p1 = REAL(x)[1];
  out.p2 = REAL(x)[2];

  if(out.family < PRIOR_IG || out.family > PRIOR_BETA)
    Rf_error("%s has unknown prior family", name);
  if(!R_FINITE(out.p1) || !R_FINITE(out.p2))
    Rf_error("%s has non-finite prior hyperparameters", name);

  return out;
}

void init_graph_state_from_backend(GraphState *g, SEXP graph_r)
{
  const char *graphTypeName;
  SEXP coords_r;
  SEXP nnIndx_r;
  SEXP nnIndxLU_r;
  SEXP degree_r, edge_i_r, edge_j_r, edge_w_r;
  SEXP time_model_r, time_delta_r;
  SEXP parent_index_r, parent_start_r, parent_count_r;
  int *nnLU;
  int i;

  std::memset(g, 0, sizeof(GraphState));

  graphTypeName = as_char_scalar(getListElement(graph_r, "graph_type"), "graph$graph_type");
  g->nNode = as_int_scalar(getListElement(graph_r, "n"), "graph$n");

  if(std::strcmp(graphTypeName, "nngp") == 0){
    g->type = GRAPH_NNGP;
    
    g->dim = as_int_scalar(getListElement(graph_r, "dim"), "graph$dim");
    coords_r = getListElement(graph_r, "coords_ord");
    require_matrix_real(coords_r, g->nNode, g->dim, "graph$coords_ord");
    g->coords = REAL(coords_r);

    nnIndx_r = getListElement(graph_r, "nnIndx");
    nnIndxLU_r = getListElement(graph_r, "nnIndxLU");
    if(TYPEOF(nnIndx_r) != INTSXP)
      Rf_error("graph$nnIndx must be integer");
    if(TYPEOF(nnIndxLU_r) != INTSXP)
      Rf_error("graph$nnIndxLU must be integer");
    
    g->nnIndx = INTEGER(nnIndx_r);
    nnLU = INTEGER(nnIndxLU_r);
    g->nnStart = (int*)R_alloc(g->nNode, sizeof(int));
    g->nnCount = (int*)R_alloc(g->nNode, sizeof(int));
    
    g->totalNnbr = 0;
    for(i = 0; i < g->nNode; i++){
      g->nnStart[i] = nnLU[i];
      g->nnCount[i] = nnLU[g->nNode + i];
      if(g->nnStart[i] < 0)
        Rf_error("negative neighbor start in graph$nnIndxLU");
      if(g->nnCount[i] < 0)
        Rf_error("negative neighbor count in graph$nnIndxLU");
      if(g->nnStart[i] + g->nnCount[i] > LENGTH(nnIndx_r))
        Rf_error("neighbor start/count exceeds graph$nnIndx length");
      g->totalNnbr += g->nnCount[i];
    }
    if(g->totalNnbr != LENGTH(nnIndx_r))
      Rf_error("graph$nnIndx length does not match graph$nnIndxLU counts");
    for(i = 0; i < g->nNode; i++){
      int k;
      for(k = 0; k < g->nnCount[i]; k++){
        int idx = g->nnIndx[g->nnStart[i] + k];
        if(idx < 0 || idx >= g->nNode)
          Rf_error("graph$nnIndx contains out-of-bounds neighbor index");
        if(idx >= i)
          Rf_error("graph$nnIndx violates NNGP history ordering");
      }
    }
    
  } else if(std::strcmp(graphTypeName, "gp") == 0){
    g->type = GRAPH_GP;
    
    g->dim = as_int_scalar(getListElement(graph_r, "dim"), "graph$dim");
    coords_r = getListElement(graph_r, "coords_ord");
    require_matrix_real(coords_r, g->nNode, g->dim, "graph$coords_ord");
    g->coords = REAL(coords_r);
    
  } else if(std::strcmp(graphTypeName, "ar1") == 0){
    g->type = GRAPH_AR1;
  } else if(std::strcmp(graphTypeName, "car") == 0 || std::strcmp(graphTypeName, "car_time") == 0){
    SEXP car_model_r;
    const char *carModelName;

    if(std::strcmp(graphTypeName, "car") == 0){
      g->type = GRAPH_CAR;
      g->nSpace = g->nNode;
      g->nTime = 1;
    } else {
      const char *timeModelName;
      g->type = GRAPH_CAR_TIME;
      g->nSpace = as_int_scalar(getListElement(graph_r, "n_space"), "graph$n_space");
      g->nTime = as_int_scalar(getListElement(graph_r, "n_time"), "graph$n_time");
      if(g->nSpace < 1 || g->nTime < 1 || g->nSpace * g->nTime != g->nNode)
        Rf_error("graph$n_space and graph$n_time malformed for car_time graph");

      time_model_r = getListElement(graph_r, "time_model");
      if(time_model_r == R_NilValue){
        g->timeModel = TIME_MODEL_AR1;
      } else {
        timeModelName = as_char_scalar(time_model_r, "graph$time_model");
        if(std::strcmp(timeModelName, "ar1") == 0){
          g->timeModel = TIME_MODEL_AR1;
        } else if(std::strcmp(timeModelName, "exp") == 0){
          int d;
          g->timeModel = TIME_MODEL_EXP;
          time_delta_r = getListElement(graph_r, "time_delta");
          if(TYPEOF(time_delta_r) != REALSXP || LENGTH(time_delta_r) != g->nTime - 1)
            Rf_error("graph$time_delta malformed for exp car_time graph");
          g->timeDelta = REAL(time_delta_r);
          for(d = 0; d < g->nTime - 1; d++){
            if(!std::isfinite(g->timeDelta[d]) || g->timeDelta[d] <= 0.0)
              Rf_error("graph$time_delta must contain positive finite gaps for exp car_time graph");
          }
        } else {
          Rf_error("unsupported graph$time_model '%s' for car_time graph", timeModelName);
        }
      }
    }

    car_model_r = getListElement(graph_r, "car_model");
    if(car_model_r == R_NilValue){
      g->carModel = CAR_MODEL_PROPER;
    } else {
      carModelName = as_char_scalar(car_model_r, "graph$car_model");
      if(std::strcmp(carModelName, "proper") == 0)
        g->carModel = CAR_MODEL_PROPER;
      else if(std::strcmp(carModelName, "leroux") == 0)
        g->carModel = CAR_MODEL_LEROUX;
      else
        Rf_error("unsupported graph$car_model '%s' for car/car_time graph", carModelName);
    }

    degree_r = getListElement(graph_r, "degree");
    edge_i_r = getListElement(graph_r, "edge_i");
    edge_j_r = getListElement(graph_r, "edge_j");
    edge_w_r = getListElement(graph_r, "edge_w");
    g->nEdge = as_int_scalar(getListElement(graph_r, "n_edge"), "graph$n_edge");

    if(TYPEOF(degree_r) != REALSXP || LENGTH(degree_r) != g->nSpace)
      Rf_error("graph$degree malformed for car/car_time graph");
    if(TYPEOF(edge_i_r) != INTSXP || LENGTH(edge_i_r) != g->nEdge)
      Rf_error("graph$edge_i malformed for car graph");
    if(TYPEOF(edge_j_r) != INTSXP || LENGTH(edge_j_r) != g->nEdge)
      Rf_error("graph$edge_j malformed for car graph");
    if(TYPEOF(edge_w_r) != REALSXP || LENGTH(edge_w_r) != g->nEdge)
      Rf_error("graph$edge_w malformed for car graph");

    g->degree = REAL(degree_r);
    g->edgeI = (int*)R_alloc(g->nEdge > 0 ? g->nEdge : 1, sizeof(int));
    g->edgeJ = (int*)R_alloc(g->nEdge > 0 ? g->nEdge : 1, sizeof(int));
    g->edgeW = REAL(edge_w_r);
    for(i = 0; i < g->nEdge; i++){
      g->edgeI[i] = INTEGER(edge_i_r)[i] - 1;
      g->edgeJ[i] = INTEGER(edge_j_r)[i] - 1;
      if(g->edgeI[i] < 0 || g->edgeI[i] >= g->nSpace ||
         g->edgeJ[i] < 0 || g->edgeJ[i] >= g->nSpace)
        Rf_error("car graph edge index out of bounds");
    }
  } else if(std::strcmp(graphTypeName, "dagar") == 0 || std::strcmp(graphTypeName, "dagar_time") == 0){
    if(std::strcmp(graphTypeName, "dagar") == 0){
      g->type = GRAPH_DAGAR;
      g->nSpace = g->nNode;
      g->nTime = 1;
    } else {
      const char *timeModelName;
      g->type = GRAPH_DAGAR_TIME;
      g->nSpace = as_int_scalar(getListElement(graph_r, "n_space"), "graph$n_space");
      g->nTime = as_int_scalar(getListElement(graph_r, "n_time"), "graph$n_time");
      if(g->nSpace < 1 || g->nTime < 1 || g->nSpace * g->nTime != g->nNode)
        Rf_error("graph$n_space and graph$n_time malformed for dagar_time graph");

      time_model_r = getListElement(graph_r, "time_model");
      if(time_model_r == R_NilValue){
        g->timeModel = TIME_MODEL_AR1;
      } else {
        timeModelName = as_char_scalar(time_model_r, "graph$time_model");
        if(std::strcmp(timeModelName, "ar1") == 0){
          g->timeModel = TIME_MODEL_AR1;
        } else if(std::strcmp(timeModelName, "exp") == 0){
          int d;
          g->timeModel = TIME_MODEL_EXP;
          time_delta_r = getListElement(graph_r, "time_delta");
          if(TYPEOF(time_delta_r) != REALSXP || LENGTH(time_delta_r) != g->nTime - 1)
            Rf_error("graph$time_delta malformed for exp dagar_time graph");
          g->timeDelta = REAL(time_delta_r);
          for(d = 0; d < g->nTime - 1; d++){
            if(!std::isfinite(g->timeDelta[d]) || g->timeDelta[d] <= 0.0)
              Rf_error("graph$time_delta must contain positive finite gaps for exp dagar_time graph");
          }
        } else {
          Rf_error("unsupported graph$time_model '%s' for dagar_time graph", timeModelName);
        }
      }
    }

    parent_index_r = getListElement(graph_r, "parent_index");
    parent_start_r = getListElement(graph_r, "parent_start");
    parent_count_r = getListElement(graph_r, "parent_count");

    if(TYPEOF(parent_index_r) != INTSXP)
      Rf_error("graph$parent_index must be integer for dagar graph");
    if(TYPEOF(parent_start_r) != INTSXP || LENGTH(parent_start_r) != g->nSpace)
      Rf_error("graph$parent_start malformed for dagar graph");
    if(TYPEOF(parent_count_r) != INTSXP || LENGTH(parent_count_r) != g->nSpace)
      Rf_error("graph$parent_count malformed for dagar graph");

    g->parentIndx = INTEGER(parent_index_r);
    g->parentStart = (int*)R_alloc(g->nSpace, sizeof(int));
    g->parentCount = (int*)R_alloc(g->nSpace, sizeof(int));
    g->totalParent = 0;

    for(i = 0; i < g->nSpace; i++){
      int k;
      g->parentStart[i] = INTEGER(parent_start_r)[i];
      g->parentCount[i] = INTEGER(parent_count_r)[i];
      if(g->parentStart[i] < 0)
        Rf_error("negative parent start in dagar graph");
      if(g->parentCount[i] < 0)
        Rf_error("negative parent count in dagar graph");
      if(g->parentStart[i] + g->parentCount[i] > LENGTH(parent_index_r))
        Rf_error("parent start/count exceeds graph$parent_index length");
      g->totalParent += g->parentCount[i];
      for(k = 0; k < g->parentCount[i]; k++){
        int idx = g->parentIndx[g->parentStart[i] + k];
        if(idx < 0 || idx >= g->nSpace)
          Rf_error("dagar graph parent index out of bounds");
        if(idx >= i)
          Rf_error("dagar graph parent index violates DAG ordering");
      }
    }
    if(g->totalParent != LENGTH(parent_index_r))
      Rf_error("graph$parent_index length does not match parent counts");
  } else {
    Rf_error("unsupported graph_type '%s'", graphTypeName);
  }
}

void init_term_state_from_backend(TermState *term, GraphState *graphs, int nGraphs, SEXP term_r, int n)
{
  SEXP map_r, obsIndx_r, node_nobs_r, obsLU_r, x_r;
  SEXP sigma_start_r, sigma_tune_r, sigma_ig_r, sigma_prior_r;
  SEXP theta_start_r, theta_tune_r, theta_bounds_r, theta_bounds_dim_r;
  SEXP theta_prior_r, theta_prior_dim_r;
  const char *termTypeName, *covModelName;
  GraphState *g;
  int *obsLU;
  int i;

  std::memset(term, 0, sizeof(TermState));

  term->name = as_char_scalar(getListElement(term_r, "name"), "term$name");
  term->label = as_char_scalar(getListElement(term_r, "label"), "term$label");

  termTypeName = as_char_scalar(getListElement(term_r, "term_type"), "term$term_type");
  if(std::strcmp(termTypeName, "nngp") != 0 &&
     std::strcmp(termTypeName, "gp")   != 0 &&
     std::strcmp(termTypeName, "ar1")  != 0 &&
     std::strcmp(termTypeName, "car")  != 0 &&
     std::strcmp(termTypeName, "dagar")  != 0 &&
     std::strcmp(termTypeName, "dagar_time")  != 0 &&
     std::strcmp(termTypeName, "car_time") != 0)
    Rf_error("unsupported term_type '%s'", termTypeName);

  term->graphIndex = as_int_scalar(getListElement(term_r, "graph_index"), "term$graph_index") - 1;
  if(term->graphIndex < 0 || term->graphIndex >= nGraphs)
    Rf_error("term$graph_index out of bounds");
  
  term->nObs = as_int_scalar(getListElement(term_r, "n_obs"), "term$n_obs");
  term->nNode = as_int_scalar(getListElement(term_r, "n_node"), "term$n_node");
  term->qLat = term->nNode;
  term->wOffset = 0;
  
  g = graphs + term->graphIndex;
  
  if(term->nNode != g->nNode)
    Rf_error("term$n_node does not match linked graph node count");

  term->sigmaSq = 1.0;
  term->sigmaSqTune = 0.1;
  term->sigmaSqShape = 2.0;
  term->sigmaSqScale = 1.0;
  term->sigmaSqPrior.family = PRIOR_IG;
  term->sigmaSqPrior.p1 = 2.0;
  term->sigmaSqPrior.p2 = 1.0;

  /* indexing */
  map_r = getListElement(term_r, "map");
  obsIndx_r = getListElement(term_r, "obsIndx");
  node_nobs_r = getListElement(term_r, "node_nobs");
  obsLU_r = getListElement(term_r, "obsIndxLU");

  if(TYPEOF(map_r) != INTSXP || LENGTH(map_r) != term->nObs)
    Rf_error("term$map malformed");
  if(TYPEOF(obsIndx_r) != INTSXP || LENGTH(obsIndx_r) != term->nObs)
    Rf_error("term$obsIndx malformed");
  if(TYPEOF(node_nobs_r) != INTSXP || LENGTH(node_nobs_r) != term->nNode)
    Rf_error("term$node_nobs malformed");

  require_matrix_int(obsLU_r, term->nNode, 2, "term$obsIndxLU");
  obsLU = INTEGER(obsLU_r);

  term->map = (int*)R_alloc(term->nObs, sizeof(int));
  term->obsIndx = (int*)R_alloc(term->nObs, sizeof(int));
  term->nodeNobs = (int*)R_alloc(term->nNode, sizeof(int));
  term->obsStart = (int*)R_alloc(term->nNode, sizeof(int));
  term->obsEnd = (int*)R_alloc(term->nNode, sizeof(int));
  term->scale = (double*)R_alloc(n, sizeof(double));

  for(i = 0; i < term->nObs; i++){
    term->map[i] = INTEGER(map_r)[i] - 1;
    term->obsIndx[i] = INTEGER(obsIndx_r)[i] - 1;
  }

  for(i = 0; i < term->nNode; i++){
    term->nodeNobs[i] = INTEGER(node_nobs_r)[i];
    term->obsStart[i] = obsLU[i] - 1;
    term->obsEnd[i] = obsLU[i + term->nNode] - 1;
  }

  /* covariate */
  x_r = getListElement(term_r, "x");
  if(x_r != R_NilValue){
    if(TYPEOF(x_r) != REALSXP || LENGTH(x_r) != n)
      Rf_error("term$x malformed");
    term->x = REAL(x_r);
    term->type = TERM_SVC;
  } else {
    term->x = NULL;
    term->type = TERM_INTERCEPT;
  }

  for(i = 0; i < n; i++)
    term->scale[i] = (term->x != NULL) ? term->x[i] : 1.0;

  /* sigmaSq controls */
  sigma_start_r = getListElement(term_r, "sigma_sq_starting");
  sigma_tune_r  = getListElement(term_r, "sigma_sq_tuning");
  sigma_ig_r    = getListElement(term_r, "sigma_sq_IG");
  sigma_prior_r = getListElement(term_r, "sigma_sq_prior");

  if(sigma_start_r != R_NilValue)
    term->sigmaSq = as_real_scalar_strict(sigma_start_r, "term$sigma_sq_starting");
  if(!(term->sigmaSq > 0.0))
    Rf_error("term$sigma_sq_starting must be positive");

  if(sigma_tune_r != R_NilValue)
    term->sigmaSqTune = as_real_scalar_strict(sigma_tune_r, "term$sigma_sq_tuning");
  if(term->sigmaSqTune < 0.0)
    Rf_error("term$sigma_sq_tuning must be nonnegative");

  if(sigma_ig_r != R_NilValue){
    if(TYPEOF(sigma_ig_r) != REALSXP || LENGTH(sigma_ig_r) != 2)
      Rf_error("term$sigma_sq_IG malformed");
    term->sigmaSqShape = REAL(sigma_ig_r)[0];
    term->sigmaSqScale = REAL(sigma_ig_r)[1];
  }
  if(!(term->sigmaSqShape > 0.0 && term->sigmaSqScale > 0.0))
    Rf_error("term$sigma_sq_IG parameters must be positive");
  term->sigmaSqPrior = parse_prior_spec(sigma_prior_r, "term$sigma_sq_prior",
                                        PRIOR_IG, term->sigmaSqShape, term->sigmaSqScale);

  /* NNGP and GP share the same correlation-model / theta parsing */
  if(g->type == GRAPH_NNGP || g->type == GRAPH_GP){

    covModelName = as_char_scalar(getListElement(term_r, "cov_model"), "term$cov_model");
    term->covModelIndex = get_cor_model_index(covModelName);
    if(term->covModelIndex < 0)
      Rf_error("unknown correlation model '%s'", covModelName);

    term->thetaDim = get_cor_model_nTheta(term->covModelIndex);
    term->corFun = get_cor_fun(term->covModelIndex);
    term->distanceMode = get_cor_model_info(term->covModelIndex)->distanceMode;

    term->theta = (double*)R_alloc(term->thetaDim, sizeof(double));
    term->thetaTune = (double*)R_alloc(term->thetaDim, sizeof(double));
    term->thetaLower = (double*)R_alloc(term->thetaDim, sizeof(double));
    term->thetaUpper = (double*)R_alloc(term->thetaDim, sizeof(double));
    term->thetaPrior = (PriorSpec*)R_alloc(term->thetaDim, sizeof(PriorSpec));

    theta_start_r  = getListElement(term_r, "theta_starting");
    theta_tune_r   = getListElement(term_r, "theta_tuning");
    theta_bounds_r = getListElement(term_r, "theta_bounds");
    theta_prior_r  = getListElement(term_r, "theta_prior");

    if(theta_tune_r == R_NilValue || TYPEOF(theta_tune_r) != REALSXP || LENGTH(theta_tune_r) != term->thetaDim)
      Rf_error(g->type == GRAPH_NNGP ?
               "term$theta_tuning malformed for nngp term" :
               "term$theta_tuning malformed for gp term");

    if(theta_bounds_r == R_NilValue || TYPEOF(theta_bounds_r) != REALSXP)
      Rf_error(g->type == GRAPH_NNGP ?
               "term$theta_bounds malformed for nngp term" :
               "term$theta_bounds malformed for gp term");

    theta_bounds_dim_r = Rf_getAttrib(theta_bounds_r, R_DimSymbol);
    if(theta_bounds_dim_r == R_NilValue || TYPEOF(theta_bounds_dim_r) != INTSXP || LENGTH(theta_bounds_dim_r) != 2)
      Rf_error(g->type == GRAPH_NNGP ?
               "term$theta_bounds must be a matrix for nngp term" :
               "term$theta_bounds must be a matrix for gp term");

    if(INTEGER(theta_bounds_dim_r)[0] != term->thetaDim || INTEGER(theta_bounds_dim_r)[1] != 2)
      Rf_error(g->type == GRAPH_NNGP ?
               "term$theta_bounds has wrong dimensions for nngp term" :
               "term$theta_bounds has wrong dimensions for gp term");

    if(theta_prior_r == R_NilValue || TYPEOF(theta_prior_r) != REALSXP)
      Rf_error(g->type == GRAPH_NNGP ?
               "term$theta_prior malformed for nngp term" :
               "term$theta_prior malformed for gp term");

    theta_prior_dim_r = Rf_getAttrib(theta_prior_r, R_DimSymbol);
    if(theta_prior_dim_r == R_NilValue || TYPEOF(theta_prior_dim_r) != INTSXP || LENGTH(theta_prior_dim_r) != 2)
      Rf_error(g->type == GRAPH_NNGP ?
               "term$theta_prior must be a matrix for nngp term" :
               "term$theta_prior must be a matrix for gp term");

    if(INTEGER(theta_prior_dim_r)[0] != term->thetaDim || INTEGER(theta_prior_dim_r)[1] != 6)
      Rf_error(g->type == GRAPH_NNGP ?
               "term$theta_prior has wrong dimensions for nngp term" :
               "term$theta_prior has wrong dimensions for gp term");

    for(i = 0; i < term->thetaDim; i++){
      term->thetaLower[i] = REAL(theta_bounds_r)[i];
      term->thetaUpper[i] = REAL(theta_bounds_r)[i + term->thetaDim];
      term->thetaPrior[i].family = (int)REAL(theta_prior_r)[i];
      term->thetaPrior[i].p1 = REAL(theta_prior_r)[i + term->thetaDim];
      term->thetaPrior[i].p2 = REAL(theta_prior_r)[i + 2 * term->thetaDim];
    }

    if(theta_start_r == R_NilValue || TYPEOF(theta_start_r) != REALSXP || LENGTH(theta_start_r) != term->thetaDim)
      Rf_error(g->type == GRAPH_NNGP ?
               "term$theta_starting malformed for nngp term" :
               "term$theta_starting malformed for gp term");

    for(i = 0; i < term->thetaDim; i++){
      term->theta[i] = REAL(theta_start_r)[i];
      term->thetaTune[i] = REAL(theta_tune_r)[i];

      if(term->thetaTune[i] > 0.0 &&
         !(term->thetaLower[i] < term->theta[i] && term->theta[i] < term->thetaUpper[i])){
        Rf_error(g->type == GRAPH_NNGP ?
                 "term$theta_starting must lie strictly inside bounds for nngp term" :
                 "term$theta_starting must lie strictly inside bounds for gp term");
      }
    }

    if(g->type == GRAPH_NNGP){
      term->B = (double*)R_alloc(g->totalNnbr > 0 ? g->totalNnbr : 1, sizeof(double));
      term->F = (double*)R_alloc(term->nNode > 0 ? term->nNode : 1, sizeof(double));

      for(i = 0; i < g->totalNnbr; i++)
        term->B[i] = 0.0;
      for(i = 0; i < term->nNode; i++)
        term->F[i] = 1.0;
    }

  } else if(g->type == GRAPH_AR1){

    theta_start_r  = getListElement(term_r, "theta_starting");
    theta_tune_r   = getListElement(term_r, "theta_tuning");
    theta_bounds_r = getListElement(term_r, "theta_bounds");
    theta_prior_r  = getListElement(term_r, "theta_prior");

    term->covModelIndex = -1;
    term->thetaDim = 1;
    term->distanceMode = 0;

    term->theta = (double*)R_alloc(1, sizeof(double));
    term->thetaTune = (double*)R_alloc(1, sizeof(double));
    term->thetaLower = (double*)R_alloc(1, sizeof(double));
    term->thetaUpper = (double*)R_alloc(1, sizeof(double));
    term->thetaPrior = (PriorSpec*)R_alloc(1, sizeof(PriorSpec));

    term->theta[0] = as_real_scalar_strict(theta_start_r, "term$theta_starting");
    term->thetaTune[0] = as_real_scalar_strict(theta_tune_r, "term$theta_tuning");

    if(theta_bounds_r == R_NilValue || TYPEOF(theta_bounds_r) != REALSXP)
      Rf_error("term$theta_bounds malformed for ar1 term");
    theta_bounds_dim_r = Rf_getAttrib(theta_bounds_r, R_DimSymbol);
    if(theta_bounds_dim_r == R_NilValue || TYPEOF(theta_bounds_dim_r) != INTSXP || LENGTH(theta_bounds_dim_r) != 2)
      Rf_error("term$theta_bounds must be a matrix for ar1 term");
    if(INTEGER(theta_bounds_dim_r)[0] != 1 || INTEGER(theta_bounds_dim_r)[1] != 2)
      Rf_error("term$theta_bounds has wrong dimensions for ar1 term");

    term->thetaLower[0] = REAL(theta_bounds_r)[0];
    term->thetaUpper[0] = REAL(theta_bounds_r)[1];
    term->thetaPrior[0] = parse_prior_spec(theta_prior_r, "term$theta_prior",
                                           PRIOR_UNIFORM, 0.0, 0.0);

    if(term->thetaTune[0] > 0.0 &&
       !(term->thetaLower[0] < term->theta[0] && term->theta[0] < term->thetaUpper[0]))
      Rf_error("term$theta_starting must lie strictly inside bounds for ar1 term");
  }
  else if(g->type == GRAPH_CAR || g->type == GRAPH_DAGAR){
    theta_start_r  = getListElement(term_r, "theta_starting");
    theta_tune_r   = getListElement(term_r, "theta_tuning");
    theta_bounds_r = getListElement(term_r, "theta_bounds");
    theta_prior_r  = getListElement(term_r, "theta_prior");

    term->thetaDim = 1;
    term->theta = (double*)R_alloc(1, sizeof(double));
    term->thetaTune = (double*)R_alloc(1, sizeof(double));
    term->thetaLower = (double*)R_alloc(1, sizeof(double));
    term->thetaUpper = (double*)R_alloc(1, sizeof(double));
    term->thetaPrior = (PriorSpec*)R_alloc(1, sizeof(PriorSpec));

    term->theta[0] = as_real_scalar_strict(theta_start_r, "term$theta_starting");
    term->thetaTune[0] = as_real_scalar_strict(theta_tune_r, "term$theta_tuning");

    if(theta_bounds_r == R_NilValue || TYPEOF(theta_bounds_r) != REALSXP)
      Rf_error(g->type == GRAPH_DAGAR ?
               "term$theta_bounds malformed for dagar term" :
               "term$theta_bounds malformed for car term");
    theta_bounds_dim_r = Rf_getAttrib(theta_bounds_r, R_DimSymbol);
    if(theta_bounds_dim_r == R_NilValue || TYPEOF(theta_bounds_dim_r) != INTSXP || LENGTH(theta_bounds_dim_r) != 2)
      Rf_error(g->type == GRAPH_DAGAR ?
               "term$theta_bounds must be a matrix for dagar term" :
               "term$theta_bounds must be a matrix for car term");
    if(INTEGER(theta_bounds_dim_r)[0] != 1 || INTEGER(theta_bounds_dim_r)[1] != 2)
      Rf_error(g->type == GRAPH_DAGAR ?
               "term$theta_bounds has wrong dimensions for dagar term" :
               "term$theta_bounds has wrong dimensions for car term");

    term->thetaLower[0] = REAL(theta_bounds_r)[0];
    term->thetaUpper[0] = REAL(theta_bounds_r)[1];
    term->thetaPrior[0] = parse_prior_spec(theta_prior_r, "term$theta_prior",
                                           PRIOR_UNIFORM, 0.0, 0.0);

    if(term->thetaTune[0] > 0.0 &&
       !(term->thetaLower[0] < term->theta[0] && term->theta[0] < term->thetaUpper[0]))
      Rf_error(g->type == GRAPH_DAGAR ?
               "term$theta_starting must lie strictly inside bounds for dagar term" :
               "term$theta_starting must lie strictly inside bounds for car term");
    if(term->theta[0] <= 0.0 || term->theta[0] >= 1.0)
      Rf_error(g->type == GRAPH_DAGAR ?
               "term$theta_starting rho must lie in (0,1) for dagar term" :
               "term$theta_starting rho must lie in (0,1) for car term");
  }
  else if(g->type == GRAPH_CAR_TIME || g->type == GRAPH_DAGAR_TIME){
    theta_start_r  = getListElement(term_r, "theta_starting");
    theta_tune_r   = getListElement(term_r, "theta_tuning");
    theta_bounds_r = getListElement(term_r, "theta_bounds");
    theta_prior_r  = getListElement(term_r, "theta_prior");
    const char *time_term_name = (g->type == GRAPH_DAGAR_TIME) ? "dagar_time" : "car_time";

    term->covModelIndex = -1;
    term->thetaDim = 2;
    term->distanceMode = 0;

    term->theta = (double*)R_alloc(2, sizeof(double));
    term->thetaTune = (double*)R_alloc(2, sizeof(double));
    term->thetaLower = (double*)R_alloc(2, sizeof(double));
    term->thetaUpper = (double*)R_alloc(2, sizeof(double));
    term->thetaPrior = (PriorSpec*)R_alloc(2, sizeof(PriorSpec));

    if(theta_start_r == R_NilValue || TYPEOF(theta_start_r) != REALSXP || LENGTH(theta_start_r) != 2)
      Rf_error("term$theta_starting malformed for %s term", time_term_name);
    if(theta_tune_r == R_NilValue || TYPEOF(theta_tune_r) != REALSXP || LENGTH(theta_tune_r) != 2)
      Rf_error("term$theta_tuning malformed for %s term", time_term_name);
    if(theta_bounds_r == R_NilValue || TYPEOF(theta_bounds_r) != REALSXP)
      Rf_error("term$theta_bounds malformed for %s term", time_term_name);
    theta_bounds_dim_r = Rf_getAttrib(theta_bounds_r, R_DimSymbol);
    if(theta_bounds_dim_r == R_NilValue || TYPEOF(theta_bounds_dim_r) != INTSXP || LENGTH(theta_bounds_dim_r) != 2)
      Rf_error("term$theta_bounds must be a matrix for %s term", time_term_name);
    if(INTEGER(theta_bounds_dim_r)[0] != 2 || INTEGER(theta_bounds_dim_r)[1] != 2)
      Rf_error("term$theta_bounds has wrong dimensions for %s term", time_term_name);
    if(theta_prior_r == R_NilValue || TYPEOF(theta_prior_r) != REALSXP)
      Rf_error("term$theta_prior malformed for %s term", time_term_name);
    theta_prior_dim_r = Rf_getAttrib(theta_prior_r, R_DimSymbol);
    if(theta_prior_dim_r == R_NilValue || TYPEOF(theta_prior_dim_r) != INTSXP || LENGTH(theta_prior_dim_r) != 2)
      Rf_error("term$theta_prior must be a matrix for %s term", time_term_name);
    if(INTEGER(theta_prior_dim_r)[0] != 2 || INTEGER(theta_prior_dim_r)[1] != 6)
      Rf_error("term$theta_prior has wrong dimensions for %s term", time_term_name);

    for(i = 0; i < 2; i++){
      term->theta[i] = REAL(theta_start_r)[i];
      term->thetaTune[i] = REAL(theta_tune_r)[i];
      term->thetaLower[i] = REAL(theta_bounds_r)[i];
      term->thetaUpper[i] = REAL(theta_bounds_r)[i + 2];
      term->thetaPrior[i].family = (int)REAL(theta_prior_r)[i];
      term->thetaPrior[i].p1 = REAL(theta_prior_r)[i + 2];
      term->thetaPrior[i].p2 = REAL(theta_prior_r)[i + 4];
      if(term->thetaTune[i] > 0.0 &&
         !(term->thetaLower[i] < term->theta[i] && term->theta[i] < term->thetaUpper[i]))
        Rf_error("term$theta_starting must lie strictly inside bounds for %s term", time_term_name);
    }
    if(term->theta[0] <= 0.0 || term->theta[0] >= 1.0)
      Rf_error("term$theta_starting rho must lie in (0,1) for %s term", time_term_name);
    if(g->timeModel == TIME_MODEL_AR1){
      if(std::fabs(term->theta[1]) >= 1.0)
        Rf_error("term$theta_starting phi must lie in (-1,1) for %s term", time_term_name);
    } else if(g->timeModel == TIME_MODEL_EXP){
      if(term->theta[1] <= 0.0)
        Rf_error("term$theta_starting lambda must be positive for exp %s term", time_term_name);
    }
  }
}

void init_sampler_state(SamplerState *s, SEXP backend_r)
{
  SEXP graphs_r;
  SEXP terms_r;
  SEXP Z_r;
  SEXP y_r;
  SEXP X_r;
  SEXP n_obs_r;
  SEXP family_r;
  int i;

  std::memset(s, 0, sizeof(SamplerState));

  n_obs_r = getListElement(backend_r, "n_obs");
  if(n_obs_r == R_NilValue)
    s->n = as_int_scalar(getListElement(backend_r, "n"), "backend$n");
  else
    s->n = as_int_scalar(n_obs_r, "backend$n_obs");
  s->p = as_int_scalar(getListElement(backend_r, "p"), "backend$p");
  s->q = as_int_scalar(getListElement(backend_r, "q"), "backend$q");
  s->qLatTotal = 0;
  s->nOmpThreads = 1;
  if(getListElement(backend_r, "n_omp_threads") != R_NilValue){
    s->nOmpThreads = as_int_scalar(getListElement(backend_r, "n_omp_threads"), "backend$n_omp_threads");
    if(s->nOmpThreads < 1)
      Rf_error("backend$n_omp_threads must be >= 1");
  }
#ifndef _OPENMP
  if(s->nOmpThreads > 1){
    Rf_warning("n_omp_threads > 1, but source was not compiled with OpenMP support");
    s->nOmpThreads = 1;
  }
#endif
  s->likelihoodFamily = LIKELIHOOD_GAUSSIAN;
  s->yObserved = NULL;
  s->offset = NULL;
  s->hasOffset = 0;
  s->nTrial = NULL;
  s->nbSize = 0.0;

  family_r = getListElement(backend_r, "family");
  if(family_r != R_NilValue){
    const char *familyName = as_char_scalar(family_r, "backend$family");
    if(std::strcmp(familyName, "gaussian") == 0){
      s->likelihoodFamily = LIKELIHOOD_GAUSSIAN;
    } else if(std::strcmp(familyName, "binomial") == 0){
      s->likelihoodFamily = LIKELIHOOD_BINOMIAL;
    } else if(std::strcmp(familyName, "negative_binomial") == 0){
      s->likelihoodFamily = LIKELIHOOD_NEGATIVE_BINOMIAL;
    } else {
      Rf_error("unsupported backend$family");
    }
  }

  if(getListElement(backend_r, "n_samples") == R_NilValue)
    s->nSamples = 1;
  else
    s->nSamples = as_int_scalar(getListElement(backend_r, "n_samples"), "backend$n_samples");


  {
    SEXP recover_process_r;
    SEXP recover_start_r;
    SEXP recover_thin_r;

    recover_process_r = getListElement(backend_r, "recover_process");
    recover_start_r   = getListElement(backend_r, "recover_start");
    recover_thin_r    = getListElement(backend_r, "recover_thin");

    if(recover_process_r == R_NilValue)
      s->recoverProcess = 0;
    else
      s->recoverProcess = as_flag_scalar(recover_process_r, "backend$recover_process");

    if(recover_start_r == R_NilValue)
      s->recoverStart = 1;
    else
      s->recoverStart = as_int_scalar(recover_start_r, "backend$recover_start");

    if(recover_thin_r == R_NilValue)
      s->recoverThin = 1;
    else
      s->recoverThin = as_int_scalar(recover_thin_r, "backend$recover_thin");

    if(s->recoverStart < 1)
      Rf_error("backend$recover_start must be >= 1");

    if(s->recoverThin < 1)
      Rf_error("backend$recover_thin must be >= 1");

    s->nRecover = compute_n_recover(s->nSamples,
                                    s->recoverProcess,
                                    s->recoverStart,
                                    s->recoverThin);
  }

  {
    SEXP warmup_r, enabled_r, batch_length_r, min_batches_r, max_batches_r, target_r, near_zero_r;

    s->warmupEnabled = 1;
    s->warmupBatchLength = 25;
    s->warmupMinBatches = 0;
    s->warmupMaxBatches = 20;
    s->warmupTargetLower = 0.15;
    s->warmupTargetUpper = 0.45;
    s->warmupNearZero = 0.02;

    warmup_r = getListElement(backend_r, "warmup");
    if(warmup_r != R_NilValue){
      if(TYPEOF(warmup_r) != VECSXP)
        Rf_error("backend$warmup must be a list");

      enabled_r = getListElement(warmup_r, "enabled");
      batch_length_r = getListElement(warmup_r, "batch_length");
      min_batches_r = getListElement(warmup_r, "min_batches");
      max_batches_r = getListElement(warmup_r, "max_batches");
      target_r = getListElement(warmup_r, "target");
      near_zero_r = getListElement(warmup_r, "near_zero");

      if(enabled_r != R_NilValue)
        s->warmupEnabled = as_flag_scalar(enabled_r, "backend$warmup$enabled");
      if(batch_length_r != R_NilValue)
        s->warmupBatchLength = as_int_scalar(batch_length_r, "backend$warmup$batch_length");
      if(min_batches_r != R_NilValue)
        s->warmupMinBatches = as_nonneg_int_scalar(min_batches_r, "backend$warmup$min_batches");
      if(max_batches_r != R_NilValue)
        s->warmupMaxBatches = as_nonneg_int_scalar(max_batches_r, "backend$warmup$max_batches");
      if(target_r != R_NilValue){
        if(TYPEOF(target_r) != REALSXP || LENGTH(target_r) != 2)
          Rf_error("backend$warmup$target must be numeric length 2");
        s->warmupTargetLower = REAL(target_r)[0];
        s->warmupTargetUpper = REAL(target_r)[1];
      }
      if(near_zero_r != R_NilValue)
        s->warmupNearZero = as_real_scalar_strict(near_zero_r, "backend$warmup$near_zero");
    }

    if(s->warmupBatchLength < 1)
      Rf_error("backend$warmup$batch_length must be positive");
    if(s->warmupMinBatches > s->warmupMaxBatches)
      Rf_error("backend$warmup$min_batches must be less than or equal to backend$warmup$max_batches");
    if(!(s->warmupTargetLower > 0.0 && s->warmupTargetLower < s->warmupTargetUpper &&
         s->warmupTargetUpper < 1.0))
      Rf_error("backend$warmup$target must satisfy 0 < lower < upper < 1");
    if(!(s->warmupNearZero >= 0.0 && s->warmupNearZero < s->warmupTargetLower))
      Rf_error("backend$warmup$near_zero must be nonnegative and less than target lower");
  }

  {
    SEXP blocking_r, batch_length_r, target_accept_r;

    s->metropolisBlocking = 0;
    s->metropolisBatchLength = 25;
    s->metropolisTargetAccept = 0.234;
    blocking_r = getListElement(backend_r, "metropolis_blocking");
    if(blocking_r != R_NilValue)
      s->metropolisBlocking = as_int_scalar(blocking_r, "backend$metropolis_blocking");
    if(s->metropolisBlocking < 0 || s->metropolisBlocking > 5)
      Rf_error("backend$metropolis_blocking has unknown blocking code");
    if(s->metropolisBlocking == 5)
      s->metropolisTargetAccept = 0.44;
    batch_length_r = getListElement(backend_r, "metropolis_batch_length");
    if(batch_length_r != R_NilValue)
      s->metropolisBatchLength = as_int_scalar(batch_length_r, "backend$metropolis_batch_length");
    target_accept_r = getListElement(backend_r, "metropolis_target_accept");
    if(target_accept_r != R_NilValue)
      s->metropolisTargetAccept = as_real_scalar_strict(target_accept_r, "backend$metropolis_target_accept");
    if(s->metropolisBatchLength < 1)
      Rf_error("backend$metropolis_batch_length must be positive");
    if(!(s->metropolisTargetAccept > 0.0 && s->metropolisTargetAccept < 1.0))
      Rf_error("backend$metropolis_target_accept must satisfy 0 < target < 1");
  }

  {
    SEXP offset_r = getListElement(backend_r, "offset_obs");
    SEXP has_offset_r = getListElement(backend_r, "has_offset");
    if(offset_r != R_NilValue){
      if(TYPEOF(offset_r) != REALSXP || LENGTH(offset_r) != s->n)
        Rf_error("backend$offset_obs malformed");
      s->offset = REAL(offset_r);
      if(has_offset_r != R_NilValue)
        s->hasOffset = as_flag_scalar(has_offset_r, "backend$has_offset");
      else
        s->hasOffset = 1;
    }
  }

  y_r = getListElement(backend_r, "y_model_obs");
  if(y_r == R_NilValue)
    y_r = getListElement(backend_r, "y_obs");
  if(y_r == R_NilValue)
    y_r = getListElement(backend_r, "y");
  if(TYPEOF(y_r) != REALSXP || LENGTH(y_r) != s->n)
    Rf_error("backend$y_obs malformed");
  if(s->likelihoodFamily == LIKELIHOOD_BINOMIAL){
    SEXP trials_r = getListElement(backend_r, "trials_obs");
    if(TYPEOF(trials_r) != INTSXP || LENGTH(trials_r) != s->n)
      Rf_error("backend$trials_obs malformed for binomial likelihood");

    s->yObserved = REAL(y_r);
    s->y = (double*)R_alloc(s->n > 0 ? s->n : 1, sizeof(double));
    s->nTrial = (int*)R_alloc(s->n > 0 ? s->n : 1, sizeof(int));

    /*
      Initialize the mutable working response to kappa. The PG update will
      replace it with kappa / omega before any binomial Gaussian working-model
      update is used. Keeping this storage separate from yObserved avoids
      mutating the user's success counts.
    */
    for(i = 0; i < s->n; i++){
      s->nTrial[i] = INTEGER(trials_r)[i];
      if(s->nTrial[i] < 1)
        Rf_error("binomial trials must be positive");
      if(!R_FINITE(s->yObserved[i]) || s->yObserved[i] < 0.0 ||
         s->yObserved[i] > (double)s->nTrial[i])
        Rf_error("binomial successes must lie between 0 and trials");
      s->y[i] = s->yObserved[i] - 0.5 * (double)s->nTrial[i];
    }
  } else if(s->likelihoodFamily == LIKELIHOOD_NEGATIVE_BINOMIAL){
    SEXP size_r = getListElement(backend_r, "nb_size");
    if(size_r == R_NilValue)
      size_r = getListElement(getListElement(backend_r, "likelihood"), "size");
    s->nbSize = as_real_scalar_strict(size_r, "backend$nb_size");
    if(!(s->nbSize > 0.0))
      Rf_error("negative-binomial size must be positive");

    s->yObserved = REAL(y_r);
    s->y = (double*)R_alloc(s->n > 0 ? s->n : 1, sizeof(double));
    for(i = 0; i < s->n; i++){
      if(!R_FINITE(s->yObserved[i]) || s->yObserved[i] < 0.0)
        Rf_error("negative-binomial counts must be nonnegative");
      s->y[i] = 0.5 * (s->yObserved[i] - s->nbSize) + std::log(s->nbSize);
    }
  } else {
    s->y = REAL(y_r);
    s->yObserved = s->y;
  }

  X_r = getListElement(backend_r, "X_obs");
  if(X_r == R_NilValue)
    X_r = getListElement(backend_r, "X");
  require_matrix_real(X_r, s->n, s->p, "backend$X_obs");
  s->X = REAL(X_r);

  {
    SEXP beta_prior_type_r = getListElement(backend_r, "beta_prior_type");
    SEXP beta_prior_mean_r = getListElement(backend_r, "beta_prior_mean");
    SEXP beta_prior_precision_r = getListElement(backend_r, "beta_prior_precision");

    s->betaPriorType = is_pg_likelihood(s) ? 1 : 0;
    s->betaPriorMean = (double*)R_alloc(s->p > 0 ? s->p : 1, sizeof(double));
    s->betaPriorPrecision = (double*)R_alloc(s->p > 0 ? s->p : 1, sizeof(double));
    for(i = 0; i < s->p; i++){
      s->betaPriorMean[i] = 0.0;
      s->betaPriorPrecision[i] =
        is_pg_likelihood(s) ? 0.01 : 0.0;
    }

    if(beta_prior_type_r != R_NilValue)
      s->betaPriorType = as_int_scalar(beta_prior_type_r, "backend$beta_prior_type");
    if(s->betaPriorType < 0 || s->betaPriorType > 1)
      Rf_error("backend$beta_prior_type must be 0 for flat or 1 for normal");

    if(beta_prior_mean_r != R_NilValue){
      if(TYPEOF(beta_prior_mean_r) != REALSXP || LENGTH(beta_prior_mean_r) != s->p)
        Rf_error("backend$beta_prior_mean malformed");
      for(i = 0; i < s->p; i++){
        s->betaPriorMean[i] = REAL(beta_prior_mean_r)[i];
        if(!R_FINITE(s->betaPriorMean[i]))
          Rf_error("backend$beta_prior_mean entries must be finite");
      }
    }
    if(beta_prior_precision_r != R_NilValue){
      if(TYPEOF(beta_prior_precision_r) != REALSXP || LENGTH(beta_prior_precision_r) != s->p)
        Rf_error("backend$beta_prior_precision malformed");
      for(i = 0; i < s->p; i++){
        s->betaPriorPrecision[i] = REAL(beta_prior_precision_r)[i];
        if(!R_FINITE(s->betaPriorPrecision[i]) || s->betaPriorPrecision[i] < 0.0)
          Rf_error("backend$beta_prior_precision entries must be finite and nonnegative");
      }
    }
    if(s->betaPriorType == 0){
      for(i = 0; i < s->p; i++)
        s->betaPriorPrecision[i] = 0.0;
    }
  }

  Z_r = getListElement(backend_r, "Z_obs");
  if(Z_r == R_NilValue)
    Z_r = getListElement(backend_r, "Z");
  s->Zp = NULL;
  s->Zi = NULL;
  s->Zx = NULL;
  s->Z_nnz = 0;
  if(Z_r != R_NilValue){
    SEXP zp_r, zi_r, zx_r, zdim_r;

    if(!Rf_inherits(Z_r, "dgCMatrix"))
      Rf_error("backend$Z_obs must be dgCMatrix");
    
    zp_r   = R_do_slot(Z_r, Rf_install("p"));
    zi_r   = R_do_slot(Z_r, Rf_install("i"));
    zx_r   = R_do_slot(Z_r, Rf_install("x"));
    zdim_r = R_do_slot(Z_r, Rf_install("Dim"));

    if(TYPEOF(zp_r) != INTSXP || LENGTH(zp_r) != s->q + 1)
      Rf_error("backend$Z@p malformed");
    if(TYPEOF(zi_r) != INTSXP)
      Rf_error("backend$Z@i malformed");
    if(TYPEOF(zx_r) != REALSXP)
      Rf_error("backend$Z@x malformed");
    if(TYPEOF(zdim_r) != INTSXP || LENGTH(zdim_r) != 2)
      Rf_error("backend$Z@Dim malformed");
    if(INTEGER(zdim_r)[0] != s->n || INTEGER(zdim_r)[1] != s->q)
      Rf_error("backend$Z_obs has wrong dimensions");

    s->Zp = INTEGER(zp_r);
    s->Zi = INTEGER(zi_r);
    s->Zx = REAL(zx_r);
    s->Z_nnz = LENGTH(zi_r);
    if(LENGTH(zx_r) != s->Z_nnz)
      Rf_error("backend$Z slots i/x length mismatch");
  } else if(s->q > 0){
    Rf_error("backend$q > 0 but backend$Z is missing");
  }

  /*
    tau^2 starting value, IG prior parameters, and MH tuning.
    All have safe defaults so the backend contract can omit them.
  */
  {
    SEXP tau_start_r = getListElement(backend_r, "tau_sq_starting");
    SEXP tau_ig_r    = getListElement(backend_r, "tau_sq_IG");
    SEXP tau_prior_r = getListElement(backend_r, "tau_sq_prior");
    SEXP tau_tune_r  = getListElement(backend_r, "tau_sq_tuning");

    if(tau_start_r == R_NilValue || LENGTH(tau_start_r) != 1)
      s->tauSq = 1.0;
    else if(TYPEOF(tau_start_r) == REALSXP)
      s->tauSq = REAL(tau_start_r)[0];
    else if(TYPEOF(tau_start_r) == INTSXP)
      s->tauSq = (double)INTEGER(tau_start_r)[0];
    else
      Rf_error("backend$tau_sq_starting malformed");

    if(s->tauSq <= 0.0)
      Rf_error("backend$tau_sq_starting must be positive");

    if(tau_ig_r == R_NilValue || LENGTH(tau_ig_r) != 2){
      s->tauSqShape = 2.0;
      s->tauSqScale = 1.0;
    } else if(TYPEOF(tau_ig_r) == REALSXP){
      s->tauSqShape = REAL(tau_ig_r)[0];
      s->tauSqScale  = REAL(tau_ig_r)[1];
    } else if(TYPEOF(tau_ig_r) == INTSXP){
      s->tauSqShape = (double)INTEGER(tau_ig_r)[0];
      s->tauSqScale  = (double)INTEGER(tau_ig_r)[1];
    } else {
      Rf_error("backend$tau_sq_IG malformed");
    }

    if(s->tauSqShape <= 0.0 || s->tauSqScale <= 0.0)
      Rf_error("backend$tau_sq_IG parameters must be positive");
    s->tauSqPrior = parse_prior_spec(tau_prior_r, "backend$tau_sq_prior",
                                     PRIOR_IG, s->tauSqShape, s->tauSqScale);

    if(tau_tune_r == R_NilValue || LENGTH(tau_tune_r) != 1)
      s->tauTune = 0.2;
    else if(TYPEOF(tau_tune_r) == REALSXP)
      s->tauTune = REAL(tau_tune_r)[0];
    else if(TYPEOF(tau_tune_r) == INTSXP)
      s->tauTune = (double)INTEGER(tau_tune_r)[0];
    else
      Rf_error("backend$tau_sq_tuning malformed");

    if(s->tauTune < 0.0)
      Rf_error("backend$tau_sq_tuning must be nonnegative");
    if(is_pg_likelihood(s))
      s->tauTune = 0.0;
  }

  s->residualModel = 0;
  s->obsPrecisionFixed = NULL;
  s->nResidualGroup = 0;
  s->residualGroupIndex = NULL;
  s->residualVariance = NULL;
  s->residualVarianceShape = NULL;
  s->residualVarianceScale = NULL;
  s->residualVariancePrior = NULL;
  s->residualVarianceTune = NULL;
  s->residualVarianceMeanLog = NULL;
  s->residualVarianceSdLog = NULL;
  s->residualScaledVhat = NULL;
  s->residualScaledWeight = NULL;
  {
    SEXP residual_r = getListElement(backend_r, "residual_model");
    SEXP type_r, prec_r, group_r, start_r, shape_r, scale_r, tune_r, meanlog_r, sdlog_r, vhat_r, weight_r;
    const char *typeName;
    int i;

    if(residual_r != R_NilValue){
      type_r = getListElement(residual_r, "type");
      if(type_r != R_NilValue){
        typeName = as_char_scalar(type_r, "backend$residual_model$type");
        if(std::strcmp(typeName, "global_tau") == 0){
          s->residualModel = 0;
        } else if(std::strcmp(typeName, "fixed_variance") == 0){
          prec_r = getListElement(residual_r, "obs_precision_obs");
          if(TYPEOF(prec_r) != REALSXP || LENGTH(prec_r) != s->n)
            Rf_error("backend$residual_model$obs_precision_obs malformed");
          s->residualModel = 1;
          s->obsPrecisionFixed = REAL(prec_r);
          s->tauTune = 0.0;
          for(i = 0; i < s->n; i++){
            if(!R_FINITE(s->obsPrecisionFixed[i]) || !(s->obsPrecisionFixed[i] > 0.0))
              Rf_error("backend fixed residual precision must be finite and positive");
          }
        } else if(std::strcmp(typeName, "group_ig_variance") == 0){
          group_r = getListElement(residual_r, "group_index_obs");
          start_r = getListElement(residual_r, "starting");
          shape_r = getListElement(residual_r, "shape");
          scale_r = getListElement(residual_r, "scale");
          tune_r = getListElement(residual_r, "tuning");
          SEXP prior_r = getListElement(residual_r, "prior");

          if(TYPEOF(group_r) != INTSXP || LENGTH(group_r) != s->n)
            Rf_error("backend$residual_model$group_index_obs malformed");
          if(TYPEOF(start_r) != REALSXP)
            Rf_error("backend$residual_model$starting malformed");
          s->nResidualGroup = LENGTH(start_r);
          if(TYPEOF(shape_r) != REALSXP || LENGTH(shape_r) != s->nResidualGroup ||
             TYPEOF(scale_r) != REALSXP || LENGTH(scale_r) != s->nResidualGroup ||
             TYPEOF(tune_r) != REALSXP || LENGTH(tune_r) != s->nResidualGroup)
            Rf_error("backend residual variance controls malformed");
          if(prior_r != R_NilValue){
            SEXP prior_dim_r = Rf_getAttrib(prior_r, R_DimSymbol);
            if(TYPEOF(prior_r) != REALSXP || prior_dim_r == R_NilValue ||
               TYPEOF(prior_dim_r) != INTSXP || LENGTH(prior_dim_r) != 2 ||
               INTEGER(prior_dim_r)[0] != s->nResidualGroup ||
               INTEGER(prior_dim_r)[1] != 6)
              Rf_error("backend residual variance prior malformed");
          }

          s->residualModel = 2;
          s->tauTune = 0.0;
          s->residualGroupIndex = (int*)R_alloc(s->n > 0 ? s->n : 1, sizeof(int));
          s->residualVariance = (double*)R_alloc(s->nResidualGroup > 0 ? s->nResidualGroup : 1, sizeof(double));
          s->residualVarianceShape = (double*)R_alloc(s->nResidualGroup > 0 ? s->nResidualGroup : 1, sizeof(double));
          s->residualVarianceScale = (double*)R_alloc(s->nResidualGroup > 0 ? s->nResidualGroup : 1, sizeof(double));
          s->residualVariancePrior = (PriorSpec*)R_alloc(s->nResidualGroup > 0 ? s->nResidualGroup : 1, sizeof(PriorSpec));
          s->residualVarianceTune = (double*)R_alloc(s->nResidualGroup > 0 ? s->nResidualGroup : 1, sizeof(double));

          for(i = 0; i < s->n; i++){
            s->residualGroupIndex[i] = INTEGER(group_r)[i] - 1;
            if(s->residualGroupIndex[i] < 0 || s->residualGroupIndex[i] >= s->nResidualGroup)
              Rf_error("backend residual group index out of bounds");
          }
          for(i = 0; i < s->nResidualGroup; i++){
            s->residualVariance[i] = REAL(start_r)[i];
            s->residualVarianceShape[i] = REAL(shape_r)[i];
            s->residualVarianceScale[i] = REAL(scale_r)[i];
            if(prior_r == R_NilValue){
              s->residualVariancePrior[i].family = PRIOR_IG;
              s->residualVariancePrior[i].p1 = s->residualVarianceShape[i];
              s->residualVariancePrior[i].p2 = s->residualVarianceScale[i];
            } else {
              s->residualVariancePrior[i].family = (int)REAL(prior_r)[i];
              s->residualVariancePrior[i].p1 = REAL(prior_r)[i + s->nResidualGroup];
              s->residualVariancePrior[i].p2 = REAL(prior_r)[i + 2 * s->nResidualGroup];
            }
            s->residualVarianceTune[i] = REAL(tune_r)[i];
            if(!R_FINITE(s->residualVariance[i]) || !(s->residualVariance[i] > 0.0))
              Rf_error("backend residual starting variance must be finite and positive");
            if(s->residualVarianceTune[i] < 0.0)
              Rf_error("backend residual variance tuning must be nonnegative");
          }
        } else if(std::strcmp(typeName, "scaled_variance") == 0){
          start_r = getListElement(residual_r, "starting");
          tune_r = getListElement(residual_r, "tuning");
          meanlog_r = getListElement(residual_r, "prior_meanlog");
          sdlog_r = getListElement(residual_r, "prior_sdlog");
          vhat_r = getListElement(residual_r, "vhat_obs");
          weight_r = getListElement(residual_r, "weight_obs");

          if(TYPEOF(start_r) != REALSXP)
            Rf_error("backend$residual_model$starting malformed");
          s->nResidualGroup = LENGTH(start_r);
          if(!(s->nResidualGroup == 1 || s->nResidualGroup == 2))
            Rf_error("backend scaled residual model must have one or two parameters");
          if(TYPEOF(tune_r) != REALSXP || LENGTH(tune_r) != s->nResidualGroup ||
             TYPEOF(meanlog_r) != REALSXP || LENGTH(meanlog_r) != s->nResidualGroup ||
             TYPEOF(sdlog_r) != REALSXP || LENGTH(sdlog_r) != s->nResidualGroup ||
             TYPEOF(vhat_r) != REALSXP || LENGTH(vhat_r) != s->n ||
             TYPEOF(weight_r) != REALSXP || LENGTH(weight_r) != s->n)
            Rf_error("backend scaled residual controls malformed");

          s->residualModel = 3;
          s->tauTune = 0.0;
          s->residualVariance = (double*)R_alloc(s->nResidualGroup, sizeof(double));
          s->residualVarianceTune = (double*)R_alloc(s->nResidualGroup, sizeof(double));
          s->residualVarianceMeanLog = (double*)R_alloc(s->nResidualGroup, sizeof(double));
          s->residualVarianceSdLog = (double*)R_alloc(s->nResidualGroup, sizeof(double));
          s->residualScaledVhat = (double*)R_alloc(s->n > 0 ? s->n : 1, sizeof(double));
          s->residualScaledWeight = (double*)R_alloc(s->n > 0 ? s->n : 1, sizeof(double));

          for(i = 0; i < s->nResidualGroup; i++){
            s->residualVariance[i] = REAL(start_r)[i];
            s->residualVarianceTune[i] = REAL(tune_r)[i];
            s->residualVarianceMeanLog[i] = REAL(meanlog_r)[i];
            s->residualVarianceSdLog[i] = REAL(sdlog_r)[i];
            if(!R_FINITE(s->residualVariance[i]) || !(s->residualVariance[i] > 0.0))
              Rf_error("backend scaled residual starting values must be finite and positive");
            if(!R_FINITE(s->residualVarianceTune[i]) || s->residualVarianceTune[i] < 0.0)
              Rf_error("backend scaled residual tuning must be nonnegative");
            if(!R_FINITE(s->residualVarianceMeanLog[i]) ||
               !R_FINITE(s->residualVarianceSdLog[i]) || !(s->residualVarianceSdLog[i] > 0.0))
              Rf_error("backend scaled residual log-prior parameters are invalid");
          }
          for(i = 0; i < s->n; i++){
            s->residualScaledVhat[i] = REAL(vhat_r)[i];
            s->residualScaledWeight[i] = REAL(weight_r)[i];
            if(!R_FINITE(s->residualScaledVhat[i]) || !(s->residualScaledVhat[i] > 0.0))
              Rf_error("backend scaled residual vhat must be finite and positive");
            if(!R_FINITE(s->residualScaledWeight[i]) ||
               s->residualScaledWeight[i] < 0.0 || s->residualScaledWeight[i] > 1.0)
              Rf_error("backend scaled residual weights must be in [0, 1]");
          }
        } else {
          Rf_error("unsupported backend$residual_model$type");
        }
      }
    }
  }

  s->obsPrecision = (double*)R_alloc(s->n > 0 ? s->n : 1, sizeof(double));
  if(s->likelihoodFamily == LIKELIHOOD_BINOMIAL){
    for(i = 0; i < s->n; i++)
      s->obsPrecision[i] = 0.25 * (double)s->nTrial[i];
  } else if(s->likelihoodFamily == LIKELIHOOD_NEGATIVE_BINOMIAL){
    for(i = 0; i < s->n; i++)
      s->obsPrecision[i] = 0.25 * (s->yObserved[i] + s->nbSize);
  } else {
    refresh_observation_precision(s);
  }

  /*
    Random-effects block metadata, starting variances, and IG priors for alpha.
    Each alpha column j belongs to block reBlockID[j], and the prior
    precision contribution is 1 / sigmaSqRE[reBlockID[j]].
  */
  s->nRE = 0;
  s->maxBlockSize = 0;
  s->reBlockID = (int*)R_alloc(s->q > 0 ? s->q : 1, sizeof(int));
  s->reBlockSize = (int*)R_alloc(1, sizeof(int));
  s->sigmaSqRE = (double*)R_alloc(1, sizeof(double));
  s->sigmaSqREShape = (double*)R_alloc(1, sizeof(double));
  s->sigmaSqREScale = (double*)R_alloc(1, sizeof(double));
  s->alphaBlockBuf = (double*)R_alloc(1, sizeof(double));

  {
    SEXP re_r = getListElement(backend_r, "re");
    SEXP sigma_start_r = getListElement(backend_r, "sigma_sq_re_starting");
    SEXP sigma_ig_r = getListElement(backend_r, "sigma_sq_re_IG");

    if(s->q > 0){
      SEXP re_q_r, block_id_r, sigma_ig_dim_r;

      if(re_r == R_NilValue || TYPEOF(re_r) != VECSXP)
        Rf_error("backend$re must be a list when backend$q > 0");

      re_q_r = getListElement(re_r, "q");
      block_id_r = getListElement(re_r, "block_id");

      if(block_id_r == R_NilValue || TYPEOF(block_id_r) != INTSXP || LENGTH(block_id_r) != s->q)
        Rf_error("backend$re$block_id malformed");

      if(re_q_r == R_NilValue)
        Rf_error("backend$re$q missing");

      s->nRE = LENGTH(re_q_r);
      if(s->nRE < 1)
        Rf_error("backend$re$q must have positive length when backend$q > 0");

      s->reBlockSize = (int*)R_alloc(s->nRE, sizeof(int));
      s->sigmaSqRE = (double*)R_alloc(s->nRE, sizeof(double));
      s->sigmaSqREShape = (double*)R_alloc(s->nRE, sizeof(double));
      s->sigmaSqREScale = (double*)R_alloc(s->nRE, sizeof(double));

      for(i = 0; i < s->nRE; i++)
        s->reBlockSize[i] = 0;

      if(sigma_start_r == R_NilValue)
        Rf_error("backend$sigma_sq_re_starting missing when backend$q > 0");

      if(TYPEOF(sigma_start_r) == REALSXP){
        if(LENGTH(sigma_start_r) != s->nRE)
          Rf_error("backend$sigma_sq_re_starting has wrong length");
        for(i = 0; i < s->nRE; i++)
          s->sigmaSqRE[i] = REAL(sigma_start_r)[i];
      } else if(TYPEOF(sigma_start_r) == INTSXP){
        if(LENGTH(sigma_start_r) != s->nRE)
          Rf_error("backend$sigma_sq_re_starting has wrong length");
        for(i = 0; i < s->nRE; i++)
          s->sigmaSqRE[i] = (double)INTEGER(sigma_start_r)[i];
      } else {
        Rf_error("backend$sigma_sq_re_starting malformed");
      }

      for(i = 0; i < s->nRE; i++){
        if(s->sigmaSqRE[i] <= 0.0)
          Rf_error("backend$sigma_sq_re_starting entries must be positive");
      }

      if(sigma_ig_r == R_NilValue || TYPEOF(sigma_ig_r) != REALSXP)
        Rf_error("backend$sigma_sq_re_IG malformed");

      sigma_ig_dim_r = Rf_getAttrib(sigma_ig_r, R_DimSymbol);
      if(sigma_ig_dim_r == R_NilValue || TYPEOF(sigma_ig_dim_r) != INTSXP || LENGTH(sigma_ig_dim_r) != 2)
        Rf_error("backend$sigma_sq_re_IG must be a matrix");
      if(INTEGER(sigma_ig_dim_r)[0] != s->nRE || INTEGER(sigma_ig_dim_r)[1] != 2)
        Rf_error("backend$sigma_sq_re_IG has wrong dimensions");

      for(i = 0; i < s->nRE; i++){
        s->sigmaSqREShape[i] = REAL(sigma_ig_r)[i];
        s->sigmaSqREScale[i]  = REAL(sigma_ig_r)[i + s->nRE];
        if(s->sigmaSqREShape[i] <= 0.0 || s->sigmaSqREScale[i] <= 0.0)
          Rf_error("backend$sigma_sq_re_IG entries must be positive");
      }

      for(i = 0; i < s->q; i++){
        s->reBlockID[i] = INTEGER(block_id_r)[i] - 1;
        if(s->reBlockID[i] < 0 || s->reBlockID[i] >= s->nRE)
          Rf_error("backend$re$block_id contains invalid block index");
        s->reBlockSize[s->reBlockID[i]]++;
      }

      for(i = 0; i < s->nRE; i++){
        if(s->reBlockSize[i] > s->maxBlockSize)
          s->maxBlockSize = s->reBlockSize[i];
      }

      s->alphaBlockBuf = (double*)R_alloc(s->maxBlockSize > 0 ? s->maxBlockSize : 1, sizeof(double));
    }
  }

  graphs_r = getListElement(backend_r, "graphs");
  terms_r = getListElement(backend_r, "process_terms_obs");
  if(terms_r == R_NilValue)
    terms_r = getListElement(backend_r, "process_terms");

  if(TYPEOF(graphs_r) != VECSXP)
    Rf_error("backend$graphs must be a list");
  if(TYPEOF(terms_r) != VECSXP)
    Rf_error("backend$process_terms must be a list");

  s->nGraphs = LENGTH(graphs_r);
  s->nTerms = LENGTH(terms_r);

  s->graphs = (GraphState*)R_alloc(s->nGraphs, sizeof(GraphState));
  s->terms = (TermState*)R_alloc(s->nTerms, sizeof(TermState));
  s->nThetaTotal = 0;
  s->covParamAccept = 0;
  s->covParamAttempts = 0;

  M_R_cholmod_start(&s->cm);
  configure_cholmod_control(&s->cm, backend_r);
  s->cm.supernodal = CHOLMOD_SUPERNODAL;
  s->cm.final_ll = 1;
  s->cm.final_super = 1;

  for(i = 0; i < s->nGraphs; i++)
    init_graph_state_from_backend(s->graphs + i, VECTOR_ELT(graphs_r, i));
  
  for(i = 0; i < s->nTerms; i++){
    init_term_state_from_backend(s->terms + i, s->graphs, s->nGraphs, VECTOR_ELT(terms_r, i), s->n);
    s->terms[i].wOffset = s->qLatTotal;
    s->qLatTotal += s->terms[i].qLat;
    s->nThetaTotal += s->terms[i].thetaDim;
  }

  s->scratch_BF_m = 0;
  for(i = 0; i < s->nGraphs; i++){
    int j;
    if(s->graphs[i].type == GRAPH_NNGP){
      for(j = 0; j < s->graphs[i].nNode; j++){
        if(s->graphs[i].nnCount[j] > s->scratch_BF_m)
          s->scratch_BF_m = s->graphs[i].nnCount[j];
      }
    }
  }
  s->scratch_BF_bk_n = 1;
  for(i = 0; i < s->nTerms; i++){
    GraphState *g = s->graphs + s->terms[i].graphIndex;
    const CorModelInfo *info;

    if(g->type != GRAPH_NNGP || s->terms[i].covModelIndex < 0)
      continue;

    info = get_cor_model_info(s->terms[i].covModelIndex);
    if(std::strcmp(info->name, "matern") == 0){
      double nuUpper = s->terms[i].thetaUpper[1];
      double nuWorkMax = nuUpper;
      int nb;

      if(!R_FINITE(nuWorkMax) || s->terms[i].theta[1] > nuWorkMax)
        nuWorkMax = s->terms[i].theta[1];

      if(!R_FINITE(nuWorkMax) || nuWorkMax <= 0.0)
        Rf_error("matern nngp term requires finite positive nu values");

      nb = 1 + (int)std::floor(nuWorkMax);
      if(nb > s->scratch_BF_bk_n)
        s->scratch_BF_bk_n = nb;
    }
  }
  s->scratch_BF_c = (double*)R_alloc((s->scratch_BF_m > 0 ? s->scratch_BF_m : 1) * s->nOmpThreads, sizeof(double));
  s->scratch_BF_C = (double*)R_alloc((s->scratch_BF_m > 0 ? s->scratch_BF_m * s->scratch_BF_m : 1) * s->nOmpThreads, sizeof(double));
  s->scratch_BF_bk = (double*)R_alloc(s->scratch_BF_bk_n * s->nOmpThreads, sizeof(double));

  s->r = (double*)R_alloc(s->n, sizeof(double));
  s->Za = (double*)R_alloc(s->n, sizeof(double));
  s->alpha = (double*)R_alloc(s->q > 0 ? s->q : 1, sizeof(double));
  s->nWork1 = (double*)R_alloc(s->n, sizeof(double));
  s->nWork2 = (double*)R_alloc(s->n, sizeof(double));
  s->latentSolveRhs = (double*)R_alloc(s->qLatTotal > 0 ? s->qLatTotal : 1, sizeof(double));
  s->latentNoiseZ = (double*)R_alloc(s->qLatTotal > 0 ? s->qLatTotal : 1, sizeof(double));

  for(i = 0; i < s->n; i++){
    s->Za[i] = 0.0;
    s->nWork1[i] = 0.0;
    s->nWork2[i] = 0.0;
  }
  for(i = 0; i < (s->q > 0 ? s->q : 1); i++)
    s->alpha[i] = 0.0;
}

void free_cov_proposal_state(SamplerState *s)
{
  int t;

  if(!s->covProp.initialized)
    return;

  for(t = 0; t < s->nTerms; t++){
    if(s->covProp.Q != NULL && s->covProp.Q[t] != NULL)
      R_Free(s->covProp.Q[t]);
  }

  if(s->covProp.M_lat_fac != NULL)
    M_cholmod_free_factor(&s->covProp.M_lat_fac, &s->cm);
  if(s->covProp.M_lat != NULL)
    M_cholmod_free_sparse(&s->covProp.M_lat, &s->cm);

  s->covProp.initialized = 0;
}

void free_sampler_state(SamplerState *s)
{
  int t;
  
  free_cov_proposal_state(s);

  for(t = 0; t < s->nTerms; t++){
    if(s->terms[t].Q != NULL)
      R_Free(s->terms[t].Q);
    if(s->terms[t].carQDiagCacheIdx != NULL)
      R_Free(s->terms[t].carQDiagCacheIdx);
    if(s->terms[t].carQOffCacheIdx != NULL)
      R_Free(s->terms[t].carQOffCacheIdx);
    if(s->terms[t].carQFac != NULL)
      M_cholmod_free_factor(&s->terms[t].carQFac, &s->cm);
    if(s->terms[t].carQ != NULL)
      M_cholmod_free_sparse(&s->terms[t].carQ, &s->cm);
  }
  
  if(s->M_lat_fac != NULL)
    M_cholmod_free_factor(&s->M_lat_fac, &s->cm);
  if(s->M_lat_sym != NULL)
    M_cholmod_free_factor(&s->M_lat_sym, &s->cm);
  if(s->M_lat != NULL)
    M_cholmod_free_sparse(&s->M_lat, &s->cm);
  M_cholmod_finish(&s->cm);
}

static void init_cov_proposal_block(CovProposalState *p,
                                    CovProposalBlock *b,
                                    int dim,
                                    int batchLength,
                                    int maxBatchHistory,
                                    double targetAccept,
                                    int scalarMode,
                                    int warmupEnabled,
                                    int warmupBatchLength,
                                    int warmupMinBatches,
                                    int warmupMaxBatches,
                                    double warmupTargetLower,
                                    double warmupTargetUpper,
                                    double warmupNearZero)
{
  int j;

  std::memset(b, 0, sizeof(CovProposalBlock));
  b->dim = dim;
  b->batchLength = batchLength;
  b->batchPos = 0;
  b->batchIndex = 0;
  b->batchAccept = 0;
  b->acceptCount = 0;
  b->maxBatchHistory = maxBatchHistory;
  b->scalarMode = scalarMode;
  b->targetAccept = targetAccept;
  b->warmupEnabled = warmupEnabled;
  b->warmupBatchLength = warmupBatchLength;
  b->warmupMinBatches = warmupMinBatches;
  b->warmupMaxBatches = warmupMaxBatches;
  b->warmupTargetLower = warmupTargetLower;
  b->warmupTargetUpper = warmupTargetUpper;
  b->warmupNearZero = warmupNearZero;
  b->warmupBatches = 0;
  b->warmupNAttempted = 0;
  b->warmupNAccepted = 0;
  b->warmupStoppedReason = 0;
  b->c0 = 1.0;
  b->c1 = 0.8;
  b->sigmaSqM = dim > 0 ? 2.4 * 2.4 / (double)dim : NA_REAL;
  if(scalarMode)
    b->sigmaSqM = 1.0;
  b->baseSigmaSqM = b->sigmaSqM;
  b->lastBatchAccept = NA_REAL;

  b->paramIndex = (int*)R_alloc(dim > 0 ? dim : 1, sizeof(int));
  b->eta = (double*)R_alloc(dim > 0 ? dim : 1, sizeof(double));
  b->etaProp = (double*)R_alloc(dim > 0 ? dim : 1, sizeof(double));
  b->z = (double*)R_alloc(dim > 0 ? dim : 1, sizeof(double));
  b->Sigma0 = (double*)R_alloc(dim > 0 ? dim * dim : 1, sizeof(double));
  b->proposalCov = (double*)R_alloc(dim > 0 ? dim * dim : 1, sizeof(double));
  b->proposalChol = (double*)R_alloc(dim > 0 ? dim * dim : 1, sizeof(double));
  b->warmupStartingEta = (double*)R_alloc(dim > 0 ? dim : 1, sizeof(double));
  b->warmupEndingEta = (double*)R_alloc(dim > 0 ? dim : 1, sizeof(double));
  b->warmupStartingProposalCov = (double*)R_alloc(dim > 0 ? dim * dim : 1, sizeof(double));
  b->warmupEndingProposalCov = (double*)R_alloc(dim > 0 ? dim * dim : 1, sizeof(double));
  b->warmupBatchAcceptHistory = (double*)R_alloc(warmupMaxBatches > 0 ? warmupMaxBatches : 1, sizeof(double));
  b->warmupProposalScaleHistory = (double*)R_alloc(warmupMaxBatches > 0 ? warmupMaxBatches : 1, sizeof(double));
  b->batchSamples = (double*)R_alloc(dim > 0 ? batchLength * dim : 1, sizeof(double));
  b->batchAcceptHistory = (double*)R_alloc(maxBatchHistory > 0 ? maxBatchHistory : 1, sizeof(double));
  b->proposalScaleHistory = (double*)R_alloc(maxBatchHistory > 0 ? maxBatchHistory : 1, sizeof(double));

  for(j = 0; j < (maxBatchHistory > 0 ? maxBatchHistory : 1); j++){
    b->batchAcceptHistory[j] = NA_REAL;
    b->proposalScaleHistory[j] = NA_REAL;
  }
  for(j = 0; j < (warmupMaxBatches > 0 ? warmupMaxBatches : 1); j++){
    b->warmupBatchAcceptHistory[j] = NA_REAL;
    b->warmupProposalScaleHistory[j] = NA_REAL;
  }
  for(j = 0; j < dim * dim; j++){
    b->Sigma0[j] = 0.0;
    b->proposalCov[j] = 0.0;
    b->proposalChol[j] = 0.0;
    b->warmupStartingProposalCov[j] = 0.0;
    b->warmupEndingProposalCov[j] = 0.0;
  }
  for(j = 0; j < (dim > 0 ? dim : 1); j++){
    b->paramIndex[j] = -1;
    b->warmupStartingEta[j] = NA_REAL;
    b->warmupEndingEta[j] = NA_REAL;
  }

  (void)p;
}

static int raw_metropolis_block_id(SamplerState *s, CovProposalState *p, int param)
{
  int type, term;

  type = p->paramType[param];
  term = p->paramTerm[param];

  if(s->metropolisBlocking == 5)
    return param;

  if(s->metropolisBlocking == 0)
    return 0;

  if(s->metropolisBlocking == 1){
    if(type == 0 || type == 3)
      return 0;
    return term + 1;
  }

  if(s->metropolisBlocking == 2){
    if(type == 0 || type == 3)
      return 0;
    return 1;
  }

  if(s->metropolisBlocking == 3){
    if(type == 2)
      return 1;
    return 0;
  }

  if(s->metropolisBlocking == 4){
    if(type == 0 || type == 3)
      return 0;
    if(type == 1)
      return 1;
    return 2;
  }

  return 0;
}

static void assign_cov_proposal_blocks(SamplerState *s, CovProposalState *p)
{
  int i, j, b, nRaw, nBlock, raw, dim, maxBatchHistory;
  int *rawToBlock, *blockDim, *blockPos;

  if(p->dim <= 0){
    p->nBlocks = 0;
    p->blocks = (CovProposalBlock*)R_alloc(1, sizeof(CovProposalBlock));
    return;
  }

  nRaw = s->metropolisBlocking == 5 ? p->dim : s->nTerms + 3;
  rawToBlock = (int*)R_alloc(nRaw, sizeof(int));
  blockDim = (int*)R_alloc(nRaw, sizeof(int));
  blockPos = (int*)R_alloc(nRaw, sizeof(int));
  for(i = 0; i < nRaw; i++){
    rawToBlock[i] = -1;
    blockDim[i] = 0;
    blockPos[i] = 0;
  }

  for(i = 0; i < p->dim; i++){
    raw = raw_metropolis_block_id(s, p, i);
    if(raw < 0 || raw >= nRaw)
      Rf_error("internal error: invalid Metropolis block id");
    if(rawToBlock[raw] < 0)
      rawToBlock[raw] = 0;
    blockDim[raw]++;
  }

  nBlock = 0;
  for(i = 0; i < nRaw; i++){
    if(rawToBlock[i] >= 0)
      rawToBlock[i] = nBlock++;
  }

  p->nBlocks = nBlock;
  p->blocks = (CovProposalBlock*)R_alloc(nBlock > 0 ? nBlock : 1, sizeof(CovProposalBlock));
  maxBatchHistory = s->nSamples > 0 ? (s->nSamples / s->metropolisBatchLength) : 0;

  for(i = 0; i < nRaw; i++){
    if(rawToBlock[i] >= 0){
      dim = blockDim[i];
      init_cov_proposal_block(p, p->blocks + rawToBlock[i],
                              dim, s->metropolisBatchLength, maxBatchHistory,
                              s->metropolisTargetAccept,
                              s->metropolisBlocking == 5,
                              s->warmupEnabled, s->warmupBatchLength,
                              s->warmupMinBatches, s->warmupMaxBatches, s->warmupTargetLower,
                              s->warmupTargetUpper, s->warmupNearZero);
    }
  }

  for(i = 0; i < p->dim; i++){
    raw = raw_metropolis_block_id(s, p, i);
    b = rawToBlock[raw];
    p->paramBlock[i] = b;
    j = blockPos[raw]++;
    p->blocks[b].paramIndex[j] = i;
  }

  for(b = 0; b < p->nBlocks; b++){
    CovProposalBlock *blk = p->blocks + b;
    for(i = 0; i < blk->dim; i++){
      int gi = blk->paramIndex[i];
      for(j = 0; j < blk->dim; j++){
        int gj = blk->paramIndex[j];
        blk->Sigma0[i + blk->dim * j] = p->Sigma0[gi + p->dim * gj];
      }
    }
  }
}

void init_cov_proposal_state(SamplerState *s)
{
  int t, j, k, info;
  double scale;
  char lower = 'L';
  CovProposalState *p;

  p = &s->covProp;
  std::memset(p, 0, sizeof(CovProposalState));

  p->sigmaSq = (double*)R_alloc(s->nTerms > 0 ? s->nTerms : 1, sizeof(double));
  p->theta = (double**)R_alloc(s->nTerms > 0 ? s->nTerms : 1, sizeof(double*));
  p->B = (double**)R_alloc(s->nTerms > 0 ? s->nTerms : 1, sizeof(double*));
  p->F = (double**)R_alloc(s->nTerms > 0 ? s->nTerms : 1, sizeof(double*));
  p->Q = (double**)R_alloc(s->nTerms > 0 ? s->nTerms : 1, sizeof(double*));
  p->logDetQ = (double*)R_alloc(s->nTerms > 0 ? s->nTerms : 1, sizeof(double));
  p->residualVariance = (double*)R_alloc(s->nResidualGroup > 0 ? s->nResidualGroup : 1, sizeof(double));

  p->dim = 0;
  if(s->residualModel == 0 && s->tauTune > 0.0)
    p->dim++;
  if(s->residualModel == 2 || s->residualModel == 3){
    for(j = 0; j < s->nResidualGroup; j++)
      if(s->residualVarianceTune[j] > 0.0)
        p->dim++;
  }
  for(t = 0; t < s->nTerms; t++){
    TermState *term = s->terms + t;
    if(term->sigmaSqTune > 0.0)
      p->dim++;
    for(j = 0; j < term->thetaDim; j++)
      if(term->thetaTune[j] > 0.0)
        p->dim++;
  }
  p->paramType = (int*)R_alloc(p->dim > 0 ? p->dim : 1, sizeof(int));
  p->paramTerm = (int*)R_alloc(p->dim > 0 ? p->dim : 1, sizeof(int));
  p->paramTheta = (int*)R_alloc(p->dim > 0 ? p->dim : 1, sizeof(int));
  p->paramBlock = (int*)R_alloc(p->dim > 0 ? p->dim : 1, sizeof(int));
  p->nBlocks = 0;
  p->totalBlockAccept = 0;
  p->batchLength = s->metropolisBatchLength;
  p->batchPos = 0;
  p->batchIndex = 0;
  p->batchAccept = 0;
  p->maxBatchHistory = s->nSamples > 0 ? (s->nSamples / p->batchLength) : 0;
  p->targetAccept = s->metropolisTargetAccept;
  p->warmupEnabled = s->warmupEnabled;
  p->warmupBatchLength = s->warmupBatchLength;
  p->warmupMinBatches = s->warmupMinBatches;
  p->warmupMaxBatches = s->warmupMaxBatches;
  p->warmupTargetLower = s->warmupTargetLower;
  p->warmupTargetUpper = s->warmupTargetUpper;
  p->warmupNearZero = s->warmupNearZero;
  p->warmupBatches = 0;
  p->warmupNAttempted = 0;
  p->warmupNAccepted = 0;
  p->warmupStoppedReason = 0;
  p->c0 = 1.0;
  p->c1 = 0.8;
  p->sigmaSqM = p->dim > 0 ? 2.4 * 2.4 / (double)p->dim : NA_REAL;
  if(s->metropolisBlocking == 5)
    p->sigmaSqM = 1.0;
  p->baseSigmaSqM = p->sigmaSqM;
  p->lastBatchAccept = NA_REAL;

  p->eta = (double*)R_alloc(p->dim > 0 ? p->dim : 1, sizeof(double));
  p->etaProp = (double*)R_alloc(p->dim > 0 ? p->dim : 1, sizeof(double));
  p->z = (double*)R_alloc(p->dim > 0 ? p->dim : 1, sizeof(double));
  p->Sigma0 = (double*)R_alloc(p->dim > 0 ? p->dim * p->dim : 1, sizeof(double));
  p->proposalCov = (double*)R_alloc(p->dim > 0 ? p->dim * p->dim : 1, sizeof(double));
  p->proposalChol = (double*)R_alloc(p->dim > 0 ? p->dim * p->dim : 1, sizeof(double));
  p->warmupStartingEta = (double*)R_alloc(p->dim > 0 ? p->dim : 1, sizeof(double));
  p->warmupEndingEta = (double*)R_alloc(p->dim > 0 ? p->dim : 1, sizeof(double));
  p->warmupStartingProposalCov = (double*)R_alloc(p->dim > 0 ? p->dim * p->dim : 1, sizeof(double));
  p->warmupEndingProposalCov = (double*)R_alloc(p->dim > 0 ? p->dim * p->dim : 1, sizeof(double));
  p->warmupBatchAcceptHistory = (double*)R_alloc(p->warmupMaxBatches > 0 ? p->warmupMaxBatches : 1, sizeof(double));
  p->warmupProposalScaleHistory = (double*)R_alloc(p->warmupMaxBatches > 0 ? p->warmupMaxBatches : 1, sizeof(double));
  p->batchSamples = (double*)R_alloc(p->dim > 0 ? p->batchLength * p->dim : 1, sizeof(double));
  p->batchAcceptHistory = (double*)R_alloc(p->maxBatchHistory > 0 ? p->maxBatchHistory : 1, sizeof(double));
  p->proposalScaleHistory = (double*)R_alloc(p->maxBatchHistory > 0 ? p->maxBatchHistory : 1, sizeof(double));

  for(j = 0; j < (p->maxBatchHistory > 0 ? p->maxBatchHistory : 1); j++){
    p->batchAcceptHistory[j] = NA_REAL;
    p->proposalScaleHistory[j] = NA_REAL;
  }
  for(j = 0; j < (p->warmupMaxBatches > 0 ? p->warmupMaxBatches : 1); j++){
    p->warmupBatchAcceptHistory[j] = NA_REAL;
    p->warmupProposalScaleHistory[j] = NA_REAL;
  }

  for(j = 0; j < p->dim * p->dim; j++){
    p->Sigma0[j] = 0.0;
    p->proposalCov[j] = 0.0;
    p->proposalChol[j] = 0.0;
    p->warmupStartingProposalCov[j] = 0.0;
    p->warmupEndingProposalCov[j] = 0.0;
  }
  for(j = 0; j < (p->dim > 0 ? p->dim : 1); j++){
    p->warmupStartingEta[j] = NA_REAL;
    p->warmupEndingEta[j] = NA_REAL;
  }

  k = 0;
  if(s->residualModel == 0 && s->tauTune > 0.0){
    p->Sigma0[k + p->dim * k] = s->tauTune * s->tauTune;
    p->paramType[k] = 0;
    p->paramTerm[k] = -1;
    p->paramTheta[k] = -1;
    k++;
  }

  if(s->residualModel == 2 || s->residualModel == 3){
    for(j = 0; j < s->nResidualGroup; j++){
      p->residualVariance[j] = s->residualVariance[j];
      if(s->residualVarianceTune[j] > 0.0){
        p->Sigma0[k + p->dim * k] = s->residualVarianceTune[j] * s->residualVarianceTune[j];
        p->paramType[k] = 3;
        p->paramTerm[k] = j;
        p->paramTheta[k] = -1;
        k++;
      }
    }
  }

  for(t = 0; t < s->nTerms; t++){
    TermState *term = s->terms + t;
    GraphState *g = s->graphs + term->graphIndex;

    p->sigmaSq[t] = term->sigmaSq;
    if(term->sigmaSqTune > 0.0){
      p->Sigma0[k + p->dim * k] = term->sigmaSqTune * term->sigmaSqTune;
      p->paramType[k] = 1;
      p->paramTerm[k] = t;
      p->paramTheta[k] = -1;
      k++;
    }

    p->theta[t] = (double*)R_alloc(term->thetaDim > 0 ? term->thetaDim : 1, sizeof(double));
    for(j = 0; j < term->thetaDim; j++){
      p->theta[t][j] = term->theta[j];
      if(term->thetaTune[j] > 0.0){
        p->Sigma0[k + p->dim * k] = term->thetaTune[j] * term->thetaTune[j];
        p->paramType[k] = 2;
        p->paramTerm[k] = t;
        p->paramTheta[k] = j;
        k++;
      }
    }

    p->B[t] = NULL;
    p->F[t] = NULL;
    p->Q[t] = NULL;
    p->logDetQ[t] = term->logDetQ;

    if(g->type == GRAPH_NNGP){
      p->B[t] = (double*)R_alloc(g->totalNnbr > 0 ? g->totalNnbr : 1, sizeof(double));
      p->F[t] = (double*)R_alloc(term->nNode > 0 ? term->nNode : 1, sizeof(double));
      for(j = 0; j < g->totalNnbr; j++) p->B[t][j] = term->B[j];
      for(j = 0; j < term->nNode; j++) p->F[t][j] = term->F[j];
    } else if(g->type == GRAPH_GP){
      p->Q[t] = R_Calloc((size_t)(term->nNode > 0 ? term->nNode : 1) *
                         (size_t)(term->nNode > 0 ? term->nNode : 1), double);
      if(term->Q != NULL){
        for(j = 0; j < term->nNode * term->nNode; j++)
          p->Q[t][j] = term->Q[j];
      } else {
        for(j = 0; j < term->nNode * term->nNode; j++)
          p->Q[t][j] = 0.0;
      }
    }
  }

  p->M_lat = M_cholmod_copy_sparse(s->M_lat, &s->cm);
  if(p->M_lat == NULL || s->cm.status != CHOLMOD_OK)
    Rf_error("init_cov_proposal_state: failed to copy M_lat");

  p->M_lat_fac = M_cholmod_copy_factor(s->M_lat_sym, &s->cm);
  if(p->M_lat_fac == NULL || s->cm.status != CHOLMOD_OK)
    Rf_error("init_cov_proposal_state: failed to copy symbolic M factor");

  p->tauSq = s->tauSq;

  if(p->dim > 0){
    assign_cov_proposal_blocks(s, p);

    for(t = 0; t < p->nBlocks; t++){
      CovProposalBlock *blk = p->blocks + t;
      scale = blk->sigmaSqM;
      if(blk->scalarMode)
        scale = scale * scale;
      for(j = 0; j < blk->dim * blk->dim; j++){
        blk->proposalCov[j] = scale * blk->Sigma0[j];
        blk->proposalChol[j] = blk->proposalCov[j];
      }

      F77_CALL(dpotrf)(&lower, &blk->dim, blk->proposalChol, &blk->dim, &info FCONE);
      if(info != 0)
        Rf_error("init_cov_proposal_state: adaptive proposal covariance is not positive definite");

      for(j = 0; j < blk->dim; j++)
        for(k = j + 1; k < blk->dim; k++)
          blk->proposalChol[j + blk->dim * k] = 0.0;
    }

    /* Preserve the historical top-level diagnostics for the default joint block. */
    scale = p->sigmaSqM;
    if(s->metropolisBlocking == 5)
      scale = scale * scale;
    for(j = 0; j < p->dim * p->dim; j++){
      p->proposalCov[j] = scale * p->Sigma0[j];
      p->proposalChol[j] = p->proposalCov[j];
    }

    F77_CALL(dpotrf)(&lower, &p->dim, p->proposalChol, &p->dim, &info FCONE);
    if(info != 0)
      Rf_error("init_cov_proposal_state: adaptive proposal covariance is not positive definite");

    for(j = 0; j < p->dim; j++)
      for(k = j + 1; k < p->dim; k++)
        p->proposalChol[j + p->dim * k] = 0.0;
  }

  p->initialized = 1;
}
