#ifndef STLMM_INTERNAL_H
#define STLMM_INTERNAL_H

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Memory.h>
#include <R_ext/Utils.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif

#include <Matrix.h>
#include <cholmod.h>
#include <BayesLogit.h>

#include <cmath>
#include <cstdio>
#include <cstring>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "corModels.h"
#include "utils.h"

enum GraphType {
  GRAPH_NNGP,
  GRAPH_GP,
  GRAPH_AR1,
  GRAPH_CAR,
  GRAPH_CAR_TIME,
  GRAPH_DAGAR,
  GRAPH_DAGAR_TIME
};

enum TimeModel {
  TIME_MODEL_AR1,
  TIME_MODEL_EXP
};

enum CarModel {
  CAR_MODEL_PROPER,
  CAR_MODEL_LEROUX
};

enum TermType {
  TERM_INTERCEPT,
  TERM_SVC
};

enum PriorFamily {
  PRIOR_IG,
  PRIOR_UNIFORM,
  PRIOR_LOG_NORMAL,
  PRIOR_GAMMA,
  PRIOR_HALF_NORMAL,
  PRIOR_HALF_T,
  PRIOR_BETA
};

enum LikelihoodFamily {
  LIKELIHOOD_GAUSSIAN,
  LIKELIHOOD_BINOMIAL,
  LIKELIHOOD_NEGATIVE_BINOMIAL
};

typedef struct {
  int family;
  double p1;
  double p2;
} PriorSpec;

typedef struct {
  int type;
  int nNode;
  int nSpace;
  int nTime;
  int timeModel;
  int carModel;
  double *timeDelta;

  int dim;
  double *coords;

  int *nnIndx;
  int *nnCount;
  int *nnStart;
  int totalNnbr;

  double *degree;
  int *edgeI;
  int *edgeJ;
  double *edgeW;
  int nEdge;

  int *parentIndx;
  int *parentCount;
  int *parentStart;
  int totalParent;

} GraphState;

typedef struct {
  const char *name;
  const char *label;
  int type;
  int graphIndex;
  int nObs;
  int nNode;
  int qLat;
  int wOffset;

  int covModelIndex;
  corFunPtr corFun;
  int thetaDim;
  int distanceMode;
  double *theta;
  double *thetaTune;
  double *thetaLower;
  double *thetaUpper;
  double *B;
  double *F;
  double *Q;
  double logDetQ;

  int *map;
  int *obsIndx;
  int *obsStart;
  int *obsEnd;
  int *nodeNobs;

  double *x;
  double *scale;

  double sigmaSq;
  double sigmaSqTune;
  double sigmaSqShape;
  double sigmaSqScale;
  PriorSpec sigmaSqPrior;
  PriorSpec *thetaPrior;

  /*
    Cached indices into the global M_lat@x array.

    These caches are built once after the symbolic pattern of M_lat is known.
    Numeric assembly then becomes direct accumulation into cached M@x slots,
    with no repeated search during MCMC iterations.
  */
  int *nngpCacheStart;
  int *nngpCacheIdx;
  int nngpCacheN;

  int *gpCacheIdx;
  int gpCacheN;

  int *ar1DiagCacheIdx;
  int *ar1OffCacheIdx;

  int *carDiagCacheIdx;
  int *carOffCacheIdx;

  cholmod_sparse *carQ;
  cholmod_factor *carQFac;
  int *carQDiagCacheIdx;
  int *carQOffCacheIdx;

  int *carTimeCacheIdx;
  int carTimeCacheN;

  int *dagarCacheStart;
  int *dagarCacheIdx;
  int dagarCacheN;
} TermState;

typedef struct {
  int dim;
  int *paramIndex;
  int batchLength;
  int batchPos;
  int batchIndex;
  int batchAccept;
  int acceptCount;
  int maxBatchHistory;
  int scalarMode;
  double targetAccept;
  double warmupTargetLower;
  double warmupTargetUpper;
  double warmupNearZero;
  double c0;
  double c1;
  double sigmaSqM;
  double baseSigmaSqM;
  double lastBatchAccept;
  int warmupEnabled;
  int warmupBatchLength;
  int warmupMinBatches;
  int warmupMaxBatches;
  int warmupBatches;
  int warmupNAttempted;
  int warmupNAccepted;
  int warmupStoppedReason;
  double *eta;
  double *etaProp;
  double *z;
  double *Sigma0;
  double *proposalCov;
  double *proposalChol;
  double *warmupStartingEta;
  double *warmupEndingEta;
  double *warmupStartingProposalCov;
  double *warmupEndingProposalCov;
  double *warmupBatchAcceptHistory;
  double *warmupProposalScaleHistory;
  double *batchSamples;
  double *batchAcceptHistory;
  double *proposalScaleHistory;
} CovProposalBlock;

typedef struct {
  int initialized;
  int dim;
  int *paramType;
  int *paramTerm;
  int *paramTheta;
  int *paramBlock;
  int nBlocks;
  CovProposalBlock *blocks;
  int totalBlockAccept;
  int batchLength;
  int batchPos;
  int batchIndex;
  int batchAccept;
  int maxBatchHistory;
  double targetAccept;
  double warmupTargetLower;
  double warmupTargetUpper;
  double warmupNearZero;
  double c0;
  double c1;
  double sigmaSqM;
  double baseSigmaSqM;
  double lastBatchAccept;
  int warmupEnabled;
  int warmupBatchLength;
  int warmupMinBatches;
  int warmupMaxBatches;
  int warmupBatches;
  int warmupNAttempted;
  int warmupNAccepted;
  int warmupStoppedReason;
  double *eta;
  double *etaProp;
  double *z;
  double *Sigma0;
  double *proposalCov;
  double *proposalChol;
  double *warmupStartingEta;
  double *warmupEndingEta;
  double *warmupStartingProposalCov;
  double *warmupEndingProposalCov;
  double *warmupBatchAcceptHistory;
  double *warmupProposalScaleHistory;
  double *batchSamples;
  double *batchAcceptHistory;
  double *proposalScaleHistory;
  double tauSq;
  double *residualVariance;
  double *sigmaSq;
  double **theta;
  double **B;
  double **F;
  double **Q;
  double *logDetQ;
  cholmod_sparse *M_lat;
  cholmod_factor *M_lat_fac;
} CovProposalState;

typedef struct {
  int *rows;
  int len;
  int cap;
} PatternCol;

typedef struct {
  int n;
  int p;
  int q;
  int qLatTotal;
  int nGraphs;
  int nTerms;
  int nSamples;
  int nThetaTotal;

  /* optional in-chain recovery of the stacked latent process */
  int recoverProcess;
  int recoverStart;
  int recoverThin;
  int nRecover;

  int covParamAccept;
  int covParamAttempts;
  int metropolisBatchLength;
  double metropolisTargetAccept;
  int warmupEnabled;
  int warmupBatchLength;
  int warmupMinBatches;
  int warmupMaxBatches;
  double warmupTargetLower;
  double warmupTargetUpper;
  double warmupNearZero;

  /*
    y is the response used by the current Gaussian working model. For the
    ordinary Gaussian likelihood it points at the observed response from R.
    For Polya-Gamma likelihoods it is sampler-owned mutable storage containing
    the pseudo-response for the current iteration. yObserved keeps the original
    observed successes/counts so the PG weights can be refreshed without losing
    the data.
  */
  double *y;
  double *yObserved;
  double *offset;
  int hasOffset;
  int likelihoodFamily;
  int *nTrial;
  double nbSize;
  double *X;
  int betaPriorType; /* 0 = flat, 1 = independent normal */
  double *betaPriorMean;
  double *betaPriorPrecision;

  int *Zp;
  int *Zi;
  double *Zx;
  int Z_nnz;

  GraphState *graphs;
  TermState *terms;

  /*
    The old dense qLatTotal x qLatTotal lookup map is intentionally not used.
    We keep the pointer in the state only to minimize intrusive changes to the
    surrounding code, but it will remain NULL in the sparse/cached workflow.
  */
  int *latentLowerMap;
  cholmod_sparse *M_lat;
  cholmod_factor *M_lat_sym;
  cholmod_factor *M_lat_fac;
  CovProposalState covProp;
  int metropolisBlocking;

  double tauSq;
  double tauSqShape;
  double tauSqScale;
  PriorSpec tauSqPrior;
  double tauTune;
  int residualModel; /* 0 = global tau_sq, 1 = fixed variance, 2 = group IG, 3 = scaled */
  double *obsPrecision;
  double *obsPrecisionFixed;
  int nResidualGroup;
  int *residualGroupIndex;
  double *residualVariance;
  double *residualVarianceShape;
  double *residualVarianceScale;
  PriorSpec *residualVariancePrior;
  double *residualVarianceTune;
  double *residualVarianceMeanLog;
  double *residualVarianceSdLog;
  double *residualScaledVhat;
  double *residualScaledWeight;

  int nRE;
  int *reBlockID;
  int *reBlockSize;
  int maxBlockSize;
  double *sigmaSqRE;
  double *sigmaSqREShape;
  double *sigmaSqREScale;
  double *alphaBlockBuf;

  double *r;
  double *Za;
  double *alpha;
  double *nWork1;
  double *nWork2;

  double *scratch_BF_c;
  double *scratch_BF_C;
  double *scratch_BF_bk;
  int scratch_BF_m;
  int scratch_BF_bk_n;
  int nOmpThreads;

  /*
    Reusable dense work buffers for sparse solves involving M_lat.
    These avoid repeated RHS/z allocation and copy-in overhead when using
    the Matrix-exposed CHOLMOD solve interface.
  */
  double *latentSolveRhs;
  double *latentNoiseZ;

  /*
    Cached indices for the A' W A observation contribution.

    We enumerate the lower term pairs (t1, t2) with t2 <= t1 once, then for
    each observation i cache the corresponding global M@x index.
  */
  int nTermPairs;
  int *termPair1;
  int *termPair2;
  int *AtAIdx;

  cholmod_common cm;
} SamplerState;

/*==========================================================================*/
/* graph helpers                                                            */
/*==========================================================================*/

SEXP getListElement(SEXP list, const char *name);
int getListIndex(SEXP list, const char *name);
void setListElementByName(SEXP list, const char *name, SEXP value);
int as_flag_scalar(SEXP x, const char *where);
int as_nonneg_int_scalar(SEXP x, const char *where);
const char *ordering_label_from_common(cholmod_common *cm);
void configure_cholmod_control(cholmod_common *cm, SEXP backend_r);
int as_int_scalar(SEXP x, const char *where);
double as_real_scalar_strict(SEXP x, const char *where);
const char *as_char_scalar(SEXP x, const char *where);
int compute_n_recover(int nSamples, int recoverProcess, int recoverStart, int recoverThin);
void require_matrix_real(SEXP x, int nr, int nc, const char *where);
void require_matrix_int(SEXP x, int nr, int nc, const char *where);
int small_chol_lower(double *A, int n);
int small_chol_solve_lower(const double *L, const double *b, double *x, int n);
double logDetFactor(cholmod_factor *L);
int cmp_int_asc(const void *a, const void *b);
void pattern_col_append(PatternCol *cols, int col, int row);
void pattern_append_lower(PatternCol *cols, int i, int j);

double graph_distance(GraphState *g, TermState *term, int i, int j, double *u);
void update_nngp_BF(SamplerState *s, GraphState *g, TermState *term);
void update_gp_Q(SamplerState *s, GraphState *g, TermState *term);
double theta_forward(double x, double lower, double upper);
double theta_inverse(double z, double lower, double upper);
double theta_log_jacobian(double x, double lower, double upper);
double logdet_Qw_blocks(SamplerState *s);
int logdet_Qw_blocks_try(SamplerState *s, double *out);
SEXP char_vector_from_term_labels(SamplerState *s);

void build_M_lat_pattern_sparse(SamplerState *s);
void build_M_lat_index_caches(SamplerState *s);
void refresh_observation_precision(SamplerState *s);
void assemble_M_lat_numeric(SamplerState *s);
void apply_A_transpose(SamplerState *s, const double *vObs, double *outLat);
void apply_A_transpose_weighted(SamplerState *s, const double *vObs, double *outLat);
void apply_A(SamplerState *s, const double *vLat, double *outObs);
void solve_M_lat(SamplerState *s, const double *b, double *x);
void apply_Vinv(SamplerState *s, const double *vObs, double *outObs);
double compute_logDetV(SamplerState *s);
double quadform_Vinv(SamplerState *s, const double *vObs);
void apply_A_transpose_multiple(SamplerState *s, const double *vObs, int ncol, double *outLat);
void apply_A_multiple(SamplerState *s, const double *vLat, int ncol, double *outObs);
void solve_M_lat_multiple(SamplerState *s, const double *B, int ncol, double *X);
void apply_Vinv_multiple(SamplerState *s, const double *Vobs, int ncol, double *outObs);
void form_XtVinvX_and_rhs(SamplerState *s, const double *ytilde, double *XtVinvX, double *XtVinvy);
void form_ZtVinvZ_and_rhs(SamplerState *s, const double *ytilde, double *ZtVinvZ, double *ZtVinvy);
void refactor_current_M(SamplerState *s, const char *where);
void update_linear_predictor(SamplerState *s,
                             const double *beta,
                             const double *alpha,
                             const double *w,
                             double *eta);
int is_pg_likelihood(const SamplerState *s);
void update_pg_working_model(SamplerState *s,
                             BayesLogit_rpg_hybrid_t pg,
                             const double *beta,
                             const double *alpha,
                             const double *w);

void init_graph_state_from_backend(GraphState *g, SEXP graph_r);
void init_term_state_from_backend(TermState *term, GraphState *graphs, int nGraphs, SEXP term_r, int n);
void init_sampler_state(SamplerState *s, SEXP backend_r);
void free_cov_proposal_state(SamplerState *s);
void free_sampler_state(SamplerState *s);
void init_cov_proposal_state(SamplerState *s);

void scatter_Z_col_to_dense(SamplerState *s, int col, double *out);
void sparse_Z_mult(SamplerState *s, const double *alpha, double *out);
void sparse_Zt_mult(SamplerState *s, const double *v, double *out);
double sparse_Z_col_dot_dense(SamplerState *s, int col, const double *v);
void update_mean_without_alpha(SamplerState *s, const double *beta, double *out);
void update_sigmaSqRE_gibbs(SamplerState *s);
double log_ig_kernel_collapsed(double x, double shape, double scale);
double log_cov_params_target_collapsed(SamplerState *s, const double *resid);
int log_cov_params_target_collapsed_try(SamplerState *s, const double *resid, double *out);
int update_cov_params_mh_collapsed(SamplerState *s, const double *resid);
void draw_gaussian_from_precision(double *Prec, double *rhs, int dim, double *draw);
void update_beta_draw_collapsed(SamplerState *s, const double *alpha, double *betaDraw);
void update_alpha_draw_collapsed(SamplerState *s, const double *beta, double *alphaDraw);
SEXP list_B_by_term(SamplerState *s);
SEXP list_F_by_term(SamplerState *s);
int should_recover_iteration(SamplerState *s, int iter1);
void recover_w_draw_collapsed(SamplerState *s, const double *resid, double *wDraw);
int term_prior_block_nnz(TermState *term, GraphState *g);
SEXP build_term_description_object(SamplerState *s, SEXP backend_r, int nSamplesDone);
void print_rule(char ch, int n);
void print_term_description(SEXP td_r);

extern "C" SEXP stLMM_collapsed_sampler(SEXP backend_r);
extern "C" SEXP stLMM_recover_w(SEXP backend_r,
                                SEXP beta_samples_r,
                                SEXP alpha_samples_r,
                                SEXP tau_samples_r,
                                SEXP residual_variance_samples_r,
                                SEXP process_sigma_samples_r,
                                SEXP theta_samples_r,
                                SEXP draw_index_r,
                                SEXP pg_iter_r);
extern "C" SEXP stLMM_predict_nngp_joint_false(SEXP support_r,
                                               SEXP new_coords_r,
                                               SEXP neighbor_index_r,
                                               SEXP w_fit_r,
                                               SEXP sigma_sq_r,
                                               SEXP theta_r,
                                               SEXP cov_model_r,
                                               SEXP n_omp_threads_r);
extern "C" SEXP stLMM_predict_nngp_vecchia_joint(SEXP coords_all_r,
                                                 SEXP n_fit_r,
                                                 SEXP neighbor_index_r,
                                                 SEXP neighbor_count_r,
                                                 SEXP w_fit_r,
                                                 SEXP sigma_sq_r,
                                                 SEXP theta_r,
                                                 SEXP cov_model_r,
                                                 SEXP n_omp_threads_r);

#endif
