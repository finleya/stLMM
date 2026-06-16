#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP get_cor_models(void);
extern SEXP mkNNIndx(SEXP coords_r,
                     SEXP n_r,
                     SEXP m_r,
                     SEXP r_r,
                     SEXP nThreads_r);
extern SEXP mkNNIndxBrute(SEXP coords_r,
                          SEXP n_r,
                          SEXP m_r,
                          SEXP r_r,
                          SEXP nThreads_r);
extern SEXP stLMM_collapsed_sampler(SEXP backend_r);
extern SEXP stLMM_recover_w(SEXP backend_r,
                           SEXP beta_samples_r,
                           SEXP alpha_samples_r,
                           SEXP tau_samples_r,
                           SEXP residual_variance_samples_r,
                           SEXP process_sigma_samples_r,
                           SEXP theta_samples_r,
                           SEXP draw_index_r,
                           SEXP pg_iter_r);
extern SEXP stLMM_predict_nngp_joint_false(SEXP support_r,
                                          SEXP new_coords_r,
                                          SEXP neighbor_index_r,
                                          SEXP w_fit_r,
                                          SEXP sigma_sq_r,
                                          SEXP theta_r,
                                          SEXP cov_model_r,
                                          SEXP n_omp_threads_r);
extern SEXP stLMM_predict_nngp_vecchia_joint(SEXP coords_all_r,
                                            SEXP n_fit_r,
                                            SEXP neighbor_index_r,
                                            SEXP neighbor_count_r,
                                            SEXP w_fit_r,
                                            SEXP sigma_sq_r,
                                            SEXP theta_r,
                                            SEXP cov_model_r,
                                            SEXP n_omp_threads_r);
extern SEXP stLMM_nngp_prediction_neighbors(SEXP support_r,
                                           SEXP new_coords_r,
                                           SEXP m_r,
                                           SEXP cov_model_r,
                                           SEXP st_scale_r,
                                           SEXP n_omp_threads_r);

static const R_CallMethodDef CallEntries[] = {
    {"get_cor_models", (DL_FUNC) &get_cor_models, 0},
    {"mkNNIndx", (DL_FUNC) &mkNNIndx, 5},
    {"mkNNIndxBrute", (DL_FUNC) &mkNNIndxBrute, 5},
    {"stLMM_collapsed_sampler", (DL_FUNC) &stLMM_collapsed_sampler, 1},
    {"stLMM_recover_w", (DL_FUNC) &stLMM_recover_w, 9},
    {"stLMM_predict_nngp_joint_false", (DL_FUNC) &stLMM_predict_nngp_joint_false, 8},
    {"stLMM_predict_nngp_vecchia_joint", (DL_FUNC) &stLMM_predict_nngp_vecchia_joint, 9},
    {"stLMM_nngp_prediction_neighbors", (DL_FUNC) &stLMM_nngp_prediction_neighbors, 6},
    {NULL, NULL, 0}
};

void R_init_stLMM(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
