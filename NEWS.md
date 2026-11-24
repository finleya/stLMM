# stLMM (development version)

- Added an Additional Article comparing small-area estimates with `sae` and
  `survey`, covering Fay–Herriot models, unit-level population prediction,
  and prediction for unsampled counties.

- Gaussian fits with `collapse_iid = TRUE` now defer IID and structured effect
  draws to `recover()`, including IID-only fits. Recovery samples both blocks
  jointly for the requested iterations. Call it before prediction, fitted
  values, or combined coefficient summaries. This avoids latent draws and IID
  sample storage during fitting; recovery costs depend on the selected draws.
- Recovered IID coefficients, predictions, sample extraction, and diagnostics
  align with the original parameter iterations through `recover_iter`.
  Separate IID sampling and augmented non-Gaussian updates are unchanged;
  recovery of previously saved fits with in-chain IID samples remains supported.
- Recovery releases native resources on errors and interrupts. Recovered
  summaries report the draw count actually used, and IID-only recovery plots
  apply `burnin` to recovered rows. Gaussian multi-chain recovery objects can
  be recovered again on a different iteration grid.
- Streamlined the grouped-effects comparison with lme4, with consistent
  posterior intervals, a Gaussian combined-coefficient diagnostic example,
  and optional detailed comparisons. Standard-deviation summaries transform
  individual variance draws before averaging. The getting-started vignette
  now explains recovery for collapsed Gaussian IID effects.

# stLMM 0.0.7

## New features

- Added experimental joint IID collapsing through `collapse_iid = TRUE`, with
  the existing adaptive Metropolis controls.
  `FALSE` preserves separate IID updates and remains the default. `NULL`
  selects collapsing for IID-only models and separate updates when structured
  effects are present. Both strategies use the same output and prediction API.

- Scalar IID variances now accept `ig()`, `log_normal()`, `gamma_dist()`,
  `half_normal()`, and `half_t()` under either sampling strategy. Half-normal
  and half-t priors are on standard deviations. Separate IID effects retain
  Gibbs updates for IG variances and IW covariance matrices; other scalar
  priors use conditional Metropolis updates. An explicit prior is required.

## Documentation and validation

- Updated manual pages and the technical document for both IID strategies,
  covariance transformations, tuning, and scalar variance prior choices.
- Expanded the grouped-effects article with a three-way comparison of separate
  IID sampling, collapsed IID sampling, and lme4. The example reports parameter
  recovery, credible intervals, combined-coefficient diagnostics, elapsed time,
  and effective samples per second.
- Added numerical reference and synthetic-data checks for Gaussian, binomial,
  probit, and negative-binomial models, including IID effects combined with
  structured processes, missing responses, recovery, and prediction.

# stLMM 0.0.6

## New features

- Added correlated grouped coefficients through `iid(group, ~ 1 + x)` with
  inverse-Wishart covariance priors. The `cov_model = "independent"` form
  shares the scalar sampler with separate IID intercept and slope terms.
- Added parent-specific IID variances and covariance matrices through `by`,
  and explicit new-group prediction through `allow_new_levels`.
- Added combined fixed-plus-IID coefficient summaries and optional
  rank-normalized R-hat, bulk ESS, and tail ESS diagnostics.
- Added estimation of negative-binomial size with an explicit Gamma prior
  and a scalar adaptive Metropolis update. Prediction and log likelihoods
  use the corresponding posterior size draws. Fixed-size models remain
  supported.

## Fixes and documentation

- Prediction no longer requires the response column in `newdata` and retains
  fitted fixed-effect bases for prediction.
- Updated manual pages, the technical document, and articles for grouped
  effects, diagnostics, prediction, and estimated negative-binomial size.

# stLMM 0.0.5

## Changes

- Use a size_t counter when iterating over CHOLMOD supernodes, removing a signed/unsigned comparison warning.

- CRAN source builds exclude local article caches and the data downloads used only by website articles.

- CAR terms now retain every node in a supplied graph even when some areas
  have no rows in the fitting data. Missing-response placeholder rows are no
  longer required for these areas; recovery and prediction retain the full
  graph support. This applies to proper and Leroux CAR terms, including
  spatially varying coefficients.

- Fits now record monotonic elapsed seconds for compiled initialization,
  proposal-calibration warmup, and the post-warmup MCMC loop separately, in
  `fit$timing$phases`. Sampling time includes all work in the loop, including
  repeated factorizations and any requested in-chain latent-effect recovery.
- The former `timing$sampler` field is renamed `timing$backend` to describe
  its existing scope accurately: the entire compiled call. Multi-chain fits
  use `backend_by_chain`/`backend_total` and add
  `phases_by_chain`/`phases_total`. The old timing names are not aliases.
  Saved fits from builds without phase timing cannot supply sampling-only
  timings. `Config/stLMM/Timing: phases-v1` identifies this capability in
  DESCRIPTION; the package version remains 0.0.5 during development.
- Updated printed timing labels and the runtime, blocking-diagnostics, and
  NIMBLE-comparison articles to distinguish sampling from broader elapsed times.

- CAR graph connectivity checks now preallocate their traversal queue, avoiding
  repeated vector copying on large graphs while preserving component labels
  and graph validation. This also benefits DAGAR and space-time models that
  use the same graph constructor.

# stLMM 0.0.4

## New features

- Added `family = "probit"` for binary probit mixed models using Albert-Chib latent-normal augmentation.
- Added fitted values, prediction, log-likelihood, recovery, tests, and documentation for probit models.
- Added support for structured probit process models using saved in-chain latent process draws.

## Changes

- NNGP prediction with `joint = TRUE` now defaults to
  `joint_method = "vecchia"` instead of the dense full method. This avoids
  building a dense covariance over new prediction nodes for moderate and large
  prediction sets. Use `joint_method = "full"` to request the previous dense
  behavior for small prediction sets.
- `stats::binomial(link = "probit")` now dispatches to the probit likelihood.
- Unsupported `stats::binomial()` links such as `"cloglog"` and `"cauchit"` now error instead of silently using the logit/Polya-Gamma path.

# stLMM 0.0.3

- Existing release.
