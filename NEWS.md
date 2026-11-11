# stLMM 0.0.5

## Changes

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
