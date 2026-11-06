# stLMM 0.0.5

## Changes

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
