Native latent-system validation
===============================

From the package source directory, run:
  Rscript tests/native/run-latent-checks.R

This separate validation suite requires the package build dependencies. It
builds and installs a temporary copy, adding native probes only to that copy.
No diagnostic native entry point is registered in the distributed package.
Logs and results are retained in the temporary directory printed by the runner.

The suite compares joint IID/structured sparse calculations with independent
dense observation-space covariance calculations. It checks mappings, prior
precision, marginal precision and log determinants, fixed-effect information,
conditional means and covariance of simulated latent draws, and numeric cache
restoration after changing covariance blocks. It covers independent, correlated,
crossed, nested, and parent-specific IID effects; unused groups; missing
responses; zero and collinear coefficient columns; unequal residual variances;
known offsets; and models without fixed effects. Structured terms include AR(1),
NNGP, GP, CAR, DAGAR, and the space-time CAR/DAGAR variants.

Covariance parameters are fixed in the latent-reference checks. Each case also verifies the
joint Gaussian distribution of fixed and random coefficients using 3,000 joint
schedule iterations and an independent dense posterior precision. The Gaussian
chain probe can also compare both coefficient update schedules. These checks do
not establish correctness, convergence, or speed of covariance-parameter
Metropolis updates by themselves. The scalar covariance suite separately checks
IG posterior moments against one-dimensional quadrature and checks accepted and
rejected proposal state restoration under all blocking choices, with warmup on
and off. The IW checks compare the transformed prior density with independent R
calculations and finite-difference Jacobians, and verify prior-only draws against
analytical covariance, precision, and log-determinant moments. Singular and
indefinite starting/scale matrices must error; no covariance jitter or ridging is
used. Extreme proposals must be rejected without modifying the current state.
Zero covariance tuning must preserve the supplied matrix exactly.

The posterior comparison suite simulates 18 Gaussian datasets covering scalar,
independent, correlated, crossed, nested, and parent-specific IID terms, without
structured effects and with AR(1) or NNGP effects. Both schedules estimate the
covariance parameters under matched IG/IW priors. Three chains per schedule
provide posterior means, 95% intervals, R-hat, bulk/tail ESS, and timings including
warmup and numeric refactorization. It compares means using their combined Monte
Carlo standard errors, including fixed-plus-IID contributions and full predictors.
These are small-system checks, not a guarantee of faster sampling for all models.
The covariance suites require the suggested posterior package for diagnostics.

Augmented-model checks cover binomial logit, probit, and negative-binomial
likelihoods with fixed or estimated size. Ten ordinary-likelihood quadrature
problems independently check posterior coefficient and covariance moments,
including scalar IG, correlated IW, estimated size, and nonzero normal fixed-effect
priors. The schedule
checks exercise all six blocking choices with warmup on/off, offsets, missing
responses, unused groups, AR(1), and NNGP effects. Saving and thinning structured
draws must preserve both the parameter samples and RNG state. A size-only warmup
case holds all IID variances fixed.

The same posterior comparison harness covers 24 simulated augmented models,
with 120 observations and 20 main grouping levels. Each uses three chains per
mode, 6,000 iterations per chain, and discards the first 1,500. Parent-specific
covariance cases use by_term blocking and 12,000 iterations, discarding 3,000;
the joint-proposal short runs can mix slowly for these weakly identified models.
These are test settings, not changes to the software defaults. To run one family
separately, set STLMM_TEST_FAMILY to binomial, probit, negative_binomial, or
estimated_size when running iid-augmented-posterior-reference.R with the library
directory and output RDS path as arguments. STLMM_TEST_IID optionally selects an
IID case (scalar, independent, covariance, crossed, nested, or by).

The public fitting/recovery reference (iid-output-reference.R) checks the joint
Gaussian posterior under fixed covariance parameters in 20 cases: ten IID
formulations, with equal/unequal residual variances and no process, AR(1), NNGP,
or GP effects. It includes offsets, missing responses, unused groups, collinear
columns, and a nonzero normal fixed-effect prior. It compares the full retained
(beta, w, alpha) posterior to independent dense calculations after recover().

Controlled performance comparisons
---------------------------------
The optional iid-benchmark.R worker is separate from correctness checks and is
not run by run-latent-checks.R or testthat. From the source root, prepare the
synthetic datasets with an installed library and an output directory:
  Rscript tests/native/iid-benchmark.R LIBRARY OUTPUT prepare
Then run individual chains in fresh processes, sequentially, for example:
  Rscript tests/native/iid-benchmark.R LIBRARY OUTPUT run correlated300 separate 1 6000
  Rscript tests/native/iid-benchmark.R LIBRARY OUTPUT run correlated300 collapsed 1 6000
Repeat for chain IDs 1, 2, and 3. A final optional argument selects Metropolis
blocking (default joint). The cases.csv file lists the twelve scenarios.
Use one BLAS/OpenMP thread, and on Linux use the same CPU affinity and wrap each
worker with /usr/bin/time -v to record whole-process peak RSS. Alternate strategy
order across chains/cases. Keep competing workloads out of timed runs.

The worker records public fitting time, native initialization/warmup/sampling
phases, post-fit recovery time, sparse counts, approximate R object sizes, and
selected posterior draws. All covariance parameters are estimated. Starts are
common between strategies and dispersed across chains. The first quarter of
post-warmup samples is excluded from diagnostics; Gaussian spatial recovery uses
every remaining draw under both strategies. Include that cost when comparing
usable predictors. Summarize all three chains together with rank-normalized
R-hat, bulk/tail ESS, and posterior means/intervals under matched priors.
Covariance, raw coefficients, combined coefficients, and predictor ESS/sec should
be reported separately. A single simulated dataset per case does not establish
universal performance or interval coverage. Whole-process RSS includes libraries,
retained samples and extraction; R object.size can count shared components more
than once, and sparse counts alone are not a total-memory measure.
