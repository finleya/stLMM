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
