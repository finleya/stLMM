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

Covariance parameters are fixed in these checks. Each case also verifies the
joint Gaussian distribution of fixed and random coefficients using 3,000 joint
schedule iterations and an independent dense posterior precision. The Gaussian
chain probe can also compare both coefficient update schedules. These checks do
not establish correctness, convergence, or speed of covariance-parameter
Metropolis updates.
