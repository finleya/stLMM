# Reuse the Gaussian comparison harness for the four augmented configurations.
Sys.setenv(STLMM_TEST_AUGMENTED='1')
source('tests/native/iid-posterior-reference.R')
