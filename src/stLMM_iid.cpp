#include "stLMM_internal.h"

/* Density on log variance: the IG density plus its transformation Jacobian.
   There is one prior/Jacobian contribution per variance, not per group. */
double log_iid_covariance_target(SamplerState *s)
{
  int b;
  double out;

  out = 0.0;
  for(b = 0; b < s->nRE; b++)
    if(s->reCovDim[b] == 0)
      out += log_ig_kernel_collapsed(s->sigmaSqRE[b], s->sigmaSqREShape[b],
                                      s->sigmaSqREScale[b]) + std::log(s->sigmaSqRE[b]);
  return out;
}
