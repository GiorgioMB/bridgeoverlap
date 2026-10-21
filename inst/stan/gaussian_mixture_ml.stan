// Two-component Gaussian mixture in two dimensions with a known normalizing
// constant. For m = 3, s_narrow = 0.5 and w = 0.5 the target is
//   f(x) = exp(log_Z) * { w N(x | (-m, 0), I) + (1 - w) N(x | (m, 0), diag(s_narrow^2, 1)) }.
// For large m Stan's sampler cannot move between the modes, so the fraction of
// draws in each mode is set by the initial values: a real sampler failure with
// a known answer.
data {
  real<lower=0> m;
  real<lower=0> s_narrow;
  real<lower=0, upper=1> w;
  real log_Z;
}
parameters {
  vector[2] x;
}
model {
  target += log_Z
            + log_mix(w,
                      normal_lpdf(x[1] | -m, 1) + std_normal_lpdf(x[2]),
                      normal_lpdf(x[1] | m, s_narrow) + std_normal_lpdf(x[2]));
}
