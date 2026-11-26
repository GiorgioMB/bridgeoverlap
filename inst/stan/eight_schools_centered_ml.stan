// Eight schools, centred parameterization (posteriordb: eight_schools_centered),
// written with `target +=` so that cmdstanr's log_prob() keeps every constant.
// The half-Cauchy prior is normalized (+ log 2), so the marginal likelihood is
// exactly log Z = -31.3113 (bridgeoverlap::eight_schools_log_evidence()).
data {
  int<lower=0> J;                 // number of schools
  array[J] real y;                // estimated treatment effects
  array[J] real<lower=0> sigma;   // standard errors
}
parameters {
  array[J] real theta;            // school effects
  real mu;                        // population mean
  real<lower=0> tau;              // population sd
}
model {
  target += cauchy_lpdf(tau | 0, 5) + log2();  // half-Cauchy(0, 5)
  target += normal_lpdf(theta | mu, tau);
  target += normal_lpdf(y | theta, sigma);
  target += normal_lpdf(mu | 0, 5);
}
