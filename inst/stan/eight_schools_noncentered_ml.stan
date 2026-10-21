// Eight schools, non-centred parameterization (posteriordb:
// eight_schools_noncentered), written with `target +=` so that every constant
// is kept. Same model and evidence as the centred version: log Z = -31.3113.
data {
  int<lower=0> J;
  array[J] real y;
  array[J] real<lower=0> sigma;
}
parameters {
  vector[J] theta_trans;          // standardized school effects
  real mu;
  real<lower=0> tau;
}
transformed parameters {
  vector[J] theta = theta_trans * tau + mu;
}
model {
  target += std_normal_lpdf(theta_trans);
  target += normal_lpdf(y | theta, sigma);
  target += normal_lpdf(mu | 0, 5);
  target += cauchy_lpdf(tau | 0, 5) + log2();  // half-Cauchy(0, 5)
}
