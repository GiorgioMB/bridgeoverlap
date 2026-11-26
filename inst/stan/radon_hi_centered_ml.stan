// Radon, varying intercept with the county uranium predictor, centred (posteriordb:
// radon_hierarchical_intercept_centered), with every normalizing constant.
// Unconstrained dimension J + 5.
data {
  int<lower=0> J;
  int<lower=0> N;
  array[N] int<lower=1, upper=J> county_idx;
  vector[N] log_uppm;
  vector[N] floor_measure;
  vector[N] log_radon;
}
parameters {
  vector[J] alpha;
  vector[2] beta;
  real mu_alpha;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_y;
}
model {
  target += normal_lpdf(sigma_alpha | 0, 1) + normal_lpdf(sigma_y | 0, 1) + 2 * log2();
  target += normal_lpdf(mu_alpha | 0, 10) + normal_lpdf(beta | 0, 10);
  target += normal_lpdf(alpha | mu_alpha, sigma_alpha);
  target += normal_lpdf(log_radon | alpha[county_idx] + log_uppm * beta[1] + floor_measure * beta[2], sigma_y);
}
