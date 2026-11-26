// Radon, varying intercept and slope, centred (posteriordb:
// radon_variable_intercept_slope_centered), with every normalizing constant.
// Unconstrained dimension 2 J + 5.
data {
  int<lower=0> N;
  int<lower=0> J;
  array[N] int<lower=1, upper=J> county_idx;
  vector[N] floor_measure;
  vector[N] log_radon;
}
parameters {
  real<lower=0> sigma_y;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
  vector[J] alpha;
  vector[J] beta;
  real mu_alpha;
  real mu_beta;
}
model {
  target += normal_lpdf(sigma_y | 0, 1) + normal_lpdf(sigma_beta | 0, 1) + normal_lpdf(sigma_alpha | 0, 1) + 3 * log2();
  target += normal_lpdf(mu_alpha | 0, 10) + normal_lpdf(mu_beta | 0, 10);
  target += normal_lpdf(alpha | mu_alpha, sigma_alpha);
  target += normal_lpdf(beta | mu_beta, sigma_beta);
  target += normal_lpdf(log_radon | alpha[county_idx] + floor_measure .* beta[county_idx], sigma_y);
}
