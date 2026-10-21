// Radon, varying intercept, non-centred (posteriordb: radon_variable_intercept_noncentered),
// with every normalizing constant. Unconstrained dimension J + 4.
data {
  int<lower=0> J;
  int<lower=0> N;
  array[N] int<lower=1, upper=J> county_idx;
  vector[N] floor_measure;
  vector[N] log_radon;
}
parameters {
  vector[J] alpha_raw;
  real beta;
  real mu_alpha;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_y;
}
model {
  vector[J] alpha = mu_alpha + sigma_alpha * alpha_raw;
  target += normal_lpdf(sigma_y | 0, 1) + normal_lpdf(sigma_alpha | 0, 1) + 2 * log2();
  target += normal_lpdf(mu_alpha | 0, 10) + normal_lpdf(beta | 0, 10);
  target += std_normal_lpdf(alpha_raw);
  target += normal_lpdf(log_radon | alpha[county_idx] + floor_measure * beta, sigma_y);
}
