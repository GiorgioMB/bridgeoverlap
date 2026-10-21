// Radon, complete pooling (posteriordb: radon_pooled), written with `target +=` and
// every normalizing constant (half-normal priors normalized by + log 2), so that the
// marginal likelihood is exact (bridgeoverlap analysis/R/radon.R computes it by quadrature).
data {
  int<lower=0> N;
  vector[N] floor_measure;
  vector[N] log_radon;
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma_y;
}
model {
  target += normal_lpdf(sigma_y | 0, 1) + log2();
  target += normal_lpdf(alpha | 0, 10) + normal_lpdf(beta | 0, 10);
  target += normal_lpdf(log_radon | alpha + beta * floor_measure, sigma_y);
}
