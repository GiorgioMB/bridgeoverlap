// Linear regression with the conjugate normal-inverse-gamma prior
//   beta | sigma^2 ~ N(0, sigma^2 v0 I_K),   sigma^2 ~ Inv-Gamma(a0, b0),
//   y | beta, sigma ~ N(X beta, sigma^2 I_n).
// The marginal likelihood has a closed form (bridgeoverlap::conjugate_log_evidence()),
// so bridge sampling estimates can be checked against the truth in any
// dimension K + 1.
data {
  int<lower=1> n;
  int<lower=1> K;
  matrix[n, K] X;
  vector[n] y;
  real<lower=0> v0;
  real<lower=0> a0;
  real<lower=0> b0;
}
parameters {
  vector[K] beta;
  real<lower=0> sigma;
}
model {
  // prior on sigma^2 expressed in sigma: |d sigma^2 / d sigma| = 2 sigma
  target += inv_gamma_lpdf(square(sigma) | a0, b0) + log(2 * sigma);
  target += normal_lpdf(beta | 0, sigma * sqrt(v0));
  target += normal_lpdf(y | X * beta, sigma);
}
