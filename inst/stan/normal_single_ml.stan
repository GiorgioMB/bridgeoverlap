// One-component counterpart of normal_mixture_ml.stan, for Bayes factors: y ~ N(mu, 1),
// mu ~ N(0, 10^2); log Z is available in closed form.
data {
  int<lower=0> N;
  array[N] real y;
}
parameters {
  real mu;
}
model {
  target += normal_lpdf(mu | 0, 10);
  target += normal_lpdf(y | mu, 1.0);
}
