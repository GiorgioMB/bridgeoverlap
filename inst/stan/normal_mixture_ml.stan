// Two-component normal mixture with unit variances (posteriordb: normal_2-normal_mixture),
// unknown weight and means, exchangeable priors, so the posterior is invariant under
// relabelling: (theta, mu1, mu2) -> (1 - theta, mu2, mu1). Written with `target +=` and
// every constant; the uniform prior on theta has density 1.
data {
  int<lower=0> N;
  array[N] real y;
}
parameters {
  real<lower=0, upper=1> theta;
  array[2] real mu;
}
model {
  target += normal_lpdf(mu | 0, 10);
  for (n in 1 : N) {
    target += log_mix(theta, normal_lpdf(y[n] | mu[1], 1.0), normal_lpdf(y[n] | mu[2], 1.0));
  }
}
