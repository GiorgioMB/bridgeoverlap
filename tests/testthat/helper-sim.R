# Gaussian target f = Z N(0, I_d) and a fixed Gaussian proposal N(m, s^2 I_d).
sim_gauss <- function(S1 = 2000, S2 = S1, d = 3, m = 0.3, s = 1.3, logZ = 1.7) {
  x <- matrix(stats::rnorm(S1 * d), S1, d)
  y <- matrix(stats::rnorm(S2 * d, m, s), S2, d)
  la <- function(t) logZ + rowSums(stats::dnorm(t, log = TRUE)) -
    rowSums(stats::dnorm(t, m, s, log = TRUE))
  list(a_post = la(x), a_prop = la(y), logZ = logZ, x = x, y = y)
}
