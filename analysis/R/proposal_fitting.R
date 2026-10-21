# Proposal fitting with a t_7 target (Section 7.3, Tables E4, E4c and prop) ---------------
# Target: independent coordinates 0.5 + s_k t_7, s_k equally spaced in [1, 2].
# Proposal: diagonal Gaussian N(m, diag(v)) fitted by moments, k = 2d parameters.

t7_target <- function(d, logZ = 1.1, nu = 7) {
  s <- if (d == 1) 1.5 else seq(1, 2, length.out = d)
  list(d = d, s = s, nu = nu, logZ = logZ,
       r = function(n) sweep(matrix(stats::rt(n * d, nu), n, d), 2, s, "*") + 0.5,
       lf = function(X) logZ + rowSums(stats::dt(sweep(X - 0.5, 2, s, "/"), nu, log = TRUE)) - sum(log(s)),
       score = function(X) {                      # gradient of log f, n x d
         Z <- sweep(X - 0.5, 2, s, "/")
         sweep(-(nu + 1) * Z / (nu + Z^2), 2, s, "/")
       },
       m0 = rep(0.5, d), v0 = s^2 * nu / (nu - 2))
}
moments <- function(X) list(m = colMeans(X), v = apply(X, 2, stats::var))
prop_from <- function(E, m, v) sweep(sweep(E, 2, sqrt(v), "*"), 2, m, "+")
ratios <- function(tg, X, Y, m, v) list(ap = tg$lf(X) - ldiag(X, m, v), aq = tg$lf(Y) - ldiag(Y, m, v))

# Skewed target for the residual second-order term of cross-fitting: independent
# coordinates 0.5 + s_k (G_k - alpha) / sqrt(alpha), G_k ~ Gamma(alpha, 1), so that the
# moment-matched Gaussian is not the overlap-maximizing one in the mean direction
# (the gradient of H at theta_0 is non-zero).
gamma_target <- function(d, logZ = 1.1, alpha = 2) {
  s <- if (d == 1) 1.5 else seq(1, 2, length.out = d)
  list(d = d, s = s, logZ = logZ,
       r = function(n) sweep(matrix(stats::rgamma(n * d, alpha, 1) - alpha, n, d) / sqrt(alpha), 2, s, "*") + 0.5,
       lf = function(X) {
         G <- alpha + sqrt(alpha) * sweep(X - 0.5, 2, s, "/")
         out <- logZ + rowSums(stats::dgamma(pmax(G, 0), alpha, 1, log = TRUE)) + d * log(sqrt(alpha)) - sum(log(s))
         out[rowSums(G <= 0) > 0] <- -Inf
         out
       },
       m0 = rep(0.5, d), v0 = s^2)
}
