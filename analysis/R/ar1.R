# Stationary Gaussian AR(1) posterior chains (Section 7.2) --------------------------------
#
# Each coordinate follows x_t = rho x_{t-1} + sqrt(1 - rho^2) e_t, so the chain is
# stationary for p = N(0, I_d) and has spectral gap gamma = 1 - rho. The exact
# integrated autocorrelation time of a function of the chain follows from the
# spectral (Mehler) expansion of the AR(1) kernel:
#   * for F2 (d = 1), Hermite polynomials h_k are eigenfunctions with eigenvalue
#     rho^k;
#   * for the radial family F1, the radial part of the 2k-th Wiener chaos is
#     spanned by the generalized Laguerre polynomial L_k^{(d/2-1)}(|x|^2 / 2),
#     with eigenvalue rho^(2k).
# Hence sigma_P^2(u) = sum_k c_k^2 (1 + lambda_k) / (1 - lambda_k) with c_k the
# orthonormal expansion coefficients of u.

gauss_hermite_prob <- function(n) {
  off <- sqrt(seq_len(n - 1L))
  J <- matrix(0, n, n)
  J[cbind(seq_len(n - 1L), 2:n)] <- off
  J[cbind(2:n, seq_len(n - 1L))] <- off
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = e$vectors[1L, ]^2)
}

gauss_laguerre_gen <- function(n, alpha) {
  k <- seq_len(n - 1L)
  J <- diag(2 * (0:(n - 1L)) + alpha + 1)
  off <- sqrt(k * (k + alpha))
  J[cbind(k, k + 1L)] <- off
  J[cbind(k + 1L, k)] <- off
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = e$vectors[1L, ]^2)
}

# Orthonormal Hermite coefficients E[u(X) h_k(X)], X ~ N(0, 1), k = 1..K.
hermite_coefs <- function(u, K = 160L, n = 320L) {
  q <- gauss_hermite_prob(n)
  x <- q$x; w <- q$w
  ux <- u(x)
  h_prev <- rep(1, n); h <- x
  out <- numeric(K)
  out[1L] <- sum(w * ux * h)
  for (k in seq_len(K - 1L)) {
    h_new <- (x * h - sqrt(k) * h_prev) / sqrt(k + 1)
    h_prev <- h; h <- h_new
    out[k + 1L] <- sum(w * ux * h)
  }
  out
}

# Orthonormal Laguerre coefficients E[u(Y) L~_k(Y)], Y ~ Gamma(alpha + 1, 1).
laguerre_coefs <- function(u, alpha, K = 160L, n = 320L) {
  q <- gauss_laguerre_gen(n, alpha)
  keep <- q$w > 1e-290
  y <- q$x[keep]; w <- q$w[keep]
  uy <- u(y)
  L_prev <- rep(1, length(y))
  L <- (alpha + 1 - y) / sqrt(alpha + 1)
  out <- numeric(K)
  out[1L] <- sum(w * uy * L)
  for (k in seq_len(K - 1L)) {
    L_new <- ((2 * k + alpha + 1 - y) * L - sqrt(k * (k + alpha)) * L_prev) /
      sqrt((k + 1) * (k + alpha + 1))
    L_prev <- L; L <- L_new
    out[k + 1L] <- sum(w * uy * L)
  }
  out
}

# Exact tau_P(pi*) for the posterior chain and tau_W (Definition 2) for S1 = S2.
ar1_tau <- function(fam, pop, rho) {
  if (fam$id == "F2") {
    cf <- hermite_coefs(function(x) pop$pistar(x))
    lam <- rho^seq_along(cf)
  } else if (fam$id == "F1") {
    alpha <- fam$d / 2 - 1
    cf <- laguerre_coefs(function(y) pop$pistar(2 * y), alpha)
    lam <- rho^(2 * seq_along(cf))
  } else {
    stop("Spectral expansions are implemented for F1 and F2 only.")
  }
  var_series <- sum(cf^2)
  lrv <- sum(cf^2 * (1 + lam) / (1 - lam))
  tau_P <- if (rho == 0) 1 else lrv / var_series
  tau_W <- (pop$var_g + pop$var_p * tau_P) / (pop$var_g + pop$var_p)
  list(tau_P = tau_P, tau_W = tau_W, var_series = var_series, var_quad = pop$var_p)
}

# u-values of R independent stationary AR(1) chains of length n in d dimensions
# (u = theta for d = 1, u = |theta|^2 otherwise); returns an n x R matrix.
ar1_u <- function(n, d, rho, R) {
  E <- matrix(stats::rnorm(n * d * R), n, d * R)
  if (rho > 0) {
    E[-1L, ] <- sqrt(1 - rho^2) * E[-1L, ]
    E <- stats::filter(E, rho, method = "recursive")
  }
  E <- unclass(E)
  if (d == 1L) return(matrix(E, n, R))
  U <- matrix(0, n, R)
  for (j in seq_len(d)) U <- U + E[, (seq_len(R) - 1L) * d + j]^2
  U
}
