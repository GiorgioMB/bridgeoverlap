# Finite-sample confidence interval ----------------------------
#
# For c in R let V1(c), V2(c) be the unbiased sample variances of the
# responsibilities of the posterior and of the proposal draws, L = log(8/delta)
# and
#   rho(c) = sum_k { sqrt(2 S_k V_k(c) L) + 7 S_k L / (3 (S_k - 1)) }.
# The interval is the smallest closed interval containing
# {c : |W(c)| <= rho(c)}; under independent sampling it covers log Z with
# probability at least 1 - delta for every sample size and every (p, g).

# With K folds (cross-fitting), each fold's estimating function is, conditionally on
# the draws its proposal was fitted on, a difference of independent sums of
# [0, 1]-valued variables with mean zero at log Z, so the bound is applied to every
# fold at level delta / K and the radii are added.
.eb_F <- function(c, a_post, a_prop, lr, L, const, fold_post = NULL, fold_prop = NULL) {   # vectorized over c
  K <- if (is.null(fold_post)) 1L else max(fold_post)
  fp <- if (K == 1L) rep(1L, length(a_post)) else fold_post; fq <- if (K == 1L) rep(1L, length(a_prop)) else fold_prop
  out <- numeric(length(c))
  for (s in split(seq_along(c), ceiling(seq_along(c) / 50))) {
    P <- 1 / (1 + exp(outer(lr + a_post, c[s], `-`))); Q <- 1 - 1 / (1 + exp(outer(lr + a_prop, c[s], `-`)))
    W <- colSums(Q) - colSums(P); rho <- 0
    for (k in seq_len(K)) {
      i1 <- fp == k; i2 <- fq == k; S1 <- sum(i1); S2 <- sum(i2)
      v1 <- (colSums(P[i1, , drop = FALSE]^2) - colSums(P[i1, , drop = FALSE])^2 / S1) / (S1 - 1)
      v2 <- (colSums(Q[i2, , drop = FALSE]^2) - colSums(Q[i2, , drop = FALSE])^2 / S2) / (S2 - 1)
      rho <- rho + sqrt(2 * S1 * pmax(v1, 0) * L) + sqrt(2 * S2 * pmax(v2, 0) * L)
    }
    out[s] <- abs(W) - rho - const
  }
  out
}

# Sharper finite-sample interval (Bennett's inequality) -----------------
#
# At c = log Z, W(c) is a sum of N independent terms with |term - mean| <= 1 and
# total variance n_ov. Bennett's inequality bounds |W(log Z)| by the root rho of
# v h(rho / v) = log(2 / d_W), h(u) = (1 + u) log(1 + u) - u, for any v >= n_ov;
# v+(c) = min(A, B) is such a v with probability >= 1 - d_A - d_B, where
#   A = sum_k (sqrt(S_k Vhat_k(c)) + sqrt(2 S_k log(2 / d_A) / (S_k - 1)))^2  [Maurer-Pontil, Thm 10]
#   B = ((sqrt(2 L_B) + sqrt(2 L_B + 4 M(c))) / 2)^2,  M(c) = sum_j pi_c(theta*_j), L_B = log(1 / d_B)
#       [n_ov <= N H, and M(log Z) is a sum of [0, 1] terms with mean N H: Chernoff lower tail].
# The interval is the hull of {c : |W(c)| <= rho(c)}; d_W = 0.8 delta, d_A = d_B = 0.1 delta.
.bennett_radius <- function(v, L) {           # vectorized: root x >= 0 of v h(x / v) = L
  lo <- rep(0, length(v)); hi <- sqrt(2 * v * L) + 2 * L     # Bernstein's root bounds Bennett's
  h <- function(u) (1 + u) * log1p(u) - u
  for (it in 1:60) { mid <- (lo + hi) / 2; f <- v * h(mid / v) - L; lo <- ifelse(f < 0, mid, lo); hi <- ifelse(f < 0, hi, mid) }
  hi
}
.fs_F <- function(c, a_post, a_prop, lr, delta, fold_post = NULL, fold_prop = NULL) {
  K <- if (is.null(fold_post)) 1L else max(fold_post)
  fp <- if (K == 1L) rep(1L, length(a_post)) else fold_post; fq <- if (K == 1L) rep(1L, length(a_prop)) else fold_prop
  dW <- 0.8 * delta / K; dA <- 0.1 * delta / K; dB <- 0.1 * delta / K; LB <- log(1 / dB)
  out <- numeric(length(c))
  for (s in split(seq_along(c), ceiling(seq_along(c) / 50))) {           # blocks of 50 grid points
    P <- 1 / (1 + exp(outer(lr + a_post, c[s], `-`)))                      # posterior responsibilities, N1 x G
    Q <- 1 - 1 / (1 + exp(outer(lr + a_prop, c[s], `-`)))                  # 1 - responsibilities of proposal draws
    W <- colSums(Q) - colSums(P); rho <- 0
    for (k in seq_len(K)) {
      i1 <- fp == k; i2 <- fq == k; S1 <- sum(i1); S2 <- sum(i2)
      v1 <- (colSums(P[i1, , drop = FALSE]^2) - colSums(P[i1, , drop = FALSE])^2 / S1) * S1 / (S1 - 1)   # S1 times the unbiased variance
      v2 <- (colSums(Q[i2, , drop = FALSE]^2) - colSums(Q[i2, , drop = FALSE])^2 / S2) * S2 / (S2 - 1)
      A <- (sqrt(pmax(v1, 0)) + sqrt(2 * S1 * log(2 / dA) / (S1 - 1)))^2 + (sqrt(pmax(v2, 0)) + sqrt(2 * S2 * log(2 / dA) / (S2 - 1)))^2
      B <- ((sqrt(2 * LB) + sqrt(2 * LB + 4 * colSums(P[i1, , drop = FALSE]))) / 2)^2
      rho <- rho + .bennett_radius(pmin(A, B), log(2 / dW))
    }
    out[s] <- abs(W) - rho
  }
  out
}
# Smallest closed interval containing {c : F(c) <= 0}, scanning from far outside the
# range of the log ratios towards chat (F(chat) < 0 since W(chat) = 0).
.hull_interval <- function(F, chat, a_post, a_prop, lr, ngrid) {
  fin <- c(a_post[is.finite(a_post)], a_prop[is.finite(a_prop)])
  B <- log(length(a_post) + length(a_prop)) + 40
  side <- function(far) {
    if (F(far) <= 0) return(sign(far - chat) * Inf)
    grid <- seq(far, chat, length.out = ngrid); Fg <- F(grid)
    k <- which(Fg <= 0)[1L]
    if (is.na(k)) return(chat)
    if (k == 1L) return(grid[1L])
    stats::uniroot(F, c(grid[k - 1L], grid[k]), tol = 1e-10)$root
  }
  c(lower = side(min(fin) + lr - B), upper = side(max(fin) + lr + B))
}
.fs_interval <- function(a_post, a_prop, lr, chat, delta, ngrid = 150L, fold_post = NULL, fold_prop = NULL) {
  .hull_interval(function(cc) .fs_F(cc, a_post, a_prop, lr, delta, fold_post, fold_prop), chat, a_post, a_prop, lr, ngrid)
}

# Finite-sample interval for Markov chain output (Proposition S-fsmcmc) --------
#
# For stationary reversible chains whose right spectral gap is at least gamma, write
# W(log Z) = U - V, with U the centred proposal sum and V the centred chain sum. U is
# sub-gamma with variance factor N2 Var_g(pi*) and scale 1/3, V with variance factor
# c0 N1 Var_p(pi*) and scale 1 / gamma (Lemma S-lazy, b(gamma) <= 1 / gamma), where
# c0 = (2 - gamma) / gamma. The two variances are bounded from the responsibilities with
# probability 1 - 3 eps, eps = 0.05 delta:
#   N2 Var_g  <= min(A2, B): A2 empirical Bernstein (Maurer-Pontil, Thm 10);
#   N1 Var_p  <= min(A1, B): A1 from a chain-Bernstein lower tail for the sum of squared
#                deviations of the posterior responsibilities from a0 = M2 / N1, which is
#                computed from the proposal draws and so independent of the chain;
#                B from the overlap counts.
# u(m; a, b, L) is the largest x with x <= m + sqrt(2 a x L) + b L. With probability
# 1 - d_W, d_W = 0.85 delta, |W(log Z)| <= rho(c) = sqrt(2 v(c) L_W) + L_W / gamma,
# v = min(A2, B) + c0 min(A1, B), L_W = log(2 / d_W). The interval is the hull of
# {c : |W(c)| <= rho(c)}.
.fs_F_mcmc <- function(c, a_post, a_prop, lr, delta, gamma) {
  c0 <- (2 - gamma) / gamma; b <- 1 / gamma; S1 <- length(a_post); S2 <- length(a_prop)
  LW <- log(2 / (0.85 * delta)); eps <- 0.05 * delta; l1 <- log(1 / eps); l2 <- log(2 / eps)
  u <- function(m, a, bb, L) ((sqrt(2 * a * L) + sqrt(2 * a * L + 4 * m + 4 * bb * L)) / 2)^2
  out <- numeric(length(c))
  for (s in split(seq_along(c), ceiling(seq_along(c) / 50))) {
    P <- 1 / (1 + exp(outer(lr + a_post, c[s], `-`)))                      # responsibilities of the posterior draws
    Pq <- 1 / (1 + exp(outer(lr + a_prop, c[s], `-`)))                     # responsibilities of the proposal draws
    M1 <- colSums(P); M2 <- S2 - colSums(Pq); W <- M2 - M1
    B <- pmin(u(M2, 1, 0, l2), u(M1, c0, b, l2))                           # N H, from either overlap count
    V2 <- (colSums(Pq^2) - colSums(Pq)^2 / S2) / (S2 - 1)
    A2 <- S2 * (sqrt(pmax(V2, 0)) + sqrt(2 * l1 / (S2 - 1)))^2
    a0 <- pmin(pmax(M2 / S1, 0), 1)                                        # E_p pi estimated from the proposal draws:
    Q <- pmax(colSums(P^2) - 2 * a0 * M1 + S1 * a0^2, 0)                   # independent of the chain, so sum_j (pi_j - a0)^2
    mq <- pmax(a0, 1 - a0)^2                                               # (pi - a0)^2 <= mq, so its variance is <= mq E(.)
    A1 <- u(Q, c0 * mq, b * mq, l1)                                        # bounds N1 E_p (pi - a0)^2 >= N1 Var_p(pi)
    v <- pmin(A2, B) + c0 * pmin(A1, B)
    out[s] <- abs(W) - (sqrt(2 * v * LW) + LW / gamma)
  }
  out
}

#' Finite-sample confidence interval for log Z from Markov chain output
#'
#' For posterior draws from one or several independent, stationary, reversible Markov
#' chains whose spectral gap is at least `gamma`, inverts a Bernstein inequality for the
#' estimating function with a variance bound computed from the responsibilities. The
#' interval covers \eqn{\log Z} with probability at least \eqn{1-\delta} for every number
#' of draws, without thinning. It is valid for every chain whose spectral gap is at least
#' `gamma`, so a smaller `gamma` gives a wider interval that is still valid. With
#' `gamma = 1` it is a more conservative alternative to [bridge_fs_interval()].
#'
#' @inheritParams bridge_logml
#' @param gamma Lower bound on the spectral gap of the posterior sampler, in (0, 1].
#' @param delta Error probability.
#' @param ngrid Number of grid points used to bracket the end points before they are
#'   refined by root finding.
#' @return Named numeric vector with elements `lower` and `upper`.
#' @export
bridge_fs_interval_mcmc <- function(a_post, a_prop, gamma, delta = 0.05, ngrid = 150L) {
  .check_scalar(delta, "delta", 0, 1)
  if (!is.numeric(gamma) || length(gamma) != 1L || !is.finite(gamma) || gamma <= 0 || gamma > 1)
    stop("`gamma` must be a number in (0, 1].", call. = FALSE)
  a_post <- .clean_log_ratios(a_post, "posterior"); a_prop <- .clean_log_ratios(a_prop, "proposal")
  lr <- log(length(a_post) / length(a_prop))
  chat <- .bridge_root(a_post, a_prop, lr)$c
  if (!is.finite(chat)) return(c(lower = -Inf, upper = Inf))
  .hull_interval(function(cc) .fs_F_mcmc(cc, a_post, a_prop, lr, delta, gamma), chat, a_post, a_prop, lr, ngrid)
}

#' Sharper finite-sample confidence interval for log Z (independent or cross-fitted draws)
#'
#' Inverts Bennett's inequality for the estimating function at each candidate value,
#' with a variance bound valid with high probability (the smaller of an empirical
#' bound and the overlap count). Covers \eqn{\log Z} with probability at least
#' \eqn{1-\delta} for every sample size and every target/proposal pair under
#' independent sampling; about 1.8 times the delta-method interval at
#' \eqn{n_{ov}\approx50} and 1.35 times at \eqn{n_{ov}\approx500}.
#'
#' @inheritParams bridge_logml
#' @param delta Error probability.
#' @param ngrid Number of grid points used to bracket the end points before they
#'   are refined by root finding.
#' @param folds_post,folds_prop Optional fold labels (cross-fitting); radii are computed
#'   per fold at level `delta / K` and added.
#' @return Named numeric vector with elements `lower` and `upper`.
#' @export
bridge_fs_interval <- function(a_post, a_prop, delta = 0.05, ngrid = 150L, folds_post = NULL, folds_prop = NULL) {
  .check_scalar(delta, "delta", 0, 1)
  a_post <- .clean_log_ratios(a_post, "posterior"); a_prop <- .clean_log_ratios(a_prop, "proposal")
  lr <- log(length(a_post) / length(a_prop))
  chat <- .bridge_root(a_post, a_prop, lr)$c
  if (!is.finite(chat)) return(c(lower = -Inf, upper = Inf))
  .fs_interval(a_post, a_prop, lr, chat, delta, ngrid, folds_post, folds_prop)
}
