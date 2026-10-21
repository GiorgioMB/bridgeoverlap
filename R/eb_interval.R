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

.eb_interval <- function(a_post, a_prop, lr, chat, delta, ngrid = 200L,
                         fold_post = NULL, fold_prop = NULL) {
  K <- if (is.null(fold_post)) 1L else max(fold_post)
  L <- log(8 * K / delta)
  sizes1 <- if (K == 1L) length(a_post) else tabulate(fold_post, K)
  sizes2 <- if (K == 1L) length(a_prop) else tabulate(fold_prop, K)
  const <- 7 * L / 3 * sum(sizes1 / (sizes1 - 1) + sizes2 / (sizes2 - 1))
  # rho(c) <= rho_max because the unbiased variance of a [0, 1] variable is at
  # most S / (4 (S - 1)); outside {|W| <= rho_max} the set is empty.
  rho_max <- sum(sqrt(2 * sizes1^2 / (4 * (sizes1 - 1)) * L)) +
    sum(sqrt(2 * sizes2^2 / (4 * (sizes2 - 1)) * L)) + const
  fin <- c(a_post[is.finite(a_post)], a_prop[is.finite(a_prop)])
  B <- log(length(a_post) + length(a_prop)) + 40
  c_far_lo <- min(fin) + lr - B
  c_far_hi <- max(fin) + lr + B
  F <- function(cc) .eb_F(cc, a_post, a_prop, lr, L, const, fold_post, fold_prop)
  Wf <- function(cc) bridge_estimating_function(cc, a_post, a_prop)

  side <- function(far, dir) {
    # dir = -1: search below chat; dir = +1: above.
    if (F(far) <= 0) return(dir * Inf)
    target <- if (dir < 0) rho_max else -rho_max
    # point beyond which |W| > rho_max (monotone W): bracket for the grid scan
    edge <- if ((Wf(far) - target) * (Wf(chat) - target) < 0) {
      stats::uniroot(function(cc) Wf(cc) - target, sort(c(far, chat)),
                     tol = 1e-10)$root
    } else {
      far
    }
    grid <- seq(edge, chat, length.out = ngrid)
    Fg <- F(grid)
    k <- which(Fg <= 0)[1L]
    if (is.na(k)) return(chat)
    if (k == 1L) return(grid[1L])
    stats::uniroot(F, c(grid[k - 1L], grid[k]), tol = 1e-10)$root
  }
  lower <- side(c_far_lo, -1)
  upper <- side(c_far_hi, +1)
  c(lower = lower, upper = upper)
}

#' Finite-sample confidence interval for log Z (independent draws)
#'
#' Computes the finite-sample interval by inverting an
#' empirical Bernstein bound for the estimating function. Under independent
#' sampling it covers \eqn{\log Z} with probability at least \eqn{1-\delta}
#' for every sample size and every target/proposal pair; it is unbounded when
#' the two samples barely overlap. It is not valid for Markov chain output or
#' when the proposal is fitted on the same draws.
#'
#' @inheritParams bridge_logml
#' @param delta Error probability.
#' @param ngrid Number of grid points used to locate the end points before
#'   refining them with [stats::uniroot()].
#' @return Named numeric vector with elements `lower` and `upper`.
#' @export
bridge_eb_interval <- function(a_post, a_prop, delta = 0.05, ngrid = 200L) {
  .check_scalar(delta, "delta", 0, 1)
  a_post <- .clean_log_ratios(a_post, "posterior")
  a_prop <- .clean_log_ratios(a_prop, "proposal")
  lr <- log(length(a_post) / length(a_prop))
  chat <- .bridge_root(a_post, a_prop, lr)$c
  if (!is.finite(chat)) return(c(lower = -Inf, upper = Inf))
  .eb_interval(a_post, a_prop, lr, chat, delta, ngrid)
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
.fs_interval <- function(a_post, a_prop, lr, chat, delta, ngrid = 150L, fold_post = NULL, fold_prop = NULL) {
  F <- function(cc) .fs_F(cc, a_post, a_prop, lr, delta, fold_post, fold_prop)
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

#' Sharper finite-sample confidence interval for log Z (independent or cross-fitted draws)
#'
#' Inverts Bennett's inequality for the estimating function at each candidate value,
#' with a variance bound valid with high probability (the smaller of an empirical
#' bound and the overlap count). Covers \eqn{\log Z} with probability at least
#' \eqn{1-\delta} for every sample size and every target/proposal pair under
#' independent sampling; about 1.8 times the delta-method interval at
#' \eqn{n_{ov}\approx50} and 1.35 times at \eqn{n_{ov}\approx500}.
#'
#' @inheritParams bridge_eb_interval
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
