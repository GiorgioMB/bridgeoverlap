# Exact computations for the two-point family (Theorems 2 and 7) ----------------------
#
# Theta = {0, 1, 2}, p = (1 - eps, eps, 0), g = (0, eps, 1 - eps), S1 = S2 = S.
# Only draws at the overlap point theta = 1 are ambiguous. With n_post and
# n_prop draws there, Delta = log(n_prop / n_post) when both are positive,
# Z-hat = 0 (Delta = -Inf) when n_prop = 0 and Z-hat = Inf when n_prop > 0 = n_post.
# H = eps / 4, beta = 1 - eps, n_ov = S eps (1 - eps) / 2.

two_point_constants <- function(S, eps) {
  H <- eps / 4
  beta <- 1 - eps
  list(H = H, beta = beta, n_ov = S * eps * (1 - eps) / 2,
       sigma_N = sqrt(beta / (2 * S * H)))
}

.support <- function(p, tail) {
  cp <- cumsum(p)
  lo <- max(1L, which(cp >= tail)[1L] - 1L)            # index of count, 0-based shift below
  hi <- which(cp >= 1 - tail)[1L] - 1L
  if (is.na(hi)) hi <- length(p) - 1L
  seq.int(lo, max(lo, hi))
}

# Exact law of Delta, the delta-method MCSE and their summaries.
# pk: pmf of n_prop on 0..S, pm: pmf of n_post on 0..S.
two_point_law <- function(S, pk, pm, tail = 1e-11, max_pairs = 4e6) {
  ks <- .support(pk, tail)
  ms <- .support(pm, tail)
  p_minus <- pk[1L]                         # n_prop = 0
  p_plus <- (1 - pk[1L]) * pm[1L]           # n_prop > 0, n_post = 0
  if (length(ks) * length(ms) <= max_pairs) {
    k <- rep(ks, times = length(ms))
    m <- rep(ms, each = length(ks))
    list(k = k, m = m, w = pk[k + 1L] * pm[m + 1L], Delta = log(k / m),
         p_minus = p_minus, p_plus = p_plus, ks = ks, ms = ms)
  } else {
    list(k = NULL, ks = ks, ms = ms, p_minus = p_minus, p_plus = p_plus)
  }
}

# Probability that |Delta| <= width (exact; never covered at +-Inf).
two_point_cover_fixed <- function(law, pk, pm, width) {
  Fk <- cumsum(pk)                          # Fk[j + 1] = P(n_prop <= j)
  m <- law$ms
  hi <- pmin(floor(m * exp(width) + 1e-9), length(pk) - 1L)
  lo <- pmax(ceiling(m * exp(-width) - 1e-9), 1L)
  pr <- ifelse(hi >= lo, Fk[hi + 1L] - Fk[lo], 0)
  sum(pm[m + 1L] * pr)
}

# Kolmogorov distance of Delta / scale to N(0, 1).
two_point_ks <- function(law, pk, pm, scale, ngrid = 8001L) {
  if (!is.null(law$k)) {
    return(ks_discrete(law$Delta / scale, law$w, law$p_minus, law$p_plus))
  }
  # Dense laws: evaluate the CDF on a fine grid (atoms are then tiny).
  Fk <- cumsum(pk)
  m <- law$ms
  wm <- pm[m + 1L]
  x <- seq(-7, 7, length.out = ngrid)
  Fx <- vapply(x, function(xx) {
    j <- pmin(floor(m * exp(xx * scale) + 1e-9), length(pk) - 1L)
    law$p_minus + sum(wm * (Fk[j + 1L] - Fk[1L]))
  }, numeric(1))
  max(abs(Fx - stats::pnorm(x)), law$p_minus, law$p_plus)
}

# Independent draws: coverage of the delta-method interval with the estimated
# MCSE, Kolmogorov distance of Delta / sigma_N and median |log(sigma-hat/sigma_N)|.
two_point_iid <- function(S, eps, z = stats::qnorm(0.975)) {
  cst <- two_point_constants(S, eps)
  pk <- stats::dbinom(0:S, S, eps)
  law <- two_point_law(S, pk, pk)
  k <- law$k; m <- law$m; w <- law$w
  q <- k / (k + m)                          # fitted responsibility at the overlap point
  n_ov_hat <- ((1 - q)^2 * k * (S - k) + q^2 * m * (S - m)) / (S - 1)
  sig_hat <- sqrt(n_ov_hat) / (m * q)
  cover <- sum(w[abs(law$Delta) <= z * sig_hat])
  dev <- abs(log(sig_hat / cst$sigma_N))
  med_dev <- weighted_median(c(dev, Inf), c(w, law$p_minus + law$p_plus))
  data.frame(family = "two-point", eps = eps, S = S, n_ov = cst$n_ov,
             coverage = cover, ks = ks_discrete(law$Delta / cst$sigma_N, w, law$p_minus, law$p_plus),
             med_abs_log_mcse = med_dev, p_zero = law$p_minus)
}

# Number of visits to state 1 in S steps of the stationary two-state chain with
# P(0,1) = eps * gamma and P(1,0) = (1 - eps) * gamma, by evaluating the
# probability generating function on the unit circle (2 x 2 matrix powers by
# repeated squaring) and inverting it with the FFT.
markov_count_pmf <- function(S, eps, gamma) {
  L <- stats::nextn(S + 1L, factors = 2L)
  z <- exp(2i * pi * (0:(L - 1L)) / L)
  p00 <- 1 - eps * gamma; p01 <- eps * gamma
  p10 <- (1 - eps) * gamma; p11 <- 1 - (1 - eps) * gamma
  # A(z) = P %*% diag(1, z)
  a <- rep(p00 + 0i, L); b <- p01 * z; c <- rep(p10 + 0i, L); d <- p11 * z
  Ra <- rep(1 + 0i, L); Rb <- rep(0i, L); Rc <- rep(0i, L); Rd <- rep(1 + 0i, L)
  e <- S - 1L
  while (e > 0L) {
    if (e %% 2L == 1L) {
      na <- Ra * a + Rb * c; nb <- Ra * b + Rb * d
      nc <- Rc * a + Rd * c; nd <- Rc * b + Rd * d
      Ra <- na; Rb <- nb; Rc <- nc; Rd <- nd
    }
    na <- a * a + b * c; nb <- a * b + b * d
    nc <- c * a + d * c; nd <- c * b + d * d
    a <- na; b <- nb; c <- nc; d <- nd
    e <- e %/% 2L
  }
  # G(z) = (1 - eps, eps z) R(z) 1
  G <- (1 - eps) * (Ra + Rb) + eps * z * (Rc + Rd)
  p <- Re(stats::fft(G)) / L
  p <- pmax(p[seq_len(S + 1L)], 0)
  p / sum(p)
}

# MCMC: coverage of intervals based on the true sigma_NP and on sigma_N
# (which ignores autocorrelation), and Kolmogorov distance of Delta / sigma_NP.
two_state_chain <- function(S, eps, gamma, z = stats::qnorm(0.975)) {
  cst <- two_point_constants(S, eps)
  tau_W <- 1 / gamma                        # Theorem 7
  sigma_NP <- sqrt(tau_W) * cst$sigma_N
  pk <- stats::dbinom(0:S, S, eps)
  pm <- markov_count_pmf(S, eps, gamma)
  law <- two_point_law(S, pk, pm)
  data.frame(family = "two-state chain", eps = eps, gamma = gamma, S = S,
             n_ov = cst$n_ov, n_ovP = gamma * cst$n_ov,
             coverage_sigmaNP = two_point_cover_fixed(law, pk, pm, z * sigma_NP),
             coverage_sigmaN = two_point_cover_fixed(law, pk, pm, z * cst$sigma_N),
             ks = two_point_ks(law, pk, pm, sigma_NP))
}
