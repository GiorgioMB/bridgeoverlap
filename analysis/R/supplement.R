# Helpers for the supplementary tables ------------------------------------------------

z975 <- stats::qnorm(0.975)
lse <- function(v) { m <- max(v); m + log(sum(exp(v - m))) }
lse2 <- function(a, b) pmax(a, b) + log1p(exp(-abs(a - b)))

# Rejection indicators of the failure diagnostics at a nominal level (Table E2):
# I_max with the corrected Bonferroni threshold Phi^{-1}(1 - alpha/2G)/sqrt(G - 1),
# the naive one with sqrt(G), and the simulated pivotal quantile; G_max with the
# Grubbs quantile; t_T with t_{G-1}; V_jack with chi^2_{G-1} / (G - 1).
rejections <- function(d, alpha = 0.05) {
  G <- d$settings$G
  q <- jackknife_null_quantiles(G, probs = 1 - alpha)
  zb <- stats::qnorm(1 - alpha / (2 * G))
  c(V_jack = isTRUE(d$V_jack > stats::qchisq(1 - alpha, G - 1) / (G - 1)),
    I_max_bonf = isTRUE(d$I_max > zb / sqrt(G - 1)),
    I_max_naive = isTRUE(d$I_max > zb / sqrt(G)),
    I_max_exact = isTRUE(d$I_max > q["I_max", 1]),
    G_max = isTRUE(d$G_max > q["G_max", 1]),
    t_T = isTRUE(abs(d$t_T) > stats::qt(1 - alpha / 2, G - 1)))
}

# Distributions ------------------------------------------------------------------------
lmvn_std <- function(X) -0.5 * rowSums(X^2) - ncol(X) / 2 * log(2 * pi)
rmvt_std <- function(n, d, nu) matrix(stats::rnorm(n * d), n, d) / sqrt(stats::rchisq(n, nu) / nu)
lmvt_std <- function(X, nu) {
  d <- ncol(X)
  lgamma((nu + d) / 2) - lgamma(nu / 2) - d / 2 * log(nu * pi) - (nu + d) / 2 * log1p(rowSums(X^2) / nu)
}
# diagonal Gaussian N(m, diag(v)) and full-covariance Gaussian fitted by moments
ldiag <- function(X, m, v) {
  rowSums(stats::dnorm(X, matrix(m, nrow(X), length(m), byrow = TRUE),
                       matrix(sqrt(v), nrow(X), length(v), byrow = TRUE), log = TRUE))
}
fit_gauss <- function(X) {
  m <- colMeans(X)
  S <- stats::cov(X)
  list(mean = m, cov = S, L = t(chol(S)))
}
rgauss <- function(n, fit) {
  d <- length(fit$mean)
  sweep(matrix(stats::rnorm(n * d), n, d) %*% t(fit$L), 2, fit$mean, "+")
}
lgauss <- function(X, fit) mvtnorm::dmvnorm(X, fit$mean, fit$cov, log = TRUE)

# Population overlap of two one-dimensional densities (log densities lp, lg)
# for sample sizes S1, S2, by quadrature split at `breaks`.
overlap_1d <- function(lp, lg, S1, S2, breaks) {
  N <- S1 + S2; s1 <- S1 / N; s2 <- S2 / N; lr0 <- log(s1 / s2)
  pistar <- function(x) stats::plogis(lg(x) - lp(x) - lr0)
  E <- function(ld, h) integrate_pieces(function(x) h(x) * exp(ld(x)), breaks, -Inf)
  Ep1 <- E(lp, pistar); Ep2 <- E(lp, function(x) pistar(x)^2)
  Eg1 <- E(lg, pistar); Eg2 <- E(lg, function(x) pistar(x)^2)
  H <- s1 * Ep1
  I <- H / (s1 * s2)
  list(H = H, NH = N * H, I = I, beta = 1 - I,
       n_ov = S2 * (Eg2 - Eg1^2) + S1 * (Ep2 - Ep1^2), sigma_N = sqrt((1 - I) / (N * H)),
       pistar = pistar, s1 = s1, s2 = s2, lp = lp, lg = lg, breaks = breaks, N = N)
}

# Moments under the overlap-weighted law m_*(d theta) = pi*(1 - pi*) m(d theta) / H.
mstar_moments <- function(ov, xi = NULL) {
  w <- function(x) {
    p <- ov$pistar(x)
    p * (1 - p) * (ov$s1 * exp(ov$lp(x)) + ov$s2 * exp(ov$lg(x)))
  }
  I <- function(h) integrate_pieces(function(x) h(x) * w(x), ov$breaks, -Inf)
  H <- I(function(x) 1)
  e1 <- I(ov$pistar) / H
  e2 <- I(function(x) ov$pistar(x)^2) / H
  out <- list(H = H, mean = e1, var = e2 - e1^2)
  if (!is.null(xi)) {
    out$E_xi <- I(xi) / H
    out$cov_xi <- I(function(x) ov$pistar(x) * xi(x)) / H - e1 * out$E_xi
  }
  out
}

# Plug-in Lyapunov ratio of Theorem 3(b): the largest ratio over |t| <= tbar,
# with c* replaced by c-hat and expectations by sample averages.
lyapunov_plugin <- function(a_post, a_prop, chat, tbar = 0.2312, ngrid = 25L) {
  S1 <- length(a_post); S2 <- length(a_prop); lr <- log(S1 / S2)
  max(vapply(seq(-tbar, tbar, length.out = ngrid), function(t) {
    pp <- stats::plogis(chat + t - a_post - lr)
    pq <- stats::plogis(chat + t - a_prop - lr)
    m3 <- S2 * mean(abs(pq - mean(pq))^3) + S1 * mean(abs(pp - mean(pp))^3)
    v <- S2 * mean((pq - mean(pq))^2) + S1 * mean((pp - mean(pp))^2)
    m3 / v^1.5
  }, numeric(1)))
}

kappa_fun <- function(t) (1 - exp(-t))^2 / (2 * cosh(t / 2)^2 + 2 / 3 * (1 - exp(-t)))

# Monte Carlo standard error of a standard deviation (normal approximation).
se_sd <- function(x) stats::sd(x) / sqrt(2 * (length(x) - 1))

# Output --------------------------------------------------------------------------------
fnum <- function(x, digits) ifelse(is.na(x), "--", formatC(x, digits = digits, format = "f"))
fval_se <- function(x, se, digits) ifelse(is.na(x), "--", sprintf("%s (%s)", fnum(x, digits), fnum(se, digits)))

supplement_root <- function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(dirname(normalizePath(f))) else if (dir.exists("analysis")) "analysis" else ".."
}

# Stationary AR(1) chain in d dimensions with N(0, I_d) marginal (n x d matrix).
ar1_matrix <- function(n, d, rho) {
  E <- matrix(stats::rnorm(n * d), n, d)
  if (rho > 0) {
    E[-1L, ] <- sqrt(1 - rho^2) * E[-1L, ]
    E <- unclass(stats::filter(E, rho, method = "recursive"))
  }
  matrix(E, n, d)
}
