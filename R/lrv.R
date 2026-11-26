# Long-run variance of the posterior responsibilities ---------------------------
#
# For Markov chain output the MCSE needs the long-run variance
# omega^2 = sigma_P^2(pi*) of the posterior responsibilities. As in Algorithm 1
# (step 4) we write omega^2 = V_post * tau, where V_post is the sample variance
# and tau the integrated autocorrelation time, and estimate tau
#   * "within":       Geyer's initial monotone sequence estimator applied within
#                     each chain and averaged over chains (the default);
#   * "multichain":   the multi-chain estimator of Vehtari et al. (2021) as
#                     implemented in posterior::ess_basic();
#   * "concat_geyer": Geyer's estimator on the concatenated chains;
#   * "concat_ar":    the autoregressive spectral estimator of coda on the
#                     concatenated chains, as used by bridgesampling's
#                     mcse_logml (shown for comparison only).
# Concatenating chains creates spurious transitions at the junctions, so the
# autocorrelation is better estimated within chains or with a multi-chain estimator.

# Biased (divisor n) autocovariances at lags 0, ..., n - 1, computed by FFT.
.autocov <- function(x) {
  n <- length(x)
  x <- x - mean(x)
  m <- stats::nextn(2L * n)
  fx <- stats::fft(c(x, numeric(m - n)))
  # as.numeric(): m * n overflows the integer range for long chains
  Re(stats::fft(Mod(fx)^2, inverse = TRUE))[seq_len(n)] / (as.numeric(m) * as.numeric(n))
}

# Geyer's (1992) initial monotone sequence estimator for a single sequence.
# Returns the long-run variance (asymptotic variance of sqrt(n) * mean) and the
# lag-0 autocovariance. As in posterior/Stan, the implied autocorrelation time
# is bounded below by 1/log10(n) to guard against antithetic artefacts.
.geyer_lrv <- function(x) {
  n <- length(x)
  if (n < 4L) return(c(lrv = NA_real_, gamma0 = NA_real_))
  g <- .autocov(x)
  g0 <- g[1L]
  if (!is.finite(g0) || g0 <= 0) return(c(lrv = 0, gamma0 = 0))
  npairs <- floor(n / 2)
  k <- seq_len(npairs)
  Gam <- g[2L * k - 1L] + g[2L * k]  # Gamma_k = gamma(2k) + gamma(2k + 1)
  neg <- which(Gam <= 0)
  m <- if (length(neg)) neg[1L] - 1L else npairs
  if (m < 1L) return(c(lrv = g0 / log10(max(n, 10)), gamma0 = g0))
  Gam <- cummin(Gam[seq_len(m)])     # initial monotone sequence
  lrv <- -g0 + 2 * sum(Gam)
  lrv <- max(lrv, g0 / log10(max(n, 10)))
  c(lrv = lrv, gamma0 = g0)
}

#' Integrated autocorrelation time of a (multi-chain) sequence
#'
#' Estimates the integrated autocorrelation time \eqn{\tau} of a scalar
#' sequence made of one or more concatenated chains. In
#' [bridge_diagnostics()] it is applied to the posterior responsibilities.
#'
#' @param x Numeric vector: the concatenation of the chains.
#' @param chain_lengths Integer vector of chain lengths summing to `length(x)`.
#' @param method One of `"within"` (Geyer's initial monotone sequence estimator
#'   within each chain, averaged over chains with chain-length weights; the
#'   default), `"multichain"`
#'   ([posterior::ess_basic()]), `"concat_geyer"` (Geyer's estimator on the
#'   concatenated sequence) or `"concat_ar"` (autoregressive spectral estimator
#'   [coda::spectrum0.ar()] on the concatenated sequence, as in the
#'   `mcse_logml` of \pkg{bridgesampling}). The last two ignore chain
#'   boundaries and are provided for comparison.
#' @return A single number, \eqn{\hat\tau}. It is `1` for a constant sequence.
#' @references Geyer, C. J. (1992). Practical Markov chain Monte Carlo.
#'   *Statistical Science*, 7, 473--483.
#'
#'   Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Bürkner, P.-C.
#'   (2021). Rank-normalization, folding, and localization: an improved
#'   \eqn{\hat R} for assessing convergence of MCMC. *Bayesian Analysis*, 16,
#'   667--718.
#' @examples
#' set.seed(1)
#' x <- as.numeric(arima.sim(list(ar = 0.8), n = 1e4))
#' iat(x)          # close to (1 + 0.8) / (1 - 0.8) = 9
#' @export
iat <- function(x, chain_lengths = length(x),
                method = c("within", "multichain", "concat_geyer", "concat_ar")) {
  method <- match.arg(method)
  x <- as.numeric(x)
  if (sum(chain_lengths) != length(x)) {
    stop("`chain_lengths` must sum to length(x).", call. = FALSE)
  }
  if (length(x) < 4L || stats::var(x) <= 0) return(1)
  ends <- cumsum(chain_lengths)
  starts <- ends - chain_lengths + 1L
  switch(method,
    within = {
      # one Geyer estimate per chain, averaged with chain-length weights; chains
      # whose responsibilities are constant (a stuck sampler) carry no
      # information about tau and are skipped (bridge_diagnostics() reports them)
      parts <- vapply(seq_along(chain_lengths), function(m) {
        .geyer_lrv(x[starts[m]:ends[m]])
      }, numeric(2))
      ok <- is.finite(parts[1, ]) & is.finite(parts[2, ]) & parts[2, ] > 0
      if (!any(ok)) return(1)
      tau_m <- parts[1, ok] / parts[2, ok]
      sum(chain_lengths[ok] * tau_m) / sum(chain_lengths[ok])
    },
    multichain = {
      if (length(unique(chain_lengths)) != 1L) {
        n <- min(chain_lengths)
        x_mat <- vapply(seq_along(chain_lengths), function(m) x[starts[m] - 1L + seq_len(n)],
                        numeric(n))
      } else {
        x_mat <- matrix(x, ncol = length(chain_lengths))
      }
      ess <- posterior::ess_basic(x_mat)
      if (!is.finite(ess) || ess <= 0) return(1)
      length(x_mat) / ess
    },
    concat_geyer = {
      p <- .geyer_lrv(x)
      if (!is.finite(p[["lrv"]]) || p[["gamma0"]] <= 0) 1 else p[["lrv"]] / p[["gamma0"]]
    },
    concat_ar = {
      spec <- coda::spectrum0.ar(x)$spec
      as.numeric(spec / stats::var(x))
    }
  )
}
