# Numerical core ---------------------------------------------------------------
#
# For a draw theta with unnormalized posterior
# f(theta) and proposal density g(theta), a = log f(theta) - log g(theta) is
# its log importance ratio. With c = log z, s1 = S1/N, s2 = S2/N and
# lr = log(s1/s2), the responsibility of the draw is
#
#   pi_c(a) = s2 z g / (s1 f + s2 z g) = 1 / (1 + exp(lr + a - c)),
#
# the probability, under the pooled law, that the draw came from the proposal.
# The estimating function is
#
#   W(c) = sum_i {1 - pi_c(a_i)} - sum_j pi_c(a_j),
#
# summing over proposal draws i and posterior draws j. Its root is the optimal
# bridge sampling estimate of log Z. All computations are on the log
# scale and use plogis()/dlogis(), which are exact at +-Inf.

.resp <- function(a, c, lr) stats::plogis(c - a - lr)          # pi_c(a)
.one_minus_resp <- function(a, c, lr) stats::plogis(a + lr - c)  # 1 - pi_c(a)
.resp_var <- function(a, c, lr) stats::dlogis(c - a - lr)      # pi (1 - pi)

.clean_log_ratios <- function(a, what) {
  if (!is.numeric(a)) {
    stop(sprintf("The %s log ratios must be numeric.", what), call. = FALSE)
  }
  a <- as.numeric(a)
  bad <- is.na(a)
  if (any(bad)) {
    warning(sprintf(
      "%d %s log ratio(s) were NA/NaN and have been set to -Inf.",
      sum(bad), what), call. = FALSE)
    a[bad] <- -Inf
  }
  a
}

# Safeguarded Newton iteration for the root of W, which is continuous and
# non-increasing. Draws with infinite log ratios contribute
# constants: a proposal draw with a = -Inf (f = 0) contributes 0, one with
# a = +Inf contributes 1; a posterior draw with a = +Inf (g = 0) contributes 0,
# one with a = -Inf contributes -1. The estimate is
# sup{c : W(c) > 0}, which is -Inf or +Inf when the samples do not overlap.
.bridge_root <- function(a_post, a_prop, lr, c0 = NULL, tol = 1e-12,
                         maxit = 500L) {
  fin_post <- is.finite(a_post)
  fin_prop <- is.finite(a_prop)
  K <- sum(a_prop == Inf) - sum(a_post == -Inf)
  ap <- a_post[fin_post]
  aq <- a_prop[fin_prop]
  W_lower <- K + length(aq)  # limit of W as c -> -Inf
  W_upper <- K - length(ap)  # limit of W as c -> +Inf
  if (W_lower <= 0) {
    return(list(c = -Inf, iterations = 0L, converged = TRUE,
                status = "no proposal draw in the overlap region (Z-hat = 0)"))
  }
  if (W_upper >= 0) {
    return(list(c = Inf, iterations = 0L, converged = TRUE,
                status = "no posterior draw in the overlap region (Z-hat = Inf)"))
  }
  W <- function(c) K + sum(stats::plogis(aq + lr - c)) - sum(stats::plogis(c - ap - lr))
  dW <- function(c) -(sum(stats::dlogis(c - ap - lr)) + sum(stats::dlogis(c - aq - lr)))

  # Bracket: beyond B = log(N) + 20 log-units outside the data every
  # responsibility is within N * exp(-B) of 0 or 1, so W has the sign of its
  # (integer-valued) limit there.
  rng <- range(c(ap, aq)) + lr
  B <- log(length(ap) + length(aq)) + 20
  lo <- rng[1] - B
  hi <- rng[2] + B
  x <- if (!is.null(c0) && is.finite(c0)) c0 else stats::median(if (length(ap)) ap else aq)
  x <- min(max(x, lo), hi)
  Wx <- W(x)
  it <- 0L
  converged <- FALSE
  while (it < maxit) {
    it <- it + 1L
    if (Wx > 0) {
      lo <- x
    } else if (Wx < 0) {
      hi <- x
    } else {
      converged <- TRUE
      break
    }
    d <- dW(x)
    xn <- if (is.finite(d) && d < 0) x - Wx / d else NA_real_
    if (!is.finite(xn) || xn <= lo || xn >= hi) xn <- 0.5 * (lo + hi)
    step <- abs(xn - x)
    x <- xn
    Wx <- W(x)
    scale <- max(1, abs(x))
    if (step <= tol * scale || (hi - lo) <= tol * scale) {
      converged <- TRUE
      break
    }
  }
  list(c = x, iterations = it, converged = converged,
       status = if (converged) "converged" else "maximum number of iterations reached")
}

#' Optimal bridge sampling estimate from log importance ratios
#'
#' Computes the optimal bridge sampling estimate of Meng and Wong (1996) as the
#' root of the estimating function
#' \deqn{W(c) = \sum_{i=1}^{S_2}\{1-\pi_c(\tilde\theta_i)\} - \sum_{j=1}^{S_1}\pi_c(\theta^*_j),}
#' where \eqn{\pi_c} is the responsibility (the fitted probability, in an
#' intercept-only logistic regression, that a draw came from the proposal).
#' The root is found by safeguarded Newton iterations on the log scale, which
#' converge from any starting value. This is the same estimator as the
#' fixed-point iteration of \pkg{bridgesampling} (with `use_neff = FALSE`), but
#' it is exact at the solution and handles non-overlapping samples.
#'
#' @param a_post Numeric vector of log ratios \eqn{\log f - \log g} at the
#'   \eqn{S_1} posterior draws, where `f` is the unnormalized posterior and `g`
#'   the (normalized) proposal density.
#' @param a_prop Numeric vector of log ratios at the \eqn{S_2} proposal draws.
#' @param c0 Optional starting value on the log scale.
#' @param tol Relative tolerance on the root.
#' @param maxit Maximum number of iterations.
#' @return A list with the estimate `logml` (which is `-Inf` or `Inf` when the
#'   two samples do not overlap), the number of `iterations`, a logical
#'   `converged` and a character `status`.
#' @seealso [bridge_diagnostics()] for the estimate together with all diagnostics.
#' @examples
#' set.seed(1)
#' # target f = Z * N(0, 1) with log Z = 2, proposal g = N(1, 1.5^2)
#' x <- rnorm(2000); y <- rnorm(2000, 1, 1.5)
#' a <- function(t) 2 + dnorm(t, log = TRUE) - dnorm(t, 1, 1.5, log = TRUE)
#' bridge_logml(a(x), a(y))$logml
#' @export
bridge_logml <- function(a_post, a_prop, c0 = NULL, tol = 1e-12, maxit = 500L) {
  a_post <- .clean_log_ratios(a_post, "posterior")
  a_prop <- .clean_log_ratios(a_prop, "proposal")
  if (length(a_post) < 1L || length(a_prop) < 1L) {
    stop("Both samples must contain at least one draw.", call. = FALSE)
  }
  lr <- log(length(a_post) / length(a_prop))
  r <- .bridge_root(a_post, a_prop, lr, c0 = c0, tol = tol, maxit = maxit)
  list(logml = r$c, iterations = r$iterations, converged = r$converged,
       status = r$status)
}

#' Estimating function of the optimal bridge sampling estimator
#'
#' Evaluates \eqn{W(c)} on a grid of values of
#' \eqn{c = \log z}. `W` is non-increasing and its root is the optimal bridge
#' sampling estimate of \eqn{\log Z}.
#'
#' @inheritParams bridge_logml
#' @param c Numeric vector of values of \eqn{\log z}.
#' @return Numeric vector of the same length as `c`.
#' @export
bridge_estimating_function <- function(c, a_post, a_prop) {
  a_post <- .clean_log_ratios(a_post, "posterior")
  a_prop <- .clean_log_ratios(a_prop, "proposal")
  lr <- log(length(a_post) / length(a_prop))
  vapply(c, function(cc) {
    sum(.one_minus_resp(a_prop, cc, lr)) - sum(.resp(a_post, cc, lr))
  }, numeric(1))
}
