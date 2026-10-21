#' Eight schools data
#'
#' The data of Rubin (1981) in the format of the `posteriordb` entries
#' `eight_schools_centered` and `eight_schools_noncentered`. With the priors
#' \eqn{\mu \sim N(0, 5^2)} and \eqn{\tau \sim} half-Cauchy(0, 5) the exact log
#' evidence is \eqn{\log Z = -31.3113}, see [eight_schools_log_evidence()].
#'
#' @return A list with elements `J`, `y` and `sigma`.
#' @export
eight_schools_data <- function() {
  list(J = 8L, y = c(28, 8, -3, 7, -1, 1, 18, 12),
       sigma = c(15, 10, 16, 11, 9, 11, 10, 18))
}

#' Exact log evidence of the eight schools model
#'
#' Integrates \eqn{\theta} and \eqn{\mu} analytically and the remaining
#' one-dimensional integral over \eqn{\log\tau} by adaptive quadrature.
#'
#' @param data A list as returned by [eight_schools_data()].
#' @param mu_sd Prior standard deviation of \eqn{\mu}.
#' @param tau_scale Scale of the half-Cauchy prior on \eqn{\tau}.
#' @return The log evidence (about -31.311347 for the default data).
#' @export
eight_schools_log_evidence <- function(data = eight_schools_data(), mu_sd = 5,
                                       tau_scale = 5) {
  log_marg <- function(tau) {
    S <- diag(data$sigma^2 + tau^2) + mu_sd^2
    R <- chol(S)
    z <- backsolve(R, data$y, transpose = TRUE)
    -0.5 * sum(z^2) - sum(log(diag(R))) - 0.5 * length(data$y) * log(2 * pi)
  }
  f <- function(lt, shift) {
    vapply(lt, function(l) {
      tau <- exp(l)
      exp(log_marg(tau) + stats::dcauchy(tau, 0, tau_scale, log = TRUE) + log(2) + l - shift)
    }, numeric(1))
  }
  grid <- seq(-15, 10, length.out = 2001)
  shift <- max(log(f(grid, 0)))
  val <- stats::integrate(f, -40, 15, shift = shift, rel.tol = 1e-12,
                          subdivisions = 1000L)$value
  log(val) + shift
}

#' Exact log evidence of the conjugate normal-inverse-gamma regression
#'
#' Closed-form marginal likelihood of the model in
#' `system.file("stan", "conjugate_regression_ml.stan", package = "bridgeoverlap")`:
#' \eqn{\beta \mid \sigma^2 \sim N(0, \sigma^2 v_0 I_K)},
#' \eqn{\sigma^2 \sim \mathrm{Inv\text{-}Gamma}(a_0, b_0)} and
#' \eqn{y \mid \beta, \sigma \sim N(X\beta, \sigma^2 I_n)}.
#'
#' @param X Design matrix (n x K).
#' @param y Response vector of length n.
#' @param v0 Prior variance scale of the coefficients.
#' @param a0,b0 Shape and scale of the inverse-gamma prior on \eqn{\sigma^2}.
#' @return The log marginal likelihood.
#' @export
conjugate_log_evidence <- function(X, y, v0, a0, b0) {
  X <- as.matrix(X)
  n <- nrow(X)
  K <- ncol(X)
  R <- chol(crossprod(X) + diag(1 / v0, K))
  b <- crossprod(X, y)
  beta_n <- backsolve(R, backsolve(R, b, transpose = TRUE))
  a_n <- a0 + n / 2
  b_n <- b0 + 0.5 * (sum(y^2) - sum(b * beta_n))
  -n / 2 * log(2 * pi) - sum(log(diag(R))) - K / 2 * log(v0) +
    a0 * log(b0) - a_n * log(b_n) + lgamma(a_n) - lgamma(a0)
}

#' Path to a Stan program shipped with the package
#'
#' @param name File name without extension: one of
#'   `"eight_schools_centered_ml"`, `"eight_schools_noncentered_ml"`,
#'   `"conjugate_regression_ml"` or `"gaussian_mixture_ml"`.
#' @return The full path to the `.stan` file.
#' @export
#' @examples
#' basename(bridgeoverlap_stan_file("eight_schools_centered_ml"))
bridgeoverlap_stan_file <- function(name) {
  f <- system.file("stan", paste0(name, ".stan"), package = "bridgeoverlap")
  if (!nzchar(f)) stop("Unknown Stan program: ", name, call. = FALSE)
  f
}
