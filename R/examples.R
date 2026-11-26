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

#' Path to a Stan program shipped with the package
#'
#' @param name File name without extension: one of
#'   `"eight_schools_centered_ml"`, `"eight_schools_noncentered_ml"`,
#'   `"normal_mixture_ml"`.
#' @return The full path to the `.stan` file.
#' @export
#' @examples
#' basename(bridgeoverlap_stan_file("eight_schools_centered_ml"))
bridgeoverlap_stan_file <- function(name) {
  f <- system.file("stan", paste0(name, ".stan"), package = "bridgeoverlap")
  if (!nzchar(f)) stop("Unknown Stan program: ", name, call. = FALSE)
  f
}
