#' bridgeoverlap: overlap certificates and failure diagnostics for bridge sampling
#'
#' The optimal bridge sampling estimator of Meng and Wong (1996) is the
#' maximum-likelihood estimate of the intercept of a logistic regression that
#' separates posterior draws from proposal draws. Its estimating function is a
#' sum of bounded terms, so its accuracy and the accuracy of its Monte Carlo
#' standard error are governed by the *overlap* of the two samples rather than
#' by moments of the importance ratios. This package reports the estimate, its MCSE, the effective-overlap certificate, a
#' finite-sample interval, a calibration test for sampler failure and
#' delete-a-group jackknife diagnostics, all computed from quantities that
#' [bridgesampling::bridge_sampler()] already produces.
#'
#' Main entry points:
#' * [bridge_diagnostics()] -- the estimate and every diagnostic, for `bridge` objects or raw log ratios.
#' * [bridge_logml()] -- the estimate alone (safeguarded Newton on the log scale).
#' * [bridge_sampler_cmdstanr()] -- run `bridgesampling` on a `cmdstanr` fit.
#' * [bridge_crossfit()] -- two-fold cross-fitted bridge sampling.
#'
#' @references
#' Meng, X.-L. and Wong, W. H. (1996). Simulating ratios of normalizing
#' constants via a simple identity: a theoretical exploration.
#' *Statistica Sinica*, 6, 831--860.
#'
#' Micaletto, G. and Vehtari, A. (2025). Bridge sampling diagnostics.
#' arXiv:2508.14487.
#' @keywords internal
"_PACKAGE"

`%||%` <- function(x, y) if (is.null(x)) y else x

# Evaluate `expr` with a fixed seed without disturbing the caller's RNG stream.
.with_seed <- function(seed, expr) {
  genv <- globalenv()
  had_seed <- exists(".Random.seed", envir = genv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = genv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = genv)
    } else if (exists(".Random.seed", envir = genv, inherits = FALSE)) {
      rm(".Random.seed", envir = genv)
    }
  }, add = TRUE)
  set.seed(seed)
  expr
}

.check_scalar <- function(x, name, lower = -Inf, upper = Inf) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < lower || x > upper) {
    stop(sprintf("`%s` must be a single number in [%s, %s].", name, lower, upper),
         call. = FALSE)
  }
  invisible(x)
}

.fmt <- function(x, digits = 4) {
  if (length(x) == 0L || is.null(x)) return("NA")
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "g", flag = "#"))
}
