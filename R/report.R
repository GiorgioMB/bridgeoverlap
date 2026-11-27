# Worked printout and Bayes factors ------------------------------------------------

#' Tabular report of the diagnostics for one run
#'
#' Arranges the output of [bridge_diagnostics()] as a table: the estimate, the MCSE, the
#' two checks with their rules, the finite-sample interval and the decision.
#'
#' @param x A `bridge_diagnostics` object.
#' @param digits Significant digits.
#' @return A data frame with columns `quantity`, `value`, `rule` and `status`.
#' @export
bridge_report <- function(x, digits = 3) {
  stopifnot(inherits(x, "bridge_diagnostics"))
  st <- x$settings
  f <- function(v) if (is.null(v) || !is.finite(v)) "--" else format(signif(v, digits), scientific = FALSE)
  ok <- function(flag) if (is.na(flag)) "--" else if (flag) "flag" else "pass"
  rows <- list(
    c("log Z-hat", format(round(x$logml, 4), nsmall = 4), "", ""),
    c("MCSE", f(x$mcse), if (st$iid) "independent draws" else sprintf("long-run variance (%s)", st$tau_method), ""),
    c(if (st$iid) "overlap check: n_ov-hat" else "overlap check: n_ovP-hat", f(if (st$iid) x$n_ov else x$n_ovP),
      sprintf(">= %s", f(st$n0)), ok(x$flags[["overlap"]])),
    c("variance check: V_jack", f(x$V_jack), sprintf("<= %s", f(x$critical[["V_jack"]])), ok(x$flags[["variance"]])),
    c("finite-sample interval", sprintf("[%s, %s]", f(x$interval[["lower"]]), f(x$interval[["upper"]])),
      sprintf("level %g, any sample size", 1 - st$delta), ""),
    c("decision", x$decision, "", "")
  )
  out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(out) <- c("quantity", "value", "rule", "status")
  out
}

#' Log Bayes factor from two independent bridge runs
#'
#' Combines the diagnostics of two independent runs, one per model, into an
#' estimate of \eqn{\log B_{12} = \log Z_1 - \log Z_2}. The MCSE is
#' \eqn{(\hat\sigma_1^2 + \hat\sigma_2^2)^{1/2}}. Because the Kolmogorov
#' distance to normality is subadditive under independent convolution, the
#' standardized error of the difference is within \eqn{D_1 + D_2} of normal,
#' so the difference passes when both runs pass; its
#' finite-sample interval, for independent draws, is
#' \eqn{[\ell_1 - u_2, u_1 - \ell_2]} with level \eqn{1 - \delta_1 - \delta_2}.
#'
#' @param d1,d2 `bridge_diagnostics` objects for models 1 and 2, from
#'   independent runs.
#' @param level Level of the normal interval.
#' @return An object of class `bridge_bayes_factor`.
#' @export
bridge_bayes_factor <- function(d1, d2, level = 0.95) {
  stopifnot(inherits(d1, "bridge_diagnostics"), inherits(d2, "bridge_diagnostics"))
  z <- stats::qnorm(1 - (1 - level) / 2)
  est <- d1$logml - d2$logml
  se <- sqrt(d1$mcse^2 + d2$mcse^2)
  pass1 <- !any(d1$flags %in% TRUE)
  pass2 <- !any(d2$flags %in% TRUE)
  fsi <- c(lower = d1$interval[["lower"]] - d2$interval[["upper"]], upper = d1$interval[["upper"]] - d2$interval[["lower"]])
  out <- list(log_bf = est, mcse = se, interval = c(lower = est - z * se, upper = est + z * se),
              fs_interval = fsi, fs_level = 1 - d1$settings$delta - d2$settings$delta,
              passes = pass1 && pass2, pass = c(model1 = pass1, model2 = pass2),
              share_model1 = d1$mcse^2 / se^2, level = level,
              actions = list(model1 = d1$actions, model2 = d2$actions))
  class(out) <- "bridge_bayes_factor"
  out
}

#' @export
print.bridge_bayes_factor <- function(x, digits = 4, ...) {
  f <- function(v) format(signif(v, digits))
  cat(sprintf("log Bayes factor (model 1 vs model 2): %s   MCSE %s\n", f(x$log_bf), f(x$mcse)))
  cat(sprintf("  %g%% interval [%s, %s]; model 1 contributes %s%% of the variance\n", 100 * x$level,
              f(x$interval[["lower"]]), f(x$interval[["upper"]]), format(round(100 * x$share_model1))))
  if (all(is.finite(x$fs_interval))) cat(sprintf("  finite-sample interval [%s, %s] (level %g)\n", f(x$fs_interval[["lower"]]), f(x$fs_interval[["upper"]]), x$fs_level))
  if (x$passes) {
    cat("  Both runs pass: report the interval.\n")
  } else {
    cat("  Flagged: run(s)", paste(names(x$pass)[!x$pass], collapse = " and "), "raised a flag.\n")
    for (m in names(x$actions)) for (a in x$actions[[m]]) cat("    - ", m, ": ", a, "\n", sep = "")
  }
  invisible(x)
}
