# Worked printout and Bayes factors ------------------------------------------------

#' Tabular report of the diagnostics for one run
#'
#' Arranges the output of [bridge_diagnostics()] as a table: each reported
#' quantity, its value, the rule applied to it and whether the rule is met.
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
  th <- x$thresholds
  rule <- st$rule %||% "gate"
  rows <- list(
    c("log Z-hat", format(round(x$logml, 4), nsmall = 4), "", ""),
    c("MCSE", f(x$mcse), if (st$iid) "independent draws" else sprintf("long-run variance (%s)", st$tau_method), ""),
    c("n_ov-hat", f(x$n_ov), if (st$iid) sprintf(">= %s (%s rule)", f(x$n_ov_threshold), rule) else "", if (st$iid) ok(x$flags[["overlap"]]) else ""),
    c("I-hat", f(x$I), if (rule == "plain") "near 1: near-perfect, report the finite-sample interval" else sprintf("threshold relaxed above %g", th$I_gate), ""),
    c("N H0-hat", f(x$NH0), "", ""),
    c("lower bound for N H", f(x$NH_lower), sprintf(">= %s for the accuracy bound", f(x$NH_required)), "")
  )
  if (!st$iid) {
    rows <- c(rows, list(c("tau-hat", f(x$tau), "within chains", ""),
                         c("n_ovP-hat", f(x$n_ovP), sprintf(">= %s (%s rule)", f(x$n_ov_threshold), rule), ok(x$flags[["overlap"]]))))
  }
  if (isTRUE(st$eb)) {
    rows <- c(rows, list(c("finite-sample interval", sprintf("[%s, %s]", f(x$eb[["lower"]]), f(x$eb[["upper"]])),
                           sprintf("level %g, any sample size", 1 - st$delta), "")))
  }
  if (st$iid) rows <- c(rows, list(c("accuracy bound", f(x$accuracy_bound), sprintf("level %g", 1 - 2 * st$delta), "")))
  if (isTRUE(st$jackknife)) {
    rows <- c(rows, list(
      c("T (t_T)", sprintf("%s (%s)", f(x$T), f(x$t_T)), sprintf("|t_T| <= %s", f(th$t_T)), ok(x$flags[["calibration"]])),
      c("V_jack", f(x$V_jack), sprintf("<= %s", f(th$V_jack)), ok(x$flags[["jackknife_variance"]])),
      c("I_max", f(x$I_max), sprintf("<= %s", f(th$I_max)), ok(x$flags[["influential_block_imax"]])),
      c("G_max", f(x$G_max), sprintf("<= %s", f(th$G_max)), ok(x$flags[["influential_block_grubbs"]]))))
  }
  rows <- c(rows, list(c("decision", x$decision %||% "", "", "")))
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
#' so the difference is certified when both runs are; its
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
  pass1 <- !any(d1$flags %in% TRUE) && !isTRUE(d1$stuck_chains > 0)
  pass2 <- !any(d2$flags %in% TRUE) && !isTRUE(d2$stuck_chains > 0)
  eb <- c(lower = d1$eb[["lower"]] - d2$eb[["upper"]], upper = d1$eb[["upper"]] - d2$eb[["lower"]])
  out <- list(log_bf = est, mcse = se, interval = c(lower = est - z * se, upper = est + z * se),
              eb = eb, eb_level = 1 - d1$settings$delta - d2$settings$delta,
              certified = pass1 && pass2, pass = c(model1 = pass1, model2 = pass2),
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
  if (all(is.finite(x$eb))) cat(sprintf("  finite-sample interval [%s, %s] (level %g)\n", f(x$eb[["lower"]]), f(x$eb[["upper"]]), x$eb_level))
  if (x$certified) {
    cat("  Both runs pass: report the interval.\n")
  } else {
    cat("  Not certified: run(s)", paste(names(x$pass)[!x$pass], collapse = " and "), "raised a flag.\n")
    for (m in names(x$actions)) for (a in x$actions[[m]]) cat("    - ", m, ": ", a, "\n", sep = "")
  }
  invisible(x)
}
