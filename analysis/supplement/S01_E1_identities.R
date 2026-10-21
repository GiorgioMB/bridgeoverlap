#!/usr/bin/env Rscript
# Supplement Table E1: identities of Lemma 1 and Theorem S1 (the Meng-Wong iteration).
# Top: iterations of the MW map to reach |x_{t+1} - x_t| <= 1e-13 from c-hat +- 20 and
# c-hat +- 200 (worst case over the two signs), final error, sup |phi'| over 4001 points
# in [c-hat - 60, c-hat + 60], and T'(Z-hat) by formula (a), finite differences and
# 1 - H/H0 (S1 = S2 = 2000). Bottom: T' over 1000 replications for F2 (mu_g = 2);
# sqrt(n_ov) sd(T') should be constant at sqrt(beta) sd_{m*}(pi*) = 0.202.
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(dirname(normalizePath(f))) else if (dir.exists("analysis")) "analysis" else ".."})
for (f in c("common.R", "families.R", "exact.R", "ar1.R", "supplement.R", "paper.R")) source(file.path(BO_ROOT, "R", f))

S <- 2000L
logZ <- 0.3
problems <- list("F1: N(0,I_10) vs N(0,1.5^2 I_10)" = family_F1(10, 1.5),
                 "F2: N(0,1) vs N(4,1)" = family_F2(4),
                 "F3: t_3(0,I_20) vs N(0,3 I_20)" = family_F3(20))
lmean <- function(lv) lse(lv) - log(length(lv))
# The MW map on the log scale, phi(x) = log T(e^x), computed with log-probabilities.
mw_phi <- function(x, ap, aq, lr) {
  x + lmean(stats::plogis(x - aq - lr, lower.tail = FALSE, log.p = TRUE)) -
    lmean(stats::plogis(x - ap - lr, log.p = TRUE)) - lr
}
# Formula (a) of Theorem S1: phi'(x) = A(x) - B(x).
mw_dphi <- function(x, ap, aq, lr) {
  lpp <- stats::plogis(x - ap - lr, log.p = TRUE)
  lq <- stats::plogis(x - aq - lr, log.p = TRUE)
  l1q <- stats::plogis(x - aq - lr, lower.tail = FALSE, log.p = TRUE)
  sum(exp(2 * lpp - lse(lpp))) - sum(exp(lq + l1q - lse(l1q)))
}
iterate_mw <- function(x0, ap, aq, lr, tol = 1e-13, maxit = 10000L) {
  x <- x0
  for (t in seq_len(maxit)) {
    xn <- mw_phi(x, ap, aq, lr)
    if (abs(xn - x) <= tol) return(c(iter = t, x = xn))
    x <- xn
  }
  c(iter = NA_real_, x = x)
}

set.seed(BO_SEED)
top <- do.call(rbind, lapply(names(problems), function(nm) {
  fam <- problems[[nm]]
  ap <- logZ + fam$lr(fam$rpost(S))
  aq <- logZ + fam$lr(fam$rprop(S))
  lr <- 0
  ch <- bridge_logml(ap, aq, tol = 1e-15)$logml
  runs <- lapply(c(-20, 20, -200, 200), function(o) iterate_mw(ch + o, ap, aq, lr))
  grid <- seq(ch - 60, ch + 60, length.out = 4001)
  h <- 1e-5
  d <- bridge_diagnostics(ap, aq, iid = TRUE, eb = FALSE, jackknife = FALSE)
  data.frame(problem = nm,
             iter_pm20 = max(runs[[1]][["iter"]], runs[[2]][["iter"]]),
             iter_pm200 = max(runs[[3]][["iter"]], runs[[4]][["iter"]]),
             final_error = max(abs(vapply(runs, function(r) r[["x"]], 0) - ch)),
             sup_abs_dphi = max(abs(vapply(grid, mw_dphi, 0, ap = ap, aq = aq, lr = lr))),
             T_formula = mw_dphi(ch, ap, aq, lr),
             T_finite_diff = (mw_phi(ch + h, ap, aq, lr) - mw_phi(ch - h, ap, aq, lr)) / (2 * h),
             T_one_minus_H = d$T)
}))
save_table(top, "tab_E1_top")
print(top, digits = 6)

fam <- family_F2(2)
pop <- overlap_population(fam)
ov <- overlap_1d(function(x) stats::dnorm(x, log = TRUE), function(x) stats::dnorm(x, 2, log = TRUE),
                 1000, 1000, c(-5, -2, 0, 1, 2, 4, 7))
limit <- sqrt(ov$beta) * sqrt(mstar_moments(ov)$var)
bottom <- do.call(rbind, lapply(quick_pick(c(250L, 1000L, 4000L, 16000L)), function(S) {
  msg("E1 bottom: S =", S)
  Tv <- run_reps(n_reps(REPS$e1), function(r) {
    d <- bridge_diagnostics(logZ + fam$lr(fam$rpost(S)), logZ + fam$lr(fam$rprop(S)),
                            iid = TRUE, eb = FALSE, jackknife = FALSE)
    data.frame(T = d$T)
  })$T
  n_ov <- overlap_at(pop, S)$n_ov
  data.frame(S = S, n_ov = n_ov, mean_T = mean(Tv), sd_T = stats::sd(Tv),
             sqrt_n_ov_sd_T = sqrt(n_ov) * stats::sd(Tv), limit = limit)
}))
save_table(bottom, "tab_E1_bottom")
print(bottom, digits = 4)
lab <- c("F1: $\\mathcal N(0,I_{10})$ vs $\\mathcal N(0,1.5^2I_{10})$", "F2: $\\mathcal N(0,1)$ vs $\\mathcal N(4,1)$",
         "F3: $t_3(0,I_{20})$ vs $\\mathcal N(0,3I_{20})$")
emit_table("tab_E1_top", cbind(lab[seq_len(nrow(top))], top$iter_pm20, top$iter_pm200, sprintf("%.0e", top$final_error),
                               fx(top$sup_abs_dphi, 3), fx(top$T_formula, 6), fx(top$T_finite_diff, 6), fx(top$T_one_minus_H, 6)))
emit_table("tab_E1_bottom", cbind(bottom$S, fx(bottom$n_ov, 1), fsign(bottom$mean_T, 5), fx(bottom$sd_T, 5),
                                  fx(bottom$sqrt_n_ov_sd_T, 3)))
res_put("e1.iter", frange(c(top$iter_pm20, top$iter_pm200), 0))
res_put("e1.limit", fx(limit, 3)); res_put("e1.obs", frange(bottom$sqrt_n_ov_sd_T, 3)); res_put("e1.R", fint(n_reps(REPS$e1)))
res_write("supplement/S01_E1_identities.R")
