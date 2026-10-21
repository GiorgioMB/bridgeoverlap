#!/usr/bin/env Rscript
# Supplement Table E7: exact computations in the two-point family of Theorem 2
# (S1 = S2 = S): Kolmogorov distance of Delta/sigma_N to N(0,1) with the bounds
# 0.04/sqrt(n_ov) (Theorem 2(b), n_ov >= 2) and 3.1/sqrt(n_ov) (Theorem 3), log10 P(Z-hat in
# {0, Inf}) with the bound e^{-4 N H}, delta-method coverage, and P(|Delta| >= 0.25)
# with the bound 2 exp{-N H kappa(0.25)} of Theorem 1.
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(dirname(normalizePath(f))) else if (dir.exists("analysis")) "analysis" else ".."})
for (f in c("common.R", "families.R", "exact.R", "ar1.R", "supplement.R", "paper.R")) source(file.path(BO_ROOT, "R", f))

rows <- data.frame(S = c(50, 100, 200, 400, 1000, 2000, 4000, 10000, 20000, 50000, 1e5, 1e6, 1e6),
                   eps = c(0.02, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01, 0.005, 0.005, 0.004, 0.004, 0.001, 0.004))
if (BO_QUICK) rows <- rows[c(1, 5, 9), ]
k25 <- kappa_fun(0.25)
tab <- do.call(rbind, lapply(seq_len(nrow(rows)), function(i) {
  S <- rows$S[i]; eps <- rows$eps[i]
  msg("E7: S =", S, "eps =", eps)
  cst <- two_point_constants(S, eps)
  NH <- 2 * S * cst$H
  pk <- stats::dbinom(0:S, S, eps)
  law <- two_point_law(S, pk, pk)
  k <- law$k; m <- law$m; w <- law$w
  q <- k / (k + m)
  n_ov_hat <- ((1 - q)^2 * k * (S - k) + q^2 * m * (S - m)) / (S - 1)
  sig_hat <- sqrt(n_ov_hat) / (m * q)
  q0 <- (1 - eps)^S
  data.frame(S = S, eps = eps, n_ov = cst$n_ov,
             kolmogorov = ks_discrete(law$Delta / cst$sigma_N, w, law$p_minus, law$p_plus),
             lower_bound = if (cst$n_ov >= 2) 0.04 / sqrt(cst$n_ov) else NA_real_,
             upper_bound = min(1, 3.1 / sqrt(cst$n_ov)),
             log10_p_nonfinite = (S * log1p(-eps) + log(2 - q0)) / log(10),
             log10_bound = -4 * NH / log(10),
             coverage = sum(w[abs(law$Delta) <= z975 * sig_hat]),
             p_abs_delta_ge_025 = law$p_minus + law$p_plus + sum(w[abs(law$Delta) >= 0.25]),
             bound_025 = min(1, 2 * exp(-NH * k25)))
}))
save_table(tab, "tab_E7")
print(tab, digits = 4)
emit_table("tab_E7", cbind(fint(tab$S), format(tab$eps), fx(tab$n_ov, 2), fx(tab$kolmogorov, 4), fx(tab$lower_bound, 4),
                           fx(tab$upper_bound, 3), fx(tab$log10_p_nonfinite, 1), fx(tab$log10_bound, 1), fx(tab$coverage, 4),
                           fx(tab$p_abs_delta_ge_025, 4), fx(tab$bound_025, 4)))
big <- tab$n_ov >= 2
res_put("e7.ksconst", frange(sqrt(tab$n_ov[big]) * tab$kolmogorov[big], 2))
res_write("supplement/S02_E7_two_point.R")
