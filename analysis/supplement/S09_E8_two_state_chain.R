#!/usr/bin/env Rscript
# Supplement Table E8: the two-point family with eps = 0.05, S1 = S2 = S and posterior
# draws from the stationary two-state chain of Theorem 7 with spectral gap gamma
# (gamma = 1: independent draws). Exact computations: log10 P(Z-hat in {0, Inf}) and its
# lower bound log10{(1 - eps)(1 - eps gamma)^(S-1)}, Kolmogorov distance of Delta/sigma_NP
# to N(0, 1) and the same times gamma sqrt(v_P) (= sqrt(n_ovP)), and exact coverage of
# log Z-hat +- 1.96 sigma_NP and +- 1.96 sigma_N.
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(dirname(normalizePath(f))) else if (dir.exists("analysis")) "analysis" else ".."})
for (f in c("common.R", "families.R", "exact.R", "ar1.R", "supplement.R", "paper.R")) source(file.path(BO_ROOT, "R", f))

eps <- 0.05
grid <- expand.grid(n_ovP = c(1, 2, 4, 8, 16, 32), gamma = c(1, 0.5, 0.2, 0.05))
if (BO_QUICK) grid <- grid[c(1, 6, 19, 24), ]
tab <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  gam <- grid$gamma[i]
  S <- round(2 * grid$n_ovP[i] / (gam * eps * (1 - eps)))
  msg("E8: gamma =", gam, "S =", S)
  res <- two_state_chain(S, eps, gam)
  lqk <- S * log1p(-eps)                                   # log P(n_prop = 0)
  lqm <- log1p(-eps) + (S - 1) * log1p(-eps * gam)         # log P(n_post = 0)
  data.frame(gamma = gam, S = S, n_ov = res$n_ov, n_ovP = res$n_ovP,
             log10_p_nonfinite = lse2(lqk, log1p(-exp(lqk)) + lqm) / log(10),
             log10_lower_bound = lqm / log(10), kolmogorov = res$ks,
             scaled_kolmogorov = sqrt(gam * res$n_ov) * res$ks,
             coverage_sigmaNP = res$coverage_sigmaNP, coverage_sigmaN = res$coverage_sigmaN)
}))
save_table(tab, "tab_E8")
print(tab, digits = 4)
emit_table("tab_E8", cbind(format(tab$gamma), fint(tab$S), fx(tab$n_ov, 1), fx(tab$n_ovP, 0), fx(tab$log10_p_nonfinite, 1),
                           fx(tab$log10_lower_bound, 1), fx(tab$kolmogorov, 4), fx(tab$scaled_kolmogorov, 3),
                           fx(tab$coverage_sigmaNP, 4), fx(tab$coverage_sigmaN, 4)), groups = group_ends(tab$gamma))
spread <- function(m) max(vapply(split(tab$coverage_sigmaNP[tab$n_ovP >= m], round(tab$n_ovP[tab$n_ovP >= m])), function(v) diff(range(v)), 0))
res_put("e8.spread2", fx(spread(1.5), 3)); res_put("e8.spread16", fx(spread(15.5), 3))
res_put("e8.c", frange(tab$scaled_kolmogorov[tab$n_ovP >= 3.5], 2))
res_put("e8.ign.005", fx(stats::median(tab$coverage_sigmaN[tab$gamma == 0.05]), 2))
res_write("supplement/S09_E8_two_state_chain.R")
