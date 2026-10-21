#!/usr/bin/env Rscript
# Supplement Table FA: joint false-alarm rate (%) of the failure diagnostics under
# correct sampling: p = N(0, I_5), g = N(0.2 * 1, 1.3^2 I_5), S1 = S2 = S, G = 20, 2000
# replications. Marginal rejection rates at nominal 5%, union of the three jackknife
# diagnostics, union of all four, and union of all four at alpha/4 each (Algorithm 1).
# The AR(1) rows use Geyer's estimator for the MCSE.
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(dirname(normalizePath(f))) else if (dir.exists("analysis")) "analysis" else ".."})
for (f in c("common.R", "families.R", "exact.R", "ar1.R", "supplement.R", "paper.R")) source(file.path(BO_ROOT, "R", f))

logZ <- 0.15
R <- n_reps(REPS$fa)
G <- 20L
invisible(jackknife_null_quantiles(G, probs = c(0.95, 1 - 0.05 / 4)))
gm <- rep(0.2, 5); gv <- rep(1.3^2, 5)
rows <- data.frame(setting = c("independent, S = 1000", "independent, S = 4000",
                               "AR(1), rho = 0.5, S = 8000", "AR(1), rho = 0.9, S = 8000"),
                   S = c(1000L, 4000L, 8000L, 8000L), rho = c(0, 0, 0.5, 0.9), iid = c(TRUE, TRUE, FALSE, FALSE))
if (BO_QUICK) rows <- rows[c(1, 3), ]
tab <- do.call(rbind, lapply(seq_len(nrow(rows)), function(i) {
  rw <- rows[i, ]
  msg("FA:", rw$setting)
  res <- run_reps(R, function(r) {
    X <- if (rw$iid) matrix(stats::rnorm(rw$S * 5), rw$S, 5) else ar1_matrix(rw$S, 5L, rw$rho)
    Y <- sweep(matrix(stats::rnorm(rw$S * 5), rw$S, 5) * 1.3, 2, gm, "+")
    d <- bridge_diagnostics(logZ + lmvn_std(X) - ldiag(X, gm, gv), logZ + lmvn_std(Y) - ldiag(Y, gm, gv),
                            iid = rw$iid, G = G, eb = FALSE, rule = "plain")
    rj <- rejections(d)
    data.frame(t_T = rj[["t_T"]], V_jack = rj[["V_jack"]], I_max = rj[["I_max_exact"]], G_max = rj[["G_max"]],
               union_a4 = any(d$flags[-1L] %in% TRUE), n_ov = d$n_ov, I = d$I, f_ov = d$flags[["overlap"]])
  })
  jk <- res$V_jack | res$I_max | res$G_max
  data.frame(setting = rw$setting, R = nrow(res), t_T = 100 * mean(res$t_T), V_jack = 100 * mean(res$V_jack),
             I_max = 100 * mean(res$I_max), G_max = 100 * mean(res$G_max),
             jackknife_union = 100 * mean(jk), union_of_four = 100 * mean(jk | res$t_T),
             union_alpha4_each = 100 * mean(res$union_a4), med_n_ov = stats::median(res$n_ov),
             med_I = stats::median(res$I), overlap = mean(res$f_ov))
}))
save_table(tab, "tab_FA")
print(tab, digits = 3)
lab <- c("independent, S = 1000" = "independent, $S=1000$", "independent, S = 4000" = "independent, $S=4000$",
         "AR(1), rho = 0.5, S = 8000" = "AR(1), $\\rho=0.5$, $S=8000$", "AR(1), rho = 0.9, S = 8000" = "AR(1), $\\rho=0.9$, $S=8000$")
emit_table("tab_FA", cbind(lab[tab$setting], fx(tab$t_T, 1), fx(tab$V_jack, 1), fx(tab$I_max, 1), fx(tab$G_max, 1),
                           fx(tab$jackknife_union, 1), fx(tab$union_of_four, 1), fx(tab$union_alpha4_each, 1)))
# union probabilities of the jackknife diagnostics under their joint fixed-G limit
npiv <- if (BO_QUICK) 4e4 else 4e5
q5 <- jackknife_null_quantiles(G, probs = 0.95)
set.seed(BO_SEED + 777)
Z <- matrix(stats::rnorm(npiv * G), npiv, G); Zc <- Z - rowMeans(Z); A <- abs(Zc)
D <- A[cbind(seq_len(npiv), max.col(A, ties.method = "first"))]; ss <- rowSums(Zc^2)
V <- ss / (G - 1); Im <- D * sqrt(G) / (G - 1); Gm <- D / sqrt(ss / (G - 1))
rI <- Im > q5["I_max", 1]; rG <- Gm > q5["G_max", 1]; rV <- V > stats::qchisq(0.95, G - 1) / (G - 1)
res_put("fa.R", fint(R)); res_put("fa.u4", frange(tab$union_of_four, 1)); res_put("fa.ua4", frange(tab$union_alpha4_each, 1))
res_put("fa.piv.jk", fx(mean(rV | rI | rG), 3)); res_put("fa.piv.IG", fx(mean(rI | rG), 3)); res_put("fa.piv.n", fsci(npiv, 0))
ind <- grepl("^independent", tab$setting)
res_put("fa.nov.1000", fint(tab$med_n_ov[tab$setting == "independent, S = 1000"])); res_put("fa.I", fx(stats::median(tab$med_I[ind]), 2))
res_put("fa.overlap.phrase", if (all(tab$overlap == 0)) "The overlap flag never fired in these runs" else
  sprintf("The overlap flag fired in at most %s\\%% of these runs", fx(100 * max(tab$overlap), 1)))
res_write("supplement/S07_FA_false_alarm.R")
