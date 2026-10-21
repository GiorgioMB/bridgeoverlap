#!/usr/bin/env Rscript
# Supplement Table NP: the near-perfect regime (Remark 1). p = N(0, I_d) and a
# full-covariance Gaussian proposal fitted by moments on an independent sample of size S
# (S1 = S2 = S, 500 replications): predicted beta = (d^2 + 3d)/(8S) and n_ov = (d^2 + 3d)/16
# against medians of 1 - I-hat and n_ov-hat, RMSE, median MCSE, coverage, the overlap flag
# n_ov-hat < 10 without and with the gate I-hat <= 1/2, and the plug-in Lyapunov ratio.
# Bottom: F1 (d = 10, sigma = 1.5) with S1 = 200 and S2 in {200, 4000, 40000}.
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(dirname(normalizePath(f))) else if (dir.exists("analysis")) "analysis" else ".."})
for (f in c("common.R", "families.R", "exact.R", "ar1.R", "supplement.R", "paper.R")) source(file.path(BO_ROOT, "R", f))

logZ <- 0.9
R <- n_reps(REPS$np)
summ <- function(res, extra) {
  data.frame(extra, R = nrow(res), med_one_minus_I = stats::median(1 - res$I),
             med_I = stats::median(res$I), med_n_ov = stats::median(res$n_ov),
             rmse = sqrt(mean(res$err^2)), med_mcse = stats::median(res$mcse),
             coverage = mean(abs(res$err) <= z975 * res$mcse), flag_ungated = mean(res$ungated),
             fs_half = stats::median(res$fs_half), fs_cover = mean(res$fs_cover), med_lyapunov = stats::median(res$lyap), row.names = NULL)
}
one <- function(ap, aq) {
  d <- bridge_diagnostics(ap, aq, iid = TRUE, eb = TRUE, jackknife = FALSE, rule = "plain", eb_method = "bennett")   # eb: Proposition 2
  data.frame(I = d$I, n_ov = d$n_ov, err = d$logml - logZ, mcse = d$mcse,
             ungated = d$flags[["overlap"]], fs_half = (d$eb[["upper"]] - d$eb[["lower"]]) / 2,
             fs_cover = d$eb[["lower"]] <= logZ && logZ <= d$eb[["upper"]], lyap = lyapunov_plugin(ap, aq, d$logml))
}
cfg <- expand.grid(S = c(1000L, 4000L), d = c(2L, 4L, 8L, 11L, 16L))
if (BO_QUICK) cfg <- cfg[c(1, 6), ]
top <- do.call(rbind, lapply(seq_len(nrow(cfg)), function(i) {
  d <- cfg$d[i]; S <- cfg$S[i]
  msg("NP: d =", d, "S =", S)
  res <- run_reps(R, function(r) {
    fit <- fit_gauss(matrix(stats::rnorm(S * d), S, d))   # independent fitting sample
    X <- matrix(stats::rnorm(S * d), S, d)
    Y <- rgauss(S, fit)
    one(logZ + lmvn_std(X) - lgauss(X, fit), logZ + lmvn_std(Y) - lgauss(Y, fit))
  })
  summ(res, data.frame(d = d, S = S, beta_pred = (d^2 + 3 * d) / (8 * S), n_ov_pred = (d^2 + 3 * d) / 16))
}))
save_table(top, "tab_NP_top")
fam <- family_F1(10, 1.5)
bottom <- do.call(rbind, lapply(quick_pick(c(200L, 4000L, 40000L)), function(S2) {
  msg("NP bottom: S1 = 200, S2 =", S2)
  res <- run_reps(R, function(r) one(logZ + fam$lr(fam$rpost(200L)), logZ + fam$lr(fam$rprop(S2))))
  summ(res, data.frame(S1 = 200L, S2 = S2))
}))
save_table(bottom, "tab_NP_bottom")
print(top, digits = 3); print(bottom, digits = 3)
sub_hdr <- "\\multicolumn{2}{l}{$\\Sone$, $\\Stwo$} & & med $\\Ihat$ & & med $\\neffhat$ & & & & & & & \\\\"
bot_rows <- sprintf("\\multicolumn{2}{l}{$%d$, $%d$} & & %s & & %s & %s & %s & %s & %s & %s & %s & %s\\\\", bottom$S1, bottom$S2,
                    fx(bottom$med_I, 3), fx(bottom$med_n_ov, 1), fx(bottom$rmse, 4), fx(bottom$med_mcse, 4),
                    fx(bottom$coverage, 3), fx(bottom$flag_ungated, 2), fx(bottom$fs_half, 4), fx(bottom$fs_cover, 3), fx(bottom$med_lyapunov, 3))
emit_table("tab_NP", cbind(top$d, top$S, fx(top$beta_pred, 4), fx(top$med_one_minus_I, 4), fx(top$n_ov_pred, 1),
                           fx(top$med_n_ov, 1), fx(top$rmse, 4), fx(top$med_mcse, 4), fx(top$coverage, 3),
                           fx(top$flag_ungated, 2), fx(top$fs_half, 4), fx(top$fs_cover, 3), fx(top$med_lyapunov, 3)),
           extra = c("\\midrule", sub_hdr, bot_rows))
dev <- 100 * abs(c(top$med_one_minus_I / top$beta_pred, top$med_n_ov / top$n_ov_pred) - 1)
res_put("np.R", fint(R)); res_put("np.dev", frange(dev, 0)); res_put("np.rmse.max", fx(max(top$rmse), 3))
res_put("np.rmse.d11", fx(max(top$rmse[top$d <= 11]), 3))
res_put("np.cov", frange(top$coverage, 2))
u8 <- top$flag_ungated[top$d <= 8]
res_put("np.ungated.d8", if (length(u8) && min(u8) == 1) "every run" else sprintf("%s\\%% of the runs", fx(100 * min(u8), 0)))
res_put("np.ungated.d11", fx(100 * mean(top$flag_ungated[top$d == 11]), 0))
res_put("np.fs.half", frange(top$fs_half, 3)); res_put("np.fs.cov", frange(top$fs_cover, 3)); res_put("np.fs.ratio", frange(top$fs_half / (1.96 * top$med_mcse), 0))
res_put("np.I.min", fx(min(top$med_I), 2)); res_put("np.lyap", frange(top$med_lyapunov, 2))
ub <- bottom[bottom$S2 == max(bottom$S2), ]
res_put("np.unbal.S2", fint(ub$S2)); res_put("np.unbal.I", fx(ub$med_I, 3)); res_put("np.unbal.nov", fx(ub$med_n_ov, 1))
res_write("supplement/S04_NP_near_perfect.R")
