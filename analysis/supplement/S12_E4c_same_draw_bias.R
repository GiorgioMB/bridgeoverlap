#!/usr/bin/env Rscript
# Supplement Table E4c: scaling of the same-draw bias (t_7 target, diagonal Gaussian
# proposal fitted on the same draws, S1 = S2 = S). Bias in MCSE units, medians of n_ov-hat
# and I-hat, and C-hat = (bias/MCSE) sqrt(n_ov-hat) / (k s2 I-hat), which is constant if
# the heuristic of Section S2.5(i) holds.
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(dirname(normalizePath(f))) else if (dir.exists("analysis")) "analysis" else ".."})
for (f in c("common.R", "families.R", "exact.R", "ar1.R", "supplement.R", "paper.R")) source(file.path(BO_ROOT, "R", f))
source(file.path(BO_ROOT, 'R', 'proposal_fitting.R'))

cfg <- data.frame(d = c(4L, 16L, 8L, 8L, 8L, 8L, 2L, 8L, 32L),
                  S = c(2000L, 2000L, 500L, 1000L, 4000L, 8000L, 2000L, 2000L, 2000L),
                  R = REPS$e4c)
if (BO_QUICK) cfg <- cfg[c(3, 9), ]
tab <- do.call(rbind, lapply(seq_len(nrow(cfg)), function(i) {
  dd <- cfg$d[i]; S <- cfg$S[i]; tg <- t7_target(dd)
  msg("E4c: d =", dd, "S =", S)
  res <- run_reps(n_reps(cfg$R[i], key = "REPS$e4c"), function(r) {
    X <- tg$r(S); ms <- moments(X)
    a <- ratios(tg, X, prop_from(matrix(stats::rnorm(S * dd), S, dd), ms$m, ms$v), ms$m, ms$v)
    d <- bridge_diagnostics(a$ap, a$aq, iid = TRUE, eb = FALSE, jackknife = FALSE)
    data.frame(D = d$logml - tg$logZ, mcse = d$mcse, n_ov = d$n_ov, I = d$I)
  })
  b <- mean(res$D) / stats::median(res$mcse)
  k <- 2L * dd
  data.frame(d = dd, k = k, S = S, R = nrow(res), bias_in_mcse = b,
             bias_se = stats::sd(res$D) / sqrt(nrow(res)) / stats::median(res$mcse),
             med_n_ov = stats::median(res$n_ov), med_I = stats::median(res$I),
             C_hat = b * sqrt(stats::median(res$n_ov)) / (k * 0.5 * stats::median(res$I)))
}))
save_table(tab, "tab_E4c")
print(tab, digits = 3)
emit_table("tab_E4c", cbind(tab$d, tab$k, tab$S, sprintf("$%s$ (%s)", fsign(tab$bias_in_mcse, 2), fx(tab$bias_se, 2)),
                            fint(tab$med_n_ov), fx(tab$med_I, 3), fx(tab$C_hat, 2)))
res_put("e4c.C", if (diff(range(tab$C_hat)) == 0) sprintf("at $%s$", fx(tab$C_hat[1], 2)) else sprintf("between $%s$ and $%s$", fx(min(tab$C_hat), 2), fx(max(tab$C_hat), 2))); res_put("e4c.absbias", frange(abs(tab$bias_in_mcse), 2))
res_write("supplement/S12_E4c_same_draw_bias.R")
