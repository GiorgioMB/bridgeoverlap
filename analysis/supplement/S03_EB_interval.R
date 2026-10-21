#!/usr/bin/env Rscript
# Supplement Table EB: the finite-sample interval of Proposition 2 (delta = 0.05,
# S1 = S2 = S, 1000 replications per row): medians of n_ov-hat and I-hat, coverage of the
# delta-method and empirical Bernstein intervals, proportion of bounded intervals and
# median width ratio (EB / delta method).
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(dirname(normalizePath(f))) else if (dir.exists("analysis")) "analysis" else ".."})
for (f in c("common.R", "families.R", "exact.R", "ar1.R", "supplement.R", "paper.R")) source(file.path(BO_ROOT, "R", f))

logZ <- -0.6
R <- n_reps(REPS$eb)
cfgs <- list(list(fam = family_F2(2), S = c(30, 100, 300, 1000, 3000)),
             list(fam = family_F2(4), S = c(100, 300, 1000, 3000)),
             list(fam = family_F1(10, 1.5), S = c(30, 100, 300, 1000)),
             list(fam = family_F3(20), S = c(30, 100, 300, 1000)))
tab <- do.call(rbind, lapply(quick_pick(cfgs, c(1L, 4L)), function(cc) {
  do.call(rbind, lapply(quick_pick(cc$S), function(S) {
    msg("EB:", cc$fam$label, "S =", S)
    res <- run_reps(R, function(r) {
      ap <- logZ + cc$fam$lr(cc$fam$rpost(S)); aq <- logZ + cc$fam$lr(cc$fam$rprop(S))
      d <- bridge_diagnostics(ap, aq, iid = TRUE, eb = TRUE, jackknife = FALSE, eb_method = "bennett")
      lo <- d$eb[["lower"]]; hi <- d$eb[["upper"]]                         # the interval of Proposition 2
      data.frame(n_ov = d$n_ov, I = d$I,
                 cov_delta = is.finite(d$logml) && isTRUE(abs(d$logml - logZ) <= z975 * d$mcse),
                 cov_eb = lo <= logZ && logZ <= hi, bounded = is.finite(lo) && is.finite(hi),
                 ratio = (hi - lo) / (2 * z975 * d$mcse))
    })
    data.frame(family = cc$fam$label, S = S, R = nrow(res), med_n_ov = stats::median(res$n_ov),
               med_I = stats::median(res$I), coverage_delta = mean(res$cov_delta),
               coverage_eb = mean(res$cov_eb), bounded = mean(res$bounded),
               width_ratio = stats::median(res$ratio))
  }))
}))
save_table(tab, "tab_EB")
print(tab, digits = 3)
lab <- c("F2 (mu_g=2)" = "F2, $\\mu_g=2$", "F2 (mu_g=4)" = "F2, $\\mu_g=4$",
         "F1 (d=10, sigma=1.5)" = "F1, $d=10$, $\\sigma=1.5$", "F3 (d=20)" = "F3, $d=20$")
emit_table("tab_EB", cbind(lab[tab$family], tab$S, ifelse(tab$med_n_ov < 100, fx(tab$med_n_ov, 1), fx(tab$med_n_ov, 0)),
                           fx(tab$med_I, 2), fx(tab$coverage_delta, 3), fx(tab$coverage_eb, 3), fx(tab$bounded, 2),
                           fx(tab$width_ratio, 1)))
res_put("eb.R", fint(R))
res_put("eb.cover.phrase", if (min(tab$coverage_eb) == 1) "covered in every replication of every configuration" else
  sprintf("covered in at least a proportion %s of the replications in every configuration", fx(min(tab$coverage_eb), 3)))
w <- which.min(tab$bounded)
res_put("eb.bounded.phrase", if (min(tab$bounded) == 1) "was bounded in every run" else
  sprintf("was bounded in all but %s\\%% of the runs of %s at $S=%d$", fx(100 * (1 - tab$bounded[w]), 0),
          sub(",.*", "", lab[tab$family[w]]), tab$S[w]))
hi <- tab$med_n_ov >= 90; md <- tab$med_n_ov >= 10 & tab$med_n_ov <= 40
res_put("eb.w.hi", frange(tab$width_ratio[hi], 1)); res_put("eb.nov.hi", frange(tab$med_n_ov[hi], 0))
res_put("eb.w.mid", frange(tab$width_ratio[md], 0)); res_put("eb.nov.mid", frange(tab$med_n_ov[md], 0))
res_put("eb.w.max", fx(max(tab$width_ratio), 0)); res_put("eb.w.max.nov", fx(tab$med_n_ov[which.max(tab$width_ratio)], 0))
res_write("supplement/S03_EB_interval.R")
