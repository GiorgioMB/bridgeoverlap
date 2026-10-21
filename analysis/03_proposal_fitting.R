#!/usr/bin/env Rscript
# Section 7.3 and Supplement Table prop: proposal fitting at a fixed posterior budget.
# t_7 target with independent coordinates 0.5 + s_k t_7 (s_k equally spaced in [1, 2]),
# diagonal Gaussian proposal fitted by moments (k = 2d parameters); S = 2000 posterior and
# 2000 proposal draws; 1000 replications for d = 2, 8 and 500 for d = 32. Arms: fixed at
# theta_0; same draws; split with an extra independent sample; split within the same
# budget (fit on the first half, bridge the second half against 2000 proposal draws);
# two-fold cross-fitting. Expected: same-draw bias -1.17 MCSEs at d = 32 (coverage 0.766,
# calibration-test rejections 10.4%); cross-fitting costs at most about 5% in sd.
# Tables E4 and E4c are in supplement/S11_E4_proposal_fitting.R and S12_E4c_same_draw_bias.R.
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else if (dir.exists("analysis")) "analysis" else "."})
for (f in c("common.R", "supplement.R", "proposal_fitting.R", "paper.R")) source(file.path(BO_ROOT, "R", f))

S <- 2000L
G <- 20L
invisible(jackknife_null_quantiles(G, probs = c(0.95, 1 - 0.05 / 4)))
arms <- c("fixed at theta_0", "same draws", "split, extra sample", "split, same budget", "cross-fitted")
rows <- list()
targets <- data.frame(target = c("t7", "t7", "t7", "skew", "skew"), d = c(2L, 8L, 32L, 8L, 32L),
                      rkey = c("prop", "prop", "prop_d32", "prop_skew", "prop_skew"))
if (BO_QUICK) targets <- targets[c(1, 4), ]
for (ti in seq_len(nrow(targets))) {
  dd <- targets$d[ti]
  tg <- if (targets$target[ti] == "t7") t7_target(dd) else gamma_target(dd)
  msg("prop:", targets$target[ti], "d =", dd)
  res <- run_reps(n_reps(REPS[[targets$rkey[ti]]], key = paste0("REPS$", targets$rkey[ti])),
                  cell = paste("prop", targets$target[ti], dd), fun = function(r) {
    X <- tg$r(S); Xextra <- tg$r(S); E <- matrix(stats::rnorm(S * dd), S, dd)   # common base variates
    h <- S / 2; A <- X[seq_len(h), , drop = FALSE]; B <- X[h + seq_len(h), , drop = FALSE]
    run <- function(arm, ap, aq, jk = FALSE, folds = NULL) {
      d <- bridge_diagnostics(ap, aq, iid = TRUE, eb = FALSE, jackknife = jk, G = G, folds = folds)
      data.frame(arm = arm, D = d$logml - tg$logZ, mcse = d$mcse, n_ov = d$n_ov,
                 t_rej = if (jk) isTRUE(abs(d$t_T) > stats::qt(0.975, G - 1)) else NA)
    }
    a0 <- ratios(tg, X, prop_from(E, tg$m0, tg$v0), tg$m0, tg$v0)
    ms <- moments(X); a1 <- ratios(tg, X, prop_from(E, ms$m, ms$v), ms$m, ms$v)
    me <- moments(Xextra); a2 <- ratios(tg, X, prop_from(E, me$m, me$v), me$m, me$v)
    mA <- moments(A); a3 <- ratios(tg, B, prop_from(E, mA$m, mA$v), mA$m, mA$v)
    mB <- moments(B)
    aA <- ratios(tg, B, prop_from(E[seq_len(h), , drop = FALSE], mA$m, mA$v), mA$m, mA$v)
    aB <- ratios(tg, A, prop_from(E[h + seq_len(h), , drop = FALSE], mB$m, mB$v), mB$m, mB$v)
    rbind(run(arms[1], a0$ap, a0$aq), run(arms[2], a1$ap, a1$aq, TRUE), run(arms[3], a2$ap, a2$aq, TRUE),
          run(arms[4], a3$ap, a3$aq),
          run(arms[5], c(aA$ap, aB$ap), c(aA$aq, aB$aq), TRUE, list(post = rep(1:2, each = h), prop = rep(1:2, each = h))))
  })
  for (a in arms) {
    x <- res[res$arm == a, ]
    mm <- stats::median(x$mcse)
    rows[[length(rows) + 1L]] <- data.frame(
      target = targets$target[ti], d = dd, proposal = a, R = nrow(x), sd = stats::sd(x$D), sd_se = se_sd(x$D),
      bias_in_mcse = mean(x$D) / mm, bias_se = stats::sd(x$D) / sqrt(nrow(x)) / mm,
      rmse = sqrt(mean(x$D^2)), coverage = mean(abs(x$D) <= z975 * x$mcse),
      med_n_ov = stats::median(x$n_ov), t_T_rejections_pct = 100 * mean(x$t_rej))
  }
}
tab <- do.call(rbind, rows)
save_table(tab, "tab_prop")
lab <- c("fixed at theta_0" = "fixed at $\\vartheta_0$", "same draws" = "same draws",
         "split, extra sample" = "split, extra sample", "split, same budget" = "split, same budget",
         "cross-fitted" = "cross-fitted")
tlab <- c(t7 = "$t_7$", skew = "Gamma")
tab$sd_vs_fixed <- tab$sd / ave(ifelse(tab$proposal == "fixed at theta_0", tab$sd, NA), tab$target, tab$d, FUN = function(v) max(v, na.rm = TRUE))
tab$sd_vs_same <- tab$sd / ave(ifelse(tab$proposal == "same draws", tab$sd, NA), tab$target, tab$d, FUN = function(v) max(v, na.rm = TRUE))
emit_table("tab_prop", cbind(ifelse(!duplicated(paste(tab$target, tab$d)), paste(tlab[tab$target], tab$d), ""), lab[tab$proposal],
                             fvse(tab$sd, tab$sd_se, 4), fx(tab$sd_vs_fixed, 3), fx(tab$sd_vs_same, 3),
                             fsvse(tab$bias_in_mcse, tab$bias_se, 2), fx(tab$rmse, 4), fx(tab$coverage, 3),
                             fint(tab$med_n_ov), fx(tab$t_T_rejections_pct, 1)), groups = group_ends(paste(tab$target, tab$d)),
           spec = "llcccccccr",
           header = "Target, $d$ & Proposal & sd & sd/fixed & sd/same & bias/MCSE & RMSE & coverage & med $\\neffhat$ & $t_\\Tp$ rej.\\ (\\%)\\\\")
skew <- tab[tab$target == "skew", ]
tab <- tab[tab$target == "t7", ]
gv <- function(d, arm, col) tab[[col]][tab$d == d & tab$proposal == arm]
ds <- unique(tab$d); dmax <- max(ds)
res_put("prop.S", fint(S))
res_put("prop.dmax", dmax)
res_put("prop.same.absbias.dmax", fx(abs(gv(dmax, "same draws", "bias_in_mcse")), 2))
res_put("prop.same.cov.dmax", fx(gv(dmax, "same draws", "coverage"), 2))
res_put("prop.same.trej.dmax", fx(gv(dmax, "same draws", "t_T_rejections_pct"), 1))
res_put("prop.split.infl", frange(100 * (vapply(ds, function(d) gv(d, "split, same budget", "sd") / gv(d, "same draws", "sd"), 0) - 1), 0))
res_put("prop.cross.cost", fx(max(0, 100 * (vapply(ds, function(d) gv(d, "cross-fitted", "sd") / gv(d, "same draws", "sd"), 0) - 1)), 0))
res_put("prop.cross.cost.fixed", frange(100 * (vapply(ds, function(d) gv(d, "cross-fitted", "sd") / gv(d, "fixed at theta_0", "sd"), 0) - 1), 0))
res_put("prop.split.cost.fixed", frange(100 * (vapply(ds, function(d) gv(d, "split, same budget", "sd") / gv(d, "fixed at theta_0", "sd"), 0) - 1), 0))
for (a in c("same draws", "cross-fitted", "split, same budget")) {
  k <- c("same draws" = "same", "cross-fitted" = "cross", "split, same budget" = "split")[[a]]
  x <- skew[skew$proposal == a, ]
  res_put(sprintf("prop.skew.%s.bias", k), if (nrow(x)) paste(sprintf("$%s$", fsign(x$bias_in_mcse, 2)), collapse = " and ") else "--")
}
res_put("prop.skew.se", if (nrow(skew)) fx(max(skew$bias_se), 2) else "--")
res_put("prop.skew.d", if (nrow(skew)) paste(unique(skew$d), collapse = " and ") else "--")
better <- ds[vapply(ds, function(d) gv(d, "split, same budget", "rmse") < gv(d, "same draws", "rmse"), TRUE)]
res_put("prop.split.rmse.phrase", if (!length(better)) "at none of the dimensions considered" else
  if (length(better) == 1) sprintf("only at $d=%d$", better) else sprintf("at $d\\in\\{%s\\}$", paste(better, collapse = ",")))
best <- all(vapply(ds, function(d) gv(d, "cross-fitted", "rmse") <= min(gv(d, "same draws", "rmse"), gv(d, "split, same budget", "rmse")), TRUE))
res_put("prop.cross.best.phrase", if (best) ", and it has the smallest root mean squared error among the arms that use only the posterior draws at hand" else "")
res_put("prop.sd.dmax.fixed", fx(gv(dmax, "fixed at theta_0", "sd"), 4))
res_put("prop.sd.dmax.same", fx(gv(dmax, "same draws", "sd"), 4))
res_put("prop.sd.dmax.extra", fx(gv(dmax, "split, extra sample", "sd"), 4))
res_put("prop.sdse.dmax", fx(max(tab$sd_se[tab$d == dmax]), 4))
res_put("prop.R.dmax", fint(gv(dmax, "fixed at theta_0", "R")))
res_write("03_proposal_fitting.R")
print(tab, digits = 3)
