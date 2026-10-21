#!/usr/bin/env Rscript
# Supplement Table LRV: long-run variance estimators in use. AR(1) posterior chains
# (gamma = 1 - rho), independent proposals, S1 = S2 = S, 500 replications per row, with
# the sample sizes. Coverage with the oracle sigma_NP, with Geyer's
# estimator, with the autoregressive estimator of coda::effectiveSize, and with the
# formula of bridgesampling::error_measures() as implemented; medians of n_ovP-hat. Bottom:
# coverage pooled over rows and binned by n_ovP-hat. Also: error_measures() for
# independent draws (F2, mu_g = 2, S = 1000).
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(dirname(normalizePath(f))) else if (dir.exists("analysis")) "analysis" else ".."})
for (f in c("common.R", "families.R", "exact.R", "ar1.R", "supplement.R", "paper.R")) source(file.path(BO_ROOT, "R", f))

logZ <- 0.25
R <- n_reps(REPS$lrv)
rows <- rbind(
  data.frame(fam = "F2", rho = 0.5, n_ovP = c(2, 4, 8, 16, 32, 64), S = c(31, 62, 123, 245, 489, 977)),
  data.frame(fam = "F2", rho = 0.9, n_ovP = c(2, 4, 8, 16, 32, 64), S = c(143, 285, 569, 1137, 2273, 4546)),
  data.frame(fam = "F1", rho = 0.5, n_ovP = c(2, 4, 8, 16, 32, 64), S = c(20, 40, 80, 159, 317, 634)),
  data.frame(fam = "F1", rho = 0.9, n_ovP = c(2, 4, 8, 16, 32, 64), S = c(65, 130, 259, 518, 1035, 2069)))
if (BO_QUICK) rows <- rows[c(1, 12, 19), ]
fams <- list(F2 = family_F2(2), F1 = family_F1(10, 1.5))
pops <- lapply(fams, overlap_population)
runs <- do.call(rbind, lapply(seq_len(nrow(rows)), function(i) {
  rw <- rows[i, ]; fam <- fams[[rw$fam]]; pop <- pops[[rw$fam]]; S <- as.integer(rw$S)
  tw <- ar1_tau(fam, pop, rw$rho)
  sigma_NP <- sqrt(S * (pop$var_g + pop$var_p * tw$tau_P)) / (2 * S * pop$H)
  msg("LRV:", fam$label, "rho =", rw$rho, "S =", S)
  res <- run_reps(R, function(r) {
    ap <- logZ + fam$lr(ar1_u(S, fam$d, rw$rho, 1L)[, 1]); aq <- logZ + fam$lr(fam$rprop(S))
    d <- bridge_diagnostics(ap, aq, iid = FALSE, jackknife = FALSE, keep_responsibilities = TRUE)
    tau_ar <- iat(d$responsibilities$posterior, S, "concat_ar")
    n_P_ar <- S * d$V_prop + S * d$V_post * tau_ar
    data.frame(Delta = d$logml - logZ, mcse = d$mcse, mcse_ar = sqrt(n_P_ar) / d$NH0,
               cv_em = sqrt(error_measures_re2(ap, aq, d$logml)), n_ovP_geyer = d$n_ovP, n_ovP_ar = d$n_ov^2 / n_P_ar)
  })
  res$fam <- rw$fam; res$rho <- rw$rho; res$n_ovP <- rw$n_ovP; res$S <- rw$S
  res$tau_W <- tw$tau_W; res$sigma_NP <- sigma_NP
  res
}))
cv <- function(e, s) mean(is.finite(e) & abs(e) <= z975 * s)
tab <- do.call(rbind, lapply(split(runs, interaction(runs$fam, runs$rho, runs$S, drop = TRUE)), function(x) {
  data.frame(family = x$fam[1], rho = x$rho[1], tau_W = x$tau_W[1], n_ovP = x$n_ovP[1], S = x$S[1], R = nrow(x),
             oracle = cv(x$Delta, x$sigma_NP), geyer = cv(x$Delta, x$mcse), ar = cv(x$Delta, x$mcse_ar),
             error_measures = cv(x$Delta, x$cv_em), med_n_ovP_geyer = stats::median(x$n_ovP_geyer),
             med_n_ovP_ar = stats::median(x$n_ovP_ar))
}))
tab <- tab[order(tab$family != "F2", tab$rho, tab$S), ]
save_table(tab, "tab_LRV")
bins <- c(0, 5, 10, 20, 30, 60, Inf)
binned <- function(n, e, s) {
  b <- cut(n, bins, right = FALSE)
  data.frame(bin = levels(b), coverage = as.numeric(tapply(abs(e) <= z975 * s, b, mean)),
             runs = as.integer(table(b)))
}
bottom <- rbind(cbind(estimator = "Geyer", binned(runs$n_ovP_geyer, runs$Delta, runs$mcse)),
                cbind(estimator = "AR", binned(runs$n_ovP_ar, runs$Delta, runs$mcse_ar)))
save_table(bottom, "tab_LRV_binned")
fam <- fams$F2
iid <- run_reps(n_reps(REPS$lrv_iid), function(r) {
  ap <- logZ + fam$lr(fam$rpost(1000)); aq <- logZ + fam$lr(fam$rprop(1000))
  d <- bridge_diagnostics(ap, aq, iid = TRUE, eb = FALSE, jackknife = FALSE)
  data.frame(Delta = d$logml - logZ, mcse = d$mcse, cv_em = sqrt(error_measures_re2(ap, aq, d$logml)), var_f2 = d$V_post / 0.25)
})
extra <- data.frame(case = "independent draws, F2 (mu_g = 2), S = 1000", R = nrow(iid),
                    coverage_delta = cv(iid$Delta, iid$mcse), coverage_error_measures = cv(iid$Delta, iid$cv_em),
                    median_ratio_em_over_mcse = stats::median(iid$cv_em / iid$mcse), median_var_f2 = stats::median(iid$var_f2))
save_table(extra, "tab_LRV_independent")
print(tab, digits = 3); print(bottom, digits = 3); print(extra, digits = 3)
flab <- c(F2 = "F2, $\\mu_g=2$", F1 = "F1")
bin_row <- function(est) {
  x <- bottom[bottom$estimator == est, ]
  sprintf("\\multicolumn{2}{l}{%s, binned} & \\multicolumn{9}{l}{%s}\\\\", est,
          paste(sprintf("%s (%d)", fx(x$coverage, 3), x$runs), collapse = ",\\quad "))
}
emit_table("tab_LRV", cbind(flab[tab$family], tab$rho, fx(tab$tau_W, 2), tab$n_ovP, tab$S, fx(tab$oracle, 3), fx(tab$geyer, 3),
                            fx(tab$ar, 3), fx(tab$error_measures, 3), fx(tab$med_n_ovP_geyer, 1), fx(tab$med_n_ovP_ar, 1)),
           extra = c("\\midrule", bin_row("Geyer"), bin_row("AR")))
pool <- function(lo, hi) {
  vapply(c("Geyer", "AR"), function(est) {
    if (est == "Geyer") { n <- runs$n_ovP_geyer; s <- runs$mcse } else { n <- runs$n_ovP_ar; s <- runs$mcse_ar }
    k <- n >= lo & n < hi
    if (!any(k)) NA_real_ else mean(abs(runs$Delta[k]) <= z975 * s[k])
  }, 0)
}
res_put("lrv.R", fint(R)); res_put("lrv.lt10", frange(pool(0, 10), 2)); res_put("lrv.10to20", frange(pool(10, 20), 2))
res_put("lrv.ge20", frange(c(bottom$coverage[bottom$bin %in% c("[20,30)", "[30,60)", "[60,Inf)")]), 2))
res_put("lrv.ge10", frange(pool(10, Inf), 2))
res_put("lrv.em", frange(tab$error_measures, 2))
res_put("lrv.iid.em", fx(extra$coverage_error_measures, 3)); res_put("lrv.iid.delta", fx(extra$coverage_delta, 3))
res_put("lrv.iid.ratio", fx(extra$median_ratio_em_over_mcse, 2)); res_put("lrv.iid.varf2", fx(extra$median_var_f2, 2))
res_write("supplement/S10_LRV_estimators.R")
