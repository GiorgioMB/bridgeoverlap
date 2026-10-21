#!/usr/bin/env Rscript
# Supplement Table R2: Theorem 10 (null calibration of the calibration test) under
# correct sampling, 2000 replications, G = 20. Ratio of sd(T) to sd_{m*}(pi*)/sqrt(N H),
# rejection rate of |t_T| > t_{19, 0.975}, and sd(t_T) (1.057 for t_19).
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(dirname(normalizePath(f))) else if (dir.exists("analysis")) "analysis" else ".."})
for (f in c("common.R", "families.R", "exact.R", "ar1.R", "supplement.R", "paper.R")) source(file.path(BO_ROOT, "R", f))

R <- n_reps(REPS$r2)
G <- 20L
dn <- function(m, s) list(l = function(x) stats::dnorm(x, m, s, log = TRUE), r = function(n) stats::rnorm(n, m, s))
dc <- list(l = function(x) stats::dcauchy(x, log = TRUE), r = function(n) stats::rcauchy(n))
cfgs <- list(
  list(name = "N(0,1), N(2,1)", p = dn(0, 1), g = dn(2, 1), S1 = 2000, S2 = 2000),
  list(name = "N(0,1), N(2,1)", p = dn(0, 1), g = dn(2, 1), S1 = 1000, S2 = 3000),
  list(name = "N(0,1), N(4,1)", p = dn(0, 1), g = dn(4, 1), S1 = 2000, S2 = 2000),
  list(name = "Cauchy, N(0,2^2)", p = dc, g = dn(0, 2), S1 = 2000, S2 = 2000))
breaks <- c(-60, -20, -8, -4, -2, -1, 0, 1, 2, 3, 4, 6, 8, 20, 60)
tab <- do.call(rbind, lapply(quick_pick(cfgs, c(1L, 4L)), function(cc) {
  ov <- overlap_1d(cc$p$l, cc$g$l, cc$S1, cc$S2, breaks)
  ms <- mstar_moments(ov)
  msg("R2:", cc$name, cc$S1, cc$S2)
  res <- run_reps(R, function(r) {
    x <- cc$p$r(cc$S1); y <- cc$g$r(cc$S2)
    d <- bridge_diagnostics(cc$p$l(x) - cc$g$l(x), cc$p$l(y) - cc$g$l(y), iid = TRUE, eb = FALSE, G = G)
    data.frame(T = d$T, t_T = d$t_T)
  })
  data.frame(setting = cc$name, S1 = cc$S1, S2 = cc$S2, R = nrow(res), n_ov = ov$n_ov,
             sd_ratio = stats::sd(res$T) / (sqrt(ms$var) / sqrt(ov$NH)),
             rejection = mean(abs(res$t_T) > stats::qt(0.975, G - 1), na.rm = TRUE),
             sd_t_T = stats::sd(res$t_T, na.rm = TRUE))
}))
save_table(tab, "tab_R2")
print(tab, digits = 3)
lab <- function(s) gsub("N\\(([^)]*)\\)", "$\\\\mathcal N(\\1)$", s)
emit_table("tab_R2", cbind(lab(tab$setting), tab$S1, tab$S2, fx(tab$n_ov, 1), fx(tab$sd_ratio, 2), fx(tab$rejection, 3), fx(tab$sd_t_T, 2)))
res_put("r2.R", fint(R)); res_put("r2.sddev", fx(100 * max(abs(tab$sd_ratio - 1)), 0))
res_put("r2.rej", frange(tab$rejection, 3)); res_put("r2.sdt", frange(tab$sd_t_T, 2))
res_write("supplement/S14_R2_calibration_null.R")
