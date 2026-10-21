#!/usr/bin/env Rscript
# Supplement Table R1: Proposition 1 under poor overlap (10^4 replications per row).
# Largest ratio n_ov-hat / (varsigma max{M_post, M_prop}) (at most 1 by part (a)); frequency
# of n_ov-hat >= 10 and its bound from part (b); 99% quantile of n_ov-hat and the bound of
# part (b) with delta = 0.01. M_post and M_prop use the population responsibilities.
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(dirname(normalizePath(f))) else if (dir.exists("analysis")) "analysis" else ".."})
for (f in c("common.R", "families.R", "exact.R", "ar1.R", "supplement.R", "paper.R")) source(file.path(BO_ROOT, "R", f))

R <- n_reps(REPS$r1)
dn <- function(m, s) list(l = function(x) stats::dnorm(x, m, s, log = TRUE), r = function(n) stats::rnorm(n, m, s))
dc <- list(l = function(x) stats::dcauchy(x, log = TRUE), r = function(n) stats::rcauchy(n))
q_bound <- function(mu, x) if (x >= mu) min(1, exp(log(2) - mu + x * (1 + log(mu) - log(x)))) else 1
cfgs <- list(
  list(name = "N(0,1), N(5,1)", p = dn(0, 1), g = dn(5, 1), S1 = 1000, S2 = 1000),
  list(name = "N(0,1), N(6,1)", p = dn(0, 1), g = dn(6, 1), S1 = 1000, S2 = 1000),
  list(name = "N(0,1), N(7,1)", p = dn(0, 1), g = dn(7, 1), S1 = 1000, S2 = 1000),
  list(name = "N(0,1), N(6,1)", p = dn(0, 1), g = dn(6, 1), S1 = 400, S2 = 1600),
  list(name = "Cauchy, N(8,0.5^2)", p = dc, g = dn(8, 0.5), S1 = 150, S2 = 150),
  list(name = "N(0,1), N(0,0.05^2)", p = dn(0, 1), g = dn(0, 0.05), S1 = 40, S2 = 40),
  list(name = "N(0,1), N(0,10^2)", p = dn(0, 1), g = dn(0, 10), S1 = 30, S2 = 30))
breaks <- c(-60, -20, -8, -4, -2, -1, -0.2, 0, 0.2, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 20, 60)
tab <- do.call(rbind, lapply(quick_pick(cfgs, c(1L, 5L)), function(cc) {
  ov <- overlap_1d(cc$p$l, cc$g$l, cc$S1, cc$S2, breaks)
  vs <- cc$S1 / (cc$S1 - 1) + cc$S2 / (cc$S2 - 1)
  msg("R1:", cc$name, cc$S1, cc$S2)
  res <- run_reps(R, function(r) {
    x <- cc$p$r(cc$S1); y <- cc$g$r(cc$S2)
    d <- bridge_diagnostics(cc$p$l(x) - cc$g$l(x), cc$p$l(y) - cc$g$l(y), iid = TRUE, eb = FALSE, jackknife = FALSE)
    M <- max(sum(ov$pistar(x)), sum(1 - ov$pistar(y)))
    data.frame(n_ov = d$n_ov, ratio = d$n_ov / (vs * M))
  })
  data.frame(setting = cc$name, S1 = cc$S1, S2 = cc$S2, R = nrow(res), n_ov = ov$n_ov, NH = ov$NH,
             max_ratio = max(res$ratio), p_n_ov_ge_10 = mean(res$n_ov >= 10), bound = q_bound(ov$NH, 10 / vs),
             q99_n_ov = stats::quantile(res$n_ov, 0.99, names = FALSE),
             bound_q99 = vs * (sqrt(ov$NH) + sqrt(2 / 3 * log(2 / 0.01)))^2)
}))
save_table(tab, "tab_R1")
print(tab, digits = 3)
lab <- function(s) { s <- gsub("N\\(([^)]*)\\)", "$\\\\mathcal N(\\1)$", s); s }
fb <- function(x) ifelse(x >= 1, "1", ifelse(x >= 0.01, fx(x, 2), fsig(x, 1)))
emit_table("tab_R1", cbind(lab(tab$setting), tab$S1, tab$S2, fx(tab$n_ov, 2), fx(tab$NH, 2), fx(tab$max_ratio, 2),
                           ifelse(tab$p_n_ov_ge_10 == 0, "0", fx(tab$p_n_ov_ge_10, 3)), fb(tab$bound),
                           fx(tab$q99_n_ov, 1), fx(tab$bound_q99, 1)))
res_put("r1.R", fint(R)); res_put("r1.maxratio", fx(max(tab$max_ratio), 2))
small <- tab$n_ov <= 3.2
res_put("r1.never.phrase", if (all(tab$p_n_ov_ge_10[small] == 0)) "$\\neffhat\\ge10$ never occurred when $\\neff\\le3.2$" else
  sprintf("$\\neffhat\\ge10$ occurred in at most a proportion %s of the replications with $\\neff\\le3.2$", fx(max(tab$p_n_ov_ge_10[small]), 4)))
res_write("supplement/S13_R1_certificate_bounds.R")
