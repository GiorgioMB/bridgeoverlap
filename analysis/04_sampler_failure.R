#!/usr/bin/env Rscript
# Section 7.4 and Supplement Table E5: sampler failure and metastability (5000
# replications per row; rejection rates in % at nominal 5%).
# Target 0.5 N(-3 e1, I) + 0.5 N(3 e1, diag(1/4, 1)); proposal N(0, diag(10, 1)).
# (a) independent draws that put weight w on the narrow component (S = 4000);
# (b) a metastable chain whose stationary law is the target: label switching with
#     probability p_sw per step and AR(1) dependence within components (S = 8000,
#     spectral MCSE). The within-component AR coefficient is a choice made
#     here: rho = 0.5.
# Expected: (a) bias 2.42/4.80/9.72 MCSEs for w = 0.6/0.7/0.9 with t_T rejections
# 45.9/96.9/100%; (b) Var/MCSE^2 1.11/1.63/4.72, mean V_jack 1.04/0.64/0.63, t_T 5.5/23.5/59.8%.
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else if (dir.exists("analysis")) "analysis" else "."})
for (f in c("common.R", "supplement.R", "paper.R")) source(file.path(BO_ROOT, "R", f))

R <- n_reps(REPS$e5)
logZ <- 0.3
rho_within <- 0.5
invisible(jackknife_null_quantiles(20L, probs = 0.95))
log_f <- function(X) logZ + lse2(log(0.5) + stats::dnorm(X[, 1], -3, 1, log = TRUE),
                                log(0.5) + stats::dnorm(X[, 1], 3, 0.5, log = TRUE)) +
  stats::dnorm(X[, 2], log = TRUE)
log_g <- function(X) stats::dnorm(X[, 1], 0, sqrt(10), log = TRUE) + stats::dnorm(X[, 2], log = TRUE)
r_g <- function(n) cbind(stats::rnorm(n, 0, sqrt(10)), stats::rnorm(n))
place <- function(z, U) cbind(ifelse(z == 1L, 3 + 0.5 * U[, 1], -3 + U[, 1]), U[, 2])
one <- function(X, Y, iid) {
  d <- bridge_diagnostics(log_f(X) - log_g(X), log_f(Y) - log_g(Y), iid = iid, eb = FALSE)
  rj <- rejections(d)
  data.frame(n_ov = d$n_ov, D = d$logml - logZ, mcse = d$mcse, V_jack = d$V_jack, rV = rj[["V_jack"]],
             rI = rj[["I_max_exact"]], rG = rj[["G_max"]], rT = rj[["t_T"]])
}
summ <- function(res, scenario) {
  data.frame(scenario = scenario, R = nrow(res), med_n_ov = stats::median(res$n_ov),
             bias_in_mcse = mean(res$D) / stats::median(res$mcse),
             var_over_mcse2 = stats::var(res$D) / mean(res$mcse^2),
             coverage = mean(abs(res$D) <= z975 * res$mcse), mean_V_jack = mean(res$V_jack, na.rm = TRUE),
             V_jack = 100 * mean(res$rV), I_max = 100 * mean(res$rI), G_max = 100 * mean(res$rG),
             t_T = 100 * mean(res$rT))
}
S <- 4000L
tab_a <- do.call(rbind, lapply(quick_pick(c(0.5, 0.6, 0.7, 0.9)), function(w) {
  msg("(a) weight on the narrow component w =", w)
  res <- run_reps(R, function(r) one(place(stats::rbinom(S, 1, w), matrix(stats::rnorm(2 * S), S, 2)), r_g(S), TRUE))
  cbind(summ(res, sprintf("(a) imbalance, $w=%g$%s", w, if (w == 0.5) " (correct)" else "")), w = w, p_sw = NA)
}))
S <- 8000L
tab_b <- do.call(rbind, lapply(quick_pick(c(1e-2, 1e-3, 2e-4)), function(p_sw) {
  msg("(b) switching probability", p_sw)
  res <- run_reps(R, function(r) {
    z <- (stats::rbinom(1, 1, 0.5) + c(0L, cumsum(stats::rbinom(S - 1L, 1, p_sw)))) %% 2L
    E <- matrix(stats::rnorm(2 * S), S, 2)
    E[-1L, ] <- sqrt(1 - rho_within^2) * E[-1L, ]
    one(place(z, unclass(stats::filter(E, rho_within, method = "recursive"))), r_g(S), FALSE)
  })
  cbind(summ(res, sprintf("(b) metastable, $p_{\\mathrm{sw}}=%s$", format(p_sw, scientific = FALSE))), w = NA, p_sw = p_sw)
}))
tab <- rbind(tab_a, tab_b)
save_table(tab, "tab_E5")
emit_table("tab_E5", cbind(tab$scenario, fint(tab$med_n_ov), fsign(tab$bias_in_mcse, 2), fx(tab$var_over_mcse2, 2),
                           fx(tab$coverage, 3), fx(tab$mean_V_jack, 3), fx(tab$V_jack, 1), fx(tab$I_max, 1),
                           fx(tab$G_max, 1), fx(tab$t_T, 1)), groups = nrow(tab_a))
a <- tab[!is.na(tab$w) & tab$w != 0.5, ]; b <- tab[!is.na(tab$p_sw), ]
res_put("e5.R", fint(R)); res_put("e5.rho", rho_within)
res_put("e5.a.bias", frange(a$bias_in_mcse, 1))
for (w in c(0.6, 0.7, 0.9)) res_put(sprintf("e5.a.t.%g", w), fx(tab$t_T[!is.na(tab$w) & tab$w == w], 0))
res_put("e5.b.var.max", fx(max(b$var_over_mcse2), 1))
res_put("e5.b.vj.slow", frange(b$mean_V_jack[b$p_sw <= 0.001], 2))
for (p_sw in c(0.001, 2e-4)) res_put(sprintf("e5.b.t.%s", format(p_sw, scientific = FALSE)), fx(b$t_T[b$p_sw == p_sw], 0))
res_write("04_sampler_failure.R")
print(tab, digits = 3)
