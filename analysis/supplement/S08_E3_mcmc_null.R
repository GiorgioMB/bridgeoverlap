#!/usr/bin/env Rscript
# Supplement Table E3: Markov chain output. Stationary AR(1) posterior chains with
# coefficient rho in {0.5, 0.9} and marginal N(0, I_5), S1 = S2 = 8000, G = 20, 5000
# replications, proposal g = N(0.2 * 1, 1.3^2 I_5). Left: spectral MCSE (Geyer's initial
# monotone sequence); right: naive MCSE that ignores autocorrelation. Rejections at 5%.
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(dirname(normalizePath(f))) else if (dir.exists("analysis")) "analysis" else ".."})
for (f in c("common.R", "families.R", "exact.R", "ar1.R", "supplement.R", "paper.R")) source(file.path(BO_ROOT, "R", f))

logZ <- 0.15
R <- n_reps(REPS$e3)
S <- 8000L
G <- 20L
invisible(jackknife_null_quantiles(G, probs = 0.95))
gm <- rep(0.2, 5); gv <- rep(1.3^2, 5)
thr_v <- stats::qchisq(0.95, G - 1) / (G - 1)
tab <- do.call(rbind, lapply(quick_pick(c(0.5, 0.9), 1L), function(rho) {
  msg("E3: rho =", rho)
  res <- run_reps(R, function(r) {
    X <- ar1_matrix(S, 5L, rho)
    Y <- sweep(matrix(stats::rnorm(S * 5), S, 5) * 1.3, 2, gm, "+")
    d <- bridge_diagnostics(logZ + lmvn_std(X) - ldiag(X, gm, gv), logZ + lmvn_std(Y) - ldiag(Y, gm, gv),
                            iid = FALSE, G = G, eb = FALSE)
    rj <- rejections(d)
    data.frame(Delta = d$logml - logZ, mcse = d$mcse, mcse_naive = d$mcse_iid, V_jack = d$V_jack,
               V_jack_naive = d$V_jack * (d$mcse / d$mcse_iid)^2, rV = rj[["V_jack"]],
               rI = rj[["I_max_exact"]], rG = rj[["G_max"]], rT = rj[["t_T"]])
  })
  data.frame(rho = rho, R = nrow(res),
             var_over_mcse2 = stats::var(res$Delta) / mean(res$mcse^2),
             coverage = mean(abs(res$Delta) <= z975 * res$mcse), mean_V_jack = mean(res$V_jack),
             V_jack = 100 * mean(res$rV), I_max = 100 * mean(res$rI), G_max = 100 * mean(res$rG),
             t_T = 100 * mean(res$rT),
             naive_var_over_mcse2 = stats::var(res$Delta) / mean(res$mcse_naive^2),
             naive_coverage = mean(abs(res$Delta) <= z975 * res$mcse_naive),
             naive_mean_V_jack = mean(res$V_jack_naive), naive_V_jack = 100 * mean(res$V_jack_naive > thr_v))
}))
save_table(tab, "tab_E3")
print(tab, digits = 3)
emit_table("tab_E3", cbind(tab$rho, fx(tab$var_over_mcse2, 3), fx(tab$coverage, 3), fx(tab$mean_V_jack, 3), fx(tab$V_jack, 1),
                           fx(tab$I_max, 1), fx(tab$G_max, 1), fx(tab$t_T, 1), fx(tab$naive_var_over_mcse2, 3),
                           fx(tab$naive_coverage, 3), fx(tab$naive_mean_V_jack, 3), fx(tab$naive_V_jack, 1)))
res_put("e3.R", fint(R))
for (rho in tab$rho) {
  x <- tab[tab$rho == rho, ]
  res_put(sprintf("e3.naive.cov.%g", rho), fx(x$naive_coverage, 2)); res_put(sprintf("e3.naive.vj.%g", rho), fx(x$naive_V_jack, 1))
  res_put(sprintf("e3.naive.factor.%g", rho), fx(x$naive_var_over_mcse2, 2))
}
res_write("supplement/S08_E3_mcmc_null.R")
