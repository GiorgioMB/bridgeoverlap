#!/usr/bin/env Rscript
# Supplement Section S5 (Proposition S-fsmcmc): the finite-sample interval for Markov chain output,
# which needs a lower bound on the spectral gap but no thinning, against the interval of
# Theorem 7 (transfer by thinning). Posterior draws: stationary Gaussian AR(1) chains with
# coefficient rho (spectral gap 1 - rho) targeting N(0, 1); proposal N(1, 1) (family F2).
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(dirname(normalizePath(f))) else if (dir.exists("analysis")) "analysis" else ".."})
for (f in c("common.R", "paper.R", "families.R")) source(file.path(BO_ROOT, "R", f))
fam <- family_F2(1)
ar1 <- function(n, rho) as.numeric(stats::filter(sqrt(1 - rho^2) * stats::rnorm(n), rho, method = "recursive", init = stats::rnorm(1)))
thin_lag <- function(n, rho, budget = 0.01) {      # smallest L with (n/L - 1) beta(L) <= budget (Corollary S-coupling, d = 1)
  for (L in seq_len(n)) if ((floor(n / L) - 1) * sqrt(-0.25 * log1p(-rho^(2 * L))) <= budget) return(L)
  n
}
RHO <- if (BO_QUICK) c(0.8) else c(0.5, 0.8, 0.95)
SS <- if (BO_QUICK) c(1000L) else c(1000L, 4000L, 16000L)
R <- n_reps(REPS$mcfs, key = "mcfs")
cells <- do.call(rbind, lapply(RHO, function(rho) do.call(rbind, lapply(SS, function(S) {
  L <- thin_lag(S, rho); msg("mcfs: rho =", rho, "S =", S, "thinning lag", L)
  out <- run_reps(R, function(r) {
    ap <- fam$lr(ar1(S, rho)); aq <- fam$lr(fam$rprop(S))           # log Z = 0
    m <- bridge_fs_interval_mcmc(ap, aq, gamma = 1 - rho, delta = 0.05)
    t <- bridge_fs_interval(ap[seq(1L, S, by = L)], aq, delta = 0.04)  # total error 0.04 + 0.01
    data.frame(chat = bridge_logml(ap, aq)$logml, m_lo = m[["lower"]], m_hi = m[["upper"]], t_lo = t[["lower"]], t_hi = t[["upper"]])
  }, cell = paste("mcfs", rho, S))
  oracle <- 2 * stats::qnorm(0.975) * stats::sd(out$chat)            # width of the interval with the true standard deviation
  data.frame(rho = rho, S = S, L = L, cov_mcmc = mean(out$m_lo <= 0 & 0 <= out$m_hi), cov_thin = mean(out$t_lo <= 0 & 0 <= out$t_hi),
             w_mcmc = stats::median(out$m_hi - out$m_lo) / oracle, w_thin = stats::median(out$t_hi - out$t_lo) / oracle)
}))))
save_table(cells, "mcmc_interval")
res_put("mcfs.R", fint(R)); res_put("mcfs.cov.min", fx(min(cells$cov_mcmc), 3))
res_put("mcfs.w.lo", fx(min(cells$w_mcmc), 1)); res_put("mcfs.w.hi", fx(max(cells$w_mcmc), 1))
big <- cells$S == max(cells$S)
res_put("mcfs.wbig.lo", fx(min(cells$w_mcmc[big]), 1)); res_put("mcfs.wbig.hi", fx(max(cells$w_mcmc[big]), 1)); res_put("mcfs.Sbig", fint(max(cells$S)))
res_put("mcfs.thin.w.lo", fx(min(cells$w_thin[is.finite(cells$w_thin)]), 1)); res_put("mcfs.thin.w.hi", if (all(is.finite(cells$w_thin))) fx(max(cells$w_thin), 1) else "$\\infty$")
res_put("mcfs.thin.L.lo", fint(min(cells$L))); res_put("mcfs.thin.L.hi", fint(max(cells$L)))
res_write("supplement/S25_mcmc_interval.R")
