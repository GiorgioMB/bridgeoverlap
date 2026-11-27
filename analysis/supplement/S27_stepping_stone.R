#!/usr/bin/env Rscript
# Supplement Table S-ss: the optimal bridge against stepping-stone sampling at an equal number of draws.
# Stepping-stone sampling (Xie et al., 2011) estimates log Z along the geometric path
# pi_beta ~ g^(1 - beta) f^beta from g to p, with K stages and N/K draws per stage. For F1 and F2 the path
# is Gaussian, so its draws are exact: the comparison isolates the estimators, at no cost for the K
# samplers that stepping-stone sampling needs in practice. The bridge uses N/2 draws from each of p and g.
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(dirname(normalizePath(f))) else if (dir.exists("analysis")) "analysis" else ".."})
for (f in c("common.R", "paper.R", "families.R")) source(file.path(BO_ROOT, "R", f))
z975 <- stats::qnorm(0.975); logZ <- 0.25; K <- 8L
beta <- ((0:K) / K)^3                                        # schedule concentrated near the proposal
lme <- function(x) { m <- max(x); m + log(mean(exp(x - m))) }
CONF <- list("F2 (mu_g=2)" = list(fam = family_F2(2), path = function(n, b) stats::rnorm(n, (1 - b) * 2)),
             "F2 (mu_g=4)" = list(fam = family_F2(4), path = function(n, b) stats::rnorm(n, (1 - b) * 4)),
             "F2 (mu_g=6)" = list(fam = family_F2(6), path = function(n, b) stats::rnorm(n, (1 - b) * 6)),
             "F1 (d=10, sigma=1.5)" = list(fam = family_F1(10, 1.5), path = function(n, b) stats::rchisq(n, 10) / ((1 - b) / 1.5^2 + b)),
             "F1 (d=40, sigma=1.3)" = list(fam = family_F1(40, 1.3), path = function(n, b) stats::rchisq(n, 40) / ((1 - b) / 1.3^2 + b)))
if (BO_QUICK) CONF <- CONF[c(1, 3, 5)]
NN <- if (BO_QUICK) c(2000L) else c(2000L, 20000L)
R <- n_reps(REPS$ss, key = "ss")
tab <- do.call(rbind, lapply(names(CONF), function(cf) do.call(rbind, lapply(NN, function(N) {
  cc <- CONF[[cf]]; fam <- cc$fam; n <- N %/% K; S <- N %/% 2; msg("ss:", cf, "N =", N)
  out <- run_reps(R, function(r) {
    ap <- logZ + fam$lr(fam$rpost(S)); aq <- logZ + fam$lr(fam$rprop(S))
    d <- bridge_diagnostics(ap, aq, iid = TRUE, interval = FALSE, jackknife = FALSE)
    lw <- lapply(seq_len(K), function(k) (beta[k + 1] - beta[k]) * (logZ + fam$lr(cc$path(n, beta[k]))))
    ss <- sum(vapply(lw, lme, 0))
    ss_se <- sqrt(sum(vapply(lw, function(x) { w <- exp(x - max(x)); stats::var(w) / (length(w) * mean(w)^2) }, 0)))
    data.frame(b_err = d$logml - logZ, b_cov = is.finite(d$logml) && abs(d$logml - logZ) <= z975 * d$mcse,
               b_pass = !isTRUE(d$flags[["overlap"]]), ss_err = ss - logZ, ss_cov = abs(ss - logZ) <= z975 * ss_se)
  }, cell = paste("ss", cf, N))
  fin <- is.finite(out$b_err)
  data.frame(config = cf, N = N, b_rmse = if (any(fin)) sqrt(mean(out$b_err[fin]^2)) else NA, b_fin = 100 * mean(fin),
             b_cov = mean(out$b_cov), b_pass = 100 * mean(out$b_pass), ss_rmse = sqrt(mean(out$ss_err^2)), ss_cov = mean(out$ss_cov))
}))))
tab$ratio <- tab$ss_rmse / tab$b_rmse
save_table(tab, "tab_ss")
lab <- function(x) { x <- gsub("mu_g=([0-9.]+)", "$\\\\mu_g=\\1$", x); x <- gsub("sigma=([0-9.]+)", "$\\\\sigma=\\1$", x); gsub("d=([0-9]+)", "$d=\\1$", x) }
rows <- sprintf("%s & %s & %s & %s & %s & %s & %s & %s\\\\", ifelse(duplicated(tab$config), "", lab(tab$config)), fint(tab$N),
                fsig(tab$b_rmse, 2), fx(tab$b_cov, 3), fx(tab$b_pass, 0), fsig(tab$ss_rmse, 2), fx(tab$ss_cov, 3), fx(tab$ratio, 2))
writeLines(c("\\begin{tabular}{@{}lrcccccc@{}}\\toprule",
             " & & \\multicolumn{3}{c}{optimal bridge} & \\multicolumn{2}{c}{stepping-stone} & \\\\",
             "\\cmidrule(lr){3-5}\\cmidrule(lr){6-7}",
             "Configuration & $N$ & RMSE & coverage & passes (\\%) & RMSE & coverage & RMSE ratio\\\\\\midrule",
             rows, "\\bottomrule\\end{tabular}"), file.path(gen_dir, "tab_ss.tex"))
poor <- tab$config == "F2 (mu_g=6)"                       # the configuration with poor overlap
res_put("ss.R", fint(R)); res_put("ss.K", fint(K))
res_put("ss.ratio.rest", frange(tab$ratio[!poor], 2)); res_put("ss.ratio.poor", frange(tab$ratio[poor], 2))
res_put("ss.pass.poor.small", fx(tab$b_pass[poor & tab$N == min(tab$N)], 0))
res_put("ss.cov", frange(tab$ss_cov, 2))
res_write("supplement/S27_stepping_stone.R")
