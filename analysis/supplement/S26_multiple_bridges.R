#!/usr/bin/env Rscript
# Supplement Section S-multi (Table S-mbar): multiple bridge sampling (MBAR) and the effective
# resistance. K = 5 unnormalized Gaussians q_k = Z_k N(2(k - 1), 1), n draws from each. The ratio
# Z_5 / Z_1 is estimated by MBAR from all draws and by the two-sample optimal bridge from the draws
# of the end states only. Theorem S-multi says that the clock of MBAR is N / R_15, the number of
# draws divided by the effective resistance between the end states in the overlap graph, whose
# edge weights are the pairwise overlaps E_m(pi_k pi_l); the two-sample bridge needs N of order 1 / H_15.
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(dirname(normalizePath(f))) else if (dir.exists("analysis")) "analysis" else ".."})
for (f in c("common.R", "paper.R")) source(file.path(BO_ROOT, "R", f))
z975 <- stats::qnorm(0.975)
K <- 5L; mu <- 2 * (0:(K - 1)); logZ <- 0:(K - 1)                    # arbitrary known normalizing constants
truth <- logZ[K] - logZ[1]
lq <- function(x) outer(x, mu, function(x, m) stats::dnorm(x, m, log = TRUE)) + rep(logZ, each = length(x))
## population overlap graph (equal allocations): edge weights w_kl = E_m(pi_k pi_l), Laplacian, resistances
pk <- function(x) exp(outer(x, mu, function(x, m) stats::dnorm(x, m, log = TRUE)))
wkl <- function(k, l) stats::integrate(function(x) { P <- pk(x) / K; m <- rowSums(P); P[, k] * P[, l] / m },
                                       min(mu) - 12, max(mu) + 12, subdivisions = 2000L, rel.tol = 1e-10)$value
Wg <- matrix(0, K, K); for (k in 1:(K - 1)) for (l in (k + 1):K) Wg[k, l] <- Wg[l, k] <- wkl(k, l)
Lap <- diag(rowSums(Wg)) - Wg
res_eff <- function(k, l) { e <- numeric(K); e[k] <- 1; e[l] <- -1; Lg <- Lap[-1, -1]; v <- e[-1]; drop(t(v) %*% solve(Lg, v)) }
R15 <- res_eff(1, K); Rseries <- sum(1 / Wg[cbind(1:(K - 1), 2:K)])
H15 <- stats::integrate(function(x) { p1 <- stats::dnorm(x, mu[1]); p5 <- stats::dnorm(x, mu[K]); 0.25 * p1 * p5 / ((p1 + p5) / 2) },
                        mu[1] - 12, mu[K] + 12, subdivisions = 2000L, rel.tol = 1e-12)$value
## MBAR: maximum-likelihood intercepts of the multinomial logistic regression, by safeguarded Newton steps
mbar <- function(LQ, lab) {
  N <- nrow(LQ); nk <- tabulate(lab, K); ls <- log(nk / N); cc <- rep(0, K)
  state <- function(cc) {
    A <- LQ + rep(ls - cc, each = N); mx <- do.call(pmax, as.data.frame(A)); lse <- mx + log(rowSums(exp(A - mx)))
    P <- exp(A - lse); list(L = sum(A[cbind(seq_len(N), lab)] - lse), P = P, W = colSums(P) - nk)
  }
  st <- state(cc)
  for (it in 1:100) {
    J <- diag(colSums(st$P)) - crossprod(st$P)
    step <- c(0, solve(J[-1, -1], st$W[-1])); tt <- 1
    repeat { new <- state(cc + tt * step); if (new$L >= st$L - 1e-10 || tt < 1e-8) break; tt <- tt / 2 }
    cc <- cc + tt * step; st <- new
    if (max(abs(st$W)) < 1e-9 * N) break
  }
  J <- diag(colSums(st$P)) - crossprod(st$P)
  V <- Reduce(`+`, lapply(1:K, function(k) { Pk <- st$P[lab == k, , drop = FALSE]; nrow(Pk) * stats::cov.wt(Pk, method = "ML")$cov }))
  Ji <- solve(J[-1, -1]); S <- Ji %*% V[-1, -1] %*% Ji               # c_1 = 0: the contrast c_5 - c_1 is c_5
  list(est = cc[K] - cc[1], se = sqrt(S[K - 1, K - 1]), n_a = S[K - 1, K - 1] / Ji[K - 1, K - 1]^2, NoverR = 1 / Ji[K - 1, K - 1])
}
NN <- if (BO_QUICK) c(50L, 500L) else c(20L, 50L, 200L, 1000L)
R <- n_reps(REPS$mbar, key = "mbar")
tab <- do.call(rbind, lapply(NN, function(n) {
  msg("mbar: n =", n)
  out <- run_reps(R, function(r) {
    th <- unlist(lapply(mu, function(m) stats::rnorm(n, m))); lab <- rep(1:K, each = n); LQ <- lq(th)
    mb <- mbar(LQ, lab)
    i1 <- lab == 1; i5 <- lab == K                                     # two-sample bridge: "posterior" p_5, "proposal" q_1
    d <- bridge_diagnostics(LQ[i5, K] - LQ[i5, 1], LQ[i1, K] - LQ[i1, 1], iid = TRUE, interval = FALSE, jackknife = FALSE)
    data.frame(mb_err = mb$est - truth, mb_cov = abs(mb$est - truth) <= z975 * mb$se, mb_na = mb$n_a, mb_NR = mb$NoverR,
               d_err = d$logml - truth, d_cov = is.finite(d$logml) && abs(d$logml - truth) <= z975 * d$mcse, d_nov = d$n_ov)
  }, cell = paste("mbar", n))
  data.frame(n = n, N = K * n, NoverR = K * n / R15, mb_rmse = sqrt(mean(out$mb_err^2)), mb_cov = mean(out$mb_cov),
             mb_na = stats::median(out$mb_na), d_fin = mean(is.finite(out$d_err)), d_cov = mean(out$d_cov), d_nov = stats::median(out$d_nov),
             mb_scaled = sqrt(mean(out$mb_err^2)) * sqrt(K * n / R15))
}))
save_table(tab, "tab_mbar")
rows <- sprintf("%s & %s & %s & %s & %s & %s & %s & %s\\\\", fint(tab$n), fx(tab$NoverR, 1), fx(tab$mb_rmse, 3), fx(tab$mb_cov, 3),
                fx(tab$mb_na, 1), fx(100 * tab$d_fin, 0), fx(tab$d_cov, 3), fx(tab$d_nov, 2))
writeLines(c("\\begin{tabular}{@{}rrcccrcc@{}}\\toprule",
             " & & \\multicolumn{3}{c}{MBAR, all draws} & \\multicolumn{3}{c}{two-sample bridge, end states}\\\\",
             "\\cmidrule(lr){3-5}\\cmidrule(l){6-8}",
             "$n$ & $N/R_{15}$ & RMSE & coverage & med $\\widehat n_a$ & finite (\\%) & coverage & med $\\neffhat$\\\\\\midrule",
             rows, "\\bottomrule\\end{tabular}"), file.path(gen_dir, "tab_mbar.tex"))
res_put("mb.R", fint(R)); res_put("mb.R15", fx(R15, 1)); res_put("mb.Rseries", fx(Rseries, 1))
res_put("mb.H15", fsci(H15, 1)); res_put("mb.scaled.lo", fx(min(tab$mb_scaled), 2)); res_put("mb.scaled.hi", fx(max(tab$mb_scaled), 2))
res_put("mb.cov.big", fx(tab$mb_cov[nrow(tab)], 3)); res_put("mb.n.big", fint(tab$n[nrow(tab)])); res_put("mb.dcov.big", fx(tab$d_cov[nrow(tab)], 3))
res_write("supplement/S26_multiple_bridges.R")
