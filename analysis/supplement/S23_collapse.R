#!/usr/bin/env Rscript
# Supplement Section S5.1 (Table S-collapse): why coverage under MCMC is a function of n^P_ov.
#  (i)  The two-state chain of Theorem S-mcmc-sharp in the small-overlap limit of
#       Proposition S-collapse, computed without Monte Carlo error: the posterior count is
#       compound Poisson with geometric sojourns, the proposal count Poisson, and both laws are
#       obtained by inverting probability generating functions with the FFT.
#  (ii) The relaxation ratio r_P = 1/(gamma tau_W) and the constant zeta_W for stationary
#       AR(1) chains, from Hermite (F2) and Laguerre (F1, d = 5) expansions.
#  (iii) The multi-scale example of Remark S-multiscale (Monte Carlo, REPS$collapse_ms draws
#       per value of kappa).
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(dirname(normalizePath(f))) else if (dir.exists("analysis")) "analysis" else ".."})
for (f in c("common.R", "paper.R")) source(file.path(BO_ROOT, "R", f))
z975 <- stats::qnorm(0.975)
XS <- seq(-5, 5, length.out = if (BO_QUICK) 401L else 1601L)
## (i) two-state chain, small-overlap limit ----------------------------------------------
pmf_fft <- function(pgf, M) pmax(Re(stats::fft(pgf(exp(2i * pi * (0:(M - 1)) / M)))) / M, 0)
cov_and_D <- function(p_post, cdf_prop, p0, sig) {     # Delta = log(n_prop / n_post)
  m <- seq_along(p_post)[-1] - 1; pm <- p_post[-1]; keep <- pm > 1e-16; m <- m[keep]; pm <- pm[keep]
  K <- length(cdf_prop) - 1; at <- function(k) cdf_prop[pmin(k, K) + 1]    # cdf_prop[k + 1] = P(n_prop <= k)
  h <- z975 * sig
  cov <- sum(pm * (at(floor(m * exp(h))) - at(ceiling(m * exp(-h)) - 1)))
  D <- 0
  for (xc in split(XS, ceiling(seq_along(XS) / 100))) {
    y <- outer(exp(xc * sig), m)
    Fv <- p0 + as.vector((matrix(at(floor(y)), nrow(y)) - cdf_prop[1]) %*% pm)
    Fm <- p0 + as.vector((matrix(at(pmax(ceiling(y) - 1, 0)), nrow(y)) - cdf_prop[1]) %*% pm)
    D <- max(D, abs(Fv - stats::pnorm(xc)), abs(Fm - stats::pnorm(xc)))
  }
  c(cov = cov, D = D)
}
limit_model <- function(nP, gam) {                     # lambda = 2 n^P_ov, sigma_{N,P}^2 = 2 / lambda
  lam <- 2 * nP; sig <- sqrt(2 / lam)
  if (gam == 0) {                                      # gamma -> 0: Delta = log(lambda / Y), Y ~ CP(lambda, Exp(1))
    kk <- seq_len(ceiling(lam + 15 * sqrt(lam) + 30)); w <- stats::dpois(kk, lam)
    FY <- function(y) exp(-lam) + vapply(y, function(v) sum(w * stats::pgamma(v, kk)), 0)
    cov <- FY(lam * exp(z975 * sig)) - FY(lam * exp(-z975 * sig))
    return(c(cov = cov, D = max(abs(1 - FY(lam * exp(-XS * sig)) - stats::pnorm(XS)))))
  }
  mu <- lam / gam; M <- 2^ceiling(log2(mu + 14 * sqrt(lam * (2 - gam)) / gam + 64))
  G <- function(z) gam * z / (1 - (1 - gam) * z)
  p_post <- pmf_fft(function(z) exp(lam * (G(z) - 1)), M)
  cov_and_D(p_post, stats::ppois(0:(ceiling(M * exp(5 * sig)) + 10), mu), stats::dpois(0, mu), sig)
}
NPS <- if (BO_QUICK) c(2, 16) else c(1, 2, 4, 8, 16, 32)
GF <- if (BO_QUICK) c(1, 0.5, 0.05, 0) else c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01, 0)
SHOW <- c(1, 0.5, 0.2, 0.05, 0)
lim <- lapply(NPS, function(nP) { r <- t(vapply(GF, function(g) limit_model(nP, g), c(cov = 0, D = 0))); rownames(r) <- GF; r })
names(lim) <- NPS
rows <- vapply(seq_along(NPS), function(i) {
  r <- lim[[i]]; nP <- NPS[i]; show <- vapply(SHOW, function(g) { j <- match(g, GF); if (is.na(j)) NA_real_ else r[j, "cov"] }, 0)
  paste(c(nP, fx(show, 3), fx(max(r[, "cov"]) - min(r[, "cov"]), 3), frange(sqrt(nP) * r[, "D"], 2)), collapse = " & ")
}, "")
writeLines(c("\\begin{tabular}{@{}rcccccrc@{}}\\toprule",
             " & \\multicolumn{5}{c}{coverage of nominal $95\\%$ intervals, $\\gamma=$} & spread over & $\\sqrt{\\neffP}\\times$\\\\",
             "$\\neffP$ & $1$ & $0.5$ & $0.2$ & $0.05$ & $\\to0$ & $\\gamma\\in[0,1]$ & Kolmogorov\\\\\\midrule",
             paste0(rows, "\\\\"), "\\bottomrule\\end{tabular}"), file.path(gen_dir, "tab_collapse.tex"))
for (nP in c(1, 2, 4, 8, 16, 32)) {
  r <- lim[[as.character(nP)]]
  res_put(sprintf("col.lim.spread%d", nP), if (is.null(r)) NA else fx(max(r[, "cov"]) - min(r[, "cov"]), 3))
}
ks <- unlist(lapply(seq_along(NPS), function(i) sqrt(NPS[i]) * lim[[i]][, "D"]))
res_put("col.lim.kslo", fx(min(ks), 2)); res_put("col.lim.kshi", fx(max(ks), 2))
r16 <- lim[["16"]]
res_put("col.lim.cov16lo", fx(min(r16[, "cov"]), 3)); res_put("col.lim.cov16hi", fx(max(r16[, "cov"]), 3))
res_put("col.lim.ks16lo", fx(min(4 * r16[, "D"]), 2)); res_put("col.lim.ks16hi", fx(max(4 * r16[, "D"]), 2))
## (ii) AR(1) chains: r_P and zeta_W ------------------------------------------------------
golub_welsch <- function(a, b) {                       # nodes and weights of the probability measure
  n <- length(a); J <- diag(a, n); J[cbind(1:(n - 1), 2:n)] <- b; J[cbind(2:n, 1:(n - 1))] <- b
  e <- eigen(J, symmetric = TRUE); list(x = e$values, w = e$vectors[1, ]^2)
}
gh <- golub_welsch(rep(0, 180), sqrt(1:179))          # N(0, 1)
alpha <- 1.5                                          # d = 5: |theta|^2 / 2 ~ Gamma(alpha + 1, 1)
gl <- golub_welsch(2 * (0:159) + alpha + 1, sqrt((1:159) * ((1:159) + alpha)))
herm <- function(u, K = 90) {                         # orthonormal Hermite coefficients under N(0, 1)
  h0 <- rep(1, length(gh$x)); h1 <- gh$x; cf <- c(sum(gh$w * u * h0), sum(gh$w * u * h1))
  for (k in 1:(K - 1)) { h2 <- (gh$x * h1 - sqrt(k) * h0) / sqrt(k + 1); h0 <- h1; h1 <- h2; cf <- c(cf, sum(gh$w * u * h1)) }
  cf
}
lag <- function(u, K = 60) {                          # orthonormal Laguerre coefficients under Gamma(alpha + 1, 1)
  y <- gl$x; L0 <- rep(1, length(y)); L1 <- 1 + alpha - y; cf <- c(sum(gl$w * u), sum(gl$w * u * L1) / sqrt(alpha + 1))
  for (k in 1:(K - 1)) {
    L2 <- ((2 * k + 1 + alpha - y) * L1 - (k + alpha) * L0) / (k + 1); L0 <- L1; L1 <- L2
    cf <- c(cf, sum(gl$w * u * L1) / exp(0.5 * (lgamma(k + 2 + alpha) - lgamma(k + 2) - lgamma(alpha + 1))))
  }
  cf
}
vlrv <- function(cf, eig) { k <- seq_along(cf)[-1] - 1; c(var = sum(cf[-1]^2), lrv = sum(cf[-1]^2 * (1 + eig(k)) / (1 - eig(k)))) }
vg_of <- function(w, u) sum(w * u^2) - sum(w * u)^2
F2 <- function(m, rho) {                              # p = N(0,1), g = N(m,1), s1 = s2; AR(1) posterior chain
  parts <- function(t) {
    pp <- stats::plogis(m * gh$x - m^2 / 2 + t); pg <- stats::plogis(m * gh$x + m^2 / 2 + t)
    a <- vlrv(herm(pp), function(k) rho^k); b <- vlrv(herm(pp * (1 - pp)), function(k) rho^k)
    c(vg = vg_of(gh$w, pg), vp = a[["var"]], lp = a[["lrv"]], vwg = vg_of(gh$w, pg * (1 - pg)), lw = b[["lrv"]])
  }
  p0 <- parts(0); tauW <- (p0[["vg"]] + p0[["lp"]]) / (p0[["vg"]] + p0[["vp"]])
  z <- max(vapply(seq(-0.25, 0.25, length.out = 11), function(t) { q <- parts(t); sqrt((q[["vwg"]] + q[["lw"]]) / (q[["vg"]] + q[["lp"]])) }, 0))
  c(r = 1 / ((1 - rho) * tauW), z = z)
}
F1 <- function(s2, rho, d = 5) {                      # p = N(0, I_5), g = N(0, s2 I_5); radial factor chain, eigenvalues rho^(2k)
  a1 <- 1 - 1 / s2; off <- -(d / 2) * log(s2)
  parts <- function(t) {
    pp <- stats::plogis(off + a1 * gl$x + t); pg <- stats::plogis(off + a1 * s2 * gl$x + t)
    a <- vlrv(lag(pp), function(k) rho^(2 * k)); b <- vlrv(lag(pp * (1 - pp)), function(k) rho^(2 * k))
    c(vg = vg_of(gl$w, pg), vp = a[["var"]], lp = a[["lrv"]], vwg = vg_of(gl$w, pg * (1 - pg)), lw = b[["lrv"]])
  }
  p0 <- parts(0); tauW <- (p0[["vg"]] + p0[["lp"]]) / (p0[["vg"]] + p0[["vp"]])
  z <- max(vapply(seq(-0.25, 0.25, length.out = 11), function(t) { q <- parts(t); sqrt((q[["vwg"]] + q[["lw"]]) / (q[["vg"]] + q[["lp"]])) }, 0))
  c(r = 1 / ((1 - rho) * tauW), rphi = 1 / ((1 - rho^2) * tauW), z = z)
}
RHO <- c(0.5, 0.8, 0.95, 0.99)
f2 <- do.call(rbind, lapply(c(1, 2, 3), function(m) t(vapply(RHO, function(r) F2(m, r), c(r = 0, z = 0)))))
f1 <- do.call(rbind, lapply(c(1.5, 2.5), function(s) t(vapply(RHO, function(r) F1(s, r), c(r = 0, rphi = 0, z = 0)))))
res_put("col.ar.F2.r", fx(max(f2[, "r"]), 2)); res_put("col.ar.F1.r", fx(max(f1[, "r"]), 2))
res_put("col.ar.F1.rphi", fx(max(f1[, "rphi"]), 2)); res_put("col.ar.z", fx(max(c(f2[, "z"], f1[, "z"])), 2))
res_put("col.ar.zP", fx(sqrt((1 + max(RHO)) / (1 - max(RHO))), 1))
## (iii) a weakly seen slow mode ------------------------------------------------------------
occupation <- function(kappa, n) {                    # fraction of [0, 1] spent in A by a stationary telegraph process
  rate <- kappa / 2; t <- numeric(n); inA <- stats::runif(n) < 0.5; occ <- numeric(n); alive <- rep(TRUE, n)
  while (any(alive)) {
    end <- pmin(t + stats::rexp(n, rate), 1)
    occ <- occ + ifelse(alive & inA, end - t, 0); t <- ifelse(alive, end, t)
    alive <- alive & t < 1; inA <- !inA
  }
  occ
}
Rm <- n_reps(REPS$collapse_ms, key = "collapse_ms")
set.seed(cell_seed("collapse_ms"))
nP <- 16; tau <- 3; mu <- 2 * tau * nP; sg <- sqrt(2 * tau / mu)
for (i in seq_along(KAPPA <- c(0.25, 1, 4, 16))) {
  kappa <- KAPPA[i]; dl <- sqrt((tau - 1) * kappa / mu); U <- 2 * occupation(kappa, Rm) - 1
  a <- stats::rpois(Rm, mu * (1 + dl * U)); b <- stats::rpois(Rm, mu)
  D <- ifelse(a > 0 & b > 0, log(pmax(b, 1) / pmax(a, 1)), ifelse(b == 0, -Inf, Inf)) / sg
  x <- sort(D); Fe <- seq_len(Rm) / Rm
  ks <- max(abs(Fe - stats::pnorm(x)), abs(Fe - 1 / Rm - stats::pnorm(x)))
  res_put(sprintf("col.ms.cov%d", i), fx(mean(abs(D) <= z975), 3)); res_put(sprintf("col.ms.ks%d", i), fx(sqrt(nP) * ks, 2))
}
res_write("supplement/S23_collapse.R")
