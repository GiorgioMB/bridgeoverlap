#!/usr/bin/env Rscript
# Computer-assisted parts of the supplement, re-done in R.
#  1. The bounds on C'_P: the two one-variable inequalities of the proof of Theorem 6,
#     gamma C'_P <= 90 - 1.13 (two-state chain) and gamma^{3/2} C'_P <= 83 y + 11 - 1.09
#     (positive chains, y = gamma^{1/2}), evaluated on 4000 equal cells of [0, 1] in interval
#     arithmetic. R has no directed rounding, so every operation is widened outwards by a
#     few ulps; this emulates outward rounding but is not a certified proof (the stated
#     verification uses outward-rounded interval arithmetic).
#  2. Lemma S3 (Bennett-type inequality for reversible chains), checked exactly on 300
#     random reversible chains with 2 to 5 states.
#  3. The conjectured Bernstein constant b(gamma) = (6 - 6 gamma + gamma^2)/(3 gamma (2 - gamma)):
#     the supremum over 3 <= k <= 3000 is attained at k = 3. Double precision resolves the
#     margin (about gamma^4/288) only for gamma >= 1e-2; smaller gamma use Rmpfr if installed.
#  4. Numerical constants quoted in the text (kappa, Remark on Theorem 5).
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(dirname(normalizePath(f))) else if (dir.exists("analysis")) "analysis" else ".."})
for (f in c("common.R", "families.R", "exact.R", "ar1.R", "supplement.R", "paper.R")) source(file.path(BO_ROOT, "R", f))

## 1. Interval verification of the C'_P bounds ------------------------------------------------
# With Rmpfr, every operation is carried out by MPFR, which rounds correctly (error at most
# half an ulp), and the result is then widened outwards by two ulps at 128 bits; decimal
# constants and cell end points are enclosed the same way. The enclosure is therefore
# rigorous, given the correctness of MPFR. Without Rmpfr the same code runs in double
# precision with an 8-ulp widening, which emulates outward rounding but is not a proof
# (glibc's exp and cosh are not guaranteed to be correctly rounded).
use_mpfr_cp <- requireNamespace("Rmpfr", quietly = TRUE)
num <- if (use_mpfr_cp) function(x) Rmpfr::mpfr(x, 128) else function(x) x
rel <- if (use_mpfr_cp) num(2)^(-126) else 8 * .Machine$double.eps
tiny <- if (use_mpfr_cp) num(2)^(-2000) else 1e-300
iv <- function(lo, hi = lo) {
  lo <- num(lo); hi <- num(hi)
  structure(list(lo = lo - abs(lo) * rel - tiny, hi = hi + abs(hi) * rel + tiny), class = "iv")
}
as_iv <- function(x) if (inherits(x, "iv")) x else iv(x, x)
# elementwise min / max by selection (no arithmetic, hence no rounding)
pmin2 <- function(a, b) { r <- a; k <- which(as.logical(b < a)); if (length(k)) r[k] <- b[k]; r }
pmax2 <- function(a, b) { r <- a; k <- which(as.logical(b > a)); if (length(k)) r[k] <- b[k]; r }
rep_num <- function(x, n) if (length(x) == n) x else num(rep(as.numeric(x), n))
`+.iv` <- function(a, b) { a <- as_iv(a); b <- as_iv(b); iv(a$lo + b$lo, a$hi + b$hi) }
`-.iv` <- function(a, b) { if (missing(b)) return(iv(-a$hi, -a$lo)); a <- as_iv(a); b <- as_iv(b); iv(a$lo - b$hi, a$hi - b$lo) }
`*.iv` <- function(a, b) { a <- as_iv(a); b <- as_iv(b)
  p1 <- a$lo * b$lo; p2 <- a$lo * b$hi; p3 <- a$hi * b$lo; p4 <- a$hi * b$hi
  iv(pmin2(pmin2(p1, p2), pmin2(p3, p4)), pmax2(pmax2(p1, p2), pmax2(p3, p4))) }
`/.iv` <- function(a, b) { a <- as_iv(a); b <- as_iv(b); stopifnot(all(as.logical(b$lo > 0))); a * iv(1 / b$hi, 1 / b$lo) }
iv_sqrt <- function(a) { z <- a$lo; k <- which(as.logical(z < 0)); if (length(k)) z[k] <- 0; iv(sqrt(z), sqrt(a$hi)) }
iv_exp <- function(a) iv(exp(a$lo), exp(a$hi))
iv_max <- function(a, b) { a <- as_iv(a); b <- as_iv(b); n <- max(length(a$lo), length(b$lo))
  iv(pmax2(rep_num(a$lo, n), rep_num(b$lo, n)), pmax2(rep_num(a$hi, n), rep_num(b$hi, n))) }
iv_h <- function(t) { z <- t$lo; k <- which(as.logical(z < 0)); if (length(k)) z[k] <- 0; iv(cosh(z / 2)^2, cosh(t$hi / 2)^2) }
q8 <- 1 - iv_exp(iv(-1) / 8); h8 <- iv_h(iv(1) / 8); e1 <- iv_exp(iv(1))
bound_two_state <- function(g) {           # gamma C'_P for the chain of Theorem 7
  ga <- 1.64 + 2.58 * g / (2 - g)
  1.29 * (3.05 * iv_max(16, 4 * ga) + 2.27 * ga + 2.59 + 0.015 * g) + 1.26 * (2 - g) +
    iv_sqrt((2 - g) * (2 * (2 - g) * h8 + 2 * q8) / (2 * e1 * q8 * q8))
}
bound_positive <- function(y) {            # gamma^{3/2} C'_P for positive chains, y = gamma^{1/2}
  r <- iv_sqrt(2 - y * y)
  tb <- y / (4 * y + 4 * r)
  82.12 * y + 0.02 * y * y * y + 0.63 * (y + r) * (2 - y * y) +
    (4 * y + 4 * r) * iv_exp(tb * 0.5) * iv_sqrt((2 - y * y) * (2 * (2 - y * y) * iv_h(tb) + 2 * tb) / (2 * e1))
}
ncell <- 4000L
cells <- iv(num(0:(ncell - 1)) / ncell, num(1:ncell) / ncell)        # all cells at once
m1 <- min(as.numeric(90 - bound_two_state(cells)$hi))
m2 <- min(as.numeric((83 * cells$lo + 11) - bound_positive(cells)$hi))
cp <- data.frame(check = c("gamma C'_P <= 90 (two-state chain)", "gamma^{3/2} C'_P <= 83 y + 11 (positive chains)"),
                 cells = ncell, arithmetic = if (use_mpfr_cp) "MPFR 128 bits, outward widening" else "double, 8-ulp widening (not certified)",
                 min_margin = c(m1, m2), reported_margin = c(1.13, 1.09), holds = c(m1 > 0, m2 > 0))
save_table(cp, "constants_CP_interval")
print(cp, digits = 4)

## 2. Lemma S3 on random reversible chains -----------------------------------------------------
set.seed(BO_SEED)
lazy_bound <- function(n, s2, c, gam, th) n * s2 / c^2 * (gam^2 * exp(th * c) / (1 - (1 - gam) * exp(th * c)) - gam - th * c)
chk <- do.call(rbind, lapply(seq_len(300), function(r) {
  m <- sample(2:5, 1)
  W <- matrix(stats::rexp(m * m), m, m); W <- W + t(W)
  P <- W / rowSums(W); pi0 <- rowSums(W) / sum(W)
  Dh <- diag(sqrt(pi0)); A <- Dh %*% P %*% solve(Dh)            # symmetric, same spectrum
  ev <- sort(eigen((A + t(A)) / 2, symmetric = TRUE, only.values = TRUE)$values, decreasing = TRUE)
  gam <- min(1, 1 - ev[2])
  u <- stats::rnorm(m); u <- u - sum(pi0 * u)
  c <- max(u); s2 <- sum(pi0 * u^2)
  th_max <- if (gam < 1) log(1 / (1 - gam)) / c else 5 / c
  worst <- -Inf
  for (n in c(1, 2, 5, 20, 100)) for (th in th_max * c(0.05, 0.25, 0.5, 0.75, 0.95)) {
    D <- exp(th * u)
    v <- pi0 * D; lmgf <- 0
    for (j in seq_len(n - 1)) { v <- as.vector(v %*% P) * D; s <- sum(v); lmgf <- lmgf + log(s); v <- v / s }
    lmgf <- lmgf + log(sum(v))
    worst <- max(worst, lmgf - lazy_bound(n, s2, c, gam, th))
  }
  data.frame(chain = r, states = m, gamma = gam, max_lhs_minus_rhs = worst)
}))
lz <- data.frame(chains = nrow(chk), states = paste(range(chk$states), collapse = "-"),
                 max_lhs_minus_rhs = max(chk$max_lhs_minus_rhs), holds = all(chk$max_lhs_minus_rhs <= 1e-10))
save_table(lz, "constants_lazy_bernstein")
print(lz, digits = 4)

## 3. The Bernstein constant b(gamma) -----------------------------------------------------------
# b(gamma) = sup_{k >= 3} {2 E L^k / (k! E L^2)}^{1/(k-2)}, L geometric on {1, 2, ...} with
# success probability gamma. With lambda = 1 - gamma, s0 = log(1/lambda) and the
# Mittag-Leffler expansion, the k-th term is t_k = s0^{-1} {R (1 + c_k)}^{1/(k-2)} with
# R = 2 gamma^3 / {lambda (1 + lambda) s0^3}, c_k = 2 sum_{n >= 1} Re(1 + 2 pi i n / s0)^{-(k+1)}.
# log t_3 - log t_k for k = 4..kmax, computed in the working precision (double or mpfr)
# so that margins far below double precision can be resolved; returns the minimal
# margin, the maximizing k and log t_3 - log b_closed.
b_margins <- function(gam, kmax, use_mpfr = FALSE, bits = 256) {
  num <- if (use_mpfr) function(x) Rmpfr::mpfr(x, bits) else function(x) x
  g <- num(gam); lam <- 1 - g; s0 <- -log1p(-g)
  R <- 2 * g^3 / (lam * (1 + lam) * s0^3)
  logt <- function(kk) {
    nmax <- min(2e5, max(1, ceiling(as.numeric(s0) / (2 * pi) * 10^(30 / (kk + 1)))))
    x <- 2 * pi * num(seq_len(nmax)) / s0
    ck <- 2 * sum((1 + x^2)^(-(kk + 1) / 2) * cos((kk + 1) * atan(x)))
    log(R * (1 + ck)) / (kk - 2) - log(s0)
  }
  l3 <- logt(3L)
  diffs <- vapply(4:kmax, function(kk) as.numeric(l3 - logt(kk)), numeric(1))
  bc <- (6 - 6 * g + g^2) / (3 * g * (2 - g))
  list(min_margin = min(diffs), argmax_k = if (min(diffs) > 0) 3L else which.min(diffs) + 3L,
       closed = as.numeric(l3 - log(bc)), margin_vs_k4 = diffs[1])
}
have_mpfr <- requireNamespace("Rmpfr", quietly = TRUE)
gams <- quick_pick(c(1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 0.03, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99), c(1L, 5L, 8L, 11L))
kmax <- if (BO_QUICK) 200L else 3000L
bt <- do.call(rbind, lapply(gams, function(gm) {
  mp <- gm < 1e-2
  if (mp && !have_mpfr) {
    return(data.frame(gamma = gm, precision = "skipped (install Rmpfr)", kmax = NA, argmax_k = NA,
                      log_t3_minus_closed_form = NA, relative_margin_vs_k4 = NA, gamma4_over_288 = gm^4 / 288,
                      min_margin_over_k = NA))
  }
  km <- kmax
  msg("b(gamma): gamma =", gm, if (mp) "(Rmpfr, 256 bits)" else "(double)", "k <=", km)
  r <- b_margins(gm, km, use_mpfr = mp)
  data.frame(gamma = gm, precision = if (mp) "mpfr 256 bits" else "double", kmax = km, argmax_k = r$argmax_k,
             log_t3_minus_closed_form = r$closed, relative_margin_vs_k4 = expm1(r$margin_vs_k4),
             gamma4_over_288 = gm^4 / 288, min_margin_over_k = r$min_margin)
}))
save_table(bt, "constants_b_gamma")
print(bt, digits = 4)

## 4. Constants quoted in the text --------------------------------------------------------------
tt <- seq(1e-6, 0.5, length.out = 200001)
consts <- data.frame(
  quantity = c("kappa(t)/t^2 as t -> 0", "kappa(1/2)/(1/2)^2", "kappa(1)", "min (1 - e^-t)/t on (0, 1/2]",
               "max (1 - e^-t) on (0, 1/2]", "max cosh^2(t/2) on (0, 1/2]"),
  value = c(kappa_fun(1e-6) / 1e-12, kappa_fun(0.5) / 0.25, kappa_fun(1), min((1 - exp(-tt)) / tt),
            max(1 - exp(-tt)), max(cosh(tt / 2)^2)),
  reported = c(0.5, 0.259, 0.135, 0.7869, 0.3935, 1.0639))
save_table(consts, "constants_text")
print(consts, digits = 5)

## 5. The MBAR overlap matrix for two states (Section S2.8) ----------------------------------------
set.seed(BO_SEED + 5)
fam <- family_F2(2)
mb <- do.call(rbind, lapply(list(c(1000, 1000), c(500, 1500), c(1500, 300)), function(SS) {
  ap <- fam$lr(fam$rpost(SS[1])); aq <- fam$lr(fam$rprop(SS[2]))
  d <- bridge_diagnostics(ap, aq, iid = TRUE, eb = FALSE, jackknife = FALSE, keep_responsibilities = TRUE)
  p <- c(d$responsibilities$posterior, d$responsibilities$proposal)           # proposal responsibilities
  W1 <- (1 - p) / SS[1]; W2 <- p / SS[2]                                      # MBAR weights
  O <- matrix(c(SS[1] * sum(W1 * W1), SS[2] * sum(W1 * W2), SS[1] * sum(W2 * W1), SS[2] * sum(W2 * W2)), 2, byrow = TRUE)
  ev <- sort(eigen(O, only.values = TRUE)$values, decreasing = TRUE)
  data.frame(S1 = SS[1], S2 = SS[2], overlap_scalar = 1 - ev[2], I_hat = d$I, abs_difference = abs(1 - ev[2] - d$I))
}))
save_table(mb, "constants_mbar")
print(mb, digits = 6)

## LaTeX macros ---------------------------------------------------------------------------------
res_put("cp.m1", fx(cp$min_margin[1], 2)); res_put("cp.m2", fx(cp$min_margin[2], 2)); res_put("cp.cells", fint(ncell))
res_put("cp.method", if (use_mpfr_cp) "outward-rounded interval arithmetic (MPFR at 128 bits)" else "interval arithmetic emulated in double precision")
res_put("lazy.n", fint(lz$chains)); res_put("lazy.states", sub("-", "--", lz$states))
res_put("lazy.phrase", if (isTRUE(lz$holds)) "held" else "failed")
res_put("b.kmax", fint(kmax)); res_put("b.gmin", fsci(min(bt$gamma[!is.na(bt$argmax_k)]), 0))
res_put("b.phrase", if (all(bt$argmax_k[!is.na(bt$argmax_k)] == 3)) "the supremum is attained at $k=3$" else "the supremum is not always attained at $k=3$")
small <- bt[!is.na(bt$relative_margin_vs_k4), ]; small <- small[which.min(small$gamma), ]
res_put("b.margin.ratio", fx(small$relative_margin_vs_k4 / small$gamma4_over_288, 3))
res_put("mbar.err", fsci(max(mb$abs_difference), 0))
res_write("supplement/S16_constants.R")
