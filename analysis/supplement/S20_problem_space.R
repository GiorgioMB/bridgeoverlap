#!/usr/bin/env Rscript
# Section 6.1: every problem is a law on [0,1] (Lemma "problem space"). With independent draws, the
# bridge output depends on (p, g, dimension) only through the law Q of pi* under m = s1 p + s2 g
# (and on S1, S2), and every Q on [0,1] with mean s2 is a problem. A draw then costs nothing: for a
# law with k atoms, a run is a pair of multinomial count vectors. This script searches the problem
# space for the worst case of the delta-method interval and of the certificate rule:
#   unconditional coverage against the true effective overlap n_ov;
#   P(miss and certified), the quantity that a procedure-level guarantee needs;
#   coverage among certified runs, which no certificate computed from the same draws can guarantee.
# The worst laws are re-evaluated with REPS$ps_refine fresh replications (no selection bias).
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(dirname(normalizePath(f))) else if (dir.exists("analysis")) "analysis" else ".."})
for (f in c("common.R", "paper.R")) source(file.path(BO_ROOT, "R", f))
z <- stats::qnorm(0.975)
CS <- c(1, 2, 3, 5, 7, 10, 15, 20, 30)

ps_run <- function(x, a, w0, w1, S1, S2, R) {
  N <- S1 + S2; s1 <- S1 / N; s2 <- S2 / N; k <- length(x)
  qpost <- pmax(c(w0, a * (1 - x)) / s1, 0); qprop <- pmax(c(a * x, w1) / s2, 0)
  np <- t(stats::rmultinom(R, S1, qpost))[, -1, drop = FALSE]          # pi* = 0 posterior draws carry no overlap
  nq <- t(stats::rmultinom(R, S2, qprop))[, -(k + 1), drop = FALSE]    # pi* = 1 proposal draws carry none either
  odds <- (1 - x) / x
  resp <- function(u) 1 / (1 + outer(exp(-u), odds))                  # pi_u = e^u pi* / (1 - pi* + e^u pi*)
  # safeguarded Newton steps on the monotone estimating function (bisection inside the bracket)
  lo <- rep(-60, R); hi <- rep(60, R); u <- rep(0, R); nt <- nq + np
  for (it in 1:40) {
    pu <- resp(u); Wu <- rowSums(nq * (1 - pu)) - rowSums(np * pu); dW <- -rowSums(nt * pu * (1 - pu))
    lo <- ifelse(Wu > 0, u, lo); hi <- ifelse(Wu > 0, hi, u)
    un <- u - Wu / dW; bad <- !is.finite(un) | un <= lo | un >= hi; un[bad] <- (lo[bad] + hi[bad]) / 2
    if (max(abs(un - u)) < 1e-10) { u <- un; break }
    u <- un }
 ok <- rowSums(np) > 0 & rowSums(nq) > 0 & abs(u) < 59
  pu <- resp(u); NH0 <- rowSums(np * pu)
  vpost <- (rowSums(np * pu^2) - NH0^2 / S1) / (S1 - 1); vprop <- (rowSums(nq * (1 - pu)^2) - rowSums(nq * (1 - pu))^2 / S2) / (S2 - 1)
  nov_hat <- ifelse(ok, S2 * vprop + S1 * vpost, 0); I_hat <- ifelse(ok, rowSums((np + nq) * pu * (1 - pu)) / N / (s1 * s2), 0)
  mcse <- ifelse(ok, sqrt(pmax(nov_hat, 0)) / pmax(NH0, 1e-300), Inf)
  H <- sum(a * x * (1 - x)); I <- H / (s1 * s2)
  list(cover = ok & abs(u) <= z * mcse, nov_hat = nov_hat, I_hat = I_hat, NH0 = ifelse(ok, NH0, 0), mcse = mcse, nov = N * H * (1 - I), I = I, N = N)
}
summarise_run <- function(r) {
  out <- c(cov = mean(r$cover))
  for (sg in c(0.1, 0.3)) {                                                    # the rule "MCSE below sg"
    m <- r$mcse < sg; k <- sub("[.]", "", format(sg))
    out[paste0(c("mpass", "mcovpass", "mjoint"), k)] <- c(mean(m), if (sum(m) >= 400) mean(r$cover[m]) else NA, mean(m & !r$cover))
  }
  for (c0 in CS) {
    p <- r$nov_hat >= c0 * pmin(1, 2 * (1 - r$I_hat))                          # continuous rule of Algorithm 1
    g <- (r$I_hat <= 0.5 & r$nov_hat >= c0) | (r$I_hat > 0.5 & (sqrt(r$NH0) - sqrt(2 / 3 * log(40)))^2 >= 16 * log(40))   # gated rule
    q <- r$nov_hat >= c0                                                         # plain rule (Algorithm 1)
    out[paste0(c("pass", "covpass", "joint", "gjoint", "qpass", "qcovpass", "qjoint"), c0)] <-
      c(mean(p), if (sum(p) >= 400) mean(r$cover[p]) else NA, mean(p & !r$cover), mean(g & !r$cover),
        mean(q), if (sum(q) >= 400) mean(r$cover[q]) else NA, mean(q & !r$cover))
  }
  out
}
## Laws: k <= 3 atoms anywhere, atoms at extreme responsibilities, smooth (discretized Beta) shapes,
## near-perfect proposals, unequal allocations; overlap scaled to a target n_ov in [1, 80].
set.seed(cell_seed("ps laws"))
n_laws <- n_reps(REPS$ps_laws); R <- n_reps(REPS$ps); kinds <- c("atoms", "edge", "mixed", "beta", "near")
laws <- vector("list", n_laws); i <- 0L; tries <- 0L
while (i < n_laws && tries < 50 * n_laws) {
  tries <- tries + 1L; kind <- kinds[i %% 5 + 1]
  if (kind == "beta") { ab <- 10^stats::runif(2, -0.7, 0.7); x <- (1:20 - 0.5) / 20; w <- x^(ab[1] - 1) * (1 - x)^(ab[2] - 1) }
  else { k <- sample(1:3, 1)
    x <- switch(kind, atoms = stats::runif(k, 0.02, 0.98), near = stats::runif(k, 0.35, 0.65),
                edge = { e <- 10^stats::runif(k, -4, log10(0.05)); ifelse(stats::runif(k) < 0.5, e, 1 - e) },
                mixed = { k <- max(k, 2); x <- stats::runif(k, 0.02, 0.98); e <- 10^stats::runif(1, -4, -1.3); x[1] <- if (stats::runif(1) < 0.5) e else 1 - e; x })
    w <- stats::rexp(length(x)) }
  w <- w / sum(w); s2 <- sample(c(0.5, 0.5, 0.5, 0.2, 0.8), 1)
  N <- floor(10^stats::runif(1, log10(60), log10(20000))); S2 <- max(2L, as.integer(round(N * s2))); S1 <- max(2L, N - S2); N <- S1 + S2
  s2 <- S2 / N; s1 <- S1 / N; nov_t <- 10^stats::runif(1, 0, log10(80)); disc <- 1 - 4 * nov_t / (N * s1 * s2)
  if (disc <= 0) next
  H <- s1 * s2 * (1 + (if (kind == "near") 1 else -1) * sqrt(disc)) / 2; atot <- H / sum(w * x * (1 - x))
  w1 <- s2 - atot * sum(w * x); w0 <- s1 - atot * sum(w * (1 - x))
  if (atot > 1 || w1 < -1e-12 || w0 < -1e-12) next
  i <- i + 1L; laws[[i]] <- list(id = i, kind = kind, x = x, a = atot * w, w0 = max(w0, 0), w1 = max(w1, 0), S1 = S1, S2 = S2)
}
laws <- laws[seq_len(i)]
res <- run_reps(length(laws), cell = "ps search", fun = function(j) {
  L <- laws[[j]]; r <- ps_run(L$x, L$a, L$w0, L$w1, L$S1, L$S2, R)
  c(id = j, nov = r$nov, I = r$I, N = r$N, summarise_run(r))
})
res <- as.data.frame(res); res$kind <- vapply(laws, `[[`, "", "kind")[res$id]
## re-evaluate the worst laws with fresh, much longer runs
R2 <- n_reps(REPS$ps_refine)
refine <- function(ids) run_reps(length(ids), cell = "ps refine", fun = function(j) {
  L <- laws[[ids[j]]]; r <- ps_run(L$x, L$a, L$w0, L$w1, L$S1, L$S2, R2); c(id = ids[j], summarise_run(r)) })
w_unc <- utils::head(res$id[res$nov >= 10][order(res$cov[res$nov >= 10])], 10)
w_joint <- utils::head(res$id[order(-res$joint10)], 10); w_gjoint <- utils::head(res$id[order(-res$gjoint10)], 10)
w_cond <- utils::head(res$id[!is.na(res$covpass10)][order(res$covpass10[!is.na(res$covpass10)])], 10)
w_qjoint <- utils::head(res$id[order(-res$qjoint10)], 10); w_mjoint <- utils::head(res$id[order(-res$mjoint03)], 10)
ref <- as.data.frame(refine(unique(c(w_unc, w_joint, w_gjoint, w_cond, w_qjoint, w_mjoint))))
save_table(res, "ps_laws"); save_table(ref, "ps_refined"); saveRDS(laws, file.path(results_dir, "ps_law_defs.rds"))
rr <- function(ids, col) ref[[col]][match(ids, ref$id)]
res_put("ps.laws", fint(nrow(res))); res_put("ps.R", fint(R)); res_put("ps.R.refine", fint(R2))
res_put("ps.unc.worst10", fx(min(rr(w_unc, "cov")), 3)); res_put("ps.unc.worst5", fx(min(res$cov[res$nov >= 5]), 3))
res_put("ps.unc.worst3", fx(min(res$cov[res$nov >= 3]), 3))
res_put("ps.joint.max10", fx(max(rr(w_joint, "joint10")), 3)); res_put("ps.gjoint.max10", fx(max(rr(w_gjoint, "gjoint10")), 3))
res_put("ps.qjoint.max10", fx(max(rr(w_qjoint, "qjoint10")), 3)); res_put("ps.mjoint.max", fx(max(rr(w_mjoint, "mjoint03")), 3))
jm <- w_mjoint[which.max(rr(w_mjoint, "mjoint03"))]
res_put("ps.mjoint.nov", fx(res$nov[res$id == jm], 1)); res_put("ps.mjoint.I", fx(res$I[res$id == jm], 2)); res_put("ps.mjoint.kind", res$kind[res$id == jm])
res_put("ps.mjoint.pass", fx(100 * rr(jm, "mpass03"), 0)); res_put("ps.mjoint.cov", fx(rr(jm, "cov"), 2))
res_put("ps.cond.worst10", if (length(w_cond)) fx(min(rr(w_cond, "covpass10"), na.rm = TRUE), 2) else "--")
j0 <- if (length(w_cond) && any(is.finite(rr(w_cond, "covpass10")))) w_cond[which.min(rr(w_cond, "covpass10"))] else res$id[1]; L0 <- laws[[j0]]
res_put("ps.sel.nov", fx(res$nov[res$id == j0], 0)); res_put("ps.sel.cov", fx(rr(j0, "cov"), 2))
res_put("ps.sel.pass", fx(100 * rr(j0, "pass10"), 0)); res_put("ps.sel.covpass", fx(rr(j0, "covpass10"), 2))
res_put("ps.sel.atom", fx(L0$x[which.max(L0$a)], 3)); res_put("ps.sel.kind", L0$kind)
print(utils::head(res[order(res$cov), c("kind", "nov", "I", "N", "cov")], 5)); print(ref[1:5, ])
res_write("supplement/S20_problem_space.R")
if (have_ggplot()) {
  st <- system2(file.path(R.home("bin"), "Rscript"), shQuote(file.path(BO_ROOT, "plot_problem_space.R")))
  if (!identical(as.integer(st), 0L)) stop("plot_problem_space.R failed (exit status ", st, ")")
}
