#!/usr/bin/env Rscript
# Supplement Section S-space (Figure S-space): calibration of the overlap check over the problem space.
# With independent draws the bridge output depends on (p, g, dimension) only through the law Q of pi*
# under m = s1 p + s2 g, and on (N1, N2) (Lemma S-space); every law Q on [0,1] with mean s2 is a
# problem. For a law with k atoms a run is a pair of multinomial count vectors, so a search over laws
# covers problems of every dimension at once. For each law we compute the coverage of the delta-method
# interval and, for the overlap check at threshold n0 and for the rule "MCSE below 0.3", the probability
# that a run passes and its interval misses. Coverage is a function of the clock n_ov (Theorem 3), so the
# threshold of Algorithm 1 is calibrated on the clock: n0 is the smallest grid value such that the worst
# coverage over all searched laws with n_ov >= n0 is at least 1 - alpha - 0.01 = 0.94. The worst laws are
# re-evaluated with fresh runs, so the reported worst cases carry no selection bias.
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(dirname(normalizePath(f))) else if (dir.exists("analysis")) "analysis" else ".."})
for (f in c("common.R", "paper.R")) source(file.path(BO_ROOT, "R", f))
z <- stats::qnorm(0.975)
CS <- c(1, 2, 3, 5, 7, 10, 15, 20, 30)
CRIT <- 0.94                                    # worst-case coverage required on the clock at the calibrated threshold

ps_run <- function(x, a, w0, w1, S1, S2, R) {
  N <- S1 + S2; s1 <- S1 / N; s2 <- S2 / N; k <- length(x)
  qpost <- pmax(c(w0, a * (1 - x)) / s1, 0); qprop <- pmax(c(a * x, w1) / s2, 0)
  np <- t(stats::rmultinom(R, S1, qpost))[, -1, drop = FALSE]          # pi* = 0 posterior draws carry no overlap
  nq <- t(stats::rmultinom(R, S2, qprop))[, -(k + 1), drop = FALSE]    # pi* = 1 proposal draws carry none either
  odds <- (1 - x) / x
  resp <- function(u) 1 / (1 + outer(exp(-u), odds))                  # pi_u = e^u pi* / (1 - pi* + e^u pi*)
  lo <- rep(-60, R); hi <- rep(60, R); u <- rep(0, R); nt <- nq + np   # safeguarded Newton steps
  for (it in 1:40) {
    pu <- resp(u); Wu <- rowSums(nq * (1 - pu)) - rowSums(np * pu); dW <- -rowSums(nt * pu * (1 - pu))
    lo <- ifelse(Wu > 0, u, lo); hi <- ifelse(Wu > 0, hi, u)
    un <- u - Wu / dW; bad <- !is.finite(un) | un <= lo | un >= hi; un[bad] <- (lo[bad] + hi[bad]) / 2
    if (max(abs(un - u)) < 1e-10) { u <- un; break }
    u <- un
  }
  ok <- rowSums(np) > 0 & rowSums(nq) > 0 & abs(u) < 59
  pu <- resp(u); NH0 <- rowSums(np * pu)
  vpost <- (rowSums(np * pu^2) - NH0^2 / S1) / (S1 - 1); vprop <- (rowSums(nq * (1 - pu)^2) - rowSums(nq * (1 - pu))^2 / S2) / (S2 - 1)
  nov_hat <- ifelse(ok, S2 * vprop + S1 * vpost, 0)
  mcse <- ifelse(ok, sqrt(pmax(nov_hat, 0)) / pmax(NH0, 1e-300), Inf)
  H <- sum(a * x * (1 - x)); I <- H / (s1 * s2)
  list(cover = ok & abs(u) <= z * mcse, nov_hat = nov_hat, mcse = mcse, nov = N * H * (1 - I), I = I, N = N)
}
summarise_run <- function(r) {
  m <- r$mcse < 0.3
  out <- c(cov = mean(r$cover), mpass = mean(m), mjoint = mean(m & !r$cover))
  for (c0 in CS) {
    q <- r$nov_hat >= c0
    out[paste0(c("qpass", "qcovpass", "qjoint"), c0)] <- c(mean(q), if (sum(q) >= 400) mean(r$cover[q]) else NA, mean(q & !r$cover))
  }
  out
}
## Laws: 1-3 atoms anywhere, atoms at extreme responsibilities, mixed, smooth (discretized Beta) and
## near-perfect laws, unequal allocations; overlap scaled to a target n_ov in [1, 80].
set.seed(cell_seed("ps laws"))
n_laws <- n_reps(REPS$ps_laws, key = "ps_laws"); R <- n_reps(REPS$ps, key = "ps"); kinds <- c("atoms", "edge", "mixed", "beta", "near")
laws <- vector("list", n_laws); i <- 0L; tries <- 0L
while (i < n_laws && tries < 50 * n_laws) {
  tries <- tries + 1L; kind <- kinds[i %% 5 + 1]
  if (kind == "beta") { ab <- 10^stats::runif(2, -0.7, 0.7); x <- (1:20 - 0.5) / 20; w <- x^(ab[1] - 1) * (1 - x)^(ab[2] - 1) }
  else {
    k <- sample(1:3, 1)
    x <- switch(kind, atoms = stats::runif(k, 0.02, 0.98), near = stats::runif(k, 0.35, 0.65),
                edge = { e <- 10^stats::runif(k, -4, log10(0.05)); ifelse(stats::runif(k) < 0.5, e, 1 - e) },
                mixed = { k <- max(k, 2); x <- stats::runif(k, 0.02, 0.98); e <- 10^stats::runif(1, -4, -1.3); x[1] <- if (stats::runif(1) < 0.5) e else 1 - e; x })
    w <- stats::rexp(length(x))
  }
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
res <- as.data.frame(run_reps(length(laws), cell = "ps search", fun = function(j) {
  L <- laws[[j]]; r <- ps_run(L$x, L$a, L$w0, L$w1, L$S1, L$S2, R)
  c(id = j, nov = r$nov, I = r$I, N = r$N, summarise_run(r))
}))
res$kind <- vapply(laws, `[[`, "", "kind")[res$id]
## re-evaluate the worst laws with fresh, much longer runs
R2 <- n_reps(REPS$ps_refine, key = "ps_refine")
top <- function(col, n = 5) utils::head(res$id[order(-res[[col]])], n)
w_q <- unlist(lapply(CS, function(c0) top(paste0("qjoint", c0))))
w_unc <- unlist(lapply(CS, function(c0) { k <- res$nov >= c0; utils::head(res$id[k][order(res$cov[k])], 5) }))
w_m <- top("mjoint", 10)
w_cond <- utils::head(res$id[!is.na(res$qcovpass10)][order(res$qcovpass10[!is.na(res$qcovpass10)])], 10)
ids <- unique(c(w_q, w_unc, w_m, w_cond))
ref <- as.data.frame(run_reps(length(ids), cell = "ps refine", fun = function(j) {
  L <- laws[[ids[j]]]; r <- ps_run(L$x, L$a, L$w0, L$w1, L$S1, L$S2, R2); c(id = ids[j], summarise_run(r))
}))
worst <- vapply(CS, function(c0) max(ref[[paste0("qjoint", c0)]], na.rm = TRUE), 0)
nov_of <- res$nov[match(ref$id, res$id)]
worst_cov <- vapply(CS, function(c0) { k <- nov_of >= c0; if (any(k)) min(ref$cov[k]) else NA_real_ }, 0)
curve <- data.frame(n0 = CS, worst_cov = worst_cov, worst_joint = worst, worst_mcse_rule = max(ref$mjoint, na.rm = TRUE))
n0 <- CS[which(worst_cov >= CRIT)[1]]
save_table(res, "ps_laws"); save_table(ref, "ps_refined"); save_table(curve, "ps_curve"); saveRDS(laws, file.path(results_dir, "ps_law_defs.rds"))
rr <- function(id, col) ref[[col]][match(id, ref$id)]
res_put("ps.laws", fint(nrow(res))); res_put("ps.R", fint(R)); res_put("ps.R.refine", fint(R2))
res_put("ps.crit", fx(CRIT, 2)); res_put("ps.nclock", if (length(n0) && !is.na(n0)) fint(n0) else "--")
res_put("ps.qjoint.max10", fx(worst[CS == 10], 3)); res_put("ps.qjoint.max5", fx(worst[CS == 5], 3)); res_put("ps.qjoint.max20", fx(worst[CS == 20], 3))
res_put("ps.mjoint.max", fx(max(ref$mjoint, na.rm = TRUE), 3))
for (c0 in c(3, 5, 7, 10, 15)) res_put(sprintf("ps.unc.worst%d", c0), fx(worst_cov[CS == c0], 3))
jc <- if (length(w_cond)) w_cond[which.min(rr(w_cond, "qcovpass10"))] else NA
res_put("ps.cond.worst10", if (!is.na(jc)) fx(rr(jc, "qcovpass10"), 2) else "--")
res_put("ps.sel.nov", if (!is.na(jc)) fx(res$nov[res$id == jc], 0) else "--"); res_put("ps.sel.cov", if (!is.na(jc)) fx(rr(jc, "cov"), 2) else "--")
res_put("ps.sel.pass", if (!is.na(jc)) fx(100 * rr(jc, "qpass10"), 0) else "--"); res_put("ps.sel.atom", if (!is.na(jc)) fx(laws[[jc]]$x[which.max(laws[[jc]]$a)], 3) else "--")
print(curve)
res_write("supplement/S20_problem_space.R")
if (have_ggplot() && requireNamespace("patchwork", quietly = TRUE)) {
  st <- system2(file.path(R.home("bin"), "Rscript"), shQuote(file.path(BO_ROOT, "plot_problem_space.R")))
  if (!identical(as.integer(st), 0L)) stop("plot_problem_space.R failed (exit status ", st, ")")
}
