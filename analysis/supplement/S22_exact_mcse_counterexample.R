#!/usr/bin/env Rscript
# Section 7.2: an exact counterexample to the rule "MCSE below 0.3" under poor overlap.
# A law Q with three atoms (Lemma "problem space"): responsibilities xA (posterior bulk), xA2 (a small
# posterior component, so that the MCSE is small but not zero) and xB (proposal bulk). A run is a pair
# of multinomial count vectors; every probability below is a finite sum over counts, computed exactly
# on a window whose omitted mass is reported (no simulation).
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(dirname(normalizePath(f))) else if (dir.exists("analysis")) "analysis" else ".."})
for (f in c("common.R", "paper.R")) source(file.path(BO_ROOT, "R", f))
z <- stats::qnorm(0.975)
exact_law <- function(xA, xA2, xB, mu2, S) {
  a2 <- (mu2 / S) * 0.5 / (1 - xA2)                      # posterior P(A2) = mu2 / S
  aB <- (0.5 - a2 * xA2 - (1 - a2) * xA) / (xB - xA); aA <- 1 - a2 - aB
  x <- c(xA, xA2, xB)
  pp <- c(aA * (1 - xA), a2 * (1 - xA2), aB * (1 - xB)) / 0.5; pq <- c(aA * xA, a2 * xA2, aB * xB) / 0.5
  H <- sum(c(aA, a2, aB) * x * (1 - x)); I <- H / 0.25; nov <- 2 * S * H * (1 - I)
  win <- function(p) { m <- S * p; sd <- sqrt(S * p) + 1; seq(max(0, floor(m - 7 * sd - 8)), ceiling(m + 7 * sd + 14)) }
  g <- expand.grid(k2 = win(pp[2]), kB = win(pp[3]), qA = win(pq[1]), q2 = win(pq[2]))
  g <- g[g$k2 + g$kB <= S & g$qA + g$q2 <= S, ]
  NP <- cbind(S - g$k2 - g$kB, g$k2, g$kB); NQ <- cbind(g$qA, g$q2, S - g$qA - g$q2)
  lw <- lfactorial(S) - rowSums(lfactorial(NP)) + NP %*% log(pp) + lfactorial(S) - rowSums(lfactorial(NQ)) + NQ %*% log(pq)
  w <- exp(drop(lw))
  resp <- function(u) { e <- exp(u); cbind(e * x[1] / (1 - x[1] + e * x[1]), e * x[2] / (1 - x[2] + e * x[2]), e * x[3] / (1 - x[3] + e * x[3])) }
  lo <- rep(-60, nrow(g)); hi <- rep(60, nrow(g))
  for (it in 1:70) { u <- (lo + hi) / 2; R <- resp(u); Wu <- rowSums(NQ * (1 - R)) - rowSums(NP * R); lo <- ifelse(Wu > 0, u, lo); hi <- ifelse(Wu > 0, hi, u) }
  u <- (lo + hi) / 2; R <- resp(u)
  NH0 <- rowSums(NP * R); Vp <- (rowSums(NP * R^2) - NH0^2 / S) / (S - 1)
  Qr <- 1 - R; sq <- rowSums(NQ * Qr); Vq <- (rowSums(NQ * Qr^2) - sq^2 / S) / (S - 1)
  nov_hat <- S * pmax(Vp, 0) + S * pmax(Vq, 0); mcse <- sqrt(nov_hat) / NH0
  NH <- rowSums((NP + NQ) * R * (1 - R)); Tcal <- 1 - NH / NH0
  cover <- abs(u) <= z * mcse; pass <- mcse < 0.3; miss_pass <- pass & !cover
  o <- order(mcse[miss_pass]); cw <- cumsum(w[miss_pass][o])
  list(nov = nov, I = I, omitted = 1 - sum(w), cov = sum(w * cover), mpass = sum(w * pass), mmiss = sum(w * miss_pass),
       med_mcse = mcse[miss_pass][o][which(cw >= cw[length(cw)] / 2)[1]], cpass = sum(w * (nov_hat >= 10)),
       med_T = stats::weighted.mean(Tcal[miss_pass], w[miss_pass]))
}
r <- exact_law(1e-3, 3e-3, 1 - 1e-4, 5, 1000L)
print(unlist(r), digits = 4)
res_put("ex.S", fint(1000)); res_put("ex.nov", fx(r$nov, 1)); res_put("ex.I", fx(r$I, 3)); res_put("ex.omitted", fsci(max(abs(r$omitted), 1e-12), 0))
res_put("ex.cov", fx(r$cov, 3)); res_put("ex.mpass", fx(100 * r$mpass, 1)); res_put("ex.mmiss", fx(100 * r$mmiss, 1))
res_put("ex.medmcse", fx(r$med_mcse, 4)); res_put("ex.cpass", fx(100 * r$cpass, 0)); res_put("ex.cpass.p", fsci(r$cpass, 0)); res_put("ex.T", fx(r$med_T, 2))
res_write("supplement/S22_exact_mcse_counterexample.R")
