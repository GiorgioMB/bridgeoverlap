#!/usr/bin/env Rscript
# Supplement Table REGEN: the jackknife variance compared with batch-means and
# regeneration-based standard errors. Two-point family of Theorem 2 (eps = 0.05) with the
# two-state posterior chain of Theorem 7, which has an atom (state 0), so that regenerative
# standard errors (Mykland et al. 1995; Doss and Tan 2014) are available exactly. For each
# estimator of the posterior long-run variance -- Geyer's initial monotone sequence,
# batch means with 20 batches and with batch size sqrt(S1), regeneration at visits to the
# atom -- and for the jackknife standard error itself (t_19 quantiles), the table reports
# the median of the standard error relative to the true sd of log Z-hat and the coverage
# of nominal 95% intervals.
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(dirname(normalizePath(f))) else if (dir.exists("analysis")) "analysis" else ".."})
for (f in c("common.R", "paper.R")) source(file.path(BO_ROOT, "R", f))
eps <- 0.05
G <- 20L
two_state_path <- function(S, gam) {           # stationary start; sojourns are geometric
  out <- integer(0); state <- as.integer(stats::runif(1) < eps)
  while (length(out) < S) {
    len <- stats::rgeom(1, if (state == 0L) eps * gam else (1 - eps) * gam) + 1L
    out <- c(out, rep(state, len)); state <- 1L - state
  }
  out[seq_len(S)]
}
bm_lrv <- function(y, nb) { b <- length(y) %/% nb; m <- colMeans(matrix(y[seq_len(nb * b)], b)); b * stats::var(m) }
regen_lrv <- function(y, x) {                  # tours start at each entry into the atom {0}
  starts <- which(x == 0L & c(TRUE, x[-length(x)] != 0L))
  if (length(starts) < 3L) return(NA_real_)
  keep <- starts[1]:(starts[length(starts)] - 1L)
  id <- cumsum(seq_along(x) %in% starts)[keep]
  Ysum <- tapply(y[keep], id, sum); Tl <- tapply(y[keep] * 0 + 1, id, sum)
  mu <- sum(Ysum) / sum(Tl)
  sum((Ysum - mu * Tl)^2) / sum(Tl)
}
cfg <- expand.grid(gamma = c(0.5, 0.1), n_ovP = c(10, 30, 100))
cfg$S <- round(2 * cfg$n_ovP / (cfg$gamma * eps * (1 - eps)))
if (BO_QUICK) cfg <- cfg[c(1, 6), ]
z <- stats::qnorm(0.975); tq <- stats::qt(0.975, G - 1)
rows <- list()
for (i in seq_len(nrow(cfg))) {
  cc <- cfg[i, ]; S <- cc$S
  msg("REGEN: gamma =", cc$gamma, "S =", S)
  res <- run_reps(n_reps(REPS$regen), cell = paste("regen", cc$gamma, S), fun = function(r) {
    x <- two_state_path(S, cc$gamma)
    yq <- sample(c(1L, 2L), S, replace = TRUE, prob = c(eps, 1 - eps))
    ap <- ifelse(x == 0L, Inf, 0); aq <- ifelse(yq == 1L, 0, -Inf)   # log Z = 0
    d <- bridge_diagnostics(ap, aq, chains = 1, G = G, keep_responsibilities = TRUE, nsim = 2e4)
    if (!is.finite(d$logml)) return(data.frame(err = d$logml, se_geyer = NA, se_bm20 = NA, se_bmsq = NA, se_regen = NA, se_jack = NA))
    pp <- d$responsibilities$posterior
    se_of <- function(lrv) sqrt(S * d$V_prop + S * lrv) / d$NH0
    data.frame(err = d$logml, se_geyer = d$mcse, se_bm20 = se_of(bm_lrv(pp, G)),
               se_bmsq = se_of(bm_lrv(pp, floor(sqrt(S)))), se_regen = se_of(regen_lrv(pp, x)),
               se_jack = sqrt(d$V_jack) * d$mcse)
  })
  ok <- is.finite(res$err)
  sdt <- stats::sd(res$err[ok])
  for (est in c("geyer", "bm20", "bmsq", "regen", "jack")) {
    se <- res[[paste0("se_", est)]]; q <- if (est %in% c("bm20", "jack")) tq else z
    rows[[length(rows) + 1L]] <- data.frame(gamma = cc$gamma, n_ovP = cc$n_ovP, S = S, R = sum(ok), estimator = est,
                                            med_ratio = stats::median(se[ok] / sdt, na.rm = TRUE),
                                            coverage = mean(ok & abs(res$err) <= q * se, na.rm = TRUE))
  }
}
tab <- do.call(rbind, rows)
save_table(tab, "tab_regen")
w <- reshape(tab[, c("gamma", "n_ovP", "S", "estimator", "med_ratio", "coverage")], idvar = c("gamma", "n_ovP", "S"),
             timevar = "estimator", direction = "wide")
ests <- c("geyer", "bm20", "bmsq", "regen", "jack")
emit_table("tab_regen", cbind(fx(w$gamma, 1), fint(w$n_ovP), fint(w$S),
                              do.call(cbind, lapply(ests, function(e) sprintf("%s (%s)", fx(w[[paste0("coverage.", e)]], 3), fx(w[[paste0("med_ratio.", e)]], 2))))),
           spec = "rrrccccc",
           header = "$\\gamma$ & $\\neffP$ & $S$ & Geyer & BM, 20 batches & BM, $\\sqrt{S}$ & regeneration & jackknife ($t_{19}$)\\\\")
res_put("regen.R", fint(n_reps(REPS$regen)))
for (e in ests) res_put(paste0("regen.cov.", e), frange(tab$coverage[tab$estimator == e], 3))
res_write("supplement/S18_LRV_regeneration.R")
