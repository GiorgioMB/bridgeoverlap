#!/usr/bin/env Rscript
# Supplement Table ADV: adversarial continuous examples (S1 = S2 = S). Smooth analogue of
# the two-point family with p/g constant (mu0 = 0) or varying (mu0 = 0.5) on the overlap;
# overlap carried by a rare mode; funnels with a full-covariance Gaussian proposal
# fitted on an independent sample in every replication. True n_ov by quadrature where
# available, medians of n_ov-hat and I-hat, coverage, KS distance of the standardized
# error Delta / MCSE to N(0,1), and coverage and boundedness of the Prop. 2 interval.
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(dirname(normalizePath(f))) else if (dir.exists("analysis")) "analysis" else ".."})
for (f in c("common.R", "families.R", "exact.R", "ar1.R", "supplement.R", "paper.R")) source(file.path(BO_ROOT, "R", f))

logZ <- 0.4
eps <- 0.02
dn <- function(x, m, s = 1) stats::dnorm(x, m, s, log = TRUE)
mix2 <- function(w, m1, s1, m2, s2) list(
  ld = function(x) lse2(log(w) + dn(x, m1, s1), log(1 - w) + dn(x, m2, s2)),
  r = function(n) { k <- stats::runif(n) < w; ifelse(k, stats::rnorm(n, m1, s1), stats::rnorm(n, m2, s2)) })
# The example definitions are also used by S17_ADV_replication_check.R.
STREAM <- if (exists("ADV_STREAM")) ADV_STREAM else "main"
one_dim <- function(name, post, prop, S_vals, breaks, R_of = function(S) n_reps(REPS$adv)) {
  list(name = name, S = S_vals, R_of = R_of, n_ov = function(S) overlap_1d(post$ld, prop$ld, S, S, breaks)$n_ov,
       draw = function(S) {
         x <- post$r(S); y <- prop$r(S)
         list(ap = logZ + post$ld(x) - prop$ld(x), aq = logZ + post$ld(y) - prop$ld(y))
       })
}
funnel <- function(d, sv, S_vals) {
  rf <- function(n) { v <- stats::rnorm(n, 0, sv); cbind(v, matrix(stats::rnorm(n * (d - 1)), n, d - 1) * exp(v / 2)) }
  lf <- function(X) stats::dnorm(X[, 1], 0, sv, log = TRUE) + rowSums(stats::dnorm(X[, -1, drop = FALSE], 0, exp(X[, 1] / 2), log = TRUE))
  list(name = sprintf("funnel, d = %d, sigma_v = %g", d, sv), S = S_vals,
       R_of = function(S) if (S >= 4000) n_reps(REPS$adv_funnel_big) else n_reps(REPS$adv_funnel),
       n_ov = function(S) NA_real_,
       draw = function(S) {
         fit <- fit_gauss(rf(S)); X <- rf(S); Y <- rgauss(S, fit)
         list(ap = logZ + lf(X) - lgauss(X, fit), aq = logZ + lf(Y) - lgauss(Y, fit))
       })
}
br <- c(-20, -14, -12, -10, -6, -3, -1, 0, 0.5, 1, 3, 6, 10, 12, 14, 20)
examples <- list(
  one_dim("smooth two-point, mu0 = 0", mix2(1 - eps, -12, 1, 0, 1), mix2(eps, 0, 1, 12, 1),
          c(50, 100, 200, 400, 800, 3200), br),
  one_dim("smooth two-point, mu0 = 0.5", mix2(1 - eps, -12, 1, 0, 1), mix2(eps, 0.5, 1, 12, 1),
          c(50, 100, 200, 400, 800, 3200), br),
  one_dim("rare mode carries overlap", mix2(0.98, 0, 1, 8, 1),
          list(ld = function(x) dn(x, 8, 1.3), r = function(n) stats::rnorm(n, 8, 1.3)),
          c(100, 300, 1000, 3000, 10000), c(-6, -3, 0, 3, 5, 8, 11, 15)),
  funnel(5, 1.5, c(60, 250, 1000, 4000, 16000)),
  funnel(10, 3, c(100, 300, 1000, 4000)))
if (exists("ADV_DEFINE_ONLY")) stop_here <- TRUE else stop_here <- FALSE
if (!stop_here) {
adv_runs <- list()
tab <- do.call(rbind, lapply(quick_pick(examples, c(1L, 3L, 5L)), function(ex) {
  do.call(rbind, lapply(quick_pick(ex$S), function(S) {
    msg("ADV:", ex$name, "S =", S)
    res <- run_reps(ex$R_of(S), cell = paste("adv", STREAM, ex$name, S), fun = function(r) {
      a <- ex$draw(S)
      d <- bridge_diagnostics(a$ap, a$aq, iid = TRUE, eb = FALSE, jackknife = FALSE, rule = "plain")
      Delta <- d$logml - logZ
      data.frame(n_ov = d$n_ov, I = d$I, z = if (is.finite(Delta)) Delta / d$mcse else Delta,
                 cover = is.finite(Delta) && isTRUE(abs(Delta) <= z975 * d$mcse),
                 flag = d$flags[["overlap"]])
    })
    adv_runs[[paste(ex$name, S)]] <<- cbind(example = ex$name, S = S, res)
    data.frame(example = ex$name, S = S, R = nrow(res), n_ov = ex$n_ov(S),
               med_n_ov_hat = stats::median(res$n_ov), med_I = stats::median(res$I),
               coverage = mean(res$cover), ks = ks_normal(res$z))
  }))
}))
saveRDS(do.call(rbind, adv_runs), file.path(results_dir, "runs_adv.rds"))
save_table(tab, "tab_ADV")
print(tab, digits = 3)
lab <- c("smooth two-point, mu0 = 0" = "smooth two-point, $\\mu_0=0$", "smooth two-point, mu0 = 0.5" = "smooth two-point, $\\mu_0=0.5$",
         "rare mode carries overlap" = "rare mode carries overlap",
         "funnel, d = 5, sigma_v = 1.5" = "funnel, $d=5$, $\\sigma_v=1.5$", "funnel, d = 10, sigma_v = 3" = "funnel, $d=10$, $\\sigma_v=3$")
first <- !duplicated(tab$example)
big <- function(x) ifelse(x >= 1000, fint(x), fx(x, 1))
emit_table("tab_ADV", cbind(ifelse(first, lab[tab$example], ""), tab$S, fx(tab$n_ov, 1), big(tab$med_n_ov_hat), fx(tab$med_I, 3),
                            fx(tab$coverage, 3), fx(tab$ks, 3)),
           groups = group_ends(tab$example))
tp0 <- tab[tab$example == "smooth two-point, mu0 = 0", ]
res_put("adv.R", fint(n_reps(REPS$adv))); res_put("adv.R.fun", sprintf("%s to %s", fint(n_reps(REPS$adv_funnel_big)), fint(n_reps(REPS$adv_funnel))))
for (S in c(50, 100, 200)) if (S %in% tp0$S) {
  res_put(sprintf("adv.tp.cov.%d", S), fx(tp0$coverage[tp0$S == S], 2)); res_put(sprintf("adv.tp.nov.%d", S), fx(tp0$n_ov[tp0$S == S], 1))
}
rm_ <- tab[tab$example == "rare mode carries overlap", ]
if (nrow(rm_)) {
  res_put("adv.rare.cov1", fx(rm_$coverage[1], 2)); res_put("adv.rare.nov1", fx(rm_$n_ov[1], 1))
  o <- order(rm_$n_ov); bad <- abs(rm_$coverage[o] - 0.95) > 0.01
  res_put("adv.rare.nominal", if (!any(bad)) fx(min(rm_$n_ov), 0) else if (max(which(bad)) < nrow(rm_)) fx(rm_$n_ov[o][max(which(bad)) + 1], 0) else "--")
}
fu <- tab[grepl("^funnel", tab$example), ]
res_put("adv.fun.small", frange(fu$coverage[fu$med_n_ov_hat < 15], 2)); res_put("adv.fun.big", frange(fu$coverage[fu$med_n_ov_hat >= 15], 2))
res_put("adv.shortfall", fx(max(0, 0.95 - tab$coverage[tab$med_n_ov_hat >= 5]), 2))
res_write("supplement/S05_ADV_adversarial.R")
}
