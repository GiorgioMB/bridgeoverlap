#!/usr/bin/env Rscript
# Proposition C (made model) and Proposition 4 (found model): checks that pass while the
# interval misses. Independent exact draws, d = 2, S1 = S2 = S.
#   bound, fitted     p = N(0, I) restricted to theta_1 > -a (a bound the bridge is not told
#                     about); g = Gaussian fitted by moments on an independent sample of size
#                     S_fit, so g spills mass ~Phi(-a) below the bound, where f = 0;
#   bound, declared   the same with the bound declared (g restricted to theta_1 > -a);
#   missed mode       p = 0.9 N(0, I) + 0.1 N((8, 0), I); the posterior draws come from the major
#                     mode only (a sampler stuck there), g fitted to them.
# Reports coverage of the delta-method interval, the flags of Algorithm 1, the rule
# "MCSE >= 0.3", the plain rule n_ov-hat < 10, and the two finite-sample intervals.
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(dirname(normalizePath(f))) else if (dir.exists("analysis")) "analysis" else ".."})
for (f in c("common.R", "paper.R")) source(file.path(BO_ROOT, "R", f))
suppressPackageStartupMessages(library(bridgeoverlap))
z <- stats::qnorm(0.975); d <- 2L; a <- 3.5
rtrunc1 <- function(n) stats::qnorm(stats::runif(n, stats::pnorm(-a), 1))          # N(0,1) on (-a, Inf)
logp_bound <- function(X) ifelse(X[, 1] > -a, rowSums(stats::dnorm(X, log = TRUE)) - log(1 - stats::pnorm(-a)), -Inf)
gfit <- function(X) list(m = colMeans(X), V = stats::cov(X))
dmvn <- function(X, g) { R <- chol(g$V); Q <- backsolve(R, t(X) - g$m, transpose = TRUE); -colSums(Q^2) / 2 - sum(log(diag(R))) - d / 2 * log(2 * pi) }
rmvn <- function(n, g) sweep(matrix(stats::rnorm(n * d), n) %*% chol(g$V), 2, g$m, `+`)
one <- function(cfg, S) {
  if (cfg$model == "bound") {
    post <- cbind(rtrunc1(S), stats::rnorm(S)); fitX <- cbind(rtrunc1(cfg$Sfit), stats::rnorm(cfg$Sfit)); g <- gfit(fitX)
    prop <- rmvn(if (cfg$declared) 3L * S else S, g)
    if (cfg$declared) { prop <- prop[prop[, 1] > -a, , drop = FALSE][seq_len(S), , drop = FALSE] }
    lg <- function(X) dmvn(X, g) - if (cfg$declared) log(stats::pnorm((g$m[1] + a) / sqrt(g$V[1, 1]))) else 0
    ap <- logp_bound(post) - lg(post); aq <- logp_bound(prop) - lg(prop); logZ <- 0
  } else {
    lp <- function(X) log(0.9 * exp(rowSums(stats::dnorm(X, log = TRUE))) + 0.1 * exp(stats::dnorm(X[, 1], 8, log = TRUE) + stats::dnorm(X[, 2], log = TRUE)))
    post <- matrix(stats::rnorm(S * d), S); g <- gfit(matrix(stats::rnorm(cfg$Sfit * d), cfg$Sfit)); prop <- rmvn(S, g)
    ap <- lp(post) - dmvn(post, g); aq <- lp(prop) - dmvn(prop, g); logZ <- 0
  }
  # f_plain: the overlap flag of Algorithm 1; pass_all: every other check (failure flags, MCSE below 0.3);
  # the interval of Proposition 2 is bridge_fs_interval().
  dg <- bridge_diagnostics(ap, aq, iid = TRUE, eb = FALSE, rule = "plain")
  fs <- bridge_fs_interval(ap, aq, delta = 0.05, ngrid = 150L)
  err <- dg$logml - logZ
  c(err = err, mcse = dg$mcse, n_ov = dg$n_ov, I = dg$I, cover = abs(err) <= z * dg$mcse,
    f_overlap = isTRUE(dg$flags[["overlap"]]), f_fail = any(unlist(dg$flags[setdiff(names(dg$flags), "overlap")]), na.rm = TRUE),
    f_mcse = !(dg$mcse < 0.3), f_plain = !(dg$n_ov >= 10), dead = sum(!is.finite(aq)),
    pass_all = !(any(unlist(dg$flags[setdiff(names(dg$flags), "overlap")]), na.rm = TRUE) || !(dg$mcse < 0.3)),
    cover_fs = fs[["lower"]] <= logZ && logZ <= fs[["upper"]], fs_half = (fs[["upper"]] - fs[["lower"]]) / 2,
    cover_eb = if (!is.null(dg$eb)) dg$eb[["lower"]] <= logZ && logZ <= dg$eb[["upper"]] else NA)
}
cfgs <- list(list(name = "bound, fitted on 4,000", model = "bound", declared = FALSE, Sfit = 4000L),
             list(name = "bound, fitted on 1,000,000", model = "bound", declared = FALSE, Sfit = 1000000L),
             list(name = "bound declared", model = "bound", declared = TRUE, Sfit = 1000000L),
             list(name = "missed mode (sampler stuck)", model = "mode", declared = FALSE, Sfit = 1000000L))
S <- 4000L; Rr <- n_reps(REPS$rare)
tab <- do.call(rbind, lapply(cfgs, function(cf) {
  m <- run_reps(Rr, cell = paste("rare", cf$name), fun = function(r) one(cf, S))
  m <- as.data.frame(m)
  data.frame(config = cf$name, R = nrow(m), coverage = mean(m$cover), med_abs_err_mcse = stats::median(abs(m$err) / m$mcse),
             mean_abs_err = mean(abs(m$err)), med_mcse = stats::median(m$mcse), med_n_ov = stats::median(m$n_ov), med_I = stats::median(m$I),
             f_overlap = 100 * mean(m$f_overlap), f_fail = 100 * mean(m$f_fail), f_mcse = 100 * mean(m$f_mcse), f_plain = 100 * mean(m$f_plain),
             no_dead = 100 * mean(m$dead == 0), pass_all = 100 * mean(m$pass_all), pass_miss = 100 * mean(m$pass_all & !m$cover),
             cov_pass = if (any(m$pass_all == 1)) mean(m$cover[m$pass_all == 1]) else NA, cov_nodead = if (any(m$dead == 0)) mean(m$cover[m$dead == 0]) else NA,
             cov_fs = mean(m$cover_fs), cov_eb = mean(m$cover_eb, na.rm = TRUE),
             fs_half = stats::median(m$fs_half), fs_ratio = stats::median(m$fs_half / (z * m$mcse)))
}))
print(tab, digits = 3)
save_table(tab, "tab_rare")
emit_table("tab_rare", cbind(tab$config, fx(tab$coverage, 3), fx(tab$med_abs_err_mcse, 1), vapply(tab$mean_abs_err, fsci, "", 1), vapply(tab$med_mcse, fsci, "", 1),
                             fx(tab$med_n_ov, 1), fx(tab$med_I, 3), fx(tab$f_fail, 0), fx(tab$f_mcse, 0), fx(tab$f_plain, 0),
                             fx(tab$no_dead, 0), fx(tab$pass_all, 0), fx(tab$pass_miss, 0), fx(tab$cov_fs, 3)),
           spec = "lrrrrrrrrrrrrr",
           header = "Model & coverage & med $|\\Delta|/\\sighat$ & mean $|\\Delta|$ & med MCSE & med $\\neffhat$ & med $\\Ihat$ & failure flags & MCSE${}\\ge0.3$ & $\\neffhat<10$ & no spill draw & pass all & pass and miss & \\cref{prop:ebci}\\\\")
for (k in seq_len(nrow(tab))) {
  key <- c("fit4k", "fit1m", "decl", "mode")[k]
  res_put(paste0("rare.", key, ".cov"), fx(tab$coverage[k], 2)); res_put(paste0("rare.", key, ".ratio"), fx(tab$med_abs_err_mcse[k], 1))
  res_put(paste0("rare.", key, ".flag"), fx(max(tab$f_fail[k], tab$f_mcse[k]), 0)); res_put(paste0("rare.", key, ".fs"), fx(tab$cov_fs[k], 3))
  res_put(paste0("rare.", key, ".plain"), fx(tab$f_plain[k], 0)); res_put(paste0("rare.", key, ".nodead"), fx(tab$no_dead[k], 0))
  res_put(paste0("rare.", key, ".passall"), fx(tab$pass_all[k], 0)); res_put(paste0("rare.", key, ".passmiss"), fx(tab$pass_miss[k], 0))
  res_put(paste0("rare.", key, ".covpass"), fx(tab$cov_pass[k], 2)); res_put(paste0("rare.", key, ".covnodead"), fx(tab$cov_nodead[k], 2))
  res_put(paste0("rare.", key, ".fail"), fx(tab$f_fail[k], 0)); res_put(paste0("rare.", key, ".fshalf"), fsig(tab$fs_half[k], 1))
}
res_put("rare.R", fint(Rr)); res_put("rare.S", fint(S)); res_put("rare.a", a)
res_write("supplement/S21_rare_proposal_mass.R")
