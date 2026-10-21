#!/usr/bin/env Rscript
# Section 7.5 / Table ES: eight schools with the default workflow of
# bridgesampling (first half of the draws fits the proposal, S/2 proposal draws).
# Top: exact posterior draws, centred and non-centred. Bottom: random-walk
# Metropolis on the centred parameterization, four chains of length S/4 started
# at stationarity, increments with covariance 2.38^2 Sigma / 10, three
# estimators of the autocorrelation time.
# Expected (top, S = 8000): centred coverage 0.920, RMSE 0.034, median MCSE 0.030,
# median n_ov 456, I 0.35; non-centred RMSE 0.007, n_ov 159, I 0.91.
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else if (dir.exists("analysis")) "analysis" else "."})
source(file.path(BO_ROOT, "R", "common.R"))
source(file.path(BO_ROOT, "R", "eight_schools.R"))
source(file.path(BO_ROOT, "R", "paper.R"))
suppressPackageStartupMessages(library(bridgesampling))

qa <- jackknife_null_quantiles(20L, probs = 1 - 0.05 / 4)
z <- stats::qnorm(0.975)

## Top: exact draws -----------------------------------------------------------------
cfg <- data.frame(param = c("centred", "centred", "centred", "non-centred", "non-centred"),
                  S = c(2000L, 8000L, 32000L, 2000L, 8000L))
cfg$R <- n_reps(REPS$es_exact)
if (BO_QUICK) cfg <- cfg[c(2, 5), ]
exact_runs <- list()
top <- do.call(rbind, lapply(seq_len(nrow(cfg)), function(i) {
  cc <- cfg[i, ]
  msg("exact draws:", cc$param, "S =", cc$S)
  res <- run_reps(cc$R, cell = paste("es exact", cc$param, cc$S), fun = function(r) {
    X <- es_exact_draws(cc$S)
    if (cc$param == "non-centred") X <- es_to_noncentered(X)
    b <- es_bounds(colnames(X))
    lp <- if (cc$param == "centred") es_lp_centered else es_lp_noncentered
    bs <- bridge_sampler(X, log_posterior = lp, data = NULL, lb = b$lb, ub = b$ub,
                         use_neff = FALSE, silent = TRUE)
    d <- bridge_diagnostics(bs, iid = TRUE, rule = "plain", eb_method = "mp")
    data.frame(err = d$logml - es$logZ, mcse = d$mcse, n_ov = d$n_ov, I = d$I, NH0 = d$NH0, T = d$T,
               t_T = d$t_T, V_jack = d$V_jack, I_max = d$I_max, G_max = d$G_max,
               eb_lo = d$eb[["lower"]] - es$logZ, eb_hi = d$eb[["upper"]] - es$logZ,
               f_overlap = d$flags[["overlap"]], f_any = any(d$flags[-1L] %in% TRUE),
               n_flags = sum(d$flags %in% TRUE),
               eb_cover = d$eb[["lower"]] <= es$logZ && es$logZ <= d$eb[["upper"]],
               eb_width = (d$eb[["upper"]] - d$eb[["lower"]]) / (2 * z * d$mcse),
               flag_gate = if (d$I <= 0.5) d$n_ov < 10 else d$NH_lower < 16 * log(40))
  })
  exact_runs[[paste(cc$param, cc$S)]] <<- res
  data.frame(parameterization = cc$param, draws = "exact", S = cc$S, R = nrow(res),
             coverage = mean(abs(res$err) <= z * res$mcse), rmse = sqrt(mean(res$err^2)),
             med_mcse = stats::median(res$mcse), med_n_ov = stats::median(res$n_ov),
             med_I = stats::median(res$I), overlap_flag_pct = 100 * mean(res$f_overlap),
             failure_flags_pct = 100 * mean(res$f_any))
}))
save_table(top, "tab_eight_schools_exact")
print(top, digits = 3)
saveRDS(do.call(rbind, lapply(names(exact_runs), function(k) cbind(cell = k, exact_runs[[k]]))),
        file.path(results_dir, "runs_eight_schools_exact.rds"))

## Worked printout (Figure 1 and tab_printout): replication 1 of the centred cell with S = 8000,
## regenerated from its seed and checked against the stored run, in printout.R.
st <- system2(file.path(R.home("bin"), "Rscript"), shQuote(file.path(BO_ROOT, "printout.R")))
if (!identical(as.integer(st), 0L)) stop("printout.R failed (exit status ", st, ")")

## Middle: diagnosis of the undercoverage of the centred runs ------------------------------
# Arms (centred parameterization, S posterior draws, the bridgesampling estimator on
# vectorized log densities):
#   split      fit on the first half, bridge the second half against S/2 proposal draws
#              (the bridgesampling default), log tau drawn on the inverse-CDF grid;
#   split*     the same with log tau drawn by exact rejection sampling;
#   fixed      the Gaussian with the exact posterior mean and covariance, S/2 + S/2 draws;
#   cross-fit  both halves bridged with the proposal fitted on the other (Algorithm 2).
set.seed(BO_SEED + 1L)
fixed_fit <- es_gfit(es_exact_draws(if (BO_QUICK) 1e5 else 1e6, "rejection"))
es_arm <- function(X, arm) {
  S <- nrow(X); h <- S %/% 2L; A <- X[seq_len(h), , drop = FALSE]; B <- X[h + seq_len(h), , drop = FALSE]
  lr_ <- function(P, fit) es_lp_centered_rows(P) - es_glog(P, fit)
  folds <- NULL
  if (arm %in% c("split", "split*")) {
    fit <- es_gfit(A); Y <- es_gdraw(h, fit); ap <- lr_(B, fit); aq <- lr_(Y, fit)
  } else if (arm == "fixed") {
    Y <- es_gdraw(h, fixed_fit); ap <- lr_(B, fixed_fit); aq <- lr_(Y, fixed_fit)
  } else {
    fA <- es_gfit(A); fB <- es_gfit(B)
    ap <- c(lr_(B, fA), lr_(A, fB)); aq <- c(lr_(es_gdraw(h, fA), fA), lr_(es_gdraw(h, fB), fB))
    folds <- list(post = rep(1:2, each = h), prop = rep(1:2, each = h))
  }
  d <- bridge_diagnostics(ap, aq, iid = TRUE, folds = folds, eb = FALSE, jackknife = FALSE, rule = "plain", eb_method = "mp")
  data.frame(arm = arm, err = d$logml - es$logZ, mcse = d$mcse, n_ov = d$n_ov, I = d$I)
}
diag_S <- quick_pick(c(2000L, 8000L, 32000L), c(1L, 2L))
diag_runs <- do.call(rbind, lapply(diag_S, function(S) {
  msg("diagnosis: S =", S)
  out <- lapply(c("split", "split*", "fixed", "cross-fit"), function(arm) {
    res <- run_reps(n_reps(REPS$es_diag), cell = paste("es diag", arm, S), fun = function(r) {
      es_arm(es_exact_draws(S, if (arm == "split") "grid" else "rejection"), arm)
    })
    cbind(S = S, res)
  })
  do.call(rbind, out)
}))
diag_runs$z <- diag_runs$err / diag_runs$mcse
saveRDS(diag_runs, file.path(results_dir, "runs_eight_schools_diag.rds"))
skew <- function(x) mean((x - mean(x))^3) / stats::sd(x)^3
diag <- do.call(rbind, lapply(split(diag_runs, list(diag_runs$arm, diag_runs$S), drop = TRUE), function(x) {
  data.frame(arm = x$arm[1], S = x$S[1], R = nrow(x), coverage = mean(abs(x$z) <= z),
             cov_se = sqrt(mean(abs(x$z) <= z) * (1 - mean(abs(x$z) <= z)) / nrow(x)),
             bias_mcse = mean(x$err) / stats::median(x$mcse), bias_se = stats::sd(x$err) / sqrt(nrow(x)) / stats::median(x$mcse),
             sd_ratio = stats::sd(x$err) / stats::median(x$mcse), rms_ratio = sqrt(mean(x$err^2) / mean(x$mcse^2)),
             mean_z = mean(x$z), sd_z = stats::sd(x$z), skew_err = skew(x$err),
             med_n_ov = stats::median(x$n_ov), med_I = stats::median(x$I))
}))
diag <- diag[order(match(diag$arm, c("split", "split*", "fixed", "cross-fit")), diag$S), ]
save_table(diag, "tab_eight_schools_diag")
arm_lab <- c("split" = "default split, grid", "split*" = "default split, rejection", "fixed" = "fixed at exact moments", "cross-fit" = "cross-fitted")
emit_table("tab_ES_diag", cbind(ifelse(!duplicated(diag$arm), arm_lab[diag$arm], ""), fint(diag$S), fint(diag$R),
                                fvse(diag$coverage, diag$cov_se, 3), fsvse(diag$bias_mcse, diag$bias_se, 2), fx(diag$sd_ratio, 2),
                                fx(diag$rms_ratio, 2), fx(diag$skew_err, 2), fint(diag$med_n_ov), fx(diag$med_I, 2)),
           groups = group_ends(diag$arm), spec = "lrrccccccc",
           header = "Proposal & $S$ & runs & coverage & bias/MCSE & sd/MCSE & rms ratio & skewness & med $\\neffhat$ & med $\\Ihat$\\\\")
print(diag, digits = 3)
if (have_ggplot()) {
  library(ggplot2)
  qq <- do.call(rbind, lapply(split(diag_runs, list(diag_runs$arm, diag_runs$S), drop = TRUE), function(x) {
    zz <- sort(x$z); n <- length(zz)
    data.frame(arm = x$arm[1], S = paste("S =", x$S[1]), theo = stats::qnorm((seq_len(n) - 0.5) / n), z = zz)
  }))
  qq$S <- factor(qq$S, levels = paste("S =", diag_S))
  pq <- ggplot(qq, aes(theo, z, colour = arm)) + geom_abline(linetype = 2, colour = "grey50") +
    geom_line() + facet_wrap(~S, nrow = 1) + coord_cartesian(xlim = c(-3.5, 3.5), ylim = c(-5, 5)) +
    theme_bw(base_size = 10) + labs(x = "standard normal quantile", y = "quantile of (log Z-hat - log Z) / MCSE", colour = NULL)
  ggsave(fig_file("fig_es_qq"), pq, width = 10, height = 3.4)
  msg("wrote generated/fig_es_qq.pdf")
}

## Bottom: random-walk Metropolis ------------------------------------------------------
set.seed(BO_SEED)
Sigma <- stats::cov(es_exact_draws(200000))
bottom_rows <- list()
rwm_runs <- list()
for (S in quick_pick(c(8000L, 32000L), 1L)) {
  msg("random-walk Metropolis: S =", S)
  res <- run_reps(n_reps(REPS$es_rwm), cell = paste("es rwm", S), fun = function(r) {
    rw <- es_rwm(S %/% 4L, 4L, Sigma)
    arr <- array(unlist(rw$chains), c(S %/% 4L, 10L, 4L))
    arr <- aperm(arr, c(1L, 3L, 2L))   # iterations x chains x variables
    rhat <- max(apply(arr, 3, posterior::rhat)); ess <- min(apply(arr, 3, posterior::ess_bulk))
    b <- es_bounds(es_names)
    ml <- coda::mcmc.list(lapply(rw$chains, coda::mcmc))
    bs <- bridge_sampler(ml, log_posterior = es_lp_centered, data = NULL, lb = b$lb, ub = b$ub,
                         use_neff = FALSE, silent = TRUE)
    chains <- bridge_chain_ids(S %/% 4L, 4L, "per_chain")
    out <- lapply(c(within = "within", concatenated = "concat_geyer", "AR, concatenated" = "concat_ar"),
                  function(tm) {
      d <- bridge_diagnostics(bs, chains = chains, tau_method = tm, jackknife = tm == "within", rule = "plain", eb_method = "mp")
      data.frame(err = d$logml - es$logZ, tau = d$tau, mcse = d$mcse, n_ovP = d$n_ovP,
                 flag = d$flags[["overlap"]], f_cal = d$flags[["calibration"]],
                 f_gmax = d$flags[["influential_block_grubbs"]], stuck = isTRUE(d$stuck_chains > 0))
    })
    cbind(estimator = names(out), do.call(rbind, out), acceptance = rw$acceptance, rhat = rhat, ess = ess)
  })
  rwm_runs[[as.character(S)]] <- res
  for (est in unique(res$estimator)) {
    x <- res[res$estimator == est, ]
    bottom_rows[[length(bottom_rows) + 1L]] <- data.frame(
      S = S, estimator = est, R = nrow(x), med_tau = stats::median(x$tau),
      med_mcse = stats::median(x$mcse), med_n_ovP = stats::median(x$n_ovP),
      flag_pct = 100 * mean(x$flag), coverage = mean(abs(x$err) <= z * x$mcse),
      pct_error_gt_1 = 100 * mean(abs(x$err) > 1), mean_acceptance = mean(x$acceptance),
      calibration_flag_pct = if (est == "within") 100 * mean(x$f_cal, na.rm = TRUE) else NA,
      G_max_flag_pct = if (est == "within") 100 * mean(x$f_gmax, na.rm = TRUE) else NA,
      std_pass_pct = 100 * mean(x$rhat < 1.01 & x$ess > 400))
  }
}
# Beyond R-hat and ESS: coverage by standard checks (max R-hat < 1.01 and min bulk-ESS > 400)
# and by the certificate (within-chain tau), pooled over S.
w_all <- do.call(rbind, lapply(names(rwm_runs), function(S) cbind(S = as.integer(S), rwm_runs[[S]][rwm_runs[[S]]$estimator == "within", ])))
std_ok <- w_all$rhat < 1.01 & w_all$ess > 400
ours_ok <- !w_all$flag & !(w_all$f_cal %in% TRUE) & !(w_all$f_gmax %in% TRUE) & !w_all$stuck
cov_ <- function(i) if (any(i)) sprintf("%s (%d)", fx(mean(abs(w_all$err[i]) <= z * w_all$mcse[i]), 2), sum(i)) else "-- (0)"
rwm2x2 <- data.frame(standard = c("pass", "pass", "fail", "fail"), ours = c("pass", "flag", "pass", "flag"),
                     coverage = c(cov_(std_ok & ours_ok), cov_(std_ok & !ours_ok), cov_(!std_ok & ours_ok), cov_(!std_ok & !ours_ok)))
save_table(rwm2x2, "tab_eight_schools_rwm_2x2")
emit_table("tab_ES_rwm2x2", cbind(rwm2x2$standard, rwm2x2$ours, rwm2x2$coverage), spec = "llr",
           header = "$\\widehat R$ and ESS & certificate and failure flags & coverage (runs)\\\\")
bottom <- do.call(rbind, bottom_rows)
save_table(bottom, "tab_eight_schools_rwm")
print(bottom, digits = 3)

## LaTeX outputs ---------------------------------------------------------------------------------
emit_table("tab_ES_top", cbind(top$parameterization, top$draws, top$S, fx(top$coverage, 3), fx(top$rmse, 3),
                               fx(top$med_mcse, 3), fint(top$med_n_ov), fx(top$med_I, 2), fx(top$overlap_flag_pct, 0),
                               fx(top$failure_flags_pct, 1)))
est_lab <- c(within = "within chains", concatenated = "concatenated", "AR, concatenated" = "AR, concatenated")
emit_table("tab_ES_bottom", cbind(bottom$S, est_lab[bottom$estimator], fx(bottom$med_tau, 0), fx(bottom$med_mcse, 3),
                                  fx(bottom$med_n_ovP, 1), fx(bottom$flag_pct, 0), fx(bottom$coverage, 3)))
# exact log Z by a second quadrature route (trapezoid rule on a fine grid in log tau)
lt <- seq(log(1e-8), log(1e4), length.out = 400001)
ld <- es_log_marg_tau(exp(lt)) + stats::dcauchy(exp(lt), 0, 5, log = TRUE) + log(2) + lt
logZ2 <- max(ld) + log(sum((exp(ld[-1] - max(ld)) + exp(ld[-length(ld)] - max(ld))) / 2 * diff(lt)))
res_put("es.logZ", fx(es$logZ, 4)); res_put("es.quad.diff", fsci(abs(logZ2 - es$logZ), 0))
res_put("es.R.exact", fint(n_reps(REPS$es_exact))); res_put("es.R.rwm", fint(n_reps(REPS$es_rwm))); res_put("es.R.diag", fint(n_reps(REPS$es_diag)))
tc <- top[top$parameterization == "centred", ]; tn <- top[top$parameterization == "non-centred", ]
res_put("es.c.I", fx(stats::median(tc$med_I), 2)); res_put("es.nc.I", frange(tn$med_I, 2))
common <- intersect(tc$S, tn$S)
res_put("es.rmse.ratio", frange(tc$rmse[match(common, tc$S)] / tn$rmse[match(common, tn$S)], 1))
res_put("es.cov", frange(top$coverage, 2))
res_put("es.overlap.phrase", if (all(top$overlap_flag_pct == 0)) "no overlap flag is raised" else
  sprintf("the overlap flag is raised in %s\\%% of runs", frange(top$overlap_flag_pct, 1)))
res_put("es.fail", frange(top$failure_flags_pct, 1))
res_put("es.cov.se", fx(sqrt(0.95 * 0.05 / n_reps(REPS$es_exact)), 3))
tcc <- top[top$parameterization == "centred", ]
pooled <- sum(tcc$coverage * tcc$R) / sum(tcc$R)
res_put("es.c.cov.pooled", fx(pooled, 3)); res_put("es.c.R.pooled", fint(sum(tcc$R)))
res_put("es.c.cov.z", fx((pooled - 0.95) / sqrt(0.95 * 0.05 / sum(tcc$R)), 1))
# diagnosis: per arm, coverage range over S and pooled bias in MCSE units
for (a in unique(diag$arm)) {
  x <- diag[diag$arm == a, ]; k <- c("split" = "split", "split*" = "splitrej", "fixed" = "fixed", "cross-fit" = "cross")[[a]]
  res_put(sprintf("es.diag.%s.cov", k), frange(x$coverage, 3)); res_put(sprintf("es.diag.%s.bias", k), frange(x$bias_mcse, 2))
  res_put(sprintf("es.diag.%s.sd", k), frange(x$sd_ratio, 2)); res_put(sprintf("es.diag.%s.skew", k), frange(x$skew_err, 2))
}
res_put("es.diag.covse", fx(max(diag$cov_se), 3))
std <- w_all$rhat < 1.01 & w_all$ess > 400
res_put("es.rwm.std.pass", fx(100 * mean(std), 0)); res_put("es.rwm.std.cov", if (any(std)) fx(mean(abs(w_all$err[std]) <= z * w_all$mcse[std]), 2) else "--")
res_put("es.rwm.std.flagged", if (any(std)) fx(100 * mean(!ours_ok[std]), 0) else "--")
for (S in tc$S) res_put(sprintf("es.c.ratio.%d", S), fx(tc$rmse[tc$S == S] / tc$med_mcse[tc$S == S], 2))
fr <- exact_runs[["centred 8000"]]
if (!is.null(fr)) {
  f1 <- fr[1, ]
  res_put("es.first.err", fx(f1$err, 3)); res_put("es.first.mcse", fx(f1$mcse, 3)); res_put("es.first.nov", fint(f1$n_ov))
  res_put("es.first.I", fx(f1$I, 2)); res_put("es.first.NH0", fint(f1$NH0)); res_put("es.first.T", fx(f1$T, 4))
  res_put("es.first.tT", fx(f1$t_T, 2)); res_put("es.first.Vjack", fx(f1$V_jack, 2)); res_put("es.first.Imax", fx(f1$I_max, 2))
  res_put("es.first.Gmax", fx(f1$G_max, 2)); res_put("es.first.eb", sprintf("[%s,%s]", fx(f1$eb_lo, 3), fx(f1$eb_hi, 3)))
  res_put("es.first.flags", if (f1$n_flags == 0) "no flag is raised" else sprintf("%d flag%s %s raised", f1$n_flags,
                                                                                     if (f1$n_flags > 1) "s" else "", if (f1$n_flags > 1) "are" else "is"))
}
res_put("es.rwm.acc", fx(mean(unlist(lapply(rwm_runs, function(x) x$acceptance[x$estimator == "within"]))), 2))
res_put("es.rwm.tau", frange(bottom$med_tau, 0))
res_put("es.rwm.cov", frange(bottom$coverage, 2))
res_put("es.rwm.errmax", fx(max(unlist(lapply(rwm_runs, function(x) abs(x$err[x$estimator == "within"])))), 0))
for (S in names(rwm_runs)) {
  x <- rwm_runs[[S]]; w <- x[x$estimator == "within", ]
  res_put(sprintf("es.rwm.n1.%s", S), fint(sum(abs(w$err) > 1)))
  res_put(sprintf("es.rwm.stuck.%s", S), fx(100 * mean(w$stuck), 0))
  for (est in names(est_lab)) {
    bb <- bottom[bottom$S == as.numeric(S) & bottom$estimator == est, ]
    key <- c(within = "within", concatenated = "concat", "AR, concatenated" = "ar")[[est]]
    res_put(sprintf("es.rwm.tau.%s.%s", key, S), fx(bb$med_tau, 0))
    res_put(sprintf("es.rwm.flag.%s.%s", key, S), fx(bb$flag_pct, 0))
    res_put(sprintf("es.rwm.cov.%s.%s", key, S), fx(100 * bb$coverage, 0))
  }
  bw <- bottom[bottom$S == as.numeric(S) & bottom$estimator == "within", ]
  res_put(sprintf("es.rwm.cal.%s", S), fx(bw$calibration_flag_pct, 0))
  res_put(sprintf("es.rwm.gmax.%s", S), fx(bw$G_max_flag_pct, 0))
}
res_write("05_eight_schools.R")
