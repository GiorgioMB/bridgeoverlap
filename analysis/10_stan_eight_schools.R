#!/usr/bin/env Rscript
# cmdstanr study 1: eight schools with Stan's NUTS sampler and bridgesampling,
# both parameterizations, exact log Z = -31.3113. This replaces the exact draws
# and the hand-written random-walk Metropolis sampler of Section 7.5 by the
# workflow practitioners use (Stan + bridgesampling): does the certificate flag
# the centred funnel, and is the MCSE calibrated for the non-centred model?
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else if (dir.exists("analysis")) "analysis" else "."})
source(file.path(BO_ROOT, "R", "common.R"))
source(file.path(BO_ROOT, "R", "stan.R"))
source(file.path(BO_ROOT, "R", "paper.R"))
stan_ready()

logZ <- eight_schools_log_evidence()
dat <- eight_schools_data()
files <- c("centred" = "eight_schools_centered_ml", "non-centred" = "eight_schools_noncentered_ml")
cfg <- data.frame(model = c("centred", "centred", "centred", "centred", "non-centred", "non-centred"),
                  adapt_delta = c(0.8, 0.8, 0.99, 0.99, 0.8, 0.8),
                  S = c(2000L, 8000L, 2000L, 8000L, 2000L, 8000L))
if (BO_QUICK) cfg <- cfg[c(1, 5), ]
R <- n_reps(REPS$stan_es)
invisible(jackknife_null_quantiles(20L, probs = 1 - 0.05 / 4))

runs <- do.call(rbind, lapply(seq_len(nrow(cfg)), function(i) {
  cc <- cfg[i, ]
  mod <- stan_model(bridgeoverlap_stan_file(files[[cc$model]]))
  msg(cc$model, "adapt_delta =", cc$adapt_delta, "S =", cc$S, "R =", R)
  cell <- paste("stan es", cc$model, cc$adapt_delta, cc$S)
  run_reps(R, cell = cell, cores = as.integer(Sys.getenv("BO_STAN_CORES", "1")), reps = rep_subset(R), fun = function(r) {
    sf <- stan_fit(mod, dat, cc$S, seed = (cell_seed(cell) + r) %% .Machine$integer.max, adapt_delta = cc$adapt_delta)
    on.exit(stan_cleanup(sf$fit))
    t_b <- system.time(bs <- bridge_sampler_cmdstanr(sf$fit))[["elapsed"]]
    t_w <- system.time(bw <- tryCatch(bridge_sampler_cmdstanr(sf$fit, method = "warp3"), error = function(e) NULL))[["elapsed"]]
    row <- stan_summary_row(bs, logZ, fit = sf$fit, warp3 = if (is.null(bw)) NA_real_ else bw$logml,
                            times = list(sampling = sf$time, bridge = t_b, warp3 = t_w))
    cbind(cc, rep = r, divergent = sf$divergent, max_treedepth = sf$max_treedepth, row)
  })
}))
save_table(runs, "stan_eight_schools_runs")
summ <- summarise_stan(runs, c("model", "adapt_delta", "S"))
save_table(summ, "stan_eight_schools_summary")
print(summ[, c("model", "adapt_delta", "S", "R", "coverage", "rmse", "med_mcse", "med_n_ov",
               "med_I", "med_tau", "med_n_ovP", "flag_overlap", "flag_any_failure",
               "coverage_error_measures", "mean_divergent")], digits = 3)

## LaTeX outputs ---------------------------------------------------------------------------------
emit_table("tab_stan_es", cbind(summ$model, format(summ$adapt_delta), summ$S, fx(summ$coverage, 3), fx(summ$rmse, 3),
                                fx(summ$med_mcse, 3), fint(summ$med_n_ov), fx(summ$med_I, 2), fx(summ$med_tau, 1),
                                fx(summ$med_n_ovP, 1), fx(100 * summ$flag_overlap, 0), fx(100 * summ$flag_any_failure, 0),
                                fx(100 * summ$runs_with_stuck_chains, 0), fx(summ$mean_divergent, 1), fx(summ$coverage_error_measures, 3)),
           groups = group_ends(paste(summ$model, summ$adapt_delta)), spec = "lrrcccrcrrrrrrc",
           header = "Parameterization & $\\delta$ & $S$ & coverage & RMSE & med MCSE & med $\\neffhat$ & med $\\Ihat$ & med $\\widehat\\tau$ & med $\\neffPhat$ & overlap flag & failure flags & stuck chains & divergences & \\texttt{error\\_measures}\\\\")
key <- function(x) sprintf("stan.es.%s.%s.%d", if (x$model == "centred") "c" else "nc", sub("[.]", "", format(x$adapt_delta, nsmall = 2)), x$S)
for (i in seq_len(nrow(summ))) {
  x <- summ[i, ]; k <- key(x)
  res_put(paste0(k, ".cov"), fx(x$coverage, 2)); res_put(paste0(k, ".flag"), fx(100 * x$flag_overlap, 0))
  res_put(paste0(k, ".fail"), fx(100 * x$flag_any_failure, 0)); res_put(paste0(k, ".div"), fx(x$mean_divergent, 0))
  res_put(paste0(k, ".I"), fx(x$med_I, 2)); res_put(paste0(k, ".nov"), fint(x$med_n_ov)); res_put(paste0(k, ".em"), fx(x$coverage_error_measures, 2))
  res_put(paste0(k, ".ratio"), fx(x$rmse_over_mcse, 1))
}
c8 <- summ[summ$model == "centred" & summ$adapt_delta == 0.8, ]
nc <- summ[summ$model == "non-centred", ]
res_put("stan.es.R", fint(R)); res_put("stan.es.c.cov", frange(c8$coverage, 2)); res_put("stan.es.c.flag", frange(100 * c8$flag_overlap, 0))
res_put("stan.es.c.div", frange(c8$mean_divergent, 0)); res_put("stan.es.c.em", frange(c8$coverage_error_measures, 2))
res_put("stan.es.nc.cov", frange(nc$coverage, 2)); res_put("stan.es.nc.I", frange(nc$med_I, 2)); res_put("stan.es.nc.flag", frange(100 * nc$flag_overlap, 0))
res_put("stan.es.logZ", fx(logZ, 4))
res_write("10_stan_eight_schools.R")
