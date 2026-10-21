#!/usr/bin/env Rscript
# cmdstanr study 3: a real sampler failure with a known answer. Target: the
# two-component mixture of Section 7.4 (wide component at -m, narrow component
# at +m, equal weights), log Z known. Three configurations:
#   m = 3, random inits   the target of Section 7.4 with Stan's default initial
#                         values; NUTS crosses between the modes only
#                         occasionally in 1000 iterations, so the run is
#                         metastable (Section 7.4(b));
#   m = 5, inits 3:1      three chains start in the wide mode and one in the
#                         narrow mode; NUTS cannot cross, so 75% of the draws
#                         are in one mode (mode imbalance, Section 7.4(a));
#   m = 5, inits 4:0      all chains start in the wide mode; the proposal is
#                         fitted to one mode only. This is the blind spot of
#                         Proposition 5: log Z-hat is off by log(1/2) with a tiny
#                         MCSE, and no diagnostic based on the bridge output can
#                         see the missing mode.
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else if (dir.exists("analysis")) "analysis" else "."})
source(file.path(BO_ROOT, "R", "common.R"))
source(file.path(BO_ROOT, "R", "stan.R"))
source(file.path(BO_ROOT, "R", "paper.R"))
stan_ready()

logZ <- 1.234
S <- 4000L
R <- n_reps(REPS$stan_toy_mix)
mod <- stan_model(bridgeoverlap_stan_file("gaussian_mixture_ml"))
invisible(jackknife_null_quantiles(20L, probs = 1 - 0.05 / 4))
cfgs <- list(
  list(config = "m = 3, random inits", m = 3, x1 = NULL),
  list(config = "m = 5, inits 3:1", m = 5, x1 = c(-5, -5, -5, 5)),
  list(config = "m = 5, inits 4:0 (one mode)", m = 5, x1 = c(-5, -5, -5, -5))
)
runs <- do.call(rbind, lapply(quick_pick(cfgs, c(1L, 2L, 3L)), function(cc) {
  dat <- list(m = cc$m, s_narrow = 0.5, w = 0.5, log_Z = logZ)
  msg(cc$config, "R =", R)
  cell <- paste("stan toy mixture", cc$config)
  run_reps(R, cell = cell, reps = rep_subset(R), cores = as.integer(Sys.getenv("BO_STAN_CORES", "1")), fun = function(r) {
    init <- if (is.null(cc$x1)) NULL else
      lapply(cc$x1, function(x1) list(x = c(x1 + stats::rnorm(1, 0, 0.1), stats::rnorm(1, 0, 0.1))))
    sf <- stan_fit(mod, dat, S, seed = (cell_seed(cell) + r) %% .Machine$integer.max, init = init)
    on.exit(stan_cleanup(sf$fit))
    x1 <- posterior::extract_variable_matrix(sf$fit$draws("x[1]"), "x[1]")
    row <- stan_summary_row(bridge_sampler_cmdstanr(sf$fit), logZ)
    cbind(config = cc$config, rep = r, divergent = sf$divergent,
          frac_narrow = mean(x1 > 0), chains_switching = sum(apply(x1, 2, function(v) any(v > 0) && any(v < 0))),
          row)
  })
}))
save_table(runs, "stan_mixture_runs")
summ <- summarise_stan(runs, "config")
summ$mean_frac_narrow <- tapply(runs$frac_narrow, runs$config, mean)[summ$config]
summ$bias_in_mcse <- tapply(runs$error, runs$config, mean)[summ$config] / summ$med_mcse   # mean bias / median MCSE, as elsewhere
save_table(summ, "stan_mixture_summary")
print(summ[, c("config", "R", "mean_frac_narrow", "bias", "bias_in_mcse", "coverage", "med_n_ov",
               "med_I", "flag_overlap", "flag_calibration", "flag_jackknife_variance",
               "flag_influential", "flag_any_failure")], digits = 3)

## LaTeX outputs ---------------------------------------------------------------------------------
lab <- c("m = 3, random inits" = "$m=3$, random initial values", "m = 5, inits 3:1" = "$m=5$, initial values 3:1",
         "m = 5, inits 4:0 (one mode)" = "$m=5$, initial values 4:0")
summ <- summ[order(match(summ$config, names(lab))), ]
emit_table("tab_stan_mix", cbind(lab[summ$config], fx(summ$mean_frac_narrow, 2), fsign(summ$bias, 3), fsign(summ$bias_in_mcse, 1),
                                 fx(summ$coverage, 3), fint(summ$med_n_ov), fx(summ$med_I, 2), fx(100 * summ$flag_overlap, 0),
                                 fx(100 * summ$flag_calibration, 0), fx(100 * summ$flag_jackknife_variance, 0),
                                 fx(100 * summ$flag_influential, 0), fx(100 * summ$flag_any_failure, 0)),
           spec = "lccccrcrrrrr",
           header = "Configuration & narrow-mode fraction & bias & bias/MCSE & coverage & med $\\neffhat$ & med $\\Ihat$ & overlap flag & $t_\\Tp$ & $\\Vjack$ & $\\Imax$ or $\\Gr$ & any failure flag\\\\")
keys <- c("m = 3, random inits" = "m3", "m = 5, inits 3:1" = "m31", "m = 5, inits 4:0 (one mode)" = "m40")
for (i in seq_len(nrow(summ))) {
  x <- summ[i, ]; k <- paste0("stan.mix.", keys[[x$config]])
  res_put(paste0(k, ".frac"), fx(x$mean_frac_narrow, 2)); res_put(paste0(k, ".bias"), fsign(x$bias, 3))
  res_put(paste0(k, ".absbiasmcse"), fx(abs(x$bias_in_mcse), 0)); res_put(paste0(k, ".cov"), fx(x$coverage, 2))
  res_put(paste0(k, ".any"), fx(100 * x$flag_any_failure, 0)); res_put(paste0(k, ".cal"), fx(100 * x$flag_calibration, 0))
  res_put(paste0(k, ".infl"), fx(100 * x$flag_influential, 0)); res_put(paste0(k, ".vj"), fx(100 * x$flag_jackknife_variance, 0))
  res_put(paste0(k, ".I"), fx(x$med_I, 2)); res_put(paste0(k, ".overlap"), fx(100 * x$flag_overlap, 0))
}
res_put("stan.mix.R", fint(R)); res_put("stan.mix.logZ", fx(logZ, 3))
res_write("12_stan_mixture.R")
