#!/usr/bin/env Rscript
# cmdstanr study 2: overlap against dimension with an exact answer. Conjugate
# normal-inverse-gamma regression of the motorcycle data (MASS::mcycle,
# standardized acceleration) on a B-spline basis with K columns, so the
# unconstrained dimension is K + 1 and log Z is available in closed form.
# For small K the fitted Gaussian proposal is nearly exact (I-hat > 1/2, the
# "near-perfect" regime of Remark 1, where a small n_ov-hat is harmless); as K
# grows the overlap degrades. Each fit is bridged twice: with the default
# split and with two-fold cross-fitting (Remark 2), which uses every draw for
# both roles.
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else if (dir.exists("analysis")) "analysis" else "."})
source(file.path(BO_ROOT, "R", "common.R"))
source(file.path(BO_ROOT, "R", "stan.R"))
source(file.path(BO_ROOT, "R", "paper.R"))
stan_ready()
if (!requireNamespace("MASS", quietly = TRUE)) stop("Package MASS is needed for the mcycle data.")

mc <- MASS::mcycle
y <- as.numeric(scale(mc$accel))
v0 <- 10; a0 <- 2; b0 <- 1
S <- 4000L
R <- n_reps(REPS$stan_reg)
mod <- stan_model(bridgeoverlap_stan_file("conjugate_regression_ml"))
invisible(jackknife_null_quantiles(20L, probs = 1 - 0.05 / 4))
invisible(jackknife_null_quantiles(16L, probs = 1 - 0.05 / 4))

runs <- do.call(rbind, lapply(quick_pick(c(4L, 8L, 16L, 32L, 64L)), function(K) {
  X <- unclass(splines::bs(mc$times, df = K, intercept = TRUE))
  attributes(X) <- list(dim = dim(X))
  truth <- conjugate_log_evidence(X, y, v0, a0, b0)
  dat <- list(n = length(y), K = K, X = X, y = y, v0 = v0, a0 = a0, b0 = b0)
  msg("K =", K, "(dimension", K + 1, ") log Z =", round(truth, 4))
  cell <- paste("stan reg", K)
  run_reps(R, cell = cell, reps = rep_subset(R), cores = as.integer(Sys.getenv("BO_STAN_CORES", "1")), fun = function(r) {
    sf <- stan_fit(mod, dat, S, seed = (cell_seed(cell) + r) %% .Machine$integer.max)
    on.exit(stan_cleanup(sf$fit))
    a <- stan_summary_row(bridge_sampler_cmdstanr(sf$fit), truth)
    b <- stan_summary_row(bridge_sampler_cmdstanr(sf$fit, crossfit = TRUE), truth, G = 16L)
    rbind(cbind(K = K, dim = K + 1L, arm = "split", rep = r, divergent = sf$divergent, a),
          cbind(K = K, dim = K + 1L, arm = "cross-fit", rep = r, divergent = sf$divergent, b))
  })
}))
save_table(runs, "stan_regression_dimension_runs")
summ <- summarise_stan(runs, c("K", "dim", "arm"))
summ$S <- S
summ <- summ[order(summ$K, summ$arm != "split"), ]   # default split first
save_table(summ, "stan_regression_dimension_summary")
print(summ[, c("dim", "arm", "R", "coverage", "bias", "rmse", "med_mcse", "rmse_over_mcse",
               "med_n_ov", "med_I", "flag_overlap", "flag_any_failure")], digits = 3)

## LaTeX outputs ---------------------------------------------------------------------------------
emit_table("tab_stan_reg", cbind(summ$dim, summ$arm, fx(summ$coverage, 3), fsign(summ$bias / summ$med_mcse, 2), fx(summ$rmse, 4),
                                 fx(summ$med_mcse, 4), fx(summ$rmse_over_mcse, 2), fint(summ$med_n_ov), fx(summ$med_I, 2),
                                 fx(100 * summ$flag_overlap, 0), fx(100 * summ$flag_any_failure, 0)),
           groups = group_ends(summ$dim), spec = "rlcccccrcrr",
           header = "Dimension & bridging & coverage & bias/MCSE & RMSE & med MCSE & RMSE/MCSE & med $\\neffhat$ & med $\\Ihat$ & overlap flag & failure flags\\\\")
sp <- summ[summ$arm == "split", ]; cf <- summ[summ$arm == "cross-fit", ]
lo <- sp[which.min(sp$dim), ]; hi <- sp[which.max(sp$dim), ]
res_put("stan.reg.R", fint(R)); res_put("stan.reg.S", fint(S)); res_put("stan.reg.dims", sprintf("%d--%d", min(summ$dim), max(summ$dim)))
res_put("stan.reg.dmin", lo$dim); res_put("stan.reg.dmax", hi$dim)
res_put("stan.reg.I.dmin", fx(lo$med_I, 2)); res_put("stan.reg.I.dmax", fx(hi$med_I, 2))
res_put("stan.reg.nov.dmin", fint(lo$med_n_ov)); res_put("stan.reg.nov.dmax", fint(hi$med_n_ov))
res_put("stan.reg.cov", frange(summ$coverage, 2)); res_put("stan.reg.flag", frange(100 * summ$flag_overlap, 0))
m <- match(sp$dim, cf$dim)
res_put("stan.reg.cf.gain", frange(100 * (1 - cf$rmse[m] / sp$rmse), 0))
res_write("11_stan_regression_dimension.R")
