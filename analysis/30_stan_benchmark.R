#!/usr/bin/env Rscript
# Section 7.3: benchmark with Stan's NUTS sampler and bridgesampling on posteriordb
# posteriors whose evidence is known exactly (radon models: quadrature in R/radon.R;
# normal mixture: R/mixture.R; eight schools: script 10). Every run uses 4 chains, 1000
# warm-up iterations and the default split of bridgesampling (first half of every chain
# fits the proposal, method "normal", use_neff = FALSE). Per run the script records our
# diagnostics, the checks of current practice (R-hat, bulk and tail ESS, divergences, the
# MCSE of bridgesampling::error_measures(), warp-III, compute times) and, for the radon
# models, the cross-fitted estimate; for the mixture, the permutation-invariant proposal
# of Fruhwirth-Schnatter (2004).
#
# Configurations can be split across cluster jobs: BO_STAN_CONFIGS="id1,id2" selects
# configurations and BO_REP_RANGE="a-b" a range of replications; every job writes
# results/bench_runs_<id>_<range>.rds, combined by 31_stan_benchmark_summary.R.
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else if (dir.exists("analysis")) "analysis" else "."})
for (f in c("common.R", "stan.R", "paper.R", "radon.R", "mixture.R", "eight_schools.R")) source(file.path(BO_ROOT, "R", f))
stan_ready()

configs <- list(
  list(id = "radon_mn_pooled", file = "radon_pooled_ml", data = "mn", model = "pooled", S = 4000L, rkey = "stan_radon"),
  list(id = "radon_mn_vi_c", file = "radon_vi_centered_ml", data = "mn", model = "vi", S = 4000L, rkey = "stan_radon"),
  list(id = "radon_mn_vi_nc", file = "radon_vi_noncentered_ml", data = "mn", model = "vi", S = 4000L, rkey = "stan_radon"),
  list(id = "radon_mn_hi_c", file = "radon_hi_centered_ml", data = "mn", model = "hi", S = 4000L, rkey = "stan_radon"),
  list(id = "radon_mn_vis_c", file = "radon_vis_centered_ml", data = "mn", model = "vis", S = 4000L, rkey = "stan_radon"),
  list(id = "radon_mn_vis_nc", file = "radon_vis_noncentered_ml", data = "mn", model = "vis", S = 4000L, rkey = "stan_radon"),
  list(id = "radon_all_vi_c", file = "radon_vi_centered_ml", data = "all", model = "vi", S = 4000L, rkey = "stan_radon_big"),
  list(id = "radon_mn_vi_c_x4", file = "radon_vi_centered_ml", data = "mn", model = "vi", S = 16000L, rkey = "stan_radon_x4"),
  list(id = "radon_mn_hi_c_x4", file = "radon_hi_centered_ml", data = "mn", model = "hi", S = 16000L, rkey = "stan_radon_x4"),
  list(id = "radon_mn_vis_c_x4", file = "radon_vis_centered_ml", data = "mn", model = "vis", S = 16000L, rkey = "stan_radon_x4"),
  list(id = "radon_mn_vis_nc_x4", file = "radon_vis_noncentered_ml", data = "mn", model = "vis", S = 16000L, rkey = "stan_radon_x4"),
  list(id = "mixture", file = "normal_mixture_ml", data = "mix", model = "mixture", S = 4000L, rkey = "stan_mix"))
sel <- Sys.getenv("BO_STAN_CONFIGS", "")
if (nzchar(sel)) configs <- Filter(function(cc) cc$id %in% strsplit(sel, ",")[[1]], configs)
if (BO_QUICK) configs <- Filter(function(cc) cc$id %in% c("radon_mn_vi_c", "mixture"), configs)

dat_cache <- list(); truth_cache <- list()
get_data <- function(which) {
  if (is.null(dat_cache[[which]])) dat_cache[[which]] <<- if (which == "mix") mixture_data() else radon_data(which)
  dat_cache[[which]]
}
get_truth <- function(cc) {
  key <- paste(cc$data, cc$model)
  if (is.null(truth_cache[[key]])) {
    truth_cache[[key]] <<- if (cc$model == "mixture") mixture_log_evidence(get_data("mix")$y)$logZ else
      radon_log_evidence(cc$model, get_data(cc$data), n_y = 91L, n_r = 121L, drop = 50)$logZ
    msg("exact log Z for", key, "=", format(truth_cache[[key]], digits = 12))
  }
  truth_cache[[key]]
}
invisible(jackknife_null_quantiles(20L, probs = 1 - 0.05 / 4))
invisible(jackknife_null_quantiles(16L, probs = 1 - 0.05 / 4))
cores <- as.integer(Sys.getenv("BO_STAN_CORES", "1"))

for (cc in configs) {
  dat <- get_data(cc$data); logZ <- get_truth(cc)
  if (cc$data == "mix") dat <- list(N = dat$N, y = dat$y)
  mod <- stan_model(bridgeoverlap_stan_file(cc$file))
  R <- n_reps(REPS[[cc$rkey]], key = paste0("REPS$", cc$rkey))
  reps <- rep_subset(R)
  if (!length(reps)) next
  cell <- paste("bench", cc$id)
  msg(cc$id, ": S =", cc$S, "replications", min(reps), "to", max(reps), "of", R)
  runs <- run_reps(R, cell = cell, cores = cores, reps = reps, fun = function(r) {
    sf <- stan_fit(mod, dat, cc$S, seed = (cell_seed(cell) + r) %% .Machine$integer.max)
    on.exit(stan_cleanup(sf$fit))
    t_b <- system.time(bs <- bridge_sampler_cmdstanr(sf$fit))[["elapsed"]]
    t_w <- system.time(bw <- tryCatch(bridge_sampler_cmdstanr(sf$fit, method = "warp3"), error = function(e) NULL))[["elapsed"]]
    row <- stan_summary_row(bs, logZ, fit = sf$fit, warp3 = if (is.null(bw)) NA_real_ else bw$logml,
                            times = list(sampling = sf$time, bridge = t_b, warp3 = t_w))
    row <- cbind(config = cc$id, model = cc$model, data = cc$data, S_total = cc$S, rep = r,
                 divergent = sf$divergent, max_treedepth = sf$max_treedepth, row)
    if (cc$model == "mixture") {
      mu1 <- posterior::extract_variable_matrix(sf$fit$draws("mu[1]"), "mu[1]")
      mu2 <- posterior::extract_variable_matrix(sf$fit$draws("mu[2]"), "mu[2]")
      row$chains_label1 <- sum(colMeans(mu1 < mu2) > 0.5)                # chains with mu_1 < mu_2
      t_s <- system.time({
        sb <- mixture_symmetric_bridge(sf$fit, dat$y)
        ds <- bridge_diagnostics(sb$a_post, sb$a_prop, chains = sb$chains)
      })[["elapsed"]]
      row$sym_error <- ds$logml - logZ; row$sym_mcse <- ds$mcse; row$sym_n_ov <- ds$n_ov; row$sym_I <- ds$I
      row$sym_pass <- !any(ds$flags %in% TRUE); row$time_symmetric <- t_s
    } else {
      t_c <- system.time(bc <- bridge_sampler_cmdstanr(sf$fit, crossfit = TRUE))[["elapsed"]]
      dc <- suppressMessages(bridge_diagnostics(bc, G = 16L))
      row$cf_error <- dc$logml - logZ; row$cf_mcse <- dc$mcse; row$cf_n_ov <- dc$n_ov; row$cf_n_ovP <- dc$n_ovP
      row$cf_pass <- !any(dc$flags %in% TRUE) && !isTRUE(dc$stuck_chains > 0); row$time_crossfit <- t_c
    }
    row
  })
  tag <- if (nzchar(Sys.getenv("BO_REP_RANGE"))) gsub("[^0-9]+", "_", Sys.getenv("BO_REP_RANGE")) else "all"
  saveRDS(runs, file.path(results_dir, sprintf("bench_runs_%s_%s.rds", cc$id, tag)))
  msg("saved", nrow(runs), "runs of", cc$id)
}
