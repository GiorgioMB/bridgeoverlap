# Helpers for the cmdstanr studies -------------------------------------------------

stan_ready <- function() {
  if (identical(tolower(Sys.getenv("BO_COMBINE")), "true")) return(invisible(NULL))   # summary-only rerun from the caches
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("These scripts need cmdstanr: install.packages('cmdstanr', ",
         "repos = c('https://stan-dev.r-universe.dev', getOption('repos')))", call. = FALSE)
  }
  p <- tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)
  if (is.null(p)) stop("CmdStan not found: run cmdstanr::install_cmdstan().", call. = FALSE)
  msg("CmdStan", as.character(cmdstanr::cmdstan_version()), "at", p)
  invisible(p)
}

.stan_models <- new.env()
stan_model <- function(file) {
  if (identical(tolower(Sys.getenv("BO_COMBINE")), "true")) return(NULL)   # summary-only rerun from the caches
  key <- normalizePath(file)
  if (is.null(.stan_models[[key]])) {
    msg("compiling", basename(file), "with model methods (once per session)")
    .stan_models[[key]] <- cmdstan_model_bridge(file)
  }
  .stan_models[[key]]
}

# Sample with NUTS (sequential chains, no console output).
stan_fit <- function(mod, data, S, seed, chains = 4L, iter_warmup = 1000L,
                     init = NULL, adapt_delta = 0.8) {
  fit <- mod$sample(data = data, chains = chains, iter_warmup = iter_warmup,
                    iter_sampling = S %/% chains, seed = seed, init = init,
                    adapt_delta = adapt_delta, refresh = 0, show_messages = FALSE,
                    show_exceptions = FALSE)
  ds <- fit$diagnostic_summary(quiet = TRUE)
  tm <- tryCatch(fit$time()$total, error = function(e) NA_real_)
  list(fit = fit, divergent = sum(ds$num_divergent), max_treedepth = sum(ds$num_max_treedepth), time = tm)
}

# Delete the CSV files of a fit once everything needed has been extracted.
stan_cleanup <- function(fit) invisible(unlink(fit$output_files(), force = TRUE))

# Diagnostics of a bridge object plus the alternative autocorrelation estimators
# of Table ES (bottom) and bridgesampling's own error measures.
# Standard MCMC checks of current practice: largest rank-normalized split R-hat and
# smallest bulk and tail ESS over all parameters (Vehtari et al. 2021).
stan_mcmc_checks <- function(fit) {
  x <- fit$draws(format = "draws_array")
  x <- posterior::subset_draws(x, variable = setdiff(posterior::variables(x), "lp__"))
  sm <- posterior::summarise_draws(x, rhat = posterior::rhat, ess_bulk = posterior::ess_bulk, ess_tail = posterior::ess_tail)
  data.frame(rhat_max = max(sm$rhat, na.rm = TRUE), ess_bulk_min = min(sm$ess_bulk, na.rm = TRUE),
             ess_tail_min = min(sm$ess_tail, na.rm = TRUE))
}

stan_summary_row <- function(bs, logZ_true, G = 20L, fit = NULL, warp3 = NULL, times = NULL) {
  t_d <- system.time(d <- suppressMessages(bridge_diagnostics(bs, G = G, keep_responsibilities = TRUE)))[["elapsed"]]
  row <- as.data.frame(d)
  pp <- d$responsibilities$posterior
  for (tm in c("concat_geyer", "concat_ar")) {
    tau <- iat(pp, d$settings$chain_lengths, tm)
    n_P <- d$settings$S2 * d$V_prop + d$settings$S1 * d$V_post * tau
    row[[paste0("tau_", tm)]] <- tau
    row[[paste0("mcse_", tm)]] <- sqrt(n_P) / d$NH0
    row[[paste0("n_ovP_", tm)]] <- d$n_ov^2 / n_P
  }
  em <- tryCatch(bridgesampling::error_measures(bs)$cv, error = function(e) NA_real_)
  row$cv_error_measures <- em
  row$mcse_logml_bs <- if (is.null(bs$mcse_logml)) NA_real_ else bs$mcse_logml
  row$logml_bs <- bs$logml
  row$truth <- logZ_true
  row$error <- d$logml - logZ_true
  row$logml_warp3 <- if (is.null(warp3)) NA_real_ else warp3
  row$mcse_logml_bs <- tryCatch(mcse_logml_bs(bs, use_ess = TRUE), error = function(e) NA_real_)
  row$mcse_logml_native <- if (is.null(bs$mcse_logml)) NA_real_ else bs$mcse_logml
  row$time_diagnostics <- t_d
  for (nm in names(times)) row[[paste0("time_", nm)]] <- times[[nm]]
  if (!is.null(fit)) row <- cbind(row, stan_mcmc_checks(fit))
  row
}

summarise_stan <- function(df, by) {
  z <- stats::qnorm(0.975)
  cov_of <- function(err, se) mean(abs(err) <= z * se, na.rm = TRUE)
  gate <- function(x, n) mean(n < 30, na.rm = TRUE)                    # plain rule (Section 6.1)
  out <- lapply(split(df, df[by], drop = TRUE), function(x) {
    data.frame(
      x[1L, by, drop = FALSE], R = nrow(x),
      coverage = cov_of(x$error, x$mcse),
      coverage_concat_geyer = cov_of(x$error, x$mcse_concat_geyer),
      coverage_concat_ar = cov_of(x$error, x$mcse_concat_ar),
      coverage_error_measures = cov_of(x$error, x$cv_error_measures),
      bias = mean(x$error), rmse = sqrt(mean(x$error^2)), med_mcse = stats::median(x$mcse),
      rmse_over_mcse = sqrt(mean(x$error^2)) / stats::median(x$mcse),
      med_n_ov = stats::median(x$n_ov), med_I = stats::median(x$I),
      med_tau = stats::median(x$tau), med_tau_concat_geyer = stats::median(x$tau_concat_geyer),
      med_tau_concat_ar = stats::median(x$tau_concat_ar), med_n_ovP = stats::median(x$n_ovP),
      flag_overlap = mean(!(x$n_ovP >= 30)),
      flag_overlap_concat_geyer = gate(x, x$n_ovP_concat_geyer),
      flag_overlap_concat_ar = gate(x, x$n_ovP_concat_ar),
      flag_any_failure = mean(x$flag_any_failure, na.rm = TRUE),
      flag_calibration = mean(x$flag_calibration, na.rm = TRUE),
      flag_jackknife_variance = mean(x$flag_jackknife_variance, na.rm = TRUE),
      flag_influential = mean(x$flag_influential_imax | x$flag_influential_grubbs, na.rm = TRUE),
      runs_with_stuck_chains = mean(x$stuck_chains > 0, na.rm = TRUE),
      mean_divergent = if ("divergent" %in% names(x)) mean(x$divergent) else NA_real_,
      row.names = NULL)
  })
  do.call(rbind, out)
}

# The MCSE that bridgesampling 1.2-1 (GitHub) reports as `mcse_logml`, recomputed from the
# stored log ratios at the final estimate, so that the "MCSE below 0.3" check of current
# practice can use it whichever version is installed. `use_ess = TRUE` is its default
# (use_ess = use_neff = TRUE): the posterior term is inflated by length / ESS, with the ESS
# of coda::effectiveSize() on the concatenated chains. Both relative variances are divided
# by the number of proposal draws, which is correct here because N1 = N2 in the default
# workflow. Numerically stable forms of numi and deni avoid overflow without Brobdingnag.
mcse_logml_bs <- function(bs, use_ess = TRUE, neff = NULL) {
  l1 <- bs$q11 - bs$q12; l2 <- bs$q21 - bs$q22
  lstar <- stats::median(l1); n1 <- length(l1); n2 <- length(l2)
  if (is.null(neff)) neff <- n1          # the estimator with s1 = S1 / N (use_neff = FALSE)
  s1 <- neff / (neff + n2); s2 <- n2 / (neff + n2)
  r <- exp(bs$logml - lstar)
  numi <- 1 / (s1 + s2 * r * exp(-(l2 - lstar)))
  deni <- 1 / (s1 * exp(l1 - lstar) + s2 * r)
  var_deni <- stats::var(deni)
  if (use_ess) var_deni <- var_deni * length(deni) / mean(coda::effectiveSize(deni))
  var_r <- mean(numi)^2 / mean(deni)^2 * (stats::var(numi) / mean(numi)^2 + var_deni / mean(deni)^2) / length(numi)
  sqrt(log1p(var_r / r^2))
}
