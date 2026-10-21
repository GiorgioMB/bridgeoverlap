#!/usr/bin/env Rscript
# Checks that the cluster environment can run the reproduction (run before hpc/submit.sh,
# which calls it): R packages, CmdStan, and one Stan program compiled with model methods.
need <- c("bridgeoverlap", "bridgesampling", "coda", "posterior", "cmdstanr", "Rcpp", "RcppEigen", "jsonlite",
          "mvtnorm", "MASS", "ggplot2", "patchwork", "Rmpfr", "Brobdingnag", "testthat")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) {
  cat("Missing R packages:", paste(miss, collapse = ", "), "\n",
      "Install from CRAN with install.packages(), cmdstanr with\n",
      "  install.packages('cmdstanr', repos = c('https://stan-dev.r-universe.dev', getOption('repos')))\n",
      "and bridgeoverlap with R CMD INSTALL . from the repository root.\n")
  quit(status = 1)
}
v <- tryCatch(cmdstanr::cmdstan_version(), error = function(e) NA)
if (is.na(v)) { cat("CmdStan not found: run cmdstanr::install_cmdstan() or set CMDSTAN.\n"); quit(status = 1) }
cat("CmdStan", v, "| bridgesampling", as.character(utils::packageVersion("bridgesampling")), "\n")
d <- file.path(tempdir(), "bo_check"); dir.create(d, showWarnings = FALSE)
Sys.setenv(BO_STAN_BUILD_DIR = d)
library(bridgeoverlap)
mod <- cmdstan_model_bridge(bridgeoverlap_stan_file("eight_schools_noncentered_ml"))
fit <- mod$sample(data = eight_schools_data(), seed = 1, chains = 2, parallel_chains = 1, refresh = 0,
                  show_messages = FALSE, iter_warmup = 300, iter_sampling = 300)
d1 <- bridge_diagnostics(bridge_sampler_cmdstanr(fit))
cat(sprintf("eight schools: log Z-hat %.3f (exact %.3f), MCSE %.3f\n", d1$logml, eight_schools_log_evidence(), d1$mcse))
if (abs(d1$logml - eight_schools_log_evidence()) > 0.5) { cat("unexpected estimate\n"); quit(status = 1) }
cat("setup OK\n")
