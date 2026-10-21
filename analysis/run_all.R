#!/usr/bin/env Rscript
# Run every analysis script in a fresh R process, in order, then assemble the
# supplementary tables. Usage (from the repository root):
#   Rscript analysis/run_all.R                 # everything that needs no Stan
#   Rscript analysis/run_all.R --stan          # full run: also the cmdstanr studies 10-12, 30-31
# (20_posteriordb_suite.R is optional: its references are
# computed with the same estimator; the benchmark of 30 has exact references.)
# The last step, supplement/S99_assemble.R, writes generated/results.tex and
# generated/MANIFEST.txt; with --stan, generated/ then holds every table, figure and
# result macro.
#   BO_QUICK=true BO_SCALE=0.01 Rscript analysis/run_all.R   # smoke test (minutes)
root <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else if (dir.exists("analysis")) "analysis" else "."})
args <- commandArgs(TRUE)
main <- c("01_universality.R", "02_markov_chains.R", "03_proposal_fitting.R", "04_sampler_failure.R",
          "05_eight_schools.R")
supp <- sort(list.files(file.path(root, "supplement"), pattern = "^S[0-9]{2}_.*[.]R$"))
stan <- c("10_stan_eight_schools.R", "11_stan_regression_dimension.R", "12_stan_mixture.R", "30_stan_benchmark.R",
          "31_stan_benchmark_summary.R")
# S98 pools the per-run records of 01, 05 and S05, so it runs after them; S99 assembles.
supp <- setdiff(supp, c("S98_operating_characteristics.R", "S99_assemble.R"))
scripts <- c(file.path(root, main), file.path(root, "supplement", supp),
             file.path(root, "supplement", "S98_operating_characteristics.R"),
             if ("--stan" %in% args) file.path(root, stan), file.path(root, "supplement", "S99_assemble.R"))
status <- vapply(scripts, function(s) {
  cat("\n=====", basename(s), "=====\n")
  t0 <- Sys.time()
  code <- system2(file.path(R.home("bin"), "Rscript"), shQuote(s))
  cat(sprintf("----- %s: exit %d after %.1f min\n", basename(s), code, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  code
}, integer(1))
cat("\nSummary:\n"); print(data.frame(script = basename(scripts), exit = unname(status)), row.names = FALSE)
if (any(status != 0)) quit(status = 1)
