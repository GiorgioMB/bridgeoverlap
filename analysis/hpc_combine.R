#!/usr/bin/env Rscript
# Cluster runs: rerun the Stan studies 10-12 in combine mode, which rebuilds their tables
# from the replications cached by the array tasks (see run_reps() in R/common.R).
root <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else "analysis"})
for (s in c("10_stan_eight_schools.R", "11_stan_regression_dimension.R", "12_stan_mixture.R")) {
  cat("== combining", s, "\n")
  code <- system2(file.path(R.home("bin"), "Rscript"), shQuote(file.path(root, s)), env = c("BO_COMBINE=true", "BO_REP_RANGE="))
  if (code != 0) stop(s, " failed in combine mode")
}
