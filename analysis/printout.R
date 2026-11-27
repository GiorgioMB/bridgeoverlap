#!/usr/bin/env Rscript
# Worked printout of Algorithm 1: replication 1 of the centred eight schools cell with independent
# draws and S = 8000, regenerated from its seed and checked against the run stored by
# 05_eight_schools.R (results/runs_eight_schools_exact.rds). Writes the macros pr.* that
# annotate the flowchart of the decision rule (supplement).
#   Rscript analysis/printout.R
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else if (dir.exists("analysis")) "analysis" else "."})
source(file.path(BO_ROOT, "R", "common.R"))
source(file.path(BO_ROOT, "R", "eight_schools.R"))
source(file.path(BO_ROOT, "R", "paper.R"))
suppressPackageStartupMessages(library(bridgesampling))
stored <- readRDS(file.path(results_dir, "runs_eight_schools_exact.rds"))
exact_runs <- list("centred 8000" = if (any(stored$cell == "centred 8000")) stored[stored$cell == "centred 8000", ] else NULL)
if (!is.null(exact_runs[["centred 8000"]])) {
  set.seed(cell_seed("es exact centred 8000") + 1L)
  X <- es_exact_draws(8000L); b <- es_bounds(colnames(X))
  bs1 <- bridge_sampler(X, log_posterior = es_lp_centered, data = NULL, lb = b$lb, ub = b$ub,
                        use_neff = FALSE, silent = TRUE)
  d1 <- bridge_diagnostics(bs1, iid = TRUE)
  stopifnot(abs((d1$logml - es$logZ) - exact_runs[["centred 8000"]]$err[1]) < 1e-10)
  rep1 <- bridge_report(d1)
}
## macros for the flowchart of the decision rule (values and rules as printed by bridge_report())
if (exists("rep1")) {
  v <- setNames(rep1$value, rep1$quantity); ru <- setNames(rep1$rule, rep1$quantity)
  last_num <- function(x) { m <- regmatches(x, gregexpr("-?[0-9]+([.][0-9]+)?", x))[[1]]; if (length(m)) m[length(m)] else "--" }
  res_put("pr.logml", v[["log Z-hat"]]); res_put("pr.mcse", v[["MCSE"]])
  res_put("pr.nov", v[["overlap check: n_ov-hat"]]); res_put("pr.thr.nov", last_num(ru[["overlap check: n_ov-hat"]]))
  res_put("pr.Vjack", v[["variance check: V_jack"]]); res_put("pr.thr.Vjack", last_num(ru[["variance check: V_jack"]]))
  res_put("pr.int", v[["finite-sample interval"]])
  res_put("pr.half", fx(stats::qnorm(0.975) * d1$mcse, 4)); res_put("pr.S", fint(8000)); res_put("pr.G", fint(d1$settings$G))
  print(rep1)
}
res_write("printout.R")
