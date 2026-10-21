#!/usr/bin/env Rscript
# Worked printout of Algorithm 1: replication 1 of the centred eight schools cell with exact
# draws and S = 8000, regenerated from its seed and checked against the run stored by
# 05_eight_schools.R (results/runs_eight_schools_exact.rds). Writes tab_printout.tex and the
# macros pr.* that annotate the decision-rule flowchart (Figure 1).
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
  qlab <- c("log Z-hat" = "$\\log\\Zhat$", "MCSE" = "MCSE $\\sighat$", "n_ov-hat" = "$\\neffhat$", "I-hat" = "$\\Ihat$",
            "N H0-hat" = "$N\\Hzero$", "lower bound for N H" = "lower bound $\\underline{N\\Hov}$",
            "finite-sample interval" = "interval of \\cref{prop:ebci}", "accuracy bound" = "accuracy bound (\\cref{cor:samplesize})",
            "T (t_T)" = "$\\Tp$ ($t_\\Tp$)", "V_jack" = "$\\Vjack$", "I_max" = "$\\Imax$", "G_max" = "$\\Gr$",
            "decision" = "decision", "tau-hat" = "$\\widehat\\tau$", "n_ovP-hat" = "$\\neffPhat$")
  tex_rule <- function(x) {
    x <- gsub(">=", "$\\ge$", x, fixed = TRUE); x <- gsub("<=", "$\\le$", x, fixed = TRUE)
    x <- gsub("|t_T|", "$|t_\\Tp|$", x, fixed = TRUE); x
  }
  vals <- rep1$value
  vals[rep1$quantity == "decision"] <- if (grepl("^report", d1$decision)) "report $\\log\\Zhat\\pm z_{0.975}\\sighat$" else "do not report the normal interval"
  vals <- gsub("^(-?[0-9.]+) \\((-?[0-9.]+)\\)$", "$\\1$ ($\\2$)", vals)
  vals <- ifelse(grepl("^-?[0-9.]+$", vals), paste0("$", vals, "$"), vals)
  vals <- ifelse(grepl("^\\[", vals), paste0("$", vals, "$"), vals)
  emit_table("tab_printout", cbind(qlab[rep1$quantity], vals, tex_rule(rep1$rule), rep1$status),
             spec = "llll", header = "Quantity & Value & Rule & Status\\\\")
  res_put("es.first.decision", d1$decision)
  res_put("es.first.bslogml.diff", fsci(abs(bs1$logml - d1$logml), 0))
}


## macros for Figure 1 (values and thresholds exactly as printed by bridge_report())
if (exists("rep1")) {
  v <- setNames(rep1$value, rep1$quantity); ru <- setNames(rep1$rule, rep1$quantity)
  last_num <- function(x) { m <- regmatches(x, gregexpr("-?[0-9]+([.][0-9]+)?", x))[[1]]; if (length(m)) m[length(m)] else "--" }
  key <- c("log Z-hat" = "logml", "MCSE" = "mcse", "n_ov-hat" = "nov", "I-hat" = "I", "N H0-hat" = "NH0",
           "lower bound for N H" = "NHlow", "finite-sample interval" = "eb", "accuracy bound" = "acc",
           "V_jack" = "Vjack", "I_max" = "Imax", "G_max" = "Gmax")
  for (q in intersect(names(key), names(v))) res_put(paste0("pr.", key[[q]]), v[[q]])
  for (q in c("n_ov-hat", "V_jack", "I_max", "G_max", "T (t_T)")) if (q %in% names(ru)) res_put(paste0("pr.thr.", gsub("[^A-Za-z]", "", q)), last_num(ru[[q]]))
  tt <- regmatches(v[["T (t_T)"]], gregexpr("-?[0-9.]+(e-?[0-9]+)?", v[["T (t_T)"]]))[[1]]
  res_put("pr.T", tt[1]); res_put("pr.tT", tt[2])
  half <- stats::qnorm(0.975) * d1$mcse
  res_put("pr.half", fx(half, 4)); res_put("pr.lo", fx(d1$logml - half, 4)); res_put("pr.hi", fx(d1$logml + half, 4))
  res_put("pr.S", fint(8000)); res_put("pr.G", fint(d1$settings$G))
  print(rep1)
}
res_write("printout.R")
