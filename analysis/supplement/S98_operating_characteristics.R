#!/usr/bin/env Rscript
# Section 7.1 / Figure 3 and Supplement Table OC-pooled: operating characteristics of the
# certificate per run. From the per-run records of 01_universality.R (all groups),
# S05_ADV_adversarial.R and 05_eight_schools.R (exact draws), the script computes the
# pooled conditional coverage P(cover | the run passes at threshold c), where a run passes
# at c if n_ov-hat >= c * min(1, 2 (1 - I-hat)) (continuous rule) or, for the gated rule,
# if n_ov-hat >= c when I-hat <= 1/2; and the coverage binned by n_ov-hat. Pooled values
# depend on the mix of configurations pooled and are reported per source.
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(dirname(normalizePath(f))) else if (dir.exists("analysis")) "analysis" else ".."})
for (f in c("common.R", "paper.R")) source(file.path(BO_ROOT, "R", f))
rd <- function(f) { p <- file.path(results_dir, f); if (file.exists(p)) readRDS(p) else NULL }
u <- rd("runs_univ.rds"); a <- rd("runs_adv.rds"); e <- rd("runs_eight_schools_exact.rds")
src <- list()
if (!is.null(u)) src$iid <- data.frame(source = "F1--F3, all designs", n = u$n_ov_hat, I = u$I_hat, cover = u$cover)
if (!is.null(a)) src$adv <- data.frame(source = "adversarial examples", n = a$n_ov, I = a$I, cover = a$cover)
if (!is.null(e)) {
  e <- e[grepl("^centred", e$cell), ]
  z <- stats::qnorm(0.975)
  src$es <- data.frame(source = "eight schools, centred", n = e$n_ov, I = e$I, cover = abs(e$err) <= z * e$mcse)
}
cs <- c(1, 2, 3, 5, 7, 10, 15, 20, 30, 50, 100)
curve <- do.call(rbind, lapply(src, function(d) {
  do.call(rbind, lapply(cs, function(c0) {
    pc <- d$n >= c0                                            # certified by the rule of Algorithm 1 at threshold c0
    data.frame(source = d$source[1], c = c0, runs = nrow(d),
               pass_rate = mean(pc), cov_pass = if (any(pc)) mean(d$cover[pc]) else NA, n_pass = sum(pc),
               cov_flag = if (any(!pc)) mean(d$cover[!pc]) else NA)
  }))
}))
save_table(curve, "tab_OC_curve")
bins <- c(0, 2, 5, 10, 20, 50, 100, Inf)
binned <- do.call(rbind, lapply(src, function(d) {
  b <- cut(d$n, bins, right = FALSE)
  data.frame(source = d$source[1], bin = levels(b), coverage = as.numeric(tapply(d$cover, b, mean)),
             runs = as.numeric(table(b)))
}))
save_table(binned, "tab_OC_binned")
sel <- curve$c %in% c(3, 5, 10, 20, 50)
cl <- function(p, n) ifelse(is.na(p), "--", sprintf("%s", fx(p, 3)))
emit_table("tab_OC_pooled", cbind(ifelse(!duplicated(curve$source[sel]), curve$source[sel], ""), fint(curve$c[sel]),
                                  fx(100 * curve$pass_rate[sel], 0), cl(curve$cov_pass[sel]), cl(curve$cov_flag[sel])),
           groups = group_ends(curve$source[sel]), spec = "lrrcc",
           header = "Source & $c$ & pass (\\%) & coverage, passing & coverage, flagged\\\\")
for (k in names(src)) {
  x <- curve[curve$source == src[[k]]$source[1], ]
  res_put(sprintf("oc.%s.cov10", k), fx(x$cov_pass[x$c == 10], 3)); res_put(sprintf("oc.%s.pass10", k), fx(100 * x$pass_rate[x$c == 10], 0))
  res_put(sprintf("oc.%s.covmin5", k), fx(min(x$cov_pass[x$c >= 5], na.rm = TRUE), 3))
  res_put(sprintf("oc.%s.runs", k), fint(x$runs[1]))
}
res_write("supplement/S98_operating_characteristics.R")

if (have_ggplot()) {   # separate process: own macros, no replication log
  st <- system2(file.path(R.home("bin"), "Rscript"), shQuote(file.path(BO_ROOT, "plot_conditional.R")))
  if (!identical(as.integer(st), 0L)) stop("plot_conditional.R failed (exit status ", st, "); its figure and macros are stale or missing")
}
