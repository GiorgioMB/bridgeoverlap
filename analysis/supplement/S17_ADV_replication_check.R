#!/usr/bin/env Rscript
# Supplement Table ADV-check: every row of Table S-ADV rerun with an independent random
# number stream and the same pre-specified number of replications. The coverages of the
# two streams are compared by z = (c1 - c2) / sqrt(se1^2 + se2^2); under correct code the
# z values are approximately standard normal, so a row with |z| > 3 would point to a bug
# (e.g. replications that are not independent) rather than to chance.
ADV_STREAM <- "check"; ADV_DEFINE_ONLY <- TRUE
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(dirname(normalizePath(f))) else if (dir.exists("analysis")) "analysis" else ".."})
source(file.path(BO_ROOT, "supplement", "S05_ADV_adversarial.R"))   # defines `examples` only
main <- readRDS(file.path(results_dir, "runs_adv.rds"))
rows <- list()
for (ex in quick_pick(examples, c(1L, 3L, 5L))) for (S in quick_pick(ex$S)) {
  m <- main[main$example == ex$name & main$S == S, ]
  if (!nrow(m)) next
  msg("ADV check:", ex$name, "S =", S)
  res <- run_reps(ex$R_of(S), cell = paste("adv", ADV_STREAM, ex$name, S), fun = function(r) {
    a <- ex$draw(S)
    d <- bridge_diagnostics(a$ap, a$aq, iid = TRUE, eb = FALSE, jackknife = FALSE)
    Delta <- d$logml - logZ
    data.frame(cover = is.finite(Delta) && isTRUE(abs(Delta) <= z975 * d$mcse))
  })
  c1 <- mean(m$cover); c2 <- mean(res$cover)
  se <- sqrt(c1 * (1 - c1) / nrow(m) + c2 * (1 - c2) / nrow(res))
  rows[[length(rows) + 1L]] <- data.frame(example = ex$name, S = S, R1 = nrow(m), R2 = nrow(res),
                                          cov_main = c1, cov_check = c2, z = if (se > 0) (c1 - c2) / se else 0)
}
tab <- do.call(rbind, rows)
save_table(tab, "tab_ADV_check")
lab <- c("smooth two-point, mu0 = 0" = "smooth two-point, $\\mu_0=0$", "smooth two-point, mu0 = 0.5" = "smooth two-point, $\\mu_0=0.5$",
         "rare mode carries overlap" = "rare mode carries overlap",
         "funnel, d = 5, sigma_v = 1.5" = "funnel, $d=5$, $\\sigma_v=1.5$", "funnel, d = 10, sigma_v = 3" = "funnel, $d=10$, $\\sigma_v=3$")
emit_table("tab_ADV_check", cbind(ifelse(!duplicated(tab$example), lab[tab$example], ""), tab$S, fint(tab$R1),
                                  fx(tab$cov_main, 3), fx(tab$cov_check, 3), fsign(tab$z, 2)),
           groups = group_ends(tab$example), spec = "lrrccr",
           header = "Example & $S$ & runs per stream & coverage, stream 1 & coverage, stream 2 & $z$\\\\")
res_put("advchk.maxz", fx(max(abs(tab$z)), 2)); res_put("advchk.n", nrow(tab))
res_put("advchk.n3", sum(abs(tab$z) > 3)); res_put("advchk.sdz", fx(stats::sd(tab$z), 2))
res_write("supplement/S17_ADV_replication_check.R")
