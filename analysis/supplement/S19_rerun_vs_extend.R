#!/usr/bin/env Rscript
# Supplement Table RERUN: what to do after an overlap flag. At the threshold cells of
# 01_universality.R (n_ov ~ 9), each replication starts with S0 posterior and S0 proposal
# draws and follows one of two policies until the continuous rule passes (at most K steps):
#   fresh   a new independent run with both samples multiplied by 1.5 times the factor
#           printed by Algorithm 1 (draws_factor = threshold / n_ov-hat);
#   extend  the same factor, but the new draws are added to the flagged run (optional
#           stopping).
# Both policies use their own random-number streams. The table reports the coverage of the
# intervals that are finally reported, the share of replications that end with a report,
# and the mean total number of draws, with two baselines: report at S0 whatever the flags,
# and report at S0 only the runs that pass.
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(dirname(normalizePath(f))) else if (dir.exists("analysis")) "analysis" else ".."})
for (f in c("common.R", "families.R", "paper.R")) source(file.path(BO_ROOT, "R", f))
logZ <- 0.7; z <- stats::qnorm(0.975); K <- 6L; margin <- 1.5
cells <- list(list(key = "f2m4", lab = "F2, $\\mu_g=4$", fam = function() family_F2(4), S0 = 280L),
              list(key = "f2m6", lab = "F2, $\\mu_g=6$", fam = function() family_F2(6), S0 = NA),
              list(key = "f1hd", lab = "F1, $d=500$, $\\sigma=1.2$", fam = function() family_F1(500, 1.2), S0 = NA),
              list(key = "f3hd", lab = "F3, $d=500$, $v=15$", fam = function() family_F3(500, 3, 15), S0 = NA))
if (BO_QUICK) cells <- cells[1:2]
diag1 <- function(ap, aq) {
  d <- bridge_diagnostics(ap, aq, iid = TRUE, eb = FALSE, jackknife = FALSE, rule = "plain")
  list(pass = !d$flags[["overlap"]], cover = is.finite(d$logml) && abs(d$logml - logZ) <= z * d$mcse,
       factor = if (is.finite(d$draws_factor)) d$draws_factor else 2)
}
rows <- list()
for (cc in cells) {
  fam <- cc$fam()
  if (is.na(cc$S0)) { pop <- overlap_population(fam, 0.5); cc$S0 <- as.integer(round(9 / (2 * pop$H * pop$beta))) }
  draw <- function(n) list(ap = logZ + fam$lr(fam$rpost(n)), aq = logZ + fam$lr(fam$rprop(n)))
  R <- n_reps(REPS$rerun)
  msg("RERUN:", cc$key, "S0 =", cc$S0, "R =", R)
  res <- run_reps(R, cell = paste("rerun", cc$key), fun = function(r) {
    # baselines and the first run (shared by both policies)
    x <- draw(cc$S0); d0 <- diag1(x$ap, x$aq)
    out <- data.frame(pass0 = d0$pass, cover0 = d0$cover)
    # fresh reruns
    S <- cc$S0; d <- d0; tot <- 2 * S; k <- 1L
    while (!d$pass && k < K) { S <- ceiling(margin * d$factor * S); y <- draw(S); d <- diag1(y$ap, y$aq); tot <- tot + 2 * S; k <- k + 1L }
    out$fresh_rep <- d$pass; out$fresh_cover <- d$cover; out$fresh_draws <- tot; out$fresh_steps <- k
    # extension of the flagged run
    S <- cc$S0; d <- d0; ap <- x$ap; aq <- x$aq; k <- 1L
    while (!d$pass && k < K) {
      Snew <- ceiling(margin * d$factor * S); y <- draw(Snew - S)
      ap <- c(ap, y$ap); aq <- c(aq, y$aq); S <- Snew; d <- diag1(ap, aq); k <- k + 1L
    }
    out$ext_rep <- d$pass; out$ext_cover <- d$cover; out$ext_draws <- 2 * S; out$ext_steps <- k
    out
  })
  se <- function(p, n) sqrt(p * (1 - p) / n)
  pol <- function(nm, cover, rep, draws) {
    n <- sum(rep); p <- if (n) mean(cover[rep]) else NA
    data.frame(cell = cc$key, label = cc$lab, S0 = cc$S0, policy = nm, runs = nrow(res), reported = mean(rep),
               coverage = p, cov_se = if (n) se(p, n) else NA, draws = if (is.null(draws)) 2 * cc$S0 else mean(draws[rep]))
  }
  rows[[cc$key]] <- rbind(pol("report at $S_0$", res$cover0, rep(TRUE, nrow(res)), NULL),
                          pol("report at $S_0$ if it passes", res$cover0, res$pass0, NULL),
                          pol("fresh rerun", res$fresh_cover, res$fresh_rep, res$fresh_draws),
                          pol("extend the flagged run", res$ext_cover, res$ext_rep, res$ext_draws))
}
tab <- do.call(rbind, rows)
save_table(tab, "tab_rerun")
emit_table("tab_rerun", cbind(ifelse(!duplicated(tab$cell), tab$label, ""), ifelse(!duplicated(tab$cell), fint(tab$S0), ""),
                              tab$policy, fx(100 * tab$reported, 1), fvse(tab$coverage, tab$cov_se, 3), fint(round(tab$draws))),
           groups = group_ends(tab$cell), spec = "lrlrcr",
           header = "Configuration & $S_0$ & policy & reported (\\%) & coverage of reported intervals & mean draws\\\\")
fr <- tab[tab$policy == "fresh rerun", ]; ex <- tab[tab$policy == "extend the flagged run", ]
res_put("rr.R", fint(n_reps(REPS$rerun))); res_put("rr.cov.fresh", frange(fr$coverage, 3)); res_put("rr.cov.extend", frange(ex$coverage, 3))
dz <- (ex$coverage - fr$coverage) / sqrt(ex$cov_se^2 + fr$cov_se^2)
res_put("rr.phrase", if (all(abs(dz) <= 2)) "the same coverage as fresh reruns within Monte Carlo error" else if (all(dz < 2) && any(dz < -2)) "lower coverage than fresh reruns" else if (all(dz > -2) && any(dz > 2)) "higher coverage than fresh reruns" else "coverage that differed from that of fresh reruns in both directions")
res_put("rr.draws.ratio", frange(ex$draws / fr$draws, 2))
res_write("supplement/S19_rerun_vs_extend.R")
