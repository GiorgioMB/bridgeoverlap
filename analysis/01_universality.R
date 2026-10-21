#!/usr/bin/env Rscript
# Section 7.1 / Figure 1 and Supplement Tables E6, HD, UNBAL, OC: universality of the
# effective overlap with independent draws, and its limits.
#   base       F1-F3 (Table S-E6);
#   hd         d = 100, 500 with I << 0.1 (radial families: exact sampling of |theta|^2);
#   unbal      F2 (mu_g = 4) with S1 != S2 (unbalanced allocations);
#   threshold  cells with n_ov ~ 9 and many replications: per-run operating
#              characteristics of the certificate (pass rate, coverage of passing and
#              of flagged runs) under the continuous and the gated rule.
# Per-run records are saved to results/runs_univ.rds for the pooled conditional
# coverage curves of supplement/S98_operating_characteristics.R.
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else if (dir.exists("analysis")) "analysis" else "."})
source(file.path(BO_ROOT, "R", "common.R"))
source(file.path(BO_ROOT, "R", "families.R"))
source(file.path(BO_ROOT, "R", "exact.R"))
source(file.path(BO_ROOT, "R", "supplement.R"))
source(file.path(BO_ROOT, "R", "paper.R"))

logZ <- 0.7
z975 <- stats::qnorm(0.975)
fam_of <- list(
  "F2 (mu_g=2)" = function() family_F2(2), "F2 (mu_g=4)" = function() family_F2(4),
  "F2 (mu_g=6)" = function() family_F2(6),
  "F1 (d=10, sigma=1.5)" = function() family_F1(10, 1.5), "F1 (d=40, sigma=1.3)" = function() family_F1(40, 1.3),
  "F3 (d=20)" = function() family_F3(20), "F3 (d=100)" = function() family_F3(100),
  "F1 (d=100, sigma=1.4)" = function() family_F1(100, 1.4), "F1 (d=500, sigma=1.2)" = function() family_F1(500, 1.2),
  "F3 (d=100, v=0.2)" = function() family_F3(100, 3, 0.2), "F3 (d=500, v=15)" = function() family_F3(500, 3, 15))
cfg <- rbind(
  # base grid of Table S-E6 (S1 = S2 = S)
  data.frame(group = "base", config = rep(c("F2 (mu_g=2)", "F2 (mu_g=4)", "F2 (mu_g=6)", "F1 (d=10, sigma=1.5)",
                                            "F1 (d=40, sigma=1.3)", "F3 (d=20)", "F3 (d=100)"), c(5, 4, 4, 4, 4, 4, 4)),
             S1 = c(10, 30, 100, 1000, 10000, 20, 60, 200, 2000, 300, 1000, 3000, 30000, 10, 30, 100, 1000,
                    100, 300, 1000, 10000, 10, 30, 100, 1000, 100, 300, 1000, 3000)),
  # high-dimensional configurations with I << 0.1: S chosen for n_ov from about 1 to 1000
  data.frame(group = "hd", config = rep(c("F1 (d=100, sigma=1.4)", "F1 (d=500, sigma=1.2)", "F3 (d=100, v=0.2)",
                                          "F3 (d=500, v=15)"), each = 5), S1 = NA_real_),
  # unbalanced allocations
  data.frame(group = "unbal", config = "F2 (mu_g=4)", S1 = c(30, 3000, 100, 10000, 300, 30000)),
  # near-threshold cells
  data.frame(group = "threshold", config = c("F2 (mu_g=4)", "F2 (mu_g=6)", "F1 (d=500, sigma=1.2)", "F3 (d=500, v=15)"),
             S1 = NA_real_))
cfg$S2 <- cfg$S1
cfg$S2[cfg$group == "unbal"] <- c(3000, 30, 10000, 100, 30000, 300)
target_hd <- c(1, 3, 10, 100, 1000)
pops <- list()
pop_of <- function(config, s1) {
  key <- paste(config, signif(s1, 6))
  if (is.null(pops[[key]])) pops[[key]] <<- overlap_population(fam_of[[config]](), s1)
  pops[[key]]
}
for (cf in unique(cfg$config[cfg$group == "hd"])) {           # S for target n_ov (S1 = S2)
  per <- 2 * pop_of(cf, 0.5)$H * pop_of(cf, 0.5)$beta
  i <- which(cfg$group == "hd" & cfg$config == cf)
  cfg$S1[i] <- cfg$S2[i] <- pmax(10, round(target_hd / per))
}
for (cf in unique(cfg$config[cfg$group == "threshold"])) {    # n_ov = 9 (S = 280 for F2, mu_g = 4)
  per <- 2 * pop_of(cf, 0.5)$H * pop_of(cf, 0.5)$beta
  i <- which(cfg$group == "threshold" & cfg$config == cf)
  cfg$S1[i] <- cfg$S2[i] <- if (cf == "F2 (mu_g=4)") 280 else round(9 / per)
}
rkey <- c(base = "univ", hd = "univ_hd", unbal = "univ_unbal", threshold = "univ_threshold")
cfg$R <- unlist(REPS[rkey[cfg$group]])
if (BO_QUICK) cfg <- cfg[c(1, 5, 14, 30, 34, 50, 52, 56, 57), ]
cfg <- cfg[!is.na(cfg$S1), ]

gate_NH <- 16 * log(2 / 0.05)
runs <- list(); rows <- list()
for (i in seq_len(nrow(cfg))) {
  cc <- cfg[i, ]
  fam <- fam_of[[cc$config]]()
  N <- cc$S1 + cc$S2; s1 <- cc$S1 / N
  pop <- pop_of(cc$config, s1)
  n_ov <- N * pop$H * pop$beta; sigma_N <- sqrt(pop$beta / (N * pop$H))
  R <- n_reps(cc$R, key = paste0("REPS$", rkey[[cc$group]]))
  msg(cc$group, cc$config, "S1 =", cc$S1, "S2 =", cc$S2, "n_ov =", signif(n_ov, 3), "I =", signif(pop$I, 3), "R =", R)
  res <- run_reps(R, function(r) {
    d <- bridge_diagnostics(logZ + fam$lr(fam$rpost(cc$S1)), logZ + fam$lr(fam$rprop(cc$S2)),
                            iid = TRUE, eb = FALSE, jackknife = FALSE,
                            rule = "continuous")   # flag_cont: the continuous rule
    data.frame(Delta = d$logml - logZ, mcse = d$mcse, n_ov_hat = d$n_ov, I_hat = d$I, NH0 = d$NH0,
               NH_lower = d$NH_lower, flag_cont = d$flags[["overlap"]])
  }, cell = paste("univ", cc$group, cc$config, cc$S1, cc$S2))
  ok <- is.finite(res$Delta) & is.finite(res$mcse)
  res$cover <- ok & abs(res$Delta) <= z975 * res$mcse
  res$flag_gate <- ifelse(res$I_hat <= 0.5, !(res$n_ov_hat >= 10), res$NH_lower < gate_NH)
  res$flag_plain <- !((res$n_ov_hat >= 10) %in% TRUE)   # the rule of Algorithm 1: pass rates below
  runs[[i]] <- cbind(group = cc$group, config = cc$config, family = fam$id, d = fam$d, S1 = cc$S1, S2 = cc$S2,
                     n_ov = n_ov, I = pop$I, res)
  pass <- !res$flag_plain
  rows[[i]] <- data.frame(
    group = cc$group, family = fam$id, config = cc$config, d = fam$d, S1 = cc$S1, S2 = cc$S2, H = pop$H,
    n_ov = n_ov, I = pop$I, R = nrow(res), coverage = mean(res$cover),
    ks = ks_normal(res$Delta / sigma_N), ks_noise = 0.8687 / sqrt(nrow(res)),
    sd_ratio = stats::sd(res$Delta[ok]) / sigma_N,
    med_abs_log_mcse = stats::median(ifelse(ok, abs(log(res$mcse / sigma_N)), Inf)),
    med_n_ov_ratio = stats::median(res$n_ov_hat) / n_ov, med_I_hat = stats::median(res$I_hat),
    pass_rate = mean(pass), cov_pass = if (any(pass)) mean(res$cover[pass]) else NA,
    cov_flag = if (any(!pass)) mean(res$cover[!pass]) else NA,
    pass_rate_gate = mean(!res$flag_gate), disagree = mean(res$flag_gate != res$flag_cont),
    I_hat_above = mean(res$I_hat > 0.5))
}
tab <- do.call(rbind, rows)
tab$C <- tab$med_abs_log_mcse * sqrt(tab$n_ov)
saveRDS(do.call(rbind, runs), file.path(results_dir, "runs_univ.rds"))
save_table(tab, "tab_E6")

## LaTeX outputs -------------------------------------------------------------------------------------
par_lab <- function(x) {
  x <- sub("^F[0-9] \\((.*)\\)$", "\\1", x)
  x <- gsub("mu_g=", "$\\\\mu_g=", x); x <- gsub("sigma=", "$\\\\sigma=", x); x <- gsub("d=", "$d=", x); x <- gsub("v=", "$v=", x)
  x <- gsub("(\\$[^,]*)(,|$)", "\\1$\\2", x); x
}
base <- tab[tab$group %in% c("base", "hd"), ]
emit_table("tab_E6", cbind(base$family, par_lab(base$config), fint(base$S1), fe(base$H), fx(base$I, 3), fx(base$n_ov, 2),
                           fx(base$coverage, 3), fx(base$ks, 3), fx(base$ks / base$ks_noise, 1), fx(base$sd_ratio, 3),
                           fx(base$med_abs_log_mcse, 3), fx(base$C, 2), fx(base$med_n_ov_ratio, 3)),
           groups = group_ends(paste(base$group, base$config)), spec = "llrrrrcccccccc",
           header = "Family & Parameters & $S$ & $\\Hov$ & $\\Iov$ & $\\neff$ & coverage & KS & KS/noise & $\\mathrm{sd}(\\Delta)/\\sigN$ & med $|\\log(\\sighat/\\sigN)|$ & $C$ & med $\\neffhat/\\neff$\\\\")
ub <- tab[tab$group == "unbal", ]
emit_table("tab_unbal", cbind(fint(ub$S1), fint(ub$S2), fx(ub$I, 3), fx(ub$n_ov, 1), fx(ub$coverage, 3),
                              fx(ub$ks, 3), fx(ub$med_abs_log_mcse, 3), fx(ub$med_I_hat, 3), fx(100 * ub$pass_rate, 0)),
           spec = "rrcrccccr",
           header = "$\\Sone$ & $\\Stwo$ & $\\Iov$ & $\\neff$ & coverage & KS & med $|\\log(\\sighat/\\sigN)|$ & med $\\Ihat$ & pass (\\%)\\\\")
th <- tab[tab$group == "threshold", ]
se <- function(p, n) sqrt(pmax(p * (1 - p), 1e-12) / n)
emit_table("tab_OC", cbind(paste0(th$family, ", ", par_lab(th$config)), fint(th$S1), fx(th$n_ov, 1), fx(th$I, 3), fint(th$R), fx(100 * th$pass_rate, 0),
                           fvse(th$cov_pass, se(th$cov_pass, th$R * th$pass_rate), 3),
                           fvse(th$cov_flag, se(th$cov_flag, th$R * (1 - th$pass_rate)), 3),
                           fx(th$coverage, 3)),
           spec = "lrrcrrccc",
           header = "Configuration & $S$ & $\\neff$ & $\\Iov$ & runs & pass (\\%) & coverage, passing & coverage, flagged & coverage, all\\\\")

R0 <- n_reps(REPS$univ)
res_put("univ.R", fint(R0)); res_put("univ.R.thr", fint(n_reps(REPS$univ_threshold)))
res_put("univ.S.min", fint(min(base$S1))); res_put("univ.S.max", fint(max(base$S1)))
res_put("univ.nov.min", fx(min(base$n_ov), 1)); res_put("univ.nov.max", fsci(max(base$n_ov), 1))
o <- order(base$n_ov); bad <- abs(base$coverage[o] - 0.95) > 0.01
res_put("univ.within", if (!any(bad)) fx(min(base$n_ov), 1) else if (max(which(bad)) < length(o)) fx(base$n_ov[o][max(which(bad)) + 1], 1) else "--")
mid <- base$n_ov >= 1 & base$n_ov <= 2.1
res_put("univ.cov12", frange(base$coverage[mid], 2))
res_put("univ.cov.min", fx(min(base$coverage), 2)); res_put("univ.cov.min.nov", fx(base$n_ov[which.min(base$coverage)], 1))
big <- base$n_ov >= 3
for (fm in c("F1", "F2", "F3")) res_put(paste0("univ.C.", fm), frange(base$C[big & base$family == fm & base$group == "base"], 2))
res_put("univ.C.hd", frange(base$C[big & base$group == "hd"], 2))
res_put("univ.C", frange(base$C[big], 2))
res_put("univ.I.base", frange(tab$I[tab$group == "base"], 2))
res_put("univ.I.F1F3low", frange(tab$I[tab$group == "base" & tab$config %in% c("F1 (d=10, sigma=1.5)", "F3 (d=20)")], 2))
res_put("univ.I.hd", frange(tab$I[tab$group == "hd"], 3))
res_put("univ.ks.exceed", sum(base$ks > 3 * base$ks_noise & base$n_ov >= 1))
w <- which.max(ifelse(base$n_ov >= 1, base$ks / base$ks_noise, -Inf))
res_put("univ.ks.worst", sprintf("%s at $\\neff=%s$ (%s times the noise level)", fx(base$ks[w], 3), fx(base$n_ov[w], 2), fx(base$ks[w] / base$ks_noise[w], 1)))
hd <- tab[tab$group == "hd", ]
res_put("univ.hd.cov", frange(hd$coverage[hd$n_ov >= 10], 3)); res_put("univ.hd.cov.small", frange(hd$coverage[hd$n_ov < 10], 2))
res_put("univ.unbal.cov", frange(ub$coverage, 3)); res_put("univ.unbal.nov", frange(ub$n_ov, 0))
res_put("univ.unbal.I", frange(ub$med_I_hat, 2))
for (k in seq_len(nrow(th))) {
  key <- c("F2 (mu_g=4)" = "f2m4", "F2 (mu_g=6)" = "f2m6", "F1 (d=500, sigma=1.2)" = "f1hd", "F3 (d=500, v=15)" = "f3hd")[[th$config[k]]]
  res_put(sprintf("univ.thr.%s.pass", key), fx(100 * th$pass_rate[k], 0))
  res_put(sprintf("univ.thr.%s.covpass", key), fx(th$cov_pass[k], 3)); res_put(sprintf("univ.thr.%s.covflag", key), fx(th$cov_flag[k], 3))
}
nb <- tab[tab$group == "base" & abs(tab$I - 0.5) < 0.1, ]

# Exact values for the two-point family of Theorem 2 (lines in Figure 1): eps in {0.05, 0.5}.
msg("two-point family (exact)")
tp <- do.call(rbind, lapply(c(0.05, 0.5), function(eps) {
  n_ov_grid <- quick_pick(exp(seq(log(0.5), log(1200), length.out = 30)))
  do.call(rbind, lapply(unique(round(2 * n_ov_grid / (eps * (1 - eps)))), function(S) two_point_iid(S, eps)))
}))
save_table(tp, "fig1_universality_two_point")
res_write("01_universality.R")

if (have_ggplot()) {   # separate process: own macros, no replication log
  st <- system2(file.path(R.home("bin"), "Rscript"), shQuote(file.path(BO_ROOT, "plot_universality.R")))
  if (!identical(as.integer(st), 0L)) stop("plot_universality.R failed (exit status ", st, "); its figure and macros are stale or missing")
}
