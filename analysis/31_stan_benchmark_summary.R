#!/usr/bin/env Rscript
# Section 7.3: summary of the Stan benchmark (runs of 30_stan_benchmark.R and of the eight
# schools study 10_stan_eight_schools.R). Writes
#   tab_bench    per configuration: accuracy, our flags, checks of current practice;
#   tab_detect   which check catches which failure, pooled over configurations;
#   tab_bf       log Bayes factors: conclusions of the default workflow and after the
#                decision rule (rerunning a flagged model with four times the draws);
#   tab_cost     compute cost.
# Checks of current practice, as a bridgesampling user would apply them:
#   R-hat/ESS    max split R-hat > 1.01, min bulk ESS < 400, or any divergence;
#   MCSE >= 0.3  the rule of Micaletto and Vehtari (2025), with the MCSE of
#                bridgesampling 1.2-1 (mcse_logml, with its effective-sample-size correction;
#                recomputed by mcse_logml_bs(), which reproduces it exactly); Table S-detect
#                also applies it to bridgesampling::error_measures() and to our within-chain
#                MCSE;
#   warp-III     |log Z-hat(normal) - log Z-hat(warp3)| > 1.96 sqrt(2) MCSE(error_measures);
#   repeated     two independent runs (consecutive replications) differ by more than
#                1.96 times the root sum of their squared error_measures() MCSEs.
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else if (dir.exists("analysis")) "analysis" else "."})
for (f in c("common.R", "stan.R", "paper.R")) source(file.path(BO_ROOT, "R", f))
z <- stats::qnorm(0.975)
files <- list.files(results_dir, pattern = "^bench_runs_.*[.]rds$", full.names = TRUE)
files <- files[basename(files) != "bench_runs_all.rds"]   # this script's own output: reading it back would double every run
drop_x <- c("sym_error", "sym_mcse", "sym_n_ov", "sym_I", "sym_pass", "time_symmetric", "chains_label1",
            "cf_error", "cf_mcse", "cf_n_ov", "cf_n_ovP", "cf_pass", "time_crossfit")
rbind_fill <- function(lst) {     # union of columns; runs from different versions or studies combine
  nm <- unique(unlist(lapply(lst, names)))
  do.call(rbind, lapply(lst, function(x) { for (n in setdiff(nm, names(x))) x[[n]] <- NA; x[, nm, drop = FALSE] }))
}
runs <- rbind_fill(lapply(files, function(f) { x <- readRDS(f); x[, setdiff(names(x), drop_x), drop = FALSE] }))
xcols <- c("config", "rep", "sym_error", "sym_mcse", "sym_n_ov", "sym_I", "sym_pass", "time_symmetric", "chains_label1",
           "cf_error", "cf_mcse", "cf_n_ov", "cf_n_ovP", "cf_pass", "time_crossfit")
extra <- do.call(rbind, lapply(files, function(f) {
  x <- readRDS(f); out <- data.frame(config = x$config, rep = x$rep)
  for (nm in xcols[-(1:2)]) out[[nm]] <- if (nm %in% names(x)) x[[nm]] else NA
  out
}))
es_file <- file.path(results_dir, "stan_eight_schools_runs.rds")
if (file.exists(es_file)) {
  es <- readRDS(es_file)
  es$config <- paste0("eight_schools_", ifelse(es$model == "centred", "c", "nc"), "_", sub("[.]", "", format(es$adapt_delta, nsmall = 2)), "_", es$S)
  runs <- rbind_fill(list(runs, es[, intersect(names(es), names(runs)), drop = FALSE]))
}
runs <- runs[order(runs$config, runs$rep), ]
runs$cover <- is.finite(runs$error) & abs(runs$error) <= z * runs$mcse
runs$cover_em <- is.finite(runs$error) & abs(runs$error) <= z * runs$cv_error_measures
# Algorithm 1 with the plain rule (Section 6.1), recomputed from the stored certificate and failure flags
runs$flag_overlap <- !(runs$n_ovP >= 30)
runs$pass <- !runs$flag_overlap & !(runs$flag_any_failure %in% TRUE)
runs$f_ours <- !runs$pass
runs$f_std <- runs$rhat_max > 1.01 | runs$ess_bulk_min < 400 | runs$divergent > 0
runs$f_em <- !(runs$cv_error_measures < 0.3)
if (is.null(runs$mcse_logml_bs)) runs$mcse_logml_bs <- NA_real_
runs$f_ml <- !(runs$mcse_logml_bs < 0.3)
runs$f_mcse <- !(runs$mcse < 0.3)
runs$f_w3 <- !is.finite(runs$logml_warp3) | abs(runs$logml_bs - runs$logml_warp3) > z * sqrt(2) * runs$cv_error_measures
runs$f_rep <- NA
for (cf in unique(runs$config)) {
  i <- which(runs$config == cf); k <- floor(length(i) / 2)
  if (k < 1) next
  a <- i[2 * seq_len(k) - 1]; b <- i[2 * seq_len(k)]
  f <- abs(runs$error[a] - runs$error[b]) > z * sqrt(runs$cv_error_measures[a]^2 + runs$cv_error_measures[b]^2)
  runs$f_rep[a] <- f; runs$f_rep[b] <- f
}
runs$f_cur <- (runs$f_std %in% TRUE) | (runs$f_ml %in% TRUE) | (runs$f_w3 %in% TRUE) | (runs$f_rep %in% TRUE)
saveRDS(runs, file.path(results_dir, "bench_runs_all.rds"))

pct <- function(x) 100 * mean(x, na.rm = TRUE)
covp <- function(cv, f) if (any(!f, na.rm = TRUE)) mean(cv[!f %in% TRUE]) else NA
summ <- do.call(rbind, lapply(split(runs, runs$config), function(x) data.frame(
  config = x$config[1], d = NA, R = nrow(x), coverage = mean(x$cover), coverage_em = mean(x$cover_em),
  bias_mcse = mean(x$error) / stats::median(x$mcse), rmse = sqrt(mean(x$error^2)), med_mcse = stats::median(x$mcse),
  med_n_ov = stats::median(x$n_ov), med_I = stats::median(x$I), med_n_ovP = stats::median(x$n_ovP),
  f_overlap = pct(x$flag_overlap), f_failure = pct(x$flag_any_failure), f_ours = pct(x$f_ours),
  f_std = pct(x$f_std), f_em = pct(x$f_em), f_ml = pct(x$f_ml), f_w3 = pct(x$f_w3), f_rep = pct(x$f_rep), f_cur = pct(x$f_cur),
  cov_pass_ours = covp(x$cover, x$f_ours), cov_pass_cur = covp(x$cover_em, x$f_cur),
  cov_pass_std = covp(x$cover_em, x$f_std))))
dims <- c(radon_mn_pooled = 3, radon_mn_vi_c = 89, radon_mn_vi_nc = 89, radon_mn_hi_c = 90, radon_mn_vis_c = 175,
          radon_mn_vis_nc = 175, radon_all_vi_c = 390, radon_mn_vi_c_x4 = 89, radon_mn_hi_c_x4 = 90,
          radon_mn_vis_c_x4 = 175, radon_mn_vis_nc_x4 = 175, mixture = 3)
summ$d <- ifelse(summ$config %in% names(dims), dims[summ$config], 10)
save_table(summ, "tab_bench")
lab <- function(id) {
  id <- sub("^radon_mn_", "radon MN, ", id); id <- sub("^radon_all_", "radon all, ", id)
  id <- sub("vis_c", "slopes, centred", id); id <- sub("vis_nc", "slopes, non-centred", id)
  id <- sub("vi_c", "intercepts, centred", id); id <- sub("vi_nc", "intercepts, non-centred", id)
  id <- sub("hi_c", "uranium, centred", id); id <- sub("pooled", "pooled", id)
  id <- sub("_x4$", ", $4\\\\times$ draws", id)
  es <- regmatches(id, regexec("^eight_schools_(c|nc)_([0-9]+)_([0-9]+)$", id))
  for (k in which(lengths(es) == 4L)) {
    e <- es[[k]]
    id[k] <- sprintf("eight schools, %s, $\\delta=%s$, $S=%s$", if (e[2] == "c") "centred" else "non-centred",
                     format(as.numeric(sub("^0", "0.", e[3]))), e[4])
  }
  id <- sub("^mixture$", "normal mixture", id)
  gsub("(?<![\\\\])_", "\\\\_", id, perl = TRUE)        # no raw underscore reaches LaTeX
}
emit_table("tab_bench", cbind(lab(summ$config), fint(summ$d), fint(summ$R), fx(summ$coverage, 3), fsign(summ$bias_mcse, 1),
                              fint(summ$med_n_ov), fx(summ$med_I, 2), fx(summ$f_overlap, 0), fx(summ$f_failure, 0),
                              fx(summ$f_std, 0), fx(summ$f_ml, 0), fx(summ$f_w3, 0), fx(summ$f_rep, 0),
                              fx(summ$cov_pass_ours, 3), fx(summ$coverage_em, 3)),
           spec = "lrrcrrcrrrrrrcc",
           header = paste0("& & & & & & & \\multicolumn{2}{c}{our flags (\\%)} & \\multicolumn{4}{c}{current practice (\\%)} & \\multicolumn{2}{c}{coverage}\\\\\n",
                           "\\cmidrule(lr){8-9}\\cmidrule(lr){10-13}\\cmidrule(lr){14-15}\n",
                           "Posterior & $d$ & runs & coverage & bias/MCSE & med $\\neffhat$ & med $\\Ihat$ & overlap & failure & $\\widehat R$/ESS & MCSE${}\\ge0.3$ & warp-III & repeated & passing & \\texttt{error\\_measures}\\\\"))

# Which check catches which failure (pooled): failure = our 95% interval misses log Z.
miss <- !runs$cover
checks <- c(f_ours = "certificate and failure flags", f_std = "$\\widehat R$, ESS, divergences",
            f_ml = "MCSE${}\\ge0.3$: \\texttt{mcse\\_logml} (1.2-1)", f_em = "MCSE${}\\ge0.3$: \\texttt{error\\_measures()}",
            f_mcse = "MCSE${}\\ge0.3$: within-chain MCSE",
            f_w3 = "warp-III disagrees", f_rep = "repeated runs disagree", f_cur = "any check of current practice")
det <- do.call(rbind, lapply(names(checks), function(k) {
  f <- runs[[k]] %in% TRUE; big <- abs(runs$error) > 0.5
  data.frame(check = checks[[k]], flag = pct(f), sens = if (any(miss)) pct(f[miss]) else NA,
             false_alarm = if (any(!miss)) pct(f[!miss]) else NA, sens_big = if (any(big)) pct(f[big]) else NA,
             cov_pass = if (any(!f)) mean(runs$cover[!f]) else NA, runs_pass = sum(!f))
}))
save_table(det, "tab_detect")
emit_table("tab_detect", cbind(det$check, fx(det$flag, 0), fx(det$sens, 0), fx(det$sens_big, 0), fx(det$false_alarm, 0),
                               fx(det$cov_pass, 3), fint(det$runs_pass)),
           spec = "lrrrrcr",
           header = "Check & flagged (\\%) & of misses (\\%) & of errors ${}>0.5$ (\\%) & of covering runs (\\%) & coverage, passing & runs passing\\\\")

# Bayes factors: replication r of model A against replication r of model B (independent runs).
cat_of <- function(x) cut(x, c(-Inf, -5, -3, -1, 1, 3, 5, Inf), right = FALSE)
bf_pairs <- list(c("radon_mn_vi_c", "radon_mn_pooled"), c("radon_mn_hi_c", "radon_mn_vi_c"),
                 c("radon_mn_vis_c", "radon_mn_vi_c"), c("radon_mn_vis_nc", "radon_mn_vi_nc"))
bf_runs <- list()
bf <- do.call(rbind, lapply(bf_pairs, function(pr) {
  A <- runs[runs$config == pr[1], ]; B <- runs[runs$config == pr[2], ]
  if (!nrow(A) || !nrow(B)) return(NULL)
  r <- intersect(A$rep, B$rep); A <- A[match(r, A$rep), ]; B <- B[match(r, B$rep), ]
  A4 <- runs[runs$config == paste0(pr[1], "_x4"), ]; B4 <- runs[runs$config == paste0(pr[2], "_x4"), ]
  truth <- A$truth[1] - B$truth[1]
  est <- (A$logml) - (B$logml)
  ok <- cat_of(est) == cat_of(truth)
  passA <- A$pass; passB <- B$pass
  # decision rule: a flagged run is replaced by the run with four times the draws. The 4x runs
  # exist for the first bench.R.x4 replications only, so the rule is evaluated on the pairs whose
  # replication has them; a model without a 4x configuration (pooling) cannot be rerun.
  estA <- A$logml; estB <- B$logml; certA <- passA; certB <- passB; r4 <- r
  if (nrow(A4)) { j <- match(r, A4$rep); use <- !passA & !is.na(j); estA[use] <- A4$logml[j[use]]; certA[use] <- A4$pass[j[use]]; r4 <- intersect(r4, A4$rep) }
  if (nrow(B4)) { j <- match(r, B4$rep); use <- !passB & !is.na(j); estB[use] <- B4$logml[j[use]]; certB[use] <- B4$pass[j[use]]; r4 <- intersect(r4, B4$rep) }
  est2 <- estA - estB; cert <- certA & certB; s <- r %in% r4
  bf_runs[[length(bf_runs) + 1L]] <<- data.frame(A = pr[1], B = pr[2], truth = truth, rep = r, est = est, pass = passA & passB,
                                                 after = s, est_after = est2, cert_after = cert)
  data.frame(A = pr[1], B = pr[2], truth = truth, pairs = length(r), med_est = stats::median(est),
             q10 = stats::quantile(est, 0.1), q90 = stats::quantile(est, 0.9), wrong_naive = pct(!ok),
             flagged_of_wrong = if (any(!ok)) pct(!(passA & passB)[!ok]) else NA,
             wrong_certified = if (any(passA & passB)) pct(!ok[passA & passB]) else NA,
             pairs_after = sum(s),
             wrong_after = if (any(cert & s)) pct(cat_of(est2[cert & s]) != cat_of(truth)) else NA,
             uncertified_after = pct(!cert[s]),
             changed = pct(cat_of(est2[s]) != cat_of(est[s])))
}))
if (!is.null(bf) && nrow(bf)) {
  save_table(bf, "tab_bf"); save_table(do.call(rbind, bf_runs), "bf_pairs")
  emit_table("tab_bf", cbind(paste(lab(bf$A), "vs", lab(bf$B)), fsign(bf$truth, 2), fint(bf$pairs),
                             sprintf("$%s$ [$%s$, $%s$]", fsign(bf$med_est, 2), fsign(bf$q10, 2), fsign(bf$q90, 2)),
                             fx(bf$wrong_naive, 0), fx(bf$flagged_of_wrong, 0), fx(bf$wrong_certified, 0),
                             fint(bf$pairs_after), fx(bf$wrong_after, 0), fx(bf$uncertified_after, 0), fx(bf$changed, 0)),
             spec = "lrrcrrrrrrr",
             header = "Comparison & true $\\log B$ & pairs & estimate, median [10\\%, 90\\%] & wrong (\\%) & flagged, of wrong (\\%) & wrong, both pass (\\%) & pairs after rule & wrong after rule (\\%) & not certified after rule (\\%) & conclusion changed (\\%)\\\\")
  vn <- bf[bf$A == "radon_mn_vis_nc", ]
  if (nrow(vn)) { res_put("bf.visnc.uncert", fx(vn$uncertified_after, 0)); res_put("bf.visnc.wrong", fx(vn$wrong_naive, 0)) }
  vs <- bf[bf$A == "radon_mn_vis_c", ]
  if (nrow(vs)) {
    res_put("bf.vis.truth", fsign(vs$truth, 2)); res_put("bf.vis.wrong", fx(vs$wrong_naive, 0))
    res_put("bf.vis.flagged", fx(vs$flagged_of_wrong, 0)); res_put("bf.vis.wrongafter", fx(vs$wrong_after, 0))
    res_put("bf.vis.changed", fx(vs$changed, 0)); res_put("bf.vis.uncert", fx(vs$uncertified_after, 0))
  }
  res_put("bf.pairs.after", fint(min(bf$pairs_after)))
  res_put("bf.wrong.max", fx(max(bf$wrong_naive), 0)); res_put("bf.wrongafter.max", fx(max(bf$wrong_after, na.rm = TRUE), 0))
}

# Compute cost (median seconds per run).
tcol <- function(x, nm) if (nm %in% names(x)) stats::median(x[[nm]], na.rm = TRUE) else NA
cost <- do.call(rbind, lapply(split(runs, runs$config), function(x) data.frame(
  config = x$config[1], sampling = tcol(x, "time_sampling"), bridge = tcol(x, "time_bridge"),
  diagnostics = tcol(x, "time_diagnostics"), warp3 = tcol(x, "time_warp3"))))
cost$overhead_pct <- 100 * cost$diagnostics / (cost$sampling + cost$bridge)
save_table(cost, "tab_cost")
emit_table("tab_cost", cbind(lab(cost$config), fx(cost$sampling, 1), fx(cost$bridge, 1), fx(cost$diagnostics, 2),
                             fx(cost$warp3, 1), fx(cost$overhead_pct, 1)),
           spec = "lrrrrr", header = "Posterior & sampling (s) & bridge sampling (s) & diagnostics (s) & warp-III (s) & diagnostics / (sampling + bridge) (\\%)\\\\")

# Mixture: label switching and the permutation-invariant proposal.
mx <- merge(runs[runs$config == "mixture", ], extra[extra$config == "mixture", ], by = c("config", "rep"))
if (nrow(mx)) {
  same <- mx$chains_label1 %in% c(0, 4)
  res_put("mix.R", fint(nrow(mx))); res_put("mix.same.pct", fx(100 * mean(same), 0))
  res_put("mix.same.bias", fsign(mean(mx$error[same]), 3)); res_put("mix.same.cov", fx(mean(mx$cover[same]), 2))
  res_put("mix.same.mcse", fx(pct(mx$mcse_logml_bs[same] < 0.3), 0)); res_put("mix.same.I", fx(stats::median(mx$I[same]), 2))
  res_put("mix.same.novP", fx(stats::median(mx$n_ovP[same]), 1)); res_put("mix.same.n", fint(sum(same)))
  res_put("mix.same.flag", fx(pct(mx$f_ours[same]), 0)); res_put("mix.same.std", fx(pct(mx$f_std[same]), 0))
  res_put("mix.mixed.flag", fx(pct(mx$f_ours[!same]), 0)); res_put("mix.mixed.std", fx(pct(mx$f_std[!same]), 0))
  res_put("mix.mixed.cov", fx(mean(mx$cover[!same]), 2))
  res_put("mix.sym.cov", fx(mean(abs(mx$sym_error) <= z * mx$sym_mcse, na.rm = TRUE), 3))
  res_put("mix.sym.bias", fsign(mean(mx$sym_error, na.rm = TRUE), 4)); res_put("mix.sym.pass", fx(pct(mx$sym_pass), 0))
}
cf <- merge(runs, extra, by = c("config", "rep"))
cf <- cf[is.finite(cf$cf_error), ]
if (nrow(cf)) {
  g <- do.call(rbind, lapply(split(cf, cf$config), function(x) data.frame(config = x$config[1],
          rmse_split = sqrt(mean(x$error^2)), rmse_cf = sqrt(mean(x$cf_error^2)),
          cov_cf = mean(abs(x$cf_error) <= z * x$cf_mcse))))
  save_table(g, "tab_bench_crossfit")
  res_put("bench.cf.gain", frange(100 * (1 - g$rmse_cf / g$rmse_split), 0)); res_put("bench.cf.cov", frange(g$cov_cf, 2))
}
for (cf in unique(runs$config)) res_put(paste0("bench.logZ.", gsub("_", "", cf)), fx(runs$truth[runs$config == cf][1], 4))
res_put("bench.R", fint(n_reps(REPS$stan_radon))); res_put("bench.R.big", fint(n_reps(REPS$stan_radon_big)))
res_put("bench.R.x4", fint(n_reps(REPS$stan_radon_x4))); res_put("bench.R.mix", fint(n_reps(REPS$stan_mix)))
res_put("bench.n.configs", nrow(summ)); res_put("bench.n.runs", fint(nrow(runs)))
res_put("bench.overhead.max", fx(max(cost$overhead_pct, na.rm = TRUE), 1)); res_put("bench.overhead.med", fx(stats::median(cost$overhead_pct, na.rm = TRUE), 1))
for (k in names(checks)) {
  x <- det[det$check == checks[[k]], ]
  res_put(sprintf("det.%s.sens", sub("^f_", "", k)), fx(x$sens, 0)); res_put(sprintf("det.%s.fa", sub("^f_", "", k)), fx(x$false_alarm, 0))
  res_put(sprintf("det.%s.cov", sub("^f_", "", k)), fx(x$cov_pass, 3))
}
# Workflows: which runs a user would report, and how their intervals cover. R-hat, ESS and
# divergences come first, as in the Stan workflow; then our checks, or the single-run checks of
# current practice (the MCSE rule with mcse_logml, warp-III); repeated runs need a second run.
k_std <- !(runs$f_std %in% TRUE)
wf_keep <- list(std = k_std, both = k_std & !runs$f_ours,
                one = k_std & !(runs$f_ml %in% TRUE) & !(runs$f_w3 %in% TRUE), cur = !runs$f_cur)
res_put("wf.all.cov", fx(mean(runs$cover), 3))
for (nm in names(wf_keep)) {
  k <- wf_keep[[nm]]; res_put(sprintf("wf.%s.rep", nm), fx(pct(k), 0)); res_put(sprintf("wf.%s.cov", nm), fx(mean(runs$cover[k]), 3))
  res_put(sprintf("wf.%s.nomix.cov", nm), fx(mean(runs$cover[k & runs$config != "mixture"]), 3))
}
res_put("wf.std.miss.caught", fx(pct(runs$f_ours[k_std & !runs$cover]), 0))
res_put("wf.std.cover.flagged", fx(pct(runs$f_ours[k_std & runs$cover]), 0))
ra <- runs$config == "radon_all_vi_c"
res_put("wf.radonall.stdpass", fx(100 - pct(runs$f_std[ra]), 0)); res_put("wf.radonall.cov", fx(mean(runs$cover[ra]), 2))
res_put("wf.radonall.covpct", fx(100 * mean(runs$cover[ra]), 0))
res_put("wf.radonall.ours", fx(pct(runs$f_ours[ra]), 0)); res_put("wf.radonall.ml", fx(pct(runs$f_ml[ra]), 0))
res_write("31_stan_benchmark_summary.R")
if (have_ggplot()) {   # separate process: own macros, no replication log
  st <- system2(file.path(R.home("bin"), "Rscript"), shQuote(file.path(BO_ROOT, "plot_bench.R")))
  if (!identical(as.integer(st), 0L)) stop("plot_bench.R failed (exit status ", st, ")")
}
