#!/usr/bin/env Rscript
# Section 7.2 / Figure 2: Markov chain draws.
# Top: exact laws for the two-state chain of Theorem 7 (eps = 0.05); the table of
# these computations is supplement/S09_E8_two_state_chain.R (Table E8).
# Bottom and Supplement Table E9: stationary AR(1) posterior chains for F1
# (d = 10, sigma = 1.5) and F2 (mu_g = 2, 4) with rho in {0, 0.5, 0.8, 0.95}; S chosen
# from the exact tau_W (spectral expansions) so that n_ovP is 1, 2, 4, 8, 20, 80
# 
# Expected: true-sigma_NP intervals cover 0.937-0.971; spectral-MCSE coverage is
# within 0.015 of nominal for n_ovP >= 20; error_measures() undercovers
# (0.59-0.88); V_jack flags an MCSE that ignores autocorrelation.
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else if (dir.exists("analysis")) "analysis" else "."})
source(file.path(BO_ROOT, "R", "common.R"))
source(file.path(BO_ROOT, "R", "families.R"))
source(file.path(BO_ROOT, "R", "exact.R"))
source(file.path(BO_ROOT, "R", "ar1.R"))
source(file.path(BO_ROOT, "R", "supplement.R"))
source(file.path(BO_ROOT, "R", "paper.R"))

## Top: two-state chain, exact ------------------------------------------------------
eps <- 0.05
top <- do.call(rbind, lapply(quick_pick(c(1, 0.5, 0.2, 0.05)), function(gam) {
  n_ovP <- quick_pick(exp(seq(log(0.5), log(100), length.out = 14)))
  S_vals <- unique(round(2 * n_ovP / (gam * eps * (1 - eps))))
  msg("two-state chain, gamma =", gam)
  do.call(rbind, lapply(S_vals, function(S) two_state_chain(S, eps, gam)))
}))
top$approx_ignoring <- 2 * stats::pnorm(1.96 * sqrt(top$gamma)) - 1
save_table(top, "fig2_two_state_chain")

## Bottom: AR(1) chains ----------------------------------------------------------------
R <- n_reps(REPS$e9)
G <- 20L
alpha <- 0.05
qa <- jackknife_null_quantiles(G, probs = 1 - alpha / 4)
logZ <- -0.4
fams <- list(family_F2(2), family_F2(4), family_F1(10, 1.5))   # order of Table E9
rows <- list()
for (fam in quick_pick(fams, c(1L, 3L))) {
  pop <- overlap_population(fam)
  per_S <- pop$var_p + pop$var_g
  for (rho in quick_pick(c(0, 0.5, 0.8, 0.95), c(1L, 3L))) {
    tw <- ar1_tau(fam, pop, rho)
    for (target in quick_pick(c(1, 2, 4, 8, 20, 80), c(1L, 5L))) {
      # Table E9: n_ovP equals the header values exactly, for every family
      eff <- target
      S <- max(4L, as.integer(round(eff * tw$tau_W / per_S)))
      v_P <- S * (pop$var_g + pop$var_p * tw$tau_P)
      sigma_NP <- sqrt(v_P) / (2 * S * pop$H)
      sigma_N <- overlap_at(pop, S)$sigma_N
      msg(fam$label, " rho =", rho, " target n_ovP =", target, " S =", S)
      res <- run_reps(R, function(r) {
        u_post <- ar1_u(S, fam$d, rho, 1L)[, 1]
        ap <- logZ + fam$lr(u_post)
        aq <- logZ + fam$lr(fam$rprop(S))
        d <- bridge_diagnostics(ap, aq, iid = FALSE, tau_method = "within", rule = "plain",   # flag: n_ovP-hat < 30
                                keep_responsibilities = TRUE, jackknife = S >= 2L * G)
        tau_ar <- iat(d$responsibilities$posterior, S, "concat_ar")
        mcse_ar <- sqrt(S * d$V_prop + S * d$V_post * tau_ar) / d$NH0
        re2 <- error_measures_re2(ap, aq, d$logml)
        data.frame(Delta = d$logml - logZ, mcse = d$mcse, mcse_iid = d$mcse_iid,
                   mcse_ar = mcse_ar, cv_em = sqrt(re2), n_P = d$n_P, n_ovP_hat = d$n_ovP,
                   I_hat = d$I, V_jack = d$V_jack, V_jack_iid = d$V_jack * (d$mcse / d$mcse_iid)^2,
                   flag_overlap = d$flags[["overlap"]])
      })
      cv <- function(se) mean(is.finite(res$Delta) & abs(res$Delta) <= 1.96 * se)
      thr_v <- stats::qchisq(1 - alpha / 4, G - 1) / (G - 1)
      rows[[length(rows) + 1L]] <- data.frame(
        family = fam$label, rho = rho, gamma = 1 - rho, target_n_ovP = target, S = S,
        tau_W = tw$tau_W, n_ovP = 2 * S * pop$H * pop$beta / tw$tau_W, R = nrow(res),
        coverage_spectral = cv(res$mcse), coverage_true_sigmaNP = cv(sigma_NP),
        coverage_ignoring = cv(sigma_N), approx_ignoring = 2 * stats::pnorm(1.96 / sqrt(tw$tau_W)) - 1,
        coverage_ar = cv(res$mcse_ar), coverage_error_measures = cv(res$cv_em),
        med_nP_over_vP = stats::median(res$n_P / v_P), med_n_ovP_hat = stats::median(res$n_ovP_hat),
        flag_n_ovP_lt_30 = mean(res$flag_overlap),
        Vjack_reject = mean(res$V_jack > thr_v, na.rm = TRUE),
        Vjack_reject_ignoring = mean(res$V_jack_iid > thr_v, na.rm = TRUE))
    }
  }
}
bottom <- do.call(rbind, rows)
save_table(bottom, "fig2_ar1_chains")
fam_lab <- c("F2 (mu_g=2)" = "F2, $\\mu_g=2$", "F2 (mu_g=4)" = "F2, $\\mu_g=4$", "F1 (d=10, sigma=1.5)" = "F1, $d=10$, $\\sigma=1.5$")
targets <- c(1, 2, 4, 8, 20, 80)
cells <- function(f) {
  keys <- unique(bottom[c("family", "rho")])
  t(vapply(seq_len(nrow(keys)), function(i) {
    x <- bottom[bottom$family == keys$family[i] & bottom$rho == keys$rho[i], ]
    f(x)
  }, character(if (identical(f, top_f)) 11 else 12)))
}
top_f <- function(x) c(fam_lab[[x$family[1]]], format(x$rho[1]), fx(x$tau_W[1], 2),
                       fx(x$coverage_spectral[match(targets, x$target_n_ovP)], 3),
                       frange(x$coverage_true_sigmaNP, 3), fx(stats::median(x$coverage_ignoring), 3))
bot_f <- function(x) c(fam_lab[[x$family[1]]], format(x$rho[1]), fx(x$tau_W[1], 2),
                       fx(x$med_nP_over_vP[match(targets, x$target_n_ovP)], 2),
                       fx(x$med_n_ovP_hat[match(c(1, 2, 4), x$target_n_ovP)], 1))
g <- group_ends(unique(bottom[c("family", "rho")])$family)
emit_table("tab_E9_top", cells(top_f), groups = g)
emit_table("tab_E9_bottom", cells(bot_f), groups = g)
pos <- bottom$rho > 0
res_put("e9.R", fint(R)); res_put("e9.ncfg", fint(nrow(bottom)))
res_put("e9.oracle", frange(bottom$coverage_true_sigmaNP, 3))
res_put("e9.ign", frange(bottom$coverage_ignoring, 2))
res_put("e9.cov.small", frange(bottom$coverage_spectral[pos & bottom$target_n_ovP == 1], 2))
res_put("e9.dev20", fx(max(abs(bottom$coverage_spectral[bottom$target_n_ovP >= 20] - 0.95)), 3))
res_put("e9.dev80", fx(max(abs(bottom$coverage_spectral[bottom$target_n_ovP >= 80] - 0.95)), 3))
res_put("e9.acc.small", frange(bottom$med_nP_over_vP[pos & bottom$target_n_ovP == 1], 2))
res_put("e9.acc.ge20", fx(min(bottom$med_nP_over_vP[bottom$target_n_ovP >= 20]), 2))
sm <- bottom$target_n_ovP <= 4
res_put("e9.overstate", fx(max(bottom$med_n_ovP_hat[sm] / bottom$n_ovP[sm]), 1))
gt <- unique(bottom[pos, c("rho", "tau_W", "family")])
res_put("e9.gtw", frange((1 - gt$rho) * gt$tau_W, 2))
res_write("02_markov_chains.R")

if (have_ggplot()) {   # separate process: own macros, no replication log
  st <- system2(file.path(R.home("bin"), "Rscript"), shQuote(file.path(BO_ROOT, "plot_mcmc.R")))
  if (!identical(as.integer(st), 0L)) stop("plot_mcmc.R failed (exit status ", st, "); its figure and macros are stale or missing")
}
