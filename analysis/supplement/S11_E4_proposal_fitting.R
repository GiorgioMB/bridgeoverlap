#!/usr/bin/env Rscript
# Supplement Table E4: proposal fitting (t_7 target, diagonal Gaussian proposal with
# k = 2d fitted parameters, G = 20). Standard deviations of log Z-hat with the proposal
# fixed at its limit, fitted on the same draws, and fitted on an independent sample
# (split); bias in MCSE units; the diagnostic R_g with the proposal draws held fixed and
# its pathwise version with common random numbers (not for d = 32); mean V_jack for
# three jackknife variants of the same-draw estimator: no refitting, refitting with the
# proposal draws held fixed, and refitting with common random numbers.
#
# R_g is not defined explicitly in the supplement. Here it is the ratio of the
# delta-method standard error contributed by the fitted parameters to the MCSE,
#   R_g = sqrt(grad' Sigma grad) / MCSE,  grad = d c-hat / d theta (implicit function
#   theorem on W), Sigma = covariance of the moment estimates,
# which is O(1) when the proposal draws are held fixed and O(n_ov^(-1/2)) for the
# pathwise derivative, the two behaviours the supplement reports.
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(dirname(normalizePath(f))) else if (dir.exists("analysis")) "analysis" else ".."})
for (f in c("common.R", "families.R", "exact.R", "ar1.R", "supplement.R", "paper.R")) source(file.path(BO_ROOT, "R", f))
source(file.path(BO_ROOT, 'R', 'proposal_fitting.R'))

G <- 20L
cfg <- data.frame(d = c(2L, 8L, 32L, 2L), S = c(2000L, 2000L, 2000L, 8000L), R = REPS$e4)
if (BO_QUICK) cfg <- cfg[c(1, 4), ]
jack_refit <- function(tg, X, E, Yfixed, block, G) {
  vapply(c(fixed = 1, crn = 2), function(mode) {
    cb <- vapply(seq_len(G), function(b) {
      keep <- block != b
      mm <- moments(X[keep, , drop = FALSE])
      Yb <- if (mode == 1) Yfixed[keep, , drop = FALSE] else prop_from(E[keep, , drop = FALSE], mm$m, mm$v)
      a <- ratios(tg, X[keep, , drop = FALSE], Yb, mm$m, mm$v)
      bridge_logml(a$ap, a$aq)$logml
    }, numeric(1))
    (G - 1) / G * sum((cb - mean(cb))^2)
  }, numeric(1))
}
R_g <- function(tg, X, Y, E, m, v, d_obj, Xfit) {
  ch <- d_obj$logml
  a <- ratios(tg, X, Y, m, v)
  wp <- stats::dlogis(ch - a$ap); wq <- stats::dlogis(ch - a$aq)
  NH <- sum(wp) + sum(wq)
  dlg <- function(Z) cbind(sweep(sweep(Z, 2, m), 2, v, "/"),
                           sweep(sweep(sweep(Z, 2, m)^2, 2, 2 * v^2, "/"), 2, 1 / (2 * v), "-"))
  post_part <- -colSums(wp * dlg(X))
  fixed <- (post_part - colSums(wq * dlg(Y))) / NH
  sc <- tg$score(Y)
  path_prop <- cbind(sc, sweep(sc * E, 2, 2 * sqrt(v), "/") + matrix(1 / (2 * v), nrow(Y), length(v), byrow = TRUE))
  path <- (post_part + colSums(wq * path_prop)) / NH
  psi <- cbind(sweep(Xfit, 2, m), sweep(sweep(Xfit, 2, m)^2, 2, v))
  Sig <- stats::cov(psi) / nrow(Xfit)
  c(fixed = sqrt(drop(t(fixed) %*% Sig %*% fixed)) / d_obj$mcse,
    pathwise = sqrt(drop(t(path) %*% Sig %*% path)) / d_obj$mcse)
}
tab <- do.call(rbind, lapply(seq_len(nrow(cfg)), function(i) {
  dd <- cfg$d[i]; S <- cfg$S[i]
  tg <- t7_target(dd)
  invisible(jackknife_null_quantiles(G, probs = 1 - 0.05 / 4))
  msg("E4: d =", dd, "S =", S)
  res <- run_reps(n_reps(cfg$R[i], key = "REPS$e4"), function(r) {
    X <- tg$r(S); Xextra <- tg$r(S); E <- matrix(stats::rnorm(S * dd), S, dd)
    out <- list()
    a <- ratios(tg, X, prop_from(E, tg$m0, tg$v0), tg$m0, tg$v0)
    out$fixed <- bridge_diagnostics(a$ap, a$aq, iid = TRUE, eb = FALSE, jackknife = FALSE)
    me <- moments(Xextra); Ye <- prop_from(E, me$m, me$v)
    a <- ratios(tg, X, Ye, me$m, me$v)
    out$split <- bridge_diagnostics(a$ap, a$aq, iid = TRUE, eb = FALSE, jackknife = FALSE)
    ms <- moments(X); Ys <- prop_from(E, ms$m, ms$v)
    a <- ratios(tg, X, Ys, ms$m, ms$v)
    ds <- bridge_diagnostics(a$ap, a$aq, iid = TRUE, eb = FALSE, G = G)
    block <- bridgeoverlap:::.block_ids_prop(S, G)
    jv <- jack_refit(tg, X, E, Ys, block, G) / ds$mcse^2
    rg <- if (dd < 32) R_g(tg, X, Ys, E, ms$m, ms$v, ds, X) else c(fixed = NA, pathwise = NA)
    data.frame(D_fixed = out$fixed$logml - tg$logZ, D_same = ds$logml - tg$logZ, D_split = out$split$logml - tg$logZ,
               mcse_same = ds$mcse, mcse_split = out$split$mcse, Rg_fixed = rg[["fixed"]], Rg_path = rg[["pathwise"]],
               V_norefit = ds$V_jack, V_refit_fixed = jv[["fixed"]], V_refit_crn = jv[["crn"]])
  })
  bias <- function(D, s) mean(D) / stats::median(s)
  bse <- function(D, s) stats::sd(D) / sqrt(length(D)) / stats::median(s)
  data.frame(d = dd, S = S, R = nrow(res), sd_fixed = stats::sd(res$D_fixed), sd_same = stats::sd(res$D_same),
             sd_split = stats::sd(res$D_split), bias_same = bias(res$D_same, res$mcse_same),
             bias_same_se = bse(res$D_same, res$mcse_same), bias_split = bias(res$D_split, res$mcse_split),
             bias_split_se = bse(res$D_split, res$mcse_split), Rg_fixed = stats::median(res$Rg_fixed),
             Rg_pathwise = stats::median(res$Rg_path), V_jack_no_refit = mean(res$V_norefit),
             V_jack_refit_fixed = mean(res$V_refit_fixed), V_jack_refit_crn = mean(res$V_refit_crn))
}))
save_table(tab, "tab_E4")
print(tab, digits = 3)
bse <- function(b, s) ifelse(is.na(b), "--", sprintf("$%s\\,(%s)$", fsign(b, 2), fx(s, 2)))
emit_table("tab_E4", cbind(tab$d, tab$S, tab$R, fx(tab$sd_fixed, 4), fx(tab$sd_same, 4), fx(tab$sd_split, 4),
                           bse(tab$bias_same, tab$bias_same_se), bse(tab$bias_split, tab$bias_split_se),
                           fx(tab$Rg_fixed, 3), fx(tab$Rg_pathwise, 3), fx(tab$V_jack_no_refit, 3),
                           fx(tab$V_jack_refit_fixed, 3), fx(tab$V_jack_refit_crn, 3)))
d2 <- tab[tab$d == 2, ]
res_put("e4.Rg.fixed", frange(tab$Rg_fixed, 2))
for (S in d2$S) {
  res_put(sprintf("e4.Rg.path.%d", S), fx(d2$Rg_pathwise[d2$S == S], 3)); res_put(sprintf("e4.Rg.fixed.%d", S), fx(d2$Rg_fixed[d2$S == S], 2))
}
res_put("e4.vj.fixed.infl", frange(100 * (tab$V_jack_refit_fixed / tab$V_jack_no_refit - 1), 0))
for (S in unique(tab$S)) res_put(sprintf("e4.vj.crn.infl.%d", S), fsign(100 * mean(tab$V_jack_refit_crn[tab$S == S] / tab$V_jack_no_refit[tab$S == S] - 1), 0))
res_write("supplement/S11_E4_proposal_fitting.R")
