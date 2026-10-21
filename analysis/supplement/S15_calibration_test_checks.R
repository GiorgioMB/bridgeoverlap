#!/usr/bin/env Rscript
# Supplement Section S5.6 (text): checks of Proposition S1 (the calibration test under
# sampler failure) and of the blind spot.
# (b) p = 0.5 N(-1.5, 0.4^2) + 0.5 N(1.5, 0.8^2), g = N(0, 1.6^2), and a sampler that puts
#     weight 1/2 + eta on the narrow component: by quadrature, (c~ - c*)/(-eta E_{m*} xi)
#     and T_inf/(eta Cov_{m*}(pi*, xi)) for eta in {0.005, 0.02, 0.05, 0.1}; by simulation
#     (S1 = S2 = 2000, 1000 replications) the bias in MCSEs and the rejection rate of t_T for
#     eta in {0, 0.01, 0.02, 0.04}, against the noncentral t_19 prediction.
# Blind spot: p = 0.7 N(0,1) + 0.3 N(10,1), g = N(0, 1.5^2), sampler draws from N(0,1).
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(dirname(normalizePath(f))) else if (dir.exists("analysis")) "analysis" else ".."})
for (f in c("common.R", "families.R", "exact.R", "ar1.R", "supplement.R", "paper.R")) source(file.path(BO_ROOT, "R", f))

S <- 2000L
G <- 20L
lmix <- function(w, m1, s1, m2, s2) function(x) lse2(log(w) + stats::dnorm(x, m1, s1, log = TRUE),
                                                     log(1 - w) + stats::dnorm(x, m2, s2, log = TRUE))
rmix <- function(w, m1, s1, m2, s2) function(n) { k <- stats::runif(n) < w; ifelse(k, stats::rnorm(n, m1, s1), stats::rnorm(n, m2, s2)) }
breaks <- c(-30, -10, -5, -3, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3, 5, 8, 10, 12, 30)

# Population root c~ and limit T_inf when the posterior draws come from p~ (log density lpt).
alt_limits <- function(lp, lg, lpt, s1 = 0.5) {
  s2 <- 1 - s1; lr0 <- log(s1 / s2)
  pic <- function(c, x) stats::plogis(c + lg(x) - lp(x) - lr0)
  E <- function(ld, h) integrate_pieces(function(x) h(x) * exp(ld(x)), breaks, -Inf)
  Ft <- function(c) s2 * E(lg, function(x) 1 - pic(c, x)) - s1 * E(lpt, function(x) pic(c, x))
  ct <- stats::uniroot(Ft, c(-8, 8), tol = 1e-12)$root
  H0 <- s1 * E(lpt, function(x) pic(ct, x))
  T_inf <- (s1 * E(lpt, function(x) pic(ct, x)^2) - s2 * E(lg, function(x) pic(ct, x) * (1 - pic(ct, x)))) / H0
  c(c_tilde = ct, T_inf = T_inf)
}

## (b) local alternatives -----------------------------------------------------------------
lp <- lmix(0.5, -1.5, 0.4, 1.5, 0.8)
lg <- function(x) stats::dnorm(x, 0, 1.6, log = TRUE)
xi <- function(x) exp(stats::dnorm(x, -1.5, 0.4, log = TRUE) - lp(x)) - exp(stats::dnorm(x, 1.5, 0.8, log = TRUE) - lp(x))  # bounded in [-2, 2]
ov <- overlap_1d(lp, lg, S, S, breaks)
ms <- mstar_moments(ov, xi)
quad <- do.call(rbind, lapply(c(0.005, 0.02, 0.05, 0.1), function(eta) {
  lim <- alt_limits(lp, lg, lmix(0.5 + eta, -1.5, 0.4, 1.5, 0.8))
  data.frame(eta = eta, ratio_bias = lim[["c_tilde"]] / (-eta * ms$E_xi),
             ratio_T = lim[["T_inf"]] / (eta * ms$cov_xi))
}))
save_table(quad, "checks_calib_alt_quadrature")
print(quad, digits = 4)
invisible(jackknife_null_quantiles(G, probs = 0.95))
crit <- stats::qt(0.975, G - 1)
sims <- do.call(rbind, lapply(quick_pick(c(0, 0.01, 0.02, 0.04)), function(eta) {
  msg("calibration test under mode imbalance, eta =", eta)
  rp <- rmix(0.5 + eta, -1.5, 0.4, 1.5, 0.8)
  res <- run_reps(n_reps(REPS$ca), function(r) {
    x <- rp(S); y <- stats::rnorm(S, 0, 1.6)
    d <- bridge_diagnostics(lp(x) - lg(x), lp(y) - lg(y), iid = TRUE, eb = FALSE, G = G)
    data.frame(D = d$logml, mcse = d$mcse, t_T = d$t_T)
  })
  ncp <- eta * ms$cov_xi * sqrt(ov$NH) / sqrt(ms$var)
  data.frame(eta = eta, R = nrow(res), bias_in_mcse = mean(res$D) / stats::median(res$mcse),
             predicted_bias_in_mcse = -eta * ms$E_xi / ov$sigma_N,
             rejection = mean(abs(res$t_T) > crit, na.rm = TRUE),
             predicted_rejection = 1 - stats::pt(crit, G - 1, ncp) + stats::pt(-crit, G - 1, ncp))
}))
save_table(sims, "checks_calib_alt_simulation")
print(sims, digits = 3)

## blind spot ---------------------------------------------------------------------------------
lpb <- lmix(0.7, 0, 1, 10, 1)
lgb <- function(x) stats::dnorm(x, 0, 1.5, log = TRUE)
lim <- alt_limits(lpb, lgb, function(x) stats::dnorm(x, log = TRUE))
res <- run_reps(n_reps(REPS$ca), function(r) {
  x <- stats::rnorm(S); y <- stats::rnorm(S, 0, 1.5)
  d <- bridge_diagnostics(lpb(x) - lgb(x), lpb(y) - lgb(y), iid = TRUE, eb = FALSE, G = G)
  data.frame(D = d$logml, mcse = d$mcse, t_T = d$t_T)
})
blind <- data.frame(c_tilde_minus_c_star = lim[["c_tilde"]], log_0.7 = log(0.7), T_inf = lim[["T_inf"]],
                    R = nrow(res), bias_in_mcse = mean(res$D) / stats::median(res$mcse),
                    coverage = mean(abs(res$D) <= z975 * res$mcse),
                    rejection_t_T = mean(abs(res$t_T) > crit, na.rm = TRUE))
save_table(blind, "checks_blind_spot")
print(blind, digits = 5)

## LaTeX macros ---------------------------------------------------------------------------------
and_list <- function(v) if (length(v) == 1) v else paste0(paste(v[-length(v)], collapse = ", "), " and ", v[length(v)])
res_put("ca.ratio.bias", and_list(fx(quad$ratio_bias, 3))); res_put("ca.ratio.T", and_list(fx(quad$ratio_T, 2)))
alt <- sims[sims$eta > 0, ]
res_put("ca.eta", and_list(format(alt$eta))); res_put("ca.sim.bias", and_list(fx(alt$bias_in_mcse, 2)))
res_put("ca.sim.rej", and_list(paste0(fx(100 * alt$rejection, 1), "\\%")))
res_put("ca.pred.rej", and_list(paste0(fx(100 * alt$predicted_rejection, 1), "\\%")))
res_put("ca.sim.null", paste0(fx(100 * sims$rejection[sims$eta == 0], 1), "\\%")); res_put("ca.R", fint(n_reps(REPS$ca)))
res_put("bs.ctilde", fx(blind$c_tilde_minus_c_star, 5)); res_put("bs.log07", fx(log(0.7), 5))
res_put("bs.Tinf", fsci(blind$T_inf, 1)); res_put("bs.bias", fx(abs(blind$bias_in_mcse), 0))
res_put("bs.cov", fx(blind$coverage, 3)); res_put("bs.rej", paste0(fx(100 * blind$rejection_t_T, 1), "\\%")); res_put("bs.R", fint(blind$R))
res_write("supplement/S15_calibration_test_checks.R")
