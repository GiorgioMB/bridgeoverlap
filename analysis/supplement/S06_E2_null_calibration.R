#!/usr/bin/env Rscript
# Supplement Table E2: null calibration with independent draws (5000 replications;
# rejection rates in % at nominal 5%). Setting 1: p = N(0, I_5) with the fixed proposal
# g = N(0.2 * 1, 1.3^2 I_5), S in {1000, 4000, 16000}, G in {10, 20}. Setting 2: Gaussian
# proposal fitted on an independent posterior sample of the same size (as in
# bridgesampling) for heavy-tailed, multimodal and high-dimensional targets, G = 20.
# The mixture is taken as 0.5 N(-3.5 e1, I_8) + 0.5 N(3.5 e1, I_8).
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(dirname(normalizePath(f))) else if (dir.exists("analysis")) "analysis" else ".."})
for (f in c("common.R", "families.R", "exact.R", "ar1.R", "supplement.R", "paper.R")) source(file.path(BO_ROOT, "R", f))

logZ <- -0.2
R <- n_reps(REPS$e2)
fixed_g <- function(d) list(m = rep(0.2, d), v = rep(1.3^2, d))
gauss5 <- function(S) {
  g <- fixed_g(5L)
  X <- matrix(stats::rnorm(S * 5), S, 5)
  Y <- sweep(matrix(stats::rnorm(S * 5), S, 5) * 1.3, 2, g$m, "+")
  list(ap = logZ + lmvn_std(X) - ldiag(X, g$m, g$v), aq = logZ + lmvn_std(Y) - ldiag(Y, g$m, g$v))
}
split_fit <- function(rpost, lpost) function(S) {
  fit <- fit_gauss(rpost(S)); X <- rpost(S); Y <- rgauss(S, fit)
  list(ap = logZ + lpost(X) - lgauss(X, fit), aq = logZ + lpost(Y) - lgauss(Y, fit))
}
rn <- function(d) function(n) matrix(stats::rnorm(n * d), n, d)
rt3 <- function(d) function(n) rmvt_std(n, d, 3)
mix8 <- list(r = function(n) { X <- matrix(stats::rnorm(n * 8), n, 8); X[, 1] <- X[, 1] + ifelse(stats::runif(n) < 0.5, -3.5, 3.5); X },
             l = function(X) { Y1 <- X; Y1[, 1] <- Y1[, 1] + 3.5; Y2 <- X; Y2[, 1] <- Y2[, 1] - 3.5
                               log(0.5) + lse2(lmvn_std(Y1), lmvn_std(Y2)) })
settings <- c(
  lapply(c(1000L, 4000L, 16000L), function(S) lapply(c(10L, 20L), function(G)
    list(name = "N(0,I_5), fixed g", S = S, G = G, draw = gauss5))),
  list(list(list(name = "N(0,I_16), split", S = 4000L, G = 20L, draw = split_fit(rn(16), lmvn_std)),
            list(name = "t_3(0,I_8), split", S = 4000L, G = 20L, draw = split_fit(rt3(8), function(X) lmvt_std(X, 3))),
            list(name = "t_3(0,I_32), split", S = 4000L, G = 20L, draw = split_fit(rt3(32), function(X) lmvt_std(X, 3))),
            list(name = "t_3(0,I_64), split", S = 1000L, G = 20L, draw = split_fit(rt3(64), function(X) lmvt_std(X, 3))),
            list(name = "mixture (+-3.5 e1), d = 8, split", S = 4000L, G = 20L, draw = split_fit(mix8$r, mix8$l)),
            list(name = "N(0,I_64), split", S = 300L, G = 20L, draw = split_fit(rn(64), lmvn_std)))))
settings <- unlist(settings, recursive = FALSE)
if (BO_QUICK) settings <- settings[c(2, 8)]
tab <- do.call(rbind, lapply(settings, function(st) {
  msg("E2:", st$name, "S =", st$S, "G =", st$G)
  invisible(jackknife_null_quantiles(st$G, probs = 0.95))
  res <- run_reps(R, function(r) {
    a <- st$draw(st$S)
    d <- bridge_diagnostics(a$ap, a$aq, iid = TRUE, G = st$G, eb = FALSE)
    rj <- rejections(d)
    data.frame(n_ov = d$n_ov, Delta = d$logml - logZ, mcse = d$mcse, V_jack = d$V_jack, V_jack_rej = rj[["V_jack"]],
               I_max_bonf = rj[["I_max_bonf"]], I_max_naive = rj[["I_max_naive"]], I_max_exact = rj[["I_max_exact"]],
               G_max = rj[["G_max"]], t_T = rj[["t_T"]])
  })
  data.frame(setting = st$name, S = st$S, G = st$G, R = nrow(res), med_n_ov = stats::median(res$n_ov),
             var_over_mcse2 = stats::var(res$Delta) / mean(res$mcse^2), mean_V_jack = mean(res$V_jack),
             V_jack = 100 * mean(res$V_jack_rej), I_max_bonf = 100 * mean(res$I_max_bonf),
             I_max_naive = 100 * mean(res$I_max_naive), I_max_exact = 100 * mean(res$I_max_exact),
             G_max = 100 * mean(res$G_max), t_T = 100 * mean(res$t_T),
             coverage = mean(abs(res$Delta) <= z975 * res$mcse))
}))
save_table(tab, "tab_E2")
print(tab, digits = 3)
lab <- c("N(0,I_5), fixed g" = "$\\mathcal N(0,I_5)$, fixed $g$", "N(0,I_16), split" = "$\\mathcal N(0,I_{16})$, split",
         "t_3(0,I_8), split" = "$t_3(0,I_{8})$, split", "t_3(0,I_32), split" = "$t_3(0,I_{32})$, split",
         "t_3(0,I_64), split" = "$t_3(0,I_{64})$, split", "mixture (+-3.5 e1), d = 8, split" = "mixture ($\\pm3.5e_1$), $d=8$, split",
         "N(0,I_64), split" = "$\\mathcal N(0,I_{64})$, split")
fixed <- tab$setting == "N(0,I_5), fixed g"
emit_table("tab_E2", cbind(lab[tab$setting], tab$S, tab$G, fint(tab$med_n_ov), fx(tab$var_over_mcse2, 3), fx(tab$mean_V_jack, 3),
                           fx(tab$V_jack, 1), fx(tab$I_max_bonf, 1), fx(tab$I_max_naive, 1), fx(tab$I_max_exact, 1),
                           fx(tab$G_max, 1), fx(tab$t_T, 1), fx(tab$coverage, 3)), groups = sum(fixed))
res_put("e2.R", fint(R))
res_put("e2.s1.bonf", frange(c(tab$I_max_bonf[fixed], tab$I_max_exact[fixed]), 1)); res_put("e2.s1.naive", frange(tab$I_max_naive[fixed], 1))
hd <- tab$setting %in% c("t_3(0,I_64), split", "N(0,I_64), split")
res_put("e2.hd", fx(max(c(tab$I_max_exact[hd], tab$G_max[hd])), 0)); res_put("e2.hd.naive", fx(max(tab$I_max_naive[hd]), 0))
res_put("e2.var", frange(tab$var_over_mcse2, 2))
res_put("e2.nov16", fint(tab$med_n_ov[tab$setting == "N(0,I_16), split"]))
res_write("supplement/S06_E2_null_calibration.R")
