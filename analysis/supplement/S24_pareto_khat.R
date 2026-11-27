#!/usr/bin/env Rscript
# Supplement Table S-khat: the Pareto k-hat diagnostic of importance sampling against the overlap
# check, for the families F1-F3 with independent draws. k-hat is computed as in Pareto smoothed
# importance sampling (Vehtari et al., 2024): a generalized Pareto distribution is fitted by the
# method of Zhang and Stephens (2009), with the weakly informative prior towards 0.5, to the
# largest min(0.2 S, 3 sqrt(S)) importance ratios, in both directions of the bridge: f/g at the
# proposal draws and g/f at the posterior draws. A run "passes" k-hat when both are at most 0.7.
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(dirname(normalizePath(f))) else if (dir.exists("analysis")) "analysis" else ".."})
for (f in c("common.R", "paper.R", "families.R")) source(file.path(BO_ROOT, "R", f))
z975 <- stats::qnorm(0.975)
gpd_khat <- function(x) {                          # x > 0, sorted increasingly; Zhang-Stephens estimate of k
  n <- length(x); m <- 30L + floor(sqrt(n)); prior <- 3
  xstar <- x[floor(n / 4 + 0.5)]
  theta <- 1 / x[n] + (1 - sqrt(m / (seq_len(m) - 0.5))) / prior / xstar
  lth <- n * vapply(theta, function(t) { k <- mean(log1p(-t * x)); log(-t / k) - k - 1 }, 0)
  w <- exp(lth - max(lth)); theta_hat <- sum(theta * w / sum(w))
  k <- mean(log1p(-theta_hat * x))
  k <- k * n / (n + 10) + 10 * 0.5 / (n + 10)      # weakly informative prior (Vehtari et al., 2024)
  if (is.finite(k)) k else Inf
}
pareto_khat <- function(log_w) {                   # k-hat of the importance ratios exp(log_w)
  S <- length(log_w); M <- ceiling(min(0.2 * S, 3 * sqrt(S)))
  lw <- sort(log_w - max(log_w)); x <- exp(lw[(S - M + 1):S]) - exp(lw[S - M])
  if (!any(x > 0)) return(-Inf)
  gpd_khat(x)
}
CONF <- list("F2 (mu_g=2)" = function() family_F2(2), "F2 (mu_g=4)" = function() family_F2(4),
             "F1 (d=10, sigma=1.5)" = function() family_F1(10, 1.5), "F1 (d=40, sigma=1.3)" = function() family_F1(40, 1.3),
             "F1 (d=100, sigma=1.4)" = function() family_F1(100, 1.4), "F3 (d=20)" = function() family_F3(20),
             "F3 (d=100)" = function() family_F3(100))
if (BO_QUICK) CONF <- CONF[c(1, 5, 6)]
SS <- if (BO_QUICK) c(100L, 1000L) else c(100L, 1000L, 10000L)
R <- n_reps(REPS$khat, key = "khat")
runs <- do.call(rbind, lapply(names(CONF), function(cf) do.call(rbind, lapply(SS, function(S) {
  fam <- CONF[[cf]](); msg("khat:", cf, "S =", S)
  out <- run_reps(R, function(r) {
    ap <- fam$lr(fam$rpost(S)); aq <- fam$lr(fam$rprop(S))            # log Z = 0
    d <- bridge_diagnostics(ap, aq, iid = TRUE, eb = FALSE, jackknife = FALSE)
    data.frame(n_ov = d$n_ov, pass = !isTRUE(d$flags[["overlap"]]),
               cover = is.finite(d$logml) && abs(d$logml) <= z975 * d$mcse,
               k_fg = pareto_khat(aq), k_gf = pareto_khat(-ap))
  }, cell = paste("khat", cf, S))
  out$config <- cf; out$S <- S; out
}))))
runs$k_max <- pmax(runs$k_fg, runs$k_gf); runs$kpass <- runs$k_max <= 0.7
save_table(runs, "khat_runs")
cov_of <- function(x) if (length(x)) mean(x) else NA_real_
tab <- do.call(rbind, lapply(split(runs, list(runs$config, runs$S), drop = TRUE), function(x) data.frame(
  config = x$config[1], S = x$S[1], k_fg = stats::median(x$k_fg), k_gf = stats::median(x$k_gf),
  kpass = 100 * mean(x$kpass), opass = 100 * mean(x$pass), cov = mean(x$cover),
  cov_kpass = cov_of(x$cover[x$kpass]), cov_opass = cov_of(x$cover[x$pass]))))
tab <- tab[order(match(tab$config, names(CONF)), tab$S), ]
save_table(tab, "tab_khat")
tex_label <- function(x) {                          # "F2 (mu_g=2)" -> "F2 ($\mu_g=2$)"
  x <- gsub("mu_g=([0-9.]+)", "$\\\\mu_g=\\1$", x); x <- gsub("sigma=([0-9.]+)", "$\\\\sigma=\\1$", x)
  gsub("d=([0-9]+)", "$d=\\1$", x)
}
rows <- sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s\\\\", ifelse(duplicated(tab$config), "", tex_label(tab$config)), fint(tab$S),
                fx(tab$k_fg, 2), fx(tab$k_gf, 2), fx(tab$kpass, 0), fx(tab$opass, 0), fx(tab$cov, 3), fx(tab$cov_kpass, 3), fx(tab$cov_opass, 3))
writeLines(c("\\begin{tabular}{@{}lrcccccccc@{}}\\toprule",
             " & & \\multicolumn{2}{c}{median $\\hat k$} & \\multicolumn{2}{c}{runs that pass (\\%)} & \\multicolumn{3}{c}{coverage of nominal $95\\%$ intervals}\\\\",
             "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(l){7-9}",
             "Configuration & $S$ & $f/g$ & $g/f$ & $\\hat k\\le0.7$ & $\\neffhat\\ge10$ & all & $\\hat k\\le0.7$ & $\\neffhat\\ge10$\\\\\\midrule",
             rows, "\\bottomrule\\end{tabular}"), file.path(gen_dir, "tab_khat.tex"))
res_put("khat.R", fint(R))
res_put("khat.kflag.all", fx(100 * mean(!runs$kpass), 0))                    # share of all runs that k-hat flags
res_put("khat.kflag.opass", fx(100 * mean(!runs$kpass[runs$pass]), 0))       # ... among the runs that pass the overlap check
res_put("khat.cov.kflag.opass", fx(cov_of(runs$cover[!runs$kpass & runs$pass]), 3))
f3 <- grepl("^F3", runs$config)
res_put("khat.F3.k", fx(stats::median(runs$k_max[f3]), 2)); res_put("khat.F3.cov", fx(cov_of(runs$cover[f3 & runs$pass]), 3))
res_write("supplement/S24_pareto_khat.R")
