#!/usr/bin/env Rscript
# Figure of Section 7.4 (fig_bench.pdf) from the runs saved by 31_stan_benchmark_summary.R:
#   (a) every Stan run: error in MCSE units against the certificate;
#   (b) every pair of runs of the comparison where Monte Carlo error decides the Kass-Raftery
#       category: error of log B against the smaller certificate of the pair.
# Writes the macros fb.* quoted with the figure.
if (!exists("BO_ROOT")) {
  BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (length(f)) dirname(normalizePath(f)) else if (dir.exists("analysis")) "analysis" else "."})
  source(file.path(BO_ROOT, "R", "common.R")); source(file.path(BO_ROOT, "R", "paper.R"))
}
if (!have_ggplot() || !requireNamespace("patchwork", quietly = TRUE)) stop("ggplot2 and patchwork are needed for the figures")
suppressPackageStartupMessages(library(ggplot2))
res_get <- function(k) base::get(k, envir = .res); plain <- function(s) gsub("{,}", ",", s, fixed = TRUE)
runs <- readRDS(file.path(results_dir, "bench_runs_all.rds"))
c_mcmc <- 30
runs$cert <- runs$n_ovP                                          # certified iff cert >= c (plain rule, MCMC)
runs$zerr <- abs(runs$error) / runs$mcse
files <- list.files(results_dir, pattern = "^bench_runs_.*[.]rds$", full.names = TRUE)
files <- files[basename(files) != "bench_runs_all.rds"]
lab1 <- do.call(rbind, lapply(files, function(f) { x <- readRDS(f)
  if ("chains_label1" %in% names(x)) data.frame(config = x$config, rep = x$rep, chains_label1 = x$chains_label1) }))
runs <- merge(runs, lab1, by = c("config", "rep"), all.x = TRUE)
runs$grp <- ifelse(runs$f_std %in% TRUE, "fails R-hat, ESS or divergence checks",
                   ifelse(runs$config == "mixture" & runs$chains_label1 %in% c(0, 4), "passes them; label switching (blind spot)",
                          "passes R-hat, ESS and divergence checks"))
runs$grp <- factor(runs$grp, levels = c("passes R-hat, ESS and divergence checks", "passes them; label switching (blind spot)", "fails R-hat, ESS or divergence checks"))
ok <- runs$grp == levels(runs$grp)[1]; cert <- runs$cert >= c_mcmc
cv <- function(m) fx(mean(runs$cover[m]), 3)
res_put("fb.c", c_mcmc); res_put("fb.cert.cov", cv(ok & cert)); res_put("fb.cert.n", fint(sum(ok & cert)))
res_put("fb.cert.big", fint(sum(ok & cert & runs$zerr > 4))); res_put("fb.unc.cov", cv(ok & !cert)); res_put("fb.unc.n", fint(sum(ok & !cert)))
res_put("fb.unc.big", fx(100 * mean(runs$zerr[ok & !cert] > 4), 0))
res_put("fb.fail.cov", cv(runs$f_std %in% TRUE)); res_put("fb.fail.cert.cov", cv((runs$f_std %in% TRUE) & cert))
res_put("fb.blind.n", fint(sum(runs$grp == levels(runs$grp)[2])))
th <- theme_bw(base_size = 9.5) + theme(panel.grid.minor = element_blank(), legend.position = "bottom", legend.direction = "vertical",
                                       legend.text = element_text(size = 8.5, hjust = 0), legend.key.size = unit(0.8, "lines"),
                                       legend.margin = margin(0, 0, 0, 0), plot.title = element_text(size = 9.5, face = "bold"),
                                       plot.subtitle = element_text(size = 8.5, colour = "grey20"))
xlo <- 0.3; ylo <- 0.03; yhi <- 300
runs$y <- pmin(pmax(runs$zerr, ylo), yhi); runs$cx <- pmax(runs$cert, xlo)   # out-of-range points drawn on the edges
runs$cls <- factor(ifelse(runs$grp == levels(runs$grp)[3], "fails", ifelse(runs$grp == levels(runs$grp)[2], "blind",
                   ifelse(runs$cert >= c_mcmc, "cert", "unc"))), levels = c("cert", "unc", "blind", "fails"))
cls_lab <- c(cert = "passes the sampler checks; certified", unc = "passes the sampler checks; not certified",
             blind = "label switching: passes every check, off by log 2", fails = "fails the sampler checks")
cls_col <- c(cert = "#0072B2", unc = "#E69F00", blind = "#D55E00", fails = "grey70")
cls_shp <- c(cert = 16, unc = 16, blind = 4, fails = 1)
pa <- ggplot(runs, aes(cx, y, colour = cls, shape = cls)) +
  annotate("rect", xmin = c_mcmc, xmax = Inf, ymin = 0, ymax = Inf, fill = "#EAF2F8") +
  geom_point(data = runs[runs$cls == "fails", ], size = 0.8, alpha = 0.6) +
  geom_point(data = runs[runs$cls %in% c("cert", "unc"), ], size = 1, alpha = 0.55) +
  geom_point(data = runs[runs$cls == "blind", ], size = 2.2, stroke = 1) +
  geom_hline(yintercept = stats::qnorm(0.975), linetype = 2) + geom_vline(xintercept = c_mcmc) +
  annotate("text", x = c_mcmc * 1.1, y = yhi, label = "certified", hjust = 0, vjust = 1, size = 3.2, fontface = "bold", colour = "#0072B2") +
  scale_x_log10(breaks = c(0.3, 1, 3, 10, 30, 100, 1000), labels = c("< 0.3", "1", "3", "10", "30", "100", "1000"), limits = c(xlo, 3000)) +
  scale_y_log10(breaks = c(0.1, 1, 1.96, 10, 100, 300), labels = c("0.1", "1", "1.96", "10", "100", "> 300"), limits = c(ylo, yhi)) +
  scale_colour_manual(values = cls_col, labels = cls_lab, name = NULL, drop = FALSE) +
  scale_shape_manual(values = cls_shp, labels = cls_lab, name = NULL, drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(size = 2.4, alpha = 1))) + th +
  labs(x = expression("certificate  "*hat(n)[ov]^P), y = "|error| / MCSE", title = "(a) Every run",
       subtitle = sprintf("Runs that pass the sampler checks:\ncertified: %s runs, coverage %s\nnot certified: %s runs, coverage %s",
                          plain(res_get("fb.cert.n")), res_get("fb.cert.cov"), plain(res_get("fb.unc.n")), res_get("fb.unc.cov")))
# (b) Bayes factors near a category boundary
pairs_of <- function(A, B, nm) {
  a <- runs[runs$config == A, ]; b <- runs[runs$config == B, ]; r <- intersect(a$rep, b$rep)
  a <- a[match(r, a$rep), ]; b <- b[match(r, b$rep), ]; truth <- a$truth[1] - b$truth[1]
  est <- a$logml - b$logml; cuts <- c(-Inf, -5, -3, -1, 1, 3, 5, Inf)
  data.frame(cmp = nm, pc = pmin(a$cert, b$cert), err = est - truth, truth = truth,
             wrong = findInterval(est, cuts) != findInterval(truth, cuts))
}
bp <- rbind(pairs_of("radon_mn_vis_c", "radon_mn_vi_c", "centred models"), pairs_of("radon_mn_vis_nc", "radon_mn_vi_nc", "non-centred models"))
truth <- bp$truth[1]; cut_lo <- -1 - truth; cut_hi <- 1 - truth
res_put("fb.bf.wrong", fint(sum(bp$wrong))); res_put("fb.bf.wrong.unc", fint(sum(bp$wrong & bp$pc < c_mcmc)))
res_put("fb.bf.cert", fint(sum(bp$pc >= c_mcmc))); res_put("fb.bf.cert.right", fint(sum(!bp$wrong & bp$pc >= c_mcmc)))
res_put("fb.bf.wrong.maxcert", fx(max(bp$pc[bp$wrong]), 0)); res_put("fb.bf.truth", fx(truth, 2))
pb <- ggplot(bp, aes(pc, err)) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = cut_lo, fill = "#FBE3D6") +
  annotate("rect", xmin = c_mcmc, xmax = Inf, ymin = cut_lo, ymax = Inf, fill = "#EAF2F8") +
  geom_hline(yintercept = 0, colour = "grey40") + geom_hline(yintercept = cut_lo, linetype = 2) + geom_vline(xintercept = c_mcmc) +
  geom_point(aes(colour = cmp, shape = cmp), size = 1.6, alpha = 0.8) +
  annotate("text", x = 1.04, y = min(bp$err) - 0.02, hjust = 0, vjust = 0, size = 3.1, colour = "#9A3A12", label = "wrong Kass-Raftery\ncategory") +
  annotate("text", x = c_mcmc * 1.1, y = Inf, label = "certified", hjust = 0, vjust = 1.4, size = 3.2, fontface = "bold", colour = "#0072B2") +
  scale_x_log10(breaks = c(1, 3, 10, 30, 100), labels = c("1", "3", "10", "30", "100"), limits = c(1, 60)) +
  scale_colour_manual(values = c("#D55E00", "#0072B2"), name = NULL) + scale_shape_manual(values = c(17, 16), name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.4, alpha = 1))) +
  coord_cartesian(ylim = c(min(bp$err) - 0.05, max(0.7, max(bp$err)))) + th +
  labs(x = "smaller certificate of the two runs", y = "estimated minus exact log Bayes factor",
       title = "(b) Bayes factors, slopes against intercepts",
       subtitle = sprintf("exact log B = %s; %s wrong conclusions, all from\npairs with a certificate below 30 (at most %s);\nall %s certified pairs right",
                          res_get("fb.bf.truth"), res_get("fb.bf.wrong"), res_get("fb.bf.wrong.maxcert"), res_get("fb.bf.cert")))
fig <- patchwork::wrap_plots(pa, pb, nrow = 1, widths = c(1.15, 1))
ggsave(fig_file("fig_bench"), fig, width = 7.4, height = 5.1)
msg("wrote generated/fig_bench.pdf")
res_write("plot_bench.R")
