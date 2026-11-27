#!/usr/bin/env Rscript
# Supplement figure fig_collapse.pdf, drawn from the results saved by supplement/S23_collapse.R
# (results/collapse_curves.rds, results/collapse_multiscale.rds):
#   Rscript analysis/plot_collapse.R
# (a) Two-state chain in the small-overlap limit: coverage of nominal 95% intervals with the
# true long-run variance is nearly one curve in n_ovP for every spectral gap. (b) sqrt(n_ovP)
# times the Kolmogorov distance stays in a narrow band. (c) A weakly seen slow mode: at fixed
# n_ovP = 16 and tau_W = 3, coverage depends on how often the slow mode relaxes during the run.
if (!exists("BO_ROOT")) {
  BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (length(f)) dirname(normalizePath(f)) else if (dir.exists("analysis")) "analysis" else "."})
  source(file.path(BO_ROOT, "R", "common.R")); source(file.path(BO_ROOT, "R", "paper.R"))
}
if (!have_ggplot() || !requireNamespace("patchwork", quietly = TRUE)) stop("ggplot2 and patchwork are needed for the figure")
suppressPackageStartupMessages(library(ggplot2))
cu <- readRDS(file.path(results_dir, "collapse_curves.rds")); ms <- readRDS(file.path(results_dir, "collapse_multiscale.rds"))
lab <- function(g) ifelse(g == 0, "0 (limit)", format(g))
cu$g <- factor(lab(cu$gamma), levels = lab(sort(unique(cu$gamma), decreasing = TRUE)))
gcol <- setNames(rev(grDevices::hcl.colors(nlevels(cu$g) + 2, "Blues 3")[seq_len(nlevels(cu$g))]), levels(cu$g))   # slower = darker
th <- theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), legend.position = "bottom",
                                       plot.title = element_text(size = 9.5, face = "bold"))
xP <- list(scale_x_log10(breaks = c(0.5, 1, 2, 4, 8, 16, 32, 64), labels = c("0.5", "1", "2", "4", "8", "16", "32", "64")),
           labs(x = expression("effective overlap under MCMC  " * n[ov]^P)))
nominal <- geom_hline(yintercept = 0.95, linetype = 3, colour = "grey40")
pa <- ggplot(cu, aes(n_ovP, coverage, colour = g)) + nominal + geom_line(linewidth = 0.7) + geom_point(size = 0.9) + xP + th +
  scale_colour_manual(values = gcol, name = expression("spectral gap " * gamma)) +
  labs(y = "coverage of nominal 95% intervals", title = "(a) coverage, two-state chain")
pb <- ggplot(cu, aes(n_ovP, ks_scaled, colour = g)) + geom_line(linewidth = 0.7) + geom_point(size = 0.9) + xP + th +
  scale_colour_manual(values = gcol, name = expression("spectral gap " * gamma)) + expand_limits(y = 0) +
  labs(y = expression(sqrt(n[ov]^P) %*% "Kolmogorov distance"), title = "(b) normal approximation")
band <- range(cu$coverage[cu$n_ovP == 16])
pc <- ggplot(ms, aes(kappa, coverage)) + nominal +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = band[1], ymax = band[2], fill = "steelblue", alpha = 0.25) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.2) + th +
  scale_x_log10(breaks = c(0.25, 1, 4, 16), labels = c("0.25", "1", "4", "16")) +
  labs(x = expression("slow relaxations per run  " * kappa), y = "coverage of nominal 95% intervals",
       title = expression(bold("(c) a weakly seen slow mode, ") * bolditalic(n)[bold(ov)]^bold(P) * bold(" = 16")))
fig <- patchwork::wrap_plots(pa, pb, pc, nrow = 1, guides = "collect") & theme(legend.position = "bottom")
ggsave(fig_file("fig_collapse"), fig, width = 10, height = 3.6)
msg("wrote", fig_file("fig_collapse"))
