#!/usr/bin/env Rscript
# Supplement figure fig_mcmc.pdf, drawn from the results saved by 02_markov_chains.R
# (results/fig2_two_state_chain.rds, results/fig2_ar1_chains.rds):
#   Rscript analysis/plot_mcmc.R
# Story. Top, two-state chain of Theorem 7 computed exactly: against n_ovP, coverage with the
# true long-run variance and the Kolmogorov distance are one curve for every spectral gap,
# while intervals that ignore autocorrelation plateau below nominal however long the run.
# Bottom, AR(1) chains: an estimated long-run variance costs coverage at small n_ovP because
# the estimate is biased downward, more so for slower chains; both recover by n_ovP ~ 20-30.
# Darker colours mean slower mixing in every panel.
if (!exists("BO_ROOT")) {
  BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (length(f)) dirname(normalizePath(f)) else if (dir.exists("analysis")) "analysis" else "."})
  source(file.path(BO_ROOT, "R", "common.R")); source(file.path(BO_ROOT, "R", "paper.R"))
}
if (!have_ggplot() || !requireNamespace("patchwork", quietly = TRUE)) stop("ggplot2 and patchwork are needed for the figure")
suppressPackageStartupMessages(library(ggplot2))
tw <- readRDS(file.path(results_dir, "fig2_two_state_chain.rds")); tw <- tw[!is.na(tw$gamma), ]
ar <- readRDS(file.path(results_dir, "fig2_ar1_chains.rds"))

pal <- function(n) grDevices::hcl.colors(n + 2, "Blues 3")[seq_len(n)]           # dark to light, no near-white
th <- theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), legend.position = "bottom",
                                       plot.title = element_text(size = 9.5, face = "bold"))
xP <- list(scale_x_log10(breaks = c(1, 3, 10, 30, 100), labels = c("1", "3", "10", "30", "100")),
           labs(x = expression("effective overlap under MCMC  " * n[ov]^P)))
nominal <- geom_hline(yintercept = 0.95, linetype = 3, colour = "grey40")

# top: two-state chain, exact
tw$g <- factor(format(tw$gamma), levels = format(sort(unique(tw$gamma), decreasing = FALSE)))
gcol <- setNames(pal(nlevels(tw$g)), levels(tw$g))                            # gamma = 0.05 (slowest) darkest
long <- rbind(data.frame(tw[, c("g", "n_ovP")], cov = tw$coverage_sigmaNP, mcse = "with the long-run variance"),
              data.frame(tw[, c("g", "n_ovP")], cov = tw$coverage_sigmaN, mcse = "ignoring autocorrelation"))
long$mcse <- factor(long$mcse, levels = c("with the long-run variance", "ignoring autocorrelation"))
ends <- do.call(rbind, lapply(split(tw, tw$g), function(x) x[which.max(x$n_ovP), ]))
ends <- ends[ends$gamma < 1, ]
pa <- ggplot(long, aes(n_ovP, cov, colour = g, linetype = mcse)) + nominal + geom_line(linewidth = 0.7) +
  annotate("text", x = ends$n_ovP * 1.08, y = ends$coverage_sigmaN, label = sprintf("gamma == %s", format(ends$gamma)),
           parse = TRUE, hjust = 0, size = 2.7, colour = gcol[as.character(ends$g)]) +
  scale_colour_manual(values = gcol, name = expression("spectral gap " * gamma)) +
  scale_linetype_manual(values = c("solid", "22"), name = "intervals") +
  scale_y_continuous(limits = c(0.25, 1.01), breaks = seq(0.3, 0.9, 0.2)) + xP + th +
  coord_cartesian(xlim = c(0.45, 260)) + labs(y = "coverage of nominal 95% intervals") +
  annotate("text", x = 0.55, y = 0.99, label = '"true long-run variance (solid): one curve for every"~gamma', parse = TRUE, hjust = 0, size = 2.7, colour = "grey20") +
  annotate("text", x = 2.1, y = 0.705, label = "ignoring autocorrelation (dashed)", hjust = 0, size = 2.7, colour = "grey20") +
  guides(colour = "none", linetype = "none") +
  ggtitle(expression(bold("Two-state chain: coverage is one curve in " * n[ov]^P)))
tri <- data.frame(x = c(8, 80, 80, 8), y = c(0.2, 0.2, 0.2 / sqrt(10), 0.2))
pb <- ggplot(tw, aes(n_ovP, ks, colour = g)) + geom_line(linewidth = 0.7) +
  geom_path(data = tri, aes(x, y), inherit.aes = FALSE, colour = "grey30", linewidth = 0.4) +
  annotate("text", x = 26, y = 0.085, label = "slope~-1/2", parse = TRUE, size = 2.7, colour = "grey20", hjust = 0.5) +
  scale_colour_manual(values = gcol, name = expression("spectral gap " * gamma)) +
  scale_y_log10(breaks = c(0.01, 0.03, 0.1, 0.3), labels = c("0.01", "0.03", "0.1", "0.3")) + xP + th +
  theme(legend.position = c(0.03, 0.03), legend.justification = c(0, 0), legend.background = element_rect(fill = alpha("white", 0.85))) +
  labs(y = "Kolmogorov distance to N(0, 1)") +
  ggtitle("... and so is the normal approximation")

# bottom: AR(1) chains, mean over the three families with their range
agg <- function(v) do.call(rbind, lapply(split(ar, list(ar$rho, round(ar$n_ovP, 0)), drop = TRUE), function(x)
  data.frame(rho = x$rho[1], n_ovP = mean(x$n_ovP), m = mean(x[[v]]), lo = min(x[[v]]), hi = max(x[[v]]))))
cs <- agg("coverage_spectral"); ct <- agg("coverage_true_sigmaNP"); vr <- agg("med_nP_over_vP")
ct <- do.call(rbind, lapply(split(ct, round(ct$n_ovP)), function(x) data.frame(n_ovP = mean(x$n_ovP), m = mean(x$m), lo = min(x$lo), hi = max(x$hi))))
cs$r <- factor(format(cs$rho), levels = format(sort(unique(ar$rho)))); vr$r <- factor(format(vr$rho), levels = levels(cs$r))
rcol <- setNames(rev(pal(nlevels(cs$r))), levels(cs$r))                     # rho = 0.95 (slowest) darkest
pc <- ggplot(cs, aes(n_ovP, m, colour = r, fill = r)) + nominal +
  geom_ribbon(data = ct, aes(n_ovP, ymin = lo, ymax = hi), inherit.aes = FALSE, fill = "grey80", alpha = 0.6) +
  geom_line(data = ct, aes(n_ovP, m), inherit.aes = FALSE, colour = "grey45", linewidth = 0.7) +
  annotate("text", x = 1.05, y = 0.962, label = "true long-run variance", hjust = 0, vjust = 0, size = 2.7, colour = "grey35") +
  geom_ribbon(aes(ymin = lo, ymax = hi), colour = NA, alpha = 0.18) + geom_line(linewidth = 0.7) + geom_point(size = 1.3) +
  scale_colour_manual(values = rcol, name = expression("AR(1) coefficient " * rho)) + scale_fill_manual(values = rcol, guide = "none") +
  scale_y_continuous(limits = c(0.76, 0.98), breaks = seq(0.8, 0.95, 0.05)) + xP + th +
  theme(legend.position = c(0.98, 0.03), legend.justification = c(1, 0), legend.background = element_rect(fill = alpha("white", 0.85))) +
  labs(y = "coverage of nominal 95% intervals") +
  ggtitle("AR(1) chains: estimating the long-run variance costs coverage")
pd <- ggplot(vr, aes(n_ovP, m, colour = r, fill = r)) + geom_hline(yintercept = 1, linetype = 3, colour = "grey40") +
  annotate("text", x = 1.05, y = 1.035, label = "unbiased", hjust = 0, size = 2.7, colour = "grey35") +
  geom_ribbon(aes(ymin = lo, ymax = hi), colour = NA, alpha = 0.18) + geom_line(linewidth = 0.7) + geom_point(size = 1.3) +
  scale_colour_manual(values = rcol, guide = "none") + scale_fill_manual(values = rcol, guide = "none") +
  scale_y_continuous(limits = c(0.25, 1.08), breaks = seq(0.3, 1, 0.1)) + xP + th +
  labs(y = expression("median of " * hat(n)[P] / v[P] * ", estimated / true")) +
  ggtitle("... because the estimate is biased downward")

fig <- patchwork::wrap_plots(pa, pb, pc, pd, ncol = 2) +
  patchwork::plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
ggsave(fig_file("fig_mcmc"), fig, width = 10, height = 7.2)
msg("wrote generated/fig_mcmc.pdf")