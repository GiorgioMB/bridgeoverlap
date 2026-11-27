#!/usr/bin/env Rscript
# Supplement Figure S-space from supplement/S20_problem_space.R (results/ps_laws.rds, ps_curve.rds,
# ps_refined.rds): (a) coverage against the effective overlap over all searched laws, with the worst
# case per bin; (b) the worst case over the laws of P(the run passes and its interval misses), against
# the threshold n0, with the calibration criterion; (c) coverage among the runs that pass at n0.
if (!exists("BO_ROOT")) {
  BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (length(f)) dirname(normalizePath(f)) else if (dir.exists("analysis")) "analysis" else "."})
  source(file.path(BO_ROOT, "R", "common.R")); source(file.path(BO_ROOT, "R", "paper.R"))
}
suppressPackageStartupMessages(library(ggplot2))
d <- readRDS(file.path(results_dir, "ps_laws.rds")); cu <- readRDS(file.path(results_dir, "ps_curve.rds"))
kl <- c(atoms = "1-3 atoms", edge = "atoms at extreme responsibilities", mixed = "mixed", beta = "smooth (Beta) laws", near = "near-perfect proposals")
d$kind <- factor(kl[d$kind], levels = kl)
pal <- setNames(c("#0072B2", "#D55E00", "#CC79A7", "#009E73", "#E69F00"), kl)
th <- theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), plot.title = element_text(size = 8.5, face = "bold"))
br <- c(1, 1.5, 2, 3, 4, 5, 7, 10, 15, 20, 30, 50, 80)
env <- do.call(rbind, lapply(seq_len(length(br) - 1), function(k) { m <- d$nov >= br[k] & d$nov < br[k + 1]
  if (any(m)) data.frame(nov = sqrt(br[k] * br[k + 1]), worst = min(d$cov[m])) }))
pa <- ggplot(d, aes(nov, cov)) + annotate("rect", xmin = 10, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey93") +
  geom_hline(yintercept = 0.95, linetype = 3) + geom_hline(yintercept = 0.94, linetype = 2, colour = "grey45") +
  geom_point(aes(colour = kind), size = 0.55, alpha = 0.55) +
  geom_step(data = env, aes(nov, worst), direction = "mid", colour = "black", linewidth = 0.55) +
  scale_x_log10(breaks = c(1, 3, 10, 30, 80)) + coord_cartesian(ylim = c(0.5, 1)) + scale_colour_manual(values = pal, name = NULL) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 1.8, alpha = 1))) + th +
  theme(legend.position = c(0.99, 0.02), legend.justification = c(1, 0), legend.key.size = unit(0.62, "lines"),
        legend.text = element_text(size = 6.8), legend.background = element_rect(fill = NA, colour = NA), legend.key = element_rect(fill = NA, colour = NA)) +
  labs(x = expression("effective overlap " * n[ov]), y = "coverage of nominal 95% intervals", title = "(a) Coverage against the clock")
pb <- ggplot(cu, aes(n0, worst_joint)) + geom_hline(yintercept = 0.06, linetype = 3) +
  geom_hline(aes(yintercept = worst_mcse_rule[1]), colour = "#009E73", linetype = 2) +
  annotate("text", x = 1, y = cu$worst_mcse_rule[1], vjust = -0.5, hjust = 0, size = 2.6, colour = "#009E73", label = "rule \"MCSE below 0.3\"") +
  geom_line() + geom_point(size = 1.6) + scale_x_log10(breaks = c(1, 3, 10, 30)) +
  coord_cartesian(ylim = c(0, max(0.14, 1.15 * cu$worst_mcse_rule[1]))) + th +
  labs(x = expression("threshold " * n[0]), y = "worst case of P(pass and miss)", title = "(b) Passing runs that miss, worst law")
cs <- cu$n0
cq <- do.call(rbind, lapply(cs, function(c0) { v <- d[[paste0("qcovpass", c0)]]; v <- v[is.finite(v)]
  if (length(v)) data.frame(c = c0, what = c("median over laws", "5% quantile", "worst law"), val = c(stats::median(v), stats::quantile(v, 0.05), min(v))) }))
cq$what <- factor(cq$what, levels = c("median over laws", "5% quantile", "worst law"))
pc <- ggplot(cq, aes(c, val, linetype = what)) + geom_hline(yintercept = 0.95, linetype = 3) + geom_line() + geom_point(size = 1.2) +
  scale_x_log10(breaks = c(1, 3, 10, 30)) + coord_cartesian(ylim = c(0.4, 1)) + th + guides(linetype = guide_legend(ncol = 1)) +
  theme(legend.position = c(0.98, 0.03), legend.justification = c(1, 0), legend.key.size = unit(0.7, "lines"),
        legend.text = element_text(size = 7), legend.background = element_rect(fill = alpha("white", 0.9), colour = NA)) +
  scale_linetype_manual(values = c(1, 2, 3), name = NULL) +
  labs(x = expression("threshold " * n[0]), y = "coverage among passing runs", title = "(c) Coverage of passing runs")
fig <- patchwork::wrap_plots(pa, pb, pc, nrow = 1, widths = c(1.35, 1, 1))
ggsave(fig_file("fig_problem_space"), fig, width = 7.4, height = 3.2)
msg("wrote", fig_file("fig_problem_space"))
