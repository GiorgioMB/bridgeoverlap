#!/usr/bin/env Rscript
# Figure S1 from the search of supplement/S20_problem_space.R (results/ps_laws.rds).
if (!exists("BO_ROOT")) {
  BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (length(f)) dirname(normalizePath(f)) else if (dir.exists("analysis")) "analysis" else "."})
  source(file.path(BO_ROOT, "R", "common.R")); source(file.path(BO_ROOT, "R", "paper.R"))
}
suppressPackageStartupMessages(library(ggplot2))
d <- readRDS(file.path(results_dir, "ps_laws.rds"))
kl <- c(atoms = "1-3 atoms", edge = "atoms at extreme responsibilities", mixed = "mixed", beta = "smooth (Beta) laws", near = "near-perfect proposals")
d$kind <- factor(kl[d$kind], levels = kl)
pal <- c("#0072B2", "#D55E00", "#CC79A7", "#009E73", "#E69F00"); names(pal) <- kl
th <- theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), legend.position = "bottom", plot.title = element_text(size = 8.5, face = "bold"))
br <- c(1, 1.5, 2, 3, 4, 5, 7, 10, 15, 20, 30, 50, 80)
env <- do.call(rbind, lapply(seq_len(length(br) - 1), function(k) { m <- d$nov >= br[k] & d$nov < br[k + 1]
  if (any(m)) data.frame(nov = sqrt(br[k] * br[k + 1]), worst = min(d$cov[m])) }))
pa <- ggplot(d, aes(nov, cov)) +
  annotate("rect", xmin = 10, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey93") +
  geom_hline(yintercept = 0.95, linetype = 3) + geom_point(aes(colour = kind), size = 0.55, alpha = 0.55) +
  geom_step(data = env, aes(nov, worst), direction = "mid", colour = "black", linewidth = 0.55) +
  scale_x_log10(breaks = c(1, 3, 10, 30, 80)) + coord_cartesian(ylim = c(0.5, 1)) + scale_colour_manual(values = pal, name = NULL) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 1.8, alpha = 1))) + th +
  theme(legend.position = c(0.99, 0.02), legend.justification = c(1, 0), legend.key.size = unit(0.62, "lines"),
        legend.text = element_text(size = 6.8), legend.background = element_rect(fill = NA, colour = NA), legend.key = element_rect(fill = NA, colour = NA)) +
  labs(x = expression("effective overlap "*n[ov]), y = "coverage of nominal 95% intervals", title = "(a) Coverage against the clock")
cs <- c(1, 2, 3, 5, 7, 10, 15, 20, 30)
jt <- data.frame(c = cs, val = vapply(cs, function(c0) max(d[[paste0("qjoint", c0)]]), numeric(1)))
pb <- ggplot(jt, aes(c, val)) + geom_hline(yintercept = 0.05, linetype = 3) +
  geom_hline(yintercept = max(d$mjoint03), colour = "#009E73", linetype = 2) +
  annotate("text", x = 1, y = max(d$mjoint03), vjust = -0.5, hjust = 0, size = 2.6, colour = "#009E73", label = "rule \"MCSE below 0.3\"") +
  geom_line() + geom_point(size = 1.6) +
  annotate("text", x = 2.2, y = jt$val[jt$c == 3] + 0.03, hjust = 0, size = 2.6, label = "Algorithm 1, certified at c") +
  scale_x_log10(breaks = c(1, 3, 10, 30)) +
  coord_cartesian(ylim = c(0, max(0.14, 1.15 * max(d$mjoint03)))) + th +
  labs(x = "threshold c", y = "worst case of P(miss and certified)", title = "(b) Certified misses, worst law")
cq <- do.call(rbind, lapply(c(1, 2, 3, 5, 7, 10, 15, 20, 30), function(c0) { v <- d[[paste0("covpass", c0)]]; v <- v[is.finite(v)]
  data.frame(c = c0, what = c("median over laws", "5% quantile", "worst law"), val = c(stats::median(v), stats::quantile(v, 0.05), min(v))) }))
cq$what <- factor(cq$what, levels = c("median over laws", "5% quantile", "worst law"))
pc <- ggplot(cq, aes(c, val, linetype = what)) + geom_hline(yintercept = 0.95, linetype = 3) + geom_line() + geom_point(size = 1.2) +
  scale_x_log10(breaks = c(1, 3, 10, 30)) + coord_cartesian(ylim = c(0.4, 1)) + th + guides(linetype = guide_legend(ncol = 1)) +
  theme(legend.position = c(0.98, 0.03), legend.justification = c(1, 0), legend.key.size = unit(0.7, "lines"),
        legend.text = element_text(size = 7), legend.background = element_rect(fill = alpha("white", 0.9), colour = NA)) +
  scale_linetype_manual(values = c(1, 2, 3), name = NULL) +
  labs(x = "threshold c", y = "coverage among certified runs", title = "(c) Coverage of certified runs")
fig <- patchwork::wrap_plots(pa, pb, pc, nrow = 1, widths = c(1.35, 1, 1))
ggsave(fig_file("fig_problem_space"), fig, width = 7.4, height = 3.2)
msg("wrote generated/fig_problem_space.pdf")
