#!/usr/bin/env Rscript
# Figure 1 (fig_universality.pdf) and the two-point figure of the supplement (fig_two_point.pdf),
# drawn from the results saved by 01_universality.R (results/tab_E6.rds and
# results/fig1_universality_two_point.rds), so that they can be redrawn in seconds without
# rerunning the simulations:   Rscript analysis/plot_universality.R
#
# Figure 1 makes one point: the effective overlap, not the number of draws, is the clock.
# Coverage (top) and the median relative error of the MCSE (bottom) are plotted against the
# total number of draws N (left), where the configurations are spread over a wide range of
# sample sizes, and against n_ov (right), where they collapse; lines join the sample sizes of
# one configuration. The two-point figure shows that both sample-size results are sharp:
# exact curves inside the corridors formed by the lower and upper bounds that apply to them.
if (!exists("BO_ROOT")) {
  BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (length(f)) dirname(normalizePath(f)) else if (dir.exists("analysis")) "analysis" else "."})
  source(file.path(BO_ROOT, "R", "common.R")); source(file.path(BO_ROOT, "R", "paper.R"))
}
if (!have_ggplot() || !requireNamespace("patchwork", quietly = TRUE)) stop("ggplot2 and patchwork are needed for the figures")
suppressPackageStartupMessages(library(ggplot2))
tab <- readRDS(file.path(results_dir, "tab_E6.rds"))
tp <- readRDS(file.path(results_dir, "fig1_universality_two_point.rds"))
pd <- tab[tab$group %in% c("base", "hd"), ]
pd <- pd[order(pd$config, pd$n_ov), ]
pd$N <- pd$S1 + pd$S2
pd$design <- factor(ifelse(pd$group == "hd", "hd", "base"), levels = c("base", "hd"))
pd$family <- factor(pd$family, levels = c("F1", "F2", "F3"))
R0 <- max(pd$R)
cov_me <- 1.96 * sqrt(0.95 * 0.05 / R0)
N10 <- pd$N / pd$n_ov * 10 # draws needed for n_ov = 10, per configuration
res_put("univ.N10.lo", fint(round(min(N10), -1))); res_put("univ.N10.hi", fint(round(max(N10), -1)))
res_put("univ.N10.ratio", fx(max(N10) / min(N10), 0))

fam_col <- c(F1 = "#0072B2", F2 = "#D55E00", F3 = "#009E73")   # Okabe-Ito
design_lab <- c(base = "base configurations", hd = expression(d == 100*","~500*","~I < 0.04))
th <- theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), legend.position = "bottom",
                                       plot.title = element_text(size = 9.5, face = "bold"))
logx <- function(lab) list(scale_x_log10(breaks = 10^(0:6), labels = scales::trans_format("log10", scales::math_format(10^.x))), labs(x = lab))
logn <- list(scale_x_log10(breaks = c(1, 10, 100, 1000), labels = c("1", "10", "100", "1000")), labs(x = expression("effective overlap " * n[ov])))
layer_lines <- list(geom_line(aes(group = config, colour = family, linetype = design), linewidth = 0.55),
                    geom_point(aes(colour = family, shape = design), size = 1.6))
scales_ <- list(scale_colour_manual(values = fam_col, name = NULL),
                scale_linetype_manual(values = c(base = "solid", hd = "22"), labels = design_lab, name = NULL),
                scale_shape_manual(values = c(base = 16, hd = 17), labels = design_lab, name = NULL))
band <- annotate("rect", xmin = 0, xmax = Inf, ymin = 0.95 - cov_me, ymax = 0.95 + cov_me, fill = "grey85", alpha = 0.8)
covy <- list(scale_y_continuous(limits = c(0.69, 0.965), breaks = seq(0.7, 0.95, 0.05)), labs(y = "coverage of nominal 95% intervals"))
mcsey <- list(scale_y_log10(breaks = c(0.003, 0.01, 0.03, 0.1, 0.3), labels = c("0.003", "0.01", "0.03", "0.1", "0.3")),
              labs(y = expression("median |log(MCSE/" * sigma[N] * ")|")))

pa <- ggplot(pd, aes(N, coverage)) + band + geom_hline(yintercept = 0.95, linetype = 3, colour = "grey40") +
  layer_lines + scales_ + logx(expression("number of draws " * N)) + covy + th + ggtitle("Coverage against the number of draws")
pb <- ggplot(pd, aes(n_ov, coverage)) + band + geom_hline(yintercept = 0.95, linetype = 3, colour = "grey40") +
  layer_lines + scales_ + logn + covy + th + ggtitle("... and against the effective overlap")
pc <- ggplot(pd, aes(N, med_abs_log_mcse)) + layer_lines + scales_ + logx(expression("number of draws " * N)) + mcsey + th +
  ggtitle("Error of the MCSE against the number of draws")
tri <- data.frame(x = c(150, 1500, 1500, 150), y = c(0.2, 0.2, 0.2 / sqrt(10), 0.2))     # slope -1/2 indicator
pdd <- ggplot(pd, aes(n_ov, med_abs_log_mcse)) + layer_lines +
  geom_path(data = tri, aes(x, y), inherit.aes = FALSE, colour = "grey30", linewidth = 0.4) +
  annotate("text", x = 480, y = 0.075, label = "slope~-1/2", parse = TRUE, size = 2.7, colour = "grey20", hjust = 1) +
  scales_ + logn + mcsey + th + ggtitle("... and against the effective overlap")

fig <- patchwork::wrap_plots(pa, pb, pc, pdd, ncol = 2, guides = "collect") +
  patchwork::plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(legend.position = "bottom", legend.box = "horizontal")
ggsave(fig_file("fig_universality"), fig, width = 8.2, height = 6.6)
msg("wrote generated/fig_universality.pdf")

## Two-point family of Theorem 2, computed exactly: both sharpness results.
# (a) Failure probability P(Z-hat in {0, inf}) = 1 - {1 - (1 - eps)^S}^2 with N H = S eps / 2,
#     between the lower bound exp(-4 N H) and the upper bound exp(-N H) + exp(-N H (1 - 1/S)),
#     which equals (1 + e^{eps/2}) exp(-N H) here and is at most (1 + e^{1/4}) exp(-N H).
# (b) Exact Kolmogorov distance between 0.04/sqrt(n_ov) (this family, n_ov >= 2) and 3.1/sqrt(n_ov).
nh <- seq(0.05, 9, length.out = 400)
fail <- do.call(rbind, lapply(c(0.05, 0.5), function(e) {
  S <- 2 * nh / e; data.frame(eps = format(e), NH = nh, p = 1 - (1 - (1 - e)^S)^2)
}))
corr_a <- data.frame(NH = nh, lo = exp(-4 * nh), hi = pmin(1, (1 + exp(1 / 4)) * exp(-nh)))
eps_lab <- lapply(c("0.05", "0.5"), function(e) bquote(epsilon == .(as.numeric(e))))
lt <- scale_linetype_manual(values = c("0.05" = "solid", "0.5" = "22"), labels = eps_lab, name = NULL)
qa <- ggplot(fail, aes(NH, p)) +
  geom_ribbon(data = corr_a, aes(NH, ymin = lo, ymax = hi), inherit.aes = FALSE, fill = "grey88") +
  geom_line(data = corr_a, aes(NH, lo), inherit.aes = FALSE, colour = "grey45", linetype = "longdash") +
  geom_line(data = corr_a, aes(NH, hi), inherit.aes = FALSE, colour = "grey45", linetype = "longdash") +
  annotate("text", x = 5.6, y = 2.28 * exp(-5.6) * 2.3, label = "upper bound (every problem)", size = 2.6, colour = "grey30", hjust = 0) +
  annotate("text", x = 0.15, y = 1e-6, label = "lower bound (this family)", size = 2.6, colour = "grey30", hjust = 0) +
  geom_line(aes(linetype = eps), linewidth = 0.7) + lt +
  scale_y_log10(breaks = 10^seq(0, -15, -3), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  coord_cartesian(ylim = c(1e-15, 1)) + th +
  labs(x = expression(N * H * ", the number of ambiguous draws"), y = expression(P(hat(Z) %in% "{0, " * infinity * "}"))) +
  ggtitle("Failure probability: exponential in N H")
tk <- tp[tp$ks > 0, ]; tk$eps <- as.character(tk$eps)     # "0.05", "0.5" (format() would pad to "0.50")
nn <- 10^seq(log10(0.45), log10(1800), length.out = 200)
corr_b <- data.frame(n = nn, lo = ifelse(nn >= 2, 0.04 / sqrt(nn), NA), hi = 3.1 / sqrt(nn))
qb <- ggplot(tk, aes(n_ov, ks)) +
  geom_ribbon(data = corr_b[nn >= 2, ], aes(n, ymin = lo, ymax = pmin(hi, 1)), inherit.aes = FALSE, fill = "grey88") +
  geom_line(data = corr_b, aes(n, hi), inherit.aes = FALSE, colour = "grey45", linetype = "longdash") +
  geom_line(data = corr_b[nn >= 2, ], aes(n, lo), inherit.aes = FALSE, colour = "grey45", linetype = "longdash") +
  annotate("text", x = 25, y = 3.1 / sqrt(25) * 1.5, label = "upper bound (every problem)", size = 2.6, colour = "grey30", hjust = 0) +
  annotate("text", x = 2.3, y = 0.04 / sqrt(2.3) / 1.6, label = "lower bound (this family)", size = 2.6, colour = "grey30", hjust = 0, vjust = 1) +
  geom_line(aes(linetype = eps), linewidth = 0.7) + lt +
  scale_x_log10(breaks = c(1, 10, 100, 1000), labels = c("1", "10", "100", "1000")) +
  scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1), labels = c("0.001", "0.01", "0.1", "1")) +
  coord_cartesian(ylim = c(5e-4, 1)) + th +
  labs(x = expression("effective overlap " * n[ov]), y = "Kolmogorov distance to N(0, 1)") +
  ggtitle(expression(bold("Normal approximation: error of order " * n[ov]^{-1/2})))
fig2 <- patchwork::wrap_plots(qa, qb, nrow = 1, guides = "collect") +
  patchwork::plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") & theme(legend.position = "bottom")
ggsave(fig_file("fig_two_point"), fig2, width = 8.6, height = 3.7)
msg("wrote generated/fig_two_point.pdf")
res_write("plot_universality.R")