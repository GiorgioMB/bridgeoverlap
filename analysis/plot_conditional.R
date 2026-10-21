#!/usr/bin/env Rscript
# Figure fig_conditional.pdf, drawn from the tables saved by
# supplement/S98_operating_characteristics.R (results/tab_OC_binned.rds, tab_OC_curve.rds):
#   Rscript analysis/plot_conditional.R
# Story: a run's own certificate says whether its interval can be trusted. (a) Coverage of the
# runs whose certificate falls in each bin: low for tiny certificates, nominal from about 5 on,
# so the threshold c = 10 leaves a margin. (b) At c = 10, passing runs cover at the nominal level
# and flagged runs do not. Writes the numbers quoted in Section 7.1 (macros ocf.*).
if (!exists("BO_ROOT")) {
  BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (length(f)) dirname(normalizePath(f)) else if (dir.exists("analysis")) "analysis" else "."})
  source(file.path(BO_ROOT, "R", "common.R")); source(file.path(BO_ROOT, "R", "paper.R"))
}
if (!have_ggplot() || !requireNamespace("patchwork", quietly = TRUE)) stop("ggplot2 and patchwork are needed for the figure")
suppressPackageStartupMessages(library(ggplot2))
bi <- readRDS(file.path(results_dir, "tab_OC_binned.rds"))
cu <- readRDS(file.path(results_dir, "tab_OC_curve.rds"))
src <- c("F1--F3, all designs" = "independent-draw designs", "adversarial examples" = "adversarial examples",
         "eight schools, centred" = "eight schools, centred")
key <- c("F1--F3, all designs" = "iid", "adversarial examples" = "adv", "eight schools, centred" = "es")
scol <- c("independent-draw designs" = "#0072B2", "adversarial examples" = "#D55E00", "eight schools, centred" = "#009E73")
bins <- c("[0,2)", "[2,5)", "[5,10)", "[10,20)", "[20,50)", "[50,100)", "[100,Inf)")
blab <- c("< 2", "2-5", "5-10", "10-20", "20-50", "50-100", "> 100")
bi <- bi[bi$runs >= 30 & !is.na(bi$coverage), ]
bi$x <- match(bi$bin, bins); bi$s <- factor(src[bi$source], levels = src)
se <- function(p, n) sqrt(p * (1 - p) / n)
bi$lo <- pmax(0, bi$coverage - 1.96 * se(bi$coverage, bi$runs)); bi$hi <- pmin(1, bi$coverage + 1.96 * se(bi$coverage, bi$runs))
th <- theme_bw(base_size = 9.5) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
                                       legend.position = "bottom", plot.title = element_text(size = 9.5, face = "bold"))
dodge <- position_dodge(width = 0.45)
pa <- ggplot(bi, aes(x, coverage, colour = s)) +
  annotate("rect", xmin = 0.5, xmax = 3.5, ymin = -Inf, ymax = Inf, fill = "grey88") +
  annotate("text", x = 0.6, y = 0.985, label = "flagged: certificate below c = 10", hjust = 0, size = 2.7, colour = "grey30") +
  geom_hline(yintercept = 0.95, linetype = 3, colour = "grey40") +
  geom_line(aes(group = s), position = dodge, linewidth = 0.6) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = dodge, width = 0.25, linewidth = 0.5) +
  geom_point(position = dodge, size = 2) +
  scale_x_continuous(breaks = seq_along(bins), labels = blab, limits = c(0.5, length(bins) + 0.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.68, 1), breaks = seq(0.7, 0.95, 0.05)) +
  scale_colour_manual(values = scol, name = NULL) + th +
  labs(x = expression("the run's certificate  " * hat(n)[ov]), y = "coverage of nominal 95% intervals") +
  ggtitle("Coverage of a run, by its own certificate")

c10 <- cu[cu$c == 10, ]
c10$lab <- sprintf("%s\n(%d%% of runs pass)", src[c10$source], as.integer(round(100 * c10$pass_rate)))
c10$y <- factor(c10$lab, levels = c10$lab[match(rev(names(src)), c10$source)])
wc <- c("runs that pass" = "#0072B2", "runs flagged" = "#D55E00")
db <- rbind(data.frame(y = c10$y, who = "runs that pass", cov = c10$cov_pass, n = c10$n_pass),
            data.frame(y = c10$y, who = "runs flagged", cov = c10$cov_flag, n = c10$runs - c10$n_pass))
db <- db[is.finite(db$cov) & db$n >= 30, ]
db$who <- factor(db$who, levels = names(wc)); db$txt <- sprintf("%.3f", db$cov)
db$lo <- db$cov - 1.96 * se(db$cov, db$n); db$hi <- db$cov + 1.96 * se(db$cov, db$n)
pb <- ggplot(db, aes(cov, y)) +
  geom_vline(xintercept = 0.95, linetype = 2, colour = "grey35") +
  annotate("text", x = 0.951, y = Inf, label = "nominal", hjust = 0, vjust = 1.3, size = 2.9, colour = "grey30") +
  geom_line(aes(group = y), colour = "grey75", linewidth = 1.8) +
  geom_errorbarh(aes(xmin = lo, xmax = hi, colour = who), height = 0.2, linewidth = 0.55) +
  geom_point(aes(colour = who), size = 3.4) +
  geom_text(data = db[db$who == "runs that pass", ], aes(label = txt, colour = who), hjust = -0.35, size = 2.9, show.legend = FALSE) +
  geom_text(data = db[db$who == "runs flagged", ], aes(label = txt, colour = who), hjust = 1.35, size = 2.9, show.legend = FALSE) +
  scale_colour_manual(values = wc, name = NULL) +
  scale_x_continuous(limits = c(0.78, 1.01), breaks = seq(0.8, 1, 0.05)) + th +
  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(), axis.text.y = element_text(size = 8.5)) +
  labs(x = "coverage of nominal 95% intervals", y = NULL) +
  ggtitle("At c = 10: pass versus flagged")
fig <- patchwork::wrap_plots(pa, pb, nrow = 1, widths = c(1.25, 1), guides = "collect") +
  patchwork::plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") & theme(legend.position = "bottom", legend.box = "vertical",
                                                                                          legend.margin = margin(0, 0, 0, 0))
ggsave(fig_file("fig_conditional"), fig, width = 7.4, height = 4.0)
msg("wrote generated/fig_conditional.pdf")

# numbers quoted in Section 7.1
for (s in names(key)) {
  k <- key[[s]]; b <- bi[bi$source == s, ]; x <- c10[c10$source == s, ]
  res_put(sprintf("ocf.%s.lt2", k), if (any(b$bin == "[0,2)")) fx(b$coverage[b$bin == "[0,2)"], 3) else "--")
  res_put(sprintf("ocf.%s.ge5", k), if (any(b$x >= 3)) frange(b$coverage[b$x >= 3], 3) else "--")
  res_put(sprintf("ocf.%s.pass", k), fx(x$cov_pass, 3)); res_put(sprintf("ocf.%s.flag", k), if (is.finite(x$cov_flag)) fx(x$cov_flag, 3) else "--")
  res_put(sprintf("ocf.%s.rate", k), fx(100 * x$pass_rate, 0)); res_put(sprintf("ocf.%s.runs", k), fint(x$runs))
}
res_write("plot_conditional.R")