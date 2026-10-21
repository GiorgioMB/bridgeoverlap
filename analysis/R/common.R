# Common setup for the analysis scripts -----------------------------------------
#
# Environment variables
#   BO_SCALE  fraction of the number of replications (default 1).
#             BO_SCALE=0.05 gives a quick smoke run of every script.
#   BO_CORES  number of forked worker processes (Unix only, default 1).
#   BO_SEED   base seed (default 20260917). Replication r of every cell uses
#             seed BO_SEED + r, so results do not depend on BO_CORES.
#   BO_QUICK  "true" keeps only a few cells of every grid (smoke test).

suppressPackageStartupMessages(library(bridgeoverlap))

if (!exists("BO_ROOT")) {
  BO_ROOT <- if (dir.exists("analysis")) normalizePath("analysis") else normalizePath(".")
}
BO_SCALE <- as.numeric(Sys.getenv("BO_SCALE", "1"))
BO_CORES <- as.integer(Sys.getenv("BO_CORES", "1"))
BO_SEED <- as.integer(Sys.getenv("BO_SEED", "20260917"))
# BO_QUICK=true trims the grids of every script to a few cells (smoke test).
BO_QUICK <- identical(tolower(Sys.getenv("BO_QUICK", "false")), "true")
quick_pick <- function(x, keep = c(1L, length(x))) if (BO_QUICK) x[unique(pmin(keep, length(x)))] else x

results_dir <- file.path(BO_ROOT, "results")
figures_dir <- file.path(BO_ROOT, "figures")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

source(file.path(BO_ROOT, "config.R"))

# Replication counts come from config.R (pre-specified). BO_SCALE < 1 is for smoke
# tests only; every lookup is logged so that the assembly step can list the counts
# used and refuse to label a scaled run as final.
.reps_log <- new.env()
n_reps <- function(n, min = 20L, key = NULL) {
  if (is.null(key)) key <- paste(deparse(substitute(n)), collapse = "")
  out <- pmax(as.integer(min), as.integer(round(n * BO_SCALE)))
  assign(key, list(key = key, specified = paste(n, collapse = "/"), used = paste(out, collapse = "/")), envir = .reps_log)
  out
}
write_reps_log <- function(script) {
  keys <- ls(.reps_log)
  if (!length(keys)) return(invisible(NULL))
  df <- do.call(rbind, lapply(keys, function(k) as.data.frame(get(k, envir = .reps_log), stringsAsFactors = FALSE)))
  df$script <- script; df$BO_SCALE <- BO_SCALE; df$BO_SEED <- BO_SEED
  utils::write.csv(df, file.path(results_dir, paste0("replications_", gsub("[^A-Za-z0-9]+", "_", basename(script)), ".csv")),
                   row.names = FALSE)
}

# Independent random-number streams per cell: replication r of the cell with key
# `cell` uses set.seed(cell_seed(cell) + r). Different cells (configurations,
# arms, scripts) therefore never share draws, which matters when runs are pooled
# across cells (conditional coverage curves).
cell_seed <- function(cell) {
  if (is.null(cell)) return(BO_SEED)
  h <- 0
  for (code in utf8ToInt(as.character(cell))) h <- (h * 131 + code) %% 19997
  BO_SEED + 100000L * as.integer(h)
}

msg <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")
  utils::flush.console()
}

# Run replications 1..R of `fun(r)`, each with its own seed, and row-bind the
# data frames they return. Failed replications are reported and dropped.
BO_SCRIPT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) basename(f) else "interactive"})
.rr_counter <- new.env(); .rr_counter$k <- 0L
# Replications to run in this process: all of them, or the subset BO_REP_RANGE = "a-b"
# (used by the cluster array jobs of the Stan benchmark; results are combined later).
rep_subset <- function(R) {
  rg <- Sys.getenv("BO_REP_RANGE", "")
  if (!nzchar(rg)) return(seq_len(R))
  ab <- as.integer(strsplit(rg, "-")[[1]])
  intersect(seq_len(R), seq(ab[1], ab[2]))
}
# Cluster support: with BO_REP_RANGE set, the replications run in this process are also
# cached per cell; with BO_COMBINE=true, run_reps() returns the cached replications of all
# processes (in replication order) instead of computing, so that rerunning a script in
# combine mode writes its tables from the union of the array tasks.
.cache_file <- function(cell, tag) {
  d <- file.path(results_dir, "cache"); dir.create(d, showWarnings = FALSE, recursive = TRUE)
  file.path(d, sprintf("%s__%s__%s.rds", sub("[.]R$", "", BO_SCRIPT), gsub("[^A-Za-z0-9.]+", "_", cell), tag))
}
run_reps <- function(R, fun, cell = NULL, seed = NULL, cores = BO_CORES, reps = seq_len(R)) {
  .rr_counter$k <- .rr_counter$k + 1L
  if (is.null(cell)) cell <- paste(BO_SCRIPT, .rr_counter$k)
  if (is.null(seed)) seed <- cell_seed(cell)
  if (identical(tolower(Sys.getenv("BO_COMBINE")), "true")) {
    pat <- sub("__[^_]*[.]rds$", "", basename(.cache_file(cell, "x")))
    fs <- list.files(file.path(results_dir, "cache"), pattern = paste0("^", gsub("([.])", "\\\\\\1", pat), "__"), full.names = TRUE)
    # A cell without cached replications means an array task did not run or failed: stop
    # rather than recompute it here (the finish job is not sized for that).
    if (!length(fs)) stop("BO_COMBINE: no cached replications for cell '", cell, "' of ", BO_SCRIPT,
                          "; rerun the corresponding line of hpc/stan_tasks.txt", call. = FALSE)
    if (length(fs)) {
      parts <- lapply(fs, readRDS)
      rr <- unlist(lapply(parts, `[[`, "reps")); outs <- unlist(lapply(parts, `[[`, "out"), recursive = FALSE)
      keep <- !duplicated(rr) & rr <= R
      msg("combined", sum(keep), "replications of", cell, "from", length(fs), "file(s)")
      if (sum(keep) < length(reps)) msg("WARNING:", length(reps) - sum(keep), "of", length(reps), "replications of", cell,
                                        "are missing (failed replications or an unfinished array task)")
      return(do.call(rbind, outs[keep][order(rr[keep])]))
    }
  }
  one <- function(r) {
    set.seed(seed + r)
    tryCatch(fun(r), error = function(e) structure(conditionMessage(e), class = "rep_error"))
  }
  out <- if (cores > 1L && .Platform$OS.type == "unix") {
    parallel::mclapply(reps, one, mc.cores = cores, mc.preschedule = TRUE)
  } else {
    lapply(reps, one)
  }
  bad <- vapply(out, function(o) inherits(o, "rep_error") || inherits(o, "try-error"), logical(1))
  if (any(bad)) {
    msg(sum(bad), "of", length(reps), "replications failed; first error:", as.character(out[[which(bad)[1]]]))
  }
  rg <- Sys.getenv("BO_REP_RANGE", "")
  if (nzchar(rg)) saveRDS(list(reps = reps[!bad], out = out[!bad]), .cache_file(cell, gsub("[^0-9]+", "-", rg)))
  do.call(rbind, out[!bad])
}

save_table <- function(x, name) {
  utils::write.csv(x, file.path(results_dir, paste0(name, ".csv")), row.names = FALSE)
  saveRDS(x, file.path(results_dir, paste0(name, ".rds")))
  msg("wrote", file.path("results", paste0(name, ".csv")))
  invisible(x)
}

# Kolmogorov distance between the empirical law of z (which may contain +-Inf,
# e.g. when Z-hat is 0 or infinite) and N(0, 1).
ks_normal <- function(z) {
  n <- length(z)
  n_minus <- sum(z == -Inf)
  n_plus <- sum(z == Inf)
  x <- sort(z[is.finite(z)])
  if (!length(x)) return(max(n_minus, n_plus) / n)
  Fr <- (n_minus + seq_along(x)) / n
  Fl <- (n_minus + seq_along(x) - 1) / n
  P <- stats::pnorm(x)
  max(abs(Fr - P), abs(Fl - P), n_minus / n, n_plus / n)
}

# Kolmogorov distance between a discrete law (atoms x with weights w, plus
# masses at -Inf and +Inf) and N(0, 1). Ties are merged.
ks_discrete <- function(x, w, p_minus = 0, p_plus = 0) {
  o <- order(x)
  x <- x[o]
  w <- w[o]
  rx <- round(x, 10)
  first <- !duplicated(rx)
  grp <- cumsum(first)
  W <- as.numeric(rowsum(w, grp))
  xu <- x[first]
  Fr <- p_minus + cumsum(W)
  Fl <- c(p_minus, utils::head(Fr, -1))
  P <- stats::pnorm(xu)
  max(abs(Fr - P), abs(Fl - P), p_minus, p_plus)
}

weighted_median <- function(x, w) {
  o <- order(x)
  cw <- cumsum(w[o]) / sum(w)
  x[o][which(cw >= 0.5)[1]]
}

se_prop <- function(p, n) sqrt(p * (1 - p) / n)

# bridgesampling::error_measures() for method "normal", re-implemented on the
# log ratios (f1 = (1 - pi)/s1, f2 = pi/s2) so that it can be applied to many
# replications quickly. The posterior term multiplies the spectral density at
# zero by var(f2) a second time, which is the normalization issue discussed in
# Section 7.2. Checked against the package in 00_check_setup.R.
error_measures_re2 <- function(a_post, a_prop, logml) {
  N1 <- length(a_post)
  N2 <- length(a_prop)
  s1 <- N1 / (N1 + N2)
  s2 <- N2 / (N1 + N2)
  lr <- log(s1 / s2)
  f1 <- stats::plogis(a_prop - logml + lr) / s1
  f2 <- stats::plogis(logml - a_post - lr) / s2
  rho_f2 <- coda::spectrum0.ar(f2)$spec
  stats::var(f1) / mean(f1)^2 / N2 + rho_f2 / N1 * stats::var(f2) / mean(f2)^2
}

# Locate the analysis directory from a script run with Rscript.
script_root <- function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else BO_ROOT
}

have_ggplot <- function() requireNamespace("ggplot2", quietly = TRUE)

save_plot <- function(p, name, width = 10, height = 4) {
  f <- file.path(figures_dir, paste0(name, ".pdf"))
  ggplot2::ggsave(f, p, width = width, height = height)
  msg("wrote", file.path("figures", paste0(name, ".pdf")))
}
