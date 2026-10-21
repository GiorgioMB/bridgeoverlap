#!/usr/bin/env Rscript
# cmdstanr study 4 (optional): the diagnostics across posteriordb posteriors, in
# the style of the suite of Micaletto and Vehtari (2025). The Stan programs use
# `~` statements, whose constants cmdstanr's log_prob() drops, so log Z is the
# normalizing constant of the program's log density; this does not affect the
# comparison because the reference is computed with the same program.
#
# Reference: the mean of R_ref long cross-fitted runs (4 chains x 8000 draws).
# Evaluation: R runs at S = 4000 and S = 16000 with the default split.
# Initial values from Pathfinder, as in the earlier suite (BO_PDB_PATHFINDER=false
# to use Stan's defaults).
#
# Posteriors are read from a local clone (PDB_PATH=/path/to/posteriordb) or
# downloaded from raw.githubusercontent.com (no GitHub API calls). Choose them
# with BO_PDB_POSTERIORS="name1,name2"; BO_PDB_DRYRUN=true only fetches them.
BO_ROOT <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else if (dir.exists("analysis")) "analysis" else "."})
source(file.path(BO_ROOT, "R", "common.R"))
source(file.path(BO_ROOT, "R", "stan.R"))
source(file.path(BO_ROOT, "R", "paper.R"))

default_posteriors <- c(
  "eight_schools-eight_schools_noncentered", "earnings-logearn_interaction",
  "gp_pois_regr-gp_pois_regr", "low_dim_gauss_mix-low_dim_gauss_mix", "arK-arK",
  "garch-garch11", "nes2000-nes", "dogs-dogs", "GLMM_Poisson_data-GLMM_Poisson_model",
  "mcycle_gp-accel_gp", "Mtbh_data-Mtbh_model", "irt_2pl-irt_2pl")
posteriors <- strsplit(Sys.getenv("BO_PDB_POSTERIORS", paste(default_posteriors, collapse = ",")), ",")[[1]]
posteriors <- quick_pick(trimws(posteriors), 1L)
dry_run <- identical(tolower(Sys.getenv("BO_PDB_DRYRUN", "false")), "true")
use_pathfinder <- !identical(tolower(Sys.getenv("BO_PDB_PATHFINDER", "true")), "false")
pdb_path <- Sys.getenv("PDB_PATH", "")
cache <- file.path(BO_ROOT, "cache", "posteriordb")
dir.create(cache, showWarnings = FALSE, recursive = TRUE)

pdb_file <- function(rel) {
  dest <- file.path(cache, basename(rel))
  if (nzchar(pdb_path)) {
    src <- file.path(pdb_path, "posterior_database", rel)
    if (!file.exists(src)) stop("Not found in PDB_PATH: ", src)
    return(src)
  }
  if (!file.exists(dest)) {
    url <- paste0("https://raw.githubusercontent.com/stan-dev/posteriordb/master/posterior_database/", rel)
    utils::download.file(url, dest, quiet = TRUE, mode = "wb")
  }
  dest
}

fetch_posterior <- function(name) {
  info <- jsonlite::read_json(pdb_file(file.path("posteriors", paste0(name, ".json"))))
  code <- pdb_file(file.path("models", "stan", paste0(info$model_name, ".stan")))
  zip <- pdb_file(file.path("data", "data", paste0(info$data_name, ".json.zip")))
  exdir <- file.path(cache, paste0(info$data_name, "_data"))
  json <- utils::unzip(zip, exdir = exdir)
  stan_file <- file.path(cache, paste0(info$model_name, ".stan"))
  if (!file.exists(stan_file)) file.copy(code, stan_file)
  list(name = name, stan_file = stan_file,
       data = jsonlite::read_json(json[grepl("\\.json$", json)][1], simplifyVector = TRUE))
}

pdbs <- lapply(posteriors, function(nm) {
  p <- tryCatch(fetch_posterior(nm), error = function(e) { msg("could not fetch", nm, ":", conditionMessage(e)); NULL })
  if (!is.null(p)) msg("fetched", nm, "(data:", paste(names(p$data), collapse = ", "), ")")
  p
})
pdbs <- Filter(Negate(is.null), pdbs)
if (dry_run) quit(save = "no")

stan_ready()
R <- n_reps(100, 10)
R_ref <- n_reps(20, 5)
invisible(jackknife_null_quantiles(20L, probs = 1 - 0.05 / 4))
invisible(jackknife_null_quantiles(16L, probs = 1 - 0.05 / 4))
make_init <- function(mod, data, seed) {
  if (!use_pathfinder) return(NULL)
  tryCatch(mod$pathfinder(data = data, seed = seed, refresh = 0, show_messages = FALSE,
                          num_paths = 4), error = function(e) NULL)
}

summaries <- list()
for (p in pdbs) {
  res <- tryCatch({
    mod <- stan_model(p$stan_file)
    msg(p$name, ": reference (", R_ref, "cross-fitted runs of 4 x 8000 draws )")
    ref <- run_reps(R_ref, function(r) {
      sf <- stan_fit(mod, p$data, 32000L, seed = BO_SEED + 10000L + r,
                     init = make_init(mod, p$data, BO_SEED + 10000L + r))
      on.exit(stan_cleanup(sf$fit))
      d <- suppressMessages(bridge_diagnostics(bridge_sampler_cmdstanr(sf$fit, crossfit = TRUE), G = 16L))
      data.frame(logml = d$logml, mcse = d$mcse)
    }, cores = 1L)
    ref_logml <- mean(ref$logml)
    out <- lapply(c(4000L, 16000L), function(S) {
      msg(p$name, ": evaluation at S =", S)
      runs <- run_reps(R, function(r) {
        sf <- stan_fit(mod, p$data, S, seed = BO_SEED + r, init = make_init(mod, p$data, BO_SEED + r))
        on.exit(stan_cleanup(sf$fit))
        cbind(posterior = p$name, S = S, rep = r, divergent = sf$divergent,
              stan_summary_row(bridge_sampler_cmdstanr(sf$fit), ref_logml))
      }, cores = 1L)
      s <- summarise_stan(runs, c("posterior", "S"))
      cbind(s, ref_logml = ref_logml, ref_sd = stats::sd(ref$logml),
            ref_mean_mcse = mean(ref$mcse))
    })
    do.call(rbind, out)
  }, error = function(e) { msg("posterior", p$name, "failed:", conditionMessage(e)); NULL })
  if (!is.null(res)) {
    summaries[[p$name]] <- res
    save_table(do.call(rbind, summaries), "posteriordb_suite_summary")
  }
}
if (length(summaries)) {
  tab <- do.call(rbind, summaries)
  print(tab[, c("posterior", "S", "R", "coverage", "rmse", "med_mcse", "rmse_over_mcse", "med_n_ov",
                "med_I", "med_n_ovP", "flag_overlap", "flag_any_failure")], digits = 3)
}

## LaTeX outputs ---------------------------------------------------------------------------------
if (length(summaries)) {
  tab <- do.call(rbind, summaries)
  tt <- function(x) sprintf("\\texttt{%s}", gsub("_", "\\\\_", x))
  emit_table("tab_pdb", cbind(ifelse(!duplicated(tab$posterior), tt(tab$posterior), ""), tab$S, fx(tab$coverage, 3), fx(tab$rmse, 3),
                              fx(tab$med_mcse, 3), fx(tab$rmse_over_mcse, 2), fint(tab$med_n_ov), fx(tab$med_I, 2),
                              fx(tab$med_n_ovP, 1), fx(100 * tab$flag_overlap, 0), fx(100 * tab$flag_any_failure, 0)),
             groups = group_ends(tab$posterior), spec = "lrcccccrcrr",
             header = "Posterior & $S$ & coverage & RMSE & med MCSE & RMSE/MCSE & med $\\neffhat$ & med $\\Ihat$ & med $\\neffPhat$ & overlap flag & failure flags\\\\")
  res_put("pdb.n", length(unique(tab$posterior))); res_put("pdb.R", fint(R)); res_put("pdb.Rref", fint(R_ref))
  res_put("pdb.cov", frange(tab$coverage, 2)); res_put("pdb.ratio", frange(tab$rmse_over_mcse, 2))
  flagged <- tab$flag_overlap > 0.5 | tab$flag_any_failure > 0.5
  res_put("pdb.nflag", sum(flagged)); res_put("pdb.nrows", nrow(tab))
  ok <- !flagged
  res_put("pdb.cov.ok", frange(tab$coverage[ok], 2))
  res_write("20_posteriordb_suite.R")
}
