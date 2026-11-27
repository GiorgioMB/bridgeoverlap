# Reproducing the paper

Every number, table and figure of the paper and its supplement is written by the scripts
in this directory to `generated/` at the repository root: tables as complete `tabular`
environments (`generated/tab_*.tex`), figures as PDF (`generated/fig_*.pdf`), and every
number quoted in the text as `\resdef{key}{value}` in `generated/results.tex`, which both
LaTeX sources input.

## Running

```sh
R CMD INSTALL .                                   # the package, from the repository root
Rscript analysis/run_all.R                        # every study that needs no Stan
Rscript analysis/run_all.R --stan                 # full run: also the Stan benchmark (30, 31)
BO_QUICK=true BO_SCALE=0.01 Rscript analysis/run_all.R   # smoke test, minutes
```

Environment variables: `BO_QUICK` (reduced grids), `BO_SCALE` (fraction of the
replications), `BO_CORES` (parallel workers), `BO_OUT` (output directory). Only a run
with `BO_SCALE=1` and no `BO_QUICK` is final; `supplement/S99_assemble.R` records this in
the macro `prov.final` and lists the expected outputs in `generated/MANIFEST.txt`.

Replication counts are fixed in advance in `config.R`. Replication `r` of a cell always
uses the seed `cell_seed(cell) + r` (`R/common.R`), so a rerun reproduces every number;
the counts actually used are listed in Table S-replications. The Stan benchmark needs
`cmdstanr` and CmdStan; `hpc/` runs it as a SLURM array.

## Where each table and figure comes from

| Paper | Output | Script |
|---|---|---|
| Figure 1 | `fig_universality` | `01_universality.R`, `plot_universality.R` |
| Figure 2 | `fig_bench` | `30_stan_benchmark.R`, `31_stan_benchmark_summary.R`, `plot_bench.R` |
| Table S-E7 | `tab_E7` | `supplement/S02_E7_two_point.R` |
| Table S-EB | `tab_EB` | `supplement/S03_EB_interval.R` |
| Table S-NP | `tab_NP` | `supplement/S04_NP_near_perfect.R` |
| Figure S-conditional | `fig_conditional` | `supplement/S98_operating_characteristics.R`, `plot_conditional.R` |
| Table S-E8, Figure S-mcmc (top) | `tab_E8`, `fig_mcmc` | `supplement/S09_E8_two_state_chain.R`, `plot_mcmc.R` |
| Table S-E9, Figure S-mcmc (bottom) | `tab_E9_top`, `tab_E9_bottom` | `02_markov_chains.R`, `plot_mcmc.R` |
| Figure S-collapse | `fig_collapse` | `supplement/S23_collapse.R`, `plot_collapse.R` |
| Table S-khat | `tab_khat` | `supplement/S24_pareto_khat.R` |
| Table S-mbar | `tab_mbar` | `supplement/S26_multiple_bridges.R` |
| Table S-LRV | `tab_LRV` | `supplement/S10_LRV_estimators.R` |
| Table S-prop | `tab_prop` | `03_proposal_fitting.R` |
| Table S-E5 | `tab_E5` | `04_sampler_failure.R` |
| Tables S-bench, S-bf, S-detect | `tab_bench`, `tab_bf`, `tab_detect` | `31_stan_benchmark_summary.R` |
| Table S-replications | `tab_replications` | `supplement/S99_assemble.R` |

Numbers quoted in the text come from the same scripts and from:
`05_eight_schools.R` (eight schools with independent draws and random-walk Metropolis),
`printout.R` (the worked example of the decision rule),
`supplement/S05_ADV_adversarial.R` (adversarial examples),
`supplement/S06_E2_null_calibration.R` and `S07_FA_false_alarm.R` (null calibration and
joint false-alarm rate), `S11_E4_proposal_fitting.R` and `S12_E4c_same_draw_bias.R`
(proposal fitting), `S18_LRV_regeneration.R` (regeneration), `S19_rerun_vs_extend.R`
(fresh rerun or extension after an overlap flag), `S22_exact_mcse_counterexample.R` (the
counting example of Section 6.2), `S23_collapse.R` (Section S5.1) and `S25_mcmc_interval.R` (the
finite-sample interval for Markov chain output).
