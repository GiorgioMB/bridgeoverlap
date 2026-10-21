# Analysis scripts

**Replication counts** are fixed in `config.R` (one key per experiment), and every script
logs the counts it used (`results/replications_*.csv` and `generated/tab_replications.tex`). **Seeds:**
replication *r* of the cell with key `cell` uses `set.seed(cell_seed(cell) + r)`
(`R/common.R`), so cells never share draws and a rerun reproduces every number. **Cluster:**
see `../hpc/README.md`.

Every script writes CSV and RDS files to `analysis/results/` and its LaTeX outputs to
`generated/` (or `BO_OUT`):
- tables, each as a complete `tabular` (`tab_*.tex`);
- figures (`fig_*.pdf`);
- the numbers quoted in the text, as `\resdef{key}{value}` lines (`results_<script>.tex`).

The last step, `supplement/S99_assemble.R`, builds the main-text Stan table, collects all
result macros into `generated/results.tex`, and writes `generated/MANIFEST.txt`,
which lists the files the LaTeX sources expect. Run everything from the repository root after
installing the package (`remotes::install_local(".")`):

```sh
Rscript analysis/run_all.R --stan     # full run: every table, figure and reported number
Rscript analysis/run_all.R            # everything that needs no Stan
BO_QUICK=true BO_SCALE=0.01 Rscript analysis/run_all.R --stan   # smoke test
```


| Variable | Default | Effect |
|---|---|---|
| `BO_SCALE` | `1` | Fraction of the replication counts of `config.R` (minimum 20 per cell). |
| `BO_QUICK` | `false` | `true` keeps only a few cells of every grid (smoke test). |
| `BO_CORES` | `1` | Forked workers for the pure-R scripts (01–05 and `supplement/`; Unix only). |
| `BO_SEED` | `20260917` | Replication *r* of every cell uses seed `BO_SEED + r`. |

```sh
BO_QUICK=true BO_SCALE=0.01 Rscript analysis/04_sampler_failure.R   # seconds
BO_CORES=8 Rscript analysis/04_sampler_failure.R                    # full scale
```

## Pure-R experiments (no Stan)

| Script | Main outputs | What to compare |
|---|---|---|
| `01_universality.R` | `tab_E6.csv`, `fig1_universality_two_point.csv`, `fig1_universality.pdf` | Coverage within 0.01 of nominal for n_ov ≥ 3.7. Median relative MCSE error ≈ C/√n_ov with C in [0.11, 0.23]. Two-point Kolmogorov distance ≈ 0.1/√n_ov. Prop. 2 interval always covers at 2.9–3.5× the width. Any failure flag fires ≈ 15% at 5% each and 3.6–5.2% at α/4. |
| `02_markov_chains.R` | `tab_E9.csv`, `fig2_*.csv`, `fig2_markov_chains.pdf` | Ignoring autocorrelation covers 2Φ(1.96√γ) − 1 (0.34 at γ = 0.05). True-σ_NP coverage is 0.937–0.971. Spectral coverage is within 0.015 for n_ov^P ≥ 20. `error_measures()` undercovers. γτ_W lies in [0.38, 0.93]. |
| `03_proposal_fitting.R` | `tab_prop.csv` | Same-draw bias exceeds one MCSE at d = 32 (coverage 0.77). Cross-fitting costs at most 5% in SD. |
| `04_sampler_failure.R` | `tab_E5.csv` | Bias 2.42–9.72 MCSEs. Calibration test rejects 46/97/100%. Var / MCSE² up to 4.72. V_jack ≈ 0.63 for slow chains. |
| `03_proposal_fitting.R` (skewed target) | `tab_prop.csv` | Residual second-order term of cross-fitting; cost relative to the fixed and the same-draw proposal. |
| `05_eight_schools.R` | `tab_eight_schools_exact.csv`, `tab_eight_schools_rwm.csv` | At S = 8000: centred n_ov ≈ 456, I ≈ 0.35, RMSE 0.034; non-centred n_ov ≈ 159, I ≈ 0.91. Random-walk Metropolis acceptance ≈ 0.08. |

Two details of the set-up are choices made here:
- **Proposal fitting (`03_proposal_fitting.R`):** the split arm bridges 1000 posterior draws against 2000
  proposal draws, keeping the proposal budget fixed.
- **Sampler failure (`04_sampler_failure.R`):** the within-component AR coefficient of the metastable
  chain is 0.5.

The within-chain autocorrelation time of the random-walk Metropolis runs is
lower here than in Table ES (median about 70 rather than 117 at S = 8000). Many
of those runs contain chains that never move, so the value depends on how
estimators treat such chains. The package reports them as `stuck_chains`.

## Further experiments

| Script | Output |
|---|---|
| `01_universality.R` (groups `hd`, `unbal`, `threshold`) | d = 100, 500 with I << 0.1; unbalanced allocations; per-run operating characteristics; per-family MCSE constants |
| `supplement/S17_ADV_replication_check.R` | every ADV row rerun with an independent stream |
| `supplement/S18_LRV_regeneration.R` | jackknife vs batch means vs regeneration |
| `supplement/S98_operating_characteristics.R` | pooled P(cover given pass at c) |
| `30_stan_benchmark.R` + `31_stan_benchmark_summary.R` | posteriordb posteriors with exact evidence (`R/radon.R`, `R/mixture.R`); checks of current practice; Bayes factors |
| `hpc_combine.R` | rebuilds the tables of scripts 10-12 from cluster array tasks |

`20_posteriordb_suite.R` is no longer run by `run_all.R`: its references came from the same
estimator, and the benchmark of script 30 uses exact evidences instead.

## Numerical checks and smaller studies

Every further table, numerical check and computer-assisted step has a script in
`analysis/supplement/`. "Exact" rows are deterministic computations, reproduced to the
printed digits; all other rows are Monte Carlo and agree up to simulation noise at full scale.

| Table | Script | Kind | Smoke-test agreement |
|---|---|---|---|
| E1 | `S01_E1_identities.R` | single draws + 1000 replications | Identities hold to all digits. The limit sqrt(beta) sd = 0.2018 matches 0.202. |
| E6 | `../01_universality.R` | 4000 replications | H and n_ov exact (e.g. 1.12e-01, 1.24). |
| E7 | `S02_E7_two_point.R` | exact | All printed digits, e.g. coverage 0.4043, 0.9595 and 0.9510. |
| EB | `S03_EB_interval.R` | 1000 replications | Width ratios 2.93 / 3.27. |
| NP | `S04_NP_near_perfect.R` | 500 replications | Median Lyapunov ratio 0.080 / 0.024. |
| ADV | `S05_ADV_adversarial.R` | 300-10000 replications | True n_ov exact (0.49, 31.4, 1.9, 191.6); see note 2. |
| E2 | `S06_E2_null_calibration.R` | 5000 replications | |
| FA | `S07_FA_false_alarm.R` | 2000 replications | |
| E3 | `S08_E3_mcmc_null.R` | 5000 replications | |
| E8 | `S09_E8_two_state_chain.R` | exact | All printed digits, e.g. Kolmogorov 0.1160 and coverage 0.8094 / 0.3083. |
| E9 | `../02_markov_chains.R` | 2000 replications | tau_W exact (1.24 ... 18.02). |
| LRV | `S10_LRV_estimators.R` | 500 replications | tau_W 9.04 and 4.04 exact. |
| E4 | `S11_E4_proposal_fitting.R` | 600-2000 replications | R_g 0.34 / 0.19; see note 1. |
| prop | `../03_proposal_fitting.R` | 500-1000 replications | Median n_ov 12 / 13 / 9. |
| E4c | `S12_E4c_same_draw_bias.R` | 400-1000 replications | Median n_ov 148 and I 0.82. |
| E5 | `../04_sampler_failure.R` | 5000 replications | Median n_ov 458 / 417 / 914 and Var/MCSE^2 4.67. |
| R1 | `S13_R1_certificate_bounds.R` | 10^4 replications | n_ov, NH and both bounds exact (9.30, 9.48, 49.2). |
| R2 | `S14_R2_calibration_null.R` | 2000 replications | n_ov exact (247.5, 85.8). |
| calibration-test checks | `S15_calibration_test_checks.R` | quadrature + 1000 replications | Quadrature ratios 1.000-1.005 and 1.002-1.041 exact; see note 3. |
| computer-assisted steps | `S16_constants.R` | numerical | C'_P margins 1.140 / 1.094; Lemma S3 holds on 300 chains; b(gamma) argmax k = 3 with margin gamma^4/288; see note 4. |


Result keys and the scripts that define them:

| Prefix | Script | Prefix | Script |
|---|---|---|---|
| `univ.` | `01_universality.R` | `e2.` | `S06_E2_null_calibration.R` |
| `e9.` | `02_markov_chains.R` | `fa.` | `S07_FA_false_alarm.R` |
| `prop.` | `03_proposal_fitting.R` | `e3.` | `S08_E3_mcmc_null.R` |
| `e5.` | `04_sampler_failure.R` | `e8.` | `S09_E8_two_state_chain.R` |
| `es.` | `05_eight_schools.R` | `lrv.` | `S10_LRV_estimators.R` |
| `e1.` | `S01_E1_identities.R` | `e4.` | `S11_E4_proposal_fitting.R` |
| `e7.` | `S02_E7_two_point.R` | `e4c.` | `S12_E4c_same_draw_bias.R` |
| `eb.` | `S03_EB_interval.R` | `r1.`, `r2.` | `S13_R1_*`, `S14_R2_*` |
| `np.` | `S04_NP_near_perfect.R` | `ca.`, `bs.` | `S15_calibration_test_checks.R` |
| `adv.` | `S05_ADV_adversarial.R` | `cp.`, `lazy.`, `b.`, `mbar.` | `S16_constants.R` |
| `stan.es.` | `10_stan_eight_schools.R` | `stan.reg.` | `11_stan_regression_dimension.R` |
| `stan.mix.` | `12_stan_mixture.R` | `pdb.` | `20_posteriordb_suite.R` |

## cmdstanr studies (need CmdStan)

These run the diagnostics on Stan's NUTS output with `bridgesampling`. Every target except
those in script 20 has an exact log Z. Each model is compiled once per session.

| Script | Model(s) | Question |
|---|---|---|
| `10_stan_eight_schools.R` (Tables Stan, S-stan-es) | Eight schools, centred (adapt_delta 0.8 and 0.99) and non-centred, S = 2000 and 8000 | Does the certificate flag the funnel under NUTS, and is the MCSE calibrated for the non-centred model? Compares three autocorrelation estimators and `error_measures()`. |
| `11_stan_regression_dimension.R` (Tables Stan, S-stan-reg) | Conjugate regression of `MASS::mcycle` on B-splines, dimension 5–65 | How overlap and calibration change with dimension (near-perfect regime to poor overlap), and split versus cross-fitting on the same fits. |
| `12_stan_mixture.R` (Tables Stan, S-stan-mix) | Two-component mixture with initial values controlled per chain | Metastability (m = 3), mode imbalance (m = 5, inits 3:1) and the blind spot (m = 5, one mode found). |
| `20_posteriordb_suite.R` (Table S-pdb) | 12 posteriordb posteriors (configurable) | Calibration against long cross-fitted reference runs, in the style of Micaletto and Vehtari (2025). Data come from a local clone (`PDB_PATH`) or `raw.githubusercontent.com`. `BO_PDB_DRYRUN=true` only fetches. |

Every script runs end to end in smoke-test mode (`BO_QUICK=true`); full-scale runs take hours
on one core, see [`../hpc/README.md`](../hpc/README.md).
