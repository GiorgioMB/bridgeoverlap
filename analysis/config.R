# Pre-specified replication counts ------------------------------------------------
#
# Every Monte Carlo experiment takes its number of
# replications from this list, and nowhere else. The values were fixed before the
# final run; final outputs come only from a run with BO_SCALE = 1, and every result
# file records the count actually used (generated/tab_replications.tex lists them).
# No run is repeated or replaced: replication r of cell `cell` always uses the seed
# cell_seed(cell) + r (see R/common.R), so a rerun reproduces the same numbers.
REPS <- list(
  ps_laws = 3000, ps = 4000, ps_refine = 50000, rare = 400,   # S21: replications per configuration (S22 is exact: no replications)
    # S20: laws Q searched, replications per law, re-evaluation of the worst laws
  # Section 7.1 (independent draws)
  univ            = 4000L,  # F1-F3 universality grid (Table S-E6, Figure 1)
  univ_hd         = 4000L,  # high-dimensional configurations with I << 0.1
  univ_unbal      = 4000L,  # unbalanced allocations S1 != S2
  univ_threshold  = 6000L,  # configurations near the threshold n_ov = 10 (per-run operating characteristics)
  eb              = 1000L,  # finite-sample interval (Table S-EB)
  np              = 500L,   # near-perfect proposals (Table S-NP)
  adv             = 2000L,  # adversarial mixtures (Table S-ADV)
  adv_funnel      = 1000L,  # funnels, S < 4000
  adv_funnel_big  = 600L,   # funnels, S >= 4000
  e2              = 5000L,  # null calibration (Table S-E2)
  fa              = 2000L,  # joint false-alarm rate (Table S-FA)
  e1              = 1000L,  # contraction factor (Table S-E1, bottom)
  r1              = 10000L, # Proposition 1 checks (Table S-R1)
  r2              = 2000L,  # Theorem 6 checks (Table S-R2)
  ca              = 1000L,  # calibration-test power and blind spot (Section S7.4)
  rerun           = 6000L,  # fresh rerun vs extension after an overlap flag (Table S-rerun)
  # Section 7.2 (Markov chains)
  e9              = 2000L,  # AR(1) universality (Table S-E9, Figure 2)
  e3              = 5000L,  # MCMC null calibration (Table S-E3)
  lrv             = 500L,   # long-run variance estimators (Table S-LRV)
  lrv_iid         = 2000L,  # error_measures() with independent draws
  regen           = 4000L,  # jackknife vs batch means vs regeneration (Table S-REGEN)
  # proposal fitting and cross-fitting
  prop            = 1000L,  # fixed budget, d = 2, 8 (Table S-prop)
  prop_d32        = 500L,   # fixed budget, d = 32
  prop_skew       = 2000L,  # skewed target (residual second-order term)
  e4              = c(2000L, 2000L, 1000L, 600L),                             # Table S-E4 rows
  e4c             = c(1000L, 1000L, 1000L, 1000L, 400L, 400L, 1000L, 1000L, 500L),  # Table S-E4c rows
  # sampler failure
  e5              = 5000L,  # Table S-E5
  # eight schools, exact draws and random-walk Metropolis (Section 7.4)
  es_exact        = 2000L,  # per (parameterization, S) with the bridgesampling workflow
  es_diag         = 2000L,  # diagnosis of the undercoverage: draws x proposal arms
  es_rwm          = 200L,   # random-walk Metropolis
  # Stan benchmark (Section 7.3); each replication is a full NUTS run
  stan_es         = 200L,   # eight schools, per configuration
  stan_radon      = 200L,   # radon models, per configuration
  stan_radon_big  = 100L,   # radon, all counties (d = 390)
  stan_radon_x4   = 100L,   # radon, four times the default number of draws (the overlap action)
  stan_mix        = 200L,   # mixture with label switching
  stan_reg        = 100L,   # conjugate regression, per dimension
  stan_toy_mix    = 100L    # two-mode mixture with controlled initial values
)
