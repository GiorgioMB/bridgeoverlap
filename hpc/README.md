# Running the full reproduction on a SLURM cluster

Every table, figure and reported number is regenerated from tagged code with the replication counts of
`analysis/config.R`. Nothing needs to be edited for the final run.

```sh
git checkout <release tag>
R CMD INSTALL .                 # bridgeoverlap and its Stan programs
Rscript hpc/check_setup.R       # R packages, CmdStan, one Stan program with model methods
hpc/submit.sh                   # runs the check again, then submits three dependent jobs
```

`hpc/check_setup.R` lists any missing R package with the command to install it. It needs
CmdStan (`cmdstanr::install_cmdstan()` or `CMDSTAN`) and about 3 GB of memory to compile.

| Job | Script | Resources (defaults) | What it does |
|---|---|---|---|
| pure R | `hpc/pure_r.slurm` | 1 node, 64 cores, 64 GB, 12 h | `analysis/run_all.R` without `--stan`: scripts 01-05, supplement S01-S22 and S98; replications on `BO_CORES` forked workers |
| Stan array | `hpc/stan_array.slurm` | 31 tasks x 1 core, 6 GB, 24 h | one line of `hpc/stan_tasks.txt` per task: script, configurations and replication range |
| finish | `hpc/finish.slurm` | 2 cores, 8 GB, 2 h | `analysis/hpc_combine.R` (scripts 10-12), `31_stan_benchmark_summary.R`, `S99_assemble.R`, checks that the run is final |

Adjust the `#SBATCH` lines to your cluster (partition, account, module loads for R, a C++
compiler). Compute nodes need no internet access: the posteriordb data are in
`analysis/data/`.

**Compilation.** A new R session always rebuilds a Stan executable when it compiles the
model methods, so every task compiles its programs in its own directory
(`BO_STAN_BUILD_DIR`, under the task's `TMPDIR`); tasks never write into the package
library. This costs about a minute per program per task.

**If a task fails.** The finish job stops with a message naming the missing cell. Rerun
that line of `hpc/stan_tasks.txt` (`sbatch --array=<line - 1> hpc/stan_array.slurm`, counting
lines without `#`), then `sbatch hpc/finish.slurm`.

**Splitting work.** `hpc/stan_tasks.txt` only splits replications; their number and seeds are
fixed by `analysis/config.R` and `analysis/R/common.R`, so a different split gives identical
results. Each task caches its replications under `analysis/results/cache/` (scripts 10-12) or
writes `analysis/results/bench_runs_<config>_<range>.rds` (script 30).

**Checks.**
- `generated/MANIFEST.txt` lists every output file, and says whether
  the run is final (every count at `BO_SCALE = 1`).
- The finish job stops if the run is not final.

**Smoke test** (minutes, one core, nothing to configure):

```sh
BO_QUICK=true BO_SCALE=0.01 Rscript analysis/run_all.R
```

**Rough cost.** Single runs on one 2-3 GHz core:
- Stan runs (sampling, bridge sampling, warp-III and cross-fitting): about 10 s for the d = 89
  radon model, and several times that for d = 175 and d = 390.
- Compiling a Stan model with model methods: about 1 minute and up to 3 GB of memory.
