#!/bin/bash
# After the pure-R job and all Stan tasks: summaries and assembly of generated/.
#SBATCH --job-name=bo-finish
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=hpc/logs/finish_%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
export BO_SCALE=1 BO_QUICK=false
export TMPDIR="${TMPDIR:-/tmp}/bo_finish_${SLURM_JOB_ID}"; mkdir -p "$TMPDIR"
export BO_STAN_BUILD_DIR="$TMPDIR/stan"; mkdir -p "$BO_STAN_BUILD_DIR"
Rscript analysis/hpc_combine.R                  # combine the per-range outputs of scripts 10-12
Rscript analysis/31_stan_benchmark_summary.R
Rscript analysis/supplement/S99_assemble.R
grep -q "Final run (BO_SCALE = 1, no BO_QUICK): yes" generated/MANIFEST.txt || { echo "not a final run"; exit 1; }
