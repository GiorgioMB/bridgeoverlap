#!/bin/bash
# Stan benchmark: one array task per line of hpc/stan_tasks.txt (comment lines skipped).
# Every task compiles its Stan programs (with model methods) in its own build directory,
# BO_STAN_BUILD_DIR: a new R session always rebuilds the executable, so tasks sharing one
# directory would overwrite each other's binaries. Compiling takes about a minute per
# program and up to 3 GB of memory.
#SBATCH --job-name=bo-stan
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=48:00:00
#SBATCH --output=hpc/logs/stan_%A_%a.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
line=$(grep -v '^#' hpc/stan_tasks.txt | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")
script=$(echo "$line" | cut -f1); configs=$(echo "$line" | cut -f2); range=$(echo "$line" | cut -f3)
export BO_SCALE=1 BO_QUICK=false BO_STAN_CORES=1 BO_STAN_CONFIGS="${configs}" BO_REP_RANGE="${range}"
export TMPDIR="${TMPDIR:-/tmp}/bo_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"; mkdir -p "$TMPDIR"
export BO_STAN_BUILD_DIR="$TMPDIR/stan"; mkdir -p "$BO_STAN_BUILD_DIR"
echo "task ${SLURM_ARRAY_TASK_ID}: ${script} configs=${configs:-all} reps=${range}"
Rscript "analysis/${script}"
