#!/bin/bash
# Pure-R experiments (scripts 01-05, supplement S01-S22, S98): one node, many cores.
# Each script parallelizes its replications over BO_CORES forked workers; results do not
# depend on the number of cores (every replication has its own seed, see analysis/R/common.R).
#SBATCH --job-name=bo-pure-r
#SBATCH --cpus-per-task=64
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=hpc/logs/pure_r_%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
export BO_SCALE=1 BO_QUICK=false BO_CORES="${SLURM_CPUS_PER_TASK}"
Rscript analysis/run_all.R          # without --stan: also runs S99_assemble.R (provisional)
