#!/bin/bash
# Submit the full reproduction with SLURM dependencies. Run from the repository root,
# after `git checkout <tag>` and `R CMD INSTALL .` (see hpc/README.md).
set -euo pipefail
mkdir -p hpc/logs
# toolchain check: compile one small Stan program with model methods (the array tasks
# compile their own copies; see stan_array.slurm). Needs about 3 GB of memory.
Rscript hpc/check_setup.R || { echo "setup check failed: see hpc/README.md"; exit 1; }
n=$(grep -vc '^#' hpc/stan_tasks.txt)
j1=$(sbatch --parsable hpc/pure_r.slurm)
j2=$(sbatch --parsable --array=0-$((n - 1)) hpc/stan_array.slurm)
j3=$(sbatch --parsable --dependency=afterok:${j1}:${j2} hpc/finish.slurm)
echo "submitted: pure R ${j1}, Stan array ${j2} (${n} tasks), finish ${j3}"
