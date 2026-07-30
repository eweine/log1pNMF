#!/bin/bash
# c-grid simulation: one array task per (seed, c_true, c_fit) cell.
#
#   10 seeds x 8 c_true x 8 c_fit = 640 tasks, indexed 0..639 with c_fit
#   varying fastest (see task_spec() in c_grid_sim_common.R).
#
# Each task is single-threaded on purpose: the 500 x 500 problem is small enough
# that one core beats any thread-level parallelism, and 640 one-core tasks pack
# into the cluster far better than 640 multi-core ones. Each task regenerates
# its own data and writes its own .rds, so tasks share no state and a failure
# costs one cell.
#
# Fits run to convergence (TOL in c_grid_sim_common.R), not to a fixed iteration
# count. Cost is very uneven across the grid: c_fit >= 1 converges in a few
# hundred iterations / tens of seconds, while the smallest c_fit values are far
# slower. MAX_TIME (2 h of optimization time per fit) is the backstop, so the
# wall clock below leaves an hour of headroom for data generation, scoring and
# a slower node.
#
# Submit:   sbatch run_c_grid_sim.sbatch
# Collect:  Rscript collect_c_grid_sim.R
#
#SBATCH --job-name=cgrid
#SBATCH --output=logs/cgrid_%A_%a.out
#SBATCH --array=0-639
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=defq

module load R/4.3.2

EXPDIR=/rafalab/eweine/log1pNMF/inst/paper_figures/experiments
cd "$EXPDIR"
mkdir -p logs

## keep the math libraries single-threaded too, so one task really is one core
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "cgrid task $SLURM_ARRAY_TASK_ID on $(hostname) at $(date)"
Rscript run_c_grid_sim.R "$SLURM_ARRAY_TASK_ID"
