#!/bin/bash
# c-grid rerun: one array task per entry of c_grid_rerun_todo.rds.
#
# Workflow:
#   mkdir -p /rafalab/eweine/log1pNMF/inst/paper_figures/experiments/logs
#   Rscript plan_c_grid_rerun.R                  # builds the worklist, prints N
#   sbatch --array=0-<N-1> run_c_grid_rerun.sbatch
#   Rscript collect_c_grid_sim.R                 # unchanged; reads the same files
#
# There is deliberately NO --array directive here: the range depends on how many
# cells the planner found, so it has to come from the command line. Submitting
# without it produces a single task with SLURM_ARRAY_TASK_ID unset, and the
# worker stops with an explicit error rather than doing something arbitrary.
#
# node13 is excluded: in round one it accepted jobs and then hung, so tasks
# scheduled there neither completed nor failed cleanly.
#
# Cost is dominated by the same corner as before (small c_fit), and warm-started
# cells start from an already-good point but still get a full MAX_TIME budget,
# so the 3 h wall clock is unchanged.
#
#SBATCH --job-name=cgrid_rerun
#SBATCH --output=/rafalab/eweine/log1pNMF/inst/paper_figures/experiments/logs/cgrid_rerun_%A_%a.out
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=defq
#SBATCH --exclude=node13

module load R/4.3.2

EXPDIR=/rafalab/eweine/log1pNMF/inst/paper_figures/experiments
cd "$EXPDIR"
mkdir -p logs

## keep the math libraries single-threaded too, so one task really is one core
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "cgrid rerun task $SLURM_ARRAY_TASK_ID on $(hostname) at $(date)"
Rscript run_c_grid_rerun.R "$SLURM_ARRAY_TASK_ID"
