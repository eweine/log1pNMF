#!/bin/bash
# Exact-vs-approx simulation: one array task per (seed, cc, method) fit.
#
#   10 seeds x 9 cc x 3 methods = 270 tasks, indexed 0..269 with method
#   varying fastest, then cc (see task_spec() in approx_sim_common.R).
#
# Each task is single-threaded on purpose: the 500 x 500 problem is small
# enough that one core beats thread-level parallelism, and 270 one-core tasks
# pack into the cluster far better than multi-core ones. Each task regenerates
# its own data and writes its own .rds, so tasks share no state and a failure
# costs one fit.
#
# Fits run to convergence (TOL = 1e-6) or to MAXITER = 250k, whichever comes
# first. There is deliberately NO optimization-time cap, and the wall clock
# below is a generous backstop -- the point of this experiment is to compare
# the CONVERGED exact and approximate fits, so nothing may stop mid-descent.
#
# Submit directly (no node-probing wrapper for this round; any task lost to a
# bad node is found afterward and resubmitted):
#
#   sbatch run_approx_sim.sbatch
#   Rscript check_approx_sim.R             # afterward: what is missing/unconverged
#
# --output is an ABSOLUTE path on purpose: Slurm resolves it against the
# submission directory and opens it before the script body runs, so a
# relative path would land wherever sbatch was invoked from. The logs/
# directory must exist BEFORE submitting (mkdir -p logs); the mkdir in the
# body runs too late for --output.
#
#SBATCH --job-name=approxsim
#SBATCH --output=/rafalab/eweine/log1pNMF/inst/paper_figures/experiments/logs/approxsim_%A_%a.out
#SBATCH --array=0-269
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=defq

module load R/4.3.2

EXPDIR=/rafalab/eweine/log1pNMF/inst/paper_figures/experiments
cd "$EXPDIR"
mkdir -p logs      # belt and braces; must already exist for --output to open

## keep the math libraries single-threaded too, so one task really is one core
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "approxsim task $SLURM_ARRAY_TASK_ID on $(hostname) at $(date)"
Rscript run_approx_sim.R "$SLURM_ARRAY_TASK_ID"
