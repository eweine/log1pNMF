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
# Fits run to convergence (TOL = 1e-6) or to MAXITER = 100k, whichever comes
# first. There is deliberately NO optimization-time cap, and the wall clock
# below is set very generously -- the point of this experiment is to compare
# the CONVERGED exact and approximate fits, so nothing may stop mid-descent.
# 48h is roughly 20x the worst fit observed in the (smaller) c-grid round;
# raise it further if the partition allows and any task ever hits it.
#
# Submit via the wrapper, which probes for hung nodes first and passes them
# to --exclude (a hung node still advertises itself as idle):
#
#   ./submit_approx_sim.sh                 # probes, then submits
#   ./blame_nodes.sh                       # after the round: who ate the tasks
#
# --output is an ABSOLUTE path on purpose: Slurm resolves it against the
# submission directory and opens it before the script body runs, so a
# relative path would land wherever sbatch was invoked from. The directory
# must exist before submitting (the wrapper mkdir -p's it).
#
#SBATCH --job-name=approxsim
#SBATCH --output=/rafalab/eweine/log1pNMF/inst/paper_figures/experiments/logs/approxsim_%A_%a.out
#SBATCH --array=0-269
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=defq
## no --exclude directive: submit_approx_sim.sh supplies one from the live
## probe, and a CLI --exclude would override anything written here anyway

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
