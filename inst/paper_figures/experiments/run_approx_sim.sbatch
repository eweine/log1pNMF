#!/bin/bash
# Exact-vs-approx simulation, extended round: one array task per
# (seed, cc, method) fit.
#
#   30 seeds x 3 c_true x (9 finite c_fit x 3 methods + c_fit = Inf, exact
#   only) = 2520 tasks, indexed 0..2519 with method varying fastest, then
#   c_fit, then c_true (see the SPEC table in approx_sim_common.R).
#
# The generator and finite grid are IDENTICAL to round one, so the 810
# round-one outputs in OUT_DIR are seeds 1-10 of this grid; the worker
# skips any task whose output already exists (a skipped task costs ~2s of
# R startup). Only the 20 new seeds and the c = Inf topic-model fits run.
#
# The grid exceeds the Slurm MaxArraySize default of 1001, so submit in
# slices; the worker adds APPROX_SIM_TASK_OFFSET to its array id. Full
# grid (ids 0-3359; completed slices are skipped in seconds):
#
#   sbatch run_approx_sim.sbatch
#   sbatch --export=ALL,APPROX_SIM_TASK_OFFSET=990  run_approx_sim.sbatch
#   sbatch --export=ALL,APPROX_SIM_TASK_OFFSET=1980 --array=0-539 run_approx_sim.sbatch
#   sbatch --export=ALL,APPROX_SIM_TASK_OFFSET=2520 --array=0-839 run_approx_sim.sbatch
#
# The last slice alone covers the EXTRA c_true = 1e-3 seeds 31-60 (ids
# 2520-3359), which is all that needs submitting if the first 2520 are
# already done. Afterward:
#
#   Rscript check_approx_sim.R             # what is missing/unconverged
#
# Each task is single-threaded on purpose: the 500 x 500 problem is small
# enough that one core beats thread-level parallelism, and one-core tasks
# pack into the cluster far better than multi-core ones. Each task regenerates
# its own data and writes its own .rds, so tasks share no state and a failure
# costs one fit.
#
# Fits run to convergence (TOL = 1e-6) or to MAXITER = 2M, whichever comes
# first (the extended run's slowest fit spent ~10h reaching 1M iterations,
# so 2M projects to ~20h against the 24h wall clock below -- close, and
# any TIMEOUT is caught by the checker).
# There is deliberately NO optimization-time cap -- the point of this
# experiment is to compare the CONVERGED fits, so nothing may stop
# mid-descent.
#
# No node-probing wrapper for this round; any task lost to a bad node is
# found afterward and resubmitted.
#
# --output is an ABSOLUTE path on purpose: Slurm resolves it against the
# submission directory and opens it before the script body runs, so a
# relative path would land wherever sbatch was invoked from. The logs/
# directory must exist BEFORE submitting (mkdir -p logs); the mkdir in the
# body runs too late for --output.
#
#SBATCH --job-name=approxsim
#SBATCH --output=/rafalab/eweine/log1pNMF/inst/paper_figures/experiments/logs/approxsim_%A_%a.out
#SBATCH --array=0-989
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=defq

module load R/4.3.2

EXPDIR=/rafalab/eweine/log1pNMF/inst/paper_figures/experiments
cd "$EXPDIR"
mkdir -p logs      # belt and braces; must already exist for --output to open

## keep the math libraries in step with the Slurm allocation (1 core unless
## the submission grants more via --cpus-per-task)
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

echo "approxsim task $SLURM_ARRAY_TASK_ID on $(hostname) at $(date)"
Rscript run_approx_sim.R "$SLURM_ARRAY_TASK_ID"
