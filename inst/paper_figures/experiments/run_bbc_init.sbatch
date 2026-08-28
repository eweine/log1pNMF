#!/bin/bash
# BBC initialization-stability experiment: one array task per fit.
#
#   93 tasks = 3 c values x (30 random-init seeds + 1 rank-1 init),
#   indexed 0..92; see the SPEC table in bbc_init_common.R.
#
# The c = 0.001 tasks (ids 0-30) get 8 cores each; everything else
# (ids 31-92: c = 1 and c = Inf) gets 1 core. Slurm cannot vary
# --cpus-per-task within one array, so submit TWO arrays:
#
#   sbatch --array=0-30  --cpus-per-task=8 run_bbc_init.sbatch
#   sbatch --array=31-92 --cpus-per-task=1 run_bbc_init.sbatch
#
# The worker reads SLURM_CPUS_PER_TASK for its thread count, so the
# command-line --cpus-per-task is the ONLY place the core count is set.
# The worker skips any task whose output already exists (~2s), so
# resubmitting after failures or cancellations is safe and cheap.
#
# Fits run to convergence (TOL = 1e-6) or MAXITER = 1M (finite c; 100k for
# the c = Inf SCD fits), whichever comes first. There is deliberately no
# optimization-time cap: the point is to compare CONVERGED fits from
# different starting points, so nothing may stop mid-descent. The 96h wall
# clock below is the backstop; if the partition rejects it, lower it and
# any TIMEOUT task is found by check_bbc_init.R and resubmitted.
#
# --output is an ABSOLUTE path on purpose (Slurm opens it before the body
# runs); the logs/ directory must exist BEFORE submitting (mkdir -p logs).
#
#SBATCH --job-name=bbcinit
#SBATCH --output=/rafalab/eweine/log1pNMF/inst/paper_figures/experiments/logs/bbcinit_%A_%a.out
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=defq

module load R/4.3.2

EXPDIR=/rafalab/eweine/log1pNMF/inst/paper_figures/experiments
cd "$EXPDIR"
mkdir -p logs      # belt and braces; must already exist for --output to open

## keep the math libraries in step with the Slurm allocation
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

Rscript run_bbc_init.R
