#!/bin/bash
# c-grid simulation: one array task per (seed, c_true, c_fit, init) cell.
#
#   10 seeds x 8 c_true x 8 c_fit x 2 inits = 1280 tasks, indexed 0..1279 with
#   init varying fastest, then c_fit (see task_spec() in c_grid_sim_common.R).
#   The two initializations of a cell are adjacent task ids on purpose.
#
# Each task is single-threaded on purpose: the 300 x 300 problem is small enough
# that one core beats any thread-level parallelism, and 1280 one-core tasks pack
# into the cluster far better than 1280 multi-core ones. Each task regenerates
# its own data and writes its own .rds, so tasks share no state and a failure
# costs one cell.
#
# Fits run to convergence (TOL) or to MAXITER = 100k, whichever comes first.
# There is no optimization-time cap: an earlier round used one and left cells
# stopped mid-descent, which is indistinguishable in the output from a cell that
# genuinely plateaued. The wall clock below is the only backstop.
#
# Submit via the wrapper, which probes for hung nodes first and passes them to
# --exclude. A hung node still advertises itself as idle and will accept an
# allocation, so sinfo alone will not protect you:
#
#   ./submit_c_grid_sim.sh                 # probes, then submits
#   ./blame_nodes.sh                       # after the round: who ate the tasks
#   Rscript plan_c_grid_rerun.R            # what still needs redoing
#   Rscript collect_c_grid_sim.R
#
# --output is an ABSOLUTE path on purpose. Slurm resolves it against the
# submission directory, not against the EXPDIR this script cd's into, and it
# opens the file before the script body runs -- so a relative "logs/..." would
# land wherever sbatch happened to be invoked from, and the mkdir below would
# come too late to create it. The directory must exist before submitting.
#
#SBATCH --job-name=cgrid
#SBATCH --output=/rafalab/eweine/log1pNMF/inst/paper_figures/experiments/logs/cgrid_%A_%a.out
#SBATCH --array=0-639
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=defq
## no --exclude directive: submit_c_grid_sim.sh supplies one from the live
## probe, and a CLI --exclude would override anything written here anyway

module load R/4.3.2

EXPDIR=/rafalab/eweine/log1pNMF/inst/paper_figures/experiments
cd "$EXPDIR"
mkdir -p logs      # belt and braces; must already exist for --output to open

## keep the math libraries single-threaded too, so one task really is one core
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "cgrid task $SLURM_ARRAY_TASK_ID on $(hostname) at $(date)"
Rscript run_c_grid_sim.R "$SLURM_ARRAY_TASK_ID"
