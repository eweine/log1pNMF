#!/bin/bash
#SBATCH --account=pi-mstephens
#SBATCH --job-name=log1p_optimization
#SBATCH --output=log1p_optimization_%A_%a.out
#SBATCH --error=log1p_optimization_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=16G
#SBATCH --array=1-6

module load R/4.3.1

# Force OpenMP to use the full allocation. SLURM often exports OMP_NUM_THREADS=1
# by default, which makes the package's parallel loops run single-threaded even
# though 48 cores are allocated. Setting it here (and echoing it) fixes/confirms.
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-48}
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
env | grep -iE 'OMP|MKL|BLAS|GOMP' || true

# the job scripts source("exact_vs_approx_common.R") with a RELATIVE path,
# so we must run from the experiments directory
EXPDIR=/home/ericweine/log1pNMF/inst/paper_figures/experiments
cd "$EXPDIR"

case $SLURM_ARRAY_TASK_ID in
  1) SCRIPT=exact_vs_approx_random_approx.R      ;;
  2) SCRIPT=exact_vs_approx_random_exact.R       ;;
  3) SCRIPT=exact_vs_approx_rank1_approx.R       ;;
  4) SCRIPT=exact_vs_approx_rank1_exact.R        ;;
  5) SCRIPT=exact_vs_approx_rank1_plus1_approx.R ;;
  6) SCRIPT=exact_vs_approx_rank1_plus1_exact.R  ;;
esac

Rscript "$EXPDIR/$SCRIPT"
