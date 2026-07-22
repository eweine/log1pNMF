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
