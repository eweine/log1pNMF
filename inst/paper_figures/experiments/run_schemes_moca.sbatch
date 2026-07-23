#!/bin/bash
#SBATCH --account=pi-mstephens
#SBATCH --job-name=log1p_moca
#SBATCH --output=log1p_moca_%A_%a.out
#SBATCH --error=log1p_moca_%A_%a.err
#SBATCH --time=36:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH --array=1-2

module load R/4.3.1

# make OpenMP use the full allocation (SLURM often defaults OMP_NUM_THREADS=1)
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-48}
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"

# the job scripts source the MOCA common with a RELATIVE path, so run from here
EXPDIR=/home/ericweine/log1pNMF/inst/paper_figures/experiments
cd "$EXPDIR"

case $SLURM_ARRAY_TASK_ID in
  1) SCRIPT=exact_vs_approx_moca_rank1_plus1_approx.R ;;
  2) SCRIPT=exact_vs_approx_moca_rank1_plus1_exact.R  ;;
esac

Rscript "$EXPDIR/$SCRIPT"
