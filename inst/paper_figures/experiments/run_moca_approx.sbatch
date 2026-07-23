#!/bin/bash
#SBATCH --job-name=log1p_moca_approx
#SBATCH --output=log1p_moca_approx_%j.out
#SBATCH --error=log1p_moca_approx_%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G

module load R/4.3.2
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-32}
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"

EXPDIR=/rafalab/eweine/log1pNMF/inst/paper_figures/experiments
cd "$EXPDIR"
Rscript "$EXPDIR/exact_vs_approx_moca_rank1_approx.R"
