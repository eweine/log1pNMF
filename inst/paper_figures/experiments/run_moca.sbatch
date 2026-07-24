#!/bin/bash
# MOCA exact-vs-approx timing: one array task per (scheme, method) combination.
# Task 0 random/approx, 1 random/exact, 2 rank1/approx, 3 rank1/exact.
#SBATCH --job-name=eva_moca
#SBATCH --output=eva_moca_%A_%a.out
#SBATCH --array=0-3
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=100G
#SBATCH --exclusive
#SBATCH --partition=defq

module load R/4.3.2
EXPDIR=/rafalab/eweine/log1pNMF/inst/paper_figures/experiments
cd "$EXPDIR"

COMBOS=("random approx" "random exact" "rank1 approx" "rank1 exact")
SPEC=${COMBOS[$SLURM_ARRAY_TASK_ID]}
echo "MOCA task $SLURM_ARRAY_TASK_ID: $SPEC on $(hostname)"
Rscript run_experiment.R moca $SPEC
