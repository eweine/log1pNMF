#!/bin/bash
#SBATCH --account=pi-mstephens
#SBATCH --job-name=log1p_moca_probe
#SBATCH --output=log1p_moca_probe_%j.out
#SBATCH --error=log1p_moca_probe_%j.err
#SBATCH --time=04:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G

# Short probe: times a few approx and a few exact iterations on MOCA to estimate
# per-iteration cost before committing to a full run. Should finish well inside
# the walltime; lower --time if you want it to schedule even faster.

module load R/4.3.1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-48}
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"

EXPDIR=/home/ericweine/log1pNMF/inst/paper_figures/experiments
cd "$EXPDIR"
Rscript "$EXPDIR/probe_moca_timing.R"
