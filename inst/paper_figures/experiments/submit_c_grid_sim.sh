#!/bin/bash
# Probe the partition, then submit the c-grid array excluding any hung nodes.
#
#   ./submit_c_grid_sim.sh [partition]
#
# Slurm cannot read an exclude list from a file inside an #SBATCH directive, and
# a CLI --exclude overrides the directive anyway, so the exclusion has to be
# assembled here at submit time. Probing immediately before submitting matters:
# node health is not stable, and a list from last week is worse than useless.
#
# Skip the probe (e.g. you already ran it) with:  SKIP_PROBE=1 ./submit_c_grid_sim.sh

set -uo pipefail
PART=${1:-defq}
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"
mkdir -p logs

## The probe is advisory. If it fails we still submit -- an unusable probe is no
## reason to block a run, and bad_nodes.txt may already hold a good list from
## blame_nodes.sh on the previous round.
if [ "${SKIP_PROBE:-0}" != "1" ]; then
  ./probe_nodes.sh "$PART" || echo "probe failed (rc=$?); submitting with whatever $PWD/bad_nodes.txt already holds"
fi

EX=$(tr -d '[:space:]' < bad_nodes.txt 2>/dev/null || true)

if [ -n "$EX" ]; then
  echo "submitting with --exclude=$EX"
  sbatch --partition="$PART" --exclude="$EX" run_c_grid_sim.sbatch
else
  echo "submitting with no exclusions"
  sbatch --partition="$PART" run_c_grid_sim.sbatch
fi
