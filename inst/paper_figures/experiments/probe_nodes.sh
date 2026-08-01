#!/bin/bash
# Find nodes that Slurm believes are healthy but which cannot actually run a
# trivial R script, and write them to bad_nodes.txt for --exclude.
#
#   ./probe_nodes.sh [partition] [timeout_seconds] [parallelism]
#
# Why a probe at all: a hung node still advertises itself to Slurm as
# idle/mixed and will happily accept an allocation, then never start the step.
# `sinfo` cannot see that -- the node looks fine right up until your task
# disappears into it. The only reliable test is to actually run something.
#
# Each node gets one `srun` of Rscript, wrapped in `timeout`. Exit codes:
#
#   0          healthy.
#   124        timeout fired -- allocation granted, step never returned. THE
#              HANG. Excluded.
#   other      srun refused. Reported WITH ITS STDERR and not excluded: a busy
#              or drained node is fine and excluding it only shrinks the pool.
#
# The timeout is PER NODE; probes run in parallel (default 8 at a time) so wall
# time is roughly ceil(nodes/parallelism) * slowest_probe, not the sum.
#
# IMPORTANT: this probe is advisory. If it cannot run, submit anyway and rely on
# blame_nodes.sh, which infers bad nodes from the previous round's real logs and
# needs no allocations at all. That is the more dependable of the two.

set -uo pipefail

PART=${1:-defq}
TMOUT=${2:-120}             # PER NODE
PAR=${3:-8}                 # how many nodes to probe at once
OUT=bad_nodes.txt
LOG=node_probe.log

for need in timeout srun sinfo; do
  command -v "$need" >/dev/null 2>&1 || { echo "required command not found: $need" >&2; exit 2; }
done

: > "$LOG"
RCDIR=$(mktemp -d)
trap 'rm -rf "$RCDIR"' EXIT

# Running inside an existing allocation breaks this completely: a nested srun is
# a STEP of that allocation, so --nodelist for any node not in it fails at once.
# That produces exactly the "every node is busy" pattern, which is meaningless.
if [ -n "${SLURM_JOB_ID:-}" ]; then
  cat <<MSG | tee -a "$LOG"
WARNING: SLURM_JOB_ID=${SLURM_JOB_ID} is set, so you are inside an allocation.
         A nested srun runs as a STEP of that allocation and cannot reach nodes
         outside it, so every probe will fail for reasons that say nothing about
         node health. Run this from a login node instead.
MSG
fi

nodes=$(sinfo -h -p "$PART" -N -o "%N %T" \
        | awk '$2 ~ /^(idle|mixed|allocated)/ {print $1}' | sort -u)

if [ -z "$nodes" ]; then
  echo "no schedulable nodes found in partition '$PART'" | tee -a "$LOG"
  exit 1
fi

nnodes=$(echo "$nodes" | wc -w | tr -d ' ')
echo "probing $nnodes node(s) in partition '$PART' (${TMOUT}s per node, $PAR at a time)" | tee -a "$LOG"

probe_one() {
  local n=$1
  timeout "$TMOUT" srun --partition="$PART" --nodelist="$n" --nodes=1 --ntasks=1 \
      --cpus-per-task=1 --time=00:02:00 --job-name=cgrid_probe \
      Rscript -e 'cat("ok\n")' >"$RCDIR/$n.out" 2>"$RCDIR/$n.err"
  echo $? > "$RCDIR/$n.rc"
}

for n in $nodes; do
  while [ "$(jobs -rp | wc -l)" -ge "$PAR" ]; do sleep 1; done
  probe_one "$n" &
done
wait

bad=(); other=(); ok=0
for n in $nodes; do
  rc=$(cat "$RCDIR/$n.rc" 2>/dev/null || echo 1)
  err=$(head -1 "$RCDIR/$n.err" 2>/dev/null | cut -c1-100)
  case $rc in
    0)   printf 'OK    %s\n' "$n"                            | tee -a "$LOG"; ok=$((ok+1)) ;;
    124) printf 'HUNG  %s  (timed out after %ss)\n' "$n" "$TMOUT" | tee -a "$LOG"; bad+=("$n") ;;
    *)   printf 'FAIL  %s  (rc=%s) %s\n' "$n" "$rc" "$err"   | tee -a "$LOG"; other+=("$n") ;;
  esac
done

printf '\n%d healthy, %d hung, %d refused\n' "$ok" "${#bad[@]}" "${#other[@]}" | tee -a "$LOG"

# If nothing at all succeeded, the probe is broken, not the cluster. Writing an
# exclude list from that would be actively harmful, so refuse and say why.
if [ "$ok" -eq 0 ]; then
  cat <<MSG | tee -a "$LOG"

NOT writing $OUT: every probe failed, which almost certainly means the probe
itself is wrong rather than that all $nnodes nodes are sick. Look at the rc and
stderr above. Common causes:
  - run from inside an allocation (see the SLURM_JOB_ID warning above)
  - the partition rejects one of the requested options
  - Rscript not on PATH on the compute nodes (module load R first)
Submit without exclusions and use blame_nodes.sh on the resulting logs instead.
MSG
  exit 3
fi

if [ ${#bad[@]} -gt 0 ]; then
  (IFS=,; echo "${bad[*]}") > "$OUT"
  echo "wrote $OUT: $(cat "$OUT")" | tee -a "$LOG"
else
  : > "$OUT"
  echo "no hung nodes; $OUT is empty" | tee -a "$LOG"
fi
