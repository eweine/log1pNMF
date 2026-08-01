#!/bin/bash
# Find nodes that Slurm believes are healthy but which cannot actually run a
# trivial R script, and write them to bad_nodes.txt for --exclude.
#
#   ./probe_nodes.sh [partition] [timeout_seconds] [parallelism]
#
# The timeout is PER NODE, not for the whole sweep. Probes run concurrently
# (default 16 at a time), so wall time is roughly
# ceil(nodes / parallelism) * slowest_probe rather than the sum over nodes --
# which matters because a hung node costs the full timeout while a healthy one
# answers in seconds.
#
# Why this works: a hung node still advertises itself to Slurm as idle/mixed and
# will happily accept an allocation, then never start the step. So `sinfo` alone
# cannot tell you about it -- the node looks fine right up until your task
# disappears into it. The only reliable test is to actually run something.
#
# Each node gets one `srun` of `Rscript -e cat("ok")`, wrapped in `timeout`. The
# exit code separates the two failure modes that matter:
#
#   124        timeout fired -- the allocation was granted but the step never
#              returned. THIS is the hang. Excluded.
#   other != 0 srun refused quickly (node busy, down, reserved, wrong
#              constraints). Reported but NOT excluded -- a busy node is fine,
#              and excluding it would just shrink the pool for no reason.
#   0          healthy.
#
# --immediate=NN keeps a merely-busy node from blocking the whole probe: srun
# gives up waiting for resources and returns non-zero fast, which lands in the
# "other" bucket rather than looking like a hang.
#
# Re-run this right before each submission; node health is not stable.

set -uo pipefail

PART=${1:-defq}
TMOUT=${2:-120}             # PER NODE
PAR=${3:-16}                # how many nodes to probe at once
IMM=30                      # seconds srun will wait for an allocation
OUT=bad_nodes.txt
LOG=node_probe.log

## the whole method rests on distinguishing "timed out" (exit 124) from other
## failures, so bail early rather than silently misclassify every node
for need in timeout srun sinfo; do
  command -v "$need" >/dev/null 2>&1 || { echo "required command not found: $need" >&2; exit 2; }
done

: > "$LOG"
bad=(); busy=(); ok=0
RCDIR=$(mktemp -d)
trap 'rm -rf "$RCDIR"' EXIT

# Only probe nodes that claim to be usable. Drop the ones Slurm already knows
# are down/drained -- Slurm will not schedule to those anyway.
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
      --cpus-per-task=1 --mem=1G --time=00:02:00 --immediate=$IMM \
      --job-name=cgrid_probe \
      Rscript -e 'cat("ok\n")' >/dev/null 2>&1
  echo $? > "$RCDIR/$n.rc"
}

for n in $nodes; do
  while [ "$(jobs -rp | wc -l)" -ge "$PAR" ]; do sleep 1; done
  probe_one "$n" &
done
wait

for n in $nodes; do
  rc=$(cat "$RCDIR/$n.rc" 2>/dev/null || echo 1)
  case $rc in
    0)   printf 'OK    %s\n' "$n"                   | tee -a "$LOG"; ok=$((ok+1)) ;;
    124) printf 'HUNG  %s  (timed out)\n' "$n"      | tee -a "$LOG"; bad+=("$n") ;;
    *)   printf 'BUSY? %s  (srun rc=%s)\n' "$n" "$rc" | tee -a "$LOG"; busy+=("$n") ;;
  esac
done

printf '\n%d healthy, %d hung, %d unavailable\n' "$ok" "${#bad[@]}" "${#busy[@]}" | tee -a "$LOG"

if [ ${#bad[@]} -gt 0 ]; then
  printf '%s\n' "$(IFS=,; echo "${bad[*]}")" > "$OUT"
  echo "wrote $OUT: $(cat "$OUT")" | tee -a "$LOG"
else
  : > "$OUT"
  echo "no hung nodes; $OUT is empty" | tee -a "$LOG"
fi
