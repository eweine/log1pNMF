#!/bin/bash
# Find nodes that Slurm believes are healthy but which cannot actually run a
# trivial R script, and write them to bad_nodes.txt for --exclude.
#
#   ./probe_nodes.sh [partition] [timeout_seconds]
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
TMOUT=${2:-120}
IMM=30                      # seconds srun will wait for an allocation
OUT=bad_nodes.txt
LOG=node_probe.log

: > "$LOG"
bad=(); busy=(); ok=0

# Only probe nodes that claim to be usable. Drop the ones Slurm already knows
# are down/drained -- Slurm will not schedule to those anyway.
nodes=$(sinfo -h -p "$PART" -N -o "%N %T" \
        | awk '$2 ~ /^(idle|mixed|allocated)/ {print $1}' | sort -u)

if [ -z "$nodes" ]; then
  echo "no schedulable nodes found in partition '$PART'" | tee -a "$LOG"
  exit 1
fi

echo "probing $(echo "$nodes" | wc -w) node(s) in partition '$PART' (timeout ${TMOUT}s)" | tee -a "$LOG"

for n in $nodes; do
  timeout "$TMOUT" srun --partition="$PART" --nodelist="$n" --nodes=1 --ntasks=1 \
      --cpus-per-task=1 --mem=1G --time=00:02:00 --immediate=$IMM \
      --job-name=cgrid_probe \
      Rscript -e 'cat("ok\n")' >/dev/null 2>&1
  rc=$?
  case $rc in
    0)   printf 'OK    %s\n' "$n"                  | tee -a "$LOG"; ok=$((ok+1)) ;;
    124) printf 'HUNG  %s  (timed out)\n' "$n"     | tee -a "$LOG"; bad+=("$n") ;;
    *)   printf 'BUSY? %s  (srun rc=%d)\n' "$n" $rc | tee -a "$LOG"; busy+=("$n") ;;
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
