#!/bin/bash
# Find nodes that Slurm believes are healthy but which cannot actually run a
# trivial R script, and write them to bad_nodes.txt for --exclude.
#
#   ./probe_nodes.sh [partition] [timeout_seconds] [module_to_load]
#
# Why a probe at all: a hung node still advertises itself to Slurm as
# idle/mixed and will happily accept an allocation, then never start the step.
# `sinfo` cannot see that -- the node looks fine right up until your task
# disappears into it. The only reliable test is to actually run something.
#
# This uses sbatch, NOT srun, and that distinction is the whole point. srun
# creates a STEP INSIDE your current allocation, so if you are sitting on a
# compute node it can only ever reach the nodes already allocated to you --
# every other node fails instantly for reasons that say nothing about node
# health. sbatch queues an independent job and works from anywhere you can
# submit from, which is what makes this usable without a login shell.
#
# One tiny job per node, all queued at once, then a single wait. The timeout is
# therefore for the WHOLE sweep, not per node. Classification at the deadline:
#
#   healthy   job finished and printed "ok"
#   HUNG      job was RUNNING the whole time and never printed anything --
#             the node took the work and did nothing with it. Excluded.
#   pending   never got scheduled (node genuinely busy). Not excluded: a busy
#             node is fine, and excluding it only shrinks the pool.
#   failed    job ended without printing "ok" -- stderr is shown.
#
# The probe is advisory. If it cannot run, submit anyway and use blame_nodes.sh,
# which infers bad nodes from the previous round's real logs, needs no
# allocations, and cannot fail this way.

set -uo pipefail

PART=${1:-defq}
TMOUT=${2:-180}                 # seconds for the whole sweep
RMOD=${3:-R/4.3.2}              # module to load inside the probe job
OUT=bad_nodes.txt
LOG=node_probe.log

for need in sbatch squeue sinfo scancel; do
  command -v "$need" >/dev/null 2>&1 || { echo "required command not found: $need" >&2; exit 2; }
done

: > "$LOG"
D=$(mktemp -d)
trap 'rm -rf "$D"' EXIT

nodes=$(sinfo -h -p "$PART" -N -o "%N %T" \
        | awk '$2 ~ /^(idle|mixed|allocated)/ {print $1}' | sort -u)
[ -z "$nodes" ] && { echo "no schedulable nodes in partition '$PART'" | tee -a "$LOG"; exit 1; }

nnodes=$(echo "$nodes" | wc -w | tr -d ' ')
echo "probing $nnodes node(s) in '$PART' via sbatch (up to ${TMOUT}s for the sweep)" | tee -a "$LOG"

jids=()
for n in $nodes; do
  jid=$(sbatch --parsable --partition="$PART" --nodelist="$n" \
        --nodes=1 --ntasks=1 --cpus-per-task=1 --time=00:02:00 \
        --job-name=cgrid_probe --output="$D/$n.out" --error="$D/$n.err" \
        --wrap="module load $RMOD >/dev/null 2>&1; Rscript -e 'cat(\"ok\n\")'" 2>"$D/$n.sub")
  if [ -n "$jid" ]; then echo "$jid" > "$D/$n.jid"; jids+=("$jid")
  else printf 'SUBMIT-FAIL %s  %s\n' "$n" "$(head -1 "$D/$n.sub")" | tee -a "$LOG"; fi
done

[ ${#jids[@]} -eq 0 ] && { echo "no probe jobs were accepted; nothing to learn" | tee -a "$LOG"; exit 3; }

## wait for the queue to drain, or for the deadline
deadline=$(( $(date +%s) + TMOUT ))
list=$(IFS=,; echo "${jids[*]}")
while [ "$(date +%s)" -lt "$deadline" ]; do
  left=$(squeue -h -j "$list" -o '%i' 2>/dev/null | wc -l | tr -d ' ')
  [ "$left" = "0" ] && break
  sleep 5
done

hung=(); pend=(); fail=(); ok=0
for n in $nodes; do
  [ -f "$D/$n.jid" ] || continue
  jid=$(cat "$D/$n.jid")
  state=$(squeue -h -j "$jid" -o '%T' 2>/dev/null | head -1)
  if grep -q '^ok$' "$D/$n.out" 2>/dev/null; then
    printf 'OK       %s\n' "$n" | tee -a "$LOG"; ok=$((ok+1))
  elif [ "$state" = "RUNNING" ]; then
    printf 'HUNG     %s  (RUNNING for %ss, produced nothing)\n' "$n" "$TMOUT" | tee -a "$LOG"
    hung+=("$n")
  elif [ -n "$state" ]; then
    printf 'PENDING  %s  (%s -- node busy, not excluded)\n' "$n" "$state" | tee -a "$LOG"
    pend+=("$n")
  else
    printf 'FAILED   %s  %s\n' "$n" "$(head -1 "$D/$n.err" 2>/dev/null | cut -c1-100)" | tee -a "$LOG"
    fail+=("$n")
  fi
done

scancel --name=cgrid_probe >/dev/null 2>&1     # clean up anything still queued

printf '\n%d healthy, %d hung, %d pending, %d failed\n' \
  "$ok" "${#hung[@]}" "${#pend[@]}" "${#fail[@]}" | tee -a "$LOG"

## A probe where nothing at all worked is a broken probe, not a dead cluster.
## Writing an exclude list from that would be worse than having no probe.
if [ "$ok" -eq 0 ] && [ ${#hung[@]} -eq 0 ]; then
  cat <<MSG | tee -a "$LOG"

NOT writing $OUT: nothing succeeded and nothing hung, so this run says nothing
about node health. Check the messages above -- likely the partition rejected an
option, or '$RMOD' is not the right module (pass it as the third argument).
Submit without exclusions and use blame_nodes.sh on the resulting logs.
MSG
  exit 3
fi

if [ ${#hung[@]} -gt 0 ]; then
  (IFS=,; echo "${hung[*]}") > "$OUT"
  echo "wrote $OUT: $(cat "$OUT")" | tee -a "$LOG"
else
  : > "$OUT"
  echo "no hung nodes; $OUT is empty" | tee -a "$LOG"
fi
