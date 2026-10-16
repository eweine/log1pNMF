#!/bin/bash
# Which nodes swallowed array tasks? Tally the hostname of every task that
# started but never finished.
#
#   ./blame_nodes.sh [logdir]        # default: logs
#
# A completed task always prints "Wrote <path>" as its last action. A task that
# started (the sbatch echo names the host) but has no "Wrote" line either died
# or hung. Grouping those by host, and comparing against the same tally for
# successful tasks, points straight at a bad node -- one that fails nearly
# everything it is given while healthy nodes fail nothing.
#
# This is stronger evidence than probe_nodes.sh, because it is the real
# workload rather than a synthetic test. Run it after a round, feed the result
# into bad_nodes.txt, and the next submission avoids those nodes.

set -uo pipefail
LOGDIR=${1:-logs}

if [ ! -d "$LOGDIR" ]; then echo "no log directory: $LOGDIR" >&2; exit 1; fi

shopt -s nullglob
logs=("$LOGDIR"/cgrid*_*.out)
if [ ${#logs[@]} -eq 0 ]; then echo "no cgrid logs in $LOGDIR" >&2; exit 1; fi

tmp=$(mktemp)
for f in "${logs[@]}"; do
  host=$(sed -n 's/^cgrid[a-z_ ]*task [0-9]* on \([^ ]*\) at .*/\1/p' "$f" | head -1)
  [ -z "$host" ] && continue                      # never even started
  if grep -q '^Wrote ' "$f"; then
    printf '%s\tok\n'   "$host" >> "$tmp"
  else
    printf '%s\tfail\n' "$host" >> "$tmp"
  fi
done

echo "task outcomes by node (${#logs[@]} logs):"
printf '%-20s %8s %8s %8s\n' NODE OK FAIL "FAIL%"
awk -F'\t' '
  {tot[$1]++; if ($2=="fail") bad[$1]++}
  END {for (n in tot)
         printf "%-20s %8d %8d %7.1f%%\n", n, tot[n]-bad[n], bad[n]+0, 100*(bad[n]+0)/tot[n]}
' "$tmp" | sort -k4 -gr

echo
echo "nodes failing everything they were given (candidates for --exclude):"
awk -F'\t' '
  {tot[$1]++; if ($2=="fail") bad[$1]++}
  END {for (n in tot) if (bad[n]+0 == tot[n]) print n}
' "$tmp" | sort | paste -sd, -

rm -f "$tmp"
