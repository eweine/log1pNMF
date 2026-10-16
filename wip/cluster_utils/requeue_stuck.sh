#!/bin/bash
# Cancel the array tasks of a running job that are stuck on a given node, and
# resubmit exactly those task ids with that node excluded. Everything else in
# the array is left alone.
#
#   ./requeue_stuck.sh <array_job_id> <node> [partition]
#   ./requeue_stuck.sh 123456 node13
#
# Safe to do mid-run: each task regenerates its own data from (seed, c_true) and
# writes its own .rds, so tasks share no state and a resubmitted task simply
# redoes work that never completed. Nothing that already finished is touched.
#
# The node is also appended to bad_nodes.txt, so submit_c_grid_sim.sh keeps
# excluding it on later rounds.
#
# Prints the plan and asks before doing anything; pass -y to skip the prompt.

set -uo pipefail

YES=0
[ "${1:-}" = "-y" ] && { YES=1; shift; }

JOB=${1:-}; NODE=${2:-}; PART=${3:-defq}
if [ -z "$JOB" ] || [ -z "$NODE" ]; then
  echo "usage: $0 [-y] <array_job_id> <node> [partition]" >&2; exit 1
fi

for need in squeue scancel sbatch; do
  command -v "$need" >/dev/null 2>&1 || { echo "not found: $need" >&2; exit 2; }
done

# Array task ids currently RUNNING on that node. %K is the array index, %N the
# node. Only running tasks matter -- pending ones have no node yet and will be
# steered away by the --exclude on the resubmission anyway.
## plain read loop rather than mapfile, which is bash 4+ only
TASKS=()
while IFS= read -r t; do
  [ -n "$t" ] && TASKS+=("$t")
done < <(squeue -h -j "$JOB" -t RUNNING -o '%K %N' 2>/dev/null \
         | awk -v n="$NODE" '$2 == n {print $1}' | sort -n)

if [ ${#TASKS[@]} -eq 0 ]; then
  echo "no RUNNING tasks of job $JOB on $NODE -- nothing to do"
  squeue -h -j "$JOB" -o '%T' 2>/dev/null | sort | uniq -c | sed 's/^/  /'
  exit 0
fi

LIST=$(IFS=,; echo "${TASKS[*]}")
echo "job $JOB has ${#TASKS[@]} task(s) running on $NODE:"
echo "  $LIST"

# Accumulate the exclusion so later submissions keep avoiding this node.
EX=$(tr -d '[:space:]' < bad_nodes.txt 2>/dev/null || true)
case ",$EX," in
  *",$NODE,"*) ;;                                   # already listed
  *) EX=${EX:+$EX,}$NODE ;;
esac

echo
echo "plan:"
echo "  1. scancel ${JOB}_[$LIST]"
echo "  2. record '$NODE' in bad_nodes.txt  (list becomes: $EX)"
echo "  3. sbatch --array=$LIST --exclude=$EX run_c_grid_sim.sbatch"
echo

if [ "$YES" != "1" ]; then
  read -r -p "proceed? [y/N] " ans
  case "$ans" in [yY]*) ;; *) echo "aborted"; exit 0 ;; esac
fi

## belt and braces: an empty LIST here would become `scancel JOB_[]`, which is
## not a no-op, and `sbatch --array=` with no value
[ -n "$LIST" ] || { echo "internal error: empty task list, refusing to act" >&2; exit 4; }

scancel "${JOB}_[$LIST]" || { echo "scancel failed" >&2; exit 3; }
echo "cancelled ${#TASKS[@]} task(s)"

echo "$EX" > bad_nodes.txt

sbatch --partition="$PART" --array="$LIST" --exclude="$EX" run_c_grid_sim.sbatch
