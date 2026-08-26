#!/bin/bash
# Fire-and-forget submitter for the approx-sim array: submits each slice of
# the 2520-task grid, retrying every 5 minutes whenever the QOS per-user
# job-submit limit rejects it (the earlier slices drain as their skipped /
# finished tasks leave the queue, freeing room). Run it ONCE on the login
# node and walk away:
#
#   nohup bash submit_approx_sim.sh 990 1980 >> submit_approx_sim.log 2>&1 &
#
# The arguments are the task offsets of the slices STILL TO SUBMIT
# (no arguments = all three: 0 990 1980). If a slice is already in the
# queue, do NOT list it -- a duplicate submission would run every one of
# its unfinished tasks twice in parallel.
#
# nohup survives logout; the loop gives up on a slice after 48h of
# rejections (something else is then wrong -- read the log). Afterward:
#
#   Rscript check_approx_sim.R      # coverage over the full grid

cd "$(dirname "$0")"

TOTAL=4480     # keep in sync with N_TASKS in approx_sim_common.R
OFFSETS="${@:-0 990 1980 2970 3960}"
for OFF in $OFFSETS; do
  N=$((TOTAL - OFF)); [ $N -gt 990 ] && N=990
  RANGE=0-$((N - 1))
  echo "$(date '+%F %T')  submitting slice: offset $OFF, array $RANGE"
  tries=0
  until sbatch --export=ALL,APPROX_SIM_TASK_OFFSET=$OFF --array=$RANGE \
               run_approx_sim.sbatch; do
    tries=$((tries + 1))
    if [ $tries -ge 576 ]; then      # 576 x 300s = 48h
      echo "$(date '+%F %T')  GIVING UP on offset $OFF after 48h of rejections"
      exit 1
    fi
    echo "$(date '+%F %T')  rejected (submit limit); retry $tries in 300s"
    sleep 300
  done
done
echo "$(date '+%F %T')  all slices submitted"
