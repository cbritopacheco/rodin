#!/usr/bin/env bash
# =============================================================================
#  run_axi_queue.sh -- runs a list of run_axi_re.sh jobs, MAXJOBS at a time.
#
#     MAXJOBS=6 NP=5 bash run_axi_queue.sh jobs.txt
#
#  jobs.txt, one job per line (blank lines and # comments ignored):
#     <Re> <geometry> <fluids, comma separated> [steps per cycle] [waveform file]
#  e.g.
#     600 S75 R 1600
#     600 S75 N,CY,D 1600
#     300 A200 R 1600 /path/waveform_cca.dat
#  Every other setting (TIME_ORDER, NEWTON_TOL, OUT, EXE, CYCLES, ...) is taken
#  from the environment and passed on to run_axi_re.sh. Finished cases (DONE)
#  are skipped by run_axi_re.sh, so the queue can be relaunched.
#  Logs: $LOGDIR/axi_Re<Re>_<geo>_<fluids>[_wave].log
# =============================================================================
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
S="$HERE/run_axi_re.sh"
JOBS="${1:?usage: run_axi_queue.sh <jobs file>}"
MAXJOBS="${MAXJOBS:-6}"
LOGDIR="${LOGDIR:-${OUT:-$HOME/Results/axi}/logs}"
mkdir -p "$LOGDIR"

running() { jobs -rp | wc -l; }

while read -r re geo fluids steps wave; do
  [[ -z "${re:-}" || "$re" == \#* ]] && continue
  while (( $(running) >= MAXJOBS )); do sleep 30; done
  tag="axi_Re${re}_${geo}_${fluids//,/}"
  [[ -n "${wave:-}" ]] && tag="${tag}_wave"
  echo "[queue] $(date '+%a %d %H:%M') start $tag"
  (
    if [[ -n "${wave:-}" ]]; then export WAVEFORM="$wave"; else unset WAVEFORM; fi
    [[ -n "${steps:-}" ]] && export STEPS="$steps"
    FLUIDS="${fluids//,/ }" bash "$S" "$re" "$geo" > "$LOGDIR/$tag.log" 2>&1 < /dev/null
    echo "[queue] $(date '+%a %d %H:%M') end   $tag (exit $?)"
  ) &
  sleep 2
done < "$JOBS"
wait
echo "[queue] $(date '+%a %d %H:%M') all jobs finished"
