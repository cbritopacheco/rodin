#!/usr/bin/env bash
# =============================================================================
#  launch_pc9_200.sh -- the whole campaign on pc9 on a budget of CORES = 200:
#  validation, then 3D, then axisymmetric, in that order.
#
#     cd ~/rodin && bash examples/viscoelastic_fluids/launch_pc9_200.sh [mode]
#
#  Modes:
#     start   (default) build, check the meshes and start the queue in the
#             background (nohup, own process group)
#     plan    print the queue (cores, phase, name) without running anything
#     status  queue progress, running cases with step/ETA/Newton, DONE/FAILED
#     stop    stop the queue and every solver it started (DONE cases are kept;
#             a relaunch skips them)
#
#  The queue: jobs start strictly in this order, each as soon as its cores are
#  free (Sum of the cores of the running jobs <= CORES):
#
#   1 validation (3D driver at its default formulation, psi_rb; FORMS adds
#     others, e.g. FORMS="psi_rb psi exp"):
#       - steady sPTT and Oldroyd-B pipe, rf 0.7 1.0 1.4 2.0, all at once
#         (2 fluids x (2+3+9+24) = 76 cores), errors and orders;
#       - pulsatile Oldroyd-B pipe (Womersley), 3D: space (rf 0.7 1.0 1.4) and
#         time (BDF1/BDF2, 16..128 steps per cycle) studies, 8 runs x 8 cores;
#       - the same with the axisymmetric driver, 8 runs x 2 cores.
#   2 3D, mesh inter, BDF2, CY against R (the viscoelastic blood model and its
#     memory-free twin), NP3D = 32 ranks each (6 x 32 = 192):
#       - carotid, Re 300 Wo 4, CCA waveform, 1600 steps/cycle: S75, A200;
#       - Re 600 Wo 5.5, sinusoid A = 0.5, 1000 steps/cycle: S75.
#   3 axisymmetric, NPAXI = 5 ranks each: every geometry (S75 A200 S50 A150 H
#     A125) x Re (600, 300, 100) x fluid group (R | N,CY,D), plus the CCA
#     waveform at Re 300 for S75 and A200 (40 jobs, 200 cores).
#
#  A job never jumps ahead of an earlier one: while the validation runs (156
#  cores) the first 3D case takes 32 of the remaining 44, the other 3D cases
#  start as the validation ends, the axisymmetric ones as 3D cases end.
#  PHASES=1 instead waits until every job of a phase has finished before the
#  next phase starts (e.g. to look at the validation first).
#
#  Environment (defaults):
#     CORES=200  NP3D=32  NPAXI=5  PHASES=0  FORMS=psi_rb
#     SOLVER3D=lulag (MUMPS LU lagged under GMRES; gmres = ASM/ILU, less memory)
#     BUILD=build-gcc12  OUT=$HOME/Results/pc9_200  LAMBDA_R=4.83
#     SKIP_VERIF=0 SKIP_3D=0 SKIP_AXI=0     (1 leaves that phase out)
#  Logs: $OUT/logs/<job>.log, the queue in $OUT/logs/queue.log
#  Results: $OUT/verification/..., $OUT/3D/Re<Re>_Wo<Wo>[_waveform_cca]/<geo>/<fluid>/,
#           $OUT/Re<Re>_Wo<Wo>[_waveform_cca]/<geo>/<fluid>/ (axisymmetric)
# =============================================================================
set -u
WHAT="${1:-start}"
case "$WHAT" in start|plan|status|stop|_run) ;; *) sed -n '2,52p' "$0"; exit 1 ;; esac

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
cd "$REPO"

export PYTHON="${PYTHON:-$(command -v python3)}"
PATH="$(printf '%s' "$PATH" | tr ':' '\n' | grep -viE 'anaconda|miniconda|conda/' | paste -sd: -)"
export PATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export OMPI_MCA_hwloc_base_binding_policy=none

BUILD="${BUILD:-build-gcc12}"
export OUT="${OUT:-$HOME/Results/pc9_200}"
export LAMBDA_R="${LAMBDA_R:-4.83}"
CORES="${CORES:-200}"
NP3D="${NP3D:-32}"
NPAXI="${NPAXI:-5}"
PHASES="${PHASES:-0}"
FORMS="${FORMS:-psi_rb}"
SOLVER3D="${SOLVER3D:-lulag}"
EXE="$REPO/$BUILD/examples/viscoelastic_fluids/ArterialLesionAxi_viscous_logarithmic_implicit/ArterialLesionAxi_viscous_logarithmic_implicit"
EXE3D="$REPO/$BUILD/examples/viscoelastic_fluids/ArterialLesion3D_viscous_logarithmic_implicit/ArterialLesion3D_viscous_logarithmic_implicit"
S3="$HERE/ArterialLesion3D_viscous_logarithmic_implicit/run_3d.sh"
SAXI="$HERE/ArterialLesionAxi_viscous_logarithmic_implicit/run_axi_re.sh"
V="$HERE/ArterialLesion3D_viscous_logarithmic_implicit/verification"
MESH="$REPO/resources/examples/viscoelastic_fluids"
WF="${WAVEFORM:-$MESH/waveform_cca.dat}"
LOG="$OUT/logs"
PIDFILE="$LOG/queue.pid"
mkdir -p "$LOG"

# ---- the queue: "cores|phase|name|command", in order ------------------------------
JOBS=()
add() { JOBS+=("$1|$2|$3|$4"); }
q() { printf '%q ' "$@"; }   # a quoted command line

if [[ "${SKIP_VERIF:-0}" != 1 ]]; then
  for F in $FORMS; do
    add 76 1 "verif_sptt_$F" "$(q env OUT="$OUT/verification/sptt_$F" NP=24 PARALLEL=1 EXE3D="$EXE3D" \
      PETSC_OPTS="-al_formulation $F" bash "$V/run_sptt_pipe.sh")"
    add 64 1 "verif_womersley3d_$F" "$(q env OUT="$OUT/verification/womersley_3d_$F" DRIVER=3d EXE="$EXE3D" \
      NP=8 MAXJOBS=8 PETSC_OPTS="-al_formulation $F" bash "$V/run_womersley.sh")"
  done
  add 16 1 verif_womersley_axi "$(q env OUT="$OUT/verification/womersley_axi" DRIVER=axi EXE="$EXE" \
    NP=2 MAXJOBS=8 bash "$V/run_womersley.sh")"
fi

if [[ "${SKIP_3D:-0}" != 1 ]]; then
  for G in S75 A200; do
    for F in R CY; do
      add "$NP3D" 2 "3D_${G}_Re300_wave_${F}" "$(q env FLUIDS="$F" NP="$NP3D" TIME_ORDER=2 MESH_LEVEL=inter \
        OUTPUT_EVERY=160 WAVEFORM="$WF" HARMONICS=10 STEPS=1600 CYCLES=4 CYCLES_R=5 SOLVER="$SOLVER3D" \
        EXE3D="$EXE3D" OUT3D="$OUT/3D" bash "$S3" 300 "$G")"
    done
  done
  for F in R CY; do
    add "$NP3D" 2 "3D_S75_Re600_${F}" "$(q env FLUIDS="$F" NP="$NP3D" TIME_ORDER=2 MESH_LEVEL=inter \
      OUTPUT_EVERY=100 STEPS=1000 CYCLES=4 CYCLES_R=5 SOLVER="$SOLVER3D" \
      EXE3D="$EXE3D" OUT3D="$OUT/3D" bash "$S3" 600 S75)"
  done
fi

if [[ "${SKIP_AXI:-0}" != 1 ]]; then
  axi() {   # re geo fluids steps [waveform]
    local re="$1" geo="$2" fl="$3" st="$4" wave="${5:-}" tag
    tag="axi_Re${re}_${geo}_${fl//,/}"
    [[ -n "$wave" ]] && tag="${tag}_wave"
    add "$NPAXI" 3 "$tag" "$(q env ${wave:+WAVEFORM="$wave"} STEPS="$st" FLUIDS="${fl//,/ }" NP="$NPAXI" \
      TIME_ORDER=2 NEWTON_TOL=5e-4 HARMONICS=10 CYCLES=4 CYCLES_R=8 RETRY=1 OUT="$OUT" EXE="$EXE" \
      bash "$SAXI" "$re" "$geo")"
  }
  for RE in 600 300 100; do   # longest first
    case "$RE" in 600) ST=1200 ;; 300) ST=800 ;; 100) ST=600 ;; esac
    for G in S75 A200 S50 A150 H A125; do
      axi "$RE" "$G" R "$ST"
      axi "$RE" "$G" N,CY,D "$ST"
    done
  done
  for G in S75 A200; do
    axi 300 "$G" R 1600 "$WF"
    axi 300 "$G" N,CY,D 1600 "$WF"
  done
fi

if [[ "$WHAT" == plan ]]; then
  total=0
  printf '%5s %5s  %s\n' cores phase job
  for j in "${JOBS[@]}"; do
    IFS='|' read -r c p n _ <<< "$j"
    printf '%5s %5s  %s\n' "$c" "$p" "$n"
    total=$(( total + c ))
  done
  echo "${#JOBS[@]} jobs, $total core-slots, budget $CORES cores, PHASES=$PHASES"
  exit 0
fi

csvcol() {   # file column-name: last value of that column
  awk -F, -v c="$2" 'NR==1{for(i=1;i<=NF;i++) if($i==c) k=i; next} END{if(k) print $k; else print "-"}' "$1" 2> /dev/null
}

if [[ "$WHAT" == status ]]; then
  if [[ -f "$PIDFILE" ]] && kill -0 "$(cat "$PIDFILE")" 2> /dev/null; then
    echo "== queue running (pid $(cat "$PIDFILE"))"
  else
    echo "== queue not running"
  fi
  echo "== $(date '+%a %d %H:%M')  solver ranks: 3D $(ps -eo args | grep -c '^[^ ]*ArterialLesion3D_viscous_logarithmic_implicit -al_')" \
       " axi $(ps -eo args | grep -c '^[^ ]*ArterialLesionAxi_viscous_logarithmic_implicit -al_')"
  [[ -f "$LOG/queue.log" ]] && { echo "== queue:"; grep -E 'start|end|waiting|all jobs' "$LOG/queue.log" | tail -12; }
  echo "== running cases (run.log written in the last 15 min):"
  find "$OUT" -name run.log -mmin -15 2> /dev/null | sort | while read -r f; do
    d=$(dirname "$f")
    total=$(grep -a -m1 -- '---- Step' "$f" | sed 's/.*Step [0-9]*\/\([0-9]*\).*/\1/')
    done_=$(( $(wc -l < "$d/case.csv" 2> /dev/null || echo 1) - 1 ))
    eta=$(grep -a 'ETA' "$f" | tail -1 | sed 's/.*ETA/ETA/')
    printf '   %-50s step %s/%s  cycle %s  newtonIts %s  retries %s  %s\n' "${d#$OUT/}" "$done_" "$total" \
      "$(csvcol "$d/case.csv" cycle)" "$(csvcol "$d/case.csv" newtonIts)" "$(csvcol "$d/case.csv" retries)" "$eta"
  done
  echo "== finished: $(find "$OUT" -name DONE 2> /dev/null | wc -l)   failed:"
  find "$OUT" -name FAILED 2> /dev/null | sed "s|$OUT/|   |; s|/FAILED||"
  exit 0
fi

if [[ "$WHAT" == stop ]]; then
  if [[ -f "$PIDFILE" ]]; then
    pg="$(cat "$PIDFILE")"
    echo "== stopping the queue (process group $pg) ..."
    kill -TERM -- "-$pg" 2> /dev/null; sleep 10
    kill -KILL -- "-$pg" 2> /dev/null
    rm -f "$PIDFILE"
  else
    echo "== no queue pid file ($PIDFILE)"
  fi
  # Unfinished cases lose their FAILED mark, so that a relaunch re-runs them.
  find "$OUT" -name FAILED 2> /dev/null | while read -r f; do
    d=$(dirname "$f"); [[ -f "$d/DONE" ]] || rm -f "$f"
  done
  exit 0
fi

# ---- the scheduler (run in the background by start) -------------------------------
if [[ "$WHAT" == _run ]]; then
  echo "$$" > "$PIDFILE"   # the session (process group) leader: stop kills the group
  pids=(); cores=(); phases=(); names=()
  i=0; n=${#JOBS[@]}
  echo "[queue] $(date '+%a %d %H:%M') $n jobs, $CORES cores, PHASES=$PHASES"
  waiting=""
  while (( i < n || ${#pids[@]} > 0 )); do
    # Reap the jobs that have ended.
    np=(); nc=(); nph=(); nn=()
    for k in ${pids[@]+"${!pids[@]}"}; do
      if kill -0 "${pids[k]}" 2> /dev/null; then
        np+=("${pids[k]}"); nc+=("${cores[k]}"); nph+=("${phases[k]}"); nn+=("${names[k]}")
      else
        wait "${pids[k]}"; rc=$?
        echo "[queue] $(date '+%a %d %H:%M') end   ${names[k]} (exit $rc)"
      fi
    done
    pids=(${np[@]+"${np[@]}"}); cores=(${nc[@]+"${nc[@]}"}); phases=(${nph[@]+"${nph[@]}"}); names=(${nn[@]+"${nn[@]}"})
    used=0; minphase=99
    for k in ${pids[@]+"${!pids[@]}"}; do
      used=$(( used + cores[k] ))
      (( phases[k] < minphase )) && minphase=${phases[k]}
    done
    # Start the next jobs, in order, while they fit.
    while (( i < n )); do
      IFS='|' read -r c p name cmd <<< "${JOBS[i]}"
      why=""
      (( used + c > CORES )) && why="$c cores (in use $used)"
      [[ "$PHASES" == 1 ]] && (( minphase < p )) && why="the end of phase $minphase"
      if [[ -n "$why" ]]; then
        if [[ "$waiting" != "$name" ]]; then
          echo "[queue] $(date '+%a %d %H:%M') $name waiting for $why"
          waiting="$name"
        fi
        break
      fi
      echo "[queue] $(date '+%a %d %H:%M') start $name ($c cores, phase $p; in use $(( used + c )))"
      bash -c "$cmd" > "$LOG/$name.log" 2>&1 < /dev/null &
      pids+=("$!"); cores+=("$c"); phases+=("$p"); names+=("$name")
      used=$(( used + c )); (( p < minphase )) && minphase=$p
      i=$(( i + 1 ))
      sleep 2
    done
    (( i < n || ${#pids[@]} > 0 )) && sleep 30
  done
  echo "[queue] $(date '+%a %d %H:%M') all jobs finished"
  rm -f "$PIDFILE"
  exit 0
fi

# ---- start: checks, build, launch -----------------------------------------------------
if [[ -f "$PIDFILE" ]] && kill -0 "$(cat "$PIDFILE")" 2> /dev/null; then
  echo "A queue is already running (pid $(cat "$PIDFILE")): bash $0 status | stop"; exit 1
fi
echo "== repo $REPO ($(git --no-optional-locks log --oneline -1 2> /dev/null))"
echo "== mpirun: $(command -v mpirun)"
[[ -d "$REPO/$BUILD" ]] || { echo "Build directory $REPO/$BUILD not found (set BUILD=...)"; exit 1; }
if [[ "${SKIP_BUILD:-0}" != 1 ]]; then
  echo "== building (the 3D driver needs ~6 GB per compiler process) ..."
  if ! cmake --build "$REPO/$BUILD" -j 8 --target ArterialLesionAxi_viscous_logarithmic_implicit \
        ArterialLesion3D_viscous_logarithmic_implicit > "$LOG/build.log" 2>&1; then
    echo "BUILD FAILED, see $LOG/build.log:"; tail -20 "$LOG/build.log"; exit 1
  fi
fi
for x in "$EXE" "$EXE3D"; do
  [[ -x "$x" ]] || { echo "missing executable $x"; exit 1; }
done
grep -aq -- "psi_rb" "$EXE3D" || { echo "$EXE3D has no psi_rb formulation: wrong sources?"; exit 1; }
"$PYTHON" -c "import numpy, h5py" 2> /dev/null || { echo "$PYTHON has no numpy/h5py (set PYTHON=...)"; exit 1; }

need=()
[[ "${SKIP_VERIF:-0}" != 1 ]] && for rf in 0.7 1.0 1.4 2.0; do need+=("pipe_verification_rf$rf"); done
[[ "${SKIP_3D:-0}" != 1 ]] && need+=(S75_pipe_inter A200_pipe_inter)
[[ "${SKIP_AXI:-0}" != 1 ]] && need+=(S75_axi_medium S50_axi_medium H_axi_medium A125_axi_medium A150_axi_medium A200_axi_medium)
for m in ${need[@]+"${need[@]}"}; do
  f="$MESH/$m.mesh"
  if [[ ! -f "$f" ]] || ! head -1 "$f" | grep -q MeshVersionFormatted; then
    echo "mesh $f missing or not MEDIT (Git LFS pointer? run: git lfs pull)"; exit 1
  fi
done
[[ -f "$WF" ]] || { echo "waveform $WF not found"; exit 1; }
if [[ "${SKIP_3D:-0}" != 1 ]]; then
  echo "== checking the 3D meshes ..."
  "$PYTHON" "$HERE/ArterialLesion3D_viscous_logarithmic_implicit/check_3d_mesh.py" \
    S75 "$MESH/S75_pipe_inter.mesh" A200 "$MESH/A200_pipe_inter.mesh" | tee "$LOG/check_3d_mesh.log"
  grep -q PROBLEMS "$LOG/check_3d_mesh.log" && { echo "3D mesh check failed"; exit 1; }
fi

busy=$(ps -eo comm | grep -c '^ArterialLesion')
if (( busy > 0 )) && [[ "${FORCE:-0}" != 1 ]]; then
  echo "$busy solver processes are already running (other runs?). Stop them first, or FORCE=1"
  echo "(the queue counts only its own jobs against CORES=$CORES)."
  exit 1
fi

bash "$0" plan | tail -1
rm -f "$PIDFILE"
setsid nohup env SKIP_BUILD=1 bash "$0" _run > "$LOG/queue.log" 2>&1 < /dev/null &
sleep 5
[[ -f "$PIDFILE" ]] || { echo "the queue did not start, see $LOG/queue.log"; exit 1; }
echo "== queue started (pid $(cat "$PIDFILE")) -> $LOG/queue.log"
echo "Status:  bash $0 status"
echo "Stop:    bash $0 stop"
echo "Plan:    bash $0 plan"
