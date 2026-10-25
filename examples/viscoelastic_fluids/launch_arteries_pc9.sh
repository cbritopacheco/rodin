#!/usr/bin/env bash
# =============================================================================
#  launch_arteries_pc9.sh -- the definitive campaign on pc9 (CORES = 200):
#  idealised lesions at the Reynolds and Womersley numbers of five arteries,
#  three fluids, 3D; then the axisymmetric Re x geometry sweep.
#
#     cd ~/rodin && bash examples/viscoelastic_fluids/launch_arteries_pc9.sh [start|plan|status|stop]
#
#  Fluids (the same steady shear viscosity, Carreau-Yasuda, for all three):
#     N   Newtonian at eta_inf (sPTT at Wi 0.01)
#     CY  Carreau-Yasuda, no memory
#     O   Owens (2006) aggregating blood (run_3d.sh; LAMBDA_H, ETA_PLASMA, A0_O,
#         ACRIT_O, EPS_O, N0_O, NIN_O): lambda = n lambda_H in rouleaux of size n(x, t),
#         lambda_H at high shear; add it with FLUIDS3D="O R CY N"
#     R   gsPTT: sPTT polymer (lambda 0.1 s, eta_p 0.0070 Pa s, eps 0.2) plus the
#         solvent eta_CY - eta_p/f. Its shear stress must INCREASE with gd
#         (eta_p <= 7.2e-3 at lambda 0.1 s; the driver checks it): with eta_p
#         0.0104 eta_s fell to 1.6e-5 Pa s and its stress decreased for 1.75 <
#         gd < 16 1/s, an ill-posed momentum equation (runs 9, 12, 13). The
#         pair (lambda, eta_p) is the one with the largest 2 Hz loss viscosity
#         eta'' under that constraint (eta'' 0.0034 Pa s, G = eta_p/lambda 0.07
#         Pa, De = lambda/T 0.12 at 72 bpm; Thurston's blood: eta'' ~ 0.01,
#         G ~ 0.1 Pa, lambda 0.01-0.1 s). lambda 4.83 s (De 5.8, lambda gd ~ 5e4
#         at the throat) was abandoned: not blood, and beyond the method.
#
#  Arteries (D, mean velocity; period T = 0.83 s, 72 bpm; Re = rho U D/eta_inf,
#  Wo = (D/2) sqrt(2 pi rho/(T eta_inf)), rho 1050, eta_inf 3.45e-3), the
#  steps per cycle N = 40 U T/D in [600, 1600] (about 12 steps per shedding
#  period of the throat jet, St 0.3, at the systolic peak) and the mesh:
#     small     D 2.5 mm  U 0.08  Re  61  Wo  1.9  Wi 3.2  N 1100  layers    (radial/digital)
#     coronary  D 3.5 mm  U 0.20  Re 213  Wo  2.7  Wi 5.7  N 1600  layers    (LAD)
#     carotid   D 6.0 mm  U 0.35  Re 639  Wo  4.6  Wi 5.8  N 1600  layers    CCA waveform
#     femoral   D 7.0 mm  U 0.15  Re 320  Wo  5.3  Wi 2.1  N  700  layers
#     aorta     D 15 mm   U 0.18  Re 822  Wo 11.4  Wi 1.2  N  600  layerswo  (infrarenal)
#  (Wi = lambda U/D; at the throat wall of S75 at the systolic peak lambda gd
#  is about 100 Wi. De = lambda/T = 0.12 for all.)
#  Inflow: the CCA waveform for the carotid, the sinusoid A = 0.5 elsewhere;
#  4 cycles (ramp + 3), XDMF 20 times per cycle. Geometries S75 and A200
#  (GEOS adds S50 etc.). Meshes: <geo>_pipe_layers (rf 1, the throat wall and
#  shear layer at h 0.035 D / 0.025 D, 244k-267k tets) and, for Wo 10, the
#  Stokes-layer meshes <geo>_pipe_layerswo (parent wall at 0.04 D), NPWO ranks.
#  Order: carotid, coronary, aorta, femoral, small; R, CY, N; S75, A200.
#
#  Phase 3 (SKIP_AXI=0; off by default): the axisymmetric sweep of
#  launch_pc9_200.sh with the new R.
#
#  Modes, queue, Scotch and partition checks: as launch_pc9_200.sh.
#  Environment (defaults):
#     CORES=200  NP3D=32  NPWO=64  NPAXI=5  PHASES=0  GEOS="S75 A200"
#     ARTERIES="carotid coronary aorta femoral small"  FLUIDS3D="R CY N"
#     SOLVER3D=lulag  BUILD=build-gcc12  OUT=$HOME/Results/arteries
#     LAMBDA_R=0.1  ETAP_R=0.0070  SKIP_3D=0  SKIP_AXI=1
#  Results: $OUT/3D/<artery>/Re<Re>_Wo<Wo>[_waveform_cca]/<geo>/<fluid>/
# =============================================================================
set -u
WHAT="${1:-start}"
case "$WHAT" in start|plan|status|stop|_run) ;; *) sed -n '2,45p' "$0"; exit 1 ;; esac

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
cd "$REPO"

export PYTHON="${PYTHON:-$(command -v python3)}"
PATH="$(printf '%s' "$PATH" | tr ':' '\n' | grep -viE 'anaconda|miniconda|conda/' | paste -sd: -)"
export PATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export OMPI_MCA_hwloc_base_binding_policy=none

BUILD="${BUILD:-build-gcc12}"
export OUT="${OUT:-$HOME/Results/arteries}"
export LAMBDA_R="${LAMBDA_R:-0.1}"
export ETAP_R="${ETAP_R:-0.0070}"
export LAMBDA_H="${LAMBDA_H:-0.005}" ETA_PLASMA="${ETA_PLASMA:-1.2e-3}"
export A0_O="${A0_O:-1}" ACRIT_O="${ACRIT_O:--1}" EPS_O="${EPS_O:-0}" N0_O="${N0_O:--1}" NIN_O="${NIN_O:--1}"
NPWO="${NPWO:-64}"
GEOS="${GEOS:-S75 A200}"
ARTERIES="${ARTERIES:-carotid coronary aorta femoral small}"
FLUIDS3D="${FLUIDS3D:-R CY N}"
SKIP_AXI="${SKIP_AXI:-1}"
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
declare -A INFO
add() { JOBS+=("$1|$2|$3|$4"); INFO["$3"]="${5:-}"; }
q() { printf '%q ' "$@"; }   # a quoted command line

RHO=1050; ETAINF=3.45e-3; PERIOD=0.83
artery() {   # name -> D U mesh waveform(0/1)
  case "$1" in
    small)    echo 0.0025 0.08 layers   0 ;;
    coronary) echo 0.0035 0.20 layers   0 ;;
    carotid)  echo 0.0060 0.35 layers   1 ;;
    femoral)  echo 0.0070 0.15 layers   0 ;;
    aorta)    echo 0.0150 0.18 layerswo 0 ;;
    *) echo "unknown artery $1" >&2; exit 1 ;;
  esac
}
derive() {   # D U -> Re Wo steps  (Re = rho U D/eta_inf; Wo = D/2 sqrt(2 pi rho/(T eta_inf)))
  # LC_ALL=C: mawk prints Wo with the locale's decimal separator ("4,6" under
  # es_CL), which PETSc rejects (-al_wo 4,6: "has no numeric value").
  LC_ALL=C awk -v D="$1" -v U="$2" -v rho="$RHO" -v eta="$ETAINF" -v T="$PERIOD" 'BEGIN {
    re = rho*U*D/eta; wo = D/2*sqrt(2*3.14159265*rho/(T*eta));
    n = int(40*U*T/D/100 + 0.5)*100; if (n < 600) n = 600; if (n > 1600) n = 1600;
    printf "%d %.1f %d\n", re + 0.5, wo, n }'
}

if [[ "${SKIP_3D:-0}" != 1 ]]; then
  for A in $ARTERIES; do
    read -r D U ML WV <<< "$(artery "$A")"
    read -r RE WO ST <<< "$(derive "$D" "$U")"
    NP="$NP3D"; [[ "$ML" == layerswo ]] && NP="$NPWO"
    WAVE=""; [[ "$WV" == 1 ]] && WAVE="$WF"
    for F in $FLUIDS3D; do
      for G in $GEOS; do
        add "$NP" 2 "3D_${A}_${G}_${F}" "$(q env FLUIDS="$F" NP="$NP" TIME_ORDER=2 MESH_LEVEL="$ML" WO="$WO" \
          OUTPUT_EVERY=$(( ST / 20 )) STEPS="$ST" CYCLES=4 CYCLES_R=4 SOLVER="$SOLVER3D" HARMONICS=10 \
          WAVEFORM="$WAVE" PETSC_EXTRA="-al_diameter $D" \
          EXE3D="$EXE3D" OUT3D="$OUT/3D/$A" bash "$S3" "$RE" "$G")" \
          "Re $RE  Wo $WO  D $D m  U $U m/s  $ST steps/cycle  $ML"
      done
    done
  done
fi

if [[ "$SKIP_AXI" != 1 ]]; then
  axi() {   # re geo fluids steps [waveform]
    local re="$1" geo="$2" fl="$3" st="$4" wave="${5:-}" tag
    tag="axi_Re${re}_${geo}_${fl//,/}"
    [[ -n "$wave" ]] && tag="${tag}_wave"
    add "$NPAXI" 3 "$tag" "$(q env ${wave:+WAVEFORM="$wave"} STEPS="$st" FLUIDS="${fl//,/ }" NP="$NPAXI" \
      TIME_ORDER=2 NEWTON_TOL=5e-4 HARMONICS=10 CYCLES=4 CYCLES_R=4 RETRY=1 OUT="$OUT" EXE="$EXE" \
      bash "$SAXI" "$re" "$geo")"
  }
  for RE in 600 300 100; do
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
  printf '%5s %5s  %-24s %s\n' cores phase job ""
  for j in "${JOBS[@]}"; do
    IFS='|' read -r c p n cmd <<< "$j"
    printf '%5s %5s  %-24s %s\n' "$c" "$p" "$n" "${INFO[$n]}"
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
  echo "== partitions of the running cases:"
  find "$OUT" -name run.log -mmin -15 2> /dev/null | sort | while read -r f; do
    p=$(grep -a -m1 '\[partition\]' "$f" | sed 's/\x1b\[[0-9;]*m//g; s/.*\[partition\] //')
    printf '   %-50s %s\n' "$(dirname "${f#$OUT/}")" "${p:--}"
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
if ! grep -qiE '^RODIN_USE_SCOTCH:BOOL=(ON|YES|TRUE|Y|1)$' "$REPO/$BUILD/CMakeCache.txt" 2> /dev/null; then
  echo "$BUILD is configured without Scotch (RODIN_USE_SCOTCH is not ON): the drivers would use the"
  echo "internal partitioner, much slower. Reconfigure and rebuild:"
  echo "   cmake -S $REPO -B $REPO/$BUILD -DRODIN_USE_SCOTCH=ON [-DSCOTCH_DIR=<scotch prefix>]"
  echo "(or ALLOW_NO_SCOTCH=1 to run anyway)"
  [[ "${ALLOW_NO_SCOTCH:-0}" == 1 ]] || exit 1
fi
grep -aq -- "\[partition\]" "$EXE3D" || { echo "$EXE3D does not report its partition: wrong sources?"; exit 1; }
"$PYTHON" -c "import numpy, h5py" 2> /dev/null || { echo "$PYTHON has no numpy/h5py (set PYTHON=...)"; exit 1; }

need=()
if [[ "${SKIP_3D:-0}" != 1 ]]; then
  for G in $GEOS; do
    for A in $ARTERIES; do read -r _ _ ML _ <<< "$(artery "$A")"; need+=("${G}_pipe_$ML"); done
  done
  need=($(printf '%s\n' "${need[@]}" | sort -u))
fi
[[ "$SKIP_AXI" != 1 ]] && need+=(S75_axi_medium S50_axi_medium H_axi_medium A125_axi_medium A150_axi_medium A200_axi_medium)
for m in ${need[@]+"${need[@]}"}; do
  f="$MESH/$m.mesh"
  if [[ ! -f "$f" ]] || ! head -1 "$f" | grep -q MeshVersionFormatted; then
    echo "mesh $f missing or not MEDIT (Git LFS pointer? run: git lfs pull)"; exit 1
  fi
done
[[ -f "$WF" ]] || { echo "waveform $WF not found"; exit 1; }
if [[ "${SKIP_3D:-0}" != 1 ]]; then
  echo "== checking the 3D meshes ..."
  args=(); for m in "${need[@]}"; do [[ "$m" == *_pipe_* ]] && args+=("${m%%_*}" "$MESH/$m.mesh"); done
  "$PYTHON" "$HERE/ArterialLesion3D_viscous_logarithmic_implicit/check_3d_mesh.py" "${args[@]}" | tee "$LOG/check_3d_mesh.log"
  grep -q PROBLEMS "$LOG/check_3d_mesh.log" && { echo "3D mesh check failed"; exit 1; }
fi

busy=$(ps -eo comm | grep -c '^ArterialLesion')
if (( busy > 0 )) && [[ "${FORCE:-0}" != 1 ]]; then
  echo "$busy solver processes are already running (other runs?). Stop them first, or FORCE=1"
  echo "(the queue counts only its own jobs against CORES=$CORES)."
  exit 1
fi

# Partition check: the 3D driver on the 3D mesh with NP3D ranks, stopped as
# soon as it has printed its [partition] line.
if [[ "${SKIP_3D:-0}" != 1 && "${SKIP_PARTITION_CHECK:-0}" != 1 ]]; then
  echo "== partition check: 3D driver, S75_pipe_layers, $NP3D ranks ..."
  PC="$LOG/partition_check"; mkdir -p "$PC"
  ( cd "$PC" && exec ${MPIRUN:-mpirun} -n "$NP3D" "$EXE3D" -al_mesh "$MESH/S75_pipe_layers.mesh" \
      -al_xdmf case -al_csv case.csv -al_cycles 1 -al_steps_per_cycle 10 ) > "$PC/run.log" 2>&1 < /dev/null &
  pcpid=$!
  for _ in $(seq 1 120); do
    grep -aq '\[partition\]' "$PC/run.log" && break
    kill -0 "$pcpid" 2> /dev/null || break
    sleep 5
  done
  # SIGTERM to mpirun ends its ranks.
  kill -TERM "$pcpid" 2> /dev/null
  for _ in $(seq 1 20); do kill -0 "$pcpid" 2> /dev/null || break; sleep 1; done
  kill -KILL "$pcpid" 2> /dev/null; wait "$pcpid" 2> /dev/null
  line=$(grep -a '\[partition\]' "$PC/run.log" | head -1 | sed 's/\x1b\[[0-9;]*m//g')
  if [[ -z "$line" ]]; then
    echo "no [partition] line from the driver, see $PC/run.log:"; tail -5 "$PC/run.log"; exit 1
  fi
  echo "   $line"
  ratio=$(printf '%s' "$line" | sed -n 's/.*max\/mean = \([0-9.eE+-]*\).*/\1/p')
  if ! awk -v r="$ratio" -v m="${MAX_IMBALANCE:-1.15}" 'BEGIN{exit !(r != "" && r <= m)}'; then
    echo "partition imbalance max/mean = ${ratio:-?} > ${MAX_IMBALANCE:-1.15}: fix the partitioner first"
    echo "(MAX_IMBALANCE=... or SKIP_PARTITION_CHECK=1 to go on anyway)"
    exit 1
  fi
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
