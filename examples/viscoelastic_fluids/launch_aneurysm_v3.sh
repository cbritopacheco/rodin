#!/usr/bin/env bash
#  Campaign v3: blood rheology in an idealised fusiform aneurysm (A200).
#
#  Design (dimensionless): with the blood's relaxation time fixed, De = lambda/T
#  is fixed and Wi = lambda U/D = De UT/D, so a case is a point of the (Wo,
#  UT/D) plane, UT/D the stroke length of the pulse (pi/2 Re/Wo^2):
#
#     point  Wo    UT/D  Wi    Re   D (mm)  U (m/s)  mesh      steps/T
#     A      2.7   47    6.5   213  3.5     0.20     layers    1600
#     B      4.6   48    6.6   639  6.0     0.35     layers    1600
#     C      5.3   10    1.4   179  7.0     0.084    layers    1200
#     D      11.4  10    1.4   822  15      0.18     layerswo  1200
#     E      1.9   27    3.6    61  2.5     0.08     layers    1200
#
#  (A, B): long stroke, the sac renewed every cycle; (C, D): short stroke, a
#  trapped population; each pair at two Wo. E: the viscous limit. Re, Wo with
#  rho = 1050 kg/m^3, eta_inf = 3.45 mPa s, T = 0.83 s. Inflow: the
#  physiological CCA pulse (Holdsworth et al. 1999, 10 harmonics, no reversal)
#  for every point. Fluids N (Newtonian, eta_inf, two-field), CY (Carreau-
#  Yasuda, two-field) and R (gsPTT, lambda 0.114 s, eta_p 0.0068 Pa s, eps 0.5),
#  5 cycles (the first ramped), BDF2.
#
#  Scheduling: JOBS cases at a time, NP ranks each (default 3 x 60); a slot is
#  given to the next case as soon as one ends (R first, the longest). When a
#  case ends, Pathlines3D runs on it in the background (one core) and the
#  slot is freed; when every case is done, postprocess_3d.py and
#  pathlines_figures.py run on the whole campaign. Finished cases (DONE) are
#  skipped, so the script can be re-run after an interruption.
#
#  Usage (from the repository root, on pc9):
#     nohup bash examples/viscoelastic_fluids/launch_aneurysm_v3.sh > ~/aneurysm_v3.log 2>&1 &
#     bash examples/viscoelastic_fluids/launch_aneurysm_v3.sh status
#     bash examples/viscoelastic_fluids/launch_aneurysm_v3.sh stop
#
#  Environment (defaults):
#     OUT=$HOME/Results/aneurysm_v3  BUILD=build-gcc12  NP=60  JOBS=3
#     POINTS="A B C D E"  FLUIDS="R CY N"  GEO=A200  CYCLES=5
#     LAMBDA_R=0.114  ETAP_R=0.0068  EPS_R=0.5  SOLVER=lulag
#     PL_SEEDS=1000  PL_RELEASES=16  PL_HORIZON=8  SKIP_BUILD=0
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
OUT="${OUT:-$HOME/Results/aneurysm_v3}"
BUILD="${BUILD:-build-gcc12}"
NP="${NP:-60}"
JOBS="${JOBS:-3}"
POINTS="${POINTS:-A B C D E}"
FLUIDS="${FLUIDS:-R CY N}"
GEO="${GEO:-A200}"
CYCLES="${CYCLES:-5}"
export LAMBDA_R="${LAMBDA_R:-0.114}" ETAP_R="${ETAP_R:-0.0068}" EPS_R="${EPS_R:-0.5}"
SOLVER="${SOLVER:-lulag}"
PL_SEEDS="${PL_SEEDS:-1000}"
PL_RELEASES="${PL_RELEASES:-16}"
PL_HORIZON="${PL_HORIZON:-8}"
WAVE="$REPO/resources/examples/viscoelastic_fluids/waveform_cca.dat"
EXE3D="$REPO/$BUILD/examples/viscoelastic_fluids/ArterialLesion3D_viscous_logarithmic_implicit/ArterialLesion3D_viscous_logarithmic_implicit"
EXEPL="$REPO/$BUILD/examples/viscoelastic_fluids/Pathlines3D/Pathlines3D"
S3="$HERE/ArterialLesion3D_viscous_logarithmic_implicit/run_3d.sh"
POST="$HERE/ArterialLesion3D_viscous_logarithmic_implicit/postprocess_3d.py"
FIGS="$HERE/Pathlines3D/pathlines_figures.py"
LOG="$OUT/logs"
PIDFILE="$LOG/launcher.pid"
mkdir -p "$LOG"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export OMPI_MCA_hwloc_base_binding_policy=none

RHO=1050; ETAINF=3.45e-3; PERIOD=0.83
point() {   # name -> D U mesh steps
  case "$1" in
    A) echo 0.0035 0.20  layers   1600 ;;
    B) echo 0.0060 0.35  layers   1600 ;;
    C) echo 0.0070 0.084 layers   1200 ;;
    D) echo 0.0150 0.18  layerswo 1200 ;;
    E) echo 0.0025 0.08  layers   1200 ;;
    *) echo "unknown point $1" >&2; return 1 ;;
  esac
}
groups() {   # D U -> Re Wo UT/D (Re = rho U D/eta_inf; Wo = D/2 sqrt(2 pi rho/(T eta_inf)))
  LC_ALL=C awk -v D="$1" -v U="$2" -v rho="$RHO" -v eta="$ETAINF" -v T="$PERIOD" 'BEGIN {
    printf "%d %.1f %.1f\n", rho*U*D/eta + 0.5, D/2*sqrt(2*3.14159265*rho/(T*eta)), U*T/D }'
}
casedir() {   # point fluid -> the folder run_3d.sh writes
  read -r D U ML ST <<< "$(point "$1")"
  read -r RE WO L <<< "$(groups "$D" "$U")"
  echo "$OUT/3D/$1/Re${RE}_Wo${WO}_waveform_cca/$GEO/$2"
}

status() {
  echo "== $OUT"
  for P in $POINTS; do
    read -r D U ML ST <<< "$(point "$P")"; read -r RE WO L <<< "$(groups "$D" "$U")"
    printf "%s  Re %-4s Wo %-5s UT/D %-5s |" "$P" "$RE" "$WO" "$L"
    for F in $FLUIDS; do
      d="$(casedir "$P" "$F")"; s="-"
      [[ -f "$d/run.log" ]] && s="run $(grep -o 'Step [0-9]*/[0-9]*' "$d/run.log" | tail -1 | cut -d' ' -f2)"
      [[ -f "$d/FAILED" ]] && s="FAILED"
      [[ -f "$d/DONE" ]] && s="DONE"
      [[ -f "$d/pathlines_lesion.csv" ]] && s="$s+PL"
      printf " %s: %-16s" "$F" "$s"
    done
    echo
  done
}
case "${1:-start}" in
  status) status; exit 0 ;;
  stop)
    [[ -f "$PIDFILE" ]] && pkill -TERM -P "$(cat "$PIDFILE")" 2>/dev/null; [[ -f "$PIDFILE" ]] && kill "$(cat "$PIDFILE")" 2>/dev/null
    pkill -u "$USER" -f "$EXE3D" 2>/dev/null; rm -f "$PIDFILE"; echo stopped; exit 0 ;;
  start) ;;
  *) echo "usage: $0 [start|status|stop]"; exit 1 ;;
esac
if [[ -f "$PIDFILE" ]] && kill -0 "$(cat "$PIDFILE")" 2>/dev/null; then
  echo "a launcher is already running (pid $(cat "$PIDFILE"))"; exit 1
fi
echo $$ > "$PIDFILE"
trap 'rm -f "$PIDFILE"' EXIT

if [[ "${SKIP_BUILD:-0}" != 1 ]]; then
  cmake --build "$REPO/$BUILD" -j 16 --target ArterialLesion3D_viscous_logarithmic_implicit Pathlines3D \
    > "$LOG/build.log" 2>&1 || { echo "build failed: $LOG/build.log"; exit 1; }
fi
[[ -x "$EXE3D" && -x "$EXEPL" ]] || { echo "executables missing under $REPO/$BUILD"; exit 1; }

pathlines() {   # folder: in the background, once
  local d="$1"
  [[ -f "$d/DONE" && ! -f "$d/pathlines_lesion.csv" ]] || return 0
  nohup "$EXEPL" -pl_case "$d" -pl_seeds "$PL_SEEDS" -pl_releases "$PL_RELEASES" -pl_horizon "$PL_HORIZON" \
    > "$d/pathlines.log" 2>&1 &
}
job() {   # point fluid
  local P="$1" F="$2" D U ML ST RE WO L
  read -r D U ML ST <<< "$(point "$P")"; read -r RE WO L <<< "$(groups "$D" "$U")"
  echo "[$(date '+%a %d %H:%M')] start $P $F (Re $RE Wo $WO, $ST steps/T, $ML, $NP ranks)"
  env FLUIDS="$F" NP="$NP" TIME_ORDER=2 MESH_LEVEL="$ML" WO="$WO" STEPS="$ST" CYCLES="$CYCLES" CYCLES_R="$CYCLES" \
    OUTPUT_EVERY=$(( ST / 20 )) SOLVER="$SOLVER" WAVEFORM="$WAVE" HARMONICS=10 RETRY=1 RESIDENCE=1 \
    PETSC_EXTRA="-al_diameter $D" EXE3D="$EXE3D" OUT3D="$OUT/3D/$P" \
    bash "$S3" "$RE" "$GEO" > "$LOG/${P}_${F}.log" 2>&1
  local st=$?
  echo "[$(date '+%a %d %H:%M')] end   $P $F (exit $st)"
  pathlines "$(casedir "$P" "$F")"
  return $st
}

# The queue: R first (the longest), then CY and N, point by point.
QUEUE=()
for F in R CY N; do
  [[ " $FLUIDS " == *" $F "* ]] || continue
  for P in $POINTS; do QUEUE+=("$P $F"); done
done
echo "[$(date '+%a %d %H:%M')] ${#QUEUE[@]} cases, $JOBS at a time x $NP ranks -> $OUT"
running=0
for item in "${QUEUE[@]}"; do
  read -r P F <<< "$item"
  if [[ -f "$(casedir "$P" "$F")/DONE" ]]; then
    echo "[skip] $P $F (DONE)"; pathlines "$(casedir "$P" "$F")"; continue
  fi
  if (( running >= JOBS )); then wait -n; running=$(( running - 1 )); fi
  job "$P" "$F" &
  running=$(( running + 1 ))
done
wait
echo "[$(date '+%a %d %H:%M')] all cases ended; waiting for the pathlines"
while pgrep -u "$USER" -f "$EXEPL -pl_case $OUT" > /dev/null; do sleep 60; done
echo "[$(date '+%a %d %H:%M')] postprocessing"
python3 "$POST" "$OUT/3D" > "$LOG/postprocess.log" 2>&1 || echo "postprocess failed: $LOG/postprocess.log"
python3 "$FIGS" "$OUT/3D" --no-render > "$LOG/pathlines_figures.log" 2>&1 || echo "pathlines figures failed: $LOG/pathlines_figures.log"
status
echo "[$(date '+%a %d %H:%M')] done. Figures with ParaView renders: python3 $FIGS $OUT/3D --pvpython <pvpython>"
