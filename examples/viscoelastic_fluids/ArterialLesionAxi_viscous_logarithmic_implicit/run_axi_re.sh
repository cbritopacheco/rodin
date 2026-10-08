#!/usr/bin/env bash
# =============================================================================
#  run_axi_re.sh -- viscoelastic (sPTT) pulsatile cases of ArterialLesionAxi at
#  ONE Reynolds number, several Wi, several geometries. One terminal per Re:
#
#     ./run_axi_re.sh 100          # terminal 1
#     ./run_axi_re.sh 300          # terminal 2
#     ./run_axi_re.sh 600          # terminal 3
#
#  Each terminal runs its cases one after another with NP MPI ranks (30 by
#  default), so three terminals need 3 x NP cores: check `nproc` first and
#  lower NP if the machine has fewer, otherwise the three sweeps slow each
#  other down.
#
#  Output, one folder per case:
#     $OUT/Re<Re>/<geometry>/Wi<Wi>/
#        run.log      solver output        case.csv    one row per time step
#        case.xdmf/h5 fields               params.txt  exact command line
#        DONE | FAILED                     (DONE cases are skipped on a re-run)
#     $OUT/Re<Re>/summary.csv              status, steps/cycle and minutes per case
#
#  Options (environment variables, all optional):
#     WI_LIST="1 2 4"      Weissenberg numbers; add 0.01 for the Newtonian limit
#     GEOS="S75 S50 H A200"   geometries, run in this order
#     NP=30                MPI ranks per case
#     STEPS=<n>            steps per cycle; default 400 for Re <= 300, 800 above
#                          (same throat CFL across Re, see the paper notes)
#     CYCLES=3  WO=4  AMP=0.5  PTT_EPS=0.1  OUTPUT_EVERY=40  MESH_LEVEL=medium
#     RETRY=1              a diverged case is re-run once with 2x STEPS into
#                          <case>_dt2 (set RETRY=0 to disable)
#     OUT=$HOME/Results/axi   MPIRUN=mpirun   EXE=... MESHDIR=...
#
#  Usage: ./run_axi_re.sh <Re> [geometry ...]     (geometries override GEOS)
# =============================================================================
set -u

if [[ $# -lt 1 ]]; then
  sed -n '2,32p' "$0"; exit 1
fi
RE="$1"; shift

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"

EXE="${EXE:-$REPO/build/examples/viscoelastic_fluids/ArterialLesionAxi_viscous_logarithmic_implicit/ArterialLesionAxi_viscous_logarithmic_implicit}"
MESHDIR="${MESHDIR:-$REPO/resources/examples/viscoelastic_fluids}"
OUT="${OUT:-$HOME/Results/axi}/Re${RE}"
MPIRUN="${MPIRUN:-mpirun}"
NP="${NP:-30}"
WI_LIST="${WI_LIST:-1 2 4}"
GEOS="${GEOS:-S75 S50 H A200}"
[[ $# -gt 0 ]] && GEOS="$*"
WO="${WO:-4}"
AMP="${AMP:-0.5}"
CYCLES="${CYCLES:-3}"
if [[ -z "${STEPS:-}" ]]; then
  STEPS=$(( RE > 300 ? 800 : 400 ))
fi
OUTPUT_EVERY="${OUTPUT_EVERY:-40}"
PTT_EPS="${PTT_EPS:-0.1}"
MESH_LEVEL="${MESH_LEVEL:-medium}"
RETRY="${RETRY:-1}"
NEWTON_ITS="${NEWTON_ITS:-8}"

if [[ ! -x "$EXE" ]]; then
  echo "Executable not found: $EXE"
  echo "  cmake -S $REPO -B $REPO/build && cmake --build $REPO/build -j --target ArterialLesionAxi_viscous_logarithmic_implicit"
  exit 1
fi
if command -v nproc > /dev/null; then
  cores=$(nproc)
elif command -v sysctl > /dev/null; then
  cores=$(sysctl -n hw.ncpu)
else
  cores=0
fi
if [[ $cores -gt 0 && $NP -gt $cores ]]; then
  echo "WARNING: NP=$NP ranks on a machine with $cores cores (oversubscribed). Set NP=<cores> or less."
fi
mkdir -p "$OUT"

# Keep a laptop awake (macOS only); harmless elsewhere.
if [[ -z "${CAFFEINATED:-}" ]] && command -v caffeinate > /dev/null; then
  export CAFFEINATED=1
  exec caffeinate -ims "$0" "$RE" "$@"
fi

SUMMARY="$OUT/summary.csv"
[[ -f "$SUMMARY" ]] || echo "geometry,Re,Wi,steps_per_cycle,status,minutes,folder" > "$SUMMARY"

run_case() {   # geometry Wi steps dir
  local geo="$1" wi="$2" steps="$3" dir="$4"
  local mesh="$MESHDIR/${geo}_axi_${MESH_LEVEL}.mesh"

  if [[ -f "$dir/DONE" ]]; then
    echo "[skip] $geo Re=$RE Wi=$wi steps=$steps (DONE)"; return 0
  fi
  if [[ ! -f "$mesh" ]]; then
    echo "[miss] $mesh not found"; return 2
  fi
  mkdir -p "$dir"; rm -f "$dir/FAILED"

  local cmd=("$MPIRUN" -n "$NP" "$EXE"
    -al_mesh "$mesh" -al_xdmf case -al_csv case.csv
    -al_re "$RE" -al_wo "$WO" -al_amplitude "$AMP" -al_wi "$wi" -al_ptt_epsilon "$PTT_EPS"
    -al_cycles "$CYCLES" -al_steps_per_cycle "$steps" -al_output_every "$OUTPUT_EVERY"
    -al_conformation_its "$NEWTON_ITS")
  printf '%q ' "${cmd[@]}" > "$dir/params.txt"; echo >> "$dir/params.txt"

  echo "[run ] $(date '+%a %d %H:%M')  $geo Re=$RE Wi=$wi steps/cycle=$steps  NP=$NP  ->  $dir"
  local t0=$SECONDS
  ( cd "$dir" && "${cmd[@]}" > run.log 2>&1 )
  local status=$?
  local minutes=$(( (SECONDS - t0) / 60 ))
  if [[ $status -eq 0 ]]; then
    touch "$dir/DONE"
    echo "$geo,$RE,$wi,$steps,DONE,$minutes,$dir" >> "$SUMMARY"
    echo "[done] $geo Re=$RE Wi=$wi in $minutes min"
  else
    echo "exit $status" > "$dir/FAILED"
    echo "$geo,$RE,$wi,$steps,FAILED,$minutes,$dir" >> "$SUMMARY"
    echo "[FAIL] $geo Re=$RE Wi=$wi (exit $status, $minutes min): $(grep -m1 -E 'diverged|Fatal|Error' "$dir/run.log")"
  fi
  return $status
}

echo "Re=$RE  geometries: $GEOS  Wi: $WI_LIST  steps/cycle: $STEPS  NP=$NP  ->  $OUT"
for geo in $GEOS; do
  for wi in $WI_LIST; do
    dir="$OUT/$geo/Wi${wi}"
    run_case "$geo" "$wi" "$STEPS" "$dir"
    st=$?
    if [[ $st -ne 0 && $st -ne 2 && "$RETRY" == "1" ]]; then
      echo "[retry] $geo Wi=$wi with $(( 2 * STEPS )) steps/cycle"
      run_case "$geo" "$wi" "$(( 2 * STEPS ))" "${dir}_dt2"
    fi
  done
done
echo "Re=$RE finished $(date). Summary: $SUMMARY"
