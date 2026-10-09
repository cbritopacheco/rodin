#!/usr/bin/env bash
# =============================================================================
#  run_axi_re.sh -- pulsatile cases of ArterialLesionAxi at ONE Reynolds
#  number, for the fluids of the study, several geometries. One terminal per
#  Re (and, if wanted, one per fluid group: see run_axi_extra.sh).
#
#     ./run_axi_re.sh 100          # terminal 1   (Wo 2.5: coronary)
#     ./run_axi_re.sh 300          # terminal 2   (Wo 4:   carotid / femoral)
#     ./run_axi_re.sh 600          # terminal 3   (Wo 5.5: iliac)
#
#  Fluids (FLUIDS="D R" by default; run_axi_extra.sh does "N CY"):
#     N    Newtonian: the sPTT at Wi = 0.01                         -> <geo>/N/
#     D    weak-elasticity sPTT (beta 0.889, eps 0.1), Wi in WI_LIST -> <geo>/Wi<wi>/
#     CY   Carreau-Yasuda of whole blood (Cho & Kensey 1991), no psi -> <geo>/CY/
#     R    generalised sPTT: one mode (eta_p 50 mPa s, lambda 4.83 s, eps 0.2,
#          De 5.6 at 70 bpm) + Carreau-Yasuda solvent; its steady shear
#          viscosity is exactly CY, so CY is its memory-free twin  -> <geo>/R/
#  Re and Wo are defined with eta_ref = 3.45 mPa s (high-shear viscosity of
#  blood) for CY and R; N and D are dimensionless-equivalent (eta_0 of D).
#
#  Output: $OUT/Re<Re>_Wo<Wo>/<geometry>/<fluid>/   (an existing $OUT/Re<Re>
#  is reused when Wo = 4, so the runs of the first sweep keep their place)
#     run.log  case.csv  case.xdmf/h5  params.txt  DONE | FAILED
#     $OUT/Re<Re>_Wo<Wo>/summary.csv
#
#  Options (environment variables):
#     FLUIDS="D R"         fluids, in this order
#     WI_LIST="2"          Wi of fluid D (use "1 2 4" for the Wi sweep)
#     GEOS="S75 S50 H A200"
#     WO=<n>               Womersley; default by Re: 100 -> 2.5, 300 -> 4, 600 -> 5.5
#     NP=30                MPI ranks per case
#     STEPS=<n>            steps per cycle; default 400 (Re <= 300), 800 above
#     CYCLES=3             cycles for N, D, CY
#     CYCLES_R=8           cycles for R (lambda spans ~6 beats: check the CSV
#                          for cycle-to-cycle periodicity and extend if needed)
#     AMP=0.5  PTT_EPS=0.1 (fluid D)  OUTPUT_EVERY=40  MESH_LEVEL=medium
#     RETRY=1              a diverged case is re-run once with 2x STEPS (<case>_dt2)
#     OUT=$HOME/Results/axi   MPIRUN=mpirun   EXE=...   MESHDIR=...
#
#  Usage: ./run_axi_re.sh <Re> [geometry ...]
# =============================================================================
set -u

if [[ $# -lt 1 ]]; then
  sed -n '2,40p' "$0"; exit 1
fi
RE="$1"; shift

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"

EXE="${EXE:-$REPO/build/examples/viscoelastic_fluids/ArterialLesionAxi_viscous_logarithmic_implicit/ArterialLesionAxi_viscous_logarithmic_implicit}"
MESHDIR="${MESHDIR:-$REPO/resources/examples/viscoelastic_fluids}"
MPIRUN="${MPIRUN:-mpirun}"
NP="${NP:-30}"
FLUIDS="${FLUIDS:-D R}"
WI_LIST="${WI_LIST:-2}"
GEOS="${GEOS:-S75 S50 H A200}"
[[ $# -gt 0 ]] && GEOS="$*"
if [[ -z "${WO:-}" ]]; then
  case "$RE" in
    100) WO=2.5 ;; 300) WO=4 ;; 600) WO=5.5 ;; *) WO=4 ;;
  esac
fi
AMP="${AMP:-0.5}"
CYCLES="${CYCLES:-3}"
CYCLES_R="${CYCLES_R:-8}"
if [[ -z "${STEPS:-}" ]]; then
  STEPS=$(( RE > 300 ? 800 : 400 ))
fi
OUTPUT_EVERY="${OUTPUT_EVERY:-40}"
PTT_EPS="${PTT_EPS:-0.1}"
MESH_LEVEL="${MESH_LEVEL:-medium}"
RETRY="${RETRY:-1}"
NEWTON_ITS="${NEWTON_ITS:-8}"
ETA_REF="${ETA_REF:-0.00345}"

OUTROOT="${OUT:-$HOME/Results/axi}"
if [[ "$WO" == "4" && -d "$OUTROOT/Re${RE}" ]]; then
  OUT="$OUTROOT/Re${RE}"
else
  OUT="$OUTROOT/Re${RE}_Wo${WO}"
fi

if [[ ! -x "$EXE" ]]; then
  echo "Executable not found: $EXE"
  echo "  cmake -S $REPO -B $REPO/build && cmake --build $REPO/build -j --target ArterialLesionAxi_viscous_logarithmic_implicit"
  exit 1
fi
if command -v nproc > /dev/null; then cores=$(nproc)
elif command -v sysctl > /dev/null; then cores=$(sysctl -n hw.ncpu)
else cores=0; fi
if [[ $cores -gt 0 && $NP -gt $cores ]]; then
  echo "WARNING: NP=$NP ranks on a machine with $cores cores (oversubscribed). Set NP=<cores> or less."
fi
mkdir -p "$OUT"

if [[ -z "${CAFFEINATED:-}" ]] && command -v caffeinate > /dev/null; then
  export CAFFEINATED=1
  exec caffeinate -ims "$0" "$RE" "$@"
fi

SUMMARY="$OUT/summary.csv"
[[ -f "$SUMMARY" ]] || echo "geometry,Re,Wo,fluid,Wi,steps_per_cycle,cycles,status,minutes,folder" > "$SUMMARY"

# fluid_args <fluid> <wi>  -> the -al_ options of that fluid (echoed)
fluid_args() {
  case "$1" in
    N)  echo "-al_wi 0.01 -al_ptt_epsilon $PTT_EPS" ;;
    D)  echo "-al_wi $2 -al_ptt_epsilon $PTT_EPS" ;;
    CY) echo "-al_fluid cy -al_eta_ref $ETA_REF" ;;
    R)  echo "-al_fluid gsptt -al_eta_ref $ETA_REF -al_eta_p 0.050 -al_ptt_epsilon 0.2 -al_de 5.635 -al_wi 0 -al_conformation_its 10 -al_inlet_developed_conformation 1" ;;
    *)  echo "unknown fluid $1" >&2; return 1 ;;
  esac
}

run_case() {   # geometry fluid wi steps cycles dir
  local geo="$1" fluid="$2" wi="$3" steps="$4" cycles="$5" dir="$6"
  local mesh="$MESHDIR/${geo}_axi_${MESH_LEVEL}.mesh"

  if [[ -f "$dir/DONE" ]]; then
    echo "[skip] $geo $fluid Wi=$wi Re=$RE Wo=$WO (DONE)"; return 0
  fi
  if [[ ! -f "$mesh" ]]; then
    echo "[miss] $mesh not found"; return 2
  fi
  mkdir -p "$dir"; rm -f "$dir/FAILED"

  local fargs
  fargs=$(fluid_args "$fluid" "$wi") || return 2
  local cmd=("$MPIRUN" -n "$NP" "$EXE"
    -al_mesh "$mesh" -al_xdmf case -al_csv case.csv
    -al_re "$RE" -al_wo "$WO" -al_amplitude "$AMP"
    -al_cycles "$cycles" -al_steps_per_cycle "$steps" -al_output_every "$OUTPUT_EVERY")
  [[ "$fluid" != "R" ]] && cmd+=(-al_conformation_its "$NEWTON_ITS")
  # shellcheck disable=SC2206
  cmd+=($fargs)
  printf '%q ' "${cmd[@]}" > "$dir/params.txt"; echo >> "$dir/params.txt"

  echo "[run ] $(date '+%a %d %H:%M')  $geo $fluid Wi=$wi Re=$RE Wo=$WO steps/cycle=$steps cycles=$cycles NP=$NP -> $dir"
  local t0=$SECONDS
  ( cd "$dir" && "${cmd[@]}" > run.log 2>&1 )
  local status=$?
  local minutes=$(( (SECONDS - t0) / 60 ))
  if [[ $status -eq 0 ]]; then
    touch "$dir/DONE"
    echo "$geo,$RE,$WO,$fluid,$wi,$steps,$cycles,DONE,$minutes,$dir" >> "$SUMMARY"
    echo "[done] $geo $fluid Wi=$wi in $minutes min"
  else
    echo "exit $status" > "$dir/FAILED"
    echo "$geo,$RE,$WO,$fluid,$wi,$steps,$cycles,FAILED,$minutes,$dir" >> "$SUMMARY"
    echo "[FAIL] $geo $fluid Wi=$wi (exit $status, $minutes min): $(grep -m1 -E 'diverged|Fatal|Error' "$dir/run.log")"
  fi
  return $status
}

echo "Re=$RE Wo=$WO  fluids: $FLUIDS  geometries: $GEOS  Wi(D): $WI_LIST  steps/cycle: $STEPS  NP=$NP  ->  $OUT"
for geo in $GEOS; do
  for fluid in $FLUIDS; do
    if [[ "$fluid" == "D" ]]; then wis="$WI_LIST"; else wis="-"; fi
    for wi in $wis; do
      case "$fluid" in
        D) sub="Wi${wi}"; cycles="$CYCLES" ;;
        R) sub="R"; cycles="$CYCLES_R" ;;
        *) sub="$fluid"; cycles="$CYCLES" ;;
      esac
      dir="$OUT/$geo/$sub"
      run_case "$geo" "$fluid" "$wi" "$STEPS" "$cycles" "$dir"
      st=$?
      if [[ $st -ne 0 && $st -ne 2 && "$RETRY" == "1" ]]; then
        echo "[retry] $geo $fluid Wi=$wi with $(( 2 * STEPS )) steps/cycle"
        run_case "$geo" "$fluid" "$wi" "$(( 2 * STEPS ))" "$cycles" "${dir}_dt2"
      fi
    done
  done
done
echo "Re=$RE Wo=$WO finished $(date). Summary: $SUMMARY"
