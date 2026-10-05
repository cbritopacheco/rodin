#!/usr/bin/env bash
# =============================================================================
#  run_axi_sweep.sh -- pulsatile (S, Re, Wi) sweep of ArterialLesionAxi.
#
#  One folder per case:  $OUT/<geometry>/Re<Re>_Wi<Wi>/
#    run.log        full solver output
#    case.csv       one row per time step (q, dp, max|u|, max|sigma|, ...)
#    case.xdmf/h5   fields every $OUTPUT_EVERY steps and at each cycle end
#    params.txt     the exact command line
#    DONE / FAILED  status marker; DONE cases are skipped on a re-run, so the
#                   script can be stopped (Ctrl-C) and relaunched at any time.
#
#  Usage (from anywhere):
#    ./run_axi_sweep.sh check          # 1-2 min sanity run, straight pipe
#    ./run_axi_sweep.sh                # tiers 1 and 2 (default)
#    ./run_axi_sweep.sh 1 2 3 4        # choose tiers
#    OUT=~/Results/axi NP=5 ./run_axi_sweep.sh 1
#
#  Tiers (each loops Wi in the order below, Newtonian-like reference first):
#    1: Re = 300           S75, S50, S25, H
#    2: Re = 100, 600      S50, S75
#    3: Re = 100, 300, 600 S25
#    4: Re = 300           A125, A150, A200  (aneurysms)
# =============================================================================
set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"

EXE="${EXE:-$REPO/build/examples/viscoelastic_fluids/ArterialLesionAxi_viscous_logarithmic_implicit/ArterialLesionAxi_viscous_logarithmic_implicit}"
MESHDIR="${MESHDIR:-$REPO/resources/examples/viscoelastic_fluids}"
OUT="${OUT:-$HOME/Results/axi_sweep}"
MPIRUN="${MPIRUN:-mpirun}"
NP="${NP:-5}"

# ---- fixed parameters -------------------------------------------------------
WO="${WO:-4}"                    # Womersley number
AMP="${AMP:-0.5}"                # q(t) = 1 + AMP sin(omega t)
CYCLES="${CYCLES:-3}"            # cycle 1 is the ramp; cycles 2-3 developed
STEPS="${STEPS:-400}"            # steps per cycle
OUTPUT_EVERY="${OUTPUT_EVERY:-40}"   # XDMF every 40 steps = 10 frames/cycle
PTT_EPS="${PTT_EPS:-0.1}"        # sPTT epsilon (0 = Oldroyd-B)
MESH_LEVEL="${MESH_LEVEL:-medium}"

# Wi = lambda Ubar/D; 0.01 is the Newtonian limit (eta_0) of the same fluid.
WI_LIST="${WI_LIST:-0.01 0.5 1 2 4}"

# ---- checks -----------------------------------------------------------------
if [[ ! -x "$EXE" ]]; then
  echo "Executable not found: $EXE"
  echo "Build it first:"
  echo "  cmake -S $REPO -B $REPO/build"
  echo "  cmake --build $REPO/build -j --target ArterialLesionAxi_viscous_logarithmic_implicit"
  exit 1
fi
mkdir -p "$OUT"

# Keep the Mac awake for the whole sweep (lid open, on power).
if [[ -z "${CAFFEINATED:-}" ]] && command -v caffeinate > /dev/null; then
  export CAFFEINATED=1
  exec caffeinate -ims "$0" "$@"
fi

SUMMARY="$OUT/summary.csv"
[[ -f "$SUMMARY" ]] || echo "geometry,Re,Wi,status,minutes,folder" > "$SUMMARY"

run_case() {   # geometry Re Wi
  local geo="$1" re="$2" wi="$3"
  local mesh="$MESHDIR/${geo}_axi_${MESH_LEVEL}.mesh"
  local dir="$OUT/$geo/Re${re}_Wi${wi}"

  if [[ -f "$dir/DONE" ]]; then
    echo "[skip] $geo Re=$re Wi=$wi (DONE)"
    return 0
  fi
  if [[ ! -f "$mesh" ]]; then
    echo "[miss] $mesh not found"
    return 0
  fi
  mkdir -p "$dir"
  rm -f "$dir/FAILED"

  local cmd=("$MPIRUN" -n "$NP" "$EXE"
    -al_mesh "$mesh" -al_xdmf case -al_csv case.csv
    -al_re "$re" -al_wo "$WO" -al_amplitude "$AMP" -al_wi "$wi"
    -al_ptt_epsilon "$PTT_EPS"
    -al_cycles "$CYCLES" -al_steps_per_cycle "$STEPS" -al_output_every "$OUTPUT_EVERY")
  printf '%q ' "${cmd[@]}" > "$dir/params.txt"; echo >> "$dir/params.txt"

  echo "[run ] $(date '+%a %H:%M')  $geo Re=$re Wi=$wi  ->  $dir"
  local t0=$SECONDS
  ( cd "$dir" && "${cmd[@]}" > run.log 2>&1 )
  local status=$?
  local minutes=$(( (SECONDS - t0) / 60 ))

  if [[ $status -eq 0 ]]; then
    touch "$dir/DONE"
    echo "$geo,$re,$wi,DONE,$minutes,$dir" >> "$SUMMARY"
    echo "[done] $geo Re=$re Wi=$wi in $minutes min"
  else
    echo "exit $status" > "$dir/FAILED"
    echo "$geo,$re,$wi,FAILED,$minutes,$dir" >> "$SUMMARY"
    echo "[FAIL] $geo Re=$re Wi=$wi (exit $status, $minutes min) -- see $dir/run.log"
  fi
}

# ---- sanity check: steady Oldroyd-B Poiseuille in the straight pipe ---------
# dp over L = 35 D must be 32 eta_0 Ubar L/D^2 (up to the inlet stress
# development, < 1%), and qIn = qOut = qTarget.
if [[ "${1:-}" == "check" ]]; then
  dir="$OUT/check_H_steady"
  mkdir -p "$dir"
  ( cd "$dir" && "$MPIRUN" -n "$NP" "$EXE" \
      -al_mesh "$MESHDIR/H_axi_coarse.mesh" -al_xdmf case -al_csv case.csv \
      -al_re 50 -al_wo 4 -al_amplitude 0 -al_wi 0.5 -al_ptt_epsilon 0 \
      -al_cycles 3 -al_steps_per_cycle 50 -al_output_every 0 > run.log 2>&1 )
  echo "exit status: $?   (log: $dir/run.log)"
  grep -E "Newton [0-9]" "$dir/run.log" | head -8
  tail -1 "$dir/case.csv" | awk -F, '{
    eta0 = 3.59e-3; rho = 1060; D = 4e-3; Re = 50; L = 35 * D;
    U = Re * eta0 / (rho * D); dp = 32 * eta0 * U * L / (D * D);
    printf "qTarget=%.6e  qIn=%.6e  qOut=%.6e\n", $3, $4, $5;
    printf "dp computed=%.6f Pa  Poiseuille=%.6f Pa  rel.err=%.3f%%\n", $8, dp, 100 * ($8 - dp) / dp;
    printf "max|u|/Ubar=%.4f (expected 2)\n", $9 / U }'
  exit 0
fi

# ---- sweep --------------------------------------------------------------
TIERS=("$@")
[[ ${#TIERS[@]} -eq 0 ]] && TIERS=(1 2)

for tier in "${TIERS[@]}"; do
  case "$tier" in
    1) for geo in S75 S50 S25 H; do for wi in $WI_LIST; do run_case "$geo" 300 "$wi"; done; done ;;
    2) for re in 100 600; do for geo in S50 S75; do for wi in $WI_LIST; do run_case "$geo" "$re" "$wi"; done; done; done ;;
    3) for re in 100 300 600; do for wi in $WI_LIST; do run_case S25 "$re" "$wi"; done; done ;;
    4) for geo in A125 A150 A200; do for wi in $WI_LIST; do run_case "$geo" 300 "$wi"; done; done ;;
    *) echo "unknown tier $tier" ;;
  esac
done

echo "Sweep finished $(date). Summary: $SUMMARY"
