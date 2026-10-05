#!/usr/bin/env bash
# =============================================================================
#  run_axi_convergence.sh -- mesh convergence of ArterialLesionAxi on the most
#  demanding case of the sweep: S75 at its largest Re and Wi.
#
#  Meshes (make_lesion_mesh.py --axi, refinement factor rf, h ~ 1/rf):
#    S75_axi_coarse.mesh  rf 0.35   ~15k vertices
#    S75_axi_medium.mesh  rf 0.50   ~29k vertices   (the sweep mesh)
#    S75_axi_fine.mesh    rf 0.70   ~56k vertices
#  Same dt for the three levels (spatial error only).
#
#  Output: $OUT/S75_Re<Re>_Wi<Wi>/<level>/{run.log, case.csv, case.xdmf, DONE}
#  Usage:  ./run_axi_convergence.sh            (Re = 600, Wi = 4)
#          RE=300 WI=2 ./run_axi_convergence.sh
# =============================================================================
set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"

EXE="${EXE:-$REPO/build/examples/viscoelastic_fluids/ArterialLesionAxi_viscous_logarithmic_implicit/ArterialLesionAxi_viscous_logarithmic_implicit}"
MESHDIR="${MESHDIR:-$REPO/resources/examples/viscoelastic_fluids}"
MPIRUN="${MPIRUN:-mpirun}"
NP="${NP:-5}"
RE="${RE:-600}"
WI="${WI:-4}"
WO="${WO:-4}"
AMP="${AMP:-0.5}"
CYCLES="${CYCLES:-3}"
STEPS="${STEPS:-400}"
OUTPUT_EVERY="${OUTPUT_EVERY:-40}"
PTT_EPS="${PTT_EPS:-0.1}"
OUT="${OUT:-$HOME/Results/axi_convergence}/S75_Re${RE}_Wi${WI}"

if [[ ! -x "$EXE" ]]; then
  echo "Executable not found: $EXE (build target ArterialLesionAxi_viscous_logarithmic_implicit)"
  exit 1
fi
if [[ -z "${CAFFEINATED:-}" ]] && command -v caffeinate > /dev/null; then
  export CAFFEINATED=1
  exec caffeinate -ims "$0" "$@"
fi

for level in coarse medium fine; do
  dir="$OUT/$level"
  [[ -f "$dir/DONE" ]] && { echo "[skip] $level"; continue; }
  mkdir -p "$dir"
  echo "[run ] $(date '+%a %H:%M')  S75 $level  Re=$RE Wi=$WI"
  t0=$SECONDS
  ( cd "$dir" && "$MPIRUN" -n "$NP" "$EXE" \
      -al_mesh "$MESHDIR/S75_axi_${level}.mesh" -al_xdmf case -al_csv case.csv \
      -al_re "$RE" -al_wo "$WO" -al_amplitude "$AMP" -al_wi "$WI" -al_ptt_epsilon "$PTT_EPS" \
      -al_cycles "$CYCLES" -al_steps_per_cycle "$STEPS" -al_output_every "$OUTPUT_EVERY" \
      > run.log 2>&1 ) && touch "$dir/DONE" || echo "[FAIL] $level -- see $dir/run.log"
  echo "       $(( (SECONDS - t0) / 60 )) min"
done
