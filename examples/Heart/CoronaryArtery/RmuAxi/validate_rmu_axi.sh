#!/usr/bin/env bash
# Validation of RmuAxi against the exact solution for Carreau-Yasuda flow in a
# straight tube, and verification of the R_mu outlet on that tube.
#
#   1. steady:    full straight tube (R_1 = 0.1 mm, L = 40 R_1) at the steady
#                 pressure Pbar of the reference rate gamma_0, on three meshes;
#                 the flow and the developed wall shear stress are compared with
#                 the WRMS solution at the same pressure (validate_rmu_axi.py).
#   2. pulsatile: the same tube, full and truncated at 20 R_1 with the single
#                 (R_mu) and rest (R_0) outlets, under Pbar [1 + 1.2 sin]. The
#                 removed part has one radius, so R_mu is exact in the
#                 quasi-steady limit and its error measures the implementation.
#
# Usage (healthy blood of Table 1, a = 1.8773):
#   BUILD=~/rodin-cmbpe2/build NP_JOBS=12 ./validate_rmu_axi.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
BUILD="${BUILD:-$REPO/build}"
EXE="${EXE:-$BUILD/examples/Heart/CoronaryArtery/RmuAxi/RmuAxi}"
OUT="${OUT:-$PWD/rmu_axi_validation}"
NP_JOBS="${NP_JOBS:-9}"                 # independent runs in parallel, one core each
GAMMAS="${GAMMAS:-30 100 400}"
SCALE="${SCALE:-1e-4}"
TUBES="1:10,1:30"                       # one radius; the step at x = 10 is fictitious
TRUNC=20
REDUCED="1:20"
CY="-rmu_cy_mu0 0.0712616 -rmu_cy_muinf 0.0041 -rmu_cy_lambda 9.2780 -rmu_cy_n 0.3572 -rmu_cy_a 1.8773"

[[ -x "$EXE" ]] || { echo "RmuAxi not found at $EXE"; exit 1; }
mkdir -p "$OUT/mesh"
# Meshes h, h/2, h/4: radial intervals per radius and axial spacing / radius.
for m in "8 1.0" "16 0.5" "32 0.25"; do
  set -- $m
  python3 "$SCRIPT_DIR/make_rmu_axi_mesh.py" --tubes "$TUBES" --truncation "$TRUNC" \
    --prefix "$OUT/mesh/tube$1" --nr "$1" --nmin "$1" --cx "$2"
done

jobs_file="$OUT/jobs.txt"; : > "$jobs_file"
common="-rmu_full_tubes $TUBES -rmu_reduced_tubes $REDUCED -rmu_truncation $TRUNC -rmu_scale $SCALE -rmu_period 0.85 $CY"
for g in $GAMMAS; do
  # 1. steady, three meshes (ramp over one cycle, then one steady cycle)
  for nr in 8 16 32; do
    d="$OUT/steady/g$g/nr$nr"; mkdir -p "$d"
    echo "$EXE -rmu_mesh $OUT/mesh/tube${nr}_full.mesh -rmu_output_prefix $d -rmu_closure none -rmu_gamma_ref $g -rmu_amplitude 0 -rmu_steps_per_cycle 50 -rmu_cycles 2 $common > $d/log.txt 2>&1" >> "$jobs_file"
  done
  # 2. pulsatile, medium mesh
  for c in none single rest; do
    d="$OUT/pulsatile/g$g/$c"; mkdir -p "$d"
    mesh="$OUT/mesh/tube16_full.mesh"; [[ "$c" != none ]] && mesh="$OUT/mesh/tube16_trunc.mesh"
    echo "$EXE -rmu_mesh $mesh -rmu_output_prefix $d -rmu_closure $c -rmu_gamma_ref $g -rmu_amplitude 1.2 -rmu_steps_per_cycle 200 -rmu_cycles 4 $common > $d/log.txt 2>&1" >> "$jobs_file"
  done
done
echo "$(wc -l < "$jobs_file") runs, $NP_JOBS at a time ..."
xargs -P "$NP_JOBS" -I{} bash -c '{}' < "$jobs_file"
python3 "$SCRIPT_DIR/validate_rmu_axi.py" "$OUT"
