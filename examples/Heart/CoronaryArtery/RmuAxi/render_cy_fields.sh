#!/usr/bin/env bash
# Fields of the steady Carreau-Yasuda straight tube and their volume rendering.
#
# Runs RmuAxi on the full tube of the validation (R_1 = 0.1 mm, L = 40 R_1,
# healthy blood) at gamma_0 = 30, 100, 400 1/s with XDMF output of the final
# steady state, then renders u/U on a cut-away developed segment together with
# the computed, exact Carreau-Yasuda and Poiseuille profiles (render_cy_tube.py).
#
#   pip install --user pyvista scipy h5py          # once
#   BUILD=~/rodin-cmbpe2/build bash render_cy_fields.sh
#
# Options (environment): NR=32|16|8 mesh, OUT=folder, GAMMAS="30 100 400",
# RENDER_ARGS="--scale 3 --spacing 0.015" (passed to render_cy_tube.py),
# SKIP_RUN=1 to render existing fields only.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
BUILD="${BUILD:-$REPO/build}"
EXE="${EXE:-$BUILD/examples/Heart/CoronaryArtery/RmuAxi/RmuAxi}"
OUT="${OUT:-$PWD/rmu_axi_fields}"
NR="${NR:-32}"
GAMMAS="${GAMMAS:-30 100 400}"
RENDER_ARGS="${RENDER_ARGS:---scale 3 --spacing 0.015}"
case "$NR" in 8) CX=1.0 ;; 16) CX=0.5 ;; 32) CX=0.25 ;; *) echo "NR must be 8, 16 or 32"; exit 1 ;; esac

TUBES="1:10,1:30"
CY="-rmu_cy_mu0 0.0712616 -rmu_cy_muinf 0.0041 -rmu_cy_lambda 9.2780 -rmu_cy_n 0.3572 -rmu_cy_a 1.8773"
common="-rmu_full_tubes $TUBES -rmu_reduced_tubes 1:20 -rmu_truncation 20 -rmu_scale 1e-4 -rmu_period 0.85 $CY"

if [[ -z "${SKIP_RUN:-}" ]]; then
  [[ -x "$EXE" ]] || { echo "RmuAxi not found at $EXE (set BUILD or EXE)"; exit 1; }
  mkdir -p "$OUT/mesh"
  python3 "$SCRIPT_DIR/make_rmu_axi_mesh.py" --tubes "$TUBES" --truncation 20 \
    --prefix "$OUT/mesh/tube$NR" --nr "$NR" --nmin "$NR" --cx "$CX"
  for g in $GAMMAS; do
    d="$OUT/g$g"; mkdir -p "$d"
    # ramp over one cycle, one steady cycle; write only the last step
    "$EXE" -rmu_mesh "$OUT/mesh/tube${NR}_full.mesh" -rmu_output_prefix "$d" -rmu_closure none \
      -rmu_gamma_ref "$g" -rmu_amplitude 0 -rmu_steps_per_cycle 50 -rmu_cycles 2 \
      -rmu_output_every 100 $common > "$d/log.txt" 2>&1 &
  done
  echo "running $(echo $GAMMAS | wc -w) cases on the nr$NR mesh (about 15 min for nr32) ..."
  wait
fi

dirs=(); for g in $GAMMAS; do dirs+=("$OUT/g$g"); done
cmd=(python3 "$SCRIPT_DIR/render_cy_tube.py" "${dirs[@]}" --out "$OUT/cy_tube.png" $RENDER_ARGS)
if [[ -z "${DISPLAY:-}" ]] && command -v xvfb-run > /dev/null; then
  xvfb-run -a -s "-screen 0 1800x1560x24" "${cmd[@]}"
else
  "${cmd[@]}" || echo "rendering failed here: copy $OUT to a machine with a display and run ${cmd[*]}"
fi
echo "done: $OUT/cy_tube.png"
