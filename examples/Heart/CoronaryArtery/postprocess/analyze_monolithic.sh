#!/usr/bin/env bash
#
# Full analysis of the monolithic (Newton-GMRES, RR-5085 eq. 15) coronary FSI
# runs: 3D haemodynamics, velocity cross-sections along the branches and the
# in-depth 0D / coupling analysis.
#
# Usage:
#   examples/Heart/CoronaryArtery/postprocess/analyze_monolithic.sh [RESULTS] [FIGURES]
#     RESULTS  default: $RODIN/results/coronary_monolithic
#     FIGURES  default: $RESULTS/figures
#
# Output (PDF vector + PNG 600 dpi):
#   fig3d_wss_diastole, fig3d_tawss, fig3d_osi    WSS (instantaneous, TAWSS, OSI)
#   fig3d_streamlines                             streamlines seeded in the whole lumen
#   fig3d_wall_motion                             cyclic wall displacement
#   fig3d_sections_<geom>_<scenario>              |u| on cross-sections, main vessel
#                                                 + branches, with locator
#   fig3d_sections_main_<geom>                    main-vessel sections, scenarios as rows
#   fig_lv_poroelastic, fig_coronary_coupling,    0D poroelastic LV and 0D-3D coupling
#   fig_phase_loops, fig_scenario_changes
#   fig_numerics                                  Newton/GMRES its, interface leakage
#   metrics_poroelastic.{csv,tex}, metrics_3d.csv cycle metrics (LaTeX table)
#   diagnostics.{csv,txt}                         periodicity, 3D-0D balance, leakage
#
# Environment:
#   PERIOD=0.85       cardiac period [s] (the LAST complete cycle is analysed)
#   PYTHON=python3
#   SEEDING=volume    streamline seeds: volume (whole lumen) | inlet
#   SEEDS=1500        number of streamline seeds
#   BRANCHES=2        branches cut besides the main vessel
#   STATIONS=4        cross-sections per vessel
#   SECTION_TIME=diastole   diastole | systole (instant of the sections)
#   WIDTH=2400        pixel width of each 3D panel
#   SKIP_3D=1         only the 0D/coupling analysis
#   CMAME_USETEX=1    typeset with the local LaTeX
#
# Python packages: numpy scipy matplotlib h5py pyvista pillow
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RODIN="${RODIN:-$(cd "$HERE/../../../.." && pwd)}"
RESULTS="${1:-$RODIN/results/coronary_monolithic}"
FIGURES="${2:-$RESULTS/figures}"
PERIOD="${PERIOD:-0.85}"
PYTHON="${PYTHON:-python3}"
WIDTH="${WIDTH:-2400}"
SEEDING="${SEEDING:-volume}"
SEEDS="${SEEDS:-1500}"
BRANCHES="${BRANCHES:-2}"
STATIONS="${STATIONS:-4}"
SECTION_TIME="${SECTION_TIME:-diastole}"

[[ -d "$RESULTS" ]] || { echo "No results directory: $RESULTS (pass it as 1st argument)"; exit 1; }
mkdir -p "$FIGURES"

missing=$("$PYTHON" - <<'EOF'
import importlib.util
need = ["numpy", "scipy", "matplotlib", "h5py", "pyvista", "PIL"]
print(" ".join(m for m in need if importlib.util.find_spec(m) is None))
EOF
)
if [[ -n "$missing" ]]; then
  echo "Missing Python packages: $missing"
  echo "  $PYTHON -m pip install numpy scipy matplotlib h5py pyvista pillow"
  [[ "$missing" == *numpy* || "$missing" == *matplotlib* ]] && exit 1
  SKIP_3D=1
fi

echo "Results: $RESULTS"
echo "Figures: $FIGURES"
find "$RESULTS" -mindepth 3 -maxdepth 3 -name 'CoronaryArtery_Explicit_PoroElastic.csv' \
  -printf '  case %h\n' | sed "s#$RESULTS/##" | sort

export PYVISTA_OFF_SCREEN=true
if [[ "${SKIP_3D:-0}" != "1" ]]; then
  echo "[3D] WSS, TAWSS, OSI, streamlines ($SEEDING seeding), wall motion ..."
  "$PYTHON" "$HERE/render_3d.py" "$RESULTS" "$FIGURES" --period "$PERIOD" \
    --width "$WIDTH" --seeding "$SEEDING" --seeds "$SEEDS"
  echo "[3D] centerlines and velocity cross-sections ..."
  "$PYTHON" "$HERE/sections_3d.py" "$RESULTS" "$FIGURES" --period "$PERIOD" \
    --branches "$BRANCHES" --stations "$STATIONS" --time "$SECTION_TIME"
fi

echo "[0D] poroelastic LV, coupling, diagnostics and solver effort ..."
m3d=()
[[ -f "$FIGURES/metrics_3d.csv" ]] && m3d=(--metrics3d "$FIGURES/metrics_3d.csv")
"$PYTHON" "$HERE/plot_poroelastic.py" "$RESULTS" "$FIGURES" --period "$PERIOD" "${m3d[@]}"

echo "Done. Figures (PDF vector + PNG 600 dpi) in $FIGURES:"
ls -1 "$FIGURES"/*.pdf | sed 's#.*/#  #'
echo "Diagnostics: $FIGURES/diagnostics.txt"
