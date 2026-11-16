#!/usr/bin/env bash
#
# Post-processing of the coronary poroelastic FSI cases (healthy vs. 80%
# stenosis x scenarios): CMAME-ready figures of the poroelastic LV, the 0D-3D
# coupling and the 3D haemodynamics (WSS, TAWSS, OSI, streamlines, wall motion).
#
# Usage:
#   examples/Heart/CoronaryArtery/postprocess/compare_poroelastic_results.sh [RESULTS] [FIGURES]
#     RESULTS  default: $RODIN/results/coronary_poroelastic_mpi
#     FIGURES  default: $RESULTS/figures
#
# Environment:
#   PERIOD=0.85   cardiac period [s] (the LAST complete cycle is analysed)
#   PYTHON=python3
#   SKIP_3D=1     only the 0D/coupling figures
#   WIDTH=2400    pixel width of each 3D panel
#   EXAGGERATION  wall-motion warp factor (default: automatic)
#   CMAME_USETEX=1  typeset with the local LaTeX
#
# Python packages: numpy matplotlib h5py pyvista pillow
#   $PYTHON -m pip install numpy matplotlib h5py pyvista pillow
# 3D rendering is off-screen (no display needed).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RODIN="${RODIN:-$(cd "$HERE/../../../.." && pwd)}"
RESULTS="${1:-$RODIN/results/coronary_poroelastic_mpi}"
FIGURES="${2:-$RESULTS/figures}"
PERIOD="${PERIOD:-0.85}"
PYTHON="${PYTHON:-python3}"
WIDTH="${WIDTH:-2400}"

[[ -d "$RESULTS" ]] || { echo "No results directory: $RESULTS"; exit 1; }
mkdir -p "$FIGURES"

missing=$("$PYTHON" - <<'EOF'
import importlib.util
need = ["numpy", "matplotlib", "h5py", "pyvista", "PIL"]
print(" ".join(m for m in need if importlib.util.find_spec(m) is None))
EOF
)
if [[ -n "$missing" ]]; then
  echo "Missing Python packages: $missing"
  echo "  $PYTHON -m pip install numpy matplotlib h5py pyvista pillow"
  [[ "$missing" == *numpy* || "$missing" == *matplotlib* ]] && exit 1
  SKIP_3D=1
fi

echo "Results: $RESULTS"
echo "Figures: $FIGURES"
find "$RESULTS" -mindepth 3 -maxdepth 3 -name 'CoronaryArtery_Explicit_PoroElastic.csv' \
  -printf '  case %h\n' | sed "s#$RESULTS/##" | sort

export PYVISTA_OFF_SCREEN=true
if [[ "${SKIP_3D:-0}" != "1" ]]; then
  echo "[3D] TAWSS/OSI over the last cycle, snapshots and panels ..."
  extra=()
  [[ -n "${EXAGGERATION:-}" ]] && extra+=(--exaggeration "$EXAGGERATION")
  "$PYTHON" "$HERE/render_3d.py" "$RESULTS" "$FIGURES" --period "$PERIOD" \
    --width "$WIDTH" "${extra[@]}"
fi

echo "[0D] poroelastic LV and coupling curves, metrics ..."
m3d=()
[[ -f "$FIGURES/metrics_3d.csv" ]] && m3d=(--metrics3d "$FIGURES/metrics_3d.csv")
"$PYTHON" "$HERE/plot_poroelastic.py" "$RESULTS" "$FIGURES" --period "$PERIOD" "${m3d[@]}"

echo "Done. Figures (PDF vector + PNG 600 dpi) in $FIGURES:"
ls -1 "$FIGURES"/*.pdf | sed 's#.*/#  #'
