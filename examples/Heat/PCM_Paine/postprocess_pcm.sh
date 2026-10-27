#!/usr/bin/env bash
# Heat-flux summary, plots and temperature contours of the PCM thickness sweep.
#
#   bash postprocess_pcm.sh [OUT] [--hours 9 13 17 21] [--cases 5 25 50 | --all-contours]
#
# OUT defaults to ./pcm_paine_runs (the run_pcm_sweep.sh default). Writes
# OUT/post/{summary.csv, summary.md, flujos_dia.png, resumen_espesor.png,
# fraccion_liquida.png, T_contornos_e_XXmm.png}. Needs python3 with numpy,
# matplotlib and h5py (pip install numpy matplotlib h5py).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${1:-$PWD/pcm_paine_runs}"
[ $# -gt 0 ] && shift
python3 -c "import numpy, matplotlib, h5py" 2>/dev/null || {
  echo "Missing python modules: pip install numpy matplotlib h5py" >&2
  exit 1
}
python3 "$HERE/postprocess_pcm.py" "$OUT" "$@"
[ -f "$OUT/post/summary.md" ] && { echo; cat "$OUT/post/summary.md"; }
