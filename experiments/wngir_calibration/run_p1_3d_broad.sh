#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
runner="$root/experiments/wngir_calibration/run_p1_3d_screen.py"
executable="$root/build-p1-3d-clang19/examples/Geometry/LevelSetWNGIRReconstruction3D"
output_root="$root/tmp/wngir_p1_3d_canonical"

run_stage() {
  local n="$1"
  local lobes="$2"
  local timeout="$3"
  local output="$output_root/n$n"
  local scratch="/tmp/wngir-p1-3d-broad-n$n"
  local common=(
    --n="$n" --lobes="$lobes" --threads=6 --timeout="$timeout"
    --barrier-max-iters=10 --cg-rtol=1e-6 --exe="$executable"
    --out-dir="$output" --scratch="$scratch"
  )

  python3 "$runner" preflight "${common[@]}"
  python3 "$runner" screen "${common[@]}" --steps=20 \
    --kappa-f=1 --kappa-s=1 \
    --kappa-d=1 --mu-hat=90
}

run_stage 5 0,1 600
run_stage 10 0,1,2 1200
run_stage 20 0:4 3600
