#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="$SCRIPT_DIR/build"

JOBS="${JOBS:-16}"
MESH_FILE="${1:-/home/wjj/Project/stage/levelset_stoke/sphere_init_aligned.mesh}"
MAX_IT="${2:-10}"
I_AXIS="${3:-0}"
J_AXIS="${4:-0}"
DT="${DT:-0.01}"
ALPHA="${ALPHA:-1.0}"
HMAX="${HMAX:-0.28}"

cmake -S "$SCRIPT_DIR" -B "$BUILD_DIR"
cmake --build "$BUILD_DIR" -j "$JOBS"

cd "$BUILD_DIR"
mkdir -p out
OMP_NUM_THREADS="${OMP_NUM_THREADS:-16}" \
  ./LevelSetStokes3DObstacle "$MESH_FILE" "$MAX_IT" "$I_AXIS" "$J_AXIS" "$DT" "$ALPHA" "$HMAX" \
  -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -ksp_converged_reason
