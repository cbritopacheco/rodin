#!/usr/bin/env bash
#
# Runs CoronaryArtery_Explicit_PoroElastic for several scenarios at once.
#
# The executable is sequential (PETSc/MUMPS on one MPI rank), so the cores are
# split between concurrent cases: each case is pinned to its own block of
# NPROC / ncases cores and its BLAS/OpenMP threads (used by the MUMPS
# factorizations) are limited to that block.
#
# Usage (from anywhere):
#   examples/Heart/CoronaryArtery/run_poroelastic_scenarios.sh
#   NPROC=50 CYCLES=4 SCENARIOS="healthy aortic_stenosis hyperviscosity" \
#     examples/Heart/CoronaryArtery/run_poroelastic_scenarios.sh
#
# Environment:
#   RODIN      repository root            (default: two levels above this script's parent)
#   BUILD      build directory            (default: $RODIN/build)
#   NPROC      cores to use               (default: 50)
#   CYCLES     cardiac cycles per case    (default: 4)
#   DT         time step [s]              (default: 1e-3)
#   SCENARIOS  cases to run               (default: healthy aortic_stenosis hyperviscosity)
#              available: healthy, aortic_stenosis, hyperviscosity, lvh
#   OUT        output root                (default: $RODIN/results/coronary_poroelastic)
#   BUILD_FIRST=1  rebuild the target before running
#   EXTRA      extra PETSc options passed to every case
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RODIN="${RODIN:-$(cd "$SCRIPT_DIR/../../.." && pwd)}"
BUILD="${BUILD:-$RODIN/build}"
NPROC="${NPROC:-50}"
CYCLES="${CYCLES:-4}"
DT="${DT:-1e-3}"
SCENARIOS="${SCENARIOS:-healthy aortic_stenosis hyperviscosity}"
OUT="${OUT:-$RODIN/results/coronary_poroelastic}"
EXTRA="${EXTRA:-}"

EXE="$BUILD/examples/Heart/CoronaryArtery_Explicit_PoroElastic"
FLUID_MESH="$RODIN/resources/examples/Heart/CoronaryArtery/Stenosis80.mesh"
SOLID_MESH="$RODIN/resources/examples/Heart/CoronaryArtery/Prismatic80.mesh"

if [[ "${BUILD_FIRST:-0}" == "1" ]]; then
  cmake --build "$BUILD" --target CoronaryArtery_Explicit_PoroElastic -j "$NPROC"
fi

[[ -x "$EXE" ]] || { echo "Missing executable: $EXE (set BUILD or BUILD_FIRST=1)"; exit 1; }
for m in "$FLUID_MESH" "$SOLID_MESH"; do
  if [[ ! -f "$m" ]] || [[ $(stat -c %s "$m") -lt 10000 ]]; then
    echo "Mesh $m is missing or an unhydrated Git LFS pointer (run: git lfs pull)"; exit 1
  fi
done

read -r -a CASES <<< "$SCENARIOS"
NCASES=${#CASES[@]}
THREADS=$(( NPROC / NCASES ))
(( THREADS >= 1 )) || { echo "NPROC=$NPROC is smaller than the number of cases"; exit 1; }
NSTEPS=$(python3 -c "print(round($CYCLES * 0.85 / $DT))")
HAVE_TASKSET=0
command -v taskset > /dev/null && HAVE_TASKSET=1

mkdir -p "$OUT"
echo "Rodin: $RODIN"
echo "Cases: ${CASES[*]}  |  $THREADS core(s) each  |  $NSTEPS steps (dt=$DT, $CYCLES cycles)"
echo "Output: $OUT"

PIDS=()
for i in "${!CASES[@]}"; do
  name="${CASES[$i]}"
  dir="$OUT/$name"
  mkdir -p "$dir"
  first=$(( i * THREADS ))
  last=$(( first + THREADS - 1 ))
  pin=()
  (( HAVE_TASKSET )) && pin=(taskset -c "$first-$last")

  # The example resolves ../resources relative to its working directory; the
  # meshes are passed explicitly, so any working directory works.
  (
    cd "$dir"
    export OMP_NUM_THREADS="$THREADS" OPENBLAS_NUM_THREADS="$THREADS" \
           MKL_NUM_THREADS="$THREADS" BLIS_NUM_THREADS="$THREADS"
    # shellcheck disable=SC2086
    exec "${pin[@]}" "$EXE" \
      -coronary_scenario "$name" \
      -coronary_results_dir "$dir" \
      -coronary_fluid_mesh "$FLUID_MESH" \
      -coronary_solid_mesh "$SOLID_MESH" \
      -coronary_dt "$DT" \
      -coronary_nsteps "$NSTEPS" \
      $EXTRA
  ) > "$dir/run.log" 2>&1 &
  PIDS+=($!)
  echo "  [$name] pid ${PIDS[-1]}  cores $first-$last  log $dir/run.log"
done

status=0
for i in "${!PIDS[@]}"; do
  if wait "${PIDS[$i]}"; then
    echo "  [${CASES[$i]}] finished"
  else
    echo "  [${CASES[$i]}] FAILED (see $OUT/${CASES[$i]}/run.log)"
    status=1
  fi
done
exit $status
