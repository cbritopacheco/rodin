#!/usr/bin/env bash
#
# Runs CoronaryArtery_Explicit_PoroElastic for every (geometry, scenario) pair:
# by default 2 geometries x 3 scenarios = 6 cases, in one batch per geometry.
# The scenarios of a batch run concurrently with CORES_PER_CASE cores each
# (3 x 50 = 150 cores); the next geometry starts when the batch finishes.
#
# The meshes are read from meshes/*.mesh.gz next to this script (plain git,
# no LFS) and decompressed into $OUT/meshes.
#
# The executable is sequential (PETSc/MUMPS on one MPI rank): the cores of a
# case are used by the threaded BLAS/OpenMP of the MUMPS factorizations, and
# each case is pinned to its own core block.
#
# Usage (from anywhere):
#   examples/Heart/CoronaryArtery/run_poroelastic_scenarios.sh
#   CORES_PER_CASE=50 CYCLES=4 GEOMETRIES="nostenosis" SCENARIOS="healthy lvh" \
#     examples/Heart/CoronaryArtery/run_poroelastic_scenarios.sh
#
# Environment:
#   RODIN          repository root            (default: two levels above this script's parent)
#   BUILD          build directory            (default: $RODIN/build)
#   CORES_PER_CASE cores per case             (default: 50)
#   CYCLES         cardiac cycles per case    (default: 4)
#   DT             time step [s]              (default: 1e-3)
#   GEOMETRIES     batches, run one after another (default: nostenosis stenosis80)
#   SCENARIOS      cases of each batch, run together (default: healthy aortic_stenosis hyperviscosity)
#                  available: healthy, aortic_stenosis, hyperviscosity, lvh
#   OUT            output root                (default: $RODIN/results/coronary_poroelastic)
#   BUILD_FIRST=1  rebuild the target before running
#   EXTRA          extra PETSc options passed to every case
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RODIN="${RODIN:-$(cd "$SCRIPT_DIR/../../.." && pwd)}"
BUILD="${BUILD:-$RODIN/build}"
CORES_PER_CASE="${CORES_PER_CASE:-50}"
CYCLES="${CYCLES:-4}"
DT="${DT:-1e-3}"
GEOMETRIES="${GEOMETRIES:-nostenosis stenosis80}"
SCENARIOS="${SCENARIOS:-healthy aortic_stenosis hyperviscosity}"
OUT="${OUT:-$RODIN/results/coronary_poroelastic}"
EXTRA="${EXTRA:-}"

EXE="$BUILD/examples/Heart/CoronaryArtery_Explicit_PoroElastic"
MESH_SRC="$SCRIPT_DIR/meshes"
MESH_DIR="$OUT/meshes"

if [[ "${BUILD_FIRST:-0}" == "1" ]]; then
  cmake --build "$BUILD" --target CoronaryArtery_Explicit_PoroElastic -j "$CORES_PER_CASE"
fi

[[ -x "$EXE" ]] || { echo "Missing executable: $EXE (set BUILD or BUILD_FIRST=1)"; exit 1; }
mesh_prefix() {
  case "$1" in
    nostenosis) echo NoStenosis ;;
    stenosis80) echo Stenosis80 ;;
    *) echo "Unknown geometry '$1' (nostenosis, stenosis80)" >&2; exit 1 ;;
  esac
}

mkdir -p "$MESH_DIR"
read -r -a GEOMS <<< "$GEOMETRIES"
read -r -a SCENS <<< "$SCENARIOS"
for g in "${GEOMS[@]}"; do
  prefix=$(mesh_prefix "$g")
  for part in Fluid Solid; do
    gz="$MESH_SRC/$prefix$part.mesh.gz"
    [[ -f "$gz" ]] || { echo "Missing $gz"; exit 1; }
    gzip -dc "$gz" > "$MESH_DIR/$prefix$part.mesh"
  done
done

NSCEN=${#SCENS[@]}
THREADS=$CORES_PER_CASE
(( THREADS >= 1 )) || { echo "CORES_PER_CASE must be >= 1"; exit 1; }
NSTEPS=$(python3 -c "print(round($CYCLES * 0.85 / $DT))")
# Pin each case to its own core block only when the machine has the cores.
HAVE_TASKSET=0
if command -v taskset > /dev/null && (( $(nproc --all) >= NSCEN * THREADS )); then
  HAVE_TASKSET=1
fi

mkdir -p "$OUT"
echo "Rodin: $RODIN"
echo "Batches (one after another): ${GEOMS[*]}"
echo "Cases per batch (together):  ${SCENS[*]}  |  $THREADS core(s) each, $(( NSCEN * THREADS )) per batch"
echo "Steps: $NSTEPS (dt=$DT, $CYCLES cycles)   Output: $OUT"

status=0
for geom in "${GEOMS[@]}"; do
  prefix=$(mesh_prefix "$geom")
  echo "[$(date '+%F %T')] batch $geom"
  PIDS=()
  for i in "${!SCENS[@]}"; do
    name="${SCENS[$i]}"
    dir="$OUT/$geom/$name"
    mkdir -p "$dir"
    first=$(( i * THREADS ))
    last=$(( first + THREADS - 1 ))
    pin=()
    (( HAVE_TASKSET )) && pin=(taskset -c "$first-$last")

    (
      cd "$dir"
      export OMP_NUM_THREADS="$THREADS" OPENBLAS_NUM_THREADS="$THREADS" \
             MKL_NUM_THREADS="$THREADS" BLIS_NUM_THREADS="$THREADS"
      # shellcheck disable=SC2086
      exec "${pin[@]}" "$EXE" \
        -coronary_scenario "$name" \
        -coronary_results_dir "$dir" \
        -coronary_fluid_mesh "$MESH_DIR/${prefix}Fluid.mesh" \
        -coronary_solid_mesh "$MESH_DIR/${prefix}Solid.mesh" \
        -coronary_dt "$DT" \
        -coronary_nsteps "$NSTEPS" \
        $EXTRA
    ) > "$dir/run.log" 2>&1 &
    PIDS+=($!)
    echo "  [$geom/$name] pid ${PIDS[-1]}  cores $first-$last  log $dir/run.log"
  done

  for i in "${!PIDS[@]}"; do
    if wait "${PIDS[$i]}"; then
      echo "  [$geom/${SCENS[$i]}] finished"
    else
      echo "  [$geom/${SCENS[$i]}] FAILED (see $OUT/$geom/${SCENS[$i]}/run.log)"
      status=1
    fi
  done
  echo "[$(date '+%F %T')] batch $geom done"
done
exit $status
