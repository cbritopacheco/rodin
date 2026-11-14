#!/usr/bin/env bash
#
# MPI runner of CoronaryArtery_Explicit_PoroElastic_MPI for every
# (geometry, scenario) pair: by default 2 geometries x 3 scenarios = 6 cases,
# one batch per geometry.  The cases of a batch run together, each with
# RANKS_PER_CASE MPI ranks (3 x 50 = 150 processes); the next geometry starts
# when the batch finishes.  Every case gets its own block of cores
# (taskset on the mpirun, ranks unbound inside the block), so concurrent
# mpiruns never bind to the same cores.
#
# The meshes are read from meshes/*.mesh.gz next to this script (plain git,
# no LFS) and decompressed into $OUT/meshes.
#
# Usage (from anywhere):
#   examples/Heart/CoronaryArtery/run_poroelastic_scenarios_mpi.sh
#   RANKS_PER_CASE=32 GEOMETRIES="stenosis80" SCENARIOS="healthy lvh" \
#     examples/Heart/CoronaryArtery/run_poroelastic_scenarios_mpi.sh
#
# Environment:
#   RODIN          repository root            (default: three levels above this script)
#   BUILD          build directory            (default: $RODIN/build)
#   MPIRUN         MPI launcher               (default: mpirun)
#   RANKS_PER_CASE MPI ranks per case         (default: 50)
#   CYCLES         cardiac cycles per case    (default: 4)
#   DT             time step [s]              (default: 1e-3)
#   OUTPUT_EVERY   XDMF snapshot stride       (default: 10 steps)
#   GEOMETRIES     batches, run one after another (default: nostenosis stenosis80)
#   SCENARIOS      cases of each batch, run together (default: healthy aortic_stenosis hyperviscosity)
#                  available: healthy, aortic_stenosis, hyperviscosity, lvh
#   OUT            output root                (default: $RODIN/results/coronary_poroelastic_mpi)
#   BUILD_FIRST=1  rebuild the target before running
#   MPI_EXTRA      extra launcher options
#   EXTRA          extra PETSc options passed to every case, e.g. the direct
#                  fluid solver:
#                  EXTRA="-coronary_fluid_ksp_type preonly -coronary_fluid_pc_type lu
#                         -coronary_fluid_pc_factor_mat_solver_type mumps"
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RODIN="${RODIN:-$(cd "$SCRIPT_DIR/../../.." && pwd)}"
BUILD="${BUILD:-$RODIN/build}"
MPIRUN="${MPIRUN:-mpirun}"
RANKS_PER_CASE="${RANKS_PER_CASE:-50}"
CYCLES="${CYCLES:-4}"
DT="${DT:-1e-3}"
OUTPUT_EVERY="${OUTPUT_EVERY:-10}"
GEOMETRIES="${GEOMETRIES:-nostenosis stenosis80}"
SCENARIOS="${SCENARIOS:-healthy aortic_stenosis hyperviscosity}"
OUT="${OUT:-$RODIN/results/coronary_poroelastic_mpi}"
MPI_EXTRA="${MPI_EXTRA:-}"
EXTRA="${EXTRA:-}"

TARGET=CoronaryArtery_Explicit_PoroElastic_MPI
EXE="$BUILD/examples/Heart/$TARGET"
MESH_SRC="$SCRIPT_DIR/meshes"
MESH_DIR="$OUT/meshes"

if [[ "${BUILD_FIRST:-0}" == "1" ]]; then
  cmake --build "$BUILD" --target "$TARGET" -j 16
fi
[[ -x "$EXE" ]] || {
  echo "Missing executable: $EXE"
  echo "Configure with -DRODIN_USE_MPI=ON and build the target $TARGET (or BUILD_FIRST=1)."
  exit 1
}

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
(( RANKS_PER_CASE >= 1 )) || { echo "RANKS_PER_CASE must be >= 1"; exit 1; }
NSTEPS=$(python3 -c "print(round($CYCLES * 0.85 / $DT))")

# Launcher flavour: no binding inside the core block of each case.
if "$MPIRUN" --version 2>&1 | grep -qi "open mpi\|open-rte\|prrte"; then
  BIND=(--bind-to none)
else
  BIND=(-bind-to none)
fi
HAVE_TASKSET=0
if command -v taskset > /dev/null && (( $(nproc --all) >= NSCEN * RANKS_PER_CASE )); then
  HAVE_TASKSET=1
fi

mkdir -p "$OUT"
echo "Rodin: $RODIN"
echo "Batches (one after another): ${GEOMS[*]}"
echo "Cases per batch (together):  ${SCENS[*]}  |  $RANKS_PER_CASE rank(s) each, $(( NSCEN * RANKS_PER_CASE )) per batch"
echo "Steps: $NSTEPS (dt=$DT, $CYCLES cycles), XDMF every $OUTPUT_EVERY   Output: $OUT"

status=0
for geom in "${GEOMS[@]}"; do
  prefix=$(mesh_prefix "$geom")
  echo "[$(date '+%F %T')] batch $geom"
  PIDS=()
  for i in "${!SCENS[@]}"; do
    name="${SCENS[$i]}"
    dir="$OUT/$geom/$name"
    mkdir -p "$dir"
    first=$(( i * RANKS_PER_CASE ))
    last=$(( first + RANKS_PER_CASE - 1 ))
    pin=()
    (( HAVE_TASKSET )) && pin=(taskset -c "$first-$last")

    (
      cd "$dir"
      export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
      # shellcheck disable=SC2086
      exec "${pin[@]}" "$MPIRUN" "${BIND[@]}" -n "$RANKS_PER_CASE" $MPI_EXTRA "$EXE" \
        -coronary_scenario "$name" \
        -coronary_results_dir "$dir" \
        -coronary_fluid_mesh "$MESH_DIR/${prefix}Fluid.mesh" \
        -coronary_solid_mesh "$MESH_DIR/${prefix}Solid.mesh" \
        -coronary_dt "$DT" \
        -coronary_nsteps "$NSTEPS" \
        -coronary_output_every "$OUTPUT_EVERY" \
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
