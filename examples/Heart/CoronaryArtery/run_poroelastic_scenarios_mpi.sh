#!/usr/bin/env bash
#
# MPI runner of CoronaryArtery_Explicit_PoroElastic_MPI for every
# (geometry, scenario) pair: by default 2 geometries x 3 scenarios = 6 cases,
# all at once with RANKS_PER_CASE MPI ranks each (6 x 30 = 180 processes).
# BATCH=geometry runs one geometry at a time instead (3 cases together, the
# next geometry when they finish).  Every case gets its own block of cores
# (taskset on the mpirun, ranks unbound inside the block) when the machine
# has enough cores, so concurrent mpiruns never bind to the same cores.
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
#   RANKS_PER_CASE MPI ranks per case         (default: 30)
#   BATCH          all: every case at once; geometry: one geometry at a time
#                  (default: all)
#   CYCLES         cardiac cycles per case    (default: 4)
#   DT             time step [s]              (default: 1e-3)
#   OUTPUT_EVERY   XDMF snapshot stride       (default: 10 steps)
#   GEOMETRIES     geometries                 (default: nostenosis stenosis80)
#   SCENARIOS      scenarios of each geometry (default: healthy aortic_stenosis hyperviscosity)
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
RANKS_PER_CASE="${RANKS_PER_CASE:-30}"
BATCH="${BATCH:-all}"
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
NGEOM=${#GEOMS[@]}
(( RANKS_PER_CASE >= 1 )) || { echo "RANKS_PER_CASE must be >= 1"; exit 1; }
NSTEPS=$(python3 -c "print(round($CYCLES * 0.85 / $DT))")

# Batches: lists of "geometry/scenario" run together.
BATCHES=()
case "$BATCH" in
  all)
    list=""
    for g in "${GEOMS[@]}"; do for sc in "${SCENS[@]}"; do list+="$g/$sc "; done; done
    BATCHES+=("$list")
    PER_BATCH=$(( NGEOM * NSCEN ))
    ;;
  geometry)
    for g in "${GEOMS[@]}"; do
      list=""
      for sc in "${SCENS[@]}"; do list+="$g/$sc "; done
      BATCHES+=("$list")
    done
    PER_BATCH=$NSCEN
    ;;
  *) echo "BATCH must be 'all' or 'geometry'"; exit 1 ;;
esac

# Launcher flavour: no binding inside the core block of each case.
if "$MPIRUN" --version 2>&1 | grep -qi "open mpi\|open-rte\|prrte"; then
  BIND=(--bind-to none)
else
  BIND=(-bind-to none)
fi
NCORES=$(nproc --all)
HAVE_TASKSET=0
if command -v taskset > /dev/null && (( NCORES >= PER_BATCH * RANKS_PER_CASE )); then
  HAVE_TASKSET=1
fi

mkdir -p "$OUT"
echo "Rodin: $RODIN"
echo "Batches: ${#BATCHES[@]} ($BATCH), $PER_BATCH case(s) each x $RANKS_PER_CASE rank(s) = $(( PER_BATCH * RANKS_PER_CASE )) processes on $NCORES cores"
(( PER_BATCH * RANKS_PER_CASE <= NCORES )) || echo "Warning: more processes than cores (oversubscribed, no pinning)"
echo "Steps: $NSTEPS (dt=$DT, $CYCLES cycles), XDMF every $OUTPUT_EVERY   Output: $OUT"

status=0
for b in "${!BATCHES[@]}"; do
  read -r -a CASES <<< "${BATCHES[$b]}"
  echo "[$(date '+%F %T')] batch $(( b + 1 ))/${#BATCHES[@]}: ${CASES[*]}"
  PIDS=()
  for i in "${!CASES[@]}"; do
    geom="${CASES[$i]%%/*}"
    name="${CASES[$i]#*/}"
    prefix=$(mesh_prefix "$geom")
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
    if (( HAVE_TASKSET )); then where="cores $first-$last"; else where="unpinned"; fi
    echo "  [$geom/$name] pid ${PIDS[-1]}  $where  log $dir/run.log"
  done

  for i in "${!PIDS[@]}"; do
    if wait "${PIDS[$i]}"; then
      echo "  [${CASES[$i]}] finished"
    else
      echo "  [${CASES[$i]}] FAILED (see $OUT/${CASES[$i]}/run.log)"
      status=1
    fi
  done
  echo "[$(date '+%F %T')] batch $(( b + 1 )) done"
done
exit $status
