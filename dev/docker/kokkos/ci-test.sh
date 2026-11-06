#!/bin/bash
# Builds Rodin with Kokkos in a rodin-kokkos-ci container and runs the Kokkos
# unit tests, plus the cuDSS or RocSolverRF test of the backend. The MPI tests
# run on 1, 2 and 4 processes, all on the one GPU of the container.
#
# Usage: ci-test.sh [source-dir] [build-dir]
# Environment: JOBS (default: nproc), RODIN_KOKKOS_BACKEND (set by the image).
set -uo pipefail

SRC=$(realpath "${1:-.}")
BUILD=${2:-$SRC/build-kokkos-ci}
JOBS=${JOBS:-$(nproc)}
LOGS=$BUILD/test-logs

options=(
  -DCMAKE_BUILD_TYPE=Release
  -DRODIN_BUILD_UNIT_TESTS=ON -DRODIN_BUILD_EXAMPLES=OFF
  -DRODIN_USE_MPI=ON -DRODIN_USE_OPENMP=ON -DRODIN_USE_PETSC=ON
  -DRODIN_USE_KOKKOS=ON -DRODIN_USE_KOKKOSKERNELS=ON -DRODIN_USE_HYPRE=ON
)
case "${RODIN_KOKKOS_BACKEND:-}" in
  cuda)
    options+=(-DRODIN_USE_CUDSS=ON) ;;
  hip)
    # hipcc would add -x hip to every source; the clang++ of ROCm compiles
    # the HIP sources of Rodin itself.
    options+=(-DCMAKE_CXX_COMPILER=/opt/rocm/lib/llvm/bin/clang++
              -DRODIN_USE_ROCSOLVER=ON -DRODIN_USE_UMFPACK=ON) ;;
  *)
    echo "RODIN_KOKKOS_BACKEND must be cuda or hip" >&2; exit 2 ;;
esac

cmake -S "$SRC" -B "$BUILD" -G Ninja "${options[@]}" || exit 1

mapfile -t targets < <(cmake --build "$BUILD" --target help \
  | grep -oE '^Rodin(Kokkos[A-Za-z]*|CuDSS|RocSolverRF)Test' | sort -u)
if [[ ${#targets[@]} -eq 0 ]]; then
  echo "no Kokkos test target" >&2; exit 1
fi
cmake --build "$BUILD" -j"$JOBS" --target "${targets[@]}" || exit 1

mkdir -p "$LOGS"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4} HDF5_USE_FILE_LOCKING=FALSE
status=0
run() {
  local log=$1 dir=$2 rc; shift 2
  (cd "$dir" && timeout 900 "$@") > "$log" 2>&1
  rc=$?
  if [[ $rc -eq 0 ]]; then
    echo "PASS $(basename "$log" .log)"
  else
    echo "FAIL $(basename "$log" .log) (rc=$rc)"; tail -30 "$log"; status=1
  fi
}
for t in "${targets[@]}"; do
  exe=$(find "$BUILD/tests" -type f -executable -name "$t" | head -1)
  if [[ -z $exe ]]; then
    echo "FAIL $t (not built)"; status=1; continue
  fi
  if [[ $t == *MPI* ]]; then
    for np in 1 2 4; do run "$LOGS/$t-np$np.log" "$(dirname "$exe")" mpirun -n $np "$exe"; done
  else
    run "$LOGS/$t.log" "$(dirname "$exe")" "$exe"
  fi
done
exit $status
