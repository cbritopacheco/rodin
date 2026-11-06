#!/bin/bash
# Builds Rodin with Kokkos in a rodin-kokkos-ci container and runs the Kokkos
# unit tests, plus the cuDSS or RocSolverRF test of the backend if any. The
# MPI tests run on 1, 2 and 4 processes, all on the one GPU of the container.
# The Kokkos examples are built, not run.
#
# Usage: ci-test.sh [--build-only] [source-dir] [build-dir]
# --build-only builds the tests and examples without running them, on a host
# without GPU.
# Environment: JOBS (default: nproc), RODIN_KOKKOS_BACKEND (set by the image).
set -uo pipefail

build_only=0
if [[ ${1:-} == --build-only ]]; then
  build_only=1; shift
fi
SRC=$(realpath "${1:-.}")
BUILD=${2:-$SRC/build-kokkos-ci}
JOBS=${JOBS:-$(nproc)}
LOGS=$BUILD/test-logs

options=(
  -DCMAKE_BUILD_TYPE=Release
  -DRODIN_BUILD_UNIT_TESTS=ON -DRODIN_BUILD_EXAMPLES=ON
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
  sycl)
    options+=(-DCMAKE_CXX_COMPILER=icpx) ;;
  *)
    echo "RODIN_KOKKOS_BACKEND must be cuda, hip or sycl" >&2; exit 2 ;;
esac

# A build that does not run needs no optimization: at -O0 on the host and the
# device, the heaviest test compiles in 4 min instead of 21 on CUDA. NVVM
# (cicc) and ptxas take most of the time at -O3. On HIP, -O0 covers both.
if [[ $build_only -eq 1 ]]; then
  options+=("-DCMAKE_CXX_FLAGS_RELEASE=-O0 -DNDEBUG")
  if [[ $RODIN_KOKKOS_BACKEND == cuda ]]; then
    export NVCC_APPEND_FLAGS="${NVCC_APPEND_FLAGS:-} -Xcicc -O0 -Xptxas -O0"
  fi
fi

# Without a GPU the NVIDIA runtime does not mount the driver, and the linker
# needs libcuda.so.1 to resolve the libraries that depend on it: the stub
# stands in for it during the build only.
if [[ $build_only -eq 1 && $RODIN_KOKKOS_BACKEND == cuda ]]; then
  stub=$(mktemp -d)
  ln -s "${CUDA_HOME:-/usr/local/cuda}/lib64/stubs/libcuda.so" "$stub/libcuda.so.1"
  export LD_LIBRARY_PATH=$stub:${LD_LIBRARY_PATH:-}
fi

cmake -S "$SRC" -B "$BUILD" -G Ninja "${options[@]}" || exit 1

help=$(cmake --build "$BUILD" --target help)
mapfile -t targets < <(grep -oE '^Rodin(Kokkos[A-Za-z]*|CuDSS|RocSolverRF)Test' <<<"$help" | sort -u)
if [[ ${#targets[@]} -eq 0 ]]; then
  echo "no Kokkos test target" >&2; exit 1
fi
# The examples of examples/Kokkos that this configuration defines.
mapfile -t examples < <(grep -rhoE 'add_executable\([A-Za-z0-9_]+' "$SRC/examples/Kokkos" \
  | cut -d'(' -f2 | sort -u | grep -xFf <(grep -oE '^[A-Za-z0-9_]+' <<<"$help"))
echo "Building ${#targets[@]} tests and ${#examples[@]} examples: ${targets[*]} ${examples[*]}"
cmake --build "$BUILD" -j"$JOBS" --target "${targets[@]}" "${examples[@]}" || exit 1
[[ $build_only -eq 1 ]] && exit 0

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
