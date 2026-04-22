#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
BUILD_DIR="${BUILD_DIR:-${REPO_ROOT}/build/codex}"
RUN_BUILD="1"
RUN_TESTS="1"
RUN_INSTALL_CHECK="0"

usage() {
  cat <<USAGE
Usage: $(basename "$0") [options]

Run repeatable Codex maintenance tasks for Rodin.

Options:
  --build-dir <path>       Build directory (default: ${BUILD_DIR})
  --skip-build             Skip cmake --build
  --skip-tests             Skip ctest label suites
  --installation-check     Run tests/installation/test_installation.sh at the end
  -h, --help               Show this help
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --build-dir)
      BUILD_DIR="$2"
      shift 2
      ;;
    --skip-build)
      RUN_BUILD="0"
      shift
      ;;
    --skip-tests)
      RUN_TESTS="0"
      shift
      ;;
    --installation-check)
      RUN_INSTALL_CHECK="1"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage
      exit 1
      ;;
  esac
done

cd "${REPO_ROOT}"
echo "[maintenance] Repository root: ${REPO_ROOT}"

echo "[maintenance] Updating submodules"
git submodule update --init --recursive

echo "[maintenance] Refreshing CMake cache"
cmake -S "${REPO_ROOT}" -B "${BUILD_DIR}" \
  -DCMAKE_BUILD_TYPE="${BUILD_TYPE:-Debug}" \
  -DRODIN_BUILD_DOC=OFF \
  -DRODIN_USE_MCSS=OFF \
  -DRODIN_BUILD_SRC=ON \
  -DRODIN_BUILD_UNIT_TESTS=ON \
  -DRODIN_BUILD_MANUFACTURED_TESTS=ON \
  -DRODIN_BUILD_EXAMPLES=OFF \
  -DRODIN_BUILD_BENCHMARKS=OFF

if [[ "${RUN_BUILD}" == "1" ]]; then
  echo "[maintenance] Building project"
  cmake --build "${BUILD_DIR}" -j"$(nproc 2>/dev/null || sysctl -n hw.ncpu || echo 4)"
fi

if [[ "${RUN_TESTS}" == "1" ]]; then
  echo "[maintenance] Running unit tests (excluding slow)"
  ctest --test-dir "${BUILD_DIR}/tests" -L unit -LE slow --output-on-failure -j "${CTEST_PARALLEL_LEVEL:-2}"

  echo "[maintenance] Running manufactured tests (excluding slow)"
  ctest --test-dir "${BUILD_DIR}/tests" -L manufactured -LE slow --output-on-failure -j 1
fi

if [[ "${RUN_INSTALL_CHECK}" == "1" ]]; then
  echo "[maintenance] Running installation regression check"
  bash "${REPO_ROOT}/tests/installation/test_installation.sh"
fi

echo "[maintenance] Completed"
