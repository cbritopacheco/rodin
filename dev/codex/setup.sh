#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
BUILD_DIR="${REPO_ROOT}/build/codex"
BUILD_TYPE="${BUILD_TYPE:-RelWithDebInfo}"
INSTALL_DEPS="0"
RUN_CONFIGURE="1"

usage() {
  cat <<USAGE
Usage: $(basename "$0") [options]

Prepare a Codex development environment for Rodin.

Options:
  --install-deps           Install dependencies (Ubuntu/macOS only)
  --skip-configure         Do not run CMake configure
  --build-dir <path>       Override build directory (default: ${BUILD_DIR})
  --build-type <type>      CMake build type (default: ${BUILD_TYPE})
  -h, --help               Show this help
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --install-deps)
      INSTALL_DEPS="1"
      shift
      ;;
    --skip-configure)
      RUN_CONFIGURE="0"
      shift
      ;;
    --build-dir)
      BUILD_DIR="$2"
      shift 2
      ;;
    --build-type)
      BUILD_TYPE="$2"
      shift 2
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

install_ubuntu_deps() {
  sudo apt-get update
  sudo apt-get install -y \
    libboost1.74-all-dev \
    libsuitesparse-dev \
    libeigen3-dev \
    libscotch-dev \
    libmetis-dev \
    libhdf5-dev \
    libomp-dev \
    petsc-dev \
    mpich \
    lcov
}

install_macos_deps() {
  brew update
  brew install \
    boost boost-mpi \
    suitesparse \
    eigen \
    scotch \
    metis \
    hdf5-mpi \
    libomp \
    petsc \
    open-mpi \
    lcov
}

echo "[codex-setup] Repository root: ${REPO_ROOT}"
cd "${REPO_ROOT}"

echo "[codex-setup] Syncing submodules"
git submodule update --init --recursive

if [[ "${INSTALL_DEPS}" == "1" ]]; then
  echo "[codex-setup] Installing dependencies"
  case "$(uname -s)" in
    Linux)
      install_ubuntu_deps
      ;;
    Darwin)
      install_macos_deps
      ;;
    *)
      echo "Unsupported OS for auto-install: $(uname -s)" >&2
      exit 1
      ;;
  esac
fi

if [[ "${RUN_CONFIGURE}" == "1" ]]; then
  echo "[codex-setup] Configuring CMake in ${BUILD_DIR}"
  cmake -S "${REPO_ROOT}" -B "${BUILD_DIR}" \
    -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
    -DRODIN_BUILD_DOC=OFF \
    -DRODIN_USE_MCSS=OFF \
    -DRODIN_BUILD_SRC=ON \
    -DRODIN_BUILD_UNIT_TESTS=ON \
    -DRODIN_BUILD_MANUFACTURED_TESTS=ON \
    -DRODIN_BUILD_EXAMPLES=OFF \
    -DRODIN_BUILD_BENCHMARKS=OFF
  echo "[codex-setup] Done. Build with: cmake --build ${BUILD_DIR} -j$(nproc 2>/dev/null || echo 4)"
else
  echo "[codex-setup] Skipped CMake configure"
fi
