#!/usr/bin/env bash
# Configures (if needed) and builds the two executables of the campaign.
#   BUILD=~/rodin/build ./00_build.sh
# PETSc and MPI are required; set PETSC_DIR (and PETSC_ARCH) if CMake does
# not find PETSc by itself.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

if [[ ! -f "$BUILD/CMakeCache.txt" ]]; then
  cmake -S "$REPO" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release \
    -DRODIN_USE_MPI=ON -DRODIN_USE_PETSC=ON -DRODIN_BUILD_EXAMPLES=ON \
    -DRODIN_INSTALL_RESOURCES=OFF
fi
cmake --build "$BUILD" -j "${BUILD_JOBS:-$CORES}" --target RmuAxi CoronaryArtery
ls -la "$RMU_EXE" "$COR_EXE"
