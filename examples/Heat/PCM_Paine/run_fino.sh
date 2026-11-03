#!/usr/bin/env bash
# Fine-mesh sweep of the Al / RT42 / Al roof sandwich on another machine:
# builds NaturalConvectionPCM (branch model/HeatTransfer), runs 5-20 mm over
# 2 periodic days on the full 2 m slope and post-processes (heat-flux summary,
# plots, temperature contours and profiles).
#
# First time on the machine:
#   git clone -b model/HeatTransfer https://github.com/cbritopacheco/rodin.git ~/rodin-heat
#   cd ~/rodin-heat && git submodule update --init third-party/termcolor third-party/mmg
# Then (in a directory for the results):
#   mkdir -p ~/runs_fino && cd ~/runs_fino
#   NP=8 NJOBS=4 bash ~/rodin-heat/examples/Heat/PCM_Paine/run_fino.sh
#
# PETSc / MPI not found by CMake: pass them, e.g.
#   CMAKE_ARGS="-DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx \
#     -DPETSc_INCLUDE_DIR='$PETSC_DIR/include;$PETSC_DIR/$PETSC_ARCH/include' \
#     -DPETSc_LIBRARY=$PETSC_DIR/$PETSC_ARCH/lib/libpetsc.so" bash run_fino.sh
#
# Mesh: h_x = 5 mm (400 cells over 2 m); PCM: 0.25 mm at its top face growing
# by 1.04 (20 cells at 5 mm, 37 at 20 mm); 4 cells per 10 mm plate. Mushy
# range 4 K (39-43 C; RT42 melts over 38-43 C) for every thickness: with 2 K,
# 20 mm loses Newton when the melt flows (~4e-4 m/s); dt 30 s. Any variable
# of run_pcm_sweep.sh can be
# overridden (SITE=paine, THICKNESSES, NDAYS, ...). The sweep skips finished
# cases, so the script can be relaunched after an interruption.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
BUILD="${BUILD:-$REPO/build}"
CMAKE_ARGS="${CMAKE_ARGS:-}"
BUILD_JOBS="${BUILD_JOBS:-4}"

echo "[rodin] $(git -C "$REPO" rev-parse --abbrev-ref HEAD) $(git -C "$REPO" rev-parse --short HEAD)"
if [ ! -f "$BUILD/CMakeCache.txt" ]; then
  # shellcheck disable=SC2086
  eval cmake -S "\"$REPO\"" -B "\"$BUILD\"" -DCMAKE_BUILD_TYPE=Release \
    -DRODIN_USE_MPI=ON -DRODIN_USE_PETSC=ON $CMAKE_ARGS
fi
cmake --build "$BUILD" -j "$BUILD_JOBS" --target NaturalConvectionPCM

export BUILD
export OUT="${OUT:-$PWD/pcm_fino}"
export THICKNESSES="${THICKNESSES:-5 10 15 20}"
export SITE="${SITE:-calama}"
export LENGTH="${LENGTH:-2.0}"
export HX="${HX:-0.005}"
export HYTOP="${HYTOP:-0.00025}"
export GRADE="${GRADE:-1.04}"
export NYMIN="${NYMIN:-20}"
export AL="${AL:-0.01}"
export NAL="${NAL:-4}"
export NDAYS="${NDAYS:-2}"
export MUSHY="${MUSHY:-4}"
export EXTRA="${EXTRA:--pcm_round 1}"
export RETRIES="${RETRIES:-0}"
export NP="${NP:-8}"
export NJOBS="${NJOBS:-4}"

bash "$HERE/run_pcm_sweep.sh"
bash "$HERE/postprocess_pcm.sh" "$OUT" --hours 6 9 12 15 18 21 --all-contours
