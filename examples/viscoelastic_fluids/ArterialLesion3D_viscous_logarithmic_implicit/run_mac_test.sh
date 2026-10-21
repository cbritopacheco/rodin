#!/usr/bin/env bash
# =============================================================================
#  run_mac_test.sh -- the hardest 3D case on a laptop (macOS or Linux):
#  S75, fluid R (the gsPTT blood model), Re 600, Wo 5.5, sinusoidal inflow
#  A = 0.5, BDF2, psi_rb, on the layer-refined coarse mesh S75_pipe_layersc
#  (121k tetrahedra, 37% of inter, the throat wall and the jet shear layer
#  at the h of inter).
#
#     cd ~/rodin && bash examples/viscoelastic_fluids/ArterialLesion3D_viscous_logarithmic_implicit/run_mac_test.sh
#
#  It builds the 3D driver (SKIP_BUILD=1 to skip), checks that the build has
#  Scotch, then runs in the background and prints how to follow it. At the
#  end (or when it stops on an error) locate_blowup.py writes where |sigma|,
#  |u| and psi were extreme in the last outputs to blowup.txt.
#
#  Environment (defaults):
#     NP       performance cores (sysctl hw.perflevel0.physicalcpu on macOS)
#     STEPS=400 (steps per cycle: dt 2.5 times that of the campaign, a solver
#              test, not a resolved run)   CYCLES=2 (the ramp cycle and one more)
#     OUTPUT_EVERY=20
#     SOLVER   lulag (MUMPS LU, lagged) with >= 16 GB of RAM, else gmres
#     VARIANT  base (psi_rb) | nodiv (psi_rb without the div sigma subscale)
#              | psi (the former formulation)
#     BUILD=build   OUT=$HOME/Results/mac_test   MPIRUN=mpirun
#  The earlier 3D failures (pc9): Re 600 S75, psi, at t = 0.108 s; Re 300
#  with the CCA waveform at t = 0.28 s (S75) and 0.59 s (A200).
# =============================================================================
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
BUILD="${BUILD:-build}"
OUT="${OUT:-$HOME/Results/mac_test}"
MPIRUN="${MPIRUN:-mpirun}"
STEPS="${STEPS:-400}"
OUTPUT_EVERY="${OUTPUT_EVERY:-20}"
CYCLES="${CYCLES:-2}"
VARIANT="${VARIANT:-base}"
MESHDIR="$REPO/resources/examples/viscoelastic_fluids"
EXE="$REPO/$BUILD/examples/viscoelastic_fluids/ArterialLesion3D_viscous_logarithmic_implicit/ArterialLesion3D_viscous_logarithmic_implicit"
PY="${PYTHON:-$(command -v python3)}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

# ---- machine ------------------------------------------------------------------
if [[ "$(uname)" == Darwin ]]; then
  ncores=$(sysctl -n hw.perflevel0.physicalcpu 2> /dev/null || sysctl -n hw.physicalcpu)
  memgb=$(( $(sysctl -n hw.memsize) / 1073741824 ))
else
  ncores=$(nproc)
  memgb=$(( $(awk '/MemTotal/ {print $2}' /proc/meminfo) / 1048576 ))
fi
NP="${NP:-$ncores}"
if [[ -z "${SOLVER:-}" ]]; then
  if (( memgb >= 16 )); then SOLVER=lulag; else SOLVER=gmres; fi
fi
case "$VARIANT" in
  base)  FORM=psi_rb; EXTRA="" ;;
  nodiv) FORM=psi_rb; EXTRA="-al_stress_div_scale 0" ;;
  psi)   FORM=psi;    EXTRA="" ;;
  *) echo "VARIANT must be base, nodiv or psi"; exit 1 ;;
esac
echo "== $(uname -s), $ncores performance cores, ${memgb} GB: NP=$NP, solver $SOLVER, variant $VARIANT ($FORM $EXTRA)"

# ---- build and checks --------------------------------------------------------------
[[ -f "$REPO/$BUILD/CMakeCache.txt" ]] || { echo "no build in $REPO/$BUILD (set BUILD=...)"; exit 1; }
if ! grep -qiE '^RODIN_USE_SCOTCH:BOOL=(ON|YES|TRUE|Y|1)$' "$REPO/$BUILD/CMakeCache.txt"; then
  echo "warning: $BUILD is configured without Scotch; the partition will be poor"
  echo "         (cmake -S $REPO -B $REPO/$BUILD -DRODIN_USE_SCOTCH=ON)"
fi
if [[ "${SKIP_BUILD:-0}" != 1 ]]; then
  echo "== building the 3D driver (~6 GB of RAM per compiler process; -j 2) ..."
  cmake --build "$REPO/$BUILD" -j 2 --target ArterialLesion3D_viscous_logarithmic_implicit \
    > "$REPO/$BUILD/mac_test_build.log" 2>&1 \
    || { echo "BUILD FAILED, see $REPO/$BUILD/mac_test_build.log"; tail -20 "$REPO/$BUILD/mac_test_build.log"; exit 1; }
fi
[[ -x "$EXE" ]] || { echo "missing $EXE"; exit 1; }
grep -aq -- "psi_rb" "$EXE" || { echo "$EXE has no psi_rb: rebuild from the current sources"; exit 1; }
MESH="$MESHDIR/S75_pipe_layersc.mesh"
if [[ ! -f "$MESH" ]] || ! head -1 "$MESH" | grep -q MeshVersionFormatted; then
  echo "== making $MESH (python gmsh) ..."
  "$PY" "$HERE/make_lesion_mesh_3d.py" S75 0.7 "$MESH" -setnumber layers 1 || exit 1
fi

# ---- run ----------------------------------------------------------------------------
mkdir -p "$OUT"
DIR="$OUT/Re600_Wo5.5/S75/R"
LOG="$OUT/mac_test_$VARIANT.log"
(
  env FLUIDS=R NP="$NP" TIME_ORDER=2 MESH_LEVEL=layersc OUTPUT_EVERY="$OUTPUT_EVERY" STEPS="$STEPS" CYCLES_R="$CYCLES" \
    SOLVER="$SOLVER" FORMULATION="$FORM" PETSC_EXTRA="$EXTRA -al_trace_newton 1" EXE3D="$EXE" \
    MPIRUN="$MPIRUN" OUT3D="$OUT" bash "$HERE/run_3d.sh" 600 S75
  if "$PY" -c "import numpy, h5py" 2> /dev/null && [[ -f "$DIR/case.xdmf" ]]; then
    "$PY" "$HERE/locate_blowup.py" "$DIR" --last 6 > "$DIR/blowup.txt" 2>&1
  fi
  echo "== mac test finished: $(date)"
) > "$LOG" 2>&1 < /dev/null &
echo "== started (pid $!): log $LOG, case $DIR"
echo "   follow:   tail -f $DIR/run.log | grep -E 'Step|timing|reject|Fatal|partition'"
echo "   CSV:      $DIR/case.csv  (maxSigma, newtonIts, dt, retries per step)"
echo "   stop:     pkill -f ArterialLesion3D_viscous_logarithmic_implicit"
echo "   at the end: $DIR/blowup.txt (needs python3 with numpy and h5py)"
