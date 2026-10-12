#!/usr/bin/env bash
# =============================================================================
#  launch_pc9.sh -- builds the drivers and launches, in the background, the
#  pending cases of the JFM study on one machine (pc9, 200 cores).
#
#     cd ~/rodin && git pull && bash examples/viscoelastic_fluids/launch_pc9.sh
#
#  What it launches (each job with nohup; it survives logging out):
#    axi  aneurysms A125, A150 (A200 is done), with run_axi_re.sh:
#         Re 100, 300: fluid R (12 ranks) and fluids D N CY (8 ranks), per geometry
#         Re 600:      fluids D N CY (8 ranks), per geometry           -> 96 cores
#    3d   with run_3d.sh, 30 ranks each:
#         S75 Re 600 N, then S75 Re 300 N | S75 Re 600 D | A200 Re 600 N then D
#                                                                       -> 90 cores
#  Usage: bash launch_pc9.sh [all | axi | 3d | 3d-restart | 3d-one | 3d-a200]   (default all)
#         3d-restart stops the running 3D jobs (only those) and relaunches them
#         from scratch with the current code; finished 3D cases (DONE) are kept.
#         3d-one stops the running 3D jobs and launches only the key case,
#         S75 Re 600 viscoelastic (D, Wi 2), with NP3D=60 ranks, CYCLES=3 and
#         GMRES (falls back to MUMPS by itself if GMRES fails).
#         3d-a200 adds (stops nothing) A200 Re 600 viscoelastic, same settings
#         as 3d-one: the aneurysm counterpart of the S75 case.
#  Environment: BUILD=build-gcc12  OUT=$HOME/Results/axi_jfm  CORES=200
#               FORCE=1 launches even if the machine looks busy
#  Logs: $OUT/logs/*.log ; per case: <case folder>/run.log, case.csv
#  Stop everything: pkill -f run_axi_; pkill -f run_3d.sh; sleep 2; pkill -f ArterialLesion
# =============================================================================
set -u
WHAT="${1:-all}"
case "$WHAT" in all|axi|3d|3d-restart|3d-one|3d-a200) ;; *) sed -n '2,32p' "$0"; exit 1 ;; esac
RESTART=0; ONE=0
if [[ "$WHAT" == 3d-restart ]]; then WHAT=3d; RESTART=1; fi
if [[ "$WHAT" == 3d-one ]]; then WHAT=3d; RESTART=1; ONE=1; fi
GEO1=S75
if [[ "$WHAT" == 3d-a200 ]]; then WHAT=3d; ONE=1; GEO1=A200; fi
NP3D="${NP3D:-60}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
cd "$REPO"

# No Anaconda on the PATH (its mpirun and Boost are not the ones we built with).
PATH="$(printf '%s' "$PATH" | tr ':' '\n' | grep -viE 'anaconda|miniconda|conda/' | paste -sd: -)"
export PATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export OMPI_MCA_hwloc_base_binding_policy=none

BUILD="${BUILD:-build-gcc12}"
export OUT="${OUT:-$HOME/Results/axi_jfm}"
CORES="${CORES:-200}"
export EXE="$REPO/$BUILD/examples/viscoelastic_fluids/ArterialLesionAxi_viscous_logarithmic_implicit/ArterialLesionAxi_viscous_logarithmic_implicit"
export EXE3D="$REPO/$BUILD/examples/viscoelastic_fluids/ArterialLesion3D_viscous_logarithmic_implicit/ArterialLesion3D_viscous_logarithmic_implicit"
S="$HERE/ArterialLesionAxi_viscous_logarithmic_implicit/run_axi_re.sh"
S3="$HERE/ArterialLesion3D_viscous_logarithmic_implicit/run_3d.sh"
MESH="$REPO/resources/examples/viscoelastic_fluids"
LOG="$OUT/logs"
mkdir -p "$LOG"

echo "== repo $REPO ($(git --no-optional-locks log --oneline -1 2>/dev/null))"
echo "== mpirun: $(command -v mpirun)"

# ---- build ------------------------------------------------------------------
if [[ ! -d "$REPO/$BUILD" ]]; then
  echo "Build directory $REPO/$BUILD not found (set BUILD=...)"; exit 1
fi
echo "== building ..."
if ! cmake --build "$REPO/$BUILD" -j 32 --target ArterialLesionAxi_viscous_logarithmic_implicit \
      ArterialLesion3D_viscous_logarithmic_implicit > "$LOG/build.log" 2>&1; then
  echo "BUILD FAILED, see $LOG/build.log:"; tail -20 "$LOG/build.log"; exit 1
fi
for x in "$EXE" "$EXE3D"; do
  [[ -x "$x" ]] || { echo "missing executable $x"; exit 1; }
done

# ---- meshes -------------------------------------------------------------------
need=()
[[ "$WHAT" != 3d ]] && need+=(A125_axi_medium A150_axi_medium)
[[ "$WHAT" != axi ]] && need+=(S75_pipe_medium A200_pipe_medium)
for m in "${need[@]}"; do
  f="$MESH/$m.mesh"
  if [[ ! -f "$f" ]] || ! head -1 "$f" | grep -q MeshVersionFormatted; then
    echo "mesh $f missing or not MEDIT (Git LFS pointer?)"; exit 1
  fi
done
if [[ "$WHAT" != axi ]]; then
  if python3 -c "import numpy" 2> /dev/null; then
    echo "== checking the 3D meshes ..."
    python3 "$HERE/ArterialLesion3D_viscous_logarithmic_implicit/check_3d_mesh.py" \
      S75 "$MESH/S75_pipe_medium.mesh" A200 "$MESH/A200_pipe_medium.mesh" | tee "$LOG/check_3d_mesh.log" \
      || { echo "3D mesh check failed"; exit 1; }
  else
    echo "(numpy not available: 3D mesh check skipped)"
  fi
fi

# ---- stop the running 3D jobs (3d-restart) ------------------------------------------
if [[ $RESTART == 1 ]]; then
  echo "== stopping the running 3D jobs ..."
  pkill -f "run_3d\.sh" ; sleep 2
  pkill -f "ArterialLesion3D_viscous_logarithmic_implicit -al_" ; sleep 10
  pkill -9 -f "ArterialLesion3D_viscous_logarithmic_implicit -al_" 2> /dev/null; sleep 2
  for d in "$OUT"/3D/Re*/*/*; do
    [[ -d "$d" && ! -f "$d/DONE" ]] && rm -f "$d/FAILED" && sed -i '/,FAILED,/d' "$(dirname "$(dirname "$d")")/summary.csv" 2> /dev/null
  done
fi

# ---- free cores -----------------------------------------------------------------
want=0
[[ "$WHAT" != 3d ]] && want=$(( want + 96 ))
if [[ "$WHAT" != axi ]]; then
  if [[ $ONE == 1 ]]; then want=$(( want + NP3D )); else want=$(( want + 90 )); fi
fi
busy=$(pgrep -f "ArterialLesion(Axi|3D)_viscous_logarithmic_implicit -al_" | wc -l)
echo "== solver processes already running: $busy ; this launch adds $want ranks ; cores: $CORES"
if (( busy + want > CORES )) && [[ "${FORCE:-0}" != 1 ]]; then
  echo "Not enough free cores ($busy busy + $want > $CORES). Wait, or rerun with FORCE=1, or launch only 'axi' or '3d'."
  exit 1
fi

launch() {   # log-name, then the command (env assignments allowed via env)
  local name="$1"; shift
  nohup "$@" > "$LOG/$name.log" 2>&1 < /dev/null &
  echo "   started $name (pid $!) -> $LOG/$name.log"
}

# ---- axisymmetric aneurysms ------------------------------------------------------
if [[ "$WHAT" != 3d ]]; then
  echo "== axisymmetric aneurysms A125 A150"
  for RE in 100 300; do
    for G in A125 A150; do
      launch "Re${RE}_${G}_R"    env FLUIDS=R        NP=12 bash "$S" "$RE" "$G"
      launch "Re${RE}_${G}_DNCY" env FLUIDS="D N CY" NP=8  bash "$S" "$RE" "$G"
    done
  done
  for G in A125 A150; do
    launch "Re600_${G}_DNCY" env FLUIDS="D N CY" NP=8 bash "$S" 600 "$G"
  done
fi

# ---- 3D -------------------------------------------------------------------------
if [[ "$WHAT" != axi && $ONE == 1 ]]; then
  echo "== 3D, $GEO1 Re 600 viscoelastic: $NP3D ranks, ${CYCLES:-3} cycles, ${SOLVER:-gmres}"
  launch "3D_${GEO1}_Re600_D" env FLUIDS=D NP="$NP3D" CYCLES="${CYCLES:-3}" SOLVER="${SOLVER:-gmres}" bash "$S3" 600 "$GEO1"
elif [[ "$WHAT" != axi ]]; then
  echo "== 3D (30 ranks each)"
  launch 3D_S75_N           bash -c "FLUIDS=N NP=30 bash '$S3' 600 S75; FLUIDS=N NP=30 bash '$S3' 300 S75"
  launch 3D_S75_Re600_D     env FLUIDS=D NP=30 bash "$S3" 600 S75
  launch 3D_A200_Re600_ND   env FLUIDS="N D" NP=30 bash "$S3" 600 A200
fi

sleep 5
echo
echo "== running solver processes now: $(pgrep -f "ArterialLesion(Axi|3D)_viscous_logarithmic_implicit -al_" | wc -l)"
if [[ $ONE == 1 ]]; then
  echo "Follow:   tail -f $LOG/3D_${GEO1}_Re600_D.log"
  echo "3D ETA:   grep -E 'timing|ETA' $OUT/3D/Re600_Wo5.5/$GEO1/Wi2/run.log | tail -2     (after ~15 min)"
else
  echo "Follow:   tail -f $LOG/3D_S75_N.log"
  echo "3D ETA:   grep ETA $OUT/3D/Re600_Wo5.5/S75/N/run.log | tail -1     (after ~15 min)"
fi
echo "Summary:  cat $OUT/Re*/summary.csv $OUT/3D/Re*/summary.csv"
