#!/usr/bin/env bash
# =============================================================================
#  launch_bdf2_pc9.sh -- the BDF2 campaign of the JFM study on pc9, 150 cores:
#
#    3D   (6 x 20 = 120 cores), m = 1 seed in cycle 1, intermediate meshes,
#         always CY against R (the viscoelastic blood model and its memory-free
#         twin, same steady viscosity: the difference is the memory alone):
#         - the carotid, Re 300 Wo 4, with its physiological CCA waveform
#           (peak/mean 3.9): S75 and A200;
#         - Re 600 Wo 5.5, sinusoid A = 0.5: S75.
#    axi  (6 x 5 = 30 cores, a queue): every geometry (S75 S50 H A125 A150
#         A200) x fluid (N CY D R) x (Re, Wo) = (100, 2.5), (300, 4), (600, 5.5),
#         sinusoidal inflow; plus the CCA waveform at Re 300 Wo 4 for S75, A200.
#
#  R has lambda = 4.83 s (LAMBDA_R) at every Re: the same blood in the same
#  geometry, El = lambda nu/D^2 fixed, so De = lambda/T = (2/pi) El Wo^2 grows
#  as Wo^2 (about 3.9, 10, 19) and Wi = El Re as Re. Its solvent
#  eta_CY - eta_p/f is positive for lambda >= 2.4 s.
#
#     cd ~/rodin && bash examples/viscoelastic_fluids/launch_bdf2_pc9.sh [all|3d|axi]
#
#  Time step (BDF2): U_jet dt/d_throat at the mean flow of S75:
#     axi  600 / 800 / 1200 steps per cycle at Re 100 / 300 / 600 (0.33, 0.30,
#          0.21), 1600 for the CCA waveform (0.15 at the mean, 0.6 at the peak);
#     3D   1600 (CCA waveform, Re 300) and 1000 (sinusoid, Re 600).
#  Newton tolerance 5e-4 (quadratic: ~2 iterations per step with the BDF2
#  extrapolated start). Cycles: 4 (N, CY, D), 8 (axi R), 5 (3D R: compare
#  the growth of uyMean, uzMean with CY).
#  Rough cost (pc9 speed unknown; check the ETA lines after ~20 min): 3D
#  ~60-100 s per step at 20 ranks -> 4-8 days; axi queue ~4-6 days.
#
#  Environment: CORES=200 (machine)  NP3D=20  NPAXI=5  MAXAXI=6
#               BUILD=build-gcc12  OUT=$HOME/Results/axi_jfm_bdf2  FORCE=1
#               LAMBDA_R=4.83
#  Logs: $OUT/logs/   3D: $OUT/3D/Re<Re>_Wo<Wo>/<geo>/<fluid>/
#        axi: $OUT/Re<Re>_Wo<Wo>[_waveform_cca]/<geo>/<fluid>/
#  Cases that finished (DONE) are skipped: the script can be relaunched.
# =============================================================================
set -u
WHAT="${1:-all}"
case "$WHAT" in all|3d|axi) ;; *) sed -n '2,30p' "$0"; exit 1 ;; esac

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
cd "$REPO"

# No Anaconda on the PATH (its mpirun and Boost are not the ones we built with).
PATH="$(printf '%s' "$PATH" | tr ':' '\n' | grep -viE 'anaconda|miniconda|conda/' | paste -sd: -)"
export PATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export OMPI_MCA_hwloc_base_binding_policy=none

BUILD="${BUILD:-build-gcc12}"
export OUT="${OUT:-$HOME/Results/axi_jfm_bdf2}"
export LAMBDA_R="${LAMBDA_R:-4.83}"
CORES="${CORES:-200}"
NP3D="${NP3D:-20}"
NPAXI="${NPAXI:-5}"
MAXAXI="${MAXAXI:-6}"
export EXE="$REPO/$BUILD/examples/viscoelastic_fluids/ArterialLesionAxi_viscous_logarithmic_implicit/ArterialLesionAxi_viscous_logarithmic_implicit"
export EXE3D="$REPO/$BUILD/examples/viscoelastic_fluids/ArterialLesion3D_viscous_logarithmic_implicit/ArterialLesion3D_viscous_logarithmic_implicit"
S3="$HERE/ArterialLesion3D_viscous_logarithmic_implicit/run_3d.sh"
Q="$HERE/ArterialLesionAxi_viscous_logarithmic_implicit/run_axi_queue.sh"
MESH="$REPO/resources/examples/viscoelastic_fluids"
WF="${WAVEFORM:-$MESH/waveform_cca.dat}"
LOG="$OUT/logs"
mkdir -p "$LOG"

echo "== repo $REPO ($(git --no-optional-locks log --oneline -1 2>/dev/null))"
echo "== mpirun: $(command -v mpirun)"

# ---- build ------------------------------------------------------------------
[[ -d "$REPO/$BUILD" ]] || { echo "Build directory $REPO/$BUILD not found (set BUILD=...)"; exit 1; }
echo "== building ..."
if ! cmake --build "$REPO/$BUILD" -j 32 --target ArterialLesionAxi_viscous_logarithmic_implicit \
      ArterialLesion3D_viscous_logarithmic_implicit > "$LOG/build.log" 2>&1; then
  echo "BUILD FAILED, see $LOG/build.log:"; tail -20 "$LOG/build.log"; exit 1
fi
for x in "$EXE" "$EXE3D"; do
  [[ -x "$x" ]] || { echo "missing executable $x"; exit 1; }
  grep -aq -- "-al_time_order" "$x" || { echo "$x has no BDF2 (-al_time_order): wrong sources?"; exit 1; }
done
grep -aq -- "-al_fluid" "$EXE3D" || { echo "$EXE3D has no -al_fluid (3D CY/R): wrong sources?"; exit 1; }

# ---- meshes and waveform ----------------------------------------------------------
need=()
[[ "$WHAT" != axi ]] && need+=(S75_pipe_inter A200_pipe_inter)
[[ "$WHAT" != 3d ]] && need+=(S75_axi_medium S50_axi_medium H_axi_medium A125_axi_medium A150_axi_medium A200_axi_medium)
for m in "${need[@]}"; do
  f="$MESH/$m.mesh"
  if [[ ! -f "$f" ]] || ! head -1 "$f" | grep -q MeshVersionFormatted; then
    echo "mesh $f missing or not MEDIT (Git LFS pointer? run: git lfs pull)"; exit 1
  fi
done
if [[ ! -f "$WF" ]]; then
  echo "waveform $WF not found"; exit 1
fi
if [[ "$WHAT" != axi ]] && python3 -c "import numpy" 2> /dev/null; then
  echo "== checking the 3D meshes ..."
  python3 "$HERE/ArterialLesion3D_viscous_logarithmic_implicit/check_3d_mesh.py" \
    S75 "$MESH/S75_pipe_inter.mesh" A200 "$MESH/A200_pipe_inter.mesh" | tee "$LOG/check_3d_mesh.log"
  grep -q PROBLEMS "$LOG/check_3d_mesh.log" && { echo "3D mesh check failed"; exit 1; }
fi

# ---- free cores ---------------------------------------------------------------------
want=0
[[ "$WHAT" != axi ]] && want=$(( want + 6 * NP3D ))
[[ "$WHAT" != 3d ]] && want=$(( want + MAXAXI * NPAXI ))
busy=$(ps -eo comm | grep -c '^ArterialLesion')
echo "== solver processes already running: $busy ; this launch adds $want ranks ; cores: $CORES"
if (( busy + want > CORES )) && [[ "${FORCE:-0}" != 1 ]]; then
  echo "Not enough free cores ($busy busy + $want > $CORES). Stop the old runs, or FORCE=1."
  exit 1
fi

launch() {   # log-name, then the command
  local name="$1"; shift
  nohup "$@" > "$LOG/$name.log" 2>&1 < /dev/null &
  echo "   started $name (pid $!) -> $LOG/$name.log"
}

# ---- 3D -----------------------------------------------------------------------------
if [[ "$WHAT" != axi ]]; then
  echo "== 3D: CY and R (lambda_R $LAMBDA_R s) at Re 300 (CCA waveform, S75 A200) and Re 600 (S75); $NP3D ranks each, BDF2, mesh inter"
  for G in S75 A200; do
    for F in CY R; do
      launch "3D_${G}_Re300_wave_${F}" env FLUIDS="$F" NP="$NP3D" TIME_ORDER=2 MESH_LEVEL=inter \
        OUTPUT_EVERY=160 WAVEFORM="$WF" HARMONICS=10 STEPS=1600 CYCLES=4 CYCLES_R=5 \
        OUT3D="$OUT/3D" bash "$S3" 300 "$G"
    done
  done
  for F in CY R; do
    launch "3D_S75_Re600_${F}" env FLUIDS="$F" NP="$NP3D" TIME_ORDER=2 MESH_LEVEL=inter \
      OUTPUT_EVERY=100 STEPS=1000 CYCLES=4 CYCLES_R=5 OUT3D="$OUT/3D" bash "$S3" 600 S75
  done
fi

# ---- axisymmetric queue ---------------------------------------------------------------
if [[ "$WHAT" != 3d ]]; then
  JOBS="$OUT/axi_jobs.txt"
  {
    echo "# Re geometry fluids steps [waveform]   (longest first)"
    for RE in 600 300 100; do
      case "$RE" in 600) ST=1200 ;; 300) ST=800 ;; 100) ST=600 ;; esac
      for G in S75 A200 S50 A150 H A125; do
        echo "$RE $G R $ST"
        echo "$RE $G N,CY,D $ST"
      done
    done
    for G in S75 A200; do
      echo "300 $G R 1600 $WF"
      echo "300 $G N,CY,D 1600 $WF"
    done
  } > "$JOBS"
  echo "== axisymmetric: $(grep -vc '^#' "$JOBS") jobs, $MAXAXI at a time, $NPAXI ranks each -> $JOBS"
  launch axi_queue env MAXJOBS="$MAXAXI" NP="$NPAXI" TIME_ORDER=2 NEWTON_TOL=5e-4 \
    HARMONICS=10 CYCLES=4 CYCLES_R=8 RETRY=1 OUT="$OUT" LOGDIR="$LOG" bash "$Q" "$JOBS"
fi

sleep 5
echo
echo "== running solver processes now: $(ps -eo comm | grep -c '^ArterialLesion')"
echo "Progress:  tail -f $LOG/axi_queue.log"
echo "3D ETA:    grep ETA $OUT/3D/Re300_Wo4_waveform_cca/S75/R/run.log | tail -1     (after ~20 min)"
echo "3D m=1:    the columns uyMean, uzMean of case.csv (seeded in cycle 1): growth = unstable"
echo "Summary:   cat $OUT/Re*/summary.csv $OUT/3D/Re*/summary.csv"
