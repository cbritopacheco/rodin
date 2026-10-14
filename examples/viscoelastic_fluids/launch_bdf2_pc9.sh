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
#     cd ~/rodin && bash examples/viscoelastic_fluids/launch_bdf2_pc9.sh [mode]
#
#  Modes:
#     all      (default) 3D and the axisymmetric queue
#     3d, axi  only that part
#     status   what is running: solver ranks, each running case with its last
#              step, ETA and Newton iterations, finished and failed cases
#     verif    the verification against exact solutions (needs no other
#              mode): steady sPTT pipe (3D, rf 0.7-2.0, 76 cores), then the
#              pulsatile Oldroyd-B pipe (Womersley) with the 3D driver (8 x 8
#              cores) and the axisymmetric one (8 x 2 cores) at the same time.
#              Results in $OUT/verification/{sptt,womersley_3d,womersley_axi}
#     verif3d  stops the running 3D cases (only those: the axisymmetric queue
#              goes on), runs verif, then launches the 3D cases
#     stop3d   stops the running 3D cases (finished ones are kept)
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
case "$WHAT" in all|3d|axi|status|verif|verif3d|stop3d|_verif) ;; *) sed -n '2,50p' "$0"; exit 1 ;; esac

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
cd "$REPO"

# Python (numpy, h5py) for the verification, taken before conda leaves the PATH.
export PYTHON="${PYTHON:-$(command -v python3)}"
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

csvcol() {   # file column-name: last value of that column
  awk -F, -v c="$2" 'NR==1{for(i=1;i<=NF;i++) if($i==c) k=i; next} END{if(k) print $k; else print "-"}' "$1" 2> /dev/null
}

if [[ "$WHAT" == status ]]; then
  echo "== $(date '+%a %d %H:%M')  solver ranks running: 3D $(ps -eo args | grep -c '^[^ ]*ArterialLesion3D_viscous_logarithmic_implicit -al_')" \
       " axi $(ps -eo args | grep -c '^[^ ]*ArterialLesionAxi_viscous_logarithmic_implicit -al_')"
  echo "== running cases (run.log written in the last 15 min):"
  find "$OUT" -name run.log -mmin -15 2> /dev/null | sort | while read -r f; do
    d=$(dirname "$f")
    total=$(grep -a -m1 -- '---- Step' "$f" | sed 's/.*Step [0-9]*\/\([0-9]*\).*/\1/')
    done_=$(( $(wc -l < "$d/case.csv" 2> /dev/null || echo 1) - 1 ))
    eta=$(grep -a 'ETA' "$f" | tail -1 | sed 's/.*ETA/ETA/')
    printf '   %-45s step %s/%s  cycle %s  newtonIts %s  %s\n' "${d#$OUT/}" "$done_" "$total" \
      "$(csvcol "$d/case.csv" cycle)" "$(csvcol "$d/case.csv" newtonIts)" "$eta"
  done
  echo "== finished: $(find "$OUT" -name DONE 2> /dev/null | wc -l)   failed:"
  find "$OUT" -name FAILED 2> /dev/null | sed "s|$OUT/|   |; s|/FAILED||"
  [[ -f "$LOG/axi_queue.log" ]] && { echo "== axi queue:"; tail -4 "$LOG/axi_queue.log"; }
  [[ -f "$LOG/verif.log" ]] && { echo "== verification:"; tail -4 "$LOG/verif.log"; }
  exit 0
fi

stop3d() {
  echo "== stopping the running 3D cases ..."
  pkill -f "run_3d\.sh" ; sleep 2
  pkill -f "ArterialLesion3D_viscous_logarithmic_implicit -al_" ; sleep 10
  pkill -9 -f "ArterialLesion3D_viscous_logarithmic_implicit -al_" 2> /dev/null; sleep 2
  for d in "$OUT"/3D/Re*/*/*; do
    [[ -d "$d" && ! -f "$d/DONE" ]] && rm -f "$d/FAILED"
  done
  for f in "$OUT"/3D/Re*/summary.csv; do [[ -f "$f" ]] && sed -i '/,FAILED,/d' "$f"; done
}
if [[ "$WHAT" == stop3d ]]; then stop3d; exit 0; fi

echo "== repo $REPO ($(git --no-optional-locks log --oneline -1 2>/dev/null))"
echo "== mpirun: $(command -v mpirun)"

# ---- build ------------------------------------------------------------------
[[ -d "$REPO/$BUILD" ]] || { echo "Build directory $REPO/$BUILD not found (set BUILD=...)"; exit 1; }
if [[ "${SKIP_BUILD:-0}" != 1 ]]; then
echo "== building ..."
if ! cmake --build "$REPO/$BUILD" -j 32 --target ArterialLesionAxi_viscous_logarithmic_implicit \
      ArterialLesion3D_viscous_logarithmic_implicit > "$LOG/build.log" 2>&1; then
  echo "BUILD FAILED, see $LOG/build.log:"; tail -20 "$LOG/build.log"; exit 1
fi
fi
for x in "$EXE" "$EXE3D"; do
  [[ -x "$x" ]] || { echo "missing executable $x"; exit 1; }
  grep -aq -- "-al_time_order" "$x" || { echo "$x has no BDF2 (-al_time_order): wrong sources?"; exit 1; }
done
grep -aq -- "-al_fluid" "$EXE3D" || { echo "$EXE3D has no -al_fluid (3D CY/R): wrong sources?"; exit 1; }

# ---- verification (the hidden mode _verif runs it in the foreground) ----------
V="$HERE/ArterialLesion3D_viscous_logarithmic_implicit/verification"
if [[ "$WHAT" == _verif ]]; then
  "$PYTHON" -c "import numpy, h5py" 2> /dev/null || { echo "$PYTHON has no numpy/h5py (set PYTHON=...)"; exit 1; }
  echo "[verif] $(date '+%a %d %H:%M') steady sPTT pipe (3D) ..."
  OUT="$OUT/verification/sptt" NP=24 PARALLEL=1 EXE3D="$EXE3D" bash "$V/run_sptt_pipe.sh" < /dev/null
  echo "[verif] $(date '+%a %d %H:%M') pulsatile Oldroyd-B pipe, 3D and axisymmetric ..."
  OUT="$OUT/verification/womersley_3d" DRIVER=3d EXE="$EXE3D" NP=8 MAXJOBS=8 \
    bash "$V/run_womersley.sh" > "$LOG/verif_womersley_3d.log" 2>&1 < /dev/null &
  OUT="$OUT/verification/womersley_axi" DRIVER=axi EXE="$EXE" NP=2 MAXJOBS=8 \
    bash "$V/run_womersley.sh" > "$LOG/verif_womersley_axi.log" 2>&1 < /dev/null &
  wait
  echo "[verif] $(date '+%a %d %H:%M') finished. Tables:"
  echo "   $OUT/verification/sptt/table.txt"
  echo "   $OUT/verification/womersley_3d/tables.txt   $OUT/verification/womersley_axi/tables.txt"
  exit 0
fi

# ---- meshes and waveform ----------------------------------------------------------
need=()
case "$WHAT" in all|3d|verif3d) need+=(S75_pipe_inter A200_pipe_inter) ;; esac
case "$WHAT" in all|axi) need+=(S75_axi_medium S50_axi_medium H_axi_medium A125_axi_medium A150_axi_medium A200_axi_medium) ;; esac
case "$WHAT" in verif|verif3d)
  for rf in 0.7 1.0 1.4 2.0; do need+=("pipe_verification_rf$rf"); done ;; esac
for m in "${need[@]}"; do
  f="$MESH/$m.mesh"
  if [[ ! -f "$f" ]] || ! head -1 "$f" | grep -q MeshVersionFormatted; then
    echo "mesh $f missing or not MEDIT (Git LFS pointer? run: git lfs pull)"; exit 1
  fi
done
if [[ ! -f "$WF" ]]; then
  echo "waveform $WF not found"; exit 1
fi
if [[ "$WHAT" =~ ^(all|3d|verif3d)$ ]] && python3 -c "import numpy" 2> /dev/null; then
  echo "== checking the 3D meshes ..."
  python3 "$HERE/ArterialLesion3D_viscous_logarithmic_implicit/check_3d_mesh.py" \
    S75 "$MESH/S75_pipe_inter.mesh" A200 "$MESH/A200_pipe_inter.mesh" | tee "$LOG/check_3d_mesh.log"
  grep -q PROBLEMS "$LOG/check_3d_mesh.log" && { echo "3D mesh check failed"; exit 1; }
fi

# ---- free cores ---------------------------------------------------------------------
[[ "$WHAT" == verif3d ]] && stop3d
want=0
case "$WHAT" in all|3d) want=$(( want + 6 * NP3D )) ;; esac
case "$WHAT" in all|axi) want=$(( want + MAXAXI * NPAXI )) ;; esac
case "$WHAT" in verif|verif3d) want=$(( want + 80 )) ;; esac
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

# ---- verification, then (verif3d) the 3D cases ------------------------------------------
if [[ "$WHAT" == verif ]]; then
  echo "== verification (sPTT pipe, then Womersley 3D and axi) -> $OUT/verification"
  launch verif env SKIP_BUILD=1 bash "$0" _verif
elif [[ "$WHAT" == verif3d ]]; then
  echo "== verification, then the 3D cases -> $OUT/verification, $OUT/3D"
  launch verif bash -c "SKIP_BUILD=1 bash '$0' _verif; SKIP_BUILD=1 bash '$0' 3d"
fi

# ---- 3D -----------------------------------------------------------------------------
if [[ "$WHAT" == all || "$WHAT" == 3d ]]; then
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
if [[ "$WHAT" == all || "$WHAT" == axi ]]; then
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
echo "Status:    bash $0 status"
echo "Progress:  tail -f $LOG/axi_queue.log"
echo "3D ETA:    grep ETA $OUT/3D/Re300_Wo4_waveform_cca/S75/R/run.log | tail -1     (after ~20 min)"
echo "3D m=1:    the columns uyMean, uzMean of case.csv (seeded in cycle 1): growth = unstable"
echo "Summary:   cat $OUT/Re*/summary.csv $OUT/3D/Re*/summary.csv"
