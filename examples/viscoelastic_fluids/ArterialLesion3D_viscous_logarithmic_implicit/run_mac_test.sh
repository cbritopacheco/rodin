#!/usr/bin/env bash
# =============================================================================
#  run_mac_test.sh -- the hardest 3D case on a laptop (macOS or Linux), as an
#  A/B test of the stabilisation: S75, fluid R (the gsPTT blood model),
#  Re 600, Wo 5.5, sinusoidal inflow A = 0.5, BDF2, dynamic subscales.
#
#     cd ~/rodin && bash examples/viscoelastic_fluids/ArterialLesion3D_viscous_logarithmic_implicit/run_mac_test.sh
#
#  Variants (VARIANTS, run at the same time, the cores shared among them):
#     new     psi_rb, stress subscale residual without d_t psi_h
#             (-al_subscale_time_residual 0); refuted on xcoarse
#     old     psi_rb with d_t psi_h in it (the default again): ran the ramp
#             cycle on xcoarse, blew up at t = 0.14 s on layersc
#     nodiv   new without the div sigma subscale u3 (-al_stress_div_scale 0)
#     exp     the exp form (Moreno et al.), without d_t E_h in its subscale
#     expold  the exp form as before (-al_subscale_time_residual 1)
#     psi     the former term-by-term psi form
#     oldbdf1 old in BDF1 (TIME_ORDER=1): BDF1 in the field and in the subscale,
#             the scheme of the 2D tests (Kolmogorov, manufactured); with d_t
#             psi_h in the subscale residual the factor 1 - tau'/dt stays >= 0
#     D       new with fluid D (sPTT, Wi 2, beta 0.889): the same flow with
#             moderate elasticity (results in <variant>/Re600_Wo5.5/S75/Wi2)
#     R01     fluid R with lambda 0.1 s and eta_p 0.0104 (the largest with a
#             positive solvent): De 0.39, lambda gd ~ 1e3 at the wall, beta 0.81,
#             the same N1 at the wall (~40 Pa) as R at 4.83 s
#     O       fluid O, the Owens (2006) aggregating blood (run_3d.sh): lambda_H
#             5 ms in the jet (lambda gd ~ 25 at the throat wall), ~0.1 s in
#             the recirculation; results in <variant>/Re600_Wo5.5/S75/O
#  The flow peaks at t = T/2 = 0.128 s and the former failures came after it,
#  in the deceleration: one cycle (CYCLES=1, the ramp cycle) is enough.
#
#  It builds the 3D driver (SKIP_BUILD=1 to skip; re-runs CMake if the build
#  has no target for it), checks Scotch, runs in the background and, at the
#  end, writes summary.txt (per variant: last t, max |sigma|, retries, status)
#  and, with numpy and h5py, <variant>/blowup.txt (where |sigma|, |u| and psi
#  were extreme in the last outputs).
#
#  Environment (defaults):
#     VARIANTS="old oldbdf1 D"   MESH_LEVEL=xcoarse (22k tetrahedra;
#              layersxc 75k, layersc 121k: the throat wall and jet shear layer
#              refined)
#     NP       performance cores (sysctl hw.perflevel0.physicalcpu on macOS),
#              shared among the variants (at least 1 each)
#     STEPS=200 (steps per cycle)   CYCLES=1   OUTPUT_EVERY=5
#     SOLVER   lulag (MUMPS LU, lagged) with >= 16 GB of RAM, else gmres
#     BUILD=build   OUT=$HOME/Results/mac_test   MPIRUN=mpirun
# =============================================================================
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
BUILD="${BUILD:-build}"
OUT="${OUT:-$HOME/Results/mac_test}"
MPIRUN="${MPIRUN:-mpirun}"
STEPS="${STEPS:-200}"
CYCLES="${CYCLES:-1}"
OUTPUT_EVERY="${OUTPUT_EVERY:-5}"
VARIANTS="${VARIANTS:-old oldbdf1 D}"
MESH_LEVEL="${MESH_LEVEL:-xcoarse}"
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
nvar=$(echo $VARIANTS | wc -w | tr -d ' ')
NPV=$(( NP / nvar )); (( NPV < 1 )) && NPV=1
if [[ -z "${SOLVER:-}" ]]; then
  if (( memgb >= 16 )); then SOLVER=lulag; else SOLVER=gmres; fi
fi
echo "== $(uname -s), $ncores performance cores, ${memgb} GB: variants [$VARIANTS], $NPV rank(s) each, mesh $MESH_LEVEL, solver $SOLVER"

# ---- build and checks --------------------------------------------------------------
[[ -f "$REPO/$BUILD/CMakeCache.txt" ]] || { echo "no build in $REPO/$BUILD (set BUILD=...)"; exit 1; }
if ! grep -qiE '^RODIN_USE_SCOTCH:BOOL=(ON|YES|TRUE|Y|1)$' "$REPO/$BUILD/CMakeCache.txt"; then
  echo "warning: $BUILD is configured without Scotch; the partition will be poor"
  echo "         (cmake -S $REPO -B $REPO/$BUILD -DRODIN_USE_SCOTCH=ON)"
fi
if [[ "${SKIP_BUILD:-0}" != 1 ]]; then
  # A build configured on another branch has no target for this driver:
  # re-run CMake on the current sources first.
  if ! cmake --build "$REPO/$BUILD" --target help 2> /dev/null | grep -q ArterialLesion3D_viscous_logarithmic_implicit; then
    echo "== re-running CMake in $BUILD (it was configured without viscoelastic_fluids) ..."
    cmake -S "$REPO" -B "$REPO/$BUILD" > "$REPO/$BUILD/mac_test_cmake.log" 2>&1 \
      || { echo "CMAKE FAILED, see $REPO/$BUILD/mac_test_cmake.log"; tail -20 "$REPO/$BUILD/mac_test_cmake.log"; exit 1; }
  fi
  echo "== building the 3D driver (~6 GB of RAM per compiler process; -j 2) ..."
  cmake --build "$REPO/$BUILD" -j 2 --target ArterialLesion3D_viscous_logarithmic_implicit \
    > "$REPO/$BUILD/mac_test_build.log" 2>&1 \
    || { echo "BUILD FAILED, see $REPO/$BUILD/mac_test_build.log"; tail -20 "$REPO/$BUILD/mac_test_build.log"; exit 1; }
fi
[[ -x "$EXE" ]] || { echo "missing $EXE"; exit 1; }
grep -aq -- "-al_subscale_time_residual" "$EXE" || { echo "$EXE is older than the sources: rebuild (without SKIP_BUILD)"; exit 1; }
MESH="$MESHDIR/S75_pipe_${MESH_LEVEL}.mesh"
if [[ ! -f "$MESH" ]] || ! head -1 "$MESH" | grep -q MeshVersionFormatted; then
  case "$MESH_LEVEL" in
    layersc)  margs=(0.7 -setnumber layers 1) ;;
    layersxc) margs=(0.7 -setnumber layers 1 -setnumber hb 0.3 -setnumber hwk 0.16 -setnumber Lwk 4
                     -setnumber Lsl 3 -setnumber hwl 0.2 -setnumber hwmax 0.14 -setnumber Ld 15) ;;
    *) echo "mesh $MESH not found"; exit 1 ;;
  esac
  echo "== making $MESH (python gmsh) ..."
  "$PY" "$HERE/make_lesion_mesh_3d.py" S75 "${margs[0]}" "$MESH" "${margs[@]:1}" || exit 1
fi

# ---- run ----------------------------------------------------------------------------
mkdir -p "$OUT/$MESH_LEVEL"
PIDS=""
for V in $VARIANTS; do
  TO=2; FL=R; LR=4.83; EP=0.050
  case "$V" in
    new)    FORM=psi_rb; EXTRA="-al_subscale_time_residual 0" ;;
    old)    FORM=psi_rb; EXTRA="-al_subscale_time_residual 1" ;;
    nodiv)  FORM=psi_rb; EXTRA="-al_stress_div_scale 0 -al_subscale_time_residual 0" ;;
    exp)    FORM=exp;    EXTRA="-al_subscale_time_residual 0" ;;
    expold) FORM=exp;    EXTRA="-al_subscale_time_residual 1" ;;
    psi)    FORM=psi;    EXTRA="" ;;
    oldbdf1) FORM=psi_rb; EXTRA="-al_subscale_time_residual 1"; TO=1 ;;
    D)      FORM=psi_rb; EXTRA="-al_subscale_time_residual 1"; FL=D ;;
    R01)    FORM=psi_rb; EXTRA="-al_subscale_time_residual 1"; LR=0.1; EP=0.0104 ;;
    O)      FORM=psi_rb; EXTRA="-al_subscale_time_residual 1"; FL=O ;;
    *) echo "unknown variant $V"; exit 1 ;;
  esac
  rm -rf "$OUT/$MESH_LEVEL/$V"
  (
    env FLUIDS="$FL" WI=2 LAMBDA_R="$LR" ETAP_R="$EP" NP="$NPV" TIME_ORDER="$TO" MESH_LEVEL="$MESH_LEVEL" OUTPUT_EVERY="$OUTPUT_EVERY" STEPS="$STEPS" \
      CYCLES_R="$CYCLES" CYCLES="$CYCLES" SOLVER="$SOLVER" FORMULATION="$FORM" PETSC_EXTRA="$EXTRA -al_trace_newton 1" \
      EXE3D="$EXE" MPIRUN="$MPIRUN" OUT3D="$OUT/$MESH_LEVEL/$V" RETRY=0 bash "$HERE/run_3d.sh" 600 S75
    d="$OUT/$MESH_LEVEL/$V/Re600_Wo5.5/S75/R"; [[ "$FL" == D ]] && d="$OUT/$MESH_LEVEL/$V/Re600_Wo5.5/S75/Wi2"
    [[ "$FL" == O ]] && d="$OUT/$MESH_LEVEL/$V/Re600_Wo5.5/S75/O"
    if "$PY" -c "import numpy, h5py" 2> /dev/null && [[ -f "$d/case.xdmf" ]]; then
      "$PY" "$HERE/locate_blowup.py" "$d" --last 8 > "$OUT/$MESH_LEVEL/$V/blowup.txt" 2>&1
    fi
  ) > "$OUT/$MESH_LEVEL/$V.log" 2>&1 < /dev/null &
  PIDS="$PIDS $!"
  echo "   started $V ($FORM $EXTRA) -> $OUT/$MESH_LEVEL/$V/Re600_Wo5.5/S75/R"
  sleep 2
done
(
  for p in $PIDS; do while kill -0 "$p" 2> /dev/null; do sleep 30; done; done
  {
    echo "variant      status  t_end      max|sigma|  max|u|   retries  Newton(max)"
    for V in $VARIANTS; do
      d="$OUT/$MESH_LEVEL/$V/Re600_Wo5.5/S75/R"; [[ "$V" == D ]] && d="$OUT/$MESH_LEVEL/$V/Re600_Wo5.5/S75/Wi2"
      [[ "$V" == O ]] && d="$OUT/$MESH_LEVEL/$V/Re600_Wo5.5/S75/O"
      st=RUNNING; [[ -f "$d/DONE" ]] && st=DONE; [[ -f "$d/FAILED" ]] && st=FAILED
      awk -F, -v v="$V" -v st="$st" 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
        { t = $c["t"]; s = $c["maxSigma"] + 0; u = $c["maxU"] + 0; r += $c["retries"];
          if (s > sm) sm = s; if (u > um) um = u; if ($c["newtonIts"] + 0 > nm) nm = $c["newtonIts"] + 0 }
        END { printf "%-12s %-7s %-10.5f %-11.4g %-8.3g %-8d %d\n", v, st, t, sm, um, r, nm }' "$d/case.csv" 2> /dev/null \
        || echo "$V $st (no case.csv)"
    done
  } > "$OUT/$MESH_LEVEL/summary.txt"
  echo "== mac test finished: $(date)" >> "$OUT/$MESH_LEVEL/summary.txt"
) > /dev/null 2>&1 < /dev/null &
echo "== follow:  tail -n 3 $OUT/$MESH_LEVEL/*/Re600_Wo5.5/S75/R/case.csv"
echo "   end:     $OUT/$MESH_LEVEL/summary.txt and $OUT/$MESH_LEVEL/<variant>/blowup.txt"
echo "   stop:    pkill -f ArterialLesion3D_viscous_logarithmic_implicit"
