#!/usr/bin/env bash
# Constant-resistance twin of run_scenarios.sh.
#
# Same scenarios, same 3D Carreau-Yasuda rheology, same calibration -- the ONLY
# difference is -coronary_constant_outlet, which freezes the reduced outlet
# closure at its high-shear plateau:
#
#     mu_ap(gammadot) == mu_inf   =>   R_inf = 8 mu_inf L / (N pi R^4)
#
# i.e. "same mu_inf, zero shear-thinning". With mu constant the WRMS integral
# is I = tau_w^4/(4 mu) exactly, so every outlet becomes a linear resistor
# while the resolved 3D field keeps its full shear-dependent viscosity. Each
# run here is the fixed-state twin of the corresponding run of
# run_scenarios.sh, and the pair isolates the reduced closure -- exactly the
# comparison of section 7.5 of the manuscript.
#
# As in run_scenarios.sh the BASELINE case is not re-run (its constant-closure
# counterpart is the mu_inf run already reported); uncomment the line below to
# compute it. Output goes to <scenario>_const/ so nothing collides with the
# R_mu runs. Runs are SEQUENTIAL. Same environment knobs as run_scenarios.sh:
#   bash run_scenarios_constant.sh                       # serial, no mpirun
#   NP=8 / LAUNCH="mpiexec -n 4" / THREADS=8 / BIN=... / PETSC_OPTS / EXTRA
set -euo pipefail

BIN=${BIN:-./examples/Heart/CoronaryArtery}
NP=${NP:-1}
THREADS=${THREADS:-40}
# Same Q0 reference as the R_mu family, so both families are reported on one
# shear-rate scale.
BASELINE_CSV=${BASELINE_CSV:-hyp2/CoronaryArtery.csv}
PETSC_OPTS=${PETSC_OPTS:-"-ksp_type preonly -pc_type lu \
-pc_factor_mat_solver_type mumps -mat_mumps_icntl_20 0 -mat_mumps_icntl_21 0 \
-ksp_converged_reason -snes_monitor -snes_converged_reason"}
EXTRA=${EXTRA:-}

if [ ! -x "$BIN" ]; then
  echo "error: binary '$BIN' not found or not executable (run from the build" >&2
  echo "       directory, or set BIN=...)" >&2
  exit 1
fi

if [ -z "${LAUNCH:-}" ]; then
  if [ "$NP" -gt 1 ] && command -v mpirun >/dev/null 2>&1; then
    LAUNCH="mpirun -n $NP"
  elif [ "$NP" -gt 1 ] && command -v mpiexec >/dev/null 2>&1; then
    LAUNCH="mpiexec -n $NP"
  else
    if [ "$NP" -gt 1 ]; then
      echo "warning: no mpirun/mpiexec found; running serially" >&2
    fi
    LAUNCH=""
  fi
fi

run () {
  local name=$1; shift
  mkdir -p "$name"
  echo "==> $name  (constant outlet closure)"
  VECLIB_MAXIMUM_THREADS=$THREADS OPENBLAS_NUM_THREADS=$THREADS \
  OMP_NUM_THREADS=$THREADS \
  $LAUNCH "$BIN" -coronary_output_prefix "$name" -coronary_constant_outlet 1 \
    "$@" $PETSC_OPTS $EXTRA
}

# run baseline_const                 # already covered by the mu_inf run
run pra2200_const  -coronary_operating_pra 2200
run alpha085_const -coronary_alpha_im 0.85
run combined_const -coronary_operating_pra 2200 -coronary_alpha_im 0.85

HERE="$(dirname "$0")"

# ---- cycle metrics of the constant-closure family ---------------------------
ARGS=()
if [ -f "$BASELINE_CSV" ]; then
  ARGS+=(--baseline "$BASELINE_CSV")
else
  echo
  echo "note: BASELINE_CSV '$BASELINE_CSV' not found; pra2200_const is used as" >&2
  echo "      the Q0 reference." >&2
  ARGS+=(--baseline pra2200_const/CoronaryArtery.csv)
fi

python3 "$HERE/retrograde_metrics.py" \
  --period 0.85 --csv-out scenario_metrics_const.csv \
  "${ARGS[@]}" \
  pra2200=pra2200_const/CoronaryArtery.csv \
  alpha085=alpha085_const/CoronaryArtery.csv \
  combined=combined_const/CoronaryArtery.csv

# ---- paired R_mu vs constant comparison -------------------------------------
if [ -f pra2200/CoronaryArtery.csv ]; then
  echo
  python3 "$HERE/compare_closures.py" --period 0.85 \
    --csv-out closure_comparison.csv \
    pra2200=pra2200/CoronaryArtery.csv,pra2200_const/CoronaryArtery.csv \
    alpha085=alpha085/CoronaryArtery.csv,alpha085_const/CoronaryArtery.csv \
    combined=combined/CoronaryArtery.csv,combined_const/CoronaryArtery.csv
else
  echo
  echo "note: R_mu runs not found; run run_scenarios.sh first, then"
  echo "      python3 compare_closures.py name=rmu.csv,const.csv ..."
fi
