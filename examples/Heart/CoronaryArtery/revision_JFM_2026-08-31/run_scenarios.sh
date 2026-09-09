#!/usr/bin/env bash
# Pathological-scenario protocol (R_mu closure): venous hypertension /
# pressure overload / combined. Runs are SEQUENTIAL (each waits for the
# previous one); each run writes ALL its output -- CSV and 3D XDMF/HDF5
# fields alike -- into its own directory via -coronary_output_prefix, so
# nothing overwrites hyp2/ or another scenario. The 3D Carreau-Yasuda
# rheology is identical across the runs by design: only the operating
# drainage pressure and the intramyocardial transmission change.
#
# The BASELINE case is deliberately NOT re-run: the reference solution
# already exists (BASELINE_CSV below, default hyp2/CoronaryArtery.csv) and is
# used as the Q0 reference of the metrics so every scenario is reported on
# the same shear-rate scale. To compute it again, uncomment the `run baseline`
# line and point BASELINE_CSV at baseline/CoronaryArtery.csv.
#
# Physiological anchors:
#   pra2200  : central venous hypertension ~16.5 mmHg (right-heart failure,
#              tricuspid regurgitation, constrictive pericarditis)
#   alpha085 : pressure-overload / subendocardial territory, p_im = 0.85 p_LV
#   combined : both (biventricular failure, worst case)
#
# Run it from the build directory:
#   bash run_scenarios.sh                        # serial, no mpirun
#   NP=8 bash run_scenarios.sh                   # mpirun -n 8 (if available)
#   LAUNCH="mpiexec -n 4" bash run_scenarios.sh  # custom launcher
#   THREADS=8 bash run_scenarios.sh              # OpenMP threads (default 40)
#   BIN=./examples/Heart/CoronaryArtery ...      # binary location
#   BASELINE_CSV=hyp2/CoronaryArtery.csv ...     # existing reference run
#   PETSC_OPTS="..." / EXTRA="..."               # solver flags / extra args
set -euo pipefail

BIN=${BIN:-./examples/Heart/CoronaryArtery}
NP=${NP:-1}
THREADS=${THREADS:-40}
BASELINE_CSV=${BASELINE_CSV:-hyp2/CoronaryArtery.csv}
# The driver already installs these as defaults; passing them explicitly keeps
# the invocation identical to the one used interactively.
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
  echo "==> $name"
  VECLIB_MAXIMUM_THREADS=$THREADS OPENBLAS_NUM_THREADS=$THREADS \
  OMP_NUM_THREADS=$THREADS \
  $LAUNCH "$BIN" -coronary_output_prefix "$name" "$@" $PETSC_OPTS $EXTRA
}

# run baseline                       # already computed; see BASELINE_CSV
run pra2200  -coronary_operating_pra 2200
run alpha085 -coronary_alpha_im 0.85
run combined -coronary_operating_pra 2200 -coronary_alpha_im 0.85

# ---- cycle metrics ----------------------------------------------------------
ARGS=()
if [ -f "$BASELINE_CSV" ]; then
  ARGS+=(--baseline "$BASELINE_CSV" "baseline=$BASELINE_CSV")
else
  echo
  echo "note: BASELINE_CSV '$BASELINE_CSV' not found; pra2200 is used as the" >&2
  echo "      Q0 reference and no baseline column is reported." >&2
  ARGS+=(--baseline pra2200/CoronaryArtery.csv)
fi

python3 "$(dirname "$0")/retrograde_metrics.py" \
  --period 0.85 --csv-out scenario_metrics.csv \
  "${ARGS[@]}" \
  pra2200=pra2200/CoronaryArtery.csv \
  alpha085=alpha085/CoronaryArtery.csv \
  combined=combined/CoronaryArtery.csv
