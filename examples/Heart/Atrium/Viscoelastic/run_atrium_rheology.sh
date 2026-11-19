#!/usr/bin/env bash
# =============================================================================
#  run_atrium_rheology.sh -- blood rheology in the left atrium (fixed fluid
#  domain of the Atrium FSI mesh): the same geometry, openings and numerics
#  for four fluids that differ by one ingredient each.
#
#     N    Newtonian, eta = eta_inf of the Carreau-Yasuda law          (u, p)
#     CY   Carreau-Yasuda                                              (u, p)
#     VE   Owens with equilibrium structure, n = n_st(gd): elasticity,
#          lambda = (eta_CY(gd) - eta_plasma)/G, no structural memory  (u, p, psi)
#     TVE  Owens (2006): n transported, rates a(gd), b(gd): elasticity
#          and thixotropy                                              (u, p, psi, n)
#     R    optional: the gsPTT of the arterial campaign (constant lambda)
#
#  CY, VE and TVE have the same steady shear viscosity, eta_CY(gd), exactly,
#  and VE and TVE the same modulus G = (eta_inf - eta_plasma)/lambda_H: CY-N
#  isolates shear thinning, VE-CY elasticity, TVE-VE the memory of the
#  structure. The kinetics are those of blood (Moyers-Gonzalez, Owens & Fang,
#  JFM 617, 2008: a = 1 1/s at rest, 2 1/s at 5.78 1/s, 0 above 900 1/s;
#  lambda_H = 5 ms): no thixotropy number is swept. The relaxation time of the
#  structure is not one number: n_st/a ~ 19 s at rest, 0.5 s at 10 1/s, 0.1 s
#  at 100 1/s, against T = 0.85 s (the driver prints the table).
#
#  The cases run one after another with NP ranks each (default 60), or JOBS at
#  a time. A finished case leaves DONE and is skipped, so the script can be
#  re-run after an interruption. When every case is done postprocess_atrium.py
#  writes the comparison to $OUT/post.
#
#  Usage (from the repository root, on pc9):
#     nohup bash examples/Heart/Atrium/Viscoelastic/run_atrium_rheology.sh > ~/atrium_rheology.log 2>&1 &
#     bash examples/Heart/Atrium/Viscoelastic/run_atrium_rheology.sh status
#     bash examples/Heart/Atrium/Viscoelastic/run_atrium_rheology.sh stop
#     bash examples/Heart/Atrium/Viscoelastic/run_atrium_rheology.sh post
#
#  Environment (defaults):
#     OUT=$HOME/Results/atrium_rheology   BUILD=build-gcc12   NP=60   JOBS=1
#     MESH=resources/examples/Heart/Atrium/AtriumFine_FluidSolid_1mm.mesh
#                    the fluid mesh of the FSI case (mm; fluid 1, solid 4). It
#                    is stored in Git LFS: on pc9 copy the real file and point
#                    MESH at it (a pointer file is rejected below).
#     FLUIDS="TVE VE CY N"     the longest first; add R for the gsPTT
#     CYCLES=6  CYCLES_TVE=8   cycles (the first with the pressure ramp)
#     STEPS=850                steps per cycle (dt = 1 ms at T = 0.85 s)
#     TIME_ORDER=2  OUTPUT_EVERY=17   BDF2; 50 snapshots per cycle
#     LV=1  LV_FEEDBACK=1      lumped ventricle on the mitral plane, coupled
#                              valve (LV=0: P_OUT, R_OUT prescribed)
#     P_IN=1200 P_IN_AMP=200   venous pressure (Pa), as in Atrium/Fluid
#     LAA_CENTER=  LAA_RADIUS= ball of the appendage, metres (the coordinates
#                              of the XDMF output: pick the centre in ParaView),
#                              as LAA_CENTER=x,y,z LAA_RADIUS=r. Empty: no
#                              appendage diagnostics.
#     PROBES=                  "x,y,z;x,y,z" in metres; PROBE_RADIUS=1.5e-3
#     CY_SET=atrium            atrium (Atrium/Fluid: 0.072, 0.0049, 6.1529,
#                              0.2081, 1.4173) or chokensey (0.056, 0.00345,
#                              1.902, 0.22, 1.25, the arterial campaign)
#     LAMBDA_H=0.005 ETA_PLASMA=1.2e-3 A0=1 N_IN=1.5 RESET_CYCLE=2   Owens
#     EPS_O=0                  linear PTT factor of the Owens dumbbells (VE and
#                              TVE). 0 is the published model. A positive value
#                              bounds the extensional stress, but the steady
#                              viscosity of VE and TVE then falls below CY.
#     LAMBDA_R=0.114 ETAP_R=0.0068 EPS_R=0.5   gsPTT (R), as launch_aneurysm_v3.sh
#     NEWTON_ITS=8 NEWTON_TOL=5e-4 PSI_ORDER=2 PC_REUSE=20
#     SOLVER=lulag             lulag (GMRES + lagged MUMPS LU), gmres (ASM/ILU)
#                              or mumps (direct)
#     PETSC_EXTRA=""           extra options for the driver
#     SKIP_BUILD=0  MPIRUN=mpirun
# =============================================================================
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"
OUT="${OUT:-$HOME/Results/atrium_rheology}"
BUILD="${BUILD:-build-gcc12}"
[[ "$BUILD" = /* ]] || BUILD="$REPO/$BUILD"
NP="${NP:-60}"
JOBS="${JOBS:-1}"
MPIRUN="${MPIRUN:-mpirun}"
MESH="${MESH:-$REPO/resources/examples/Heart/Atrium/AtriumFine_FluidSolid_1mm.mesh}"
FLUIDS="${FLUIDS:-TVE VE CY N}"
CYCLES="${CYCLES:-6}"
CYCLES_TVE="${CYCLES_TVE:-8}"
STEPS="${STEPS:-850}"
TIME_ORDER="${TIME_ORDER:-2}"
OUTPUT_EVERY="${OUTPUT_EVERY:-17}"
LV="${LV:-1}"
LV_FEEDBACK="${LV_FEEDBACK:-1}"
P_IN="${P_IN:-1200}"
P_IN_AMP="${P_IN_AMP:-200}"
P_OUT="${P_OUT:-0}"
R_OUT="${R_OUT:-0}"
LAA_CENTER="${LAA_CENTER:-}"
LAA_RADIUS="${LAA_RADIUS:-}"
PROBES="${PROBES:-}"
PROBE_RADIUS="${PROBE_RADIUS:-1.5e-3}"
CY_SET="${CY_SET:-atrium}"
LAMBDA_H="${LAMBDA_H:-0.005}"
ETA_PLASMA="${ETA_PLASMA:-1.2e-3}"
A0="${A0:-1}"
N_IN="${N_IN:-1.5}"
RESET_CYCLE="${RESET_CYCLE:-2}"
EPS_O="${EPS_O:-0}"
LAMBDA_R="${LAMBDA_R:-0.114}"
ETAP_R="${ETAP_R:-0.0068}"
EPS_R="${EPS_R:-0.5}"
NEWTON_ITS="${NEWTON_ITS:-8}"
NEWTON_TOL="${NEWTON_TOL:-5e-4}"
PSI_ORDER="${PSI_ORDER:-2}"
PC_REUSE="${PC_REUSE:-20}"
SOLVER="${SOLVER:-lulag}"
PETSC_EXTRA="${PETSC_EXTRA:-}"
EXE="$BUILD/examples/Heart/Atrium/Viscoelastic/AtriumViscoelastic"
POST="$HERE/postprocess_atrium.py"
LOG="$OUT/logs"
PIDFILE="$LOG/launcher.pid"
GMRES_OPTS="-ksp_type gmres -ksp_gmres_restart 200 -ksp_max_it 2000 -ksp_rtol 1e-8 -pc_type asm -pc_asm_overlap 1 -sub_pc_type ilu -sub_pc_factor_levels 1"
LULAG_OPTS="-ksp_type gmres -ksp_gmres_restart 100 -ksp_max_it 300 -ksp_rtol 1e-8 -pc_type lu -pc_factor_mat_solver_type mumps"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export OMPI_MCA_hwloc_base_binding_policy=none

case "$CY_SET" in
  atrium)    CY=(-av_cy_eta_zero 0.0720 -av_cy_eta_inf 0.0049 -av_cy_lambda 6.1529 -av_cy_n 0.2081 -av_cy_a 1.4173); ETA_INF=0.0049 ;;
  chokensey) CY=(-av_cy_eta_zero 0.056 -av_cy_eta_inf 0.00345 -av_cy_lambda 1.902 -av_cy_n 0.22 -av_cy_a 1.25); ETA_INF=0.00345 ;;
  *) echo "CY_SET must be atrium or chokensey"; exit 1 ;;
esac

status() {
  echo "== $OUT"
  for F in $FLUIDS; do
    d="$OUT/$F"; s="-"
    [[ -f "$d/run.log" ]] && s="run $(grep -o 'Step [0-9]*/[0-9]*' "$d/run.log" | tail -1 | cut -d' ' -f2)"
    [[ -f "$d/FAILED" ]] && s="FAILED ($(cat "$d/FAILED"))"
    [[ -f "$d/DONE" ]] && s="DONE"
    printf "%-4s %s\n" "$F" "$s"
  done
}
post() {
  local args=("$OUT" --fields)
  [[ -n "$LAA_CENTER" && -n "$LAA_RADIUS" ]] && args+=(--appendage "$LAA_CENTER,$LAA_RADIUS")
  python3 "$POST" "${args[@]}" > "$LOG/postprocess.log" 2>&1 \
    && tail -n +1 "$LOG/postprocess.log" | tail -12 \
    || echo "postprocess failed: $LOG/postprocess.log"
}
mkdir -p "$LOG"
case "${1:-start}" in
  status) status; exit 0 ;;
  post) post; exit 0 ;;
  stop)
    [[ -f "$PIDFILE" ]] && pkill -TERM -P "$(cat "$PIDFILE")" 2>/dev/null
    [[ -f "$PIDFILE" ]] && kill "$(cat "$PIDFILE")" 2>/dev/null
    pkill -u "$USER" -f "$EXE" 2>/dev/null; rm -f "$PIDFILE"; echo stopped; exit 0 ;;
  start) ;;
  *) echo "usage: $0 [start|status|stop|post]"; exit 1 ;;
esac
if [[ -f "$PIDFILE" ]] && kill -0 "$(cat "$PIDFILE")" 2>/dev/null; then
  echo "a launcher is already running (pid $(cat "$PIDFILE"))"; exit 1
fi
echo $$ > "$PIDFILE"
trap 'rm -f "$PIDFILE"' EXIT

if [[ ! -f "$MESH" ]]; then
  echo "Mesh not found: $MESH"; exit 1
fi
if ! head -1 "$MESH" | grep -q MeshVersionFormatted; then
  echo "$MESH is not a MEDIT mesh (a Git LFS pointer?): copy the real file and set MESH."; exit 1
fi
MESH="$(cd "$(dirname "$MESH")" && pwd)/$(basename "$MESH")"
if [[ "${SKIP_BUILD:-0}" != 1 ]]; then
  cmake --build "$BUILD" -j 16 --target AtriumViscoelastic > "$LOG/build.log" 2>&1 \
    || { echo "build failed: $LOG/build.log"; exit 1; }
fi
[[ -x "$EXE" ]] || { echo "executable missing: $EXE"; exit 1; }

run_case() {   # fluid
  local F="$1" dir="$OUT/$1" cycles="$CYCLES" fargs
  local owens=(-av_fluid owens -av_owens_lambda_h "$LAMBDA_H" -av_owens_eta_plasma "$ETA_PLASMA"
               -av_owens_a0 "$A0" -av_owens_n_inlet "$N_IN" -av_owens_epsilon "$EPS_O")
  case "$F" in
    N)   fargs=(-av_fluid newtonian -av_eta_ref "$ETA_INF") ;;
    CY)  fargs=(-av_fluid cy) ;;
    VE)  fargs=("${owens[@]}" -av_owens_equilibrium 1) ;;
    TVE) fargs=("${owens[@]}" -av_owens_equilibrium 0 -av_owens_reset_cycle "$RESET_CYCLE"); cycles="$CYCLES_TVE" ;;
    R)   fargs=(-av_fluid gsptt -av_eta_p "$ETAP_R" -av_ptt_epsilon "$EPS_R" -av_lambda "$LAMBDA_R") ;;
    *) echo "unknown fluid $F (N, CY, VE, TVE or R)"; return 2 ;;
  esac
  if [[ -f "$dir/DONE" ]]; then
    echo "[skip] $F (DONE)"; return 0
  fi
  mkdir -p "$dir"; rm -f "$dir/FAILED"
  local cmd=("$MPIRUN" -n "$NP" "$EXE"
    -av_mesh "$MESH" -av_xdmf case -av_csv case.csv
    "${fargs[@]}" "${CY[@]}"
    -av_p_in "$P_IN" -av_p_in_amplitude "$P_IN_AMP"
    -av_lv_model "$LV" -av_lv_feedback "$LV_FEEDBACK" -av_p_out "$P_OUT" -av_r_out "$R_OUT"
    -av_time_order "$TIME_ORDER" -av_cycles "$cycles" -av_steps_per_cycle "$STEPS"
    -av_output_every "$OUTPUT_EVERY"
    -av_conformation_its "$NEWTON_ITS" -av_conformation_tol "$NEWTON_TOL" -av_psi_order "$PSI_ORDER"
    -av_formulation psi_rb -av_dynamic_subscales 1 -av_adaptive_dt 1 -av_newton_predict 1
    -av_pc_reuse_solves "$PC_REUSE" -av_residence_time 1 -av_metrics 1)
  [[ -n "$LAA_CENTER" && -n "$LAA_RADIUS" ]] && cmd+=(-av_appendage_center "$LAA_CENTER" -av_appendage_radius "$LAA_RADIUS")
  [[ -n "$PROBES" ]] && cmd+=(-av_probes "$PROBES" -av_probe_radius "$PROBE_RADIUS")
  # shellcheck disable=SC2206
  [[ "$SOLVER" == gmres ]] && cmd+=($GMRES_OPTS)
  # shellcheck disable=SC2206
  [[ "$SOLVER" == lulag ]] && cmd+=($LULAG_OPTS)
  # shellcheck disable=SC2206
  [[ -n "$PETSC_EXTRA" ]] && cmd+=($PETSC_EXTRA)
  printf '%q ' "${cmd[@]}" > "$dir/params.txt"; echo >> "$dir/params.txt"

  echo "[run ] $(date '+%a %d %H:%M')  $F  cycles=$cycles steps/cycle=$STEPS BDF$TIME_ORDER NP=$NP $SOLVER -> $dir"
  local t0=$SECONDS
  ( cd "$dir" && "${cmd[@]}" > run.log 2>&1 )
  local st=$? minutes=$(( (SECONDS - t0) / 60 ))
  if [[ $st -eq 0 ]]; then
    touch "$dir/DONE"; echo "[done] $F in $minutes min"
  else
    echo "exit $st" > "$dir/FAILED"
    echo "[FAIL] $F (exit $st, $minutes min): $(grep -m1 -E 'diverged|Fatal|Error|error' "$dir/run.log")"
  fi
  return $st
}

echo "[$(date '+%a %d %H:%M')] fluids: $FLUIDS  $JOBS at a time x $NP ranks  mesh: $MESH  ->  $OUT"
running=0
for F in $FLUIDS; do
  if (( running >= JOBS )); then wait -n; running=$(( running - 1 )); fi
  run_case "$F" &
  running=$(( running + 1 ))
done
wait
echo "[$(date '+%a %d %H:%M')] all cases ended; postprocessing"
post
status
