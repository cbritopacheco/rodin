#!/usr/bin/env bash
# =============================================================================
#  run_3d.sh -- pulsatile 3D cases of ArterialLesion3D at ONE Reynolds number
#  and ONE geometry, for the fluids N and D of the axisymmetric study, to test
#  whether the axisymmetric base state survives non-axisymmetric (m >= 1)
#  perturbations.
#
#     ./run_3d.sh 600 S75            # N then D, 30 ranks each
#     FLUIDS=D ./run_3d.sh 600 S75   # only D (run N in another terminal)
#
#  Fluids (same parameters as run_axi_re.sh):
#     N   Newtonian: the sPTT at Wi = 0.01                 -> <geo>/N/
#     D   weak-elasticity sPTT (beta 0.889, eps 0.1), Wi    -> <geo>/Wi<wi>/
#     CY  Carreau-Yasuda of whole blood (eta_ref 3.45 mPa s) -> <geo>/CY/
#     R   generalised sPTT, the viscoelastic twin of CY      -> <geo>/R/
#
#  Seed: the inlet profile is multiplied by 1 + PERTURB y/R during the first
#  cycle (flow rate unchanged), then the flow evolves freely. The CSV columns
#  uyMean, uzMean (volume means of u_y, u_z over Ubar) vanish for any
#  axisymmetric flow: if they decay after cycle 1 the axisymmetric state is
#  stable; if they grow, the m = 1 mode is unstable.
#
#  Output: $OUT/Re<Re>_Wo<Wo>/<geometry>/<fluid>/
#     run.log  case.csv  case.xdmf/h5  params.txt  DONE | FAILED
#     $OUT/Re<Re>_Wo<Wo>/summary.csv
#
#  Options (environment variables):
#     FLUIDS="CY R"  WI=2   NP=30
#     TIME_ORDER=2    BDF2 (1: BDF1, the scheme of the first campaign)
#     WO=<n>          default by Re: 100 -> 2.5, 300 -> 4, 600 -> 5.5
#     STEPS=<n>       steps per cycle; default 600 (Re 100), 1000 (Re 300, 600):
#                     U_jet dt/d_throat ~ 0.25 for S75
#     CYCLES=4        cycle 1: ramp + seed; cycles 2-4 free
#     CYCLES_R=6      cycles for R (lambda spans several periods: compare the
#                     growth of uyMean, uzMean with CY rather than periodicity)
#     ETAP_R=0.050    polymer viscosity of R (Pa s); the solvent eta_CY - eta_p/f
#                     must stay positive: eta_p <= 0.0104 at lambda 0.1 s, 0.0081
#                     at 0.05 s, 0.031 at 1 s (N1 at the wall ~40-60 Pa either way)
#     LAMBDA_R=4.83   relaxation time of R (s), a material constant, the same
#                     at every Re (De = lambda/T varies with Wo at fixed
#                     El = lambda nu/D^2); solvent positive for lambda >= 2.4 s
#     WAVEFORM=<file> physiological flow-rate waveform instead of the sinusoid
#                     (output Re<Re>_Wo<Wo>_<file name>), HARMONICS=10
#     PERTURB=0.01    seed amplitude delta (0: no seed)
#     AMP=0.5  PTT_EPS=0.1  OUTPUT_EVERY=40  MESH_LEVEL=inter  NEWTON_ITS=8
#                     (inter: S75 rf 0.85, 328k tets; A200 rf 1.1, 202k tets;
#                      about half the cost of medium)
#     PSI_ORDER=2     quadrature degree of the log-conformation integrals
#                     (2: 4 points per tetrahedron; 4: 11 points, slower)
#     NEWTON_TOL=5e-4 Newton stops when max|dpsi|/max(1,|psi|) < NEWTON_TOL;
#                     Newton is quadratic, so the error left is ~NEWTON_TOL^2
#     FORMULATION=psi_rb constitutive discretisation (-al_formulation); all
#                     three have split OSS on momentum/continuity:
#                       psi_rb psi equation, OSS of its full residual (default)
#                       exp    Moreno et al. (2019, 2020): sigma equation in
#                              E = exp(psi), OSS of its full residual; loses
#                              the small eigenvalues of psi at blood Wi
#                       psi    the former psi equation with term-by-term OSS
#     DYNAMIC=1       dynamic subscales (psi_rb, exp; 0: quasi-static)
#     ADAPTIVE_DT=1   step control: a failed step (Newton not converged or
#                     diverging, linear solver failure, velocity guard) is
#                     retried with dt/2 (up to 6 times); dt grows back x1.5
#                     after 10 easy steps, never above T/STEPS (-al_adaptive_dt)
#     NEWTON_PREDICT=1 accept the Newton iterate when the quadratic prediction
#                     C d_k^2 of the next increment is below NEWTON_TOL/10,
#                     saving the confirming iteration (-al_newton_predict)
#     PC_REUSE=20     lagged preconditioner: reuse the ASM/ILU (or LU) factors
#                     for up to PC_REUSE solves; rebuilt when a solve needs
#                     more than PC_REUSE_GROWTH=2 times the iterations of the
#                     first one after the build, when one fails (repeated
#                     fresh), or when dt changes; 1 rebuilds at every solve
#     OUTLET_PSI=1    weak inflow condition for psi where the outlet has
#                     backflow (-al_outlet_conformation_inflow)
#     PSI_SKEW=0      skew-symmetric psi advection, + (1/2)(div u) psi (-al_psi_skew)
#     PETSC_EXTRA=""  extra PETSc options, e.g. MUMPS block low-rank:
#                     "-mat_mumps_icntl_35 1 -mat_mumps_cntl_7 1e-10"
#     SOLVER=gmres    linear solver of each Newton iteration:
#                       gmres  GMRES + additive Schwarz (overlap 1) + ILU(1)
#                       lulag  GMRES + MUMPS LU, lagged over PC_REUSE solves
#                       mumps  direct LU (MUMPS), the driver's default
#     RETRY=1         a failed case is re-run: with gmres, first with mumps
#                     (<case>_lu); then with mumps and 2x STEPS (<case>_dt2)
#     OUT3D=<dir>     default $OUT/3D (OUT as for run_axi_re.sh: $HOME/Results/axi_jfm)
#     MPIRUN=mpirun   EXE3D=...   MESHDIR=...
#
#  Usage: ./run_3d.sh <Re> <geometry>
# =============================================================================
set -u

if [[ $# -lt 2 ]]; then
  sed -n '2,40p' "$0"; exit 1
fi
RE="$1"; GEO="$2"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"

EXE="${EXE3D:-$REPO/build/examples/viscoelastic_fluids/ArterialLesion3D_viscous_logarithmic_implicit/ArterialLesion3D_viscous_logarithmic_implicit}"
MESHDIR="${MESHDIR:-$REPO/resources/examples/viscoelastic_fluids}"
MPIRUN="${MPIRUN:-mpirun}"
NP="${NP:-30}"
FLUIDS="${FLUIDS:-CY R}"
WI="${WI:-2}"
if [[ -z "${WO:-}" ]]; then
  case "$RE" in
    100) WO=2.5 ;; 300) WO=4 ;; 600) WO=5.5 ;; *) WO=4 ;;
  esac
fi
if [[ -z "${STEPS:-}" ]]; then
  STEPS=$(( RE > 100 ? 1000 : 600 ))
fi
CYCLES="${CYCLES:-4}"
CYCLES_R="${CYCLES_R:-6}"
TIME_ORDER="${TIME_ORDER:-2}"
WAVEFORM="${WAVEFORM:-}"
HARMONICS="${HARMONICS:-10}"
ETA_REF="${ETA_REF:-0.00345}"
LAMBDA_R="${LAMBDA_R:-4.83}"
ETAP_R="${ETAP_R:-0.050}"
if [[ -n "$WAVEFORM" ]]; then
  [[ -f "$WAVEFORM" ]] || { echo "Waveform file not found: $WAVEFORM"; exit 1; }
  WAVEFORM="$(cd "$(dirname "$WAVEFORM")" && pwd)/$(basename "$WAVEFORM")"
fi
PERTURB="${PERTURB:-0.01}"
AMP="${AMP:-0.5}"
PTT_EPS="${PTT_EPS:-0.1}"
OUTPUT_EVERY="${OUTPUT_EVERY:-40}"
MESH_LEVEL="${MESH_LEVEL:-inter}"
NEWTON_ITS="${NEWTON_ITS:-8}"
PSI_ORDER="${PSI_ORDER:-2}"
NEWTON_TOL="${NEWTON_TOL:-5e-4}"
PETSC_EXTRA="${PETSC_EXTRA:-}"
RETRY="${RETRY:-1}"
SOLVER="${SOLVER:-gmres}"
GMRES_OPTS="-ksp_type gmres -ksp_gmres_restart 200 -ksp_max_it 2000 -ksp_rtol 1e-8 -pc_type asm -pc_asm_overlap 1 -sub_pc_type ilu -sub_pc_factor_levels 1"
LULAG_OPTS="-ksp_type gmres -ksp_gmres_restart 100 -ksp_max_it 300 -ksp_rtol 1e-8 -pc_type lu -pc_factor_mat_solver_type mumps"
FORMULATION="${FORMULATION:-psi_rb}"
DYNAMIC="${DYNAMIC:-1}"
ADAPTIVE_DT="${ADAPTIVE_DT:-1}"
NEWTON_PREDICT="${NEWTON_PREDICT:-1}"
PC_REUSE="${PC_REUSE:-20}"
PC_REUSE_GROWTH="${PC_REUSE_GROWTH:-2}"
OUTLET_PSI="${OUTLET_PSI:-1}"
PSI_SKEW="${PSI_SKEW:-0}"
OUT="${OUT3D:-${OUT:-$HOME/Results/axi_jfm}/3D}/Re${RE}_Wo${WO}"   # never inside the axisymmetric folders
[[ -n "$WAVEFORM" ]] && OUT="${OUT}_$(basename "${WAVEFORM%.*}")"
MESH="$MESHDIR/${GEO}_pipe_${MESH_LEVEL}.mesh"

if [[ ! -x "$EXE" ]]; then
  echo "Executable not found: $EXE  (set EXE3D=...)"
  echo "  cmake --build <build> -j --target ArterialLesion3D_viscous_logarithmic_implicit"
  exit 1
fi
if [[ ! -f "$MESH" ]]; then
  echo "Mesh not found: $MESH"; exit 1
fi
if ! head -1 "$MESH" | grep -q MeshVersionFormatted; then
  echo "$MESH is not a MEDIT mesh (a Git LFS pointer?)"; exit 1
fi
mkdir -p "$OUT"
SUMMARY="$OUT/summary.csv"
[[ -f "$SUMMARY" ]] || echo "geometry,Re,Wo,fluid,Wi,steps_per_cycle,cycles,status,minutes,folder" > "$SUMMARY"

run_case() {   # fluid steps dir solver
  local fluid="$1" steps="$2" dir="$3" solver="$4" wi fargs cycles="$CYCLES"
  case "$fluid" in
    N)  wi=0.01; fargs=(-al_wi 0.01 -al_ptt_epsilon "$PTT_EPS") ;;
    D)  wi="$WI"; fargs=(-al_wi "$WI" -al_ptt_epsilon "$PTT_EPS") ;;
    CY) wi=-; fargs=(-al_fluid cy -al_eta_ref "$ETA_REF") ;;
    R)  wi=-; cycles="$CYCLES_R"
        fargs=(-al_fluid gsptt -al_eta_ref "$ETA_REF" -al_eta_p "$ETAP_R" -al_ptt_epsilon 0.2
               -al_lambda "$LAMBDA_R" -al_de 0 -al_wi 0 -al_inlet_developed_conformation 1) ;;
    *) echo "unknown fluid $fluid (3D: N, D, CY or R)"; return 2 ;;
  esac
  if [[ -f "$dir/DONE" ]]; then
    echo "[skip] $GEO $fluid Re=$RE Wo=$WO (DONE)"; return 0
  fi
  mkdir -p "$dir"; rm -f "$dir/FAILED"
  local cmd=("$MPIRUN" -n "$NP" "$EXE"
    -al_mesh "$MESH" -al_xdmf case -al_csv case.csv
    -al_re "$RE" -al_wo "$WO" -al_amplitude "$AMP" "${fargs[@]}" -al_time_order "$TIME_ORDER"
    -al_cycles "$cycles" -al_steps_per_cycle "$steps" -al_output_every "$OUTPUT_EVERY"
    -al_conformation_its "$NEWTON_ITS" -al_conformation_tol "$NEWTON_TOL" -al_psi_order "$PSI_ORDER"
    -al_inlet_perturbation "$PERTURB" -al_perturbation_cycles 1
    -al_formulation "$FORMULATION" -al_dynamic_subscales "$DYNAMIC"
    -al_adaptive_dt "$ADAPTIVE_DT" -al_newton_predict "$NEWTON_PREDICT"
    -al_pc_reuse_solves "$PC_REUSE" -al_pc_reuse_growth "$PC_REUSE_GROWTH"
    -al_outlet_conformation_inflow "$OUTLET_PSI" -al_psi_skew "$PSI_SKEW")
  [[ -n "$WAVEFORM" ]] && cmd+=(-al_waveform "$WAVEFORM" -al_harmonics "$HARMONICS")
  # shellcheck disable=SC2206
  [[ "$solver" == gmres ]] && cmd+=($GMRES_OPTS)
  # shellcheck disable=SC2206
  [[ "$solver" == lulag ]] && cmd+=($LULAG_OPTS)
  # shellcheck disable=SC2206
  [[ -n "$PETSC_EXTRA" ]] && cmd+=($PETSC_EXTRA)
  printf '%q ' "${cmd[@]}" > "$dir/params.txt"; echo >> "$dir/params.txt"

  echo "[run ] $(date '+%a %d %H:%M')  3D $GEO $fluid Wi=$wi Re=$RE Wo=$WO steps/cycle=$steps cycles=$cycles BDF$TIME_ORDER NP=$NP $solver -> $dir"
  local t0=$SECONDS
  ( cd "$dir" && "${cmd[@]}" > run.log 2>&1 )
  local status=$?
  local minutes=$(( (SECONDS - t0) / 60 ))
  if [[ $status -eq 0 ]]; then
    touch "$dir/DONE"
    echo "$GEO,$RE,$WO,$fluid,$wi,$steps,$cycles,DONE,$minutes,$dir" >> "$SUMMARY"
    echo "[done] $GEO $fluid in $minutes min"
  else
    echo "exit $status" > "$dir/FAILED"
    echo "$GEO,$RE,$WO,$fluid,$wi,$steps,$cycles,FAILED,$minutes,$dir" >> "$SUMMARY"
    echo "[FAIL] $GEO $fluid (exit $status, $minutes min): $(grep -m1 -E 'diverged|Fatal|Error' "$dir/run.log")"
  fi
  return $status
}

echo "3D Re=$RE Wo=$WO $GEO  fluids: $FLUIDS  inflow: ${WAVEFORM:-sinusoid A=$AMP}  BDF$TIME_ORDER  steps/cycle: $STEPS  cycles: $CYCLES (R $CYCLES_R)  seed: $PERTURB  mesh: $MESH_LEVEL  NP=$NP  ->  $OUT"
for fluid in $FLUIDS; do
  if [[ "$fluid" == "D" ]]; then sub="Wi${WI}"; else sub="$fluid"; fi
  dir="$OUT/$GEO/$sub"
  run_case "$fluid" "$STEPS" "$dir" "$SOLVER"
  st=$?
  if [[ $st -ne 0 && $st -ne 2 && "$RETRY" == "1" && "$SOLVER" == gmres ]]; then
    echo "[retry] $GEO $fluid with the direct solver (mumps)"
    run_case "$fluid" "$STEPS" "${dir}_lu" mumps
    st=$?
  fi
  if [[ $st -ne 0 && $st -ne 2 && "$RETRY" == "1" ]]; then
    echo "[retry] $GEO $fluid with $(( 2 * STEPS )) steps/cycle (mumps)"
    run_case "$fluid" "$(( 2 * STEPS ))" "${dir}_dt2" mumps
  fi
done
echo "3D Re=$RE $GEO finished $(date). Summary: $SUMMARY"
