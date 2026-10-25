#!/usr/bin/env bash
# =============================================================================
#  run_diag_R.sh -- diagnostic case of the R fluid (gsPTT) on pc9: S75 carotid
#  (D 6 mm, Re 639, Wo 4, CCA waveform), BDF2, 1 cycle, mesh layersc (121k).
#  Runs 9, 12, 13 (eta_p 0.0104) blew up at t = 0.584 s; here eta_p = 0.0070,
#  the largest with an increasing solvent stress. XDMF dumps when max|sigma|
#  > DUMP Pa, then locate_blowup.py if it fails.
#
#     cd ~/rodin && bash examples/viscoelastic_fluids/ArterialLesion3D_viscous_logarithmic_implicit/run_diag_R.sh
#
#  Environment: ETAP_R=0.0070 LAMBDA_R=0.1 EPS_R=0.2 MESH_LEVEL=layersc NP=32 BUILD=build-gcc12
#               OUT3D=$HOME/Results/diag_etap<ETAP_R> DUMP=100 SKIP_BUILD=0
# =============================================================================
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
cd "$REPO"
BUILD="${BUILD:-build-gcc12}"
ETAP_R="${ETAP_R:-0.0070}"
LAMBDA_R="${LAMBDA_R:-0.1}"
EPS_R="${EPS_R:-0.2}"
MESH_LEVEL="${MESH_LEVEL:-layersc}"
NP="${NP:-32}"
DUMP="${DUMP:-100}"
OUT3D="${OUT3D:-$HOME/Results/diag_l${LAMBDA_R}_etap${ETAP_R}_eps${EPS_R}}"
EXE3D="$REPO/$BUILD/examples/viscoelastic_fluids/ArterialLesion3D_viscous_logarithmic_implicit/ArterialLesion3D_viscous_logarithmic_implicit"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

if [[ "${SKIP_BUILD:-0}" != 1 ]]; then
  echo "== building the 3D driver in $BUILD ..."
  cmake --build "$REPO/$BUILD" -j 8 --target ArterialLesion3D_viscous_logarithmic_implicit \
    > "$REPO/$BUILD/diag_build.log" 2>&1 \
    || { echo "BUILD FAILED, see $REPO/$BUILD/diag_build.log"; tail -20 "$REPO/$BUILD/diag_build.log"; exit 1; }
fi
grep -aq -- "-al_dump_stress" "$EXE3D" || { echo "$EXE3D has no -al_dump_stress: git pull and rebuild"; exit 1; }

echo "== R: eta_p=$ETAP_R lambda=$LAMBDA_R eps=$EPS_R, mesh $MESH_LEVEL, $NP ranks -> $OUT3D"
env FLUIDS=R NP="$NP" TIME_ORDER=2 MESH_LEVEL="$MESH_LEVEL" WO=4 STEPS=1600 CYCLES=1 CYCLES_R=1 \
  SOLVER=lulag HARMONICS=10 LAMBDA_R="$LAMBDA_R" ETAP_R="$ETAP_R" EPS_R="$EPS_R" OUTPUT_EVERY=80 RETRY=0 \
  WAVEFORM="$REPO/resources/examples/viscoelastic_fluids/waveform_cca.dat" \
  PETSC_EXTRA="-al_diameter 0.006 -al_dump_stress $DUMP" \
  EXE3D="$EXE3D" OUT3D="$OUT3D" \
  bash "$HERE/run_3d.sh" 639 S75

D="$OUT3D/Re639_Wo4_waveform_cca/S75/R"
echo "== max|sigma| around the former blow-up (t = 0.55-0.65 s):"
awk -F, 'NR==1{for(i=1;i<=NF;i++){if($i=="t")a=i; if($i=="maxSigma")b=i; if($i=="maxU")c=i}} NR>1 && $a>0.55 && $a<0.65 && NR%20==0 {printf "   t=%.4f  max|sigma|=%.4g Pa  max|u|=%.3f m/s\n",$a,$b,$c}' "$D/case.csv"
if [[ -f "$D/FAILED" ]]; then
  echo "== FAILED: where (locate_blowup.py)"
  python3 "$HERE/locate_blowup.py" "$D" --last 8 | tee "$D/blowup.txt"
else
  echo "== DONE: passed the cycle"
fi
