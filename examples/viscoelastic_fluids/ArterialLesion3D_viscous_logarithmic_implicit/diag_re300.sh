#!/usr/bin/env bash
# =============================================================================
#  diag_re300.sh -- the failure of fluid R at Re 300 with the CCA waveform,
#  on the "layersc" mesh (37% of the cells of inter, the same h in the
#  throat wall and shear layers), so that it is reached in a few hours.
#
#     cd ~/rodin && nohup bash examples/viscoelastic_fluids/ArterialLesion3D_viscous_logarithmic_implicit/diag_re300.sh [geo] > diag_re300.log 2>&1 &
#
#  Runs, at the same time (NP ranks each, default 32):
#     base    psi_rb as in the campaign
#     nodiv   psi_rb without the div sigma subscale u3 (-al_stress_div_scale 0)
#     psi     the former formulation
#  2 cycles (S75 failed at t = 0.28 s in the ramp cycle, A200 at 0.59 s in
#  the second), XDMF every 20 steps. Each
#  run stops at its first fatal error; then locate_blowup.py prints where
#  |sigma|, |u| and psi were extreme in the last outputs:
#     $OUT/<variant>/blowup.txt
#  Environment: NP=32  OUT=$HOME/Results/diag_re300  MESH_LEVEL=layersc
#               VARIANTS="base nodiv psi"  STEPS=1600  CYCLES_R=2
#  The mesh: resources/examples/viscoelastic_fluids/<geo>_pipe_layersc.mesh,
#  made if missing with
#     python3 make_lesion_mesh_3d.py <geo> 0.7 <mesh> -setnumber layers 1
# =============================================================================
set -u
GEO="${1:-S75}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
NP="${NP:-32}"
OUT="${OUT:-$HOME/Results/diag_re300}"
MESH_LEVEL="${MESH_LEVEL:-layersc}"
VARIANTS="${VARIANTS:-base nodiv psi}"
STEPS="${STEPS:-1600}"
CYCLES_R="${CYCLES_R:-2}"
MESHDIR="$REPO/resources/examples/viscoelastic_fluids"
PY="${PYTHON:-$(command -v python3)}"
MESH="$MESHDIR/${GEO}_pipe_${MESH_LEVEL}.mesh"
if [[ ! -f "$MESH" ]]; then
  case "$MESH_LEVEL" in
    layersc) rf=0.7 ;; layers) rf=1 ;; *) echo "mesh $MESH not found"; exit 1 ;;
  esac
  echo "== making $MESH ..."
  "$PY" "$HERE/make_lesion_mesh_3d.py" "$GEO" "$rf" "$MESH" -setnumber layers 1 || exit 1
fi
mkdir -p "$OUT"
for V in $VARIANTS; do
  case "$V" in
    base)  F=psi_rb; X="" ;;
    nodiv) F=psi_rb; X="-al_stress_div_scale 0" ;;
    psi)   F=psi;    X="" ;;
    *) echo "unknown variant $V"; continue ;;
  esac
  (
    env FLUIDS=R NP="$NP" TIME_ORDER=2 MESH_LEVEL="$MESH_LEVEL" OUTPUT_EVERY=20 \
      WAVEFORM="$MESHDIR/waveform_cca.dat" HARMONICS=10 STEPS="$STEPS" CYCLES_R="$CYCLES_R" SOLVER=lulag \
      FORMULATION="$F" PETSC_EXTRA="$X" OUT3D="$OUT/$V" \
      bash "$HERE/run_3d.sh" 300 "$GEO" > "$OUT/$V.log" 2>&1
    d=$(find "$OUT/$V" -name case.xdmf -printf '%h\n' | head -1)
    [[ -n "$d" ]] && "$PY" "$HERE/locate_blowup.py" "$d" --last 6 > "$OUT/$V/blowup.txt" 2>&1
  ) < /dev/null &
  echo "   started $V ($F $X) -> $OUT/$V.log"
  sleep 2
done
echo "== progress: tail -2 \$(find $OUT -name case.csv)    results: $OUT/<variant>/blowup.txt"
wait
echo "== all variants finished"
