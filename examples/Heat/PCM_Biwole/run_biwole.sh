#!/usr/bin/env bash
# Validation of NaturalConvectionPCM against the experiment without fins of
# Biwole, Eclache & Kuznik, Energy Build. 62 (2013) 59-67, Sec. 2.3:
# Plexiglas cavity 167 mm high x 30 mm wide, PCM of their Table 2, a fixed
# temperature on each vertical wall, top and bottom adiabatic, no-slip.
#
#   melt  : hot wall 40 C, cold wall 20 C, PCM solid at 20 C; front every
#           20 min up to 160 min (their Fig. 5).
#   steady: hot wall 38 C, cold wall 28 C, PCM fully liquid; velocity
#           profiles at z = 41, 82, 123 mm and at s = 15 mm (their Figs. 6-7).
#
# The cavity is the solver's rectangle at 90 deg: x (length) is the height,
# y (thickness) the width; g = -|g| e_x. Top face (y = W) = hot wall, held by a
# Robin penalty h_pen (T - T_hot) without sun nor sky radiation
# (-pcm_top_fixed); bottom face (y = 0) = cold wall, the room condition with
# h_in = h_pen (-pcm_attic 0 -pcm_hin -pcm_tin); ends (x = 0, H) adiabatic.
#
# Not given in the paper (assumed here): initial temperature of the melt case
# (20 C, the cold wall), mushy range (1 K = 2 x their Delta T of Sec. 3),
# mu = rho nu = 890 x 1.1e-5 Pa s.
#
# Mesh: h_x = 0.5 mm along the height; across, cells of ~0.1 mm at the hot
# wall growing by 1.03 (78 cells, ~1 mm at the cold wall). The first liquid
# layer at the hot wall must span several cells: with a uniform 0.5 mm, Newton
# fails at t ~ 16 s, when the melt starts to flow (front speed
# q/(rho L_F) ~ 0.03 mm/s at q ~ 5 kW/m^2).
#
#   cd ~/runs && NP=4 bash ~/rodin-heat/examples/Heat/PCM_Biwole/run_biwole.sh
#   python3 ~/rodin-heat/examples/Heat/PCM_Biwole/postprocess_biwole.py biwole_runs
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
BIN="${BIN:-$REPO/build/examples/Heat/NaturalConvectionPCM}"
OUT="${OUT:-$PWD/biwole_runs}"
MPIRUN="${MPIRUN:-mpirun}"
NP="${NP:-4}"
CASES="${CASES:-melt steady}"
DT="${DT:-2}"                 # s
NX="${NX:-334}"               # cells along the height (0.5 mm)
NY="${NY:-78}"                # cells across the width
GRADE="${GRADE:-1.03}"        # growth from the hot wall
MUSHY="${MUSHY:-1}"           # K
HPEN="${HPEN:-1e6}"           # W/m^2K, wall penalty
TEND_MELT="${TEND_MELT:-9600}"     # s, 160 min
TEND_STEADY="${TEND_STEADY:-7200}" # s, 2 h to the steady state (check maxU in the CSV)
EXTRA="${EXTRA:-}"

[ -x "$BIN" ] || { echo "Missing $BIN: cmake --build $REPO/build -j --target NaturalConvectionPCM" >&2; exit 1; }

# PCM of Biwole et al. (2013), Table 2.
PCM="-pcm_rho_s 890 -pcm_rho_l 890 -pcm_cp_s 2110 -pcm_cp_l 2110 -pcm_k_s 0.182 -pcm_k_l 0.182
     -pcm_latent 1.81e5 -pcm_mu 9.79e-3 -pcm_beta 3.09e-3 -pcm_tm 25 -pcm_cmushy 1e5 -pcm_mushy $MUSHY"
GEO="-pcm_length 0.167 -pcm_thickness 0.030 -pcm_angle 90 -pcm_al 0
     -pcm_nx $NX -pcm_ny $NY -pcm_ygrade $GRADE"
WALLS="-pcm_emissivity 0 -pcm_absorptivity 0 -pcm_attic 0 -pcm_hin $HPEN -pcm_hpen $HPEN"

run() {
  local name="$1" thot="$2" tcold="$3" tinit="$4" tend="$5" every="$6"
  local dir="$OUT/$name"
  if [ -f "$dir/DONE" ]; then echo "[skip] $dir"; return 0; fi
  mkdir -p "$dir"
  echo "[run ] $name: T_hot = $thot C, T_cold = $tcold C, T(0) = $tinit C, dt = $DT s ($(date +%H:%M))"
  # shellcheck disable=SC2086
  (cd "$dir" && $MPIRUN -n "$NP" "$BIN" $GEO $PCM $WALLS \
      -pcm_top_fixed "$thot" -pcm_tin "$tcold" -pcm_tinit "$tinit" \
      -pcm_dt "$DT" -pcm_tend "$tend" -pcm_output_every "$every" -pcm_print_every "$every" \
      -snes_max_it 50 $EXTRA > run.log 2>&1) && touch "$dir/DONE" \
    && echo "[done] $name ($(date +%H:%M))" || echo "[FAIL] $name: see $dir/run.log" >&2
}

every=$(awk -v dt="$DT" 'BEGIN { printf "%d", 1200 / dt + 0.5 }')   # 20 min
for c in $CASES; do
  case "$c" in
    melt)   run melt 40 20 20 "$TEND_MELT" "$every" ;;
    steady) run steady 38 28 33 "$TEND_STEADY" "$every" ;;
  esac
done
echo "Next: python3 $HERE/postprocess_biwole.py $OUT"
