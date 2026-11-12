#!/usr/bin/env bash
# Thickness sweep of NaturalConvectionPCM: RT42 roof layer between two
# aluminium plates of AL m (sandwich), typical clear January day in SITE
# (paine or calama), north-facing roof at ANGLE deg. Outside: solar gain,
# convection h_a(wind) and long-wave radiation to the sky. Inside: convection
# only (h_in = HIN) to a ventilated attic (entretecho: ceiling U_c, ACH).
#
#   cd ~/runs && bash ~/rodin/examples/Heat/PCM_Paine/run_pcm_sweep.sh
#   NP=8 NJOBS=10 bash run_pcm_sweep.sh            # pc9: all cases at once
#   THICKNESSES="5 20" REFINE=2 bash run_pcm_sweep.sh   # mesh check
#   AL=0 HX=0.005 GRADE=1.04 NYMIN=20 bash run_pcm_sweep.sh # PCM alone
#
# Every case runs in OUT/e_XXmm (CSV + XDMF/HDF5 of u, p, T, Y every hour);
# finished cases (file DONE) are skipped, so the script can be relaunched.
# Then: bash postprocess_pcm.sh OUT
#
# Mesh (structured, triangles): h_x = HX (10 mm: 200 cells over 2 m); NAL
# cells across each plate; across the PCM, cells of HYTOP = 0.25 mm at its top
# face growing by GRADE = 1.08 downwards,
# n_y = max(NYMIN, ln(1 + e (GRADE - 1)/HYTOP)/ln GRADE): 13 cells at 5 mm,
# 26 at 20 mm (0.25-1.7 mm). Melting and
# solidification start at the top face; the first liquid layer must span
# several cells: with a uniform h_y of 1.25 mm (or 0.5 mm) at 50 mm, Newton
# iterations grow from ~5 to 30 within minutes of the onset of melting and the
# step diverges, while 0.25 mm keeps ~9 iterations (the 5 mm run that worked
# had 0.25 mm). At 5 cm the layer cooled from above at night reaches
# Ra ~ 1e6 (dT 5 K, Pr ~ 240): cells of width ~2e = 10 cm (20 h_x) and thermal
# boundary layers of ~4 mm. REFINE=2 doubles n_x, n_y (GRADE -> sqrt(GRADE)).
#
# GRADING=both (-pcm_ygrade_both): HYTOP at BOTH faces of the PCM, growing
# geometrically toward the mid-plane by r <= GRADEMAX up to HYMAX: n_y = 2m
# with the smallest m such that HYTOP (r^m - 1)/(r - 1) = e/2 has
# r <= GRADEMAX and HYTOP r^(m-1) <= HYMAX. Lesson of the Biwole validation:
# the melt front, which crosses the whole layer, and the liquid against the
# second wall (no-slip, thin crust) need fine cells in the interior and at
# both faces; with one-sided grading (1.25 mm cells at the far wall) T
# undershot by 10 K and Newton failed down to dt = 0.03 s.
# NP_MAP="5:25 10:25 ..." gives the MPI ranks of each thickness (others: NP).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"

BUILD="${BUILD:-$REPO/build}"
BIN="${BIN:-$BUILD/examples/Heat/NaturalConvectionPCM}"
OUT="${OUT:-$PWD/pcm_runs}"
MPIRUN="${MPIRUN:-mpirun}"     # may carry options, e.g. "mpirun --oversubscribe"
NP="${NP:-4}"                 # MPI ranks per case
NJOBS="${NJOBS:-1}"           # cases running at the same time
THICKNESSES="${THICKNESSES:-5 10 15 20}"   # mm of PCM
SITE="${SITE:-calama}"        # weather: paine or calama (paine_weather.py --site), or santiago
                              # (pvgis_weather.py: PVGIS TMY of Quinta Normal, January)
TMY_JSON="${TMY_JSON:-}"      # santiago: PVGIS TMY JSON already downloaded (else downloaded)
DIA="${DIA:-promedio}"        # santiago: promedio | caluroso | despejado | day of the month
AL="${AL:-0.01}"              # m, aluminium plate below and above the PCM (0: none)
NAL="${NAL:-4}"               # cells across each plate
HIN="${HIN:-2}"               # W/m^2K, convection only on the attic side
ANGLE="${ANGLE:-25}"          # deg, roof slope (typical gable roof)
FACING="${FACING:-north}"     # roof orientation (north: most sun in Chile)
NDAYS="${NDAYS:-2}"           # periodic days; the last one is analysed
DT="${DT:-30}"                # s
MUSHY="${MUSHY:-2}"           # K, T_l - T_s (solver default 1 K; RT42 melts over 38-43 C)
MAXIT="${MAXIT:-100}"         # Newton iterations per step
RETRIES="${RETRIES:-1}"       # restarts of a failed case (see run_case)
LENGTH="${LENGTH:-2.0}"       # m, along the slope
HX="${HX:-0.01}"              # m
HYTOP="${HYTOP:-0.00025}"     # m, cell height at the top face
GRADE="${GRADE:-1.08}"        # height ratio of consecutive cells, top -> bottom
NYMIN="${NYMIN:-10}"
GRADING="${GRADING:-top}"     # top (HYTOP at the top face, GRADE) or both (see above)
HYMAX="${HYMAX:-0.001}"       # m, largest cell across the PCM (GRADING=both)
GRADEMAX="${GRADEMAX:-1.06}"  # largest growth ratio (GRADING=both)
NP_MAP="${NP_MAP:-}"          # per-thickness MPI ranks, "mm:np mm:np ..."
REFINE="${REFINE:-1}"
TIN="${TIN:-24}"              # C, room below the ceiling
UCEIL="${UCEIL:-0.38}"        # W/m^2K, insulated ceiling (OGUC 4.1.10)
ACH="${ACH:-5}"               # 1/h, attic ventilation
ATTICH="${ATTICH:-0.9}"       # m, attic volume per unit roof area
EXTRA="${EXTRA:-}"            # more options for the solver

if [ ! -x "$BIN" ]; then
  echo "Missing $BIN: cmake --build $BUILD -j --target NaturalConvectionPCM" >&2
  exit 1
fi
mkdir -p "$OUT"
WEATHER="$OUT/${SITE}_enero_${FACING}_${ANGLE}deg.csv"
if [ "$SITE" = santiago ]; then
  WEATHER="$OUT/santiago_enero_${DIA}_${FACING}_${ANGLE}deg.csv"
  if [ -n "$TMY_JSON" ]; then
    python3 "$HERE/pvgis_weather.py" --input "$TMY_JSON" --dia "$DIA" --tilt "$ANGLE" --facing "$FACING" -o "$WEATHER"
  else
    python3 "$HERE/pvgis_weather.py" --save "$OUT/tmy_santiago.json" --dia "$DIA" --tilt "$ANGLE" --facing "$FACING" -o "$WEATHER"
  fi
else
  python3 "$HERE/paine_weather.py" --site "$SITE" --tilt "$ANGLE" --facing "$FACING" -o "$WEATHER"
fi
[ -s "$WEATHER" ] || { echo "No weather file $WEATHER" >&2; exit 1; }

TEND=$(awk -v d="$NDAYS" 'BEGIN { printf "%d", d * 86400 }')
EVERY=$(awk -v dt="$DT" 'BEGIN { printf "%d", 3600 / dt + 0.5 }')
NX=$(awk -v L="$LENGTH" -v h="$HX" -v r="$REFINE" 'BEGIN { printf "%d", r * int(L / h + 0.5) }')

{
  echo "date: $(date)"
  echo "git: $(git -C "$REPO" rev-parse --short HEAD 2>/dev/null || echo '?')"
  echo "WEATHER=$WEATHER ($(head -n 1 "$WEATHER"))"
  echo "THICKNESSES=$THICKNESSES SITE=$SITE AL=$AL NAL=$NAL HIN=$HIN ANGLE=$ANGLE FACING=$FACING NDAYS=$NDAYS DT=$DT"
  echo "LENGTH=$LENGTH HX=$HX HYTOP=$HYTOP GRADE=$GRADE NYMIN=$NYMIN REFINE=$REFINE NP=$NP"
  echo "GRADING=$GRADING HYMAX=$HYMAX GRADEMAX=$GRADEMAX NP_MAP=$NP_MAP MUSHY=$MUSHY MPIRUN=$MPIRUN"
  echo "TIN=$TIN UCEIL=$UCEIL ACH=$ACH ATTICH=$ATTICH EXTRA=$EXTRA"
} > "$OUT/sweep.txt"

run_case() {
  local mm="$1"
  local e ny grade dir np both=0
  e=$(awk -v mm="$mm" 'BEGIN { printf "%.6f", mm / 1000 }')
  np=$(echo "$NP_MAP" | tr ' ' '\n' | awk -F: -v mm="$mm" '$1 == mm { print $2 }' | tail -n 1)
  np="${np:-$NP}"
  if [ "$GRADING" = both ]; then
    both=1
    # smallest half count m whose ratio r (HYTOP (r^m - 1)/(r - 1) = e/2) keeps
    # r <= GRADEMAX and the mid-plane cell HYTOP r^(m-1) <= HYMAX
    read -r ny grade < <(awk -v e="$e" -v h="$HYTOP" -v hm="$HYMAX" -v gm="$GRADEMAX" -v rf="$REFINE" '
      BEGIN {
        half = e / 2
        for (m = 1; m < 10000; ++m) {
          if (h * m >= half) r = 1
          else {
            lo = 1 + 1e-12; hi = 3
            for (k = 0; k < 200; ++k) { r = (lo + hi) / 2; if (h * (r ^ m - 1) / (r - 1) < half) lo = r; else hi = r }
            r = (lo + hi) / 2
          }
          if (r <= gm && h * r ^ (m - 1) <= hm * (1 + 1e-9)) break
        }
        printf "%d %.8g\n", rf * 2 * m, r ^ (1 / rf)
      }')
  else
    ny=$(awk -v e="$e" -v h="$HYTOP" -v g="$GRADE" -v n="$NYMIN" -v r="$REFINE" \
         'BEGIN { m = (g > 1) ? log(1 + e * (g - 1) / h) / log(g) : e / h;
                  m = (m == int(m)) ? m : int(m) + 1; if (m < n) m = n; printf "%d", r * m }')
    grade=$(awk -v g="$GRADE" -v r="$REFINE" 'BEGIN { printf "%.8g", g ^ (1 / r) }')
  fi
  dir="$OUT/$(printf 'e_%02dmm' "$mm")"
  if [ -f "$dir/DONE" ]; then
    echo "[skip] $dir"
    return 0
  fi
  mkdir -p "$dir"
  # If Newton fails (the run aborts), the case restarts with a 4 K mushy range
  # (RT42 melts over 38-43 C, so still physical). A smaller dt does not help:
  # in the tests it failed earlier, in the same state (melt-layer flow
  # ~3e-4 m/s). With 2 K and 100 iterations, 5-20 mm passed in the tests.
  local dt="$DT" every="$EVERY" mushy="$MUSHY" try=0 retryOpts=""
  while :; do
    echo "thickness_mm=$mm al_mm=$(awk -v a="$AL" 'BEGIN { printf "%g", a * 1000 }') nx=$NX ny=$ny nal=$NAL grade=$grade grading=$GRADING np=$np angle=$ANGLE dt=$dt mushy=$mushy ndays=$NDAYS" > "$dir/case.txt"
    echo "[run ] e = $mm mm  mesh $NX x $ny (grade $grade, $GRADING)  np = $np  dt = $dt s  mushy = $mushy K  ($(date +%H:%M))"
    # shellcheck disable=SC2086
    if (cd "$dir" && $MPIRUN -n "$np" "$BIN" \
          -pcm_thickness "$e" -pcm_length "$LENGTH" -pcm_nx "$NX" -pcm_ny "$ny" -pcm_ygrade "$grade" -pcm_ygrade_both "$both" \
          -pcm_angle "$ANGLE" -pcm_weather "$WEATHER" -pcm_dt "$dt" -pcm_tend "$TEND" -pcm_mushy "$mushy" \
          -pcm_output_every "$every" -pcm_print_every "$every" \
          -pcm_al "$AL" -pcm_al_n "$NAL" -pcm_hin "$HIN" \
          -pcm_attic 1 -pcm_tin "$TIN" -pcm_uceil "$UCEIL" -pcm_attic_ach "$ACH" \
          -pcm_attic_height "$ATTICH" -snes_max_it "$MAXIT" $retryOpts $EXTRA > run.log 2>&1); then
      touch "$dir/DONE"
      echo "[done] e = $mm mm  ($(date +%H:%M))"
      return 0
    fi
    try=$((try + 1))
    local failed="$dir/run_failed_dt${dt}_mushy${mushy}.log"
    cp "$dir/run.log" "$failed"
    if [ "$try" -gt "$RETRIES" ]; then
      echo "[FAIL] e = $mm mm: see $failed" >&2
      return 0
    fi
    echo "[retry] e = $mm mm: Newton failed with dt = $dt s, mushy = $mushy K ($(grep -m1 'Fatal' "$dir/run.log" || true))"
    if awk -v m="$mushy" 'BEGIN { exit !(m < 2) }'; then
      mushy=2
    else
      mushy=4
      retryOpts="-pcm_round 1"
    fi
  done
}

for mm in $THICKNESSES; do
  while [ "$(jobs -rp | wc -l | tr -d ' ')" -ge "$NJOBS" ]; do
    sleep 5
  done
  run_case "$mm" &
done
wait
echo "Results in $OUT. Next: bash $HERE/postprocess_pcm.sh $OUT"
