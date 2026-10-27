#!/usr/bin/env bash
# Thickness sweep of NaturalConvectionPCM: RT42 roof layer, typical clear
# January day in Paine (RM), north-facing roof at ANGLE deg, ventilated attic
# underneath (entretecho: h_in 6 W/m^2K, ceiling U_c, ventilation ACH).
#
#   cd ~/runs && bash ~/rodin/examples/Heat/PCM_Paine/run_pcm_sweep.sh
#   NP=8 NJOBS=10 bash run_pcm_sweep.sh            # pc9: all cases at once
#   THICKNESSES="5 25 50" REFINE=2 bash run_pcm_sweep.sh   # mesh check
#
# Every case runs in OUT/e_XXmm (CSV + XDMF/HDF5 of u, p, T, Y every hour);
# finished cases (file DONE) are skipped, so the script can be relaunched.
# Then: bash postprocess_pcm.sh OUT
#
# Mesh (structured, triangles): h_x = HX (5 mm: 400 cells over 2 m) and
# n_y = max(NYMIN, e/HYMAX) across the layer, so h_y <= 1.25 mm: >= 20 cells
# across every layer and 40 at 5 cm. At 5 cm the layer cooled from above at
# night reaches Ra = g beta dT e^3/(nu alpha) ~ 1e6 (dT 5 K, Pr ~ 240): cells
# of width ~2e = 10 cm (20 h_x) and thermal boundary layers of ~e/(2 Nu) ~
# 4 mm (3+ h_y). REFINE=2 halves both sizes for a mesh-independence check.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"

BUILD="${BUILD:-$REPO/build}"
BIN="${BIN:-$BUILD/examples/Heat/NaturalConvectionPCM}"
OUT="${OUT:-$PWD/pcm_paine_runs}"
MPIRUN="${MPIRUN:-mpirun}"     # may carry options, e.g. "mpirun --oversubscribe"
NP="${NP:-4}"                 # MPI ranks per case
NJOBS="${NJOBS:-1}"           # cases running at the same time
THICKNESSES="${THICKNESSES:-5 10 15 20 25 30 35 40 45 50}"   # mm
ANGLE="${ANGLE:-25}"          # deg, roof slope (typical gable roof)
FACING="${FACING:-north}"     # roof orientation (north: most sun in Chile)
NDAYS="${NDAYS:-3}"           # periodic days; the last one is analysed
DT="${DT:-30}"                # s
MUSHY="${MUSHY:-1}"           # K, T_l - T_s (solver default 1 K)
RETRIES="${RETRIES:-2}"       # restarts of a failed case (see run_case)
LENGTH="${LENGTH:-2.0}"       # m, along the slope
HX="${HX:-0.005}"             # m
HYMAX="${HYMAX:-0.00125}"     # m
NYMIN="${NYMIN:-20}"
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
WEATHER="$OUT/paine_enero_${FACING}_${ANGLE}deg.csv"
python3 "$HERE/paine_weather.py" --tilt "$ANGLE" --facing "$FACING" -o "$WEATHER"

TEND=$(awk -v d="$NDAYS" 'BEGIN { printf "%d", d * 86400 }')
EVERY=$(awk -v dt="$DT" 'BEGIN { printf "%d", 3600 / dt + 0.5 }')
NX=$(awk -v L="$LENGTH" -v h="$HX" -v r="$REFINE" 'BEGIN { printf "%d", r * int(L / h + 0.5) }')

{
  echo "date: $(date)"
  echo "git: $(git -C "$REPO" rev-parse --short HEAD 2>/dev/null || echo '?')"
  echo "THICKNESSES=$THICKNESSES ANGLE=$ANGLE FACING=$FACING NDAYS=$NDAYS DT=$DT"
  echo "LENGTH=$LENGTH HX=$HX HYMAX=$HYMAX NYMIN=$NYMIN REFINE=$REFINE NP=$NP"
  echo "TIN=$TIN UCEIL=$UCEIL ACH=$ACH ATTICH=$ATTICH EXTRA=$EXTRA"
} > "$OUT/sweep.txt"

run_case() {
  local mm="$1"
  local e ny dir
  e=$(awk -v mm="$mm" 'BEGIN { printf "%.6f", mm / 1000 }')
  ny=$(awk -v e="$e" -v h="$HYMAX" -v n="$NYMIN" -v r="$REFINE" \
       'BEGIN { m = e / h; m = (m == int(m)) ? m : int(m) + 1; if (m < n) m = n; printf "%d", r * m }')
  dir="$OUT/$(printf 'e_%02dmm' "$mm")"
  if [ -f "$dir/DONE" ]; then
    echo "[skip] $dir"
    return 0
  fi
  mkdir -p "$dir"
  # If Newton fails (the run aborts), the case restarts with a 2 K mushy range
  # (RT42 melts over 38-43 C, so still physical), then also with dt/2. A 1 K
  # range under ~2 cells of melt at the onset of melting can make Newton
  # diverge in thick layers (seen at 50 mm); a smaller dt alone does not help.
  # Retries also allow 100 Newton iterations (default 30).
  local dt="$DT" every="$EVERY" mushy="$MUSHY" try=0 retryOpts=""
  while :; do
    echo "thickness_mm=$mm nx=$NX ny=$ny angle=$ANGLE dt=$dt mushy=$mushy ndays=$NDAYS" > "$dir/case.txt"
    echo "[run ] e = $mm mm  mesh $NX x $ny  dt = $dt s  mushy = $mushy K  ($(date +%H:%M))"
    # shellcheck disable=SC2086
    if (cd "$dir" && $MPIRUN -n "$NP" "$BIN" \
          -pcm_thickness "$e" -pcm_length "$LENGTH" -pcm_nx "$NX" -pcm_ny "$ny" \
          -pcm_angle "$ANGLE" -pcm_weather "$WEATHER" -pcm_dt "$dt" -pcm_tend "$TEND" -pcm_mushy "$mushy" \
          -pcm_output_every "$every" -pcm_print_every "$every" \
          -pcm_attic 1 -pcm_tin "$TIN" -pcm_uceil "$UCEIL" -pcm_attic_ach "$ACH" \
          -pcm_attic_height "$ATTICH" $retryOpts $EXTRA > run.log 2>&1); then
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
    retryOpts="-snes_max_it 100"
    if awk -v m="$mushy" 'BEGIN { exit !(m < 2) }'; then
      mushy=2
    else
      dt=$(awk -v d="$dt" 'BEGIN { printf "%g", d / 2 }')
      every=$((every * 2))
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
