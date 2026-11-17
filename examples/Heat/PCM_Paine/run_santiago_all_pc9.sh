#!/usr/bin/env bash
# Both design days of Santiago on pc9 with at most MAXNP (50) MPI ranks:
#   1) summer: hottest January day of the PVGIS TMY  -> OUT_SUMMER (pcm_santiago)
#   2) winter: coldest July day (lowest daily mean T) -> OUT_WINTER (pcm_santiago_frio)
# Each season runs its 10 thicknesses at once (one group), with ranks sized to
# the mesh: 4 for 5-20 mm, 5 for 25-35 mm, 6 for 40-50 mm (49 in total).
# Finished cases (DONE) are skipped, so a completed summer sweep is not rerun.
#
#   mkdir -p ~/runs_santiago && cd ~/runs_santiago
#   nohup bash ~/rodin/examples/Heat/PCM_Paine/run_santiago_all_pc9.sh > all.log 2>&1 &
#
# Room set-point: TIN_SUMMER = 24 C (cooling), TIN_WINTER = 20 C (heating).
set -euo pipefail
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export MAXNP="${MAXNP:-50}"
export TANDAS="${TANDAS:-5 10 15 20 25 30 35 40 45 50}"
export NP_MAP="${NP_MAP:-5:4 10:4 15:4 20:4 25:5 30:5 35:5 40:6 45:6 50:6}"
OUT_SUMMER="${OUT_SUMMER:-$PWD/pcm_santiago}"
OUT_WINTER="${OUT_WINTER:-$PWD/pcm_santiago_frio}"
SEASONS="${SEASONS:-summer winter}"

for s in $SEASONS; do
  case "$s" in
    summer) out="$OUT_SUMMER"; month=1; dia=caluroso; tin="${TIN_SUMMER:-24}" ;;
    winter) out="$OUT_WINTER"; month=7; dia=frio;     tin="${TIN_WINTER:-20}" ;;
    *) echo "unknown season $s" >&2; exit 1 ;;
  esac
  # reuse the TMY downloaded by an earlier run (no internet needed afterwards)
  tmy="${TMY_JSON:-}"
  for f in "$OUT_SUMMER/tmy_santiago.json" "$OUT_WINTER/tmy_santiago.json"; do
    [ -z "$tmy" ] && [ -s "$f" ] && tmy="$f"
  done
  echo "[season] $s: month $month, $dia day, T_in = $tin C -> $out ($(date +%d/%m\ %H:%M))"
  OUT="$out" MONTH="$month" DIA="$dia" TIN="$tin" TMY_JSON="$tmy" bash "$HERE/run_santiago_pc9.sh" \
    || echo "[warn] season $s ended with errors, see above" >&2
done
