#!/usr/bin/env bash
# seat_reports.sh — seat report and force report of the seat study, for a
# chosen set of designs and openings, including runs still in progress
# (marked in yellow, provisional values).
#
#   bash examples/Valve/seat_reports.sh [RUNS] [OUTDIR]
#     RUNS    folder with the case folders (default: runs_seat)
#     OUTDIR  where the reports go (default: seat_reports)
#   env: DESIGNS  seat designs, the names of make_meshes.py --seat
#                 (default: "cad fillet rampup rampdn lip"; also chamfer, upramp)
#        OPENINGS openings in % (default: "5 10 25")
#
# Writes OUTDIR/runs/ (links to the selected case folders, nothing copied),
# OUTDIR/seat_report.{pdf,csv} (flow in the seat) and
# OUTDIR/forces_report.{pdf,_summary.csv} (forces, losses, velocity planes).
# Can be rerun at any time: the links follow the case folders.
set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNS="$(cd "${1:-runs_seat}" && pwd)" || { echo "no folder ${1:-runs_seat}"; exit 1; }
OUT="${2:-seat_reports}"
DESIGNS="${DESIGNS:-cad fillet rampup rampdn lip}"
OPENINGS="${OPENINGS:-5 10 25}"
VALVE="${VALVE:-FPC8}"

mkdir -p "$OUT/runs"
find "$OUT/runs" -maxdepth 1 -type l -delete
n=0
for d in $DESIGNS; do
  name="$VALVE"; [ "$d" != cad ] && name="$VALVE-$d"
  for o in $OPENINGS; do
    tag=$(printf "%s_open%03d" "$name" "$o")
    if [ -d "$RUNS/$tag" ]; then
      ln -sfn "$RUNS/$tag" "$OUT/runs/$tag"; n=$((n + 1))
    else
      echo "  (no $tag in $RUNS)"
    fi
  done
done
echo "$n cases linked in $OUT/runs"
[ "$n" -gt 0 ] || exit 1

python3 "$SCRIPT_DIR/seat_stats.py" "$OUT/runs" -o "$OUT/seat_report.pdf"
python3 "$SCRIPT_DIR/postprocess_valves.py" "$OUT/runs" -o "$OUT/forces_report.pdf"
echo "Reports in $OUT/"
