#!/usr/bin/env bash
# make_videos.sh — presentation videos (volume rendering of the speed, valve
# as transparent solids) of chosen cases, with valve_video.py.
#
#   bash examples/Valve/make_videos.sh [RUNS] [OUTDIR]
#   env: CASES    case folders to film (default: the best seat designs and the
#                 original at 10 and 25 %)
#        QUALITY  draft | normal | high (default high: 1920x1080, 1.5 mm voxels)
#        EXTRA    more valve_video.py options, e.g. "--orbit 80 --opacity 1.5"
#
# Needs: pip install pyvista imageio imageio-ffmpeg h5py
# First check one view quickly:  QUALITY=draft EXTRA=--preview bash make_videos.sh
set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNS="${1:-runs_seat}"
OUT="${2:-videos}"
CASES="${CASES:-FPC8_open025 FPC8-rampup_open025 FPC8-lip_open025 FPC8-rampup_open010 FPC8-lip_open010}"
QUALITY="${QUALITY:-high}"
EXTRA="${EXTRA:-}"
mkdir -p "$OUT"
for c in $CASES; do
  if [ ! -f "$RUNS/$c/$c.xdmf" ]; then echo "  (no $RUNS/$c, skipped)"; continue; fi
  ext=mp4; case " $EXTRA " in *" --preview "*) ext=png ;; esac
  echo "== $c"
  python3 "$SCRIPT_DIR/valve_video.py" "$RUNS/$c" -o "$OUT/$c.$ext" --quality "$QUALITY" $EXTRA \
    2>&1 | grep -v -e "WARN" -e "^$"
done
echo "Videos in $OUT/"
