#!/usr/bin/env bash
# run_valves_seat3.sh
#
# Third round of the seat study (FPC8, pulp as a Bingham plastic): design D
# ("lip"), bidirectional. Vertical lands 1.6 mm deep above the floor on both
# sides (the knife seals and is retained on them and on the floor whatever
# the flow direction), then 30 deg ramps 1.8 mm deep to the bore on both
# sides, ramp-bore edges rounded R 3 (tangent to the bore), floor corners
# R 0.6. Compared with the chamfer of the first round (B: 30 deg, 1.2 mm,
# vertical walls of 2.2 mm with sharp floor corners), deeper ramps, shorter
# lands, no floor corners. 5, 10 and 25 % opening; meshes
# (resources/examples/Valve/seat/):
#   python3 make_meshes.py FPC8.stp FPC16.stp --only FPC8 --openings 10,25 \
#       --seat lip --hm 0.032 --wake-end 4.5 --hs 0.006 --outdir seat
#   python3 make_meshes.py FPC8.stp FPC16.stp --only FPC8 --openings 5 \
#       --seat lip --hk 0.01 --hm 0.032 --wake-end 4.5 --hs 0.006 --outdir seat
# Same settings as run_valves_seat2.sh (5 and 10 % driven by 150 psi, times
# in jet units; 25 % at 188.8 m^3/h). Same output folder, so seat_stats.py
# compares it with the other designs.
#
# Usage, from the folder that holds runs_seat/:
#   NP=45 NJOBS=3 bash examples/Valve/run_valves_seat3.sh
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export VALVES="${VALVES:-FPC8-lip}"
export OPENINGS="${OPENINGS:-5 10 25}"
export PRESSURE_OPENINGS="${PRESSURE_OPENINGS:-5 10}"
exec bash "$SCRIPT_DIR/run_valves_seat.sh"
