#!/usr/bin/env bash
# run_valves_seat100b.sh
#
# Seat study at 100 % opening, FPC8 and FPC16: "fillet" and "lipr" (the
# "lip" seat with the floor corners rounded as "fillet" and a small rounding
# at the land-ramp edges; seat_plan.py draws both). Pulp as a Bingham
# plastic, flow rates of the specification (FPC8 188.8 m^3/h, FPC16
# 873 m^3/h), atmospheric outlet, PSPG tau without its 2/dt term
# (TAU_P_DT_FACTOR=0, see run_valves_seat100.sh). Meshes
# (resources/examples/Valve/seat/):
#   python3 make_meshes.py FPC8.stp FPC16.stp --only FPC8 --openings 100 \
#       --seat lipr --hm 0.032 --wake-end 4.5 --hs 0.006 --outdir seat
#   python3 make_meshes.py FPC8.stp FPC16.stp --only FPC16 --openings 100 \
#       --seat fillet,lipr --hm 0.032 --wake-end 4.5 --hs 0.0046 --outdir seat
# (FPC8-fillet_open100 is the one of run_valves_seat100.sh; the FPC16 seat
# zone is 0.0046 D = 1.8 mm, the same fraction of its 5 mm groove as the
# 1.2 mm of the FPC8 in its 3.4 mm one.)
#
# Usage, from the folder that holds runs_seat/ (finished cases are skipped):
#   NP=60 NJOBS=4 bash examples/Valve/run_valves_seat100b.sh
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export VALVES="${VALVES:-FPC8-fillet FPC8-lipr FPC16-fillet FPC16-lipr}"
export OPENINGS="${OPENINGS:-100}"
export PRESSURE_OPENINGS="${PRESSURE_OPENINGS:-}"
export Q_FPC8="${Q_FPC8:-188.8}"
export Q_FPC16="${Q_FPC16:-873}"
export TAU_P_DT_FACTOR="${TAU_P_DT_FACTOR:-0}"
exec bash "$SCRIPT_DIR/run_valves_seat.sh"
