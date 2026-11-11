#!/usr/bin/env bash
# run_valves_seat2.sh
#
# Second round of the seat study (FPC8, pulp as a Bingham plastic): design A
# (floor corners R 1.5 mm, "fillet") and two A + C combinations, a short ramp
# of 25 deg with rounded edges on one side of the groove and the vertical
# wall on the other:
#   FPC8-fillet   A: rectangular groove, floor corners R 1.5 mm
#   FPC8-rampup   A + C, ramp upstream (the flow enters the pocket along it),
#                 vertical wall downstream (takes the pressure load)
#   FPC8-rampdn   A + C, ramp downstream (flow against the ramp), vertical
#                 wall upstream
# at 5, 10 and 25 % opening: 9 cases. Meshes (resources/examples/Valve/seat/):
#   python3 make_meshes.py FPC8.stp FPC16.stp --only FPC8 --openings 10,25 \
#       --seat fillet,rampup,rampdn --hm 0.032 --wake-end 4.5 --hs 0.006 --outdir seat
#   python3 make_meshes.py FPC8.stp FPC16.stp --only FPC8 --openings 5 \
#       --seat fillet,rampup,rampdn --hk 0.01 --hm 0.032 --wake-end 4.5 --hs 0.006 --outdir seat
# (5 %: gap of 6.8 mm under the knife, hence 2 mm cells at the knife).
#
# 5 and 10 % are driven by the line pressure (150 psi upstream, atmospheric
# downstream; the flow rate is a result): the 188.8 m^3/h of the table would
# need several MPa across the valve. Their times are in jet units D/U_jet
# (U_jet = sqrt(2 dp/rho) = 34 m/s, D/U_jet = 6 ms): 0.24 s, average from
# 0.09 s, dt <= 0.3 ms. 25 % keeps the flow rate of the table and the times
# of run_valves_seat.sh (1.75 s, average from 0.75 s).
#
# Usage, from a scratch folder (3 cases at a time with 15 processes each):
#   NP=45 NJOBS=3 bash examples/Valve/run_valves_seat2.sh
# then
#   python3 examples/Valve/seat_stats.py runs_seat -o seat_report.pdf
#   python3 examples/Valve/postprocess_valves.py runs_seat -o seat_forces.pdf
# Same output folder as run_valves_seat.sh: cases already run (FPC8-fillet at
# 10 and 25 %) are skipped.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export VALVES="${VALVES:-FPC8-fillet FPC8-rampup FPC8-rampdn}"
export OPENINGS="${OPENINGS:-5 10 25}"
export PRESSURE_OPENINGS="${PRESSURE_OPENINGS:-5 10}"
exec bash "$SCRIPT_DIR/run_valves_seat.sh"
