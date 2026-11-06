#!/usr/bin/env bash
# run_valves_seat.sh
#
# Seat (seal) design study of the FPC8 with the pulp (Bingham, as in
# run_valves_bingham.sh): the CAD seat against three alternatives that keep
# the flat floor and the vertical walls (what seals and retains the knife)
# and leave less dead pocket for the solids, at 10 and 25 % opening. Meshes:
#   python3 make_meshes.py FPC8.stp FPC16.stp --only FPC8 --openings 10,25 \
#       --seat cad,fillet,chamfer,upramp --hm 0.032 --wake-end 4.5 --hs 0.006 --outdir seat
# (resources/examples/Valve/seat/), seat zone meshed at 1.2 mm.
#
#   FPC8          CAD seat: rectangular groove 10 mm x 3.4 mm (reference)
#   FPC8-fillet   floor corners rounded R 1.5 mm (knife tip corners alike)
#   FPC8-chamfer  lead-in chamfers 30 deg, 1.2 mm deep, on both top edges
#   FPC8-upramp   upstream wall replaced by a 15 deg ramp; downstream wall
#                 vertical (unidirectional)
#
# Same pulp and flow rate as run_valves_bingham.sh (188.8 m^3/h). Note that
# at 10 % the imposed flow rate forces a jet of ~55 m/s and a pressure drop
# of ~2.7 MPa, above the 150 psi rating: PRESSURE_OPENINGS=10 runs the 10 %
# cases driven by the line pressure instead (150 psi -> 0, the flow rate is
# then a result, ~110 m^3/h).
#
# Shorter runs than the force study: 14 D/U (1.75 s), average from 6 D/U;
# the flow in the seat settles in a few D/U. XDMF every 10 steps (the seat
# statistics are time averages of the saved fields: seat_stats.py).
#
# Usage, from a scratch folder:
#   bash examples/Valve/run_valves_seat.sh
#   NJOBS=2 bash examples/Valve/run_valves_seat.sh          # two cases at a time
#   PRESSURE_OPENINGS=10 bash examples/Valve/run_valves_seat.sh
# Then
#   python3 examples/Valve/seat_stats.py runs_seat -o seat_report.pdf
#   python3 examples/Valve/postprocess_valves.py runs_seat -o seat_forces.pdf
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RODIN_ROOT="${RODIN_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"

export NP="${NP:-8}"
export NJOBS="${NJOBS:-1}"
export VALVES="${VALVES:-FPC8 FPC8-fillet FPC8-chamfer FPC8-upramp}"
export OPENINGS="${OPENINGS:-10 25}"
export MESH_DIR="${MESH_DIR:-$RODIN_ROOT/resources/examples/Valve/seat}"
export OUT_DIR="${OUT_DIR:-$PWD/runs_seat}"

export MODEL=bingham
export RHO="${RHO:-1770}"
export TAU_Y="${TAU_Y:-7.0}"
export MU_P="${MU_P:-0.0182}"
export Q_FPC8="${Q_FPC8:-188.8}"
export POUT="${POUT:-0}"
export PLINE="${PLINE:-1.0342e6}"
export PRESSURE_OPENINGS="${PRESSURE_OPENINGS:-}"
export SOLVER="${SOLVER:-iterative}"

export T_CONV="${T_CONV:-14}"
export T_AVG_CONV="${T_AVG_CONV:-6}"
export OUTPUT_EVERY="${OUTPUT_EVERY:-10}"

exec bash "$SCRIPT_DIR/run_valves.sh"
