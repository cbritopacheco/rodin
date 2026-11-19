#!/usr/bin/env bash
# run_valves_seat100.sh
#
# Seat study (FPC8, pulp as a Bingham plastic) at 100 % opening: the CAD
# seat and the designs fillet, rampup, rampdn and lip, velocity inlet at
# 188.8 m^3/h, atmospheric outlet. Meshes (resources/examples/Valve/seat/):
#   python3 make_meshes.py FPC8.stp FPC16.stp --only FPC8 --openings 100 \
#       --seat cad,fillet,rampup,rampdn,lip --hm 0.032 --wake-end 4.5 --hs 0.006 --outdir seat
#
# TAU_P_DT_FACTOR=0: the PSPG tau without its 2/dt term. With it, tau_p ~
# dt/2 is ~7 times below h/2|u| on the coarse cells (D/10) of the inlet pipe
# (CFL ~ 0.15 at 100 %) and the P1/P1 pressure wiggles next to the velocity
# inlet (the red/blue patches of the first FPC8 100 % run, which also raised
# p_in and zeta). Straight FPC8 pipe test, same flow: in-slab pressure
# spread at the inlet 1.3 kPa -> 0.33 kPa, far field 160 -> 33 Pa, Picard
# unchanged. The convective tau_M keeps 2/dt.
#
# Same output folder as the other rounds, so seat_stats.py and
# seat_reports.sh compare them (OPENINGS="100").
#
# Usage, from the folder that holds runs_seat/:
#   NP=45 NJOBS=5 bash examples/Valve/run_valves_seat100.sh
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export VALVES="${VALVES:-FPC8 FPC8-fillet FPC8-rampup FPC8-rampdn FPC8-lip}"
export OPENINGS="${OPENINGS:-100}"
export PRESSURE_OPENINGS="${PRESSURE_OPENINGS:-}"
export TAU_P_DT_FACTOR="${TAU_P_DT_FACTOR:-0}"
exec bash "$SCRIPT_DIR/run_valves_seat.sh"
