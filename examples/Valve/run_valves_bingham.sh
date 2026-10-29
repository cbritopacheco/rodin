#!/usr/bin/env bash
# run_valves_bingham.sh
#
# FPC8 (8") with the pulp as a Bingham plastic, openings 0, 25, 75, 100 %,
# on the meshes refined in the wake (resources/examples/Valve/wake/).
# Upper values of the pulp table (8" row):
#
#   density            1770 kg/m^3   (1.770 t/m^3)
#   yield stress          7.0 Pa
#   plastic viscosity  0.0182 Pa s
#   flow rate           188.8 m^3/h  (U = Q/A = 1.62 m/s in D = 203.2 mm)
#
# Open gate: velocity inlet with flow rate Q (developed profile chosen from
# Re_B and Hanks' critical Reynolds number: Re_B ~ 3.2e4 > Re_c ~ 1.8e4, so
# the 1/7 law), outlet at atmospheric pressure (p = 0 gauge).
# Closed gate (0 %): no flow rate can cross it; pressure inlet with the line
# pressure PIN_CLOSED (150 psi, the CWP of the valve) upstream and 0
# downstream (fluid at rest, the yield stress plays no role).
#
# The runs go to t_end = T_CONV D/U and never stop at a steady state
# (-kv_steady_tol 0); averages over [T_AVG_CONV D/U, t_end].
# D/U = 0.126 s: t_end = 2.5 s, average over the last 1.5 s.
#
# Usage (from anywhere, after building KnifeGateValve):
#   bash examples/Valve/run_valves_bingham.sh
# Overrides, e.g. three cases at a time with ~10 processes each:
#   NJOBS=3 bash examples/Valve/run_valves_bingham.sh
# Every other variable of run_valves.sh can be set the same way.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RODIN_ROOT="${RODIN_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"

export NP="${NP:-32}"
export NJOBS="${NJOBS:-1}"
export VALVES="${VALVES:-FPC8}"
export OPENINGS="${OPENINGS:-0 25 75 100}"
export MESH_DIR="${MESH_DIR:-$RODIN_ROOT/resources/examples/Valve/wake}"
export OUT_DIR="${OUT_DIR:-$PWD/runs_bingham}"

export MODEL=bingham
export RHO="${RHO:-1770}"             # kg/m^3
export TAU_Y="${TAU_Y:-7.0}"          # Pa
export MU_P="${MU_P:-0.0182}"         # Pa s
export Q_FPC8="${Q_FPC8:-188.8}"      # m^3/h
export POUT="${POUT:-0}"              # Pa, atmospheric outlet
export PIN_CLOSED="${PIN_CLOSED:-1.0342e6}"   # Pa, 150 psi held by the closed gate
export POUT_CLOSED="${POUT_CLOSED:-0}"

# Iterative solver: FGMRES + Schur field split with BoomerAMG. The MUMPS LU
# of these meshes (550 k unknowns) needs 25-30 GB; this a fraction of it.
# SOLVER=direct goes back to the LU.
export SOLVER="${SOLVER:-iterative}"

export T_CONV="${T_CONV:-20}"
export T_AVG_CONV="${T_AVG_CONV:-8}"
export OUTPUT_EVERY="${OUTPUT_EVERY:-10}"

exec bash "$SCRIPT_DIR/run_valves.sh"
