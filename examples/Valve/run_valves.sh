#!/usr/bin/env bash
# run_valves.sh
#
# Runs KnifeGateValve for every valve (FPC8, FPC16) and every opening
# (0, 25, 50, 75, 100 %), one folder per case, and collects the
# time-averaged forces of all cases in one CSV.
#
#   runs/<valve>_open<NNN>/            XDMF every OUTPUT_EVERY steps, case CSV
#   runs/<valve>_open<NNN>/case.log    full log of the case
#   runs/forces_summary.csv            one row per case
#
# Open gate: velocity inlet (developed profile, U = UMEAN), outlet at the
# line pressure (150 psi). Closed gate: pressure inlet, 150 psi upstream and
# 0 downstream (hydrostatic, a few steps suffice).
#
# A case whose case.log ends with the final summary is skipped, so the
# script can be relaunched after an interruption. A failed case is reported
# and the script moves on.
#
# Usage (from anywhere):
#   bash examples/Valve/run_valves.sh
# Overrides through the environment, e.g.
#   NP=12 UMEAN=2.5 TEND=0.8 VALVES="FPC16" OPENINGS="25 50" bash examples/Valve/run_valves.sh
set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RODIN_ROOT="${RODIN_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
EXE="${EXE:-$RODIN_ROOT/build/examples/Valve/KnifeGateValve}"
MESH_DIR="${MESH_DIR:-$RODIN_ROOT/resources/examples/Valve}"
OUT_DIR="${OUT_DIR:-$PWD/runs}"

NP="${NP:-12}"                        # MPI processes
VALVES="${VALVES:-FPC8 FPC16}"
OPENINGS="${OPENINGS:-0 25 50 75 100}"
OUTPUT_EVERY="${OUTPUT_EVERY:-10}"    # XDMF every N time steps

# Water at 20 C in a 150 psi line (defaults of the executable, repeated here
# so the summary records them).
RHO="${RHO:-998.2}"
K="${K:-1.002e-3}"                    # Pa s^n (= mu for n = 1)
N="${N:-1}"
UMEAN="${UMEAN:-3.0}"                 # m/s, open gate
PLINE="${PLINE:-1.0342e6}"            # Pa, 150 psi
TEND="${TEND:-0.6}"                   # s, open gate
TAVG="${TAVG:-0.3}"                   # s, start of the time average
TEND_CLOSED="${TEND_CLOSED:-2e-3}"    # s, closed gate (hydrostatic, no ramp: a few steps)
EXTRA="${EXTRA:-}"                    # any other -kv_* or PETSc options

MPIRUN="${MPIRUN:-mpirun}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"   # MPI only: no oversubscription

[ -x "$EXE" ] || { echo "Executable not found: $EXE" >&2;
  echo "Build it with: cmake --build build -j --target KnifeGateValve" >&2; exit 1; }
mkdir -p "$OUT_DIR"

failed=()
for valve in $VALVES; do
  for o in $OPENINGS; do
    tag=$(printf "%s_open%03d" "$valve" "$o")
    mesh="$MESH_DIR/$tag.mesh"
    dir="$OUT_DIR/$tag"
    log="$dir/case.log"
    if [ ! -f "$mesh" ]; then
      echo "[$tag] mesh not found: $mesh (skipped)"; failed+=("$tag"); continue
    fi
    if [ -f "$log" ] && grep -q "average over the last\|Steady state reached\|End of the time interval" "$log"; then
      echo "[$tag] already done (skipped)"; continue
    fi
    mkdir -p "$dir"

    if [ "$o" -eq 0 ]; then
      bc=(-kv_inlet pressure -kv_pin "$PLINE" -kv_pout 0 -kv_ramp 0 -kv_tend "$TEND_CLOSED" -kv_tavg 0)
    else
      bc=(-kv_inlet velocity -kv_umean "$UMEAN" -kv_pout "$PLINE" -kv_tend "$TEND" -kv_tavg "$TAVG")
    fi

    echo "[$tag] start $(date '+%F %T')  ($NP processes)  -> $dir"
    (
      cd "$dir" &&
      $MPIRUN -n "$NP" "$EXE" \
        -kv_mesh "$mesh" -kv_out "$tag" \
        -kv_rho "$RHO" -kv_K "$K" -kv_n "$N" \
        "${bc[@]}" \
        -kv_output_every "$OUTPUT_EVERY" \
        $EXTRA
    ) > "$log" 2>&1
    status=$?
    if [ $status -ne 0 ]; then
      echo "[$tag] FAILED (exit $status), see $log"; failed+=("$tag")
    else
      echo "[$tag] done  $(date '+%F %T')"
    fi
  done
done

# ---- Summary: time averages of every case CSV -----------------------------
python3 - "$OUT_DIR" "$TAVG" "$UMEAN" "$PLINE" <<'PY'
import csv, glob, math, os, sys

out, tavg, umean, pline = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
rows = []
for d in sorted(glob.glob(os.path.join(out, "*_open*"))):
    tag = os.path.basename(d)
    f = os.path.join(d, tag + ".csv")
    if not os.path.isfile(f):
        continue
    data = list(csv.DictReader(open(f)))
    if not data:
        continue
    valve, opening = tag.split("_open")
    opening = int(opening)
    t0 = 0.0 if opening == 0 else tavg
    win = [r for r in data if float(r["t"]) >= t0 - 1e-12] or data[-1:]
    w = [float(r["dt"]) for r in win]
    W = sum(w)

    def mean(key):
        return sum(wi * float(r[key]) for wi, r in zip(w, win)) / W

    def std(key):
        m = mean(key)
        return math.sqrt(max(sum(wi * (float(r[key]) - m) ** 2 for wi, r in zip(w, win)) / W, 0.0))

    q, pin, pout = mean("Q_in"), mean("p_in_mean"), mean("p_out_mean")
    row = {
        "valve": valve, "opening_pct": opening,
        "inlet": "pressure" if opening == 0 else "velocity",
        "U_mean_m_s": 0.0 if opening == 0 else umean,
        "t_end_s": float(data[-1]["t"]), "t_avg_from_s": t0, "steps": len(data),
        "steps_averaged": len(win),
        "unconverged_picard_steps": open(os.path.join(d, "case.log"), errors="ignore").read().count("Picard did not converge")
                                     if os.path.isfile(os.path.join(d, "case.log")) else -1,
        "Q_m3_s": q, "p_in_mean_Pa": pin, "p_out_mean_Pa": pout, "dp_Pa": pin - pout,
        "zeta": mean("zeta") if opening else 0.0,
    }
    for c in "xyz":
        row[f"F{c}_N"] = mean(f"F{c}")
    for c in "xyz":
        row[f"F{c}_dyn_N"] = mean(f"F{c}_dyn")
    for c in "xyz":
        row[f"F{c}_direct_N"] = mean(f"F{c}_direct")
    for c in "xyz":
        row[f"F{c}_pressure_N"] = mean(f"F{c}_pressure")
    for c in "yz":
        row[f"F{c}_std_N"] = std(f"F{c}")
    rows.append(row)

summary = os.path.join(out, "forces_summary.csv")
if rows:
    rows.sort(key=lambda r: (r["valve"], r["opening_pct"]))
    with open(summary, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        for r in rows:
            wr.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
    print(f"\nSummary ({len(rows)} cases): {summary}")
    for r in rows:
        print(f"  {r['valve']:6s} {r['opening_pct']:3d} %   Fz = {r['Fz_N']:11.4g} N   "
              f"Fy_dyn = {r['Fy_dyn_N']:10.4g} N   dp = {r['dp_Pa']:10.4g} Pa")
else:
    print("No case CSV found in", out)
PY

if [ ${#failed[@]} -gt 0 ]; then
  echo "Failed cases (or missing meshes): ${failed[*]}"
  exit 1
fi
