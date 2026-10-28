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
# Open gate: velocity inlet (developed profile, U = UMEAN, or the flow rate
# Q_<valve> in m^3/h), outlet at POUT (default the line pressure, 150 psi).
# Closed gate: pressure inlet, PIN_CLOSED (150 psi) upstream and POUT_CLOSED
# (0) downstream (hydrostatic, a few steps suffice).
#
# Fluid: MODEL=powerlaw (RHO, K, N; default water) or MODEL=bingham (RHO,
# TAU_Y, MU_P; PAPA_M optional, see KnifeGateValve.cpp).
# run_valves_bingham.sh sets these for the pulp of the 8" valve.
#
# A case whose case.log ends with the final summary is skipped, so the
# script can be relaunched after an interruption. A failed case is reported
# and the script moves on.
#
# Usage (from anywhere):
#   bash examples/Valve/run_valves.sh
# Overrides through the environment, e.g.
#   NP=8 UMEAN=2.5 T_CONV=30 VALVES="FPC16" OPENINGS="25 50" bash examples/Valve/run_valves.sh
# NJOBS > 1 runs that many cases at the same time with NP/NJOBS processes
# each (more throughput than one case on all cores once a case has fewer
# than ~30 k cells per process).
#
# Times are set in convective units D/U (D the bore, U the mean velocity),
# so FPC8 and FPC16 get the same number of flow passes and of time steps:
#   t_end = T_CONV D/U,  average over [T_AVG_CONV D/U, t_end],  dt0 = DT_CONV D/U.
# FPC8 at 3 m/s: D/U = 0.068 s; FPC16: 0.127 s. With the defaults (20, 8,
# 0.0075): ~2700 steps per case, the first 8 D/U (the jet and the separated
# region behind the gate settling) discarded, 12 D/U averaged, i.e. some
# 10 periods of the shear-layer oscillation at the gate edge (St ~ 0.2 on
# the gap height). Forces are written to the case CSV at EVERY step; the
# summary gives the mean, the standard deviation and the standard error of
# the mean (batch means): if the standard error is not small against the
# mean, run longer (T_CONV) rather than trusting the number.
set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RODIN_ROOT="${RODIN_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
EXE="${EXE:-$RODIN_ROOT/build/examples/Valve/KnifeGateValve}"
MESH_DIR="${MESH_DIR:-$RODIN_ROOT/resources/examples/Valve}"
OUT_DIR="${OUT_DIR:-$PWD/runs}"

NP="${NP:-8}"                         # MPI processes in total (8 on the Mac)
NJOBS="${NJOBS:-1}"                   # cases run at the same time (NP/NJOBS processes each)
VALVES="${VALVES:-FPC8 FPC16}"
OPENINGS="${OPENINGS:-0 25 50 75 100}"
OUTPUT_EVERY="${OUTPUT_EVERY:-10}"    # XDMF every N time steps

# Water at 20 C in a 150 psi line (defaults of the executable, repeated here
# so the summary records them).
RHO="${RHO:-998.2}"
K="${K:-1.002e-3}"                    # Pa s^n (= mu for n = 1)
N="${N:-1}"
UMEAN="${UMEAN:-3.0}"                 # m/s, open gate
Q_FPC8="${Q_FPC8:-}"                  # m^3/h; when set, replaces UMEAN for FPC8
Q_FPC16="${Q_FPC16:-}"                # m^3/h; when set, replaces UMEAN for FPC16
MODEL="${MODEL:-powerlaw}"            # powerlaw | bingham
TAU_Y="${TAU_Y:-0}"                   # Pa, Bingham yield stress
MU_P="${MU_P:-0}"                     # Pa s, Bingham plastic viscosity
PAPA_M="${PAPA_M:-}"                  # s, Papanastasiou exponent (empty: executable default)
PLINE="${PLINE:-1.0342e6}"            # Pa, 150 psi
POUT="${POUT:-$PLINE}"                # Pa, outlet pressure with the gate open
PIN_CLOSED="${PIN_CLOSED:-$PLINE}"    # Pa, upstream pressure held by the closed gate
POUT_CLOSED="${POUT_CLOSED:-0}"       # Pa, downstream of the closed gate
T_CONV="${T_CONV:-20}"                # t_end in units of D/U (open gate)
T_AVG_CONV="${T_AVG_CONV:-8}"         # start of the time average, units of D/U
DT_CONV="${DT_CONV:-0.0075}"          # initial dt, units of D/U (adaptive afterwards)
DTMAX_CONV="${DTMAX_CONV:-0.015}"     # largest dt, units of D/U (CFL ~ 6 in the 25 % jet)
NBATCH="${NBATCH:-6}"                 # batches for the standard error of the mean
TEND_CLOSED="${TEND_CLOSED:-2e-3}"    # s, closed gate (hydrostatic, no ramp: a few steps)
EXTRA="${EXTRA:-}"                    # any other -kv_* or PETSc options

MPIRUN="${MPIRUN:-mpirun}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"   # MPI only: no oversubscription

[ -x "$EXE" ] || { echo "Executable not found: $EXE" >&2;
  echo "Build it with: cmake --build build -j --target KnifeGateValve" >&2; exit 1; }
mkdir -p "$OUT_DIR"
NPJ=$(( NP / NJOBS )); [ "$NPJ" -ge 1 ] || NPJ=1

if [ "$MODEL" = "bingham" ]; then
  fluid=(-kv_model bingham -kv_rho "$RHO" -kv_tau_y "$TAU_Y" -kv_mu_p "$MU_P")
  [ -n "$PAPA_M" ] && fluid+=(-kv_papa_m "$PAPA_M")
else
  fluid=(-kv_model powerlaw -kv_rho "$RHO" -kv_K "$K" -kv_n "$N")
fi

failed=()
pids=(); ptags=()
# Waits for the background cases until fewer than $1 are running (portable:
# the bash 3.2 of macOS has no wait -n).
throttle() {
  while :; do
    [ "${#pids[@]}" -eq 0 ] && return
    local keep_p=() keep_t=() i
    for i in "${!pids[@]}"; do
      if kill -0 "${pids[$i]}" 2>/dev/null; then
        keep_p+=("${pids[$i]}"); keep_t+=("${ptags[$i]}")
      else
        wait "${pids[$i]}" || failed+=("${ptags[$i]}")
      fi
    done
    pids=("${keep_p[@]+"${keep_p[@]}"}"); ptags=("${keep_t[@]+"${keep_t[@]}"}")
    [ "${#pids[@]}" -lt "$1" ] && return
    sleep ${THROTTLE_SLEEP:-20}
  done
}

run_case() {  # tag dir log bc...
  local tag="$1" dir="$2" log="$3"; shift 3
  echo "[$tag] start $(date '+%F %T')  ($NPJ processes)  -> $dir"
  (
    cd "$dir" &&
    $MPIRUN -n "$NPJ" "$EXE" \
      -kv_mesh "$mesh" -kv_out "$tag" \
      "${fluid[@]}" \
      "$@" \
      -kv_steady_tol 0 \
      -kv_output_every "$OUTPUT_EVERY" \
      $EXTRA
  ) > "$log" 2>&1
  local status=$?
  if [ $status -ne 0 ]; then
    echo "[$tag] FAILED (exit $status), see $log"; return 1
  fi
  echo "[$tag] done  $(date '+%F %T')"
}

for valve in $VALVES; do
  for o in $OPENINGS; do
    tag=$(printf "%s_open%03d" "$valve" "$o")
    mesh="$MESH_DIR/$tag.mesh"
    dir="$OUT_DIR/$tag"
    log="$dir/case.log"
    if [ ! -f "$mesh" ]; then
      echo "[$tag] mesh not found: $mesh (skipped)"; failed+=("$tag"); continue
    fi
    mkdir -p "$dir"
    if [ -f "$log" ] && grep -qE "average over the last|Steady state reached|End of the time interval" "$log"; then
      echo "[$tag] already done (skipped)"; continue
    fi
    case "$valve" in
      FPC8)  D=0.2032; Qh="$Q_FPC8" ;;
      FPC16) D=0.382;  Qh="$Q_FPC16" ;;
      *) echo "[$tag] unknown bore for $valve (skipped)"; failed+=("$tag"); continue ;;
    esac
    if [ -n "$Qh" ]; then
      Qs=$(python3 -c "print($Qh/3600.0)")                     # m^3/s
      U=$(python3 -c "import math; print($Qs/(math.pi*$D**2/4))")
      inflow=(-kv_q "$Qs")
    else
      U="$UMEAN"
      inflow=(-kv_umean "$UMEAN")
    fi
    tc=$(python3 -c "print($D/abs($U))")                     # D/U, s
    tend=$(python3 -c "print($T_CONV*$tc)")
    tavg=$(python3 -c "print($T_AVG_CONV*$tc)")
    dt0=$(python3 -c "print($DT_CONV*$tc)")
    dtmax=$(python3 -c "print($DTMAX_CONV*$tc)")
    printf "D=%s\nU=%s\ntc=%s\ntavg=%s\ntend=%s\n" "$D" "$U" "$tc" "$tavg" "$tend" > "$dir/case.env"

    if [ "$o" -eq 0 ]; then
      bc=(-kv_inlet pressure -kv_pin "$PIN_CLOSED" -kv_pout "$POUT_CLOSED" -kv_ramp 0
          -kv_tend "$TEND_CLOSED" -kv_tavg 0)
    else
      bc=(-kv_inlet velocity "${inflow[@]}" -kv_pout "$POUT" -kv_tend "$tend" -kv_tavg "$tavg"
          -kv_dt "$dt0" -kv_dtmax "$dtmax")
    fi

    if [ "$NJOBS" -gt 1 ]; then
      throttle "$NJOBS"
      run_case "$tag" "$dir" "$log" "${bc[@]}" &
      pids+=($!); ptags+=("$tag")
    else
      run_case "$tag" "$dir" "$log" "${bc[@]}" || failed+=("$tag")
    fi
  done
done
[ "$NJOBS" -gt 1 ] && throttle 1

# ---- Summary: time averages of every case CSV -----------------------------
python3 - "$OUT_DIR" "$NBATCH" "$UMEAN" <<'PY'
import csv, glob, math, os, sys

out, nbatch, umean = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
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
    env = {}
    if os.path.isfile(os.path.join(d, "case.env")):
        for line in open(os.path.join(d, "case.env")):
            k, _, v = line.strip().partition("=")
            env[k] = float(v)
    tc = env.get("tc", 0.0)
    t0 = 0.0 if opening == 0 else env.get("tavg", 0.5 * float(data[-1]["t"]))
    # Steps of t^{n+1} >= t0: the closed gate keeps only its last state.
    win = [r for r in data if float(r["t"]) > t0 + 1e-12] if opening else data[-1:]
    win = win or data[-1:]
    w = [float(r["dt"]) for r in win]
    W = sum(w)

    def sem(key):
        # Standard error of the mean from nbatch consecutive batches of equal
        # duration (batch means): valid when a batch is longer than the
        # correlation time of the signal; with fewer than 2 batches, nan.
        if opening == 0 or len(win) < 2 * nbatch:
            return float("nan") if opening else 0.0
        t = [float(r["t"]) for r in win]
        ta, tb = t[0] - w[0], t[-1]
        edges = [ta + (tb - ta) * i / nbatch for i in range(nbatch + 1)]
        means = []
        for i in range(nbatch):
            sel = [(wi, float(r[key])) for wi, r, ti in zip(w, win, t) if edges[i] < ti <= edges[i + 1] + 1e-12]
            ws = sum(a for a, _ in sel)
            if ws > 0:
                means.append(sum(a * b for a, b in sel) / ws)
        if len(means) < 2:
            return float("nan")
        m = sum(means) / len(means)
        return math.sqrt(sum((x - m) ** 2 for x in means) / (len(means) - 1) / len(means))

    def mean(key):
        return sum(wi * float(r[key]) for wi, r in zip(w, win)) / W

    def std(key):
        m = mean(key)
        return math.sqrt(max(sum(wi * (float(r[key]) - m) ** 2 for wi, r in zip(w, win)) / W, 0.0))

    q, pin, pout = mean("Q_in"), mean("p_in_mean"), mean("p_out_mean")
    row = {
        "valve": valve, "opening_pct": opening,
        "inlet": "pressure" if opening == 0 else "velocity",
        "U_mean_m_s": 0.0 if opening == 0 else env.get("U", umean),
        "t_end_s": float(data[-1]["t"]), "t_avg_from_s": t0, "steps": len(data),
        "steps_averaged": len(win),
        "D_over_U_s": tc, "window_D_over_U": (W / tc) if (tc > 0 and opening) else 0.0,
        "unconverged_picard_steps": open(os.path.join(d, "case.log"), errors="ignore").read().count("accepted UNCONVERGED")
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
    for c in "yz":
        row[f"F{c}_sem_N"] = sem(f"F{c}")
    row["Fy_dyn_sem_N"] = sem("Fy_dyn")
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
        print(f"  {r['valve']:6s} {r['opening_pct']:3d} %   Fz = {r['Fz_N']:11.4g} +- {r['Fz_sem_N']:.2g} N   "
              f"Fy_dyn = {r['Fy_dyn_N']:10.4g} +- {r['Fy_dyn_sem_N']:.2g} N   dp = {r['dp_Pa']:10.4g} Pa   "
              f"window = {r['window_D_over_U']:.3g} D/U")
else:
    print("No case CSV found in", out)
PY

if [ ${#failed[@]} -gt 0 ]; then
  echo "Failed cases (or missing meshes): ${failed[*]}"
  exit 1
fi
