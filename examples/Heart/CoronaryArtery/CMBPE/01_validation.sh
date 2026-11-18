#!/usr/bin/env bash
# Validation: truncation benchmark of the reduced closures (RmuAxi).
#
# For every blood, inlet radius (Womersley sweep at fixed shear rates) and
# reference shear rate, the stepped pipe is solved in full and on the
# truncated geometry with each closure. The problems are two-dimensional and
# small, so every run uses one core and they run concurrently. A convergence
# check of the reference (time step halved, mesh refined) is added for the
# most severe case.
#
#   CORES=60 ./01_validation.sh
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

OUT="${OUT:-$OUT_ROOT/validation}"
TUBES="${TUBES:-1:10,0.6:30,0.35:30}"
TRUNC="${TRUNC:-13}"
SCALES="${SCALES:-1e-4 1e-3 3e-3}"   # inlet radius (m): Wo ~ 0.05-0.1, 0.5-1, 1.5-3
GAMMAS="${GAMMAS:-100 400}"          # mean apparent shear rate of the narrowest tube (1/s)
AMPLITUDE="${AMPLITUDE:-1.2}"        # p_in = Pbar [1 + A sin wt]; A > 1 reverses the flow
PERIOD="${PERIOD:-0.85}"
STEPS="${STEPS:-200}"
CYCLES="${CYCLES:-4}"
CLOSURES="${CLOSURES:-tree single rest plateau}"
MESH_ARGS="${MESH_ARGS:---nr 16 --nmin 8 --cx 0.5}"
FINE_MESH_ARGS="${FINE_MESH_ARGS:---nr 32 --nmin 16 --cx 0.25}"
CONV_STATE="${CONV_STATE:-hyperviscous}"
CONV_SCALE="${CONV_SCALE:-1e-4}"
CONV_GAMMA="${CONV_GAMMA:-100}"

[[ -x "$RMU_EXE" ]] || { echo "RmuAxi not found at $RMU_EXE: run 00_build.sh"; exit 1; }
AXI_DIR="$REPO/examples/Heart/CoronaryArtery/RmuAxi"
mkdir -p "$OUT"

REDUCED="$(python3 - "$TUBES" "$TRUNC" <<'EOF'
import sys
tubes = [tuple(map(float, t.split(":"))) for t in sys.argv[1].split(",") if t]
xT = float(sys.argv[2]); x = 0.0; out = []
for R, L in tubes:
    if x + L > xT + 1e-12:
        out.append((R, L - max(0.0, xT - x)))
    x += L
print(",".join("%g:%g" % t for t in out))
EOF
)"

mkdir -p "$OUT/mesh" "$OUT/mesh_fine"
python3 "$AXI_DIR/make_rmu_axi_mesh.py" --tubes "$TUBES" --truncation "$TRUNC" \
  --prefix "$OUT/mesh/rmu" $MESH_ARGS
python3 "$AXI_DIR/make_rmu_axi_mesh.py" --tubes "$TUBES" --truncation "$TRUNC" \
  --prefix "$OUT/mesh_fine/rmu" $FINE_MESH_ARGS

JOBS="$OUT/jobs.tsv"
: > "$JOBS"

axi_cmd() {   # mesh closure state scale gamma steps
  echo "$RMU_EXE -rmu_mesh $1 -rmu_output_prefix . -rmu_closure $2" \
       "-rmu_full_tubes $TUBES -rmu_reduced_tubes $REDUCED -rmu_truncation $TRUNC" \
       "-rmu_scale $4 -rmu_gamma_ref $5 -rmu_amplitude $AMPLITUDE -rmu_period $PERIOD" \
       "-rmu_steps_per_cycle $6 -rmu_cycles $CYCLES $(cy_rmu "$3")"
}

add_group() {   # root meshdir state scale gamma steps tag
  local root="$1" mdir="$2" state="$3" scale="$4" gamma="$5" steps="$6" tag="$7"
  local base="$root/s${scale}_g${gamma}/$state"
  add_job "$JOBS" "$tag/full" 1 "$base/full" - "$(axi_cmd "$mdir/rmu_full.mesh" none "$state" "$scale" "$gamma" "$steps")"
  for c in $CLOSURES; do
    add_job "$JOBS" "$tag/$c" 1 "$base/$c" - "$(axi_cmd "$mdir/rmu_trunc.mesh" "$c" "$state" "$scale" "$gamma" "$steps")"
  done
}

# The full runs are listed first: they are the longest.
for scale in $SCALES; do
  for gamma in $GAMMAS; do
    for state in $STATES; do
      add_group "$OUT/main" "$OUT/mesh" "$state" "$scale" "$gamma" "$STEPS" "main/s${scale}_g${gamma}/$state"
    done
  done
done

# Convergence of the reference and of the closure errors.
add_group "$OUT/conv_dt" "$OUT/mesh" "$CONV_STATE" "$CONV_SCALE" "$CONV_GAMMA" "$((2 * STEPS))" \
  "conv_dt/$CONV_STATE"
add_group "$OUT/conv_h" "$OUT/mesh_fine" "$CONV_STATE" "$CONV_SCALE" "$CONV_GAMMA" "$STEPS" \
  "conv_h/$CONV_STATE"

cat > "$OUT/campaign.txt" <<EOF
tubes=$TUBES truncation=$TRUNC reduced=$REDUCED
scales=$SCALES gammas=$GAMMAS amplitude=$AMPLITUDE period=$PERIOD
steps=$STEPS cycles=$CYCLES closures=$CLOSURES
mesh=$MESH_ARGS fine_mesh=$FINE_MESH_ARGS
convergence: state=$CONV_STATE scale=$CONV_SCALE gamma=$CONV_GAMMA (dt/2 and fine mesh)
EOF

echo "$(wc -l < "$JOBS") jobs on $CORES cores"
python3 "$CMBPE_DIR/queue.py" "$JOBS" --cores "$CORES" || true

for root in main conv_dt conv_h; do
  echo "== $root"
  python3 "$AXI_DIR/compare_rmu_axi.py" "$OUT/$root" --cycles "$CYCLES" | tee "$OUT/$root/summary.md" || true
done
# Reference convergence: the full runs of conv_dt/conv_h against main.
python3 - "$OUT" "$CONV_STATE" "$CONV_SCALE" "$CONV_GAMMA" "$CYCLES" <<'EOF' | tee "$OUT/reference_convergence.md"
import csv, sys, os
out, state, scale, gamma, cycles = sys.argv[1:6]
cycles = int(cycles)
def stats(path):
    rows = [r for r in csv.DictReader(open(path)) if int(float(r["cycle"])) == cycles - 1]
    q = [abs(float(r["qIn"])) for r in rows]
    e = [float(r["dissipation"]) for r in rows]
    w = [abs(float(r["wssParent"])) for r in rows]
    return sum(q) / len(q), sum(e) / len(e), sum(w) / len(w)
g = "s%s_g%s/%s/full/RmuAxi.csv" % (scale, gamma, state)
ref = os.path.join(out, "main", g)
if not os.path.isfile(ref):
    sys.exit(0)
r = stats(ref)
print("| reference | mean abs Q | mean dissipation | mean abs WSS parent |")
print("|---|---|---|---|")
print("| main | %.6g | %.6g | %.6g |" % r)
for tag in ("conv_dt", "conv_h"):
    p = os.path.join(out, tag, g)
    if os.path.isfile(p):
        s = stats(p)
        print("| %s | %.6g (%+.3f%%) | %.6g (%+.3f%%) | %.6g (%+.3f%%) |" % (
            tag, s[0], 100 * (s[0] / r[0] - 1), s[1], 100 * (s[1] / r[1] - 1), s[2], 100 * (s[2] / r[2] - 1)))
EOF
