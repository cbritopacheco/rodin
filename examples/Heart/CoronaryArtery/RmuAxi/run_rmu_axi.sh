#!/usr/bin/env bash
# Truncation benchmark of the reduced outlet closures (RmuAxi).
#
# For every blood (Carreau-Yasuda fits of table 1), every inlet radius
# (Womersley sweep at fixed shear rates) and every reference shear rate, the
# stepped pipe is solved once in full (the reference) and once per closure on
# the truncated geometry. compare_rmu_axi.py then reads the last cycle of
# every run and reports the error of each closure against the full solution.
#
# Usage (from anywhere; every variable can be overridden):
#   BUILD=~/rodin/build NP=4 ./run_rmu_axi.sh
#   SCALES="1e-4 1e-3 3e-3" GAMMAS="100 400" ./run_rmu_axi.sh
#
# Geometry, in units of the inlet radius: TUBES is the full series R:L,...;
# TRUNC the truncation x_T, which must lie downstream of the first tube. The
# reduced tubes (what the truncated geometry removes) are derived from them.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
BUILD="${BUILD:-$REPO/build}"
EXE="${EXE:-$BUILD/examples/Heart/CoronaryArtery/RmuAxi/RmuAxi}"
NP="${NP:-1}"
OUT="${OUT:-$PWD/rmu_axi_runs}"

TUBES="${TUBES:-1:10,0.6:30,0.35:30}"
TRUNC="${TRUNC:-13}"
SCALES="${SCALES:-1e-4}"          # inlet radius (m)
GAMMAS="${GAMMAS:-100}"           # mean apparent shear rate of the narrowest tube (1/s)
AMPLITUDE="${AMPLITUDE:-1.2}"     # p_in = Pbar [1 + A sin(omega t)]; A > 1 reverses
PERIOD="${PERIOD:-0.85}"
STEPS="${STEPS:-200}"             # steps per cycle
CYCLES="${CYCLES:-4}"
CLOSURES="${CLOSURES:-tree single rest plateau}"
STATES="${STATES:-healthy diabetic hyperviscous}"
MESH_ARGS="${MESH_ARGS:---nr 16 --nmin 8 --cx 0.5}"
EXTRA="${EXTRA:-}"                # extra -rmu_* options for every run

[[ -x "$EXE" ]] || { echo "RmuAxi not found at $EXE (cmake --build \$BUILD --target RmuAxi)"; exit 1; }

launcher() {
  if [[ "$NP" -gt 1 ]]; then echo "mpirun -n $NP"; fi
}

# Carreau-Yasuda fits (table 1), SI units: mu0 muinf lambda n a.
cy_args() {
  case "$1" in
    healthy)      echo "-rmu_cy_mu0 0.0712616 -rmu_cy_muinf 0.0041   -rmu_cy_lambda 9.2780  -rmu_cy_n 0.3572 -rmu_cy_a 1.8873" ;;
    diabetic)     echo "-rmu_cy_mu0 0.0973528 -rmu_cy_muinf 0.0049622 -rmu_cy_lambda 15.7825 -rmu_cy_n 0.2670 -rmu_cy_a 1.0829" ;;
    hyperviscous) echo "-rmu_cy_mu0 0.3467269 -rmu_cy_muinf 0.0053264 -rmu_cy_lambda 17.4172 -rmu_cy_n 0.2200 -rmu_cy_a 0.6957" ;;
    *) echo "unknown state $1" >&2; exit 1 ;;
  esac
}

# Tubes downstream of x_T.
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

mkdir -p "$OUT/mesh"
python3 "$SCRIPT_DIR/make_rmu_axi_mesh.py" --tubes "$TUBES" --truncation "$TRUNC" \
  --prefix "$OUT/mesh/rmu" $MESH_ARGS

cat > "$OUT/campaign.txt" <<EOF
tubes=$TUBES
truncation=$TRUNC
reduced=$REDUCED
scales=$SCALES
gammas=$GAMMAS
amplitude=$AMPLITUDE
period=$PERIOD
steps=$STEPS
cycles=$CYCLES
mesh_args=$MESH_ARGS
extra=$EXTRA
EOF

expected=$((STEPS * CYCLES + 1))

run_case() {   # dir mesh closure state scale gamma
  local dir="$1" mesh="$2" closure="$3" state="$4" scale="$5" gamma="$6"
  if [[ -f "$dir/RmuAxi.csv" && "$(wc -l < "$dir/RmuAxi.csv")" -ge "$expected" ]]; then
    echo "  [skip] $dir"
    return
  fi
  mkdir -p "$dir"
  echo "  [run ] $dir"
  # shellcheck disable=SC2046
  $(launcher) "$EXE" \
    -rmu_mesh "$mesh" -rmu_output_prefix "$dir" -rmu_closure "$closure" \
    -rmu_full_tubes "$TUBES" -rmu_reduced_tubes "$REDUCED" -rmu_truncation "$TRUNC" \
    -rmu_scale "$scale" -rmu_gamma_ref "$gamma" -rmu_amplitude "$AMPLITUDE" \
    -rmu_period "$PERIOD" -rmu_steps_per_cycle "$STEPS" -rmu_cycles "$CYCLES" \
    $(cy_args "$state") $EXTRA > "$dir/log.txt" 2>&1 \
    || { echo "  [FAIL] $dir (see $dir/log.txt)"; return 0; }
}

for scale in $SCALES; do
  for gamma in $GAMMAS; do
    for state in $STATES; do
      base="$OUT/s${scale}_g${gamma}/$state"
      echo "== scale=$scale gamma=$gamma state=$state"
      run_case "$base/full" "$OUT/mesh/rmu_full.mesh" none "$state" "$scale" "$gamma"
      for closure in $CLOSURES; do
        run_case "$base/$closure" "$OUT/mesh/rmu_trunc.mesh" "$closure" "$state" "$scale" "$gamma"
      done
    done
  done
done

python3 "$SCRIPT_DIR/compare_rmu_axi.py" "$OUT" --cycles "$CYCLES"
