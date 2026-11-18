#!/usr/bin/env bash
# Coupled LV-0D / coronary-3D campaign, scheduled on a core budget.
#
# Jobs (NP cores each, listed by priority; the queue starts every job that
# fits and whose dependencies are done):
#   base     3 bloods x {rmu, rest, plateau} at P_RA = 700 Pa, and the
#            hyperviscous blood x {rmu, rest, plateau} at P_RA = 2400 Pa;
#   tree     3 bloods x {rmu, rest} with a multi-generation arteriolar limb;
#   mesh     healthy x {rmu, rest} on the coarse and fine meshes, reading the
#            calibration of the medium (base) run: fixed terminal parameters;
#   dt       hyperviscous at P_RA = 2400 Pa x {rmu, rest} at smaller dt;
#   matched  per base group, the rest closure scaled until its mean |Q_in|
#            equals that of rmu (secant; 2-4 sequential runs).
#
#   CORES=60 NP=10 ./02_coronary.sh
#
# With 60 cores and NP=10, six runs proceed at once; the fine mesh uses
# NP_FINE. Check the time per step in the first log (grep "Step") to size the
# campaign: a run is CYCLES * 0.85 / DT steps.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

OUT="${OUT:-$OUT_ROOT/coronary}"
NP="${NP:-10}"
NP_FINE="${NP_FINE:-20}"
PERIOD=0.85
DT="${DT:-1e-3}"
DTS="${DTS:-5e-4}"                 # add 2.5e-4 for a third level (4x the cost)
CYCLES="${CYCLES:-3}"
PRA="${PRA:-700}"
PRA_HIGH="${PRA_HIGH:-2400}"
ALPHA="${ALPHA:-0.75}"
OUTPUT_EVERY="${OUTPUT_EVERY:-10}" # XDMF every 10 ms; 0 disables it
MESH_MEDIUM="${MESH_MEDIUM:-$MESH_DIR/coronary_medium.mesh}"
MESH_COARSE="${MESH_COARSE:-$MESH_DIR/coronary_coarse.mesh}"
MESH_FINE="${MESH_FINE:-$MESH_DIR/coronary_fine.mesh}"
# Illustrative generations r:v:f (m, m/s, share of the arteriolar drop);
# replace with the morphometry adopted in the manuscript.
TREE="${TREE:-100e-6:15e-3:0.35,50e-6:8e-3:0.35,25e-6:5e-3:0.30}"
DO_TREE="${DO_TREE:-1}"
DO_MESH="${DO_MESH:-1}"
DO_DT="${DO_DT:-1}"
DO_MATCH="${DO_MATCH:-1}"

[[ -x "$COR_EXE" ]] || { echo "CoronaryArtery not found at $COR_EXE: run 00_build.sh"; exit 1; }
for m in "$MESH_MEDIUM" "$MESH_COARSE" "$MESH_FINE"; do
  [[ "$(head -c 7 "$m")" != "version" ]] || { echo "$m is a Git LFS pointer"; exit 1; }
done
mkdir -p "$OUT"
JOBS="$OUT/jobs.tsv"
: > "$JOBS"

steps() { python3 -c "print(int(round($CYCLES * $PERIOD / $1)))"; }

# PROBE=1: time a few steps of one run at several rank counts and stop. Use it
# to choose NP before launching the campaign.
if [[ "${PROBE:-0}" == "1" ]]; then
  PROBE_NPS="${PROBE_NPS:-10 15 20 30}"
  PROBE_STEPS="${PROBE_STEPS:-12}"
  for np in $PROBE_NPS; do
    d="$OUT/probe/np$np"; mkdir -p "$d"
    (cd "$d" && $MPIRUN -n "$np" "$COR_EXE" -coronary_mesh "$MESH_MEDIUM" -coronary_dt "$DT" \
      -coronary_nsteps "$PROBE_STEPS" -coronary_calibration reference -coronary_pra "$PRA" \
      -coronary_alpha_im "$ALPHA" -coronary_output_every "$OUTPUT_EVERY" $(cy_cor healthy) \
      -coronary_outlet_closure rmu -coronary_output_prefix . > log.txt 2>&1) || true
    python3 - "$d/log.txt" "$np" "$(steps "$DT")" "$CORES" <<'EOF2'
import re, sys
log, np, n, cores = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
t = [float(m.group(1)) for m in re.finditer(r"step = \d+\s+total = ([0-9.eE+-]+) s", open(log, errors="ignore").read())]
if len(t) < 3:
    print("np=%d: not enough steps in %s" % (np, log)); sys.exit()
per = sum(t[1:]) / len(t[1:])
print("np=%3d: %.2f s/step -> %.1f h per run of %d steps; %d runs at once on %d cores"
      % (np, per, per * n / 3600.0, n, cores // np, cores))
EOF2
  done
  exit 0
fi

# Command of one coupled run, without the closure and the output prefix.
base_cmd() {   # np mesh dt state extra...
  local np="$1" mesh="$2" dt="$3" state="$4"; shift 4
  echo "$MPIRUN -n $np $COR_EXE -coronary_mesh $mesh -coronary_dt $dt" \
       "-coronary_nsteps $(steps "$dt") -coronary_calibration reference" \
       "-coronary_pra $PRA -coronary_alpha_im $ALPHA -coronary_output_every $OUTPUT_EVERY" \
       "$(cy_cor "$state") $*"
}

closure_args() {
  case "$1" in
    rmu)     echo "-coronary_outlet_closure rmu" ;;
    rest)    echo "-coronary_outlet_closure rest" ;;
    plateau) echo "-coronary_outlet_closure plateau" ;;
  esac
}

run_job() {   # name np dir after mesh dt state closure extra...
  local name="$1" np="$2" dir="$3" after="$4" mesh="$5" dt="$6" state="$7" closure="$8"; shift 8
  add_job "$JOBS" "$name" "$np" "$dir" "$after" \
    "$(base_cmd "$np" "$mesh" "$dt" "$state" "$@") $(closure_args "$closure")" \
    "-coronary_output_prefix . -coronary_calibration_out calibration.csv"
}

groups=()
for state in $STATES; do
  for c in rmu rest plateau; do
    run_job "pra$PRA/$state/$c" "$NP" "$OUT/pra$PRA/$state/$c" - "$MESH_MEDIUM" "$DT" "$state" "$c"
  done
  groups+=("pra$PRA/$state|$state|")
done
for c in rmu rest plateau; do
  run_job "pra$PRA_HIGH/hyperviscous/$c" "$NP" "$OUT/pra$PRA_HIGH/hyperviscous/$c" - \
    "$MESH_MEDIUM" "$DT" hyperviscous "$c" -coronary_operating_pra "$PRA_HIGH"
done
groups+=("pra$PRA_HIGH/hyperviscous|hyperviscous|-coronary_operating_pra $PRA_HIGH")

if [[ "$DO_TREE" == "1" ]]; then
  for state in $STATES; do
    for c in rmu rest; do
      run_job "tree/$state/$c" "$NP" "$OUT/tree/$state/$c" - "$MESH_MEDIUM" "$DT" "$state" "$c" \
        -coronary_arteriolar_tree "$TREE"
    done
  done
fi

if [[ "$DO_MESH" == "1" ]]; then
  cal="$OUT/pra$PRA/healthy/rmu/calibration.csv"
  mkdir -p "$OUT/mesh"
  ln -sfn "../pra$PRA/healthy" "$OUT/mesh/medium"
  for c in rmu rest; do
    run_job "mesh/coarse/$c" "$NP" "$OUT/mesh/coarse/$c" "pra$PRA/healthy/rmu" \
      "$MESH_COARSE" "$DT" healthy "$c" -coronary_calibration_in "$cal"
    run_job "mesh/fine/$c" "$NP_FINE" "$OUT/mesh/fine/$c" "pra$PRA/healthy/rmu" \
      "$MESH_FINE" "$DT" healthy "$c" -coronary_calibration_in "$cal"
  done
fi

if [[ "$DO_DT" == "1" ]]; then
  mkdir -p "$OUT/dt"
  ln -sfn "../pra$PRA_HIGH/hyperviscous" "$OUT/dt/dt$DT"
  for dt in $DTS; do
    for c in rmu rest; do
      run_job "dt/dt$dt/$c" "$NP" "$OUT/dt/dt$dt/$c" - "$MESH_MEDIUM" "$dt" hyperviscous "$c" \
        -coronary_operating_pra "$PRA_HIGH"
    done
  done
fi

if [[ "$DO_MATCH" == "1" ]]; then
  for g in "${groups[@]}"; do
    IFS='|' read -r gdir state extra <<< "$g"
    add_job "$JOBS" "$gdir/match" "$NP" "$OUT/$gdir/match" "$gdir/rmu,$gdir/rest" \
      "python3 $CMBPE_DIR/matched_secant.py --group $OUT/$gdir --period $PERIOD --cycles $CYCLES" \
      "--run \"$(base_cmd "$NP" "$MESH_MEDIUM" "$DT" "$state" $extra)\""
  done
fi

cat > "$OUT/campaign.txt" <<EOF
dt=$DT dts=$DTS cycles=$CYCLES pra=$PRA pra_high=$PRA_HIGH alpha=$ALPHA
np=$NP np_fine=$NP_FINE cores=$CORES output_every=$OUTPUT_EVERY
meshes: $MESH_COARSE $MESH_MEDIUM $MESH_FINE
tree=$TREE
EOF

echo "$(wc -l < "$JOBS") jobs on $CORES cores"
python3 "$CMBPE_DIR/queue.py" "$JOBS" --cores "$CORES" || true

python3 "$REPO/examples/Heart/CoronaryArtery/CFD/compare_coronary_closures.py" "$OUT" \
  --period "$PERIOD" --cycles "$CYCLES" | tee "$OUT/summary.md" || true
