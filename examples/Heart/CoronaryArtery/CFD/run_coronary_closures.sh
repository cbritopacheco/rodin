#!/usr/bin/env bash
# Closure comparison of the coupled LV-0D / coronary-3D model for the three
# Carreau-Yasuda fits of the manuscript (table 1).
#
# Calibration: the bed resistances are those of the healthy blood at its
# resting operating point (-coronary_calibration reference), so the healthy
# R_mu reproduces the prescribed resting flow and every blood perfuses the
# same bed. At fixed blood, three closures are compared:
#   rmu      shear dependent (R_mu);
#   rest     constant, frozen at the resting operating point of that blood
#            (coincides with R_mu at rest: the fair constant comparator);
#   plateau  constant at mu_inf (the comparison of the first submission).
# With MATCH=1 a fourth one, rest_matched, scales the rest closure by secant
# iteration until its mean |Q_in| equals that of R_mu.
#
# Optional studies (0/1):
#   PRA_STUDY   hyperviscous blood at the elevated drainage pressure;
#   DT_STUDY    time-step refinement of the most severe case;
#   MESH_STUDY  mesh refinement at fixed terminal parameters (the finest mesh
#               writes the calibration, the others read it);
#   TREE_STUDY  multi-generation arteriolar limb instead of a single one.
#
# Usage:
#   BUILD=~/rodin/build NP=8 ./run_coronary_closures.sh
#   NP=8 MATCH=1 DT_STUDY=1 ./run_coronary_closures.sh
#   NP=8 MESH_STUDY=1 MESHES="fine:/p/fine.mesh medium:/p/med.mesh coarse:/p/c.mesh" ./run_coronary_closures.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
BUILD="${BUILD:-$REPO/build}"
EXE="${EXE:-$BUILD/examples/Heart/CoronaryArtery/CFD/CoronaryArtery}"
MESH="${MESH:-$REPO/resources/examples/Heart/CoronaryArtery/Output3D.mesh}"
NP="${NP:-8}"
OUT="${OUT:-$PWD/coronary_closures}"

PERIOD=0.85
DT="${DT:-1e-3}"
CYCLES="${CYCLES:-3}"
STATES="${STATES:-healthy diabetic hyperviscous}"
CLOSURES="${CLOSURES:-rmu rest plateau}"
PRA="${PRA:-700}"                  # calibration and normal drainage pressure (Pa)
PRA_HIGH="${PRA_HIGH:-2400}"       # elevated drainage pressure (Pa)
ALPHA="${ALPHA:-0.75}"             # p_im = alpha p_LV
MATCH="${MATCH:-0}"
MATCH_TOL="${MATCH_TOL:-1e-3}"     # relative tolerance on mean |Q_in|
MATCH_MAXIT="${MATCH_MAXIT:-4}"
PRA_STUDY="${PRA_STUDY:-1}"
DT_STUDY="${DT_STUDY:-0}"
DTS="${DTS:-1e-3 5e-4 2.5e-4}"
MESH_STUDY="${MESH_STUDY:-0}"
MESHES="${MESHES:-}"               # "name:path ..." finest first
TREE_STUDY="${TREE_STUDY:-0}"
# Illustrative generations r:v:f (m, m/s, share of the arteriolar drop):
# replace with the morphometry adopted in the manuscript.
TREE="${TREE:-100e-6:15e-3:0.35,50e-6:8e-3:0.35,25e-6:5e-3:0.30}"
EXTRA="${EXTRA:-}"

[[ -x "$EXE" ]] || { echo "CoronaryArtery not found at $EXE"; exit 1; }
[[ -f "$MESH" ]] || { echo "mesh not found at $MESH (git lfs pull?)"; exit 1; }

cy_args() {
  case "$1" in
    healthy)      echo "-coronary_cy_mu0 0.0712616 -coronary_cy_muinf 0.0041   -coronary_cy_lambda 9.2780  -coronary_cy_n 0.3572 -coronary_cy_a 1.8873" ;;
    diabetic)     echo "-coronary_cy_mu0 0.0973528 -coronary_cy_muinf 0.0049622 -coronary_cy_lambda 15.7825 -coronary_cy_n 0.2670 -coronary_cy_a 1.0829" ;;
    hyperviscous) echo "-coronary_cy_mu0 0.3467269 -coronary_cy_muinf 0.0053264 -coronary_cy_lambda 17.4172 -coronary_cy_n 0.2200 -coronary_cy_a 0.6957" ;;
    *) echo "unknown state $1" >&2; exit 1 ;;
  esac
}

closure_args() {   # closure [scale]
  case "$1" in
    rmu)      echo "-coronary_outlet_closure rmu" ;;
    rest)     echo "-coronary_outlet_closure rest -coronary_outlet_constant_scale ${2:-1}" ;;
    plateau)  echo "-coronary_outlet_closure plateau" ;;
    *) echo "unknown closure $1" >&2; exit 1 ;;
  esac
}

nsteps() {   # dt -> steps of CYCLES cycles
  python3 -c "print(int(round($CYCLES * $PERIOD / $1)))"
}

# Mean |Q_in| over the last cycle of a run.
mean_inflow() {
  python3 - "$1/CoronaryArtery.csv" "$PERIOD" "$CYCLES" <<'EOF'
import csv, sys
path, T, n = sys.argv[1], float(sys.argv[2]), int(sys.argv[3])
rows = list(csv.DictReader(open(path)))
t0 = (n - 1) * T
q = [abs(float(r["CoronaryInletFlux"])) for r in rows if float(r["t"]) >= t0 - 1e-12]
print(sum(q) / len(q))
EOF
}

run() {   # dir dt mesh args...
  local dir="$1" dt="$2" mesh="$3"; shift 3
  local steps; steps="$(nsteps "$dt")"
  if [[ -f "$dir/CoronaryArtery.csv" && "$(wc -l < "$dir/CoronaryArtery.csv")" -gt "$steps" ]]; then
    echo "  [skip] $dir"
    return 0
  fi
  mkdir -p "$dir"
  echo "  [run ] $dir"
  # shellcheck disable=SC2068
  mpirun -n "$NP" "$EXE" \
    -coronary_mesh "$mesh" -coronary_output_prefix "$dir" \
    -coronary_dt "$dt" -coronary_nsteps "$steps" \
    -coronary_calibration reference -coronary_pra "$PRA" -coronary_alpha_im "$ALPHA" \
    $@ $EXTRA > "$dir/log.txt" 2>&1 || { echo "  [FAIL] $dir (see $dir/log.txt)"; return 1; }
}

case_dir() {   # state pra closure
  echo "$OUT/pra$2/$1/$3"
}

state_runs() {   # state pra
  local state="$1" pra="$2"
  local op=""
  [[ "$pra" != "$PRA" ]] && op="-coronary_operating_pra $pra"
  for closure in $CLOSURES; do
    run "$(case_dir "$state" "$pra" "$closure")" "$DT" "$MESH" \
      $(cy_args "$state") $(closure_args "$closure") $op \
      -coronary_calibration_out "$(case_dir "$state" "$pra" "$closure")/calibration.csv" || true
  done

  if [[ "$MATCH" == "1" ]]; then
    local rdir; rdir="$(case_dir "$state" "$pra" rmu)"
    [[ -f "$rdir/CoronaryArtery.csv" ]] || run "$rdir" "$DT" "$MESH" $(cy_args "$state") $(closure_args rmu) $op
    local target; target="$(mean_inflow "$rdir")"
    # secant on s: m(s) is the mean |Q_in| of the rest closure scaled by s
    local s0=1 m0 s1 m1 s2 it
    local d0; d0="$(case_dir "$state" "$pra" rest)"
    [[ -f "$d0/CoronaryArtery.csv" ]] || run "$d0" "$DT" "$MESH" $(cy_args "$state") $(closure_args rest 1) $op
    m0="$(mean_inflow "$d0")"
    s1="$(python3 -c "print($s0 * $m0 / $target)")"
    for it in $(seq 1 "$MATCH_MAXIT"); do
      local d1="$OUT/pra$pra/$state/rest_matched_it$it"
      run "$d1" "$DT" "$MESH" $(cy_args "$state") $(closure_args rest "$s1") $op || break
      m1="$(mean_inflow "$d1")"
      echo "  [match] $state pra=$pra it=$it s=$s1 mean|Q|=$m1 target=$target"
      if python3 -c "import sys; sys.exit(0 if abs($m1-$target) <= $MATCH_TOL*abs($target) else 1)"; then
        rm -rf "$OUT/pra$pra/$state/rest_matched"
        cp -r "$d1" "$OUT/pra$pra/$state/rest_matched"
        echo "$s1" > "$OUT/pra$pra/$state/rest_matched/scale.txt"
        break
      fi
      s2="$(python3 -c "print($s1 - ($m1 - $target) * ($s1 - $s0) / ($m1 - $m0))")"
      s0="$s1"; m0="$m1"; s1="$s2"
    done
  fi
}

echo "== closure comparison at P_RA = $PRA Pa"
for state in $STATES; do
  echo "-- $state"
  state_runs "$state" "$PRA"
done

if [[ "$PRA_STUDY" == "1" ]]; then
  echo "== hyperviscous blood at P_RA = $PRA_HIGH Pa"
  state_runs hyperviscous "$PRA_HIGH"
fi

if [[ "$DT_STUDY" == "1" ]]; then
  echo "== time-step study (hyperviscous, P_RA = $PRA_HIGH Pa)"
  for dt in $DTS; do
    for closure in rmu rest; do
      run "$OUT/dt/dt$dt/$closure" "$dt" "$MESH" $(cy_args hyperviscous) $(closure_args "$closure") \
        -coronary_operating_pra "$PRA_HIGH" || true
    done
  done
fi

if [[ "$MESH_STUDY" == "1" ]]; then
  [[ -n "$MESHES" ]] || { echo "MESH_STUDY=1 needs MESHES=\"name:path ...\" (finest first)"; exit 1; }
  echo "== mesh study at fixed terminal parameters (healthy, rmu and rest)"
  cal=""
  for entry in $MESHES; do
    name="${entry%%:*}"; path="${entry#*:}"
    for closure in rmu rest; do
      dir="$OUT/mesh/$name/$closure"
      if [[ -z "$cal" ]]; then
        run "$dir" "$DT" "$path" $(cy_args healthy) $(closure_args "$closure") \
          -coronary_calibration_out "$OUT/mesh/calibration.csv" || true
      else
        run "$dir" "$DT" "$path" $(cy_args healthy) $(closure_args "$closure") \
          -coronary_calibration_in "$cal" || true
      fi
    done
    cal="$OUT/mesh/calibration.csv"
  done
fi

if [[ "$TREE_STUDY" == "1" ]]; then
  echo "== multi-generation arteriolar limb: $TREE"
  for state in $STATES; do
    for closure in rmu rest; do
      run "$OUT/tree/$state/$closure" "$DT" "$MESH" $(cy_args "$state") $(closure_args "$closure") \
        -coronary_arteriolar_tree "$TREE" || true
    done
  done
fi

python3 "$SCRIPT_DIR/compare_coronary_closures.py" "$OUT" --period "$PERIOD" --cycles "$CYCLES"
