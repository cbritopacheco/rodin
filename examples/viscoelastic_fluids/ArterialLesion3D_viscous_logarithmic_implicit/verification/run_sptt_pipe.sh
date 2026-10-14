#!/usr/bin/env bash
# =============================================================================
#  run_sptt_pipe.sh -- verification of ArterialLesion3D against the fully
#  developed sPTT pipe flow (sptt_pipe.py; derivation in sptt_pipe_solution.tex).
#
#     bash run_sptt_pipe.sh                       # 4 meshes x 2 fluids, concurrently
#     RF_LIST="0.7 1.0" NP=8 bash run_sptt_pipe.sh
#     nohup bash run_sptt_pipe.sh > verif.log 2>&1 &
#
#  Straight pipe, x in [-2D, 12D], quasi-uniform tetrahedra h/D = 0.1/rf:
#     rf 0.5 -> h = 0.2 D    (  2.1k vertices,   7.6k tets; two elements per radius:
#                              pre-asymptotic, not in the default list)
#     rf 0.7 -> h = 0.143 D  (  5.0k vertices,  20k tets)
#     rf 1.0 -> h = 0.1 D    ( 12k vertices,    53k tets)
#     rf 1.4 -> h = 0.071 D  ( 30k vertices,   141k tets)
#     rf 2.0 -> h = 0.05 D   ( 79k vertices,   399k tets)
#  Fluids: beta = 0.5 (eta_s = eta_p = 1.795 mPa s), Wi = lambda Ubar/D = 0.5,
#     OB    epsilon = 0     (Oldroyd-B; Poiseuille with eta_0)
#     sPTT  epsilon = 0.25  (f = 2.4 at the wall: strongly shear-thinning)
#  Re = 10, Wo = 1, steady inflow (A = 0) ramped over the first cycle; two
#  cycles of 20 steps: steady, Newton |dpsi| ~ 1e-14 at the end.
#  Errors over 5 <= x/D <= 10 (developed region), see compare_sptt_pipe.py.
#
#  All cases run at the same time; the finest mesh gets NP processes and the
#  others NP (rf/rf_max)^3 (at least 2): with the defaults 2+3+9+24 per fluid,
#  76 cores in total; the finest pair takes ~1-2 h on pc9. PARALLEL=0 runs them one after the other with NP each.
#
#  Environment: NP=24  PARALLEL=1  BUILD=build-gcc12 (else build)  EXE3D=<driver>
#               MPIRUN=mpirun  OUT=$HOME/Results/axi_jfm/verification
#               RF_LIST="0.7 1.0 1.4 2.0"  CASES="OB sPTT"  PYTHON=python3
#               TIME_ORDER=2 (BDF2; 1 = BDF1: the steady state is the same)
#               PETSC_OPTS="..." extra solver options (default: the driver's MUMPS LU;
#               the GMRES+ASM/ILU set of run_3d.sh is ~1.7x faster on the linear solve)
#  Output: $OUT/<case>_rf<rf>/ (run.log, case.csv, XDMF), $OUT/summary.csv and
#          the table of errors and observed orders (also in $OUT/table.txt).
# =============================================================================
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"
# python (numpy, h5py) is taken before conda is removed from the PATH.
PY="${PYTHON:-$(command -v python3)}"
PATH="$(printf '%s' "$PATH" | tr ':' '\n' | grep -viE 'anaconda|miniconda|conda/' | paste -sd: -)"
export PATH OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export OMPI_MCA_hwloc_base_binding_policy=none
if [[ -z "${BUILD:-}" ]]; then
  BUILD=build; [[ -d "$REPO/build-gcc12" ]] && BUILD=build-gcc12
fi
EXE="${EXE3D:-$REPO/$BUILD/examples/viscoelastic_fluids/ArterialLesion3D_viscous_logarithmic_implicit/ArterialLesion3D_viscous_logarithmic_implicit}"
MPIRUN="${MPIRUN:-mpirun}"
NP="${NP:-24}"
PARALLEL="${PARALLEL:-1}"
OUT="${OUT:-$HOME/Results/axi_jfm/verification}"
RF_LIST="${RF_LIST:-0.7 1.0 1.4 2.0}"
CASES="${CASES:-OB sPTT}"
MESHDIR="$REPO/resources/examples/viscoelastic_fluids"

ETA_S=1.795e-3; ETA_P=1.795e-3; WI=0.5; RE=10; WO=1

[[ -x "$EXE" ]] || { echo "Executable not found: $EXE (set BUILD=... or EXE3D=...)"; exit 1; }
"$PY" -c "import numpy, h5py" 2> /dev/null || { echo "$PY: numpy and h5py are required (set PYTHON=...)"; exit 1; }
echo "== driver $EXE"
echo "== mpirun $(command -v "$MPIRUN"),  python $PY"
mkdir -p "$OUT"
SUMMARY="$OUT/summary.csv"

mesh_for() {   # rf -> path; generated with gmsh if it is not in the repository
  local rf="$1" m="$MESHDIR/pipe_verification_rf$1.mesh"
  if [[ ! -f "$m" ]]; then
    m="$OUT/pipe_verification_rf$rf.mesh"
    if [[ ! -f "$m" ]]; then
      echo "   generating $m" >&2
      "$PY" "$HERE/../make_lesion_mesh_3d.py" H "$rf" "$m" -setnumber Lu 2 -setnumber Ld 12 \
        -setnumber hb 0.1 -setnumber hl 0.1 -setnumber hw 0.1 -setnumber hwf 0.1 >&2 || return 1
    fi
  fi
  echo "$m"
}

eps_for() { case "$1" in OB) echo 0 ;; sPTT) echo 0.25 ;; *) return 1 ;; esac; }

rf_max=$(printf '%s\n' $RF_LIST | sort -g | tail -1)
np_for() {
  if [[ $PARALLEL == 1 ]]; then
    "$PY" -c "import math; print(max(2, math.ceil($NP * ($1 / $rf_max) ** 3 - 1e-9)))"
  else
    echo "$NP"
  fi
}

run_case() {   # case rf mesh eps np
  local c=$1 rf=$2 mesh=$3 eps=$4 np=$5 dir="$OUT/${1}_rf$2"
  [[ -f "$dir/DONE" ]] && return 0
  mkdir -p "$dir"
  echo "[run ] $(date '+%d %H:%M') $c rf=$rf np=$np -> $dir"
  ( cd "$dir" && "$MPIRUN" -n "$np" "$EXE" -al_mesh "$mesh" \
      -al_re $RE -al_wo $WO -al_amplitude 0 -al_wi $WI -al_ptt_epsilon "$eps" \
      -al_eta_s $ETA_S -al_eta_p $ETA_P -al_cycles 2 -al_steps_per_cycle 20 \
      -al_output_every 1000 -al_conformation_its 12 -al_conformation_tol 1e-9 \
      -al_psi_order 2 -al_time_order "${TIME_ORDER:-2}" \
      -al_xdmf case -al_csv case.csv ${PETSC_OPTS:-} > run.log 2>&1 ) \
    && { touch "$dir/DONE"; echo "[end ] $(date '+%d %H:%M') $c rf=$rf"; } \
    || echo "[FAIL] $c rf=$rf, see $dir/run.log"
}

for rf in $RF_LIST; do
  mesh=$(mesh_for "$rf") || { echo "no mesh for rf=$rf"; exit 1; }
  for c in $CASES; do
    eps=$(eps_for "$c") || { echo "unknown case $c"; exit 1; }
    if [[ $PARALLEL == 1 ]]; then
      run_case "$c" "$rf" "$mesh" "$eps" "$(np_for "$rf")" &
    else
      run_case "$c" "$rf" "$mesh" "$eps" "$NP"
    fi
  done
done
wait

rm -f "$SUMMARY.tmp"
for rf in $RF_LIST; do
  for c in $CASES; do
    eps=$(eps_for "$c"); dir="$OUT/${c}_rf$rf"
    [[ -f "$dir/DONE" ]] || continue
    line=$("$PY" "$HERE/compare_sptt_pipe.py" "$dir" --eta-s $ETA_S --eta-p $ETA_P \
             --eps "$eps" --re $RE --wi $WI --profile "$dir/profile.csv") \
      || { echo "[cmp ] failed for $dir"; continue; }
    echo "case=$c rf=$rf h_over_D=$("$PY" -c "print(0.1/$rf)") $line" >> "$SUMMARY.tmp"
  done
done

# Table and observed orders (log2 of successive error ratios over log2 of h ratios).
[[ -f "$SUMMARY.tmp" ]] || exit 0
"$PY" - "$SUMMARY.tmp" "$SUMMARY" << 'EOF' | tee "$OUT/table.txt"
import sys, math
rows = []
for line in open(sys.argv[1]):
    rows.append(dict(kv.split("=", 1) for kv in line.split()))
keys = ["case", "rf", "h_over_D", "err_u", "err_sigma", "err_G", "err_tauw_mean", "err_Q_region_vs_nominal", "n_wall_zero"]
seen = {}
for r in rows:
    seen[(r["case"], float(r["rf"]))] = r
with open(sys.argv[2], "w") as f:
    f.write(",".join(keys) + "\n")
    for k in sorted(seen):
        f.write(",".join(seen[k].get(x, "") for x in keys) + "\n")
for case in sorted({k[0] for k in seen}):
    rs = [seen[k] for k in sorted(seen) if k[0] == case]
    print("\n%s" % case)
    print("%8s %10s %10s %10s %10s   %6s %6s %6s" % ("h/D", "e_u", "e_sigma", "e_G", "e_tauw", "p_u", "p_sig", "p_G"))
    prev = None
    for r in rs:
        h = float(r["h_over_D"]); eu = float(r["err_u"]); es = float(r["err_sigma"])
        eg = abs(float(r["err_G"])); et = float(r["err_tauw_mean"])
        if prev:
            lr = math.log(prev[0] / h)
            pu, ps, pg = (math.log(prev[1] / eu) / lr, math.log(prev[2] / es) / lr,
                          math.log(prev[3] / eg) / lr if eg > 0 and prev[3] > 0 else float("nan"))
            print("%8.4f %10.3e %10.3e %10.3e %10.3e   %6.2f %6.2f %6.2f" % (h, eu, es, eg, et, pu, ps, pg))
        else:
            print("%8.4f %10.3e %10.3e %10.3e %10.3e" % (h, eu, es, eg, et))
        prev = (h, eu, es, eg)

# LaTeX rows for the table of sptt_pipe_solution.tex (booktabs).
def fmt(x):
    m, e = ("%.2e" % x).split("e")
    return "$%s{\\times}10^{%d}$" % (m, int(e))
with open(sys.argv[2].replace("summary.csv", "table.tex"), "w") as f:
    f.write("\\begin{table}\n  \\centering\n  \\begin{tabular}{cccccccc}\n    \\toprule\n"
            "    $h/D$ & $e_{\\boldsymbol{u}}$ & $p$ & $e_{\\boldsymbol{\\sigma}}$ & $p$\n"
            "      & $|G_h/G-1|$ & $p$ & $\\tau_{w,h}/\\tau_w-1$ (\\%) \\\\\n    \\midrule\n")
    for case in sorted({k[0] for k in seen}):
        rs = [seen[k] for k in sorted(seen) if k[0] == case]
        f.write("    \\multicolumn{8}{l}{%s} \\\\\n" % ("Oldroyd-B ($\\varepsilon=0$)" if case == "OB" else "sPTT ($\\varepsilon=0.25$)"))
        prev = None
        for r in rs:
            h = float(r["h_over_D"]); e = [float(r["err_u"]), float(r["err_sigma"]), abs(float(r["err_G"]))]
            p = ["--"] * 3
            if prev:
                lr = math.log(prev[0] / h)
                p = ["%.2f" % (math.log(a / b) / lr) for a, b in zip(prev[1], e)]
            f.write("    %.3f & %s & %s & %s & %s & %s & %s & %s \\\\\n"
                    % (h, fmt(e[0]), p[0], fmt(e[1]), p[1], fmt(e[2]), p[2], "$%+.1f$" % (100 * float(r["err_tauw_mean"]))))
            prev = (h, e)
    f.write("    \\bottomrule\n  \\end{tabular}\n"
            "  \\caption{Fully developed pipe flow: relative $L^2$ errors of velocity and\n"
            "  polymer stress over $5\\le z/D\\le10$, error of the pressure gradient and mean\n"
            "  error of the wall shear stress, with observed orders $p$.}\n"
            "  \\label{tab:sptt-pipe}\n\\end{table}\n")
EOF
rm -f "$SUMMARY.tmp"
echo "Summary: $SUMMARY   LaTeX rows: $OUT/table.tex"
"$PY" "$HERE/plot_sptt_pipe.py" "$OUT" 2> /dev/null && echo "Figure: $OUT/profiles.pdf"
