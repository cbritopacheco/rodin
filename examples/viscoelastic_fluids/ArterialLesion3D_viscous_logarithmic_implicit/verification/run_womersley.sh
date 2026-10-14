#!/usr/bin/env bash
# =============================================================================
#  run_womersley.sh -- pulsatile Oldroyd-B pipe flow against its exact solution
#  (womersley_ob.py), for the 3D or the axisymmetric driver:
#
#    space  BDF2 at a fine dt (STEPS_SPACE steps per cycle), mesh by mesh;
#    time   on one mesh, BDF1 and BDF2 at 16, 32, 64, 128 steps per cycle.
#
#     bash run_womersley.sh                       # 3D, both studies, one run after another
#     DRIVER=axi bash run_womersley.sh            # axisymmetric driver
#     STUDY=space MESHES="0.7 1.0" bash run_womersley.sh
#     nohup bash run_womersley.sh > womersley.log 2>&1 &
#
#  Flow: Re = 10, Wo = 3, q(t) = 1 + 0.5 sin(omega t), beta = 0.5 (eta_s = eta_p
#  = 1.795 mPa s), lambda omega = 1 (De = lambda/T = 1/(2 pi)), epsilon = 0.
#  For Oldroyd-B the inflow profile of the drivers is exact, so the whole pipe
#  carries the exact solution once the transient has gone: 4 cycles, ramp over
#  half a cycle, the last cycle compared (periodicity residual ~1e-8).
#  Outputs at 8 phases per cycle (STEPS must be multiples of 8); errors over
#  5 <= x/D <= 10, see compare_womersley.py.
#
#  Meshes: 3D  pipe_verification_rf{0.7,1.0,1.4}.mesh (resources; h/D = 0.1/rf)
#          axi make_pipe_axi_mesh.py rf {1,2,4} (h/D = 0.1, 0.05, 0.025), generated
#
#  Environment (defaults):
#    DRIVER=3d | axi       STUDY=all | space | time
#    MESHES  (3d: "0.7 1.0 1.4", axi: "1 2 4")   TIME_MESH (3d: 0.7, axi: 2)
#    STEPS_LIST="16 32 64 128"   STEPS_SPACE=128   ORDERS="1 2"   ORDER_SPACE=2
#    NP=8 processes per run      MAXJOBS=1 runs at a time (1 = one after another)
#    OUT=$HOME/Results/axi_jfm/verification/womersley_<driver>
#    BUILD (build-gcc12 if present, else build)  EXE=<driver>  MPIRUN=mpirun
#    PETSC_OPTS="..."  PYTHON=python3 (numpy, h5py; scipy and matplotlib optional)
#    BDF2_EXTRAP=<mask> -al_bdf2_extrapolate of the drivers (default: theirs)
#  A run that has finished leaves DONE and is not repeated. Results:
#    $OUT/space.txt, $OUT/time.txt (tables with observed orders),
#    $OUT/womersley_table.tex, $OUT/womersley.pdf, $OUT/summary.csv
# =============================================================================
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"
PY="${PYTHON:-$(command -v python3)}"
PATH="$(printf '%s' "$PATH" | tr ':' '\n' | grep -viE 'anaconda|miniconda|conda/' | paste -sd: -)"
export PATH OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export OMPI_MCA_hwloc_base_binding_policy=none

DRIVER="${DRIVER:-3d}"
STUDY="${STUDY:-all}"
if [[ -z "${BUILD:-}" ]]; then
  BUILD=build; [[ -d "$REPO/build-gcc12" ]] && BUILD=build-gcc12
fi
case "$DRIVER" in
  3d)  TARGET=ArterialLesion3D_viscous_logarithmic_implicit
       MESHES="${MESHES:-0.7 1.0 1.4}"; TIME_MESH="${TIME_MESH:-0.7}" ;;
  axi) TARGET=ArterialLesionAxi_viscous_logarithmic_implicit
       MESHES="${MESHES:-1 2 4}"; TIME_MESH="${TIME_MESH:-2}" ;;
  *)   echo "DRIVER must be 3d or axi"; exit 1 ;;
esac
EXE="${EXE:-$REPO/$BUILD/examples/viscoelastic_fluids/$TARGET/$TARGET}"
MPIRUN="${MPIRUN:-mpirun}"
NP="${NP:-8}"
MAXJOBS="${MAXJOBS:-1}"
STEPS_LIST="${STEPS_LIST:-16 32 64 128}"
STEPS_SPACE="${STEPS_SPACE:-128}"
ORDERS="${ORDERS:-1 2}"
ORDER_SPACE="${ORDER_SPACE:-2}"
OUT="${OUT:-$HOME/Results/axi_jfm/verification/womersley_$DRIVER}"
MESHDIR="$REPO/resources/examples/viscoelastic_fluids"

ETA_S=1.795e-3; ETA_P=1.795e-3; RE=10; WO=3; AMP=0.5; DE=0.15915494309189535
CYCLES=4; RAMP=0.5

[[ -x "$EXE" ]] || { echo "Executable not found: $EXE (set BUILD=... or EXE=...)"; exit 1; }
"$PY" -c "import numpy, h5py" 2> /dev/null || { echo "$PY: numpy and h5py are required (set PYTHON=...)"; exit 1; }
grep -aq -- "-al_time_order" "$EXE" \
  || { echo "$EXE does not know -al_time_order: rebuild it with the BDF2 sources"; exit 1; }
mkdir -p "$OUT/meshes"
echo "== $(date '+%d/%m %H:%M')  driver $EXE"
echo "== study=$STUDY  meshes: $MESHES  time mesh: $TIME_MESH  steps: $STEPS_LIST  NP=$NP  MAXJOBS=$MAXJOBS"

mesh_for() {
  local rf="$1" m
  if [[ "$DRIVER" == axi ]]; then
    m="$OUT/meshes/axi_pipe_rf$rf.mesh"
    [[ -f "$m" ]] || "$PY" "$HERE/make_pipe_axi_mesh.py" "$rf" "$m" >&2 || return 1
  else
    m="$MESHDIR/pipe_verification_rf$rf.mesh"
    if [[ ! -f "$m" ]]; then
      m="$OUT/meshes/pipe_verification_rf$rf.mesh"
      [[ -f "$m" ]] || "$PY" "$HERE/../make_lesion_mesh_3d.py" H "$rf" "$m" -setnumber Lu 2 \
        -setnumber Ld 12 -setnumber hb 0.1 -setnumber hl 0.1 -setnumber hw 0.1 \
        -setnumber hwf 0.1 >&2 || return 1
    fi
  fi
  echo "$m"
}

run_case() {   # rf steps order
  local rf=$1 n=$2 o=$3 dir="$OUT/rf${1}_N${2}_o${3}" mesh
  [[ -f "$dir/DONE" ]] && return 0
  (( n % 8 == 0 )) || { echo "steps per cycle must be a multiple of 8 ($n)"; return 1; }
  mesh=$(mesh_for "$rf") || { echo "no mesh rf=$rf"; return 1; }
  mkdir -p "$dir"
  echo "[run ] $(date '+%d %H:%M') rf=$rf N=$n BDF$o np=$NP -> $dir"
  ( cd "$dir" && $MPIRUN -n "$NP" "$EXE" -al_mesh "$mesh" \
      -al_re $RE -al_wo $WO -al_amplitude $AMP -al_de $DE -al_wi 0 -al_ptt_epsilon 0 \
      -al_eta_s $ETA_S -al_eta_p $ETA_P -al_cycles $CYCLES -al_ramp_cycles $RAMP \
      -al_steps_per_cycle "$n" -al_output_every $(( n / 8 )) -al_time_order "$o" \
      ${BDF2_EXTRAP:+-al_bdf2_extrapolate $BDF2_EXTRAP} \
      -al_conformation_its 12 -al_conformation_tol 1e-9 \
      -al_xdmf case -al_csv case.csv ${PETSC_OPTS:-} > run.log 2>&1 ) \
    && { touch "$dir/DONE"; echo "[end ] $(date '+%d %H:%M') rf=$rf N=$n BDF$o"; } \
    || echo "[FAIL] rf=$rf N=$n BDF$o, see $dir/run.log"
}

slot() { while (( $(jobs -rp | wc -l) >= MAXJOBS )); do sleep 10; done; }

# Meshes first, in this shell (the runs may start together).
for rf in $( { [[ "$STUDY" != time ]] && echo $MESHES; [[ "$STUDY" != space ]] && echo $TIME_MESH; } | tr ' ' '\n' | sort -u); do
  mesh_for "$rf" > /dev/null || { echo "no mesh rf=$rf"; exit 1; }
done

queue=()
if [[ "$STUDY" == all || "$STUDY" == space ]]; then
  for rf in $MESHES; do queue+=("$rf $STEPS_SPACE $ORDER_SPACE"); done
fi
if [[ "$STUDY" == all || "$STUDY" == time ]]; then
  for o in $ORDERS; do for n in $STEPS_LIST; do queue+=("$TIME_MESH $n $o"); done; done
fi
seen=" "
for job in "${queue[@]}"; do
  key="${job// /_}"
  [[ "$seen" == *" $key "* ]] && continue
  seen="$seen$key "
  if (( MAXJOBS > 1 )); then
    slot
    run_case $job &
  else
    run_case $job
  fi
done
wait

# ---- comparisons -------------------------------------------------------------
CMP=("$PY" "$HERE/compare_womersley.py" --eta-s $ETA_S --eta-p $ETA_P --amplitude $AMP)
SUM="$OUT/summary.csv.tmp"
rm -f "$SUM"
for d in "$OUT"/rf*_N*_o*; do
  [[ -f "$d/DONE" ]] || continue
  b=$(basename "$d"); rf=${b#rf}; rf=${rf%%_N*}; n=${b#*_N}; n=${n%%_o*}; o=${b##*_o}
  ref="$OUT/rf${rf}_N$(( 2 * n ))_o$o"
  args=("$d")
  [[ -f "$ref/DONE" ]] && args+=(--ref "$ref")
  line=$("${CMP[@]}" "${args[@]}") || { echo "[cmp ] failed for $d"; continue; }
  echo "rf=$rf N=$n order=$o $line" >> "$SUM"
done
[[ -f "$SUM" ]] || exit 0

"$PY" - "$SUM" "$OUT" "$DRIVER" "$STEPS_SPACE" "$ORDER_SPACE" "$TIME_MESH" << 'EOF' | tee "$OUT/tables.txt"
import sys, math, os
rows = [dict(kv.split("=", 1) for kv in l.split()) for l in open(sys.argv[1])]
out, driver, NS, OS, TM = sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], float(sys.argv[6])
f = lambda r, k: float(r[k]) if k in r else float("nan")
keys = ["rf", "N", "order", "dt", "err_u", "err_sigma", "err_sxr", "err_sxx", "err_u_max",
        "periodicity", "dself_u", "dself_sigma"]
with open(os.path.join(out, "summary.csv"), "w") as fh:
    fh.write(",".join(keys) + "\n")
    for r in sorted(rows, key=lambda r: (float(r["rf"]), int(r["order"]), int(r["N"]))):
        fh.write(",".join(r.get(k, "") for k in keys) + "\n")
p = lambda a, b, ha, hb: math.log(a / b) / math.log(ha / hb) if a > 0 and b > 0 else float("nan")
fmt = lambda x: "--" if x != x else ("$%.2f$" % x)
sci = lambda x: "--" if x != x else "$%s{\\times}10^{%d}$" % (("%.2e" % x).split("e")[0], int(("%.2e" % x).split("e")[1]))
pr = lambda x: "" if x != x else "%.2f" % x
tex = []

# space study
sp = sorted([r for r in rows if r["N"] == NS and r["order"] == OS], key=lambda r: float(r["rf"]))
if sp:
    print("\nSPACE (%s, BDF%s, %s steps/cycle): relative space-time L2 errors over the last cycle" % (driver, OS, NS))
    print("%8s %10s %6s %10s %6s %10s %10s %10s" % ("h/D", "e_u", "p", "e_sigma", "p", "e_sxr", "e_sxx", "periodic"))
    tex += ["\\begin{table}", "  \\centering", "  \\begin{tabular}{ccccccc}", "    \\toprule",
            "    $h/D$ & $e_{\\boldsymbol{u}}$ & $p$ & $e_{\\boldsymbol{\\sigma}}$ & $p$ & $e_{\\sigma_{rz}}$ & $e_{\\sigma_{zz}}$ \\\\",
            "    \\midrule"]
    prev = None
    for r in sp:
        h = 0.1 / float(r["rf"])
        pu = p(prev[1], f(r, "err_u"), prev[0], h) if prev else float("nan")
        ps = p(prev[2], f(r, "err_sigma"), prev[0], h) if prev else float("nan")
        print("%8.4f %10.3e %6s %10.3e %6s %10.3e %10.3e %10.1e" % (h, f(r, "err_u"), pr(pu),
              f(r, "err_sigma"), pr(ps), f(r, "err_sxr"), f(r, "err_sxx"), f(r, "periodicity")))
        tex.append("    %.4f & %s & %s & %s & %s & %s & %s \\\\" % (h, sci(f(r, "err_u")), fmt(pu),
                   sci(f(r, "err_sigma")), fmt(ps), sci(f(r, "err_sxr")), sci(f(r, "err_sxx"))))
        prev = (h, f(r, "err_u"), f(r, "err_sigma"))
    tex += ["    \\bottomrule", "  \\end{tabular}",
            "  \\caption{Pulsatile Oldroyd-B pipe flow (%s driver), BDF2 with $T/\\Delta t=%s$: relative"
            " space-time $L^2$ errors against the exact solution over the last period, $5\\le z/D\\le 10$,"
            " and observed orders $p$.}" % ("3D" if driver == "3d" else "axisymmetric", NS),
            "  \\label{tab:womersley-space}", "\\end{table}", ""]

# time study
tm = [r for r in rows if abs(float(r["rf"]) - TM) < 1e-9]
orders = sorted({r["order"] for r in tm})
table = {}
for o in orders:
    rs = sorted([r for r in tm if r["order"] == o], key=lambda r: int(r["N"]))
    if len(rs) < 2:
        continue
    print("\nTIME (%s, mesh rf=%g, BDF%s): d_N = ||x_N - x_2N|| at t = C T; e = error vs exact" % (driver, TM, o))
    print("%6s %10s %10s %6s %10s %6s %10s %10s" % ("N", "dt", "d_u", "p", "d_sigma", "p", "e_u", "e_sigma"))
    prev = None
    for r in rs:
        du, ds, dt = f(r, "dself_u"), f(r, "dself_sigma"), f(r, "dt")
        pu = p(prev[1], du, prev[0], dt) if prev and du == du else float("nan")
        ps = p(prev[2], ds, prev[0], dt) if prev and ds == ds else float("nan")
        print("%6s %10.3e %10.3e %6s %10.3e %6s %10.3e %10.3e" % (r["N"], dt, du, pr(pu), ds, pr(ps),
              f(r, "err_u"), f(r, "err_sigma")))
        table.setdefault(int(r["N"]), {})[o] = (du, pu, ds, ps)
        prev = (dt, du, ds)
if table:
    cols = "".join("cccc" for _ in orders)
    head = " & ".join("\\multicolumn{4}{c}{BDF%s}" % o for o in orders)
    sub = " & ".join("$d_{\\boldsymbol{u}}$ & $p$ & $d_{\\boldsymbol{\\sigma}}$ & $p$" for _ in orders)
    tex += ["\\begin{table}", "  \\centering", "  \\begin{tabular}{c%s}" % cols, "    \\toprule",
            "    & %s \\\\" % head, "    $T/\\Delta t$ & %s \\\\" % sub, "    \\midrule"]
    for N in sorted(table):
        cells = []
        for o in orders:
            du, pu, ds, ps = table[N].get(o, (float("nan"),) * 4)
            cells += [sci(du), fmt(pu), sci(ds), fmt(ps)]
        if all(c == "--" for c in cells):
            continue
        tex.append("    %d & %s \\\\" % (N, " & ".join(cells)))
    tex += ["    \\bottomrule", "  \\end{tabular}",
            "  \\caption{Pulsatile Oldroyd-B pipe flow (%s driver, mesh $h/D=%.3g$): self-convergence in"
            " time, $d=\\|x_{\\Delta t}-x_{\\Delta t/2}\\|/\\|x\\|$ at $t=4T$ over $5\\le z/D\\le 10$, and"
            " observed orders $p$.}" % ("3D" if driver == "3d" else "axisymmetric", 0.1 / TM),
            "  \\label{tab:womersley-time}", "\\end{table}", ""]
open(os.path.join(out, "womersley_table.tex"), "w").write("\n".join(tex) + "\n")
EOF
"$PY" "$HERE/plot_womersley.py" "$OUT" 2> /dev/null && echo "Figure: $OUT/womersley.pdf"
rm -f "$SUM"
echo "Tables: $OUT/tables.txt   LaTeX rows: $OUT/womersley_table.tex   CSV: $OUT/summary.csv"
