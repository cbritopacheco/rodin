#!/usr/bin/env bash
# Thickness sweep 5-50 mm of the Al / RT42 / Al roof sandwich under the
# January weather of Santiago (pvgis_weather.py: PVGIS TMY of Quinta Normal,
# -33.445, -70.683; the hottest real January day of the TMY by default,
# DIA=promedio for the mean day), for pc9: builds NaturalConvectionPCM, runs
# the thicknesses in groups (one group after the other, the cases of a group
# at once, each with its own MPI ranks) and post-processes.
#
#   mkdir -p ~/runs_santiago && cd ~/runs_santiago
#   nohup bash ~/rodin/examples/Heat/PCM_Paine/run_santiago_pc9.sh > sweep.log 2>&1 &
#
# The TMY is downloaded from PVGIS (re.jrc.ec.europa.eu). Without internet on
# pc9, download it elsewhere and pass it:
#   curl -o tmy.json "https://re.jrc.ec.europa.eu/api/v5_3/tmy?lat=-33.445&lon=-70.683&outputformat=json"
#   TMY_JSON=$PWD/tmy.json bash run_santiago_pc9.sh
#
# TANDAS (separated by |): 5-20 mm with 25 ranks each (100), then 25-35 mm
# with 30 (90), then 40-50 mm with 50 (150); a group whose ranks exceed MAXNP
# (150) is refused. Several mpiruns at
# once under Open MPI would bind their ranks to the same cores: MPIRUN below
# uses --bind-to none and one thread per rank (Open MPI syntax).
#
# Mesh (lessons of the Biwole validation): h_x = 5 mm (400 cells over 2 m);
# across the PCM 0.25 mm at BOTH faces growing by <= 1.06 up to <= 1 mm at
# mid-thickness (GRADING=both): 18 cells at 5 mm, 28 at 10, 44 at 20, 56 at 30,
# 74 at 40, 94 at 50 mm; 4 cells per 10 mm aluminium plate. Mushy range 4 K
# (RT42 melts over 38-43 C) with rounded kinks, as the earlier fine sweep;
# dt 30 s, halved automatically where Newton fails (column `substeps` of the
# CSV, lines [dt] in run.log). Any variable of run_pcm_sweep.sh can be
# overridden; finished cases (DONE) are skipped on relaunch.
set -euo pipefail
export LC_ALL=C   # decimal point in awk/printf (es_CL locales print 1,03)

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
if [ -z "${BUILD:-}" ]; then
  BUILD="$REPO/build"
  [ -d "$REPO/build-pcm" ] && BUILD="$REPO/build-pcm"   # pc9
fi
CMAKE_ARGS="${CMAKE_ARGS:-}"
BUILD_JOBS="${BUILD_JOBS:-8}"

echo "[rodin] $(git -C "$REPO" rev-parse --abbrev-ref HEAD) $(git -C "$REPO" rev-parse --short HEAD)"
if [ ! -f "$BUILD/CMakeCache.txt" ]; then
  # shellcheck disable=SC2086
  eval cmake -S "\"$REPO\"" -B "\"$BUILD\"" -DCMAKE_BUILD_TYPE=Release \
    -DRODIN_USE_MPI=ON -DRODIN_USE_PETSC=ON $CMAKE_ARGS
fi
cmake --build "$BUILD" -j "$BUILD_JOBS" --target NaturalConvectionPCM
grep -q pcm_ygrade_both "$REPO/examples/Heat/NaturalConvectionPCM.cpp" \
  || { echo "This NaturalConvectionPCM.cpp has no -pcm_ygrade_both: update the source first" >&2; exit 1; }

export BUILD
export OUT="${OUT:-$PWD/pcm_santiago}"
TANDAS="${TANDAS:-5 10 15 20|25 30 35|40 45 50}"
MAXNP="${MAXNP:-150}"
export NP_MAP="${NP_MAP:-5:25 10:25 15:25 20:25 25:30 30:30 35:30 40:50 45:50 50:50}"
export NP="${NP:-25}"
export SITE="${SITE:-santiago}"
export DIA="${DIA:-caluroso}"
export TMY_JSON="${TMY_JSON:-}"
export ANGLE="${ANGLE:-25}"
export FACING="${FACING:-north}"
export LENGTH="${LENGTH:-2.0}"
export HX="${HX:-0.005}"
export GRADING="${GRADING:-both}"
export HYTOP="${HYTOP:-0.00025}"
export HYMAX="${HYMAX:-0.001}"
export GRADEMAX="${GRADEMAX:-1.06}"
export AL="${AL:-0.01}"
export NAL="${NAL:-4}"
export NDAYS="${NDAYS:-2}"
export DT="${DT:-30}"
export MUSHY="${MUSHY:-4}"
export EXTRA="${EXTRA:--pcm_round 1}"
export RETRIES="${RETRIES:-0}"
export MPIRUN="${MPIRUN:-mpirun --bind-to none -x OMP_NUM_THREADS=1 -x OPENBLAS_NUM_THREADS=1 -x MKL_NUM_THREADS=1}"

np_of() { echo "$NP_MAP" | tr ' ' '\n' | awk -F: -v mm="$1" -v d="$NP" '$1 == mm { n = $2 } END { print (n ? n : d) }'; }
IFS='|' read -r -a tandas <<< "$TANDAS"
for g in "${tandas[@]}"; do
  total=0; count=0
  for mm in $g; do total=$((total + $(np_of "$mm"))); count=$((count + 1)); done
  if [ "$total" -gt "$MAXNP" ]; then
    echo "[group] $g: $total ranks > MAXNP = $MAXNP, refused" >&2
    exit 1
  fi
  echo "[group] $g mm: $count cases, $total ranks ($(date +%d/%m\ %H:%M))"
  THICKNESSES="$g" NJOBS="$count" bash "$HERE/run_pcm_sweep.sh"
  # the next groups reuse the TMY downloaded by the first one
  if [ -z "$TMY_JSON" ] && [ -s "$OUT/tmy_santiago.json" ]; then export TMY_JSON="$OUT/tmy_santiago.json"; fi
done
bash "$HERE/postprocess_pcm.sh" "$OUT" --hours 6 9 12 15 18 21 --all-contours
python3 "$HERE/figures_paper.py" "$OUT"   # paper figures + CSV in $OUT/figures
