#!/usr/bin/env bash
# Comparacion de variantes del modelo en la auricula 2D de Qureshi, para pc9.
#
#   bash pc9_run.sh                 # todo: calibracion, matriz, reporte, comparacion
#   NPROC=60 bash pc9_run.sh
#   bash pc9_run.sh dose            # solo el barrido de dosis de factor tisular
#   bash pc9_run.sh report          # rehacer tablas y figuras sin recalcular
#   bash pc9_run.sh look            # figuras rapidas comparando variantes (quicklook.py)
#
# Por que un rango MPI por caso y no 60 rangos en un caso: las mallas tienen
# 10189 celdas. Con 60 rangos serian 170 celdas por rango y el tiempo se iria
# entero en comunicacion. La matriz son 40 corridas independientes, asi que la
# forma eficiente de usar 60 nucleos es correrlas todas a la vez con un rango
# cada una. Si algun dia la malla sube a ~1e5 celdas, conviene 2-4 rangos por
# caso y la mitad de trabajos simultaneos.
set -euo pipefail
cd "$(dirname "$0")"

NPROC=${NPROC:-60}
RESERVA=${RESERVA:-4}                      # nucleos que se dejan libres
JOBS=${JOBS:-$(( NPROC - RESERVA ))}
OUT=${OUT:-$PWD/results_pc9}
CYCLES=${CYCLES:-60}                       # segundos de especies despues de la semilla
WARMUP=${WARMUP:-3}                        # ciclos de solo flujo antes
KIN=${KIN:-qureshi,ptagx,ptagx_mono,ptagx_notm,ptagx_base}
SHAPES=${SHAPES:-CW,WS,BR,CA}
STAGE=${1:-all}

REPO=$(cd ../../../../.. && pwd)
export PTAGX_EXE=${PTAGX_EXE:-$REPO/build/examples/Heart/Thrombosis/LeftAtrium2DPTAGX/LeftAtrium2DPTAGX}
# --bind-to none: cada caso es un mpirun independiente; con el amarre por
# defecto de OpenMPI todos se fijan al nucleo 0 y comparten un solo nucleo.
export PTAGP_MPIRUN=${PTAGP_MPIRUN:-"mpirun -n 1 --bind-to none"}

if [ ! -x "$PTAGX_EXE" ]; then
  echo "No encuentro el ejecutable: $PTAGX_EXE" >&2
  echo "Compilar con:  cd $REPO/build && cmake . && make LeftAtrium2DPTAGX" >&2
  exit 1
fi

# Las mallas van versionadas (texto plano, ~400 kB, sin git-lfs). Si faltan, se
# regeneran desde geometry/, que no necesita red.
if [ ! -f meshes/Q_CW.mesh ]; then
  echo "== regenerando mallas"
  python3 make_qureshi_meshes.py
fi

mkdir -p "$OUT"
echo "== pc9: $NPROC nucleos, $JOBS trabajos simultaneos, salida en $OUT"
echo "== ejecutable: $PTAGX_EXE"

ncase() { python3 - "$1" "$2" <<'PY'
import sys
print(len(sys.argv[1].split(',')) * len(sys.argv[2].split(',')) * 2)
PY
}

case "$STAGE" in
  all|calib|coag)
    if [ "$STAGE" != "coag" ]; then
      echo "== calibracion de la velocidad de pared de la orejuela (CW, SR y FA)"
      python3 run_Q1.py --stage calib --jobs 2 --out "$OUT" --calib-cycles 4 \
        2>&1 | tee -a "$OUT/calib.log"
    fi
    if [ "$STAGE" != "calib" ]; then
      N=$(ncase "$KIN" "$SHAPES")
      echo "== matriz: $N corridas de $CYCLES s ($KIN x $SHAPES x SR,FA)"
      python3 run_Q1.py --stage coag --jobs "$JOBS" --out "$OUT" \
        --kinetics "$KIN" --shapes "$SHAPES" --seeds P2 \
        --cycles "$CYCLES" --warmup "$WARMUP" 2>&1 | tee -a "$OUT/coag.log"
    fi
    if [ "$STAGE" = "calib" ]; then exit 0; fi
    ;&                                       # sigue al reporte
  report)
    echo "== reporte por variante (checks contra la Tabla 1 de Qureshi)"
    python3 run_Q1.py --stage report --out "$OUT" --warmup "$WARMUP" \
      --kinetics "$KIN" --shapes "$SHAPES" 2>&1 | tail -40
    echo "== comparacion entre variantes"
    python3 compare_models.py --out "$OUT" --warmup "$WARMUP"
    ;;
  dose)
    echo "== barrido de la dosis de factor tisular, sin semilla (iniciacion espontanea)"
    python3 dose_scan.py --out "$OUT" --jobs "$JOBS" --shapes CW \
      --cycles "$CYCLES" --warmup "$WARMUP"
    ;;
  look)
    echo "== vista rapida: series de fibrina, campos y barrido de dosis"
    python3 quicklook.py --out "$OUT" --warmup "$WARMUP" "${@:2}"
    exit 0
    ;;
  *)
    echo "etapa desconocida: $STAGE (use all, calib, coag, report, dose o look)" >&2
    exit 2
    ;;
esac

echo
echo "== listo."
echo "   comparacion entre variantes : $OUT/summary.md"
echo "   veredicto por variante      : $OUT/<variante>/verdict.md"
echo "   figuras                     : $OUT/Q1_fibrin_<variante>.png"
