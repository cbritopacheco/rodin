#!/bin/bash
# Campaña final del atrio en la malla fina 2 (más fina en el centro):
#   1. recompila Rodin (driver con reinicio) y genera Fine2.mesh si no existe;
#   2. arranque: N durante SPIN ciclos, hasta un flujo desarrollado;
#   3. desde ese mismo estado, N, CY, VE y TVE durante CYCLES ciclos, y los
#      controles numéricos de TVE (devss=0.5, implicit=1) durante CTL ciclos;
#   4. estadística pareada ciclo a ciclo.
# Se puede relanzar: salta lo hecho. Variables: SPIN (8), CYCLES (8), CTL (3),
# TARGET (2600000), REGION ("x,y,z,r" en mm, p. ej. la orejuela; opcional).
set -uo pipefail
source "$HOME/gpu/env.sh"
RUNS=${RUNS:-$HOME/runs}
DIR=${DIR:-$RUNS/final}
SPIN=${SPIN:-8}
CYCLES=${CYCLES:-8}
CTL=${CTL:-3}
TARGET=${TARGET:-2600000}
REGION=${REGION:-}
IN_MESH=${IN_MESH:-$RUNS/fine/Fine.mesh}   # malla que se remalla
GUIDE=${GUIDE:-$RUNS/fine/fine_TVE}        # corrida que guía el refinamiento
SIZE_ARGS=${SIZE_ARGS:-}
EXTRA=${EXTRA:-}                           # claves extra para todas las corridas
R=$RUNS/refine_atrium.py
MMG=$HOME/gpu/mmg/build/bin/mmg3d_O3
BIN=$RODIN_BIN/KokkosAtriumRheologyPETSc
mkdir -p "$DIR" && cd "$DIR" || exit 1
log() { echo "[$(date '+%d %H:%M')] $*"; }

# 0. Detener lo que esté corriendo.
pkill -f run10.sh; pkill -f run_fine.sh; pkill -f run_fine2.sh
pkill -f KokkosAtriumRheologyPETSc; while pgrep -f KokkosAtriumRheologyPETSc > /dev/null; do sleep 5; done
FREE=$(df -BG --output=avail . | tail -1 | tr -dc 0-9)
log "GPU libre; $FREE GB libres en disco"
[ "$FREE" -ge 50 ] || log "AVISO: la campaña escribe unos 40 GB"

# 1. Driver con reinicio y malla.
if ! grep -aq '\[restart\]' "$BIN"; then
  SRC=$HOME/gpu/rodin/examples/Kokkos/Heart/AtriumRheology.cpp
  if [ -f "$RUNS/AtriumRheology.cpp" ]; then
    # Driver nuevo descargado en ~/runs: se compila solo ese ejecutable.
    log "recompilando KokkosAtriumRheologyPETSc con $RUNS/AtriumRheology.cpp"
    cp "$RUNS/AtriumRheology.cpp" "$SRC"
    CXX=$(sed -n 's/^CMAKE_CXX_COMPILER:[A-Z]*=//p' "$HOME/gpu/rodin/build-cuda/CMakeCache.txt")
    NVCC_WRAPPER_DEFAULT_COMPILER=$CXX PETSC_DIR=$HOME/gpu/petsc \
      cmake --build "$HOME/gpu/rodin/build-cuda" -j "$(nproc)" --target KokkosAtriumRheologyPETSc > build.log 2>&1 \
      || { tail -20 build.log; log "la compilación falló: ver $DIR/build.log"; exit 1; }
  else
    log "recompilando Rodin desde la rama model/AtriumKokkos"
    SKIP_TEST=1 bash "$HOME/install_rodin_gpu.sh" build > build.log 2>&1 || { log "la compilación falló: ver $DIR/build.log"; exit 1; }
    source "$HOME/gpu/env.sh"
  fi
  grep -aq '\[restart\]' "$BIN" || { log "el driver no tiene reinicio: falta ~/runs/AtriumRheology.cpp o el push de la rama"; exit 1; }
fi
log "driver con reinicio: ok"
if [ ! -f Fine2.mesh ]; then
  if [ -f "$RUNS/fine2/Fine2.mesh" ]; then
    ln -s "$RUNS/fine2/Fine2.mesh" Fine2.mesh
  else
    F=1.3
    for it in 1 2 3; do
      # shellcheck disable=SC2086
      python3 "$R" size mesh="$IN_MESH" run="$GUIDE" target="$TARGET" fudge="$F" out=size.sol $SIZE_ARGS || exit 1
      "$MMG" -in "$IN_MESH" -sol size.sol -out out.mesh -hgrad 1.3 -hausd 0.05 -v 1 > mmg.log 2>&1 \
        || { log "mmg3d falló"; tail -5 mmg.log; exit 1; }
      N=$(awk '/^ *Tetrahedra/ {getline; print $1; exit}' out.mesh)
      log "malla, pasada $it: $N celdas (objetivo $TARGET)"
      awk -v n="$N" -v t="$TARGET" 'BEGIN {r = n / t; exit !(r > 0.93 && r < 1.07)}' && break
      F=$(awk -v f="$F" -v n="$N" -v t="$TARGET" 'BEGIN {print f * n / t}')
    done
    python3 "$R" clean in=out.mesh out=Fine2.mesh || exit 1
    rm -f out.mesh out.sol
  fi
fi
if [ ! -f test.ok ]; then
  log "prueba: TVE, 20 pasos"
  nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -lms 500 > mem.txt &
  MON=$!
  # shellcheck disable=SC2086
  "$BIN" Fine2.mesh fluid=TVE t_end=0.02 output_every=100000 $EXTRA output=test > test.log 2>&1
  ok=$?
  kill $MON 2> /dev/null
  log "pico GPU $(sort -n mem.txt | tail -1) MiB; $(grep -E '^Rheology' test.log)"
  [ $ok -eq 0 ] || { tail -3 test.log; log "la prueba falló: baja TARGET y borra Fine2.mesh"; exit 1; }
  touch test.ok
fi

run() {
  local n=$1; shift
  [ -f "$n.done" ] && return
  log "$n"
  # shellcheck disable=SC2086
  if "$BIN" Fine2.mesh "$@" $EXTRA output="$n" > "$n.log" 2>&1; then touch "$n.done"; else log "$n FALLÓ"; fi
  grep -E '^Rheology|Diverged' "$n.log" | tail -1 | cut -c1-160
}

# 2. Arranque newtoniano (solo escribe el punto de reinicio).
run spin_N fluid=N cycles="$SPIN" output_every=1000000
[ -f spin_N.restart ] || { log "no hay spin_N.restart"; exit 1; }
log "energía cinética media por ciclo del arranque (mJ): $(awk -F, 'NR > 1 {e[$2] += $11; n[$2]++} END {for (c = 1; c in n; c++) printf "%.3f ", 1e3 * e[c] / n[c]}' spin_N.csv)"

# 3. Producción desde el mismo estado; TVE y sus controles primero.
run prod_TVE fluid=TVE cycles="$CYCLES" restart=spin_N.restart output_every=50
run ctl_TVE_devss05 fluid=TVE cycles="$CTL" devss=0.5 restart=spin_N.restart output_every=50
run ctl_TVE_imp1 fluid=TVE cycles="$CTL" implicit=1 restart=spin_N.restart output_every=50
python3 "$RUNS/stats_models.py" runs=TVE:prod_TVE,devss05:ctl_TVE_devss05,imp1:ctl_TVE_imp1 \
  pairs=devss05/TVE,imp1/TVE skip=1 out=controles > controles.log 2>&1 && log "controles: $DIR/controles_summary.txt"
for f in N CY VE; do run "prod_$f" fluid=$f cycles="$CYCLES" restart=spin_N.restart output_every=50; done

# 4. Estadística.
python3 "$RUNS/stats_models.py" skip=2 ${REGION:+region=$REGION} out=final > stats.log 2>&1 \
  && log "estadística: $DIR/final_summary.txt"
log "fin"
