#!/bin/bash
# Campaña del atrio en GPU:
#   1. recompila el driver si hace falta y genera la malla adaptada Atrium.mesh
#      (remallando IN_MESH guiado por la corrida GUIDE) si no existe;
#   2. arranque: N durante SPIN ciclos, hasta un flujo desarrollado;
#   3. desde ese mismo estado, TVE y sus controles numéricos (devss=0.5,
#      implicit=1; CTL ciclos), y N, CY y VE, durante CYCLES ciclos;
#   4. estadística pareada ciclo a ciclo.
# Cada corrida escribe un punto de reinicio cada CK pasos. Si se cae (p. ej. un
# Xid de la GPU), se retoma sola desde el último, hasta MAXTRY veces, como
# <nombre>_r1, _r2...; <nombre>.parts lista las partes y stats_models.py las lee
# juntas. Se puede relanzar: salta lo hecho.
# Variables: DIR ($HOME/runs/campaign), TARGET (1200000 celdas), SIZE_ARGS
# ("wall=0.85 cap=1.0": tamaño parejo, centro casi tan fino como la pared),
# SPIN (8), CYCLES (10), CTL (3), CK (170), MAXTRY (20), REGION ("x,y,z,r" mm).
set -uo pipefail
source "$HOME/gpu/env.sh"
export LC_ALL=C
RUNS=${RUNS:-$HOME/runs}
DIR=${DIR:-$RUNS/campaign}
SPIN=${SPIN:-8}
CYCLES=${CYCLES:-10}
CTL=${CTL:-3}
CK=${CK:-170}
MAXTRY=${MAXTRY:-20}
TARGET=${TARGET:-1200000}
SIZE_ARGS=${SIZE_ARGS:-wall=0.85 cap=1.0}
REGION=${REGION:-}
IN_MESH=${IN_MESH:-$RUNS/fine/Fine.mesh}   # malla que se remalla
GUIDE=${GUIDE:-$RUNS/fine/fine_TVE}        # corrida que guía el refinamiento
EXTRA=${EXTRA:-}                           # claves extra para todas las corridas
PERIOD=${PERIOD:-0.85}
STEPS=$(echo "$EXTRA" | sed -n 's/.*steps=\([0-9]*\).*/\1/p'); STEPS=${STEPS:-850}
R=$RUNS/refine_atrium.py
MMG=$HOME/gpu/mmg/build/bin/mmg3d_O3
BIN=$RODIN_BIN/KokkosAtriumRheologyPETSc
MESH=Atrium.mesh
mkdir -p "$DIR" && cd "$DIR" || exit 1
log() { echo "[$(date '+%d %H:%M')] $*"; }
ckstep() { od -An -t u8 -j 16 -N 8 "$1" | tr -d ' '; }   # paso guardado en un punto de reinicio

# 0. Detener lo que esté corriendo.
pkill -f run10.sh; pkill -f run_fine.sh; pkill -f run_fine2.sh
pkill -f KokkosAtriumRheologyPETSc; while pgrep -f KokkosAtriumRheologyPETSc > /dev/null; do sleep 5; done
log "GPU libre; $(df -BG --output=avail . | tail -1 | tr -dc 0-9) GB libres; potencia máx. $(nvidia-smi --query-gpu=power.limit --format=csv,noheader)"

# 1. Driver y malla.
if ! grep -aq 'checkpoint_every' "$BIN"; then
  log "recompilando Rodin desde la rama model/AtriumKokkos"
  SKIP_TEST=1 bash "$HOME/install_rodin_gpu.sh" build > build.log 2>&1 || { log "la compilación falló: ver $DIR/build.log"; exit 1; }
  source "$HOME/gpu/env.sh"
  grep -aq 'checkpoint_every' "$BIN" || { log "el driver no tiene checkpoint_every: ¿se hizo git push de la rama?"; exit 1; }
fi
log "driver: ok"
if [ ! -f "$MESH" ]; then
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
  python3 "$R" clean in=out.mesh out="$MESH" || exit 1
  rm -f out.mesh out.sol test.ok
fi
if [ ! -f test.ok ]; then
  log "prueba: TVE, 20 pasos"
  nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -lms 500 > mem.txt &
  MON=$!
  # shellcheck disable=SC2086
  "$BIN" "$MESH" fluid=TVE t_end=0.02 output_every=100000 $EXTRA output=test > test.log 2>&1
  ok=$?
  kill $MON 2> /dev/null
  log "pico GPU $(sort -n mem.txt | tail -1) MiB; $(grep -E '^Rheology' test.log)"
  [ $ok -eq 0 ] || { tail -3 test.log; log "la prueba falló"; exit 1; }
  touch test.ok
fi

run() {
  local n=$1; shift
  [ -f "$n.done" ] && return 0
  local base=() cycles=0 from="" at=0 a
  for a in "$@"; do
    case $a in
      cycles=*) cycles=${a#cycles=} ;;
      restart=*) from=${a#restart=}; at=$(ckstep "$from") ;;
      *) base+=("$a") ;;
    esac
  done
  local end=$((at + cycles * STEPS)) part=$n try=0 limit=("cycles=$cycles")
  echo "$n" > "$n.parts"
  while :; do
    log "$part (pasos $((at + 1))-$end)"
    # shellcheck disable=SC2086
    if "$BIN" "$MESH" "${base[@]}" ${from:+restart=$from} "${limit[@]}" checkpoint_every="$CK" \
         $EXTRA output="$part" > "$part.log" 2>&1; then
      touch "$n.done"
      grep -E '^Rheology' "$part.log" | tail -1 | cut -c1-160
      return 0
    fi
    log "$part FALLÓ: $(grep -m1 -E 'Diverged|did not converge|PETSC ERROR: [A-Z]|what\(\)' "$part.log" | cut -c1-150)"
    sudo -n dmesg -T 2> /dev/null | grep -i xid | tail -1
    try=$((try + 1))
    [ "$try" -le "$MAXTRY" ] || { log "$n: sin más reintentos"; return 1; }
    if [ -f "$part.restart" ]; then from=$part.restart; at=$(ckstep "$from"); fi
    local rem=$((end - at))
    [ "$rem" -gt 0 ] || { touch "$n.done"; return 0; }
    limit=("cycles=$(((rem + STEPS - 1) / STEPS))" "t_end=$(awk -v r="$rem" -v p="$PERIOD" -v s="$STEPS" 'BEGIN {printf "%.9f", r * p / s}')")
    part=${n}_r$try
    echo "$part" >> "$n.parts"
    sleep 30
  done
}
first() { local n=$1; shift; run "$n" "$@" cycles="$CYCLES"; }

# 2. Arranque newtoniano (solo escribe el punto de reinicio).
run spin_N fluid=N cycles="$SPIN" output_every=1000000 || exit 1
SPIN_RESTART=$(tail -1 spin_N.parts).restart
[ "$(ckstep "$SPIN_RESTART")" -eq $((SPIN * STEPS)) ] || { log "el arranque no llegó al final"; exit 1; }
cp "$SPIN_RESTART" spin.restart
log "energía cinética media por ciclo del arranque (mJ): $(sed 's/$/.csv/' spin_N.parts | xargs cat | awk -F, 'NR == 1 {for (i = 1; i <= NF; i++) if ($i == "kineticEnergy") k = i} $1 != "t" {t[$1] = $2 SUBSEP $k} END {for (x in t) {split(t[x], v, SUBSEP); e[v[1]] += v[2]; n[v[1]]++} for (c = 1; c in n; c++) printf "%.3f ", 1e3 * e[c] / n[c]}')"

# 3. Producción desde el mismo estado; TVE y sus controles primero.
first prod_TVE fluid=TVE restart=spin.restart output_every=50
run ctl_TVE_devss05 fluid=TVE cycles="$CTL" devss=0.5 restart=spin.restart output_every=50
run ctl_TVE_imp1 fluid=TVE cycles="$CTL" implicit=1 restart=spin.restart output_every=50
python3 "$RUNS/stats_models.py" runs=TVE:prod_TVE,devss05:ctl_TVE_devss05,imp1:ctl_TVE_imp1 \
  pairs=devss05/TVE,imp1/TVE skip=1 out=controles > controles.log 2>&1 && log "controles: $DIR/controles_summary.txt"
for f in N CY VE; do first "prod_$f" fluid=$f restart=spin.restart output_every=50; done

# 4. Estadística.
python3 "$RUNS/stats_models.py" skip=2 ${REGION:+region=$REGION} out=final > stats.log 2>&1 \
  && log "estadística: $DIR/final_summary.txt"
log "fin"
