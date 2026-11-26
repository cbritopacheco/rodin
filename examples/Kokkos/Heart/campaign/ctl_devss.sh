#!/bin/bash
# Control DEVSS para la comparación VE/TVE: ¿la menor variabilidad y la mayor
# estasis de TVE frente a VE sobreviven con devss=0.5?
#   - alarga ctl_TVE_devss05 (3 ciclos) hasta los mismos 10 ciclos de prod_*;
#   - corre VE con devss=0.5 desde el mismo arranque, 10 ciclos;
#   - estadística pareada en los ciclos 11-18, junto con las corridas de producción.
# Retoma sola tras una caída, como run_final.sh. Se puede relanzar: salta lo hecho.
set -uo pipefail
source "$HOME/gpu/env.sh"
export LC_ALL=C
RUNS=${RUNS:-$HOME/runs}
DIR=${DIR:-$RUNS/campaign}
CYCLES=${CYCLES:-10}
CK=${CK:-170}
MAXTRY=${MAXTRY:-20}
EXTRA=${EXTRA:-}                           # claves extra, como en run_final.sh
PERIOD=${PERIOD:-0.85}
STEPS=$(echo "$EXTRA" | sed -n 's/.*steps=\([0-9]*\).*/\1/p'); STEPS=${STEPS:-850}
BIN=$RODIN_BIN/KokkosAtriumRheologyPETSc
MESH=Atrium.mesh
cd "$DIR" || exit 1
log() { echo "[$(date '+%d %H:%M')] $*"; }
ckstep() { od -An -t u8 -j 16 -N 8 "$1" | tr -d ' '; }

pgrep -f "$BIN" > /dev/null && { log "ya hay una corrida en la GPU"; exit 1; }

run() {   # igual que en run_final.sh
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

S0=$(ckstep spin.restart)
END=$((S0 + CYCLES * STEPS))

# 1. TVE con devss=0.5: de donde quedó el control hasta el paso END.
FROM=$(tail -1 ctl_TVE_devss05.parts).restart
AT=$(ckstep "$FROM")
[ $(((END - AT) % STEPS)) -eq 0 ] || { log "$FROM está en el paso $AT, fuera de un fin de ciclo"; exit 1; }
run ctl_TVE_devss05_ext fluid=TVE devss=0.5 restart="$FROM" cycles=$(((END - AT) / STEPS)) output_every=50 || exit 1
cat ctl_TVE_devss05.parts ctl_TVE_devss05_ext.parts > ctl_TVE_d05.parts

# 2. VE con devss=0.5 desde el arranque.
run ctl_VE_devss05 fluid=VE devss=0.5 restart=spin.restart cycles="$CYCLES" output_every=50 || exit 1

# 3. Estadística: producción (devss=1) y control (devss=0.5) juntos.
python3 "$RUNS/stats_models.py" \
  runs=N:prod_N,CY:prod_CY,VE:prod_VE,TVE:prod_TVE,VEd:ctl_VE_devss05,TVEd:ctl_TVE_d05 \
  pairs=TVE/VE,TVEd/VEd,VEd/VE,TVEd/TVE,TVE/CY,TVEd/CY skip=2 out=devss > devss.log 2>&1 \
  && log "estadística: $DIR/devss_summary.txt" || log "la estadística falló: ver $DIR/devss.log"
log "fin"
