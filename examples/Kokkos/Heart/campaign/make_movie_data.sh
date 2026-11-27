#!/bin/bash
# Un ciclo más de TVE con salida densa para el video (atrium_movie.py).
# Parte del final de prod_TVE (mismo estado, misma carga) y escribe cada 5 pasos:
# 171 cuadros por ciclo, unos 8 GB, unos 40 min en la 5090.
# Espera a que termine la campaña (run_fa.sh) para no compartir la GPU.
# Variables: DIR ($HOME/runs/fa), EVERY (5), EXTRA ("p_in=900"), RUN (prod_TVE).
set -uo pipefail
source "$HOME/gpu/env.sh"
DIR=${DIR:-$HOME/runs/fa}
EVERY=${EVERY:-5}
EXTRA=${EXTRA:-p_in=900}
RUN=${RUN:-prod_TVE}
BIN=$RODIN_BIN/KokkosAtriumRheologyPETSc
cd "$DIR" || exit 1
log() { echo "[$(date '+%d %H:%M')] $*"; }

while pgrep -f run_fa.sh > /dev/null || pgrep -f "$BIN" > /dev/null; do
  log "la GPU está ocupada; espero"
  sleep 600
done
FROM=$(tail -1 "$RUN.parts").restart
[ -f "$FROM" ] || { log "no encuentro $FROM"; exit 1; }
log "movie_TVE desde $FROM, salida cada $EVERY pasos"
# shellcheck disable=SC2086
if "$BIN" Atrium.mesh fluid=TVE $EXTRA restart="$FROM" cycles=1 output_every="$EVERY" \
     output=movie_TVE > movie_TVE.log 2>&1; then
  log "listo: $(ls movie_TVE.velocity.*.h5 | wc -l) cuadros; ahora atrium_movie.py"
else
  log "falló: $(grep -m1 -E 'PETSC ERROR|what\(\)|did not converge' movie_TVE.log)"
  exit 1
fi
