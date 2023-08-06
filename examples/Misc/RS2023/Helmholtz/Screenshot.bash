#!/bin/bash

mesh="$1"
gf="$2"
n="$3"

TMPFILE="$(mktemp)"
echo "# BEGIN" >> ${TMPFILE}
echo "window 0 0 600 400" >> ${TMPFILE}
echo "solution $mesh ${gf}_0000.sol" >> ${TMPFILE}
echo "{" >> ${TMPFILE}
echo "  keys R" >> ${TMPFILE}
echo "  zoom 1.5" >> ${TMPFILE}
echo "}" >> ${TMPFILE}
echo "{" >> ${TMPFILE}
for i in $(seq -f "%04g" 0 "$n")
do
  echo "  solution ${mesh} ${gf}_$i.sol screenshot ${gf}_$i.png" >> ${TMPFILE}
done
echo "}" >> ${TMPFILE}
echo "# END" >> ${TMPFILE}

cat ${TMPFILE}

echo "-- Written to: ${TMPFILE}"

# GLVIS="(glvis -run ${TMPFILE})"

