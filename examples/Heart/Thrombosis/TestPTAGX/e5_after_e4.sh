#!/bin/bash
export PTAGP_MPIRUN="mpirun --allow-run-as-root -n 1"
cd "$(dirname "$0")"
until grep -q "STAGE E4 done" chain.log; do sleep 30; done
python3 E5_ShearThreshold/run_E5.py --jobs 1 > E5_ShearThreshold/run.log 2>&1
echo "STAGE E5 done" >> chain.log
