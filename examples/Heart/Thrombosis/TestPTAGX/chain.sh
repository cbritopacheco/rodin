#!/bin/bash
# Runs the whole TestPTAGX battery on two cores (two chains).
export PTAGP_MPIRUN="mpirun --allow-run-as-root -n 1"
HERE=$(cd "$(dirname "$0")" && pwd)
EXE=${PTAGX_EXE_FROZEN:-$HERE/../../../../build/examples/Heart/Thrombosis/LeftAtrium2DPTAGX/LeftAtrium2DPTAGX}
chainA() {
  while pgrep -f "E3_DOAC/run_E3.py" > /dev/null; do sleep 20; done
  python3 $HERE/E3_DOAC/run_E3.py --jobs 1 > $HERE/E3_DOAC/run.log 2>&1; echo "STAGE E3 done"
  python3 $HERE/E4_PatchThreshold/run_E4.py --jobs 1 > $HERE/E4_PatchThreshold/run.log 2>&1; echo "STAGE E4 done"
  while [ ! -f $HERE/Q1_Qureshi2021/results/calibration.json ]; do sleep 30; done
  python3 $HERE/Q1_Qureshi2021/run_Q1.py --stage coag --jobs 1 --exe-x $EXE --kinetics ptagx \
    --shapes CW,BR > $HERE/Q1_Qureshi2021/coag_ptagx.log 2>&1; echo "STAGE Q1 ptagx done"
}
chainB() {
  while pgrep -f "run_Q1.py --stage calib" > /dev/null; do sleep 20; done
  python3 $HERE/Q1_Qureshi2021/run_Q1.py --stage calib --jobs 1 --exe-x $EXE > $HERE/Q1_Qureshi2021/calib.log 2>&1
  echo "STAGE Q1 calib done"
  python3 $HERE/Q1_Qureshi2021/run_Q1.py --stage coag --jobs 1 --exe-x $EXE --kinetics qureshi \
    > $HERE/Q1_Qureshi2021/coag_qureshi.log 2>&1; echo "STAGE Q1 qureshi done"
}
chainA & chainB & wait
python3 $HERE/Q1_Qureshi2021/run_Q1.py --stage report --exe-x $EXE > $HERE/Q1_Qureshi2021/report.log 2>&1
echo "ALL DONE"
