#!/bin/zsh
cd "/Users/bgm/MB Work/where-atoms-come-from/LANE_EXACT_C_FORM"
while [ ! -f run_all4.done ]; do sleep 10; done
python3 -u t7_assumptions.py > t7_assumptions.txt 2>&1
echo DONE > run_all5.done
