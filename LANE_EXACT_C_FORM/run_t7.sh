#!/bin/zsh
cd "/Users/bgm/MB Work/where-atoms-come-from/LANE_EXACT_C_FORM"
while [ ! -f run_t6.done ]; do sleep 5; done
python3 -u t7_assumptions.py > t7_assumptions.txt 2>&1
echo T7DONE > run_t7.done
