#!/bin/zsh
cd "/Users/bgm/MB Work/where-atoms-come-from/LANE_EXACT_C_FORM"
while [ ! -f run_rest.done ]; do sleep 5; done
python3 -u t6_scaling.py > t6_scaling.txt 2>&1
echo T6DONE > run_t6.done
