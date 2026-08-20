#!/bin/zsh
cd "/Users/bgm/MB Work/where-atoms-come-from/LANE_EXACT_C_FORM"
while [ ! -f run_all2.done ]; do sleep 10; done
python3 -u t5_range.py > t5_range.txt 2>&1
echo DONE > run_all3.done
