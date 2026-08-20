#!/bin/zsh
cd "/Users/bgm/MB Work/where-atoms-come-from/LANE_EXACT_C_FORM"
while pgrep -f t1_additivity.py > /dev/null; do sleep 5; done
python3 -u t2_separation.py > t2_separation.txt 2>&1
python3 -u t5_range.py > t5_range.txt 2>&1
python3 -u t4_superposition.py > t4_superposition.txt 2>&1
echo ALLDONE > run_rest.done
