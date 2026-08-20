#!/bin/zsh
cd "/Users/bgm/MB Work/where-atoms-come-from/LANE_EXACT_C_FORM"
python3 -u t5_range.py   > t5_range.txt 2>&1
python3 -u t6_scaling.py > t6_scaling.txt 2>&1
python3 -u t7_assumptions.py > t7_assumptions.txt 2>&1
python3 -u t2_separation.py > t2_separation.txt 2>&1
python3 -u t4_superposition.py > t4_superposition.txt 2>&1
echo DONE > run_all4.done
