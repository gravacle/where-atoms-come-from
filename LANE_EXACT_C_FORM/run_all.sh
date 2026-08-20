#!/bin/zsh
# LANE_EXACT_C_FORM -- run order.  t0 must pass before anything else is believed.
cd "/Users/bgm/MB Work/where-atoms-come-from/LANE_EXACT_C_FORM"
for s in t0_clauses t1_additivity t2_separation t4_superposition t5_range t6_scaling t7_assumptions t8_weakfield; do
  python3 -u $s.py > $s.txt 2>&1
  echo "$s: $(grep -c '\[FAIL\]' $s.txt) failures"
done
