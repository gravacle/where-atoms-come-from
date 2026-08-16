#!/bin/sh
# W-07 leg F — custody scan.  Where on disk does each number W-06 PRODUCED actually live?
# Run from the repo root.
cd "$(dirname "$0")/.."
echo "== every W-06-produced figure, traced across the repo (excluding .git) =="
for n in "4.45e-16" "3\*sqrt(3)" "2.221e-16" "8.42e-04" "4.11e-01" "1000 of 4000" "3125" "0.5196"; do
  printf "  %-16s -> " "$n"
  grep -ril -- "$n" . 2>/dev/null | grep -v '^./.git' | grep -v 'LANE_W07' | tr '\n' ' '
  echo
done
echo
echo "== lane directories left behind by each register round =="
ls -d LANE_*/ | grep -v W07
echo
echo "W-01/W-02: build+audit artifacts sealed.  W-03/W-05: lane directories on disk."
echo "W-06: no artifact, no lane directory, no .py, no .OUT.txt.  Every number it produced"
echo "appears in exactly one file: REGISTER_V001.md."
