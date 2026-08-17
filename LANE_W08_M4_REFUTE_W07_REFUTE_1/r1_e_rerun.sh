#!/bin/sh
# W-08 / M4-REFUTE-1 leg E — REPRODUCE BOTH LANES BEFORE JUDGING EITHER.
W07="/Users/bgm/MB Work/where-atoms-come-from/LANE_W07_RECURRENCE_ISOLATION"
M4="/Users/bgm/MB Work/where-atoms-come-from/LANE_W08_M4_REFUTE_W07"
echo "== SEALS =="
cd "$W07" && printf "  W-07 : "; shasum -a 256 -c SEALS.sha256 2>&1 | grep -c OK
cd "$M4"  && printf "  M4   : "; shasum -a 256 -c SEALS.sha256 2>&1 | grep -c OK
cd "$M4"  && shasum -a 256 -c SEALS.sha256 2>&1 | grep -v ': OK' | sed 's/^/  M4 FAILED: /'
echo
echo "== RE-RUN M4's OWN PYTHON LEGS AND DIFF AGAINST ITS SEALED OUTPUT =="
for f in m4_b_rho_algebra m4_c_diophantine m4_d_forced_agreement m4_e_3sqrt3 m4_f_sweep m4_g_attained_exact m4_h_corrections; do
  python3 "$M4/$f.py" > "/tmp/r1_$f.txt" 2>&1
  if diff -q "/tmp/r1_$f.txt" "$M4/$f.OUT.txt" >/dev/null 2>&1; then
    echo "  $f.py : IDENTICAL to sealed output"
  else
    echo "  $f.py : DIFFERS ---------------"
    diff "$M4/$f.OUT.txt" "/tmp/r1_$f.txt" | head -30 | sed 's/^/      /'
  fi
done
echo
echo "== RE-RUN W-07's PYTHON LEGS =="
for f in w07_a_carrier w07_b_dressed w07_c_scaling w07_c2_exact w07_d_carrier_recur w07_e_isolation; do
  python3 "$W07/$f.py" > "/tmp/r1_$f.txt" 2>&1
  if diff -q "/tmp/r1_$f.txt" "$W07/$f.OUT.txt" >/dev/null 2>&1; then
    echo "  $f.py : IDENTICAL to sealed output"
  else
    echo "  $f.py : DIFFERS"; diff "$W07/$f.OUT.txt" "/tmp/r1_$f.txt" | head -20 | sed 's/^/      /'
  fi
done
echo
echo "  Platform: python3 $(python3 -c 'import sys;print(sys.version.split()[0])'), numpy $(python3 -c 'import numpy;print(numpy.__version__)'), float64 default."
