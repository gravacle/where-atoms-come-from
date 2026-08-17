#!/bin/sh
# W-08 / M4 leg A — REPRODUCE W-07 BEFORE ATTACKING IT.
# Verify the seals, re-run every W-07 script, diff against its sealed .OUT.txt.
# A refuter that cannot reproduce the thing it attacks is attacking something else.
W07="/Users/bgm/MB Work/where-atoms-come-from/LANE_W07_RECURRENCE_ISOLATION"
echo "== SEALS =="
cd "$W07" && shasum -a 256 -c SEALS.sha256 2>&1 | grep -c OK | sed 's/^/  files verifying OK: /'
shasum -a 256 -c SEALS.sha256 2>&1 | grep -v OK | sed 's/^/  FAILED: /'
echo
echo "== RE-RUN AND DIFF (empty diff = byte-identical reproduction) =="
for f in w07_a_carrier w07_b_dressed w07_c_scaling w07_c2_exact w07_d_carrier_recur w07_e_isolation; do
  python3 "$W07/$f.py" > "/tmp/m4_$f.txt" 2>&1
  if diff -q "/tmp/m4_$f.txt" "$W07/$f.OUT.txt" >/dev/null 2>&1; then
    echo "  $f.py : IDENTICAL to sealed output"
  else
    echo "  $f.py : DIFFERS ------------------"
    diff "$W07/$f.OUT.txt" "/tmp/m4_$f.txt" | sed 's/^/      /'
  fi
done
sh "$W07/w07_f_custody.sh" > /tmp/m4_w07_f.txt 2>&1
if diff -q /tmp/m4_w07_f.txt "$W07/w07_f_custody.OUT.txt" >/dev/null 2>&1; then
  echo "  w07_f_custody.sh : IDENTICAL to sealed output"
else
  echo "  w07_f_custody.sh : DIFFERS"; diff "$W07/w07_f_custody.OUT.txt" /tmp/m4_w07_f.txt | sed 's/^/      /'
fi
echo
echo "  Platform: python3 $(python3 -c 'import sys;print(sys.version.split()[0])'), numpy $(python3 -c 'import numpy;print(numpy.__version__)'), double precision throughout."
