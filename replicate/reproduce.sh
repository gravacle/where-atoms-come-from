#!/usr/bin/env bash
# ONE COMMAND THAT REPRODUCES THE WHOLE CORPUS.
#   ./replicate/reproduce.sh          re-run every lane script, diff against the SEALED output
#   ./replicate/reproduce.sh --seals  verify SHA-256 seals only (fast)
#   ./replicate/reproduce.sh --quick  skip the scripts known to take minutes
# Exits non-zero if any seal fails or any script's output differs from its sealed .txt.
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
FAIL=0

hdr(){ printf '\n\033[1m%s\033[0m\n%s\n' "$1" "$(printf '=%.0s' $(seq 1 78))"; }

hdr "0. ENVIRONMENT"
printf "  python  %s\n  numpy   %s\n" "$(python3 -V 2>&1)" "$(python3 -c 'import numpy;print(numpy.__version__)' 2>/dev/null || echo MISSING)"
python3 -c 'import numpy' 2>/dev/null || { echo "  numpy is required"; exit 2; }

hdr "1. SEALS  (SHA-256 over every registered document and lane)"
bad=0
for f in *.sha256; do
  case "$f" in *.seal.sha256) continue;; esac
  if ! shasum -a 256 -c "$f" >/dev/null 2>&1; then echo "  FAIL  $f"; bad=$((bad+1)); fi
done
[ "$bad" -eq 0 ] && echo "  all $(ls *.sha256 | grep -vc '\.seal\.' ) manifests verify" || FAIL=1

if [ "${1:-}" = "--seals" ]; then
  hdr "RESULT"; [ "$FAIL" -eq 0 ] && echo "  SEALS OK" || echo "  SEAL FAILURE"; exit "$FAIL"
fi

SLOW="o1_gap_a_compatibility.py o1_converse_theorem.py o4_L3.py o4_toric.py verify_o4.py o10_arbitrary_coupling.py"
hdr "2. RE-RUN EVERY LANE SCRIPT AND DIFF AGAINST ITS SEALED OUTPUT"
printf "  %-30s %-26s %s\n" "lane" "script" "result"
printf "  %-30s %-26s %s\n" "------------------------------" "--------------------------" "------"
for py in LANE_*/*.py; do
  lane="${py%%/*}"; base="$(basename "$py")"; txt="${py%.py}.txt"
  if [ ! -f "$txt" ]; then printf "  %-30s %-26s %s\n" "$lane" "$base" "no sealed output — SKIP"; continue; fi
  if [ "${1:-}" = "--quick" ] && [[ " $SLOW " == *" $base "* ]]; then
    printf "  %-30s %-26s %s\n" "$lane" "$base" "skipped (--quick)"; continue; fi
  out="$TMP/$lane.$base.out"
  ( cd "$lane" && python3 -u "$base" ) >"$out" 2>&1
  if diff -q "$txt" "$out" >/dev/null 2>&1; then
    printf "  %-30s %-26s \033[32m%s\033[0m\n" "$lane" "$base" "IDENTICAL"
  else
    n=$(diff "$txt" "$out" | grep -c '^[<>]' || true)
    printf "  %-30s %-26s \033[31m%s\033[0m\n" "$lane" "$base" "DIFFERS ($n lines)"
    FAIL=1
  fi
done

hdr "3. THE STATUS LEDGER IS GENERATED, NOT TYPED"
before=$(shasum -a 256 STATUS_LEDGER_V001.md | cut -d' ' -f1)
./ledger/status.py >/dev/null
after=$(shasum -a 256 STATUS_LEDGER_V001.md | cut -d' ' -f1)
if [ "$before" = "$after" ]; then echo "  regenerating the grid from ledger/status_ledger.tsv reproduces it byte-for-byte"
else echo "  FAIL — the committed grid does not match what the ledger generates"; FAIL=1; fi

hdr "RESULT"
[ "$FAIL" -eq 0 ] && echo "  EVERYTHING REPRODUCES" || echo "  REPRODUCTION FAILED — see above"
exit "$FAIL"
