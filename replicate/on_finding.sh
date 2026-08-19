#!/usr/bin/env bash
# THE DOCUMENTATION GATE. Run after every finding, BEFORE committing.
# It FAILS when the documentation has fallen behind the work, so the process is enforced
# rather than remembered. (D-7: a caught assumption that lives only in prose will be re-made.)
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"
FAIL=0; warn(){ echo "  FAIL  $1"; FAIL=1; }; ok(){ echo "  ok    $1"; }

echo "DOCUMENTATION GATE"; echo "=================================================================="

# 1. every lane has a sealed manifest that verifies from the repo root
for d in LANE_* model replicate; do
  [ -d "$d" ] || continue
  m="$d.sha256"; [ "$d" = "model" ] && m="model.sha256"; [ "$d" = "replicate" ] && m="replicate.sha256"
  if [ ! -f "$m" ]; then warn "$d has no sealed manifest ($m)"
  elif ! shasum -a 256 -c "$m" >/dev/null 2>&1; then warn "$m does not verify"; fi
done
[ "$FAIL" -eq 0 ] && ok "every lane is sealed and verifies from the repo root"

# 2. every lane has sealed output (three naming conventions in use across the corpus)
empty=0; thin=0
for d in LANE_*/ model/; do
  d="${d%/}"; [ -d "$d" ] || continue
  nt=$(find "$d" -name '*.txt' | wc -l | tr -d ' ')
  np=$(find "$d" -name '*.py' | wc -l | tr -d ' ')
  if [ "$nt" -eq 0 ]; then echo "        NO OUTPUT AT ALL: $d"; empty=$((empty+1))
  elif [ "$nt" -lt "$np" ]; then thin=$((thin+1)); fi
done
if [ "$empty" -eq 0 ]; then
  ok "every lane has sealed output ($thin lanes have fewer .txt than .py -- helper modules)"
else warn "$empty lanes have no output at all"; fi

# 3. the grid is GENERATED, never typed
before=$(shasum -a 256 STATUS_LEDGER_V001.md 2>/dev/null | cut -d' ' -f1)
./ledger/status.py >/dev/null 2>&1
[ "$before" = "$(shasum -a 256 STATUS_LEDGER_V001.md | cut -d' ' -f1)" ] \
  && ok "the status grid regenerates byte-for-byte from ledger/status_ledger.tsv" \
  || warn "the committed grid does not match what the ledger generates — never hand-edit it"

# 3b. THE PLAN is generated too
pb=$(shasum -a 256 THE_PLAN_V001.md 2>/dev/null | cut -d' ' -f1)
./ledger/plan.py >/dev/null 2>&1
pa=$(shasum -a 256 THE_PLAN_V001.md 2>/dev/null | cut -d' ' -f1)
if [ -z "$pb" ] || [ -z "$pa" ]; then warn "THE_PLAN could not be hashed"
elif [ "$pb" = "$pa" ]; then ok "the plan regenerates byte-for-byte from ledger/plan.tsv"
else warn "the committed plan does not match what ledger/plan.py generates"; fi

# 4. the model still validates
if ( cd model && python3 validate_model.py >/dev/null 2>&1 ); then ok "model/validate_model.py passes"
else warn "model/validate_model.py FAILS — a finding has broken the first-principles model"; fi
if ( cd model && python3 validate_formation.py >/dev/null 2>&1 ); then ok "model/validate_formation.py passes (17 formation checks)"
else warn "model/validate_formation.py FAILS -- a finding has broken the formation half"; fi
if ( cd LANE_T9_CARRIERINDEP && python3 t9_sweep.py 2>&1 | grep -q "0 FAIL" ); then ok "carrier independence holds (3 carriers, 32 checks)"
else warn "LANE_T9_CARRIERINDEP FAILS -- a finding no longer survives a second carrier"; fi
if ( cd LANE_T10_PARAMS && python3 t10_params.py >/dev/null 2>&1 ); then ok "no conclusion moves with a free parameter (28 settings)"
else warn "LANE_T10_PARAMS FAILS -- a conclusion now depends on a free parameter"; fi
if ( cd model && python3 count_law.py >/dev/null 2>&1 ); then ok "model/count_law.py passes"
else warn "model/count_law.py FAILS — the record-count law no longer holds"; fi

# 5. top-level documents are sealed
for f in REGISTER_V001.md STATUS_LEDGER_V001.md CORE_FRAMEWORK_V001.md MODEL.md REPLICATE.md; do
  [ -f "$f.sha256" ] && shasum -a 256 -c "$f.sha256" >/dev/null 2>&1 || warn "$f is unsealed or stale"
done
[ "$FAIL" -eq 0 ] && ok "anchor, register, ledger, MODEL.md and REPLICATE.md all sealed"

echo
echo "NOT MACHINE-CHECKABLE — confirm by hand before committing:"
echo "  [ ] the REGISTER has an entry for this finding, including anything WITHDRAWN by it"
echo "  [ ] the ledger row was set with ./ledger/status.py, never by editing the grid"
echo "  [ ] if a clause moved, CORE_FRAMEWORK_V001.md was amended AND re-read to confirm the patch landed"
echo "  [ ] if the model's construction or its limits changed, MODEL.md says so"
echo "  [ ] if what is claimed changed, REPLICATE.md says so"
echo "  [ ] every reported ZERO has a positive control beside it that would have registered a non-zero"
echo "  [ ] no self-check has a literal expected value, and no fit is reported without a noise floor (D-8)"
echo "  [ ] no negative verdict rests on an ensemble average alone (D-6)"
echo
[ "$FAIL" -eq 0 ] && echo "GATE PASSED" || echo "GATE FAILED"
exit "$FAIL"
