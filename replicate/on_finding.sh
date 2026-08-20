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
# The renderer's EXIT CODE must be checked. It was previously swallowed by >/dev/null 2>&1, so a row
# written with a status outside the closed vocabulary left the grid stale and the gate still PASSED
# -- the same class of bug as D-8 (two empty hashes comparing equal and reporting success).
before=$(shasum -a 256 STATUS_LEDGER_V001.md 2>/dev/null | cut -d' ' -f1)
if ! rc=$(./ledger/status.py 2>&1 >/dev/null); then
  warn "ledger/status.py REFUSED TO RENDER: $rc"
else
  [ -n "$before" ] && [ "$before" = "$(shasum -a 256 STATUS_LEDGER_V001.md | cut -d' ' -f1)" ] \
    && ok "the status grid regenerates byte-for-byte from ledger/status_ledger.tsv" \
    || warn "the committed grid does not match what the ledger generates — never hand-edit it"
fi
# every status in the ledger must be in the closed vocabulary, checked directly
badv=$(awk -F'\t' 'NR==FNR{if(FNR>1)v[$1]=1;next} FNR>1 && !($4 in v){print $1"="$4}' \
       ledger/STATUS_VOCAB.tsv ledger/status_ledger.tsv)
[ -z "$badv" ] && ok "every status is in the closed vocabulary" \
                || warn "STATUS OUTSIDE THE CLOSED VOCABULARY: $badv"

# 3a. THE GROUNDING DEBT — printed on EVERY run, because a guard that is only a checklist item is
# a guard that changes nothing. H-3 sat OPEN for the life of the program while 162 rows were marked
# PROVED about a stipulated definition. This line makes the true state impossible not to see.
ug=$(./ledger/status.py ungrounded 2>/dev/null | head -1)
case "$ug" in
  "0 of "*) ok "$ug" ;;
  "") warn "could not compute the grounding debt" ;;
  *) warn "GROUNDING DEBT — $ug (H-3). A row about the stipulated definition is not a row about the world." ;;
esac

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
echo "NOT MACHINE-CHECKABLE — and note what that has been worth. Every guard on this list was once"
echo "an instruction from the principal. Converting an instruction into a checkbox made it LOOK"
echo "handled and changed nothing; H-3 stayed open for the program's whole life. When a guard matters,"
echo "MAKE THE TOOL REFUSE, as ./ledger/status.py now does for PROVED without a grounding."
echo ""
echo "confirm by hand before committing:"
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
