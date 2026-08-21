#!/usr/bin/env python3
"""Render PROOF_V002 claim blocks from a claim-table JSON, filling STATUS and CARRIER MARK
   mechanically from the ledger and the T-9 audit. Typing a status by hand is how a proof goes
   stale; check_proof.py refuses a stale one, so it is never typed here either."""
import json, sys, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOT = "·"
EMD = "—"

def ledger():
    st = {}
    with open(os.path.join(ROOT, "ledger/status_ledger.tsv")) as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            p = line.rstrip("\n").split("\t")
            if p and p[0]:
                st[p[0]] = p[3] if len(p) > 3 else "?"
    return st

def plan():
    """Plan task IDs are citable too: a debt block rests on the TASK that owes the work, and the
       plan is as much a source of truth as the status ledger. A plan row is never a carrier mark."""
    pl = {}
    with open(os.path.join(ROOT, "ledger/plan.tsv")) as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            p = line.rstrip("\n").split("\t")
            if p and p[0]:
                pl[p[0]] = p[3] if len(p) > 3 else "?"
    return pl

def marks():
    m = {}
    with open(os.path.join(ROOT, "LANE_T9_AUDIT/T9_carrier_audit.tsv")) as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            p = line.rstrip("\n").split("\t")
            if p and p[0]:
                m[p[0]] = p[1]
    return m

def cell(text):
    """A literal pipe inside a cell breaks the row and the field vanishes; check_proof.py then
       reports a missing field, which is the right complaint about the wrong cause. Escape it."""
    return str(text).strip().replace("|", "\\|")


def render(claim, st, mk, pl=None):
    pl = pl or {}
    ids = claim["rows"]
    cells, has2 = [], False
    for r in ids:
        if r not in st and r in pl:
            cells.append("`%s` (%s, PLAN)" % (r, pl[r]))
            continue
        s = st.get(r, "MISSING-FROM-LEDGER")
        c = mk.get(r, "UNAUDITED")
        if c == "TWO-CARRIER":
            has2 = True
        cells.append("`%s` (%s, %s)" % (r, s, c))
    scope = claim["scope"].strip()
    pre = "SINGLE-CARRIER %s " % EMD
    if not has2 and not scope.startswith("SINGLE-CARRIER"):
        scope = pre + scope[0].lower() + scope[1:]
    if has2 and scope.startswith(pre):
        scope = scope[len(pre):]
    gate = "none" if claim["gate_name"].strip().lower() == "none" or claim["gate_file"].strip().lower() == "none" \
        else "`%s` :: `%s`" % (claim["gate_file"].strip(), claim["gate_name"].strip())
    mf = claim["model_fn"].strip()
    mf = "none" if mf.lower() == "none" else "`%s`" % mf
    return "\n".join([
        "### %s %s %s" % (claim["id"], EMD, claim["headline"]),
        "",
        claim["prose"].strip(),
        "",
        "| | |",
        "|---|---|",
        "| **model** | %s |" % cell(mf),
        "| **gate** | %s |" % gate,
        "| **grounding** | %s |" % cell(claim["grounding"]),
        "| **rows** | %s |" % (" %s " % DOT).join(cells),
        "| **scope** | %s |" % cell(scope),
        "",
    ])

if __name__ == "__main__":
    st, mk, pl = ledger(), marks(), plan()
    data = json.load(open(sys.argv[1]))
    claims = data["claims"] if isinstance(data, dict) else data
    for c in claims:
        print(render(c, st, mk, pl))
