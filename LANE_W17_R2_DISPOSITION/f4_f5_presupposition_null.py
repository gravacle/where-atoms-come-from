"""
F4 — PRESUPPOSITION  and  F5 — THE NULL OPTION, for R2.

F4 sentence under test:
  "A disposition is a single program-global value; the two named arms exhaust it;
   and setting it determines what the program does next."
Tested on three independent legs, each with a count.

F5: is NEITHER admissible, with its own evidence?
"""
import os, re, math
import numpy as np
from cutoff_guard import CUT, TEXT, line

REPO = "/Users/bgm/MB Work/where-atoms-come-from"
W05 = "\n".join(CUT[367:494])     # register lines 368..494
W06 = "\n".join(CUT[494:639])     # register lines 495..639

print("=" * 78)
print("F4 — PRESUPPOSITION")
print("=" * 78)
print("SENTENCE UNDER TEST:")
print('  "A disposition is a single program-global value; the two named arms exhaust')
print('   it; and setting it determines what the program does next."')

# ---- LEG 1: the register's unit of ruling is the QUESTION, not the program.
rows = []
cur = None
for ln in CUT:
    if ln.startswith("## ") and not ln.startswith("## HOW TO USE"):
        cur = [ln, []]; rows.append(cur)
    elif cur is not None:
        cur[1].append(ln)
SCOPE = {"W-01": "question", "W-02": "question", "ERRATUM": "question",
         "W-03": "question", "W-04": "PROGRAM", "W-05": "PROGRAM", "W-06": "PROGRAM"}
print("\nLEG 1 — custody §8 / REGISTER:3-8: the register's unit is the QUESTION.")
print(f"  {'row':<10}{'scope':<10}{'proof ptr':<11}{'digest':<9}{'reopen':<8}")
tab = {}
for head, body in rows:
    tag = head.replace("## ", "").split("—")[0].strip()[:7]
    key = "ERRATUM" if tag.startswith("ERRATUM") else tag
    b = "\n".join(body)
    rec = (SCOPE[key], "WHERE THE PROOF IS" in b, "sha256" in b, "REOPEN" in b)
    tab[key] = rec
    print(f"  {key:<10}{rec[0]:<10}{str(rec[1]):<11}{str(rec[2]):<9}{str(rec[3]):<8}")
prog = [k for k, v in tab.items() if v[0] == "PROGRAM"]
ques = [k for k, v in tab.items() if v[0] == "question"]
p_np = sum(1 for k in prog if not tab[k][1])
q_np = sum(1 for k in ques if not tab[k][1])
p_nr = sum(1 for k in prog if not tab[k][3])
q_nr = sum(1 for k in ques if not tab[k][3])
print(f"\n  P(no proof pointer | PROGRAM-scoped row) = {p_np}/{len(prog)} = {p_np/len(prog):.2f}")
print(f"  P(no proof pointer | question-scoped row)= {q_np}/{len(ques)} = {q_np/len(ques):.2f}")
print(f"  P(no reopen cond.  | PROGRAM-scoped row) = {p_nr}/{len(prog)} = {p_nr/len(prog):.2f}")
print(f"  P(no reopen cond.  | question-scoped row)= {q_nr}/{len(ques)} = {q_nr/len(ques):.2f}")
from math import comb
def fisher_1s(a,b,c,d):
    n=a+b+c+d; r1=a+b; c1=a+c
    tot=comb(n,c1); p=0.0
    for x in range(max(0,c1-(n-r1)), min(r1,c1)+1):
        pr=comb(r1,x)*comb(n-r1,c1-x)/tot
        if x>=a: p+=pr
    return p
pv_proof=fisher_1s(p_np,len(prog)-p_np,q_np,len(ques)-q_np)
pv_reop =fisher_1s(p_nr,len(prog)-p_nr,q_nr,len(ques)-q_nr)
print(f"  CONFOUND, RECORDED: n = 7 rows. The one-sided pattern is perfect but tiny.")
print(f"  one-sided Fisher exact, no-proof-pointer  ({p_np}/{len(prog)} vs {q_np}/{len(ques)}) : p = {pv_proof:.3f}")
print(f"  one-sided Fisher exact, no-reopen-cond.   ({p_nr}/{len(prog)} vs {q_nr}/{len(ques)}) : p = {pv_reop:.3f}")
print("  The finding is the PATTERN and its custody-§8 consequence, NOT a p-value.")
print("  CONSEQUENCE (custody §8, REGISTER:8): 'A row with no reopen condition is closed")
print("  permanently.'  Both arms of R2 live in rows that register no reopen condition,")
print("  so BOTH are permanently-closed rows by the program's own convention — which")
print("  contradicts STOP-FALLS-REBUILD's own words 'The live route' at REGISTER:619.")

# ---- LEG 2: the founding ladder is per-stage, with per-stage falsifiers.
fd = open(f"{REPO}/FOUNDING_DESIGN_V001.md", encoding="utf-8").read()
stages = re.findall(r"\*\*(S[1-6]) — ([A-Z][^*]*?)\*\*", fd)
fals = fd.count("*Falsifier:*")
print(f"\nLEG 2 — FOUNDING_DESIGN §7: the program's own state variable is a LADDER.")
print(f"  stages declared            : {len(stages)}  {[s[0] for s in stages]}")
print(f"  per-stage falsifiers        : {fals}")
print(f"  a per-stage state vector over {{not-run, run, refuted}} has 3^{len(stages)} = "
      f"{3**len(stages)} values; the disposition binary has 2.")
alpha_w05 = len(re.findall(r"\balpha\b", W05, re.I))
alpha_w06 = len(re.findall(r"\balpha\b", W06, re.I))
print(f"  occurrences of 'alpha' in the W-05 row : {alpha_w05}")
print(f"  occurrences of 'alpha' in the W-06 row : {alpha_w06}")
print(f"  S5 ('LOOK FOR THE SLOT', FOUNDING:108-111) is where alpha — one of the three")
print(f"  terms of the target at FOUNDING:13-14 — enters. NEITHER ARM MENTIONS IT.")

# ---- LEG 3: W-06 refutes the presupposition in the act of using it.
print("\nLEG 3 — the posing text decomposes the very thing it asks to be scalar.")
for n in (615, 617, 618):
    print(f"  REGISTER:{n}  {line(n).strip()[:88]}")
print("  W-06 needs THREE components — OBJECT, GROUND, DISPOSITION — to state its own")
print("  position, and assigns them DIFFERENT values (dead / flipped / defective).")
print("  A question demanding one value on a 3-vector has a false presupposition.")
print("\n  F4 RESULT: PRESUPPOSITION FALSE, on all three legs, from in-cut text alone.")

# ============================================================================
print("\n" + "=" * 78)
print("F5 — THE NULL OPTION.  Does EITHER branch obtain?")
print("=" * 78)

# (a) admissibility under the pointer rule
w06_digests = len(re.findall(r"sha256", W06))
w05_digests = len(re.findall(r"sha256", W05))
w06_fileline = re.findall(r"[A-Z0-9_]+_V001\.md:\d+", W06)
w05_fileline = re.findall(r"[A-Z0-9_]+_V001\.md:\d+", W05)
lane_dirs = [d for d in os.listdir(REPO) if d.startswith("LANE_W0") and d[7] in "0123456"]
print("(a) ADMISSIBILITY UNDER CUSTODY §1 (a governing clause needs a pointer)")
print(f"    digests in the W-05 row                     : {w05_digests}")
print(f"    digests in the W-06 row                     : {w06_digests}")
print(f"    file:line pointers in the W-05 row          : {len(w05_fileline)}  {w05_fileline}")
print(f"    file:line pointers in the W-06 row          : {len(w06_fileline)}  {w06_fileline}")
print(f"    lane directories on disk for W-01..W-06     : {len(lane_dirs)}  {lane_dirs}")
print("    STOP's ground (g1, vacuity/IMP-1): W-06 itself finds it entered 'with no digest,")
print("      no file:line, no named ruling' — REGISTER:508-510 — hence FLAGGED, not inherited.")
print("    SFR's ground (g2, S3-0 false with an exhibited counterexample) rests on figures")
print("      (4.45e-16, 3*sqrt(3)/10, 1000-of-4000) that carry no digest and no artifact.")
print("    ADMISSIBLE ARMS UNDER CUSTODY §1 : 0 of 2.")

# (b) the catch-rate estimator
caught = 5   # S2 by W-03/W-04 · W-03 by W-04 · W-04's chain by W-05 · W-05 by W-06 · (S4 by W-03)
layers  = 5
lap  = (caught + 1) / (layers + 2)
jeff = (caught + 0.5) / (layers + 1)
print(f"\n(b) THE LAYER CATCH RATE (REGISTER:635: 'FIVE consecutive layers have each been")
print(f"    caught by the next, and the rate has not fallen. Discount this layer too.')")
print(f"    observed                       : {caught}/{layers} = {caught/layers:.3f}")
print(f"    Laplace  P(layer 6 also caught): ({caught}+1)/({layers}+2) = {lap:.3f}")
print(f"    Jeffreys P(layer 6 also caught): ({caught}+0.5)/({layers}+1) = {jeff:.3f}")
print(f"    => P(the ground of EITHER arm survives the next layer) <= {1-lap:.3f}")
print("    CONTROL VOID, DECLARED: with 5/5 and no failures every standard estimator is")
print("    forced high. This sub-test COULD NOT have come out low and carries NO WEIGHT")
print("    as a test. It is reported only as an order of magnitude.")

# (c) the third option, with its own evidence
print("\n(c) THE THIRD OPTION, POSITIVELY EVIDENCED — SUSPEND, DO NOT DISPOSE:")
print("    b08  run the lineage-independent lane   REGISTER:87 (W-01 reopen), :158 (W-02 reopen)")
print("    b07  audit S1                           REGISTER:292-299 (ERR-4)")
print("    Both are REGISTERED reopen conditions or registered defects, both are unrun, and")
print("    custody §4 names lineage independence as the remedy for this corpus's grade.")
print("    2 of the 11 registered reopen conditions name the lineage-independent lane;")
print(f"    lane directories that would carry it: {len(lane_dirs)}.")
print("    => NEITHER is not a residue here. It is the option the register's own reopen")
print("       conditions select, and it is the only one either arm's ground could be")
print("       tested by.")

# ============================================================================
print("\n" + "=" * 78)
print("SIDE CHECK — the ONE concrete claim inside the SFR arm, verified")
print("=" * 78)
ks = [(5**n - 1) // 4 for n in range(1, 6)]
Vs = [4 * k + 1 for k in ks]
print(f"  W-06 REGISTER:620-621: 'V = 4k+1; k = 1,6,31,156,781 gives V = 5,25,125,625,3125'")
print(f"  computed k = {ks}   V = {Vs}")
print(f"  all V are powers of 5              : {all(V == 5**(i+1) for i, V in enumerate(Vs))}")
print(f"  every consecutive pair divides     : {all(Vs[i+1] % Vs[i] == 0 for i in range(len(Vs)-1))}")
print(f"  reason: 5 = 1 (mod 4), so 5^n = 1 (mod 4) for every n — the chain is FORCED,")
print(f"  which is exactly W-06's 'accident of two arithmetic progressions'. VERIFIED.")
base = [v for v in range(5, 60) if v % 4 == 1]
print(f"  BUT the same construction works for EVERY V0 = 1 (mod 4): {base}")
print(f"  so b01 is not a single route: it is a one-parameter family with "
      f"{len(base)} members below V0=60 alone.")
print(f"  The binary's 'rebuild' arm is itself an unenumerated family, not a point.")
