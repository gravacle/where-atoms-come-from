# LANE W11-R-CROSS  LEG G -- THE STEELMAN'S LOAD-BEARING CITATIONS, RESOLVED AT THE BYTES.
# Custody sec1's pointer rule: no term without a digest, a file:line, or a named ruling.  The
# steelman's whole Reading-A case rests on four pointers.  I resolve each one and quote it, and I
# record where the steelman UNDER-claimed as well as where it over-claimed.
import hashlib
import os

REPO = "/Users/bgm/MB Work/where-atoms-come-from"
FILES = ["FOUNDING_DESIGN_V001.md", "S1_CARRIER_K1_V001.md",
         "S2_FORMATION_CONDITION_ON_K1_V001.md", "S2_FORMATION_CONDITION_ON_K1_AUDIT_V001.md",
         "S3_THE_CROSSING_AUDIT_V001.md", "REGISTER_V001.md"]

print("== G0  DIGESTS OF EVERY FILE I QUOTE ==")
for f in FILES:
    p = os.path.join(REPO, f)
    h = hashlib.sha256(open(p, "rb").read()).hexdigest()
    side = p + ".sha256"
    ok = ""
    if os.path.exists(side):
        txt = open(side).read()
        ok = "SIDECAR MATCHES" if h[:16] in txt or h in txt else "SIDECAR MISMATCH"
    print(f"  {f:<48} {h[:32]}...  {ok}")


def show(f, lo, hi, tag):
    p = os.path.join(REPO, f)
    lines = open(p).read().split("\n")
    print(f"\n  --- {tag}   [{f}:{lo}-{hi}] ---")
    for i in range(lo - 1, min(hi, len(lines))):
        print(f"  {i+1:>5}| {lines[i]}")


print("\n== G1  THE PRE-REGISTERED CONTACT POINT.  RESOLVES, AND IT IS NAMED IN ADVANCE ==")
show("FOUNDING_DESIGN_V001.md", 114, 120, "sec8 EXTERNAL CONTACT, NAMED IN ADVANCE")
print("  VERDICT: the steelman's citation is EXACT.  'S2 the trivial-connection limit must give")
print("  the known trivial answer' is there, in a section headed NAMED IN ADVANCE, and it is")
print("  attached to STAGE S2 -- the stage that writes the FORMATION CONDITION.")

print("\n== G2  WHAT 'THE KNOWN TRIVIAL ANSWER' IS.  RESOLVES ==")
show("S2_FORMATION_CONDITION_ON_K1_V001.md", 575, 585, "S2 sec5 trivial-connection check")
print("  VERDICT: EXACT.  'No formation at trivial connection is the known trivial answer.'")

print("\n== G3  CHOICE LEDGER C4.  RESOLVES -- AND THE STEELMAN CONCEDED MORE THAN IT HAD TO ==")
show("S2_FORMATION_CONDITION_ON_K1_V001.md", 226, 233, "S2 sec3.2 CRITERION / C4")
print("  The steelman's self_flag says: 'C4's LETTER is about CONSTANCY in (W_F,W_C); the edge")
print("  convention is NOT constant, it merely fires with no field.  So the disqualification rests")
print("  on FOUNDING_DESIGN:117-118 directly and on C4 only derivatively.  THAT MAPPING IS MINE.'")
print("  AT THE BYTES THAT CONCESSION IS TOO GENEROUS.  C4 has TWO clauses and the second is a")
print("  VERDICT clause: 'A condition that returns the same verdict at the trivial connection as at")
print("  a generic one is not a formation condition.'  The edge convention returns FIRES at both.")
print("  On the second clause C4 disqualifies it directly, exactly as it disqualified readings (i)")
print("  and (ii).  THE STEELMAN'S CONCESSION HERE IS LAZY, NOT FORCED, and I record it against my")
print("  own side: its weakest self-declared seam is stronger than it says.")

print("\n== G4  COR-J, THE PREMISE MY LEG E SHOWS IS DOING THE WORK.  RESOLVES ==")
show("S3_THE_CROSSING_AUDIT_V001.md", 198, 207, "S3 audit sec2.5 -- COR-J, the body")
show("S3_THE_CROSSING_AUDIT_V001.md", 798, 798, "S3 audit sec8 correction table -- COR-J")
print("  VERDICT: EXACT.  The instruction 'Add it to the CHOICE LEDGER' was never carried out:")
print("  COR-J appears ONCE in the whole register (inside W-05's narrative, see G6) and in no")
print("  CHOICE LEDGER of any artifact.  It is the premise my leg E shows is doing the selecting.")

print("\n== G5  S2 AUDIT CHOICE LEDGER A2 -- THE CLOCK, LEDGERED WITH THE WRONG WHY ==")
show("S2_FORMATION_CONDITION_ON_K1_AUDIT_V001.md", 658, 658, "S2 audit CHOICE LEDGER A2")
print("  VERDICT: BOTH LANES ARE RIGHT ABOUT THIS AND IT IS THE HINGE.  The CHOICE is 'Time =")
print("  number of circuits'.  The ALTERNATIVES column offers only 'a real parameter t with a")
print("  Hamiltonian' and 'no time at all'.  THE WHY COLUMN JUSTIFIES IT WITH 'EDGE COUNT is")
print("  carrier-supplied combinatorics (S1 :16-22)' -- a justification for the alternative it does")
print("  not take.  The sub-choice edge -> circuit is made and never argued.")

print("\n== G6  WHAT THE REGISTER CARRIES ABOUT ALL OF THIS ==")
reg = open(os.path.join(REPO, "REGISTER_V001.md")).read()
for term in ("COR-J", "COR-F", "gauge-invariant", "edge tick", "correct trivial limit",
             "trivial-connection", "gauge covariance", "admissibility criterion"):
    print(f"  occurrences of {term!r:<26} in REGISTER_V001.md : {reg.count(term)}")
print("  -> 'correct trivial limit' appears once, as a checked PROPERTY of W-01's construction")
print("     (REGISTER:47-49).  It has never been used as a CRITERION on a convention by any row.")
print("     COR-J -- the premise that decides the convention question once leg E is run -- has")
print("     zero register rows.")
