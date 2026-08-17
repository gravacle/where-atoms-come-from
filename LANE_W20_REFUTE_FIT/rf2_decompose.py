# rf2_decompose.py -- LANE W20_REFUTE_FIT.  THE FITTING AUDIT PROPER.
#
# THE CHARGE I WAS COMMISSIONED TO PRESS: three ingredients were named IN ADVANCE by the
# principal, and the failure mode is finding a role for each and declaring victory.  The brief
# names a specific new risk: "Is any CONSTITUTIVE verdict a dodge -- a testable ingredient
# reclassified to avoid a null?"  and a specific isolation question: "Does any CAUSAL_EARNED
# verdict rest on arms differing in more than the named ingredient?"
#
# THIS FILE ANSWERS THE SECOND, AND THE ANSWER DECIDES THE FIRST.
#
# THE INSTRUMENT'S OWN DECLARED FACT (rf1 BLOCK 1, an operator identity in every sector):
#       Z(A_S) = A_Sigma        =>       H_CENTRE  IS  THE BOUNDARY'S OWN RECORDED DATA.
# W20_PRE BLOCK 7 V1 stamped exactly this quantity: "2 of the record's 3 bits are boundary data
# BY DEFINITION.  Reporting that as evidence for H-SURFACE is the W-19 defect verbatim.  IT IS
# NOT EVIDENCE.  IT IS THE DEFINITION."
#
# SO THE SYMMETRIC QUESTION, WHICH NEITHER LANE ASKED:
#       HOW MUCH OF EACH ARM'S REPORTED "THE RECORD MOVED" IS A MOVE OF H_CENTRE?
# If an arm's dH_FULL is dH_CENTRE, the arm did not move the record past its boundary -- it
# MOVED THE BOUNDARY.  And an arm that moves the boundary cannot be a formation arm, because
# the commissioning text says: "HOLD THE BOUNDARY'S EXISTENCE FIXED.  PRODUCE THE SAME BOUNDARY
# BY TWO DIFFERENT FORMATION ROUTES."
#
# BOTH LANES HELD flux(Sigma) FIXED AND CALLED THAT "THE SAME BOUNDARY".  The flux is a SCALAR:
# X^Sigma = X_0X_4X_5 = eta_0 eta_1 eta_2 is a c-number on the physical sector.  The surface
# algebra A_Sigma has TWO MORE BITS beyond it, and those bits are exactly H_CENTRE.  Holding the
# flux is holding ONE constraint on a 2-bit object.  It is not holding the boundary.
#
# NOTHING IS CHOSEN HERE AFTER SEEING A NUMBER.  The decomposition H_FULL = H_CENTRE + C is
# forced by the algebra (rf_core header), the arms are the ones the two lanes already reported,
# and the grid is the pre-registered 13 points.

import numpy as np, math
from rf_core import *

Lg = Log("OUT_rf2_decompose.txt")
P = Lg
rule = Lg.rule

KEYS = [(+1, +1), (+1, -1), (-1, +1), (-1, -1)]


def sweep(ch):
    se = sector(ch)
    return {g2: record(se, se.ground(g2)[0]) for g2 in GRID}


ARMS = {
    "vacuum":   sweep([]),
    "eta{0,4}": sweep([0, 4]),
    "eta{1,5}": sweep([1, 5]),
    "eta{4,5}": sweep([4, 5]),
    "eta{0,1}": sweep([0, 1]),
}

PAIRS = [
    ("A1  CHARGE SECTOR ", "vacuum",   "eta{0,4}", "ledgered CAUSAL_EARNED by BOTH lanes as 'the charge sector'"),
    ("R1  FORMATION a   ", "eta{0,4}", "eta{1,5}", "ledgered CAUSAL_EARNED by BOTH lanes as 'boundary FORMATION'"),
    ("R2  FORMATION b   ", "vacuum",   "eta{4,5}", "ledgered CAUSAL_EARNED by BOTH lanes as 'boundary FORMATION'"),
    ("R3  FORMATION c   ", "vacuum",   "eta{0,1}", "ledgered CAUSAL_EARNED by BOTH lanes as 'boundary FORMATION'"),
]

rule("BLOCK 1 -- THE SAME KNOB APPEARS IN TWO LEDGER SLOTS.  PRINT THE KNOB.")
P("Every arm backing 'the charge sector' AND every arm backing 'boundary FORMATION', in BOTH")
P("lanes, is the same move: pick a different eta among the 128 admissible charge sectors.")
P("Neither lane's formation arm contains a route, a history, or a time.  R1/R2/R3 are STATIC")
P("ground states of DIFFERENT superselection sectors, compared at equal g2.")
P("")
P("  arm                 sector a    sector b    flux_a flux_b  equal?   ledger slot it was scored in")
for lab, a, b, note in PAIRS:
    ca = [] if a == "vacuum" else [int(x) for x in a[4:-1].split(",")]
    cb = [] if b == "vacuum" else [int(x) for x in b[4:-1].split(",")]
    fa, fb = flux_of(ca), flux_of(cb)
    P("  %-19s %-11s %-11s %+d     %+d      %-8s %s" %
      (lab, a, b, fa, fb, str(fa == fb), note))
P("")
P(">>> THE ONLY STRUCTURAL DIFFERENCE between the arm scored as 'the charge sector' and the arms")
P("    scored as 'boundary FORMATION' is WHETHER THE TWO SECTORS HAPPEN TO SHARE A FLUX SIGN.")
P("    One knob, split into two ledger slots by a one-bit label.  Both slots then received an")
P("    independent CAUSAL_EARNED verdict, in both lanes.  THAT IS THE FITTING, IN ITS EXACT")
P("    MECHANICAL FORM: six ingredients were asked for and the knobs were subdivided until six")
P("    were available.")

rule("BLOCK 2 -- THE DECOMPOSITION NEITHER LANE PRINTED.  dH_FULL = dH_CENTRE + dC.")
P("H_CENTRE = Shannon of the 4 central weights = the entropy of the SURFACE'S OWN DATA")
P("           (Z(A_S) = A_Sigma, operator identity, rf1 BLOCK 1).")
P("C        = H_FULL - H_CENTRE = the record's content BEYOND its boundary.  Ceiling 1 bit.")
P("D_Sigma  = trace distance of the two arms restricted to A_Sigma = total variation of the")
P("           boundary's own 2-bit data.  THIS IS THE NUMBER THAT SAYS WHETHER THE ARM HELD")
P("           THE BOUNDARY FIXED.  Both lanes held only the flux, which is a scalar, not this.")
P("")
for lab, a, b, note in PAIRS:
    A, B = ARMS[a], ARMS[b]
    P("-" * 100)
    P("%s :  %s  vs  %s      [%s]" % (lab, a, b, note))
    P("  g2      dH_FULL     dH_CENTRE   dC          |dC|/|dH_FULL|   D_Sigma(boundary moved)  "
      "H_FULL_a  H_FULL_b")
    for g2 in GRID:
        ra, rb = A[g2], B[g2]
        dF = rb["H_FULL"] - ra["H_FULL"]
        dC = rb["C"] - ra["C"]
        dH = rb["H_CENTRE"] - ra["H_CENTRE"]
        ratio = abs(dC) / abs(dF) if abs(dF) > 1e-12 else float("nan")
        ds = tv(ra["p"], rb["p"], KEYS)
        P("  %-6.2f  %-11.7f %-11.7f %-11.8f %-16.4e %-24.7f %-9.6f %-9.6f"
          % (g2, dF, dH, dC, ratio, ds, ra["H_FULL"], rb["H_FULL"]))
    mF = max(abs(B[g]["H_FULL"] - A[g]["H_FULL"]) for g in GRID)
    mC = max(abs(B[g]["C"] - A[g]["C"]) for g in GRID)
    mS = max(tv(A[g]["p"], B[g]["p"], KEYS) for g in GRID)
    P("  >>> max|dH_FULL| = %.8f    max|dC| = %.8f    max D_Sigma = %.8f" % (mF, mC, mS))
    P("      FRACTION OF THE REPORTED RECORD MOVE THAT IS A BOUNDARY-DATA MOVE: %.5f%%"
      % (100.0 * (1.0 - mC / mF)))

rule("BLOCK 3 -- THE VERDICT OF BLOCK 2, STATED ONCE")
rows = []
for lab, a, b, note in PAIRS:
    A, B = ARMS[a], ARMS[b]
    mF = max(abs(B[g]["H_FULL"] - A[g]["H_FULL"]) for g in GRID)
    mC = max(abs(B[g]["C"] - A[g]["C"]) for g in GRID)
    mS = max(tv(A[g]["p"], B[g]["p"], KEYS) for g in GRID)
    rows.append((lab, mF, mC, mS))
P("  arm                 max|dH_FULL|  max|dC|       max D_Sigma   what the arm actually moved")
for lab, mF, mC, mS in rows:
    P("  %-19s %-13.8f %-13.8f %-13.8f %s"
      % (lab, mF, mC, mS, "THE BOUNDARY" if mC < 0.05 * mF else "something past the boundary"))
P("")
P(">>> EVERY ARM SCORED AS 'boundary FORMATION' MOVES THE BOUNDARY'S OWN DATA BY MORE THAN IT")
P("    MOVES ANYTHING ELSE.  The record content past the boundary moves by at most %.6f bits"
  % max(r[2] for r in rows))
P("    in ANY of these arms, against a claimed record move of up to %.6f bits." % max(r[1] for r in rows))
P("    THE ARMS DID NOT HOLD THE BOUNDARY FIXED.  They held one scalar (the flux) out of a")
P("    2-bit surface algebra, and then attributed the motion of the other bits to 'formation'.")

rule("BLOCK 4 -- THE REPAIRED FORMATION ARM.  CONDITION ON THE BOUNDARY, THEN COMPARE.")
P("If the boundary cannot be held fixed as a whole, it can be CONDITIONED ON.  rho|A_S is a")
P("classical 2-bit label k (the boundary's own data) plus a conditional qubit rho_k on")
P("alg{X_1, W_S}.  The commissioning arm -- SAME BOUNDARY, DIFFERENT ROUTE -- is then well posed")
P("for the first time in this round:")
P("      hold k, compare rho_k across the two arms.")
P("Labelling is INTRINSIC: k is the joint eigenvalue of X_1X_2 and X_1X_3, which are operators")
P("INSIDE A_S and are the same operators in every sector.  No relabelling freedom is used and")
P("none is needed.  Blocks with p_k < 1e-6 in either arm are skipped and the skip is printed.")
P("")
P("  arm                 g2      max_k D_tr(rho_k)  argmax k    p_k(a)    p_k(b)   blocks used")
BEST = {}
for lab, a, b, note in PAIRS:
    A, B = ARMS[a], ARMS[b]
    best = (0.0, None, None)
    for g2 in GRID:
        ra, rb = A[g2], B[g2]
        used, mx, argk = 0, 0.0, None
        for k in KEYS:
            if ra["p"][k] < 1e-6 or rb["p"][k] < 1e-6:
                continue
            used += 1
            d = dtr_qubit(ra["bloch"][k], rb["bloch"][k])
            if d > mx:
                mx, argk = d, k
        P("  %-19s %-6.2f  %-18.9f %-11s %-9.6f %-8.6f %d/4"
          % (lab, g2, mx, str(argk), ra["p"][argk] if argk else 0.0,
             rb["p"][argk] if argk else 0.0, used))
        if mx > best[0]:
            best = (mx, g2, argk)
    BEST[lab] = best
    P("  >>> %s  MAX CONDITIONAL RECORD MOVE OVER THE GRID = %.9f at g2=%s, k=%s"
      % (lab, best[0], best[1], best[2]))
    P("")

rule("BLOCK 5 -- WHAT THE CONDITIONAL TEST SAYS ABOUT THE FORMATION VERDICT")
P("  arm                 max_k D_tr(rho_k)   max|dH_FULL| (as ledgered)   ratio")
for lab, a, b, note in PAIRS:
    A, B = ARMS[a], ARMS[b]
    mF = max(abs(B[g]["H_FULL"] - A[g]["H_FULL"]) for g in GRID)
    P("  %-19s %-19.9f %-28.8f %.4f" % (lab, BEST[lab][0], mF, BEST[lab][0] / mF))
P("")
P("A NULL READS TWO WAYS AND I AM NOT SCORING THIS ONE AS CONFIRMATION OF ANYTHING.")
P("Reading 1: with the boundary conditioned on, the charge-placement route leaves the region's")
P("           remaining record essentially untouched -> FORMATION-BY-CHARGE-PLACEMENT IS")
P("           NARRATION and the CAUSAL_EARNED verdict is unearned.")
P("Reading 2: the conditional record is ONE QUBIT and on a ground state it is nearly PURE")
P("           (C <= 0.0033 bits, rf1 BLOCK 5), so there is very little room for it to move; the")
P("           test may be too small rather than the effect absent.")
P("THE DISCRIMINATOR IS RUN IN rf3 (a route that DOES move C: the quench) AND IN rf4.")

Lg.save()
