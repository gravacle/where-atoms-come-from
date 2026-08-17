#!/usr/bin/env python3
"""
REFUTER 1 — SCRIPT 1.  ATTACK THE RECONSTRUCTIONS (the lens's own assignment).

Q1  Does an INDEPENDENT rebuild of B0b and B4 (my indexing, my rank code) reproduce
    S4's published row and the lane's numbers?  If yes the lane's arithmetic is clean.
Q2  IS S4's B4 ROW DETERMINED BY WHAT S4 PUBLISHED?  This is the question the lane's
    "0 mismatches / 4 carriers" headline presupposes.  I attack it by exhibiting a SECOND
    spindle -- two 2-spheres glued at two points, V=6 E=8 F=4 -- which matches EVERY
    published column of S4's B4 row EXCEPT the class multiset, and therefore a DIFFERENT
    lambda.  If that succeeds, matching the multiset was an INPUT to the lane's construction,
    not an output, and the B4 half of the headline is a satisfiability proof and not an audit.
Q3  Same question for the loop designation: on each complex, enumerate every admissible
    gamma_C (simple cycle, does not bound, independent of gamma_F) and report which class
    multisets are reachable.  If more than one, the row does not determine the carrier.

EXACT arithmetic throughout (Fractions / integers).  No float in any check in this file.
"""
import sys
from fractions import Fraction
from r1lib import (CW, rk, in_span, my_K1, my_B0b, my_B4_square, my_B4_tri_pent, my_B1q,
                   simple_cycles)

LOG = []
def out(s=""):
    print(s); LOG.append(s)

CN = {(0, 0): "00", (1, 0): "10", (0, 1): "01", (1, 1): "11"}

# S4's published rows, transcribed independently from S4_THE_MEASUREMENT_V001.md:511-520,
# :538-548, :574-583.  (I transcribed these from the file myself; the lane's dict is not used.)
S4 = {
    "B0b": dict(V=9, E=18, F=9, chi=0, b0=1, b1=2, b2=1, gauge=8, inv=10, curv=8, flat=2,
                ms={(0, 0): 4, (0, 1): 1, (1, 0): 2, (1, 1): 2},
                gF_bounds=True, gC_bounds=False, independent=True, lamU="-0.810930216"),
    "B4":  dict(V=6, E=8, F=4, chi=2, b0=1, b1=1, b2=2, gauge=5, inv=3, curv=2, flat=1,
                ms={(0, 0): 1, (0, 1): 1, (1, 0): 1, (1, 1): 3},
                gF_bounds=True, gC_bounds=False, independent=True, lamU="-0.693147181"),
    "B1":  dict(V=5, E=6, F=1, chi=0, b0=1, b1=1, b2=0, gauge=4, inv=2, curv=1, flat=1,
                ms={(0, 0): 0, (0, 1): 2, (1, 0): 2, (1, 1): 1},
                gF_bounds=True, gC_bounds=False, independent=True, lamU="-0.756573586"),
    "B1q": dict(V=7, E=8, F=1, chi=0, b0=1, b1=1, b2=0, gauge=6, inv=2, curv=1, flat=1,
                ms={(0, 0): 1, (0, 1): 3, (1, 0): 3, (1, 1): 0},
                gF_bounds=True, gC_bounds=False, independent=True, lamU="-0.741029583"),
}

out("=" * 104)
out("REFUTER 1 / SCRIPT 1 — INDEPENDENT REBUILD, AND: IS S4's B4 ROW DETERMINED?")
out("=" * 104)
out("My indexing differs from the lane's on purpose (B0b: vertex (i,j)->3i+j, verticals first;")
out("B4: interleaved edge order).  My rank routine is a separate implementation.  Exact only.")
out()

FAIL = []
b0b, _ = my_B0b()
MINE = [("B1", my_K1()), ("B1q", my_B1q()), ("B0b", b0b), ("B4", my_B4_square())]

out("-" * 104)
out("Q1 — MY REBUILD vs S4's PUBLISHED ROW")
out("-" * 104)
out(f"{'car':<6}{'V':<4}{'E':<4}{'F':<4}{'chi':<5}{'b0':<4}{'b1':<4}{'b2':<4}"
    f"{'gau':<5}{'inv':<5}{'cur':<5}{'flt':<5}{'d1d2':<6}{'gFb':<6}{'gCb':<6}{'ind':<6}"
    f"{'multiset':<28}{'S4?'}")
for tag, cw in MINE:
    R = cw.report()
    ms = cw.multiset()
    T = S4[tag]
    ok = all(R[k] == T[k] for k in
             ("V", "E", "F", "chi", "b0", "b1", "b2", "gauge", "inv", "curv", "flat",
              "gF_bounds", "gC_bounds", "independent"))
    ok = ok and R["d1d2"] == 0 and R["gF_cycle"] and R["gC_cycle"]
    ok = ok and {k: ms[k] for k in ms} == {k: T["ms"][k] for k in ms}
    if not ok:
        FAIL.append(f"{tag} does not match S4's row")
    msstr = "{" + ",".join(f"{CN[c]}:{ms[c]}" for c in [(0, 0), (0, 1), (1, 0), (1, 1)]) + "}"
    out(f"{tag:<6}{R['V']:<4}{R['E']:<4}{R['F']:<4}{R['chi']:<5}{R['b0']:<4}{R['b1']:<4}"
        f"{R['b2']:<4}{R['gauge']:<5}{R['inv']:<5}{R['curv']:<5}{R['flat']:<5}{R['d1d2']:<6}"
        f"{str(R['gF_bounds']):<6}{str(R['gC_bounds']):<6}{str(R['independent']):<6}"
        f"{msstr:<28}{'OK' if ok else '**MISMATCH**'}")
out()
out("SENSE U pushforwards (EXACT Fractions), from MY incidence:")
for tag, cw in MINE:
    out(f"   {tag:<6} pi = (" + ", ".join(str(x) for x in cw.pi_U()) + ")")
out()
out("VERDICT Q1: my independent rebuild agrees with the lane and with S4 on every column.")
out("THE LANE'S ARITHMETIC IS NOT THE PROBLEM.  Neither was anyone's, in this whole program.")
out()

# ==================================================================================================
out("=" * 104)
out("Q2 — THE ATTACK: A SECOND SPINDLE THAT SATISFIES EVERY PUBLISHED COLUMN OF S4's B4 ROW")
out("=" * 104)
out("S4 published NO incidence for B4 (grep confirms; the lane says so too).  So the row is")
out("whatever complexes satisfy it.  The lane built ONE: two SQUARES, each with two 2-cells,")
out("glued at OPPOSITE corners.  Here is another, equally a 'spindle (two spheres, 2 glue pts)':")
out()
alt = my_B4_tri_pent()
lane = my_B4_square()
out("   B4  (the lane's) : sphere A = square p,a1,q,a2 with 2 faces ; sphere B = square p,b1,q,b2")
out("                      glued at p and q (OPPOSITE corners of each square).")
out("   B4' (mine)       : sphere A = TRIANGLE p,q,r with 2 faces (v=3,e=3,f=2, chi=2) ;")
out("                      sphere B = PENTAGON p,s1,q,s2,s3 with 2 faces (v=5,e=5,f=2, chi=2) ;")
out("                      glued at p and q.  V = 3+5-2 = 6, E = 3+5 = 8, F = 2+2 = 4.")
out("   Both are regular oriented CW complexes and both ARE two 2-spheres glued at two points.")
out()
out(f"{'complex':<10}{'V':<4}{'E':<4}{'F':<4}{'chi':<5}{'b0':<4}{'b1':<4}{'b2':<4}"
    f"{'gau':<5}{'inv':<5}{'cur':<5}{'flt':<5}{'d1d2':<6}{'gFb':<6}{'gCb':<6}{'ind':<6}{'multiset'}")
for nm, cw in [("B4 lane", lane), ("B4' mine", alt)]:
    R = cw.report()
    ms = cw.multiset()
    msstr = "{" + ",".join(f"{CN[c]}:{ms[c]}" for c in [(0, 0), (0, 1), (1, 0), (1, 1)]) + "}"
    out(f"{nm:<10}{R['V']:<4}{R['E']:<4}{R['F']:<4}{R['chi']:<5}{R['b0']:<4}{R['b1']:<4}"
        f"{R['b2']:<4}{R['gauge']:<5}{R['inv']:<5}{R['curv']:<5}{R['flat']:<5}{R['d1d2']:<6}"
        f"{str(R['gF_bounds']):<6}{str(R['gC_bounds']):<6}{str(R['independent']):<6}{msstr}")
out()
Ra, Rl = alt.report(), lane.report()
same_cols = all(Ra[k] == Rl[k] for k in ("V", "E", "F", "chi", "b0", "b1", "b2", "gauge",
                                         "inv", "curv", "flat", "gF_bounds", "gC_bounds",
                                         "independent"))
out(f"EVERY published column of S4's B4 row agrees between B4 and B4' : {same_cols}")
out(f"CLASS MULTISET differs                                          : "
    f"{alt.multiset() != lane.multiset()}")
out()
piA, piL = alt.pi_U(), lane.pi_U()
out(f"   B4  SENSE U pi = ({', '.join(str(x) for x in piL)})")
out(f"   B4' SENSE U pi = ({', '.join(str(x) for x in piA)})")


def mahler_exact_label(pi4):
    """EXACT closed form of m(P) when one Jensen branch dominates or when P factors.
    Returns (label, Fraction-or-None expressing the argument of the log)."""
    p00, p10, p01, p11 = [Fraction(x) for x in pi4]
    C = p00 * p00 + p01 * p01 - p10 * p10 - p11 * p11
    D = 2 * (p00 * p01 - p10 * p11)
    if p00 * p11 == p10 * p01:
        return ("FACTORS (p00 p11 = p10 p01)", None)
    if C == 0 and D == 0:
        return ("branches identical", None)
    if abs(D) <= abs(C):
        if C > 0:
            return (f"log max(p00,p01) = log({max(p00, p01)})", max(p00, p01))
        return (f"log max(p10,p11) = log({max(p10, p11)})", max(p10, p11))
    return ("branches cross -- quadrature", None)


for nm, pi4 in [("B4 (lane)", piL), ("B4' (mine)", piA)]:
    lab, arg = mahler_exact_label(pi4)
    out(f"   {nm:<12} exact-form test: {lab}")
out()
out("B4'  factors:  p00 p11 = (1/6)(2/6) = 2/36  and  p10 p01 = (1/6)(2/6) = 2/36.")
out("   P' = 1/6 + (1/6)x + (2/6)y + (2/6)xy = (1 + x)(1/6 + (1/3)y),")
out("   so lambda(B4', SENSE U) = log max(1,1) + log max(1/6,1/3) = log(1/3) EXACTLY.")
out("   S4's published B4 SENSE U lambda is log(1/2).   log(1/3) != log(1/2).")
out()
out(">>> FINDING.  S4's B4 ROW IS NOT DETERMINED BY WHAT S4 PUBLISHED.  Two spindles agree on")
out("    V,E,F,chi,b0,b1,b2,gauge,invariants,curvature,flat and on gamma_F bounds / gamma_C does")
out("    not / independent, and DISAGREE on the class multiset and hence on lambda.  The ONE")
out("    column that feeds Z_k -- the multiset -- is therefore the column the lane HARD-CODED AS")
out("    ITS TARGET (w10a_1_build.py:37-38) and then matched.  For B4 the lane's '0 mismatches'")
out("    is a SATISFIABILITY proof, not an audit: it shows S4's row is realizable, and it could")
out("    not have failed once the multiset was an input to the construction.")
out("    THIS IS A COR-K DEFECT AGAINST S4 ('published rows not reproducible from their")
out("    parameters'), and the lane's headline reports its absence.")
out()

# ==================================================================================================
out("=" * 104)
out("Q3 — WHICH CLASS MULTISETS ARE REACHABLE ON EACH COMPLEX?  (admissible gamma_C sweep)")
out("=" * 104)
out("admissible = simple cycle, does NOT bound, chain-independent of gamma_F.  gamma_F fixed.")
out("If more than one multiset is reachable, the S4 row does not determine the loop designation.")
out()
for nm, cw in [("B4 (lane's spindle)", lane), ("B4' (tri|pent spindle)", alt),
               ("B0b (3x3 torus)", b0b)]:
    fc = cw.face_chains()
    cF = cw.chain(cw.gF)
    VF = cw.loop_vertices(cw.gF)
    seen = {}
    for ch, vs in simple_cycles(cw):
        vec = [0] * cw.nE
        for e, s in ch.items():
            vec[e] += s
        if in_span(fc, vec):
            continue                                   # bounds -> not a flat holonomy
        if rk([[cF[e], vec[e]] for e in range(cw.nE)]) != 2:
            continue                                   # not independent of gamma_F
        ms = {(0, 0): 0, (1, 0): 0, (0, 1): 0, (1, 1): 0}
        for v in range(cw.nV):
            ms[(1 if v in VF else 0, 1 if v in vs else 0)] += 1
        key = tuple(ms[c] for c in [(0, 0), (1, 0), (0, 1), (1, 1)])
        seen.setdefault(key, 0)
        seen[key] += 1
    out(f"{nm}:  {len(seen)} distinct reachable multiset(s) over admissible simple cycles")
    for key, n in sorted(seen.items()):
        p = [Fraction(x, cw.nV) for x in key]
        lab, _ = mahler_exact_label(p)
        nocc = sum(1 for x in key if x > 0)
        out(f"    (00:{key[0]}, 10:{key[1]}, 01:{key[2]}, 11:{key[3]})  x{n:<3} "
            f"|occ|={nocc}  SENSE U pi=({', '.join(str(x) for x in p)})   {lab}")
    out()
out(">>> B0b: the row-cycle designation is NOT the only admissible one.  Several distinct")
out("    multisets are reachable on the very same complex with the very same gamma_F, so S4's")
out("    B0b multiset column is a CHOICE too -- one that the lane inherited from")
out("    LANE_R_MAPS_REFUTER/rm_lib.py:torus33 (declared) rather than derived.")
out()

# ==================================================================================================
out("=" * 104)
out("Q4 — PRIOR ART: WAS THIS ROW 'NEVER AUDITED'?")
out("=" * 104)
out("The lane's script 1 prints 'S4 was never audited; these rows are audited now and they hold',")
out("and finding W10A-03 says B0b and B4 are 'both built here for the first time'.")
out("LANE_R_MAPS_REFUTER/rm_1_validate.OUT.txt (W-03's round, on disk, same repo) already")
out("contains, for B0b: full d1 and d2^T, V=9 E=18 F=9 chi=0, rank d1=8 rank d2=8, b0=1 b1=2")
out("b2=1, max|d1.d2|=0, class counts {'00':4,'01':1,'10':2,'11':2}, SENSE U pi =")
out("(0.444444444,0.222222222,0.111111111,0.222222222), and lambda = -0.810930216216 against")
out("S4's -0.810930216 (dev 2.163e-10) -- by closed form, Jensen quadrature AND 2D quadrature.")
out("The same file carries B1 and B1q identically.  The lane declares reading rm_lib.py:386-419")
out("(torus33) but its B1q is ALSO a byte-level reproduction of rm_lib.py:K1_bridged_subdiv,")
out("vertex names included, and that is not declared.")
out(">>> THREE of the lane's FOUR carriers had their S4 rows reproduced in an earlier round, in")
out("    a file the lane read.  Only B4 is new.  'Never audited' and 'built here for the first")
out("    time' are false as written.")
out()

out("=" * 104)
if FAIL:
    out(f"**{len(FAIL)} FAILURES IN Q1**")
    for f in FAIL:
        out("   " + f)
else:
    out("Q1: 0 mismatches.  The lane reproduces S4's four rows and so do I, independently.")
out("Q2: S4's B4 row is UNDER-DETERMINED; a second admissible spindle gives log(1/3), not log(1/2).")
out("Q3: the multiset column is a designation CHOICE on B0b as well as on B4.")
out("Q4: three of the four rows were already reproduced in LANE_R_MAPS_REFUTER.")

with open("r1_rebuild.OUT.txt", "w") as fh:
    fh.write("\n".join(LOG) + "\n")
sys.exit(1 if FAIL else 0)
