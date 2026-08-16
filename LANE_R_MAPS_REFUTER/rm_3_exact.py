"""
rm_3_exact.py -- LANE R (MAPS REFUTER).
Independent re-derivation of every number the attack rests on, plus the
existence criterion for class-compatible maps and the rank-G non-invariance.
Nothing here reuses the closed forms of rm_2; the check route is DIRECT
SCHEDULE-B SIMULATION, which uses no Mahler identity at all.
"""
import numpy as np, math, sys, itertools
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from rm_lib import *

BAR = "=" * 78
print(BAR); print("R1  INDEPENDENT RE-DERIVATION: direct schedule-B simulation"); print(BAR)
f0, c0 = 1.0, math.sqrt(2.0)
rk, _ = relation_lattice_rank(f0, c0)
print(f"   (f,c) = (1.0, sqrt(2)), relation lattice rank = {rk} -> generic, lambda_B = m(P)")
print(f"   N = 4e6.  Convergence of a schedule-B average is O(1/N) at best; deviations")
print(f"   of order 1e-6 are the method's own error, not disagreement.\n")
CASES = [
    ("B1  SENSE U",            (0.0, 0.4, 0.4, 0.2)),
    ("B1s SENSE U",            (0.0, 5/11, 5/11, 1/11)),
    ("B1p SENSE U",            (0.0, 0.5, 0.5, 0.0)),
    ("B1q SENSE U",            (1/7, 3/7, 3/7, 0.0)),
    ("B0a SENSE U",            (2/9, 4/9, 3/9, 0.0)),
    ("B0b SENSE U",            (4/9, 2/9, 1/9, 2/9)),
    ("A2 collapse (3,4,4)/11", (0.0, 4/11, 4/11, 3/11)),
    ("A2 collapse (3,3,5)/11", (0.0, 5/11, 3/11, 3/11)),
    ("A2 collapse (2,4,5)/11", (0.0, 5/11, 4/11, 2/11)),
    ("A3 transported (0,4/7,3/7,0)", (0.0, 4/7, 3/7, 0.0)),
    ("A5 K1[2,0]",             (0.0, 4/7, 2/7, 1/7)),
    ("A5 K1[10,0]",            (0.0, 12/15, 2/15, 1/15)),
    ("A5 K1[100,0]",           (0.0, 102/105, 2/105, 1/105)),
]
print(f"   {'case':32s} {'closed form':>18s} {'direct N=4e6':>18s} {'dev':>10s}")
worst = 0.0
for nm, pi in CASES:
    lc = lambda_B_closed(pi)
    ld = lambda_B_direct(pi, f0, c0, N=4000000)
    worst = max(worst, abs(lc - ld))
    print(f"   {nm:32s} {lc:18.12f} {ld:18.12f} {abs(lc-ld):10.2e}")
print(f"   worst |closed - direct| = {worst:.2e}   (S4's own worst on this check: 3.0e-06)")

print("\n" + BAR)
print("R2  TWO INDEPENDENT CONSTRUCTIONS OF B1s AGREE")
print(BAR)
a = ALL["B1s"]().pi_uniform()
b = K1_partial_subdiv(3, 3).pi_uniform()
print(f"   explicit K1_subdivided()      pi = {tuple(round(x,12) for x in a)}")
print(f"   parametric K1_partial_subdiv(3,3) pi = {tuple(round(x,12) for x in b)}")
print(f"   max component difference = {max(abs(x-y) for x,y in zip(a,b)):.2e}")

print("\n" + BAR)
print("R3  WHEN DOES A CLASS-COMPATIBLE MAP EXIST AT ALL?")
print(BAR)
print("""
PROPOSITION (rm-2).  A class-compatible vertex map c : K -> K' exists iff
    supp(classes(K))  SUBSET OF  supp(classes(K')) ,
i.e. iff every vertex class OCCUPIED on K is also occupied on K'.
PROOF.  (=>) c(w) must have class(w), so that class is occupied on K'.
(<=) send each w to any chosen vertex of K' of class(w).  QED

COROLLARY.  The condition the claim calls 'a map is present' is really
'the two carriers already agree on which classes are occupied' -- which is the
first half of the very datum lambda is a function of.  It cannot be weakened to
the existence of a map, a homeomorphism, a homotopy equivalence, or a collapse.
""")
names = ["B1", "B1s", "B1p", "B1q", "B0a", "B0b"]
Ks = {n: ALL[n]() for n in names}
sup = {n: set(Ks[n].class_counts()) for n in names}
print(f"   {'':6s}" + "".join(f"{n:>7s}" for n in names) + "   <- target K'")
for s in names:
    print(f"   {s:6s}" + "".join(
        f"{('yes' if sup[s] <= sup[t] else '.'):>7s}" for t in names))
print("   (row = source K, column = target K'; 'yes' = a class-compatible map exists)")
print(f"\n   occupied classes: " + ", ".join(
    f"{n}={sorted(''.join(map(str,c)) for c in sup[n])}" for n in names))

print("\n" + BAR)
print("R4  EDGES ON NEITHER LOOP -- where subdivision breaks class-compatibility")
print(BAR)
print("""
Subdividing an edge on NEITHER designated loop creates a class-(0,0) vertex.  If the
carrier had none, no class-compatible map to the un-subdivided carrier can exist
(R3).  K1 is the special case where EVERY edge lies on a loop -- which is the only
reason the claim's exhibit worked.
""")
print(f"   {'carrier':8s} {'E':>4s} {'on gF or gC':>12s} {'on NEITHER':>11s} "
      f"{'has (0,0)?':>11s} {'subdiv safe?':>13s}")
for n in names:
    K = Ks[n]
    onl = set(K.gF) | set(K.gC)
    neither = K.nE - len(onl)
    has00 = (0, 0) in K.class_counts()
    safe = (neither == 0) or has00
    print(f"   {n:8s} {K.nE:4d} {len(onl):12d} {neither:11d} {str(has00):>11s} {str(safe):>13s}")
print("""
   B1p has 1 edge (the bridge) on neither loop and no (0,0) vertex.  Subdividing it --
   S4's own Control 4 -- destroys class-compatibility at the first subdivision.
   The claim generalised from the one carrier in the corpus where it cannot fail.""")

print("\n" + BAR)
print("R5  rank G IS NOT A HOMEOMORPHISM INVARIANT EITHER  (this WIDENS the box)")
print(BAR)
def rankG(pi):
    corners = [(0,0),(1,0),(0,1),(1,1)]
    S = [corners[i] for i, w in enumerate(pi) if w > 1e-14]
    if len(S) <= 1:
        return 0
    diffs = [[S[i][0]-S[0][0], S[i][1]-S[0][1]] for i in range(1, len(S))]
    return exact_rank(diffs)
print("""
S4 Control 3 asserted only that the VALUE of lambda is not a homeomorphism invariant.
The B1p/B1q pair shows something strictly stronger, which S4 filed under a different
control (Control 4, 'the spectator') and therefore never stated:
""")
for n in ["B1p", "B1q", "B1", "B1s"]:
    pi = Ks[n].pi_uniform()
    print(f"   {n:5s} pi = ({pi[0]:.4f},{pi[1]:.4f},{pi[2]:.4f},{pi[3]:.4f})  "
          f"rank G = {rankG(pi)}   lambda = {lambda_B_closed(pi):.12f}")
print("""
   B1p and B1q are the SAME TOPOLOGICAL SPACE (B1q = B1p with the bridge subdivided).
   rank G moves 1 -> 2 under a subdivision.  So the FORMATION CRITERION of W-02 /
   Theorem S4-1 -- not merely the rate -- fails to be a homeomorphism invariant.
   Under a homeomorphism, formation goes from 'sees the product W_F.W_C only' to
   'separates curvature from flat holonomy'.  That is a qualitative change and it is
   NOT in the register.""")

print("\n" + BAR)
print("R6  A CLOSING SANITY CHECK ON THE CLAIM'S OWN ARITHMETIC")
print(BAR)
print("   The claim's stated transported weights (0,5/11,5/11,1/11):")
print(f"     sum = {0+5/11+5/11+1/11:.15f}  (must be 1)")
print(f"     lambda = m(5/11 + 5/11 x + 1/11 y) = {lambda_B_closed((0.0,5/11,5/11,1/11)):.12f}")
print(f"     = S4's OWN B1s SENSE U row, -0.724759919461, to {abs(lambda_B_closed((0.0,5/11,5/11,1/11))+0.724759919461):.1e}")
print("""     -> the claim's 'both sides' number IS S4's B1s number.  The claim did not
        compute a new quantity; it evaluated S4's B1s row on B1 as well.  The B1 side
        of S4's Control 3, -0.756573586, was simply discarded and never defended.""")
print("\n" + BAR)
