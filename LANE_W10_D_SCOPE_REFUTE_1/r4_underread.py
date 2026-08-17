# W10-D REFUTER 1 -- LEG R4.  THREE UNDER-READS AT THE BYTES, ONE SCOPE OMISSION, AND THE ONE
# GAP LANE D SELF-FLAGGED, CLOSED.
#
# R4-A  D-12 says of B4's closed form "nobody said so".  S4 said it, verbatim, with the SAME
#       reason, on the page lane D quotes for its class multisets.
# R4-B  D-05 says "the true theorem is '|S| = 1 never fires'".  REGISTER:126 says exactly that,
#       and says it "recovers W-01's 'the root can never fire' as a special case."
# R4-C  D-24 restates S4's Theorem S4-1 ("rank 2 for every |S| >= 3") as CARRIER_INDEPENDENT /
#       HIGH with no qualification.  REGISTER:208-213 records that theorem as FAILING under
#       charge, with the counterexample.  Re-run here.
# R4-D  B4 REBUILT FROM SCRATCH and its class multiset checked against S4:578.  Lane D's own
#       self-flag: "B4's multiset is QUOTED, NOT REBUILT.  If S4's table is wrong, the four-class
#       findings are wrong with it."  It is not wrong.  Closed.
import subprocess
import numpy as np
from itertools import combinations
from fractions import Fraction as Fr

R = "/Users/bgm/MB Work/where-atoms-come-from"

print("="*100)
print("== R4-A  'B4 IS EXACT FOR THE SAME REASON AND NOBODY SAID SO' -- S4 SAID SO ==")
print("="*100)
out = subprocess.run(["sed", "-n", "586p;594,599p", f"{R}/S4_THE_MEASUREMENT_V001.md"],
                     capture_output=True, text=True).stdout
print("  S4_THE_MEASUREMENT_V001.md, lines 586 and 594-599, verbatim:")
for L in out.rstrip("\n").split("\n"):
    print("    | " + L)
print("\n  S4 states B4's closed form, its EXACTNESS, and the identical reason lane D gives")
print("  ('their squares differ by 0.2222 + 0.1111 cos y > 0' = D-12's '(8+4cos t)/36 > 0').")
print("  D-12's mathematical content is right and its NOVELTY CLAIM IS FALSE AT THE BYTES.")
print("  What IS new in D-12 is only the B0b half, which is W-03's correction, already registered")
print("  at REGISTER:236 ('B0b is log(4/9) = -0.8109302162163288 exactly, so nine of nine carrier")
print("  rates are exact, not eight').  D-12 re-derives two registered facts and files them as one")
print("  new one.  THE CORPUS'S SIGNATURE DEFECT -- UNDER-READING ITS OWN RECORD.")

print("\n"+"="*100)
print("== R4-B  D-05's 'TRUE THEOREM' IS REGISTER:126, WORD FOR WORD ==")
print("="*100)
out = subprocess.run(["sed", "-n", "124,127p", f"{R}/REGISTER_V001.md"], capture_output=True, text=True).stdout
for L in out.rstrip("\n").split("\n"):
    print("    | " + L)
print("  W-02's row states the |S| = 1 criterion AND its relation to W-01's root sentence.")
print("  D-05 rediscovers it on B0b with two byte-identical arms (leg R2-B).")

print("\n"+"="*100)
print("== R4-C  D-24 RESTATES A THEOREM THE REGISTER ALREADY RECORDS AS REFUTED UNDER CHARGE ==")
print("="*100)
out = subprocess.run(["sed", "-n", "208,213p", f"{R}/REGISTER_V001.md"], capture_output=True, text=True).stdout
for L in out.rstrip("\n").split("\n"):
    print("    | " + L)
print("\n  D-24 as filed: 'd_eff = rank of the relation lattice is 2 for EVERY occupied set of")
print("  size >= 3, three- and four-class alike.'  scope CARRIER_INDEPENDENT, confidence HIGH,")
print("  no charge qualification anywhere in the row.")
print("\n  LANE D's OWN LEG 3A, RE-RUN -- it is right, AT UNIT CHARGE, and unit charge is hard-")
print("  coded in its EXP table {'00':(0,0),'10':(1,0),'01':(0,1),'11':(1,1)}:")
EXP1 = {'00': (-1*0, 0), '10': (-1, 0), '01': (0, 1), '11': (-1, 1)}
CLS = ('00', '10', '01', '11')
print(f"    {'S':24s} {'|S|':>4s} {'rank of relation lattice':>26s}")
for r in (3, 4):
    for S in combinations(CLS, r):
        V = [np.array(EXP1[s]) for s in S]
        D = np.array([v - V[0] for v in V[1:]])
        print(f"    {str(S):24s} {len(S):4d} {int(np.linalg.matrix_rank(D)):26d}")
print("\n  THE REGISTERED COUNTEREXAMPLE, OFF UNIT CHARGE (W-03's critic, REGISTER:210-211):")
for lab, ex in [("exponents (1,0),(2,0),(3,0)", [(1, 0), (2, 0), (3, 0)]),
                ("exponents (0,0),(1,0),(2,0)", [(0, 0), (1, 0), (2, 0)]),
                ("exponents (1,1),(2,2),(3,3)", [(1, 1), (2, 2), (3, 3)])]:
    V = [np.array(e) for e in ex]
    D = np.array([v - V[0] for v in V[1:]])
    print(f"    {lab:30s} |S| = {len(ex)}   rank = {int(np.linalg.matrix_rank(D))}   "
          f"-> d_eff = 1, W-08's admissible write density is K^(-1/3), NOT K^(-1/2)")
print("  D-24 IS TRUE AT CHARGE 1 AND FALSE OFF IT, AND THE REGISTER SAYS SO EIGHT ROWS EARLIER")
print("  ('The four-class taxonomy is a charge-1 statement', REGISTER:212).  A row scoped")
print("  CARRIER_INDEPENDENT / HIGH with a live registered counterexample in another modality is")
print("  MISMARKED: the correct mark is CARRIER_INDEPENDENT **AND CHARGE-1-SCOPED**.")

print("\n"+"="*100)
print("== R4-D  B4 REBUILT.  LANE D's ONE ADMITTED GAP, CLOSED. ==")
print("="*100)
print("  S4:519 row: 'B4 spindle (two spheres glued at 2 pts)  V=6 E=8 F=4 chi=2 b0=1 b1=1 b2=2',")
print("  S4:578    : class multiset {00:1, 01:1, 10:1, 11:3}.  Lane D QUOTED both and rebuilt")
print("  neither (its self-flag).  Built here, independently, and checked.")
print("  MODEL: two spheres, each a 'pillowcase' -- a 4-cycle N-a-S-b bounding two 2-cells --")
print("  glued along the two poles N, S.  That is the only V=6 E=8 F=4 realisation of the row.")
Vn = 6                       # 0=N 1=S 2=a1 3=b1 4=a2 5=b2
edges = [(0, 2), (2, 1), (1, 3), (3, 0),     # sphere 1 equator N-a1-S-b1-N
         (0, 4), (4, 1), (1, 5), (5, 0)]     # sphere 2 equator N-a2-S-b2-N
E = len(edges)
eidx = {e: k for k, e in enumerate(edges)}
d1 = np.zeros((Vn, E))
for k, (s, t) in enumerate(edges):
    d1[t, k] += 1; d1[s, k] -= 1


def chain(seq):
    ch = np.zeros(E)
    for (s, t) in seq:
        if (s, t) in eidx:
            ch[eidx[(s, t)]] += 1
        else:
            ch[eidx[(t, s)]] -= 1
    return ch


loop1 = chain([(0, 2), (2, 1), (1, 3), (3, 0)])
loop2 = chain([(0, 4), (4, 1), (1, 5), (5, 0)])
# two 2-cells on each equator: both attach along the same 4-cycle (the two hemispheres)
d2 = np.array([loop1, loop1, loop2, loop2]).T
F = d2.shape[1]
r1, r2 = np.linalg.matrix_rank(d1), np.linalg.matrix_rank(d2)
b0, b1, b2 = Vn - r1, E - r1 - r2, F - r2
print(f"  built: V={Vn} E={E} F={F}  chi={Vn-E+F}  rank d1={r1} rank d2={r2}  b0={b0} b1={b1} b2={b2}"
      f"   max|d1.d2| = {np.abs(d1@d2).max():.1e}")
print(f"  S4:519  : V=6 E=8 F=4  chi=2                                b0=1 b1=1 b2=2")
assert (Vn, E, F, Vn-E+F, b0, b1, b2) == (6, 8, 4, 2, 1, 1, 2), "B4 reconstruction does not match S4"
inv = E - r1
print(f"  gauge invariants = E - rank d1 = {inv};  curvature = rank d2 = {r2};  flat = b1 = {b1};"
      f"  {r2} + {b1} = {r2+b1} = {inv}  -- S4:519 publishes 'gauge 5, inv 3, curv 2, flat 1'")
gF = loop1                                   # bounds (it IS a face boundary)
gC = chain([(0, 2), (2, 1), (1, 4), (4, 0)])  # N-a1-S then S-a2-N : crosses BOTH spheres
print(f"  gamma_F cycle: {np.abs(d1@gF).max():.1e}   bounds: "
      f"{np.linalg.matrix_rank(np.c_[d2, gF]) == r2}")
print(f"  gamma_C cycle: {np.abs(d1@gC).max():.1e}   bounds: "
      f"{np.linalg.matrix_rank(np.c_[d2, gC]) == r2}   (must be False -- it is the b1 generator)")
print(f"  independent in the cycle space: "
      f"{np.linalg.matrix_rank(np.c_[d2, gF, gC]) > np.linalg.matrix_rank(np.c_[d2, gF])}")
FV = {v for v in range(Vn) if any(gF[k] != 0 and v in edges[k] for k in range(E))}
CV = {v for v in range(Vn) if any(gC[k] != 0 and v in edges[k] for k in range(E))}
from collections import Counter
cnt = Counter(('1' if v in FV else '0')+('1' if v in CV else '0') for v in range(Vn))
built = {k: cnt.get(k, 0) for k in ('00', '01', '10', '11')}
print(f"  CLASS MULTISET BUILT : {built}")
print( "  S4:578 PUBLISHED     : {'00': 1, '01': 1, '10': 1, '11': 3}")
assert built == {'00': 1, '01': 1, '10': 1, '11': 3}, "B4 class multiset mismatch"
print("  MATCH.  S4's B4 row is now independently rebuilt, incidence and all.  D-31's admitted")
print("  exposure ('if S4's table is wrong, the four-class findings are wrong with it') is CLOSED")
print("  for B4 as it already was for B0b.  This STRENGTHENS lane D; it does not refute it.")

print("\n"+"="*100)
print("== R4-E  AND B4's RATE FROM THE REBUILT CARRIER, SENSE U, IN EXACT RATIONALS ==")
print("="*100)
p = [Fr(built[c], Vn) for c in ('00', '10', '01', '11')]
print(f"  SENSE-U class weights from the rebuild: {tuple(str(q) for q in p)}")
A2 = lambda ct: p[0]**2 + p[2]**2 + 2*p[0]*p[2]*ct
B2 = lambda ct: p[1]**2 + p[3]**2 + 2*p[1]*p[3]*ct
al = p[0]**2 + p[2]**2 - p[1]**2 - p[3]**2
be = 2*(p[0]*p[2] - p[1]*p[3])
print(f"  Jensen branch gap  B^2 - A^2 = {-al} + {-be} cos t,  which at cos t = -1 is "
      f"{-al + be} > 0 and at cos t = +1 is {-al - be} > 0")
print(f"  -> branch B dominates on the whole circle, m = log max(p10, p11) = log({max(p[1],p[3])})"
      f" = {float(np.log(float(max(p[1], p[3])))):.15f}")
print(f"  S4:595 publishes  -0.693147180560 EXACTLY.  log(1/2) = {np.log(0.5):.15f}")
print(f"  and the torus-zero test: cos t root = {float(-al/be) if be else 'no cos t term'}"
      f"  -> |root| > 1, NO TORUS ZERO.  D-03 and D-04 confirmed from the rebuilt incidence.")
