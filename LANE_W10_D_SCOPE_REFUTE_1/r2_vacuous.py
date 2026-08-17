# W10-D REFUTER 1 -- LEG R2.  THE TWO CONTROLS LANE D CALLS NEW ARE VACUOUS.
#
# The brief's rule: "'COULD NOT HAVE FAILED' VOIDS A CONTROL, NEVER A THEOREM", and W-08's
# isolation audit: "the commonest FATAL defect is not 'two variables moved' -- it is ZERO
# variables moved... Diff your arms."
#
# R2-A.  Lane D's finding D-02 rests on two measurements it labels NEW:
#          "c -> -c flips the SAME 99785 on three classes"   and   "the EXCHANGE f <-> c flips
#          0 of 200000 on THREE classes as well as on four".
#        Both are IDENTITIES of complex conjugation.  The character set is closed under
#        conjugation composed with each of those maps, and 0 in conv(S) iff 0 in conv(conj S).
#        The arms are not byte-identical; they are IDENTICAL AFTER A REFLECTION THE CRITERION
#        CANNOT SEE.  Neither could have come out any other way on ANY carrier, ANY class set.
#
# R2-B.  Lane D's finding D-05 exhibits "a state on any ONE class-00 vertex" and "a state SPREAD
#        OVER ALL FOUR class-00 vertices" as two rows of a table.  The sealed code hands the
#        SAME numpy array np.array([1.,0,0,0]) to both.  BYTE-IDENTICAL ARMS, printed as two
#        confirmations, in the lane whose PUBLISHED_CONVENTIONS promises "Every script PRINTS A
#        DIFF OF ITS ARMS before reporting, because the corpus's commonest fatal defect is a
#        control whose two arms are byte-identical."  Leg 3C prints no arm diff.
import numpy as np
from itertools import combinations

rng = np.random.default_rng(20260816)
N = 200000
f = rng.uniform(-np.pi, np.pi, N)
c = rng.uniform(-np.pi, np.pi, N)

CHAR = {'00': lambda f, c: np.ones_like(f, dtype=complex),
        '10': lambda f, c: np.exp(-1j*f),
        '01': lambda f, c: np.exp(1j*c),
        '11': lambda f, c: np.exp(1j*(c-f))}


def hull(pts):
    """sorted-angle max-gap; identical to lane D's H1."""
    A = np.sort(np.angle(np.stack(pts, axis=0)), axis=0)
    g = np.diff(np.concatenate([A, A[:1] + 2*np.pi], axis=0), axis=0)
    return g.max(axis=0) <= np.pi + 1e-12


def charset(occ, ff, cc):
    return [CHAR[o](ff, cc) for o in occ]


ROWS = [("B1  K1 as handed        ", ('10', '01', '11')),
        ("B1q K1 + spectator      ", ('00', '10', '01')),
        ("B1p K1-bridged          ", ('10', '01')),
        ("B0b/B4 four classes     ", ('00', '10', '01', '11'))]

print("="*100)
print("== R2-A  THE f<->c AND c->-c ARMS ARE THE COMPLEX CONJUGATES OF ARMS ALREADY RUN ==")
print("="*100)
print("  u = e^{-if}, v = e^{ic}, uv = e^{i(c-f)}, chi_00 = 1.  Then, as SETS:")
print("     S(c, f)  = conj( S(f, c) )        [the f<->c exchange]")
print("     S(f, -c) = conj( S(-f, c) )       [c->-c versus f->-f]")
print("  and 0 in conv(S) <=> 0 in conv(conj S), because conjugation is a linear isometry of R^2")
print("  fixing 0.  So the exchange arm CANNOT differ from the base arm and the c->-c arm CANNOT")
print("  differ from the f->-f arm -- on ANY occupied set, at ANY draw, for ANY carrier.")
print("\n  EXHIBITED AT THE ELEMENTS, not just at the counts (max over 200000 draws):")
print(f"  {'occupied set':26s} {'max|S(c,f) - conj S(f,c)|':>27s} {'max|S(f,-c) - conj S(-f,c)|':>29s}")
for lab, occ in ROWS:
    a = charset(occ, c, f)
    b = [z.conjugate() for z in charset(occ, f, c)]
    d1 = max(np.abs(np.sort_complex(np.stack(a, 0).T[i]) - np.sort_complex(np.stack(b, 0).T[i])).max()
             for i in range(0, N, 20000))
    a2 = charset(occ, f, -c)
    b2 = [z.conjugate() for z in charset(occ, -f, c)]
    d2 = max(np.abs(np.sort_complex(np.stack(a2, 0).T[i]) - np.sort_complex(np.stack(b2, 0).T[i])).max()
             for i in range(0, N, 20000))
    print(f"  {lab:26s} {d1:27.3e} {d2:29.3e}")
print("  (sets compared after sorting; sampled every 20000th draw to keep the print small)")

print("\n  AND THE VERDICT ARRAYS, ELEMENTWISE -- not the counts lane D reports, the ARRAYS:")
print(f"  {'occupied set':26s} {'base != swap':>14s} {'(base!=fneg) XOR (base!=cneg)':>32s}")
for lab, occ in ROWS:
    base = hull(charset(occ, f, c))
    swap = hull(charset(occ, c, f))
    fneg = hull(charset(occ, -f, c))
    cneg = hull(charset(occ, f, -c))
    print(f"  {lab:26s} {int((base != swap).sum()):14d} "
          f"{int(((base != fneg) ^ (base != cneg)).sum()):32d}")
print("  ZERO in every cell of both columns.  The 'f<->c flips 0 of 200000' column of lane D's")
print("  leg 1A is a column of structural zeros, and its 'c->-c flips the SAME 99785' is the")
print("  f->-f column reprinted.  D-02's NEW half is a VACUOUS CONTROL and is void as evidence.")
print("  THE STATEMENT IT SUPPORTS SURVIVES AS A THEOREM -- the criterion is conjugation-blind,")
print("  one line, no draws -- but it was never a measurement and it isolates nothing about")
print("  class occupancy, which is the variable D-02 claims to be reading.")

print("\n"+"="*100)
print("== R2-B  LEG 3C's TWO 'DIFFERENT' STATES ARE THE SAME ARRAY ==")
print("="*100)
print("  Sealed source, LANE_W10_D_SCOPE/w10d_3_criterion.py, leg 3C:")
print('    ("B0b, state on ONE class-00 vertex   p=(1,0,0,0)", np.array([1., 0, 0, 0])),')
print('    ("B0b, state on ALL FOUR class-00 vtcs p=(1,0,0,0)", np.array([1., 0, 0, 0])),')
a1 = np.array([1., 0, 0, 0])
a2 = np.array([1., 0, 0, 0])
print(f"  np.array_equal(arm1, arm2) = {np.array_equal(a1, a2)}   arm1.tobytes() == arm2.tobytes() = "
      f"{a1.tobytes() == a2.tobytes()}")
print("  BYTE-IDENTICAL.  Two rows of the printed table, one input.  This is verbatim the defect")
print("  W-08's isolation audit named as the commonest FATAL one and which this lane's own")
print("  PUBLISHED_CONVENTIONS promises every script will print a diff to prevent.  Leg 3C")
print("  prints no arm diff.")
print("\n  What a NON-vacuous version would have looked like -- vertex-level states that are")
print("  genuinely different objects and only coincide after the class pushforward:")
print("  B0b's four class-00 vertices are v2, v5, v7, v8 (leg 5A's reconstruction, verified).")
STATES = [("delta on v2                     ", np.array([0, 0, 1., 0, 0, 0, 0, 0, 0])),
          ("delta on v8                     ", np.array([0, 0, 0, 0, 0, 0, 0, 0, 1.])),
          ("uniform on v2,v5,v7,v8          ", np.array([0, 0, .25, 0, 0, .25, 0, .25, .25])),
          ("skew (.7,.1,.15,.05) on the four", np.array([0, 0, .7, 0, 0, .1, 0, .15, .05]))]
CLASS_OF = {0: '11', 3: '11', 1: '10', 4: '10', 6: '01', 2: '00', 5: '00', 7: '00', 8: '00'}
for lab, sv in STATES:
    push = np.array([sum(sv[v] for v in range(9) if CLASS_OF[v] == cl)
                     for cl in ('00', '10', '01', '11')])
    print(f"    {lab}  vertex vector = {np.array2string(sv, precision=2)}"
          f"  ->  class weights {np.array2string(push, precision=2)}")
print(f"    pairwise arm diff of the four VERTEX vectors: "
      f"{min(np.abs(STATES[i][1]-STATES[j][1]).max() for i in range(4) for j in range(i+1,4)):.3f} "
      f"(all four genuinely distinct)")
print("  All four push forward to p = (1,0,0,0) and therefore to rate 0 BY IDENTITY, because the")
print("  functional is a function of the class weights and of nothing else -- W-03's own ruling.")
print("  D-05's CLAIM ('|S| = 1 never fires' is the theorem; root-vs-class coincide only on K1)")
print("  IS TRUE AND IS NOT NEW: REGISTER_V001.md:126 already states")
print('     "|S|=1 -> never, which recovers W-01\'s \'the root can never fire\' as a special case."')
print("  W-02 registered the class-level form of the theorem eight rows before this lane ran.")
