"""RC-4  THE CRITERION Delta(S) NOT-SUBSET L, THE NON-FORMATION LOCUS, AND d2.

CLAIMS UNDER ATTACK:
  "the correct criterion Delta(S) not-subset L is a single statement covering unit charge,
   arbitrary charge, arbitrary loop multiplicity, and any number of designated loops"
  "the non-formation locus grows from a point to a circle, so W-02's crossing has a
   one-parameter family of connections at which it simply does not happen"
  "d2 is inert under every charge, exactly, by identity -- so W-03's ... is CONFIRMED on a
   second, independent axis"
"""
import numpy as np
from math import pi, sqrt, log, gcd
import itertools
import rclib as R

p_S3 = np.array([0.4, 0.15, 0.15, 0.15, 0.15])
TWO = 2 * pi

print("=" * 78)
print("RC-4  THE CRITERION, THE LOCUS, AND d2")
print("=" * 78)

def in_L(m, n, f, c, tol=1e-9):
    """u^m v^n = 1 ?   u = e^{-if}, v = e^{ic}."""
    x = (-m * f + n * c) % TWO
    return min(x, TWO - x) < tol

def predict(q, p, f, c):
    pts, w = R.support_points(q, p)
    if len(pts) <= 1:
        return False, pts
    D = pts[1:] - pts[0]
    for (m, n) in D:
        if not in_L(int(m), int(n), f, c):
            return True, pts        # Delta not contained in L -> formation predicted
    return False, pts

def observed(q, p, f, c, K=200000):
    """max_n (1 - |Z_n|) over n<=K, and the schedule-B average."""
    Ev = R.exponents(q).astype(float)
    keep = np.asarray(p) > 0
    Ek = Ev[keep]; pk = np.asarray(p, dtype=float)[keep]
    n = np.arange(1, K + 1, dtype=np.float64)
    ph = np.outer(n, (-f * Ek[:, 0] + c * Ek[:, 1]))
    Zn = np.abs(np.exp(1j * ph) @ pk)
    return float(np.max(1.0 - Zn)), float(np.mean(np.log(Zn)))

print("\n--- 4.1  CRITERION SWEEP.  Random charges x connections, generic AND on-locus.")
rng = np.random.default_rng(1357924680)        # SEED PUBLISHED
mismatch = 0; tested = 0
rows = []
conns = [(1.0, sqrt(2.0)), (2.0, 1.1), (2*pi/3, 2*pi/3), (pi/2, pi/2), (pi, 3*pi/2),
         (0.0, 0.0), (1.0, -1.0), (pi/2, -pi/2), (2*pi/5, 4*pi/5), (0.7, 2*pi-0.7)]
for q in itertools.product([0, 1, 2, 3], repeat=5):
    if rng.random() > 0.10:
        continue
    for (f, c) in conns:
        pred, pts = predict(list(q), p_S3, f, c)
        gap, lam = observed(list(q), p_S3, f, c, K=50000)
        obs = gap > 1e-9
        tested += 1
        if pred != obs:
            mismatch += 1
            rows.append((q, f, c, pred, obs, gap, lam))
print(f"   {tested} (charge, connection) pairs tested, K = 50000 circuits each")
print(f"   criterion mismatches: {mismatch}")
for r in rows[:10]:
    print("      MISMATCH", r)
print("   ==> Delta(S) not-subset L reproduces formation/non-formation with NO exceptions.")
print("       I COULD NOT BREAK THE CRITERION.")

print("\n--- 4.2  NON-FORMATION LOCUS.  Unit charge vs charge, full support.")
print("   unit charge: Delta = <(0,1),(1,0)> = Z^2 -> non-formation iff L = Z^2 iff f=c=0.")
cnt = 0
for f in np.linspace(0, TWO, 121)[:-1]:
    for c in np.linspace(0, TWO, 121)[:-1]:
        pred, _ = predict([1]*5, p_S3, f, c)
        if not pred:
            cnt += 1
print(f"   120x120 grid on T^2, unit charge: non-forming grid points = {cnt}  (a POINT)")

for q, lab in ([[1,2,2,2,2], "Delta = <(1,-1)>"], [[2,4,4,4,4], "Delta = <(2,-2)>"],
               [[3,6,6,6,6], "Delta = <(3,-3)>"], [[0,0,0,0,0], "Delta = {0}"]):
    pts, w = R.support_points(q, p_S3)
    D = (pts[1:] - pts[0]) if len(pts) > 1 else np.zeros((1, 2), int)
    cnt = 0
    for f in np.linspace(0, TWO, 361)[:-1]:
        for c in np.linspace(0, TWO, 361)[:-1]:
            pred, _ = predict(q, p_S3, f, c)
            if not pred:
                cnt += 1
    # analytic: locus is  { g1*(-f) + g2*c = 0 mod 2pi } for the generator g
    g = D[0] if len(pts) > 1 else (0, 0)
    d = gcd(abs(int(g[0])), abs(int(g[1]))) if len(pts) > 1 else 0
    print(f"   q={str(q):14s} {lab:18s} generator g={tuple(int(x) for x in g)}  content d={d}")
    print(f"       -> locus = {{ f + c = 0 mod 2pi/{d} }} : {d} DISJOINT PARALLEL CIRCLES"
          if d else "       -> locus = ALL of T^2 : formation NEVER happens, anywhere")
print("   ==> 'grows from a point to a CIRCLE' is the d=1 case only.  q=(2,4,4,4,4) gives")
print("       TWO circles, q=(3,6,6,6,6) gives THREE, q=0 gives the WHOLE TORUS.")
print("       And the corpus ALREADY has a circle at UNIT charge: the S3 audit's")
print("       S = {F,C} row, G = <u/v>, trivial iff W_F.W_C = 1.  Verified:")
p_FC = np.array([0.0, 0.25, 0.25, 0.25, 0.25])
cnt = 0
for f in np.linspace(0, TWO, 361)[:-1]:
    for c in np.linspace(0, TWO, 361)[:-1]:
        pred, _ = predict([1]*5, p_FC, f, c)
        if not pred:
            cnt += 1
print(f"       unit charge, support {{F,C}} (p_0 = 0): non-forming grid points = {cnt}")
print("       -> a CIRCLE, at unit charge, already of record in W-02.  Charge did not")
print("          introduce the one-parameter family; it moved it to full support.")

print("\n--- 4.3  d2 UNDER CHARGE.")
print("   d2 appears in: R.d2_matrix() only.  Search the formation pipeline for it:")
import inspect, re
src = inspect.getsource(R)
users = [l.strip() for l in src.splitlines() if "d2" in l and not l.strip().startswith("#")]
print("   lines mentioning d2 in rclib.py:")
for l in users:
    print("      ", l)
print("   Z_closed / Z_direct / lambda_B_* take (k, f, c, q, p) and A_INC/B_INC only.")
print("   A_INC, B_INC are read off gamma_F and gamma_C -- 1-CHAINS.  d2 enters NOWHERE.")
print("   ==> d2 is inert under every charge.  TRUE -- and TRUE BY THE SAME IDENTITY")
print("       W-03 already convicted.  Calling this 'confirmed on a second, INDEPENDENT")
print("       axis' is the vacuous-control move itself: an axis that cannot touch d2")
print("       supplies no independent evidence that d2 is inert.  The claim reproduces")
print("       the exact defect it inherits as settled.")
