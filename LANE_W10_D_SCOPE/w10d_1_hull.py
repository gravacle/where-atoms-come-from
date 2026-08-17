# W10-D leg 1 -- W-01's CONVEX HULL CRITERION AND ITS ADVERTISED VIRTUE, ON FOUR-CLASS CARRIERS.
# Verifies W-09's ruling INDEPENDENTLY (different hull algorithm, second algorithm cross-check),
# then goes past it: which sign/exchange symmetries the criterion acquires off K1, and the
# quantifier, tested on a four-class carrier's OWN published weights.
import numpy as np
from itertools import combinations
from math import comb

rng = np.random.default_rng(20260816)
N = 200000
f = rng.uniform(-np.pi, np.pi, N)
c = rng.uniform(-np.pi, np.pi, N)

CHAR = {'00': lambda f, c: np.ones_like(f, dtype=complex),
        '10': lambda f, c: np.exp(-1j*f),
        '01': lambda f, c: np.exp(1j*c),
        '11': lambda f, c: np.exp(1j*(c-f))}

def hull_H1(pts):
    """sorted-angle max-gap test"""
    A = np.sort(np.angle(np.stack(pts, axis=0)), axis=0)
    g = np.diff(np.concatenate([A, A[:1] + 2*np.pi], axis=0), axis=0)
    return g.max(axis=0) <= np.pi + 1e-12

def hull_H2(pts):
    """independent: 0 in conv(S) iff 0 on a segment joining two points, or inside a triangle."""
    n = len(pts)
    P = np.stack(pts, axis=0)
    X, Y = P.real, P.imag
    out = np.zeros(P.shape[1], dtype=bool)
    for i, j in combinations(range(n), 2):          # 0 on segment i-j
        cross = X[i]*Y[j] - Y[i]*X[j]
        dot = X[i]*X[j] + Y[i]*Y[j]
        out |= (np.abs(cross) <= 1e-12) & (dot <= 1e-12)
    for i, j, k in combinations(range(n), 3):       # 0 in triangle ijk
        d1 = X[i]*Y[j] - Y[i]*X[j]                  # cross(Pi,Pj) = cross(Pi-0, Pj-0)
        d2 = X[j]*Y[k] - Y[j]*X[k]
        d3 = X[k]*Y[i] - Y[k]*X[i]
        out |= ((d1 >= -1e-12) & (d2 >= -1e-12) & (d3 >= -1e-12)) | \
               ((d1 <= 1e-12) & (d2 <= 1e-12) & (d3 <= 1e-12))
    return out

def region(occ, ff, cc):
    pts = [CHAR[o](ff, cc) for o in occ]
    h1, h2 = hull_H1(pts), hull_H2(pts)
    assert (h1 == h2).all(), f"HULL ALGORITHMS DISAGREE on {occ}: {(h1!=h2).sum()}"
    return h1

ROWS = [("B1  K1 as handed        ", ('10', '01', '11')),
        ("B1q K1 + spectator      ", ('00', '10', '01')),
        ("B1p K1-bridged          ", ('10', '01')),
        ("B0b ring torus, meeting ", ('00', '10', '01', '11')),
        ("B4  spindle             ", ('00', '10', '01', '11'))]

print("=" * 100)
print("ARM DIFF FIRST.  The corpus's commonest FATAL defect is two byte-identical arms reported")
print("as a confirmation.  These are the character sets actually handed to the SAME evaluator:")
for lab, occ in ROWS:
    print(f"   {lab}  occupied = {occ}")
seen = {}
for lab, occ in ROWS:
    seen.setdefault(occ, []).append(lab.strip())
print("   DISTINCT ARMS:", len(seen), "of", len(ROWS),
      "-- B0b and B4 are THE SAME ARM for this criterion (same occupied set); that is reported")
print("   as one arm below, not as two confirmations.")
print("=" * 100)

print("\n== 1A  FIRING REGION, f->-f, c->-c, f<->c  (200000 draws, seed 20260816) ==")
print(f"{'carrier':26s} {'|S|':>3s} {'fire':>9s} {'f->-f flips':>12s} {'c->-c flips':>12s} {'f<->c flips':>12s}")
res = {}
for lab, occ in ROWS:
    base = region(occ, f, c)
    fneg = region(occ, -f, c)
    cneg = region(occ, f, -c)
    swap = region(occ, c, f)
    res[lab] = base
    print(f"{lab:26s} {len(occ):3d} {base.mean():9.6f} {int((base!=fneg).sum()):12d} "
          f"{int((base!=cneg).sum()):12d} {int((base!=swap).sum()):12d}")

print("\n  READ.  Three classes: the verdict moves under f->-f AND under c->-c, at ~half the")
print("  grid; four classes: 0 of 200000 under either.  The exchange f<->c moves NEITHER, on")
print("  three OR four classes -- so 'distinguishes curvature from flat holonomy' was never a")
print("  distinction between the two ROLES even on K1; it is a sensitivity to the SIGN of each,")
print("  and THAT is what four classes destroys.  Recorded because it narrows W-09's finding.")

print("\n== 1B  THE EXACT VALUES, NOT MEASUREMENTS ==")
p_half3 = 2.0**(-3+1)*sum(comb(2, k) for k in range(2))
print(f"  Wendel(1962), N=3 points, d=2:  P(0 in hull) = 1 - {p_half3} = {1-p_half3}  = 1/4 EXACTLY")
print("     (three occupied classes divide by one character to {1, e^{i a}, e^{i b}} with (a,b)")
print("      jointly uniform on T^2 -- checked below by exhibiting the reduction for each set)")
for lab, occ in ROWS[:2]:
    occ_l = list(occ)
    d = occ_l[0]
    print(f"     {lab.strip():22s} divide by chi_{d}: angles become "
          f"{[f'chi_{o}/chi_{d}' for o in occ_l]}")
cf = (np.cos(f) + np.cos(c) <= 0)
four = res["B0b ring torus, meeting "]
print(f"\n  FOUR CLASSES closed form  0 in conv{{1,u,v,uv}} <=> cos f + cos c <= 0 :"
      f"  agrees on {int((four==cf).sum())} of {N}")
print("     and (f,c) -> (pi-f, pi-c) preserves the uniform measure while flipping the sign of")
print(f"     cos f + cos c, so the region is 1/2 EXACTLY.  measured {four.mean():.6f}")
print("  BOTH VALUES ARE EXACT.  W-09's ruling is REPRODUCED from an independent implementation")
print("  with a second, structurally different hull algorithm agreeing on every draw.")

print("\n== 1C  THE QUANTIFIER, ON A FOUR-CLASS CARRIER'S OWN PUBLISHED WEIGHTS ==")
print("  W-01 is registered at REGISTER_V001.md:43 as 'vanishes IFF 0 lies in the convex hull'.")
print("  HULL is an EXISTENCE statement over ready states; VANISHES is a statement about ONE.")
W = {'B0b': np.array([4/9, 2/9, 1/9, 2/9]), 'B4': np.array([1/6, 1/6, 1/6, 3/6]),
     'B1 ': np.array([0.0, 2/5, 2/5, 1/5]), 'SENSE-C': np.array([.25, .25, .25, .25])}
for nm, p in W.items():
    Z1 = p[0] + p[1]*np.exp(-1j*f) + p[2]*np.exp(1j*c) + p[3]*np.exp(1j*(c-f))
    occ = tuple(k for k, w in zip(('00', '10', '01', '11'), p) if w > 0)
    h = region(occ, f, c)
    print(f"  {nm:8s} p={np.array2string(p, precision=4)}  hull fires {h.mean():.4f} of the grid;"
          f"  min_grid |Z_1| = {np.abs(Z1).min():.6e}")
print("  On B0b's and B4's OWN weights the hull says 'fires' on half the connections and the")
print("  actual overlap |Z_1| is bounded away from zero on the WHOLE grid.  The gap is not a")
print("  measure-zero technicality off K1: it is half the parameter space.")

print("\n== 1D  WHY -- AND IT IS EXACT, NOT NUMERICAL ==")
print("  P(x,y) = (p00 + p01 y) + x (p10 + p11 y) has a zero on T^2 iff")
print("  |p00 + p01 y| = |p10 + p11 y| for some |y|=1, i.e. iff")
print("     (p00^2+p01^2-p10^2-p11^2) + 2(p00 p01 - p10 p11) cos t = 0  has a root with |cos t|<=1.")
for nm, p in W.items():
    A = p[0]**2 + p[2]**2 - p[1]**2 - p[3]**2
    B = 2*(p[0]*p[2] - p[1]*p[3])
    if abs(B) < 1e-15:
        sol = "no cos t term; " + ("identically zero" if abs(A) < 1e-15 else "no root")
    else:
        r = -A/B
        sol = f"cos t = {r:+.6f}  ->  {'ROOT EXISTS' if abs(r)<=1 else 'NO ROOT'}"
    print(f"  {nm:8s} A={A:+.6f} B={B:+.6f}   {sol}")
print("  B0b and B4 have NO torus zero of P at their own published weights.  The four-class")
print("  'firing region 1/2' is a statement about which (f,c) admit SOME firing state, and on")
print("  these two carriers the state that would fire is NOT the one the carrier hands you.")
print("\n  NOTE AGAINST THE BRIEF I WAS GIVEN: 'a four-class P HAS torus zeros' is FALSE as a")
print("  general statement.  It is true at SENSE-C (1/4,1/4,1/4,1/4), where P factors as")
print("  (1+x)(1+y)/4 and the zero set is a CURVE, and false at both four-class carriers' own")
print("  SENSE-U weights.  Torus zeros are a property of the WEIGHTS, not of the class count.")
