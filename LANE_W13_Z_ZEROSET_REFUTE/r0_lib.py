"""LANE W-13 / Z_REFUTE -- shared library.  WRITTEN FROM SCRATCH; it does NOT import z0_lib.
Where a lane-Z routine is needed for a cross-check it is imported EXPLICITLY at the call site and
labelled there, so that no figure in this refuter silently inherits the lane's code."""
import math
import numpy as np
from fractions import Fraction as F

# ---------------------------------------------------------------- ROUTE 2: the S_4-invariant form
def pairing_triple(p):
    """The three 2+2 pairing differences, EXACT.  T_A T_B T_C = (S1^2-S2^2)(D1^2-D2^2)."""
    p00, p10, p01, p11 = p
    return (p00 + p10 - p01 - p11,
            p00 + p11 - p10 - p01,
            p00 + p01 - p10 - p11)

def has_torus_zero_invariant(p):
    """EXACT.  True iff P has a zero on T^2.  Manifestly invariant under S_4 on the weights."""
    a, b, c = pairing_triple(p)
    return a * b * c <= 0

# ---------------------------------------------------------------- ROUTE 3: sorted, pairing-free
def strat_sorted(p):
    """EXACT stratum from the SORTED multiset alone.  Shares no line with the lane's classifier."""
    w1, w2, w3, w4 = sorted([F(q) for q in p], reverse=True)
    DA = w1 + w2 - w3 - w4          # >= 0 always
    DB = w1 + w3 - w2 - w4          # >= 0 always
    DC = w1 + w4 - w2 - w3          # sign decides
    if DB == 0:                     # w1 == w2 and w3 == w4  <=> multiset {a,a,b,b}
        return 'CURVE'
    if DA == 0:
        return 'ONE'
    if DC > 0:
        return 'EMPTY'
    if DC < 0:
        return 'TWO'
    return 'ONE'

# ---------------------------------------------------------------- ROUTE 1: the OTHER Jensen group
def strat_ygrouping(p):
    """EXACT.  Same geometry read with P = A'(y) + x B'(y):  A' = p00+p01 y, B' = p10+p11 y.
    A DIFFERENT pairing of the four weights, so a different-looking predicate."""
    p00, p10, p01, p11 = p
    S1, S2 = p00 + p01, p10 + p11
    D1, D2 = abs(p00 - p01), abs(p10 - p11)
    cI = (p00 == p01) and (p10 == p11)      # A',B' share the zero y = -1
    cII = (p00 == p10) and (p01 == p11)     # P = (1+x)(p00 + p01 y)
    cIII = (p00 == p11) and (p01 == p10)    # |A'| == |B'| identically
    if cI or cII or cIII:
        return 'CURVE'
    pr = (S1 - S2) * (D1 - D2)
    return 'EMPTY' if pr > 0 else ('TWO' if pr < 0 else 'ONE')

# ---------------------------------------------------------------- closed-form zero data (mine)
def zero_angles(p):
    """Angles s0 (for x0 = e^{i s0}) of the isolated zeros, EXACT cos, float sin.
    Returns [] if none, [s0] if one, [s0, -s0] if two, None if the zero set is a curve."""
    if strat_sorted(p) == 'CURVE':
        return None
    p00, p10, p01, p11 = [F(q) for q in p]
    C = p00**2 + p10**2 - p01**2 - p11**2
    D = 2 * (p00 * p10 - p01 * p11)
    if D == 0:
        return []
    cs = F(-C, 1) / D
    if cs > 1 or cs < -1:
        return []
    if cs == 1 or cs == -1:
        return [0.0 if cs == 1 else math.pi]
    a = math.acos(float(cs))
    return [a, -a]

def zero_points(p):
    """(x0, y0) pairs as complex, from the angles.  y0 = -A(x0)/B(x0)."""
    ang = zero_angles(p)
    if ang is None or not ang:
        return ang if ang is None else []
    p00, p10, p01, p11 = [float(q) for q in p]
    out = []
    for s in ang:
        x0 = complex(math.cos(s), math.sin(s))
        A = p00 + p10 * x0
        B = p01 + p11 * x0
        out.append((x0, -A / B))
    return out

# ---------------------------------------------------------------- MY closed form for det M
def detM_closedform(p, s0):
    """THEOREM R2 of this refuter:  det M = - Delta * Delta' * sin(s0) / |B(x0)|^2, with
         Delta  = p10 p01 - p00 p11        (Delta = 0  <=>  P FACTORS)
         Delta' = p10 p11 - p00 p01
    Derived in r2(a); no SVD, no numerical differentiation."""
    p00, p10, p01, p11 = [float(q) for q in p]
    x0 = complex(math.cos(s0), math.sin(s0))
    Dl = p10 * p01 - p00 * p11
    Dp = p10 * p11 - p00 * p01
    B = p01 + p11 * x0
    return -Dl * Dp * math.sin(s0) / (abs(B) ** 2)

def detM_svd(p, x0, y0):
    """The lane's route, reproduced here ONLY as the comparison arm for r2(b)."""
    p00, p10, p01, p11 = [float(q) for q in p]
    alpha = x0 * (p10 + p11 * y0)
    beta = y0 * (p01 + p11 * x0)
    return alpha.real * beta.imag - alpha.imag * beta.real

# ---------------------------------------------------------------- evaluators
def Pabs(p, s, t):
    p00, p10, p01, p11 = [float(q) for q in p]
    x = np.exp(1j * np.asarray(s, dtype=float))
    y = np.exp(1j * np.asarray(t, dtype=float))
    return np.abs(p00 + p10 * x + p01 * y + p11 * x * y)

def jensen_branches(p, t):
    p00, p10, p01, p11 = [float(q) for q in p]
    c, s = np.cos(t), np.sin(t)
    return np.hypot(p00 + p10 * c, p10 * s), np.hypot(p01 + p11 * c, p11 * s)

def mahler_1d(p, n=1 << 20):
    t = (np.arange(n) + 0.5) * (2 * np.pi / n)
    a, b = jensen_branches(p, t)
    return float(np.mean(np.log(np.maximum(a, b))))

# ---------------------------------------------------------------- ROUTE 4: 2-D Newton, no Jensen
def newton_zeros_2d(p, n_starts=4000, seed=1, tol=1e-13, iters=80):
    """Find zeros of P on T^2 by Newton on (Re P, Im P) = 0 in the ANGLES (s,t).
    Uses NO Jensen reduction and NO closed form.  Returns the list of converged (s,t) mod 2pi.
    This is the only route in this lane that can distinguish 'two points' from 'a curve'
    without assuming the reduction."""
    p00, p10, p01, p11 = [float(q) for q in p]
    rng = np.random.default_rng(seed)
    st = rng.random((n_starts, 2)) * 2 * np.pi
    s, t = st[:, 0].copy(), st[:, 1].copy()
    for _ in range(iters):
        x = np.exp(1j * s); y = np.exp(1j * t)
        Pv = p00 + p10 * x + p01 * y + p11 * x * y
        dPs = 1j * (p10 * x + p11 * x * y)          # dP/ds
        dPt = 1j * (p01 * y + p11 * x * y)          # dP/dt
        J = np.stack([np.stack([dPs.real, dPt.real], -1),
                      np.stack([dPs.imag, dPt.imag], -1)], -2)   # 2x2 real Jacobian
        rhs = np.stack([Pv.real, Pv.imag], -1)
        det = J[:, 0, 0] * J[:, 1, 1] - J[:, 0, 1] * J[:, 1, 0]
        ok = np.abs(det) > 1e-14
        inv00 = np.where(ok, J[:, 1, 1] / np.where(ok, det, 1.0), 0.0)
        inv01 = np.where(ok, -J[:, 0, 1] / np.where(ok, det, 1.0), 0.0)
        inv10 = np.where(ok, -J[:, 1, 0] / np.where(ok, det, 1.0), 0.0)
        inv11 = np.where(ok, J[:, 0, 0] / np.where(ok, det, 1.0), 0.0)
        ds = inv00 * rhs[:, 0] + inv01 * rhs[:, 1]
        dt = inv10 * rhs[:, 0] + inv11 * rhs[:, 1]
        step = np.hypot(ds, dt)
        cap = np.minimum(1.0, 0.5 / np.maximum(step, 1e-300))
        s = s - cap * ds; t = t - cap * dt
    res = Pabs(p, s, t)
    keep = res < tol
    return np.stack([s[keep] % (2 * np.pi), t[keep] % (2 * np.pi)], -1)

def count_clusters(pts, eps=1e-6, cap=200):
    """Number of DISTINCT points of T^2 (angle pairs) at resolution eps, with MODULAR distance
    (so a zero sitting at angle 0 is not split by the branch cut).  O(n * cap): a running list
    of representatives, capped -- anything reaching the cap is treated as a positive-dimensional
    solution set.  Returns (count, capped_flag)."""
    P = np.asarray(pts, dtype=float)
    if P.size == 0:
        return 0, False
    reps = []
    TWO = 2 * np.pi
    for q in P:
        new = True
        for r in reps:
            d0 = abs((q[0] - r[0] + np.pi) % TWO - np.pi)
            d1 = abs((q[1] - r[1] + np.pi) % TWO - np.pi)
            if d0 * d0 + d1 * d1 < eps * eps:
                new = False; break
        if new:
            reps.append(q)
            if len(reps) >= cap:
                return cap, True
    return len(reps), False

# ---------------------------------------------------------------- named states
def fr(*a):
    return tuple(F(q) for q in a)

# CORPUS STATES, EACH WITH ITS SOURCE LINE.  Read off S4_THE_MEASUREMENT_V001.md:575-585
# by this refuter, independently of lane Z's NAMED table.  Order is (p00,p10,p01,p11).
S4_575 = [
    ("B0a ring torus, disjoint   {00:2, 01:3, 10:4}",        fr(F(2,9), F(4,9), F(3,9), F(0))),
    ("B0b ring torus, meeting    {00:4, 01:1, 10:2, 11:2}",  fr(F(4,9), F(2,9), F(1,9), F(2,9))),
    ("B3  horn torus             {01:2, 10:2, 11:1}",        fr(F(0), F(2,5), F(2,5), F(1,5))),
    ("B1  K1                     {01:2, 10:2, 11:1}",        fr(F(0), F(2,5), F(2,5), F(1,5))),
    ("B4  spindle                {00:1, 01:1, 10:1, 11:3}",  fr(F(1,6), F(1,6), F(1,6), F(1,2))),
    ("B2  K1 both filled         {01:2, 10:2, 11:1}",        fr(F(0), F(2,5), F(2,5), F(1,5))),
    ("B1p K1-bridged             {01:3, 10:3}",              fr(F(0), F(1,2), F(1,2), F(0))),
    ("B1q K1-bridged + spectator {00:1, 01:3, 10:3}",        fr(F(1,7), F(3,7), F(3,7), F(0))),
    ("B1s K1 subdivided          {01:5, 10:5, 11:1}",        fr(F(0), F(5,11), F(5,11), F(1,11))),
]
K1_REG   = fr(F(0), F(3,10), F(3,10), F(2,5))          # N1 as registered (registrar's brief)
S1_PUB   = fr(F(0), F(0), F(1,2), F(1,2))              # S1 sec6 ready state via M1_08 T1
SENSEC4  = fr(F(1,4), F(1,4), F(1,4), F(1,4))          # S4 sense C, four classes
B0b_LANE = fr(F(4,9), F(2,9), F(2,9), F(1,9))          # what LANE Z used for B0b
B0b_S4   = fr(F(4,9), F(2,9), F(1,9), F(2,9))          # what S4:575 actually says
CENTROID = fr(F(0), F(1,3), F(1,3), F(1,3))            # M1_06's counterexample state
