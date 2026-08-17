#!/usr/bin/env python3
"""LANE W-10 B REFUTE 2 (LENS 2, SCOPE) — shared library.  WRITTEN FROM SCRATCH.

Nothing is imported from LANE_W10_B_MULTISET/b_lib.py.  Every evaluator here is an
independent code path, so an agreement with the lane is a corroboration and not a
re-print.  Conventions are the register's and the lane's, restated:

    index 0 = class 00 (in neither loop)   character 1
    index 1 = class 10 (in gamma_F only)   character u = conj(W_F) = e^{-i f}
    index 2 = class 01 (in gamma_C only)   character v = W_C      = e^{+i c}
    index 3 = class 11 (in both loops)     character uv
    P(x,y) = p0 + p1 x + p2 y + p3 xy      Z_k = P(u^k, v^k)      lambda = m(P)

EVALUATORS
  R1  mahler_jensen   mpmath, dps configurable, Jensen in y, integrand split at the
                      branch crossings (found by SOLVING the trig equation, not by the
                      lane's closed-form acos) and at both branch minima.
  R2  mahler_lawton   Lawton's theorem: m(P) = lim_N m(P(x, x^N)) with the one-variable
                      measure from the ROOTS (Jensen's formula).  A completely different
                      algorithm -- root finding, not quadrature.  This is the evaluator
                      that makes the cross-check independent.
  R3  dominated_exact EXACT rational certificate, with the criterion DERIVED here in the
                      sorted form  w_max + w_min >= w_mid1 + w_mid2  and cross-checked
                      against the pairwise test on all three pairings.
  R4  mahler_ergodic  direct average of log|Z_k| at a DIFFERENT irrational connection
                      from the lane's (f = 0.9, c = pi/e), N configurable.
  R5  mahler_subtorus EXACT-from-roots average on the (f,c) = 0.1*(20,11) subtorus.
  R6  rate_order4     CLOSED FORM at S1's published connection (W_F=-1, W_C=-i),
                      derived on the page, in exact rationals.
"""
import itertools
from fractions import Fraction

import numpy as np
import mpmath as mp

LBL = ('00', '10', '01', '11')
PERMS = list(itertools.permutations(range(4)))


# ------------------------------------------------------------------ R1
def mahler_jensen(p, dps=30, pairing='y'):
    """m(p0+p1x+p2y+p3xy) by Jensen in y (or x).  Complex p allowed.

    m = (1/2pi) INT log max(|a+b e^{it}|,|c+d e^{it}|) dt with (a,b),(c,d) the pairing.
    Split points are found by mp.findroot on SA-SB (NOT by a closed-form acos), so the
    split rule is a different code path from the lane's.
    """
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        idx = ((0, 1), (2, 3)) if pairing == 'y' else ((0, 2), (1, 3))
        a, b = (mp.mpmathify(complex(p[i])) for i in idx[0])
        c, d = (mp.mpmathify(complex(p[i])) for i in idx[1])
        ra, rb, rc, rd = abs(a), abs(b), abs(c), abs(d)
        be = (mp.arg(b) - mp.arg(a)) if (ra > 0 and rb > 0) else mp.mpf(0)
        de = (mp.arg(d) - mp.arg(c)) if (rc > 0 and rd > 0) else mp.mpf(0)
        tp = 2 * mp.pi

        def SA(t):
            return (ra - rb) ** 2 + 4 * ra * rb * mp.cos((t + be) / 2) ** 2

        def SB(t):
            return (rc - rd) ** 2 + 4 * rc * rd * mp.cos((t + de) / 2) ** 2

        def g(t):
            return SA(t) - SB(t)

        # split points: sign changes of g on a coarse scan, refined by bisection/findroot
        pts = [mp.mpf(0), tp]
        M = 720
        prev_t = mp.mpf(0)
        prev = g(prev_t)
        for i in range(1, M + 1):
            t = tp * i / M
            cur = g(t)
            if (prev > 0 and cur < 0) or (prev < 0 and cur > 0):
                lo, hi = prev_t, t
                for _ in range(int(3.5 * dps) + 20):
                    mid = (lo + hi) / 2
                    if (g(lo) > 0) == (g(mid) > 0):
                        lo = mid
                    else:
                        hi = mid
                pts.append((lo + hi) / 2)
            prev_t, prev = t, cur
        # branch minima / zeros
        if ra > 0 and rb > 0:
            pts.append(mp.fmod(mp.pi - be, tp) % tp)
        if rc > 0 and rd > 0:
            pts.append(mp.fmod(mp.pi - de, tp) % tp)
        pts = sorted(set(float(x) for x in pts))
        pts = [mp.mpf(x) for x in pts if 0 <= x <= float(tp)]
        if pts[0] > 0:
            pts = [mp.mpf(0)] + pts
        if pts[-1] < tp:
            pts = pts + [tp]

        def f(t):
            big = SA(t)
            s = SB(t)
            if s > big:
                big = s
            return mp.log(big) / 2 if big > 0 else mp.mpf('-inf')

        tot = mp.mpf(0)
        for lo, hi in zip(pts[:-1], pts[1:]):
            if hi - lo > mp.mpf(10) ** (-dps + 4):
                tot += mp.quad(f, [lo, hi])
        return tot / tp
    finally:
        mp.mp.dps = old


# ------------------------------------------------------------------ R2
def mahler_lawton(p, N=400):
    """m(P) = lim_N m(P(x, x^N)) (Lawton 1983).  One-variable measure from ROOTS.

    P(x,x^N) = p0 + p1 x + p2 x^N + p3 x^{N+1}.  m = log|lead| + sum log max(1,|root|).
    Completely different algorithm from any quadrature.  float64 root finding, so this
    is a 1e-8-class check, never a headline digit.
    """
    coef = np.zeros(N + 2, dtype=complex)
    coef[0] += p[0]
    coef[1] += p[1]
    coef[N] += p[2]
    coef[N + 1] += p[3]
    # numpy.roots wants highest degree first
    c = coef[::-1]
    nz = np.nonzero(np.abs(c) > 0)[0]
    c = c[nz[0]:]
    lead = c[0]
    r = np.roots(c)
    return float(np.log(abs(lead)) + np.sum(np.log(np.maximum(1.0, np.abs(r)))))


# ------------------------------------------------------------------ R3
def sorted_domination(p):
    """EXACT.  For REAL NON-NEGATIVE p, decide whether EVERY one of the 24 arrangements
    has a pointwise-dominated Jensen pairing, in which case lambda = log(max p) exactly.

    DERIVATION (on the page, exact rationals).  For the pairing {X,Y}|{Z,W},
        SA - SB  =  (X^2+Y^2-Z^2-W^2) + 2(XY - ZW) cos t ,
    monotone in cos t, so one branch dominates iff its value at cos t = +1 and at
    cos t = -1 have the same sign, i.e. iff sign(X+Y-Z-W) = sign(|X-Y|-|Z-W|).
    With w1>=w2>=w3>=w4 the three pairings give, respectively,
        {w1,w2}|{w3,w4}:  w1-w2 >= w3-w4
        {w1,w3}|{w2,w4}:  w1-w3 >= w2-w4
        {w1,w4}|{w2,w3}:  w1+w4 >= w2+w3
    and ALL THREE ARE THE SAME INEQUALITY  w1 + w4 >= w2 + w3.
    Returns (all_dominated, max_weight, per_pairing_flags).
    """
    q = sorted((Fraction(x) for x in p), reverse=True)
    w1, w2, w3, w4 = q
    flags = []
    for (i, j), (k, l) in (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))):
        X, Y, Z, W = q[i], q[j], q[k], q[l]
        hi = (X + Y) ** 2 - (Z + W) ** 2
        lo = (X - Y) ** 2 - (Z - W) ** 2
        flags.append((hi >= 0 and lo >= 0) or (hi <= 0 and lo <= 0))
    return (all(flags), w1, flags, (w1 + w4 >= w2 + w3))


def dominated_exact_perm(p):
    """Exact per-arrangement test over all 24 permutations; returns (#dominated, values)."""
    n = 0
    vals = set()
    for s in PERMS:
        q = [Fraction(p[s[i]]) for i in range(4)]
        a, b, c, d = q
        hi = (a + b) ** 2 - (c + d) ** 2
        lo = (a - b) ** 2 - (c - d) ** 2
        if hi >= 0 and lo >= 0:
            n += 1
            vals.add(max(a, b))
        elif hi <= 0 and lo <= 0:
            n += 1
            vals.add(max(c, d))
    return n, sorted(vals)


# ------------------------------------------------------------------ R4
def mahler_ergodic(p, f=0.9, c=float(np.pi / np.e), N=4_000_000):
    """Direct average of log|Z_k| at a DIFFERENT irrational connection from the lane's."""
    k = np.arange(1, N + 1, dtype=np.float64)
    u = np.exp(-1j * f * k)
    v = np.exp(1j * c * k)
    Z = p[0] + p[1] * u + p[2] * v + p[3] * u * v
    return float(np.mean(np.log(np.abs(Z))))


# ------------------------------------------------------------------ R5
def mahler_subtorus_resonant(p, dps=40):
    """EXACT-from-roots average of log|Z_k| at S3/S4's headline f=2.0, c=1.1.

    (f,c) = 0.1*(20,11) so u^k = w^{-20}, v^k = w^{11} with w = e^{0.1 i k} dense in the
    circle.  w^{20} Z = p0 w^20 + p1 + p2 w^31 + p3 w^11 =: Q(w), and the average of
    log|Z| over the closure is m(Q), Jensen from the roots.  Independent of the lane's
    implementation (mpmath polyroots at 40 digits here).
    """
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        coef = [mp.mpf(0)] * 32          # index = exponent
        coef[20] += mp.mpmathify(complex(p[0]))
        coef[0] += mp.mpmathify(complex(p[1]))
        coef[31] += mp.mpmathify(complex(p[2]))
        coef[11] += mp.mpmathify(complex(p[3]))
        c = coef[::-1]                    # highest first
        while c and c[0] == 0:
            c = c[1:]
        while c and c[-1] == 0:           # factor out w^m: |w|=1 so it contributes 0
            c = c[:-1]
        if not c:
            return mp.mpf('-inf')
        if len(c) == 1:
            return mp.log(abs(c[0]))
        r = mp.polyroots(c, maxsteps=200, extraprec=400)
        tot = mp.log(abs(c[0]))
        for z in r:
            a = abs(z)
            if a > 1:
                tot += mp.log(a)
        return tot
    finally:
        mp.mp.dps = old


# ------------------------------------------------------------------ R6
def rate_order4_exact(p):
    """CLOSED FORM at S1's published connection W_F = -1, W_C = -i (S1 sec6), order 4.

    u = conj(W_F) = -1, v = W_C = -i.  (u^k, v^k) for k=1..4 is
       (-1,-i), (1,-1), (-1, i), (1, 1),
    so with a,b,c,d = p00,p10,p01,p11 the four moduli-squared are
       |a - b - i(c - d)|^2 = (a-b)^2 + (c-d)^2      (k=1 and k=3, complex conjugates)
       (a + b - c - d)^2                             (k=2)
       (a + b + c + d)^2                             (k=4)
    hence   rate = (1/4)[ log((a-b)^2+(c-d)^2) + log|a+b-c-d| + log(a+b+c+d) ].
    Exact in rationals up to the final logs.  Depends on p ONLY through the partition
    {{00,10},{01,11}} -- the Jensen pairing -- so its invariance group is that
    partition's stabiliser, of order 8.
    """
    a, b, c, d = (Fraction(x) for x in p)
    T1 = (a - b) ** 2 + (c - d) ** 2
    T2 = abs(a + b - c - d)
    T3 = a + b + c + d
    return T1, T2, T3


def rate_order4(p, dps=30):
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        T1, T2, T3 = rate_order4_exact(p)
        if T1 == 0 or T2 == 0 or T3 == 0:
            return mp.mpf('-inf')
        F = lambda q: mp.mpf(q.numerator) / mp.mpf(q.denominator)
        return (mp.log(F(T1)) + mp.log(F(T2)) + mp.log(F(T3))) / 4
    finally:
        mp.mp.dps = old


# ------------------------------------------------------------------ helpers
def apply_perm(p, s):
    return tuple(p[s[i]] for i in range(4))


def distinct_arrays(p):
    return len({apply_perm(tuple(p), s) for s in PERMS})


def cyc(s):
    seen = [False] * 4
    out = []
    for i in range(4):
        if seen[i] or s[i] == i:
            seen[i] = True
            continue
        cy, j = [], i
        while not seen[j]:
            seen[j] = True
            cy.append(LBL[j])
            j = s[j]
        out.append('(' + ' '.join(cy) + ')')
    return ''.join(out) or 'e'


def compose(s, t):
    return tuple(s[t[i]] for i in range(4))


def is_subgroup(S):
    S = set(S)
    if tuple(range(4)) not in S:
        return False
    return all(compose(s, t) in S for s in S for t in S)


def hdr(t):
    print()
    print('=' * 96)
    print(t)
    print('=' * 96)
