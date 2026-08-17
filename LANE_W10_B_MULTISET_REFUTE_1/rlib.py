#!/usr/bin/env python3
"""LANE W-10 B REFUTE 1 — shared library.  WRITTEN FROM SCRATCH, NOT IMPORTED FROM THE
LANE UNDER ATTACK.  Nothing in this file imports b_lib.py.

CLASS ORDER (same as the lane's, so the arrays are comparable at the digit):
    index 0 = class 00  character 1     index 1 = class 10  character u = conj(W_F)
    index 2 = class 01  character v     index 3 = class 11  character uv
    P(x,y) = p0 + p1 x + p2 y + p3 xy      lambda = m(P)

THE POINT OF THIS FILE.  Every number in the lane's LEG B (the complex-coefficient
experiment, on which its whole attribution finding rests) is produced by ONE reduction:
m(P) = (1/2pi) INT log max(|A|,|B|) dt.  Its four "independent evaluators" are
    E1 = that reduction in y,  E2 = that reduction in x,  E5 = the float64 twin of E1,
    E3 = exact, but only for REAL NON-NEGATIVE dominated arrays,
    E4 = ergodic, but LEG A only, i.e. real arrays only.
So the complex half of the lane has NO evaluator that does not assume the Jensen
reduction.  This file supplies two that do not:

  R1  m_grid   : the double integral itself, midpoint tensor grid, Richardson in 1/n^2.
                 No Jensen, no branch analysis, no splitting.  Complex coefficients are
                 nothing special to it.
  R2  m_erg    : the corpus's OWN direct method (S4 sec4.2 schedule-B) run on COMPLEX
                 arrays, which the lane never did:  (1/N) sum_k log|P(u^k,v^k)|.
  R3  m_split  : an independent Jensen-family evaluator built on a DIFFERENT
                 decomposition from the lane's --
                     m = (1/2)[log max(a0,a1) + log max(a2,a3)] + (1/8pi) INT |log(SA/SB)|
                 -- whose smooth half is exact (Jensen on each branch separately) and
                 whose remainder is the only quadrature.  Shares the reduction, not the
                 code.
  R4  m_finite : average of log|Z_k| over a FINITE orbit (exact in Fractions where asked).
  R5  m_subtorus: Mahler measure of a one-variable polynomial from its ROOTS (numpy.roots),
                 used for the resonant connection.

PRECISION.  R1/R2/R4/R5 are float64; R3 is mpmath at dps 30.  Every claim below 1e-10 is
R3 or exact rational arithmetic, never R1/R2.
"""
import itertools
from fractions import Fraction

import numpy as np
import mpmath as mp

LBL = ('00', '10', '01', '11')
PERMS = list(itertools.permutations(range(4)))
MATCHINGS = (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2)))
MATCH_NAME = {0: '{00,10}|{01,11}', 1: '{00,01}|{10,11}', 2: '{00,11}|{10,01}'}


# ------------------------------------------------------------------ R1: raw 2D integral
def m_grid(p, n=2048, chunk=256):
    """(1/(2pi)^2) INT INT log|P(e^{is},e^{it})| ds dt by a midpoint tensor grid.
    NO Jensen reduction anywhere.  Midpoint offset keeps the grid off the (1,1) corner
    for real non-negative arrays but NOT off a general zero; the zeros of log are
    integrable and the midpoint rule converges like C/n^2 with a log correction, so
    Richardson in 1/n^2 is applied by m_grid_rich."""
    p = np.asarray(p, dtype=complex)
    s = (np.arange(n) + 0.5) * 2 * np.pi / n
    x = np.exp(1j * s)
    tot = 0.0
    for lo in range(0, n, chunk):
        hi = min(lo + chunk, n)
        y = np.exp(1j * s[lo:hi])[:, None]
        val = (p[0] + p[1] * x[None, :]) + y * (p[2] + p[3] * x[None, :])
        tot += np.log(np.abs(val)).sum()
    return float(tot / (n * n))


def m_grid_rich(p, n=1024):
    """Richardson: two grids n and 2n, extrapolate assuming error ~ C/n^2."""
    a = m_grid(p, n)
    b = m_grid(p, 2 * n)
    return (4 * b - a) / 3.0, a, b


# ------------------------------------------------------------------ R2: ergodic, complex
def m_erg(p, f=1.0, c=float(np.sqrt(2.0)), N=2_000_000, chunk=500_000):
    """The corpus's own direct method, run on COMPLEX coefficients.
    u = e^{-i f k}, v = e^{i c k}, exactly the lane's Z_of convention."""
    p = np.asarray(p, dtype=complex)
    tot = 0.0
    done = 0
    while done < N:
        m = min(chunk, N - done)
        k = np.arange(done + 1, done + m + 1, dtype=np.float64)
        u = np.exp(-1j * f * k)
        v = np.exp(1j * c * k)
        tot += np.log(np.abs(p[0] + p[1] * u + p[2] * v + p[3] * u * v)).sum()
        done += m
    return float(tot / N)


# ------------------------------------------------------------------ R3: split Jensen
def m_split(p, dps=30):
    """m = (1/2)[log max(a0,a1) + log max(a2,a3)] + (1/8pi) INT |log(SA/SB)| dt.
    The first bracket is EXACT (Jensen on each branch alone).  Only the |log| remainder
    is quadratured, and it is split at every crossing and every branch zero."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        q = [mp.mpmathify(complex(z)) for z in p]
        a0, a1, a2, a3 = (abs(z) for z in q)
        beta = (mp.arg(q[1]) - mp.arg(q[0])) if (a0 > 0 and a1 > 0) else mp.mpf(0)
        delta = (mp.arg(q[3]) - mp.arg(q[2])) if (a2 > 0 and a3 > 0) else mp.mpf(0)
        TP = 2 * mp.pi

        def SA(t):
            return (a0 - a1) ** 2 + 4 * a0 * a1 * mp.cos((t + beta) / 2) ** 2

        def SB(t):
            return (a2 - a3) ** 2 + 4 * a2 * a3 * mp.cos((t + delta) / 2) ** 2

        # crossings of SA and SB, and the two branch minima, as split points
        C0 = a0 ** 2 + a1 ** 2 - a2 ** 2 - a3 ** 2
        K = 2 * a0 * a1 * mp.expj(beta) - 2 * a2 * a3 * mp.expj(delta)
        pts = [mp.mpf(0), TP]
        if abs(K) > 0 and abs(C0) <= abs(K):
            ph = mp.acos(-C0 / abs(K))
            ps = mp.arg(K)
            for sgn in (ph, -ph):
                t = mp.fmod(sgn - ps, TP)
                pts.append(t + TP if t < 0 else t)
        for sh in (beta, delta):
            t = mp.fmod(mp.pi - sh, TP)
            pts.append(t + TP if t < 0 else t)
        pts = sorted(pts)

        def g(t):
            sa, sb = SA(t), SB(t)
            if sa <= 0 or sb <= 0:
                return mp.mpf('+inf')
            return abs(mp.log(sa) - mp.log(sb))

        tot = mp.mpf(0)
        for lo, hi in zip(pts[:-1], pts[1:]):
            if hi - lo > mp.mpf(10) ** (-dps + 4):
                tot += mp.quad(g, [lo, hi])
        smooth = (mp.log(max(a0, a1)) + mp.log(max(a2, a3))) / 2
        return smooth + tot / (8 * mp.pi)
    finally:
        mp.mp.dps = old


# ------------------------------------------------------------------ R4: finite orbit
def m_finite(p, f, c, K):
    """(1/K) sum_{k=1..K} log|Z_k| -- exact when the connection has order K."""
    p = np.asarray(p, dtype=complex)
    k = np.arange(1, K + 1)
    u = np.exp(-1j * f * k)
    v = np.exp(1j * c * k)
    return float(np.mean(np.log(np.abs(p[0] + p[1] * u + p[2] * v + p[3] * u * v))))


# ------------------------------------------------------------------ R5: subtorus
def m_poly_roots(coeffs_by_exp):
    """Logarithmic Mahler measure of a one-variable polynomial given as {exp: coeff},
    from its roots: m = log|lead| + sum log max(1,|root|).  numpy.roots."""
    deg = max(coeffs_by_exp)
    lo = min(coeffs_by_exp)
    c = np.zeros(deg + 1, dtype=complex)
    for e, v in coeffs_by_exp.items():
        c[e] = v
    c = c[lo:]                       # strip the common w^lo factor (modulus 1 on |w|=1)
    poly = c[::-1]                   # numpy wants highest first
    nz = np.nonzero(np.abs(poly) > 0)[0]
    poly = poly[nz[0]:]
    r = np.roots(poly)
    return float(np.log(np.abs(poly[0])) + np.sum(np.log(np.maximum(1.0, np.abs(r)))))


def m_resonant(p):
    """S3/S4's exactly-resonant connection f=2.0, c=1.1 = 0.1*(20,11).
    u^k = e^{-i 2k} = w^{-20}, v^k = e^{i 1.1k} = w^{11} with w = e^{i 0.1 k}, whose orbit
    is dense in the circle.  Z = p0 + p1 w^-20 + p2 w^11 + p3 w^-9; multiply by w^20."""
    return m_poly_roots({20: complex(p[0]), 0: complex(p[1]),
                         31: complex(p[2]), 11: complex(p[3])})


# ------------------------------------------------------------------ structure helpers
def apply_perm(p, s):
    return [p[s[i]] for i in range(4)]


def fluxes(p):
    """the three matching fluxes, as absolute values in [0,pi]."""
    A = [np.angle(complex(z)) for z in p]
    raw = [A[0] + A[1] - A[2] - A[3], A[0] + A[2] - A[1] - A[3], A[0] + A[3] - A[1] - A[2]]
    out = []
    for r in raw:
        r = (r + np.pi) % (2 * np.pi) - np.pi
        out.append(abs(r))
    return out


def matching_of(s):
    """which matching of the ORIGINAL indices lands on the Newton diagonal under s."""
    d = frozenset((s[0], s[3]))
    for i, ((a, b), (c, e)) in enumerate(MATCHINGS):
        if d == frozenset((a, b)) or d == frozenset((c, e)):
            return i
    raise AssertionError


def blocks(vals, tol):
    """group 24 values into blocks by tolerance; returns (sizes, reps, labels)."""
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    lab = [-1] * len(vals)
    reps = []
    for i in order:
        placed = False
        for b, r in enumerate(reps):
            if abs(vals[i] - r) <= tol:
                lab[i] = b
                placed = True
                break
        if not placed:
            reps.append(vals[i])
            lab[i] = len(reps) - 1
    sizes = [lab.count(b) for b in range(len(reps))]
    return sizes, reps, lab


def compose(s, t):
    return tuple(s[t[i]] for i in range(4))


def is_subgroup(S):
    S = set(S)
    if tuple(range(4)) not in S:
        return False
    return all(compose(s, t) in S for s in S for t in S)


def max_collinear(p):
    """size of the largest subset of the four coefficients lying on one line through 0."""
    best = 1
    n = len(p)
    for mask in range(1, 1 << n):
        idx = [i for i in range(n) if mask >> i & 1]
        ok = True
        for i in idx:
            for j in idx:
                if i < j:
                    zi, zj = complex(p[i]), complex(p[j])
                    if abs(zi) > 0 and abs(zj) > 0:
                        if abs((zi.conjugate() * zj).imag) / (abs(zi) * abs(zj)) > 1e-12:
                            ok = False
        if ok:
            best = max(best, len(idx))
    return best


def dominated_pairing(p):
    """True if SOME Jensen pairing of the MODULI is dominated pointwise for EVERY flux,
    i.e. | |a|-|b| | >= |c|+|d| for some split.  This is the correct degeneracy test;
    'r_max >= sum of the other three' is the special case c or d = the two smallest."""
    r = sorted((abs(complex(z)) for z in p), reverse=True)
    a, b, c, d = r
    for (i, j, k, l) in ((0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2)):
        if abs(r[i] - r[j]) >= r[k] + r[l]:
            return True
    return False


def hdr(t):
    print()
    print('=' * 96)
    print(t)
    print('=' * 96)
