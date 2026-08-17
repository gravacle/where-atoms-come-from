#!/usr/bin/env python3
"""LANE W-10 B — shared library.  Everything needed to re-run is here or in CONVENTIONS.txt.

CONVENTIONS (canonical order of the four vertex classes, used EVERYWHERE below):

    index 0 = class 00  (v in neither loop)      character  1
    index 1 = class 10  (v in gamma_F only)      character  u = conj(W_F)
    index 2 = class 01  (v in gamma_C only)      character  v = W_C
    index 3 = class 11  (v in both loops)        character  uv

    P(x,y) = p[0] + p[1] x + p[2] y + p[3] xy          Z_k = P(u^k, v^k)
    lambda = m(P) = (1/(2pi)^2) INT INT log|P(e^{is},e^{it})| ds dt   (logarithmic Mahler measure)

JENSEN REDUCTION.  P is linear in y, so integrating y first (Jensen) gives
    m(P) = (1/2pi) INT log max(|A(e^{it})|, |B(e^{it})|) dt,   A = p0 + p1 x,  B = p2 + p3 x
and P is linear in x too, giving the SECOND, independent reduction
    m(P) = (1/2pi) INT log max(|A'(e^{it})|,|B'(e^{it})|) dt,  A' = p0 + p2 y, B' = p1 + p3 y.
Both are valid for COMPLEX coefficients.  They are different code paths on the same P.

EVALUATORS PROVIDED (four, deliberately independent):
    E1  m_jensen(p,'y')    mpmath, dps 30, integrand split analytically at the branch
                           crossings and at the branch zeros, tanh-sinh on each piece
    E2  m_jensen(p,'x')    same machinery, the OTHER pairing
    E3  m_dominated_exact  EXACT.  When one Jensen branch dominates the other pointwise
                           the integral collapses to log(max of that branch's two
                           coefficients); the domination test is done in exact rationals
    E4  m_ergodic          the corpus's own method (S4 sec4.2 "direct schedule-B
                           simulation"): (1/N) sum_{k<=N} log|P(u^k,v^k)| at a fixed
                           irrational connection.  No quadrature anywhere in it.

PRECISION.  E1/E2 run at mpmath dps=30 and are reported to 18 places.  E4 is numpy
float64.  E3 is exact (fractions.Fraction) up to one final log of a rational.
"""
import itertools
from fractions import Fraction

import numpy as np
import mpmath as mp

LBL = ('00', '10', '01', '11')
PAIRINGS = {'y': ((0, 1), (2, 3)), 'x': ((0, 2), (1, 3))}


# --------------------------------------------------------------------------- E1/E2
def m_jensen(p, pairing='y', dps=30, extra_split=()):
    """Logarithmic Mahler measure of p0 + p1 x + p2 y + p3 xy.  Complex p allowed."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        (i, j), (k, l) = PAIRINGS[pairing]
        a, b, c, d = (mp.mpmathify(complex(p[t])) for t in (i, j, k, l))
        ra, rb, rc, rd = (abs(z) for z in (a, b, c, d))
        beta = (mp.arg(b) - mp.arg(a)) if (ra > 0 and rb > 0) else mp.mpf(0)
        delta = (mp.arg(d) - mp.arg(c)) if (rc > 0 and rd > 0) else mp.mpf(0)
        two_pi = 2 * mp.pi

        # |a + b e^{i th}|^2 = (ra-rb)^2 + 4 ra rb cos^2(th/2).  This form is used
        # instead of ra^2+rb^2+2 ra rb cos(th) because the latter loses ALL significant
        # digits near a branch zero (1 + cos(th) with cos(th) rounded to -1), which is
        # exactly where the integrand matters.  Caught on the first run: SENSE C
        # (1/4,1/4,1/4,1/4) returned -inf under the unstable form.
        def SA(t):
            return (ra - rb) ** 2 + 4 * ra * rb * mp.cos((t + beta) / 2) ** 2

        def SB(t):
            return (rc - rd) ** 2 + 4 * rc * rd * mp.cos((t + delta) / 2) ** 2

        # SA - SB = C0 + Re(K e^{it}) = C0 + |K| cos(t + arg K)
        C0 = ra ** 2 + rb ** 2 - rc ** 2 - rd ** 2
        K = 2 * ra * rb * mp.expj(beta) - 2 * rc * rd * mp.expj(delta)
        pts = [mp.mpf(0), two_pi]
        if abs(K) > 0 and abs(C0) <= abs(K):
            phi0 = mp.acos(-C0 / abs(K))
            psi = mp.arg(K)
            for s in (phi0, -phi0):
                t = mp.fmod(s - psi, two_pi)
                if t < 0:
                    t += two_pi
                pts.append(t)
        # zeros of a branch (log singularity) -- put them on a split point
        if ra == rb and ra > 0:
            t = mp.fmod(mp.pi - beta, two_pi)
            pts.append(t + two_pi if t < 0 else t)
        if rc == rd and rc > 0:
            t = mp.fmod(mp.pi - delta, two_pi)
            pts.append(t + two_pi if t < 0 else t)
        for t in extra_split:
            pts.append(mp.mpf(t))
        pts = sorted(pts)

        def f(t):
            m2 = SA(t)
            s2 = SB(t)
            big = m2 if m2 > s2 else s2
            return mp.log(big) / 2 if big > 0 else mp.mpf('-inf')

        tot = mp.mpf(0)
        for lo, hi in zip(pts[:-1], pts[1:]):
            if hi - lo > mp.mpf(10) ** (-dps + 4):
                tot += mp.quad(f, [lo, hi])
        return tot / two_pi
    finally:
        mp.mp.dps = old


# --------------------------------------------------------------------------- E5
_GLX, _GLW = np.polynomial.legendre.leggauss(160)


def m_fast(p, pairing='y'):
    """float64 twin of m_jensen: identical analytic split, Gauss-Legendre per piece.
    Accurate to ~1e-14 when neither branch vanishes on the circle.  Used only where a
    1e-9 decision is being made (telling 0.01 apart from 0), never for a headline digit."""
    (i, j), (k, l) = PAIRINGS[pairing]
    a, b, c, d = (complex(p[t]) for t in (i, j, k, l))
    ra, rb, rc, rd = abs(a), abs(b), abs(c), abs(d)
    beta = (np.angle(b) - np.angle(a)) if (ra > 0 and rb > 0) else 0.0
    delta = (np.angle(d) - np.angle(c)) if (rc > 0 and rd > 0) else 0.0
    C0 = ra ** 2 + rb ** 2 - rc ** 2 - rd ** 2
    K = 2 * ra * rb * np.exp(1j * beta) - 2 * rc * rd * np.exp(1j * delta)
    pts = [0.0, 2 * np.pi]
    if abs(K) > 0 and abs(C0) <= abs(K):
        phi0 = np.arccos(-C0 / abs(K))
        psi = np.angle(K)
        for s in (phi0, -phi0):
            pts.append((s - psi) % (2 * np.pi))
    # always split at each branch's minimum: that is where a zero or a near-zero of
    # that branch sits, and Gauss-Legendre must not straddle it
    pts.append((np.pi - beta) % (2 * np.pi))
    pts.append((np.pi - delta) % (2 * np.pi))
    pts = sorted(set(pts))
    tot = 0.0
    for lo, hi in zip(pts[:-1], pts[1:]):
        if hi - lo <= 0:
            continue
        t = 0.5 * (hi - lo) * _GLX + 0.5 * (hi + lo)
        SA = (ra - rb) ** 2 + 4 * ra * rb * np.cos((t + beta) / 2) ** 2
        SB = (rc - rd) ** 2 + 4 * rc * rd * np.cos((t + delta) / 2) ** 2
        tot += 0.5 * (hi - lo) * np.dot(_GLW, 0.5 * np.log(np.maximum(SA, SB)))
    return tot / (2 * np.pi)


# --------------------------------------------------------------------------- E3
def m_dominated_exact(p, pairing='y'):
    """If one Jensen branch dominates pointwise, return (True, Fraction M) with
    lambda = log(M) EXACTLY; else (False, None).  Requires real non-negative rational p.

    A = p_i + p_j x, B = p_k + p_l x.  SA - SB at cos t = +1 is (S1^2 - S2^2), at
    cos t = -1 it is (D1^2 - D2^2), and SA - SB is monotone in cos t, so A dominates
    iff both are >= 0 and B dominates iff both are <= 0.  All in exact rationals.
    """
    (i, j), (k, l) = PAIRINGS[pairing]
    a, b, c, d = (Fraction(p[t]) for t in (i, j, k, l))
    if min(a, b, c, d) < 0:
        return (False, None)
    S1, S2, D1, D2 = a + b, c + d, a - b, c - d
    hi = S1 * S1 - S2 * S2
    lo = D1 * D1 - D2 * D2
    if hi >= 0 and lo >= 0:
        return (True, max(a, b))
    if hi <= 0 and lo <= 0:
        return (True, max(c, d))
    return (False, None)


# --------------------------------------------------------------------------- E4
def m_ergodic(p, f=1.0, c=float(np.sqrt(2.0)), N=2_000_000):
    """The corpus's own direct method (S4 sec4.2).  No quadrature.  float64."""
    k = np.arange(1, N + 1, dtype=np.float64)
    u = np.exp(-1j * f * k)          # u = conj(W_F);  W_F = e^{if}
    v = np.exp(1j * c * k)           # v = W_C
    Z = p[0] + p[1] * u + p[2] * v + p[3] * u * v
    return float(np.mean(np.log(np.abs(Z))))


# --------------------------------------------------------------------------- Z_k
def Z_of(p, f, c, ks):
    """Z_k = sum_ab p_ab (u^a v^b)^k  with u = e^{-if}, v = e^{ic}.  ks an array."""
    ks = np.asarray(ks)
    u = np.exp(-1j * f * ks)
    v = np.exp(1j * c * ks)
    return p[0] + p[1] * u + p[2] * v + p[3] * u * v


# --------------------------------------------------------------------------- S4
PERMS = list(itertools.permutations(range(4)))


def apply_perm(p, s):
    """position i of the new array receives the old value p[s[i]]."""
    return type(p)(p[s[i]] for i in range(4)) if isinstance(p, tuple) else np.array([p[s[i]] for i in range(4)])


def matching_of(s):
    """the D4 invariant: which two VALUES land on the Newton square's diagonal."""
    return (frozenset((s[0], s[3])), frozenset((s[1], s[2])))


def matching_key(s):
    return frozenset(matching_of(s))


def cycle_notation(s):
    """s as a permutation of the labels 00,10,01,11 (i -> s[i] reading 'position i takes value s[i]')."""
    seen = [False] * 4
    out = []
    for i in range(4):
        if seen[i] or s[i] == i:
            seen[i] = True
            continue
        cyc = []
        j = i
        while not seen[j]:
            seen[j] = True
            cyc.append(LBL[j])
            j = s[j]
        out.append('(' + ' '.join(cyc) + ')')
    return ''.join(out) if out else 'e'


def compose(s, t):
    """apply s then t:  new[i] = p[s[t[i]]]."""
    return tuple(s[t[i]] for i in range(4))


def is_subgroup(S):
    S = set(S)
    if tuple(range(4)) not in S:
        return False
    for s in S:
        for t in S:
            if compose(s, t) not in S:
                return False
    return True


def collinearity_defect(p):
    """0 iff the four coefficients are real up to ONE common phase (pairwise ratios real).
    Returns max over pairs of |Im(conj(p_i) p_j)| / (|p_i||p_j|)."""
    w = 0.0
    for i in range(4):
        for j in range(i + 1, 4):
            n = abs(p[i]) * abs(p[j])
            if n > 0:
                w = max(w, abs((np.conj(p[i]) * p[j]).imag) / n)
    return w


def flux(p, matching=((0, 3), (1, 2))):
    """arg( prod over first pair / prod over second pair ) in (-pi,pi]."""
    (i, j), (k, l) = matching
    z = p[i] * p[j] / (p[k] * p[l])
    return float(np.angle(z))


def hdr(t):
    print()
    print('=' * 96)
    print(t)
    print('=' * 96)
