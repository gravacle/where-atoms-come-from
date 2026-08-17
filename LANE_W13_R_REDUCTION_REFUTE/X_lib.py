#!/usr/bin/env python3
"""
X_lib — INDEPENDENT machinery for the W-13/R REFUTER.  Written from scratch; the target
lane's R_lib.py is NOT imported anywhere in this directory, so every number here is a
second implementation and a disagreement with the lane would be a finding.

Contents:
  * K1's registered pi, P, and m(P) by two independent quadratures (Jensen max-form, and
    the root-product form of the y-resultant) so the target's L3 choice is checked, not reused.
  * an exact three-chunk phase reducer for a Fraction rotation number (my own; same IDEA as
    M1_06's An/dA split, which is the corpus's published method, but written here).
  * the relation lattice L(omega) by exhaustive search over a box, with EXACT arithmetic.
  * m() of a one-variable complex Laurent polynomial by roots.
"""
import hashlib
import math
from fractions import Fraction
import numpy as np

PI_K1 = (0.0, 0.3, 0.3, 0.4)          # (p00, p10, p01, p11)


def P_eval(pi, x, y):
    p00, p10, p01, p11 = pi
    return p00 + p10 * x + p01 * y + p11 * x * y


def m_maxform(pi, nq):
    """m(P) = (1/2pi) INT log max(|p00+p10 e^{it}|, |p01+p11 e^{it}|) dt."""
    p00, p10, p01, p11 = pi
    t = np.arange(nq) * (2.0 * np.pi / nq)
    e = np.exp(1j * t)
    return float(np.mean(np.log(np.maximum(np.abs(p00 + p10 * e), np.abs(p01 + p11 * e)))))


def m_one_var(coef):
    """m of sum_i coef[i] z^i (complex allowed), by Jensen on the roots."""
    c = np.asarray(coef, dtype=complex)
    nz = np.nonzero(np.abs(c) > 0)[0]
    c = c[nz[0]:nz[-1] + 1]
    if len(c) == 1:
        return float(np.log(abs(c[0])))
    r = np.roots(c[::-1])
    return float(np.log(abs(c[-1])) + np.sum(np.log(np.maximum(np.abs(r), 1.0))))


# ---------------------------------------------------------------- exact high-precision reals
def frac_sqrt(n, prec=60):
    """Fraction approximation of sqrt(n) accurate to 10^-prec."""
    S = 10 ** prec
    return Fraction(math.isqrt(n * S * S), S)


def frac_pi(prec=60):
    """Machin, pure integer."""
    S = 10 ** (prec + 20)

    def atan_inv(x):
        total, term, k, xx = 0, S // x, 0, x * x
        while term:
            total += term // (2 * k + 1) if k % 2 == 0 else -(term // (2 * k + 1))
            term //= xx
            k += 1
        return total
    return Fraction(16 * atan_inv(5) - 4 * atan_inv(239), S)


def frac_theta_star(prec=60):
    """theta* = arccos(-2/3)/(2 pi), exact Fraction to 10^-prec.
       arccos(-2/3) = pi - arctan(sqrt5/2).  arctan by argument halving + series."""
    G = 30
    p = prec + G
    S = 10 ** p
    PIs = int(frac_pi(p) * S)
    s5 = math.isqrt(5 * S * S)
    # arctan(s5/(2S)) scaled by S
    t = (s5 * S) // (2 * S)
    halv = 0
    while t > S // 100:
        r = math.isqrt(S * S + t * t)
        t = (t * S) // (S + r)
        halv += 1
    total, term, k, tt = 0, t, 0, (t * t) // S
    while term:
        total += term // (2 * k + 1) if k % 2 == 0 else -(term // (2 * k + 1))
        term = (term * tt) // S
        k += 1
    at = total << halv
    phi = PIs - at
    th = (phi * S) // (2 * PIs)
    return Fraction(th // 10 ** G, 10 ** prec)


class ExactRot:
    """frac(k * alpha) for an EXACT Fraction alpha, k <= 1.6e7, residual ~1e-21.
       Three 39-bit chunks so k*chunk stays inside int64."""

    def __init__(self, alpha, D=2 ** 39):
        r = Fraction(alpha) % 1
        self.D = D
        self.A1 = (r.numerator * D) // r.denominator
        r -= Fraction(self.A1, D)
        self.A2 = (r.numerator * D * D) // r.denominator
        r -= Fraction(self.A2, D * D)
        self.A3 = (r.numerator * D ** 3) // r.denominator
        r -= Fraction(self.A3, D ** 3)
        self.rho = float(r)

    def frac(self, k):
        D = self.D
        return np.mod(((k * self.A1) % D) / D + (k * self.A2) / float(D) ** 2
                      + (k * self.A3) / float(D) ** 3 + k * self.rho, 1.0)


# ---------------------------------------------------------------- the relation lattice
def relation_lattice(alpha, beta, B=40, tol=None):
    """All (m,n) != 0 with |m|,|n| <= B and m*alpha + n*beta in Z.
       alpha, beta are exact Fractions -> the test is EXACT."""
    out = []
    a, b = Fraction(alpha), Fraction(beta)
    for m in range(-B, B + 1):
        for n in range(-B, B + 1):
            if m == 0 and n == 0:
                continue
            if (m * a + n * b).denominator == 1:
                out.append((m, n))
    return out


def arm_hash(v):
    return hashlib.sha256(np.ascontiguousarray(np.asarray(v, dtype=np.float64)).tobytes()).hexdigest()[:24]
