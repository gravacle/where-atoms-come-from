#!/usr/bin/env python3
"""
R_lib — shared machinery for LANE W-13 / R.

Contains, and nothing else:
  * K1's registered pi and the polynomial P.
  * m(P) by the Jensen reduction (continuous integrand; see CHOICE LEDGER L3).
  * the exact-integer phase reduction of M1_06 (An/dA split), reused and cited.
  * high-precision constants in pure Python integer fixed-point:
        pi (Machin), sqrt (isqrt), arctan (argument halving),
        and theta_star = arccos(-2/3)/(2 pi).
  * the arms-diff guard.

Precision policy: every reported number is produced by the pure-integer or float64 routines
here.  mpmath is imported ONLY inside cross_check_constants() and never supplies a reported
number.
"""
import hashlib
import numpy as np
from fractions import Fraction

# ------------------------------------------------------------------ K1's registered pi
PI_K1 = (0.0, 0.3, 0.3, 0.4)          # (p00, p10, p01, p11)   REGISTER: K1's ready state


def P_eval(pi, x, y):
    p00, p10, p01, p11 = pi
    return p00 + p10 * x + p01 * y + p11 * x * y


def Zk_from_phases(pi, fa, fb):
    """|Z_k| given fractional phases fa = frac(k alpha), fb = frac(k beta)."""
    x = np.exp(2j * np.pi * fa)
    y = np.exp(2j * np.pi * fb)
    return np.abs(P_eval(pi, x, y))


# ------------------------------------------------------------------ m(P), Jensen reduction
def m_jensen(pi, nq):
    """m(P) = (1/2pi) INT log max(|p00+p10 e^{it}|, |p01+p11 e^{it}|) dt, trapezoid on nq nodes.
       Integrand continuous; deterministic; no random input."""
    p00, p10, p01, p11 = pi
    t = np.arange(nq) * (2.0 * np.pi / nq)
    e = np.exp(1j * t)
    a = np.abs(p00 + p10 * e)
    b = np.abs(p01 + p11 * e)
    return float(np.mean(np.log(np.maximum(a, b))))


# ------------------------------------------------------------------ pure-integer high precision
def _isqrt_scaled(n, S):
    """floor(sqrt(n) * S) for integer n, via isqrt."""
    import math
    return math.isqrt(n * S * S)


def hp_pi(prec):
    """floor(pi * 10^prec) - 2 guard digits dropped.  Machin: pi = 16 atan(1/5) - 4 atan(1/239)."""
    S = 10 ** (prec + 20)

    def atan_inv(x):          # atan(1/x) * S, x integer
        total = 0
        term = S // x
        k = 0
        xx = x * x
        while term != 0:
            total += term // (2 * k + 1) if k % 2 == 0 else -(term // (2 * k + 1))
            term //= xx
            k += 1
        return total
    v = 16 * atan_inv(5) - 4 * atan_inv(239)
    return v // 10 ** 20


def hp_sqrt(n, prec):
    """floor(sqrt(n) * 10^prec)."""
    return _isqrt_scaled(n, 10 ** prec)


def hp_atan(num, den, prec):
    """atan(num/den) * 10^prec, num,den positive integers, by argument halving + series.
       atan(t) = 2 atan( t / (1 + sqrt(1+t^2)) )."""
    G = 25
    S = 10 ** (prec + G)
    t = (num * S) // den                       # t * S
    halvings = 0
    while t > S // 100:                        # shrink until t < 0.01
        # t <- t / (1 + sqrt(1+t^2))
        r = _isqrt_scaled((S * S + t * t) // 1, 1)          # sqrt(S^2+t^2) as integer*1
        r = __import__('math').isqrt(S * S + t * t)
        t = (t * S) // (S + r)
        halvings += 1
    # series atan(t) = t - t^3/3 + t^5/5 - ...
    total = 0
    term = t
    tt = (t * t) // S
    k = 0
    while term != 0:
        total += term // (2 * k + 1) if k % 2 == 0 else -(term // (2 * k + 1))
        term = (term * tt) // S
        k += 1
    total <<= halvings
    return total // 10 ** G


def hp_theta_star(prec):
    """theta* = arccos(-2/3) / (2 pi), as a Fraction with denominator 10^prec.
       arccos(-2/3) = pi - arctan(sqrt(5)/2)."""
    G = 30
    p = prec + G
    S = 10 ** p
    PI = hp_pi(p)
    s5 = hp_sqrt(5, p)                       # sqrt(5)*S
    # arctan(sqrt5/2): pass num/den as integers scaled: num = s5, den = 2*S
    at = hp_atan(s5, 2 * S, p)
    phi = PI - at                            # arccos(-2/3) * S
    th = (phi * S) // (2 * PI)               # theta* * S
    return Fraction(th // 10 ** G, 10 ** prec)


def cross_check_constants(prec=60):
    """mpmath cross-check.  Never supplies a reported number."""
    try:
        import mpmath as mp
    except Exception:
        return None
    mp.mp.dps = prec + 20
    ref = mp.acos(mp.mpf(-2) / 3) / (2 * mp.pi)
    mine = hp_theta_star(prec)
    return float(abs(mp.mpf(mine.numerator) / mp.mpf(mine.denominator) - ref))


# ------------------------------------------------------------------ M1_06's exact phase split
class PhaseReducer:
    """frac(k * alpha) for k up to ~1e7 with residual error ~1e-21.
       EXACTLY M1_06's construction (LANE_W08_M1_IDENTIFICATION/M1_06_liouville.py:118-124),
       reused rather than reinvented.  D = 2^39 keeps k*An inside int64 for k <= 1.6e7."""

    def __init__(self, alpha, D=2 ** 39):
        self.D = D
        a = alpha % 1.0
        self.An = int(np.floor(a * D))
        self.dA = a - self.An / D

    def frac(self, k):
        """k : int64 numpy array."""
        return np.mod(((k * self.An) % self.D) / self.D + k * self.dA, 1.0)


# ------------------------------------------------------------------ arms-diff guard
def arm_hash(vec):
    return hashlib.sha256(np.ascontiguousarray(np.asarray(vec, dtype=np.float64)).tobytes()).hexdigest()[:24]


def diff_arms(name_a, va, name_b, vb):
    """Print the guard.  Hashes OUTPUTS (W-10 N-6: a guard that hashes inputs is void)."""
    ha, hb = arm_hash(va), arm_hash(vb)
    n = min(len(va), len(vb))
    d = float(np.max(np.abs(np.asarray(va)[:n] - np.asarray(vb)[:n])))
    verdict = "VOID (arms byte-identical)" if ha == hb else "arms differ"
    print("   ARMS-DIFF GUARD  %-28s %s" % (name_a, ha))
    print("                    %-28s %s" % (name_b, hb))
    print("                    max|A-B| over %d terms = %.6e   -> %s" % (n, d, verdict))
    return ha != hb


class HighPhaseReducer:
    """frac(k*alpha) for an EXACT Fraction alpha, k <= 1.6e7, residual error ~2e-21.
       Three 39-bit chunks: alpha = A1/D + A2/D^2 + A3/D^3 + rho, D = 2^39, |rho| < D^-3.
       k*A_i stays inside int64 for k <= 1.6e7.  Removes the ~1e-10 phase drift that a
       float64 alpha would carry at k = 1e7."""

    def __init__(self, alpha_frac, D=2 ** 39):
        from fractions import Fraction as _F
        r = _F(alpha_frac) % 1
        self.D = D
        self.A1 = (r.numerator * D) // r.denominator
        r = r - _F(self.A1, D)
        self.A2 = (r.numerator * D * D) // r.denominator
        r = r - _F(self.A2, D * D)
        self.A3 = (r.numerator * D * D * D) // r.denominator
        self.rho = float(r - _F(self.A3, D * D * D))

    def frac(self, k):
        D = self.D
        t1 = ((k * self.A1) % D) / D
        t2 = (k * self.A2) / float(D) ** 2
        t3 = (k * self.A3) / float(D) ** 3
        return np.mod(t1 + t2 + t3 + k * self.rho, 1.0)
