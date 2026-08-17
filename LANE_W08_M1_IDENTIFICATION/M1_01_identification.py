#!/usr/bin/env python3
"""
M1_01 — THE IDENTIFICATION.   Z_k = P(u^k, v^k) ?

Conventions: M1_00_conventions.txt.  Precision: IEEE double unless a block says EXACT.
Seeds: numpy default_rng(20260816) for the float sweep; the exact block uses no randomness
beyond an enumeration.

Tests, in order:
  T1  float sweep over random connections x random ready states x k, both candidate
      conventions u = conj(W_F) and u = W_F.  Worst deviation reported for each.
  T2  EXACT arithmetic (Gaussian rationals): W_F, W_C in the 4th roots of unity, rational
      |s_v|^2.  Deviation is an exact rational; reported as an exact 0 or not.
  T3  gauge invariance of Z_k under the FULL gauge action (connection AND section).
  T4  dressing invariance: Z_k built from the tree-dressed section t equals Z_k built from s.
  T5  the degenerate-carrier check: p00 = 0 identically on K1 over random states.
"""
import numpy as np
from fractions import Fraction
from itertools import product

FACE_V = [0, 1, 2]
CYC_V  = [0, 3, 4]
EDGES  = [(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 0)]   # e1..e6
TREE   = {1: [0], 2: [0, 1], 3: [3], 4: [3, 4]}             # vertex -> edge indices, path v0->v


def branch_vectors(WF, WC):
    mF = np.ones(5, dtype=complex); mF[FACE_V] = WF
    mC = np.ones(5, dtype=complex); mC[CYC_V] = WC
    return mF, mC


def Z_direct(s, WF, WC, k):
    """Z_k = <M_dF^k s, M_c^k s> with <z,w> = conj(z) w  (S1 sec3)."""
    mF, mC = branch_vectors(WF, WC)
    a = (mF ** k) * s
    b = (mC ** k) * s
    return np.vdot(a, b)          # numpy vdot conjugates the FIRST argument


def pushforward(s):
    q = np.abs(s) ** 2
    p00 = 0.0
    p10 = q[1] + q[2]
    p01 = q[3] + q[4]
    p11 = q[0]
    return p00, p10, p01, p11


def P_eval(p, x, y):
    p00, p10, p01, p11 = p
    return p00 + p10 * x + p01 * y + p11 * x * y


def holonomies(a):
    WF = np.exp(1j * (a[0] + a[1] + a[2]))
    WC = np.exp(1j * (a[3] + a[4] + a[5]))
    return WF, WC


def dress(s, a):
    """t_v = W(tree path v0->v)^{-1} s_v, tree {e1,e2,e4,e5}, root v0."""
    t = np.array(s, dtype=complex)
    for v, path in TREE.items():
        ph = np.exp(1j * sum(a[e] for e in path))
        t[v] = s[v] / ph
    return t


print("=" * 78)
print("M1_01 — THE IDENTIFICATION Z_k = P(u^k, v^k) ON K1")
print("=" * 78)

# ---------------------------------------------------------------- T1 float sweep
rng = np.random.default_rng(20260816)
NDRAW, KMAX = 4000, 60
worst_conj, worst_plain = 0.0, 0.0
arg_conj = arg_plain = None
for _ in range(NDRAW):
    a = rng.uniform(0.0, 2 * np.pi, 6)
    s = rng.normal(size=5) + 1j * rng.normal(size=5)
    s = s / np.linalg.norm(s)
    WF, WC = holonomies(a)
    p = pushforward(s)
    k = int(rng.integers(1, KMAX + 1))
    Z = Z_direct(s, WF, WC, k)
    u_c, v_c = np.conj(WF), WC                 # CONVENTION A (claimed)
    u_p, v_p = WF, WC                          # CONVENTION B (the alternative named in the brief)
    d_c = abs(Z - P_eval(p, u_c ** k, v_c ** k))
    d_p = abs(Z - P_eval(p, u_p ** k, v_p ** k))
    if d_c > worst_conj:
        worst_conj, arg_conj = d_c, (a.copy(), k)
    if d_p > worst_plain:
        worst_plain, arg_plain = d_p, (a.copy(), k)

print("\nT1  float64 sweep: %d draws, random connection, random complex ready state, k<=%d" %
      (NDRAW, KMAX))
print("    CONVENTION A   u = conj(W_F), v = W_C   worst |Z_k - P(u^k,v^k)|  = %.3e" % worst_conj)
print("    CONVENTION B   u = W_F,       v = W_C   worst |Z_k - P(u^k,v^k)|  = %.3e" % worst_plain)
print("    -> CONVENTION A is the one that makes the identity exact.  B is wrong by O(1).")

# also: is |Z_k| the same under B?  (it is not, in general — record it)
rng2 = np.random.default_rng(20260816 + 1)
worst_abs_B = 0.0
for _ in range(2000):
    a = rng2.uniform(0.0, 2 * np.pi, 6)
    s = rng2.normal(size=5) + 1j * rng2.normal(size=5); s /= np.linalg.norm(s)
    WF, WC = holonomies(a); p = pushforward(s); k = int(rng2.integers(1, 60))
    Z = Z_direct(s, WF, WC, k)
    worst_abs_B = max(worst_abs_B, abs(abs(Z) - abs(P_eval(p, WF ** k, WC ** k))))
print("    and |Z_k| vs |P(W_F^k, W_C^k)| under B: worst discrepancy = %.3e  (B fails in modulus too)"
      % worst_abs_B)

# CONVENTION C: the only other consistent choice, (u,v) -> (W_F, conj W_C) = (conj u, conj v).
# P has REAL coefficients, so P(conj a, conj b) = conj P(a,b): this convention returns
# conj(Z_k) and therefore leaves |Z_k|, Omega_N and lambda unchanged.  Checked, not assumed.
rng2b = np.random.default_rng(20260816 + 2)
worst_C = 0.0
for _ in range(2000):
    a = rng2b.uniform(0.0, 2 * np.pi, 6)
    s = rng2b.normal(size=5) + 1j * rng2b.normal(size=5); s /= np.linalg.norm(s)
    WF, WC = holonomies(a); p = pushforward(s); k = int(rng2b.integers(1, 60))
    Z = Z_direct(s, WF, WC, k)
    worst_C = max(worst_C, abs(np.conj(Z) - P_eval(p, WF ** k, np.conj(WC) ** k)))
print("    CONVENTION C   u = W_F, v = conj(W_C):  worst |conj(Z_k) - P(u^k,v^k)| = %.3e"
      % worst_C)
print("    -> C is A composed with complex conjugation; it changes Z_k but not |Z_k|, not")
print("       Omega_N and not lambda.  A is the convention of record because <z,w>=conj(z)w")
print("       conjugates the FIRST slot and M_dF acts there.")

# ---------------------------------------------------------------- T2 exact arithmetic
# Gaussian rationals as (re, im) pairs of Fraction.
def gmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def gadd(a, b):
    return (a[0] + b[0], a[1] + b[1])
def gconj(a):
    return (a[0], -a[1])
def gpow(a, n):
    r = (Fraction(1), Fraction(0))
    for _ in range(n):
        r = gmul(r, a)
    return r

ROOTS4 = {1: (Fraction(1), Fraction(0)), 1j: (Fraction(0), Fraction(1)),
          -1: (Fraction(-1), Fraction(0)), -1j: (Fraction(0), Fraction(-1))}
R4 = list(ROOTS4.values())

# rational squared-moduli, summing to 1 — several ready states including S1's published one
STATES = [
    [Fraction(1, 2), Fraction(0), Fraction(0), Fraction(1, 4), Fraction(1, 4)],   # S1 published p
    [Fraction(1, 5)] * 5,
    [Fraction(2, 7), Fraction(1, 7), Fraction(1, 7), Fraction(2, 7), Fraction(1, 7)],
    [Fraction(1), Fraction(0), Fraction(0), Fraction(0), Fraction(0)],            # root only
    [Fraction(0), Fraction(1, 2), Fraction(1, 2), Fraction(0), Fraction(0)],      # FACE only
]
maxdev = Fraction(0)
ncase = 0
for WFg in R4:
    for WCg in R4:
        for q in STATES:
            p11, p10, p01 = q[0], q[1] + q[2], q[3] + q[4]
            for k in range(1, 13):
                # direct: sum_v conj((M_dF^k s)_v) (M_c^k s)_v = sum_v conj(bF_v^k) bC_v^k |s_v|^2
                Z = (Fraction(0), Fraction(0))
                for v in range(5):
                    bF = WFg if v in FACE_V else (Fraction(1), Fraction(0))
                    bC = WCg if v in CYC_V else (Fraction(1), Fraction(0))
                    term = gmul(gconj(gpow(bF, k)), gpow(bC, k))
                    Z = gadd(Z, (term[0] * q[v], term[1] * q[v]))
                u, v_ = gconj(WFg), WCg
                Pv = gadd(gadd((p10 * gpow(u, k)[0], p10 * gpow(u, k)[1]),
                               (p01 * gpow(v_, k)[0], p01 * gpow(v_, k)[1])),
                          (lambda t: (p11 * t[0], p11 * t[1]))(gmul(gpow(u, k), gpow(v_, k))))
                dev = abs(Z[0] - Pv[0]) + abs(Z[1] - Pv[1])
                maxdev = max(maxdev, dev)
                ncase += 1
print("\nT2  EXACT (Gaussian-rational) check: %d cases  [W_F,W_C in 4th roots of unity;" % ncase)
print("    5 rational ready states incl. S1's published p=(1/2,0,0,1/4,1/4); k=1..12]")
print("    max exact deviation |Z_k - P(u^k,v^k)| = %s   (Fraction, printed exactly)" % maxdev)
print("    -> the identification is an ALGEBRAIC IDENTITY, not a numerical coincidence.")

# ---------------------------------------------------------------- T3 gauge invariance (full action)
rng3 = np.random.default_rng(20260816 + 3)
worst_gauge = 0.0
for _ in range(2000):
    a = rng3.uniform(0, 2 * np.pi, 6)
    s = rng3.normal(size=5) + 1j * rng3.normal(size=5); s /= np.linalg.norm(s)
    th = rng3.uniform(0, 2 * np.pi, 5)
    a2 = a.copy()
    for e, (src, tgt) in enumerate(EDGES):
        a2[e] = a[e] + th[tgt] - th[src]          # S1:63, on the CONNECTION
    s2 = np.exp(1j * th) * s                       # and on the SECTION
    k = int(rng3.integers(1, 40))
    Z1 = Z_direct(s, *holonomies(a), k)
    Z2 = Z_direct(s2, *holonomies(a2), k)
    worst_gauge = max(worst_gauge, abs(Z1 - Z2))
print("\nT3  full gauge action (connection AND section): worst |Z_k - Z_k^gauge| = %.3e  (float64)"
      % worst_gauge)

# ---------------------------------------------------------------- T4 dressing invariance
rng4 = np.random.default_rng(20260816 + 4)
worst_dress = 0.0
for _ in range(2000):
    a = rng4.uniform(0, 2 * np.pi, 6)
    s = rng4.normal(size=5) + 1j * rng4.normal(size=5); s /= np.linalg.norm(s)
    t = dress(s, a)
    k = int(rng4.integers(1, 40))
    worst_dress = max(worst_dress, abs(Z_direct(s, *holonomies(a), k) -
                                       Z_direct(t, *holonomies(a), k)))
print("T4  tree dressing t_v = W(path)^-1 s_v: worst |Z_k[s] - Z_k[t]| = %.3e" % worst_dress)
print("    (expected 0: |t_v| = |s_v| for a U(1) dressing, and Z_k sees only |s_v|^2)")

# ---------------------------------------------------------------- T5 p00 = 0 on K1
rng5 = np.random.default_rng(20260816 + 5)
mx = 0.0
for _ in range(5000):
    s = rng5.normal(size=5) + 1j * rng5.normal(size=5); s /= np.linalg.norm(s)
    p = pushforward(s)
    mx = max(mx, abs(p[0]), abs(sum(p) - 1.0))
print("\nT5  max(|p00|, |sum p - 1|) over 5000 random ready states = %.3e" % mx)
print("    p00 = 0 is forced by K1's incidence (no vertex outside both loops), not chosen.")

# ---------------------------------------------------------------- corollaries used downstream
print("\nCOROLLARIES (checked here, used in M1_02/03/04):")
print("  Z_0 = P(1,1) = sum p = 1;  |Z_k| <= 1 for every k and every connection.")
rng6 = np.random.default_rng(20260816 + 6)
mx1 = 0.0
for _ in range(20000):
    a = rng6.uniform(0, 2 * np.pi, 6)
    s = rng6.normal(size=5) + 1j * rng6.normal(size=5); s /= np.linalg.norm(s)
    k = int(rng6.integers(1, 200))
    mx1 = max(mx1, abs(Z_direct(s, *holonomies(a), k)))
print("  max |Z_k| over 20000 random (connection, state, k) draws = %.15f  (<= 1)" % mx1)
print("\nDONE M1_01")
