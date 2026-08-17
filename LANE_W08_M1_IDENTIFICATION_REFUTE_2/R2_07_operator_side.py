#!/usr/bin/env python3
"""
R2_07 — THE THREE-CONNECTION COMPARISON REBUILT FROM THE 5x5 BRANCH OPERATORS, NOT FROM THE
POLYNOMIAL.  A fresh implementation: no import from the target lane, no use of P at all in
the Z_k evaluation.  This is the check that the whole of M1_03 is a statement about
Z_k = <M_dF^k s, M_c^k s> and not about a substitute object.

ISOLATION LEDGER: HELD FIXED -- carrier K1, the ready state s (S3/S4's own
p = (0.4,0.15,0.15,0.15,0.15) as amplitudes, so pi = (0,0.3,0.3,0.4)), the k-grid, float64.
THE ONE THING THAT MOVES: the implementation route (5x5 diagonal operators + vdot  vs  the
3-term polynomial).  If they disagree the identification is not what carries M1_03.
Precision: float64.  Also an EXACT integer check for the order-4 connection.
"""
import numpy as np
from fractions import Fraction

FACE_V = [0, 1, 2]
CYC_V = [0, 3, 4]
q = np.array([0.4, 0.15, 0.15, 0.15, 0.15])
s = np.sqrt(q).astype(complex)                       # amplitudes; any phases give the same Z_k


def Zk_operator(WF, WC, K):
    mF = np.ones(5, dtype=complex); mF[FACE_V] = WF
    mC = np.ones(5, dtype=complex); mC[CYC_V] = WC
    out = np.empty(K, dtype=complex)
    aF = np.ones(5, dtype=complex)
    aC = np.ones(5, dtype=complex)
    for k in range(K):
        aF = aF * mF
        aC = aC * mC
        out[k] = np.vdot(aF * s, aC * s)             # <M_dF^k s, M_c^k s>
    return out


def Zk_poly(u, v, K):
    k = np.arange(1, K + 1)
    return 0.3 * u ** k + 0.3 * v ** k + 0.4 * (u * v) ** k


print("=" * 78)
print("R2_07 — M1_03 REBUILT FROM THE OPERATORS.  fresh implementation, no target imports.")
print("=" * 78)
K = 100000
CASES = [
    ("A  Diophantine  W_F=e^{i*2^(1/3)*2pi}, W_C=e^{-i*4^(1/3)*2pi}",
     np.exp(2j * np.pi * (2.0 ** (1 / 3))), np.exp(2j * np.pi * (4.0 ** (1 / 3)))),
    ("B  S1 published  W_F=-1, W_C=-i", -1 + 0j, -1j),
    ("C  S3/S4 headline  W_F=e^{2i}, W_C=e^{1.1i}", np.exp(2j), np.exp(1.1j)),
]
print("\n   ready state s = sqrt(p), p = (0.4,0.15,0.15,0.15,0.15)  ->  pi = (0,0.3,0.3,0.4)")
print("   K = %d, float64.\n" % K)
print("   case                                              worst |Z_op - P(u^k,v^k)|   "
      "(1/K)sum log|Z_op|")
for (name, WF, WC) in CASES:
    Zo = Zk_operator(WF, WC, K)
    u, v = np.conj(WF), WC
    Zp = Zk_poly(u, v, K)
    d = float(np.max(np.abs(Zo - Zp)))
    lam = float(np.mean(np.log(np.abs(Zo))))
    print("   %-50s %.3e         %.9f" % (name, d, lam))

print("\n   EXACT check for case B (Gaussian integers, no float):")
WF = (Fraction(-1), Fraction(0)); WC = (Fraction(0), Fraction(-1))
qf = [Fraction(2, 5), Fraction(3, 20), Fraction(3, 20), Fraction(3, 20), Fraction(3, 20)]


def gmul(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def gcj(a): return (a[0], -a[1])
def gpow(a, n):
    r = (Fraction(1), Fraction(0))
    for _ in range(n): r = gmul(r, a)
    return r


worst = Fraction(0)
mods = []
for k in range(1, 13):
    Z = (Fraction(0), Fraction(0))
    for vtx in range(5):
        bF = WF if vtx in FACE_V else (Fraction(1), Fraction(0))
        bC = WC if vtx in CYC_V else (Fraction(1), Fraction(0))
        t = gmul(gcj(gpow(bF, k)), gpow(bC, k))
        Z = (Z[0] + t[0] * qf[vtx], Z[1] + t[1] * qf[vtx])
    u = gcj(WF); v = WC
    Pv = (Fraction(3, 10) * gpow(u, k)[0] + Fraction(3, 10) * gpow(v, k)[0]
          + Fraction(2, 5) * gmul(gpow(u, k), gpow(v, k))[0],
          Fraction(3, 10) * gpow(u, k)[1] + Fraction(3, 10) * gpow(v, k)[1]
          + Fraction(2, 5) * gmul(gpow(u, k), gpow(v, k))[1])
    worst = max(worst, abs(Z[0] - Pv[0]) + abs(Z[1] - Pv[1]))
    mods.append(Z[0] * Z[0] + Z[1] * Z[1])
print("     max exact |Z_k - P(u^k,v^k)|, k=1..12 = %s" % worst)
print("     |Z_k|^2 for k=1..4 (exact) = %s" % [str(m) for m in mods[:4]])
print("     product over one period = %s  ->  lambda = (1/4) log(that) = %.15f"
      % (mods[0] * mods[1] * mods[2] * mods[3],
         0.125 * float(np.log(float(mods[0] * mods[1] * mods[2] * mods[3])))))
print("     -(1/2) log 5 = %.15f" % (-0.5 * np.log(5)))
print("\nDONE R2_07")
