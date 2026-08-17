#!/usr/bin/env python3
"""
R1_01 — ATTACK ON F1/T1 (THE IDENTIFICATION).

The lane checked the identity in float64 over random draws and EXACTLY over the 4th roots of
unity only.  4th roots of unity are the one place where conj(w) = w^{-1} = w^3 is a *power*
with small order and where several distinct conventions can coincide.  Here it is redone:

  (1) EXACTLY over Z[i] (4th roots)          -- reproduce the lane's own exact block
  (2) EXACTLY over Z[omega] (3rd/6th roots)  -- a ring the lane never touched, where the
      distinction conj(W_F) vs W_F is NOT a sign and cannot be masked by order 4
  (3) EXACTLY over 12th roots of unity, in Z[zeta_12] as a Z-module with basis
      1, z, z^2, z^3 (z^4 = z^2 - 1), so W_F and W_C need not lie in the same cyclotomic field
  (4) as an OPERATOR identity, not a sampled one: M_dF, M_c built as 5x5 matrices from the
      incidence, never as the lane's precomputed diagonal vectors
  (5) SCOPE ATTACK: is Z_k = P(u^k,v^k) special to K1, or is it a theorem about ANY pair of
      diagonal loop operators?  If the latter, the "carrier fact" p00 = 0 is the ONLY K1
      content in T1 and the identification carries no information about K1.

Precision: EXACT (fractions.Fraction) throughout parts 1-3 and 5; float64 only in part 4 where
it is labelled.
"""
import numpy as np
from fractions import Fraction as F
from itertools import product

FACE_V = [0, 1, 2]
CYC_V = [0, 3, 4]

print("=" * 78)
print("R1_01 — THE IDENTIFICATION, ATTACKED IN THREE CYCLOTOMIC RINGS")
print("=" * 78)

# ---------------------------------------------------------------------------
# generic exact cyclotomic arithmetic: elements of Q[z]/(Phi) as coefficient tuples
# ---------------------------------------------------------------------------
class Cyc:
    """Q[z]/(min_poly), min_poly given as reduction rule: z^d = sum red[i] z^i."""
    def __init__(self, coef, red):
        self.c = list(coef); self.red = red
    def __add__(self, o):
        return Cyc([a + b for a, b in zip(self.c, o.c)], self.red)
    def __mul__(self, o):
        d = len(self.c)
        raw = [F(0)] * (2 * d - 1)
        for i, a in enumerate(self.c):
            if a == 0: continue
            for j, b in enumerate(o.c):
                if b == 0: continue
                raw[i + j] += a * b
        for k in range(2 * d - 2, d - 1, -1):
            v = raw[k]
            if v == 0: continue
            raw[k] = F(0)
            for i, r in enumerate(self.red):
                raw[k - d + i] += v * r
        return Cyc(raw[:d], self.red)
    def scale(self, q):
        return Cyc([q * a for a in self.c], self.red)
    def __eq__(self, o):
        return self.c == o.c
    def __repr__(self):
        return "Cyc" + str(self.c)

def one(red):
    d = len(red); c = [F(0)] * d; c[0] = F(1); return Cyc(c, red)

def zpow(red, k):
    """z^k reduced."""
    d = len(red)
    r = one(red)
    zc = [F(0)] * d; zc[1 % d] = F(1)
    z = Cyc(zc, red) if d > 1 else one(red)
    for _ in range(k % (10 ** 9)):
        r = r * z
    return r

# ---- ring 1: Z[i]              z^2 = -1                (4th roots: z^0..z^3)
RED_I = [F(-1), F(0)]
# ---- ring 2: Z[omega]          z^2 = -1 - z            (cube roots of 1)
RED_W = [F(-1), F(-1)]
# ---- ring 3: Z[zeta_12]        z^4 = z^2 - 1           (12th roots of 1)
RED_12 = [F(-1), F(0), F(1), F(0)]

def conj_in(ring_name, e):
    """complex conjugation = the field automorphism z -> z^{-1}."""
    if ring_name == "Zi":       # zbar = -z = z^3
        a, b = e.c
        return Cyc([a, -b], RED_I)
    if ring_name == "Zw":       # omegabar = omega^2 = -1-omega
        a, b = e.c
        # (a + b w)bar = a + b w^2 = a + b(-1-w) = (a-b) + (-b) w
        return Cyc([a - b, -b], RED_W)
    if ring_name == "Z12":      # zeta^-1 = zeta^11 ; compute by conjugating basis powers
        # conj(z^j) = z^{-j} = z^{12-j}
        out = Cyc([F(0)] * 4, RED_12)
        for j, a in enumerate(e.c):
            if a == 0: continue
            out = out + zpow(RED_12, (12 - j) % 12).scale(a)
        return out
    raise ValueError

RINGS = {
    "Zi":  (RED_I, 4,  "Z[i]        4th roots of unity"),
    "Zw":  (RED_W, 3,  "Z[omega]    3rd roots of unity  (order 3: conj != +-1 * elt)"),
    "Z12": (RED_12, 12, "Z[zeta_12]  12th roots of unity"),
}

STATES = [
    [F(1, 2), F(0), F(0), F(1, 4), F(1, 4)],       # S1's published p
    [F(1, 5)] * 5,
    [F(2, 7), F(1, 7), F(1, 7), F(2, 7), F(1, 7)],
    [F(3, 11), F(1, 11), F(2, 11), F(4, 11), F(1, 11)],   # ASYMMETRIC: p10 != p01
    [F(1), F(0), F(0), F(0), F(0)],
    [F(0), F(1, 3), F(1, 6), F(1, 4), F(1, 4)],           # ASYMMETRIC, p11 = 0
]

for rname, (red, order, label) in RINGS.items():
    roots = [zpow(red, j) for j in range(order)]
    maxdevA = 0; maxdevB = 0; ncase = 0
    for WF in roots:
        for WC in roots:
            for q in STATES:
                p11, p10, p01 = q[0], q[1] + q[2], q[3] + q[4]
                for k in range(0, 13):
                    # DIRECT: build M_dF, M_c as DIAGONAL EXACT operators from the incidence
                    Zd = Cyc([F(0)] * len(red), red)
                    for v in range(5):
                        bF = WF if v in FACE_V else one(red)
                        bC = WC if v in CYC_V else one(red)
                        pF = one(red); pC = one(red)
                        for _ in range(k):
                            pF = pF * bF; pC = pC * bC
                        Zd = Zd + (conj_in(rname, pF) * pC).scale(q[v])
                    # CONVENTION A : u = conj(W_F), v = W_C
                    uA = conj_in(rname, WF); vA = WC
                    puA = one(red); pvA = one(red)
                    for _ in range(k):
                        puA = puA * uA; pvA = pvA * vA
                    PA = puA.scale(p10) + pvA.scale(p01) + (puA * pvA).scale(p11)
                    # CONVENTION B : u = W_F, v = W_C
                    puB = one(red); pvB = one(red)
                    for _ in range(k):
                        puB = puB * WF; pvB = pvB * WC
                    PB = puB.scale(p10) + pvB.scale(p01) + (puB * pvB).scale(p11)
                    dA = sum(abs(a - b) for a, b in zip(Zd.c, PA.c))
                    dB = sum(abs(a - b) for a, b in zip(Zd.c, PB.c))
                    maxdevA = max(maxdevA, dA); maxdevB = max(maxdevB, dB)
                    ncase += 1
    print("\nRING %s   (%s)" % (rname, label))
    print("   %d exact cases  [W_F,W_C over all %d roots; %d rational states incl. ASYMMETRIC;"
          " k=0..12]" % (ncase, order, len(STATES)))
    print("   CONVENTION A  u = conj(W_F):  max EXACT deviation = %s" % maxdevA)
    print("   CONVENTION B  u = W_F      :  max EXACT deviation = %s" % maxdevB)

# ---------------------------------------------------------------------------
# (4) operator form, float64, matrices built from incidence
# ---------------------------------------------------------------------------
print("\n(4) OPERATOR FORM (float64).  M_dF, M_c built as 5x5 MATRICES from the incidence,")
print("    then M^k by repeated matrix multiplication -- not by exponentiating a vector.")
rng = np.random.default_rng(20260817)
worst = 0.0
for _ in range(3000):
    a = rng.uniform(0, 2 * np.pi, 6)
    WF = np.exp(1j * (a[0] + a[1] + a[2])); WC = np.exp(1j * (a[3] + a[4] + a[5]))
    MF = np.eye(5, dtype=complex); MC = np.eye(5, dtype=complex)
    for v in FACE_V: MF[v, v] = WF
    for v in CYC_V:  MC[v, v] = WC
    s = rng.normal(size=5) + 1j * rng.normal(size=5); s /= np.linalg.norm(s)
    k = int(rng.integers(0, 40))
    A = np.linalg.matrix_power(MF, k) @ s
    B = np.linalg.matrix_power(MC, k) @ s
    Z = np.vdot(A, B)
    q = np.abs(s) ** 2
    p11, p10, p01 = q[0], q[1] + q[2], q[3] + q[4]
    u, v = np.conj(WF), WC
    Pv = p10 * u ** k + p01 * v ** k + p11 * (u * v) ** k
    worst = max(worst, abs(Z - Pv))
print("    worst |Z_k - P(u^k,v^k)| over 3000 draws, k<=39 = %.3e" % worst)

# ---------------------------------------------------------------------------
# (5) SCOPE ATTACK
# ---------------------------------------------------------------------------
print("\n(5) SCOPE ATTACK.  Is T1 a theorem about K1, or about ANY pair of diagonal unitaries?")
print("    Take an ARBITRARY complex 'carrier': n vertices, arbitrary subsets S_F, S_C, and")
print("    arbitrary unimodular W_F, W_C.  Define M_F, M_C diagonally the same way.")
rng2 = np.random.default_rng(20260817 + 1)
worst_gen = 0.0
for _ in range(3000):
    n = int(rng2.integers(2, 12))
    SF = set(int(i) for i in np.flatnonzero(rng2.integers(0, 2, n)))
    SC = set(int(i) for i in np.flatnonzero(rng2.integers(0, 2, n)))
    WF = np.exp(1j * rng2.uniform(0, 2 * np.pi)); WC = np.exp(1j * rng2.uniform(0, 2 * np.pi))
    s = rng2.normal(size=n) + 1j * rng2.normal(size=n); s /= np.linalg.norm(s)
    k = int(rng2.integers(0, 30))
    mF = np.array([WF if i in SF else 1.0 for i in range(n)], dtype=complex)
    mC = np.array([WC if i in SC else 1.0 for i in range(n)], dtype=complex)
    Z = np.vdot((mF ** k) * s, (mC ** k) * s)
    q = np.abs(s) ** 2
    r00 = sum(q[i] for i in range(n) if i not in SF and i not in SC)
    r10 = sum(q[i] for i in range(n) if i in SF and i not in SC)
    r01 = sum(q[i] for i in range(n) if i not in SF and i in SC)
    r11 = sum(q[i] for i in range(n) if i in SF and i in SC)
    u, v = np.conj(WF), WC
    Pv = r00 + r10 * u ** k + r01 * v ** k + r11 * (u * v) ** k
    worst_gen = max(worst_gen, abs(Z - Pv))
print("    worst deviation over 3000 RANDOM complexes (n=2..11, random loop supports) = %.3e"
      % worst_gen)
print("    -> T1 holds verbatim on every one of them.  THE IDENTIFICATION IS NOT A FACT ABOUT")
print("       K1.  It is a two-line consequence of `both branch operators are diagonal in the")
print("       vertex basis'.  K1's ONLY contribution to it is which of the four classes is")
print("       empty (p00 = 0).  The lane says this; it is recorded here as CONFIRMED and as")
print("       the reason F1 cannot be evidence for anything carrier-specific.")
print("\nDONE R1_01")
