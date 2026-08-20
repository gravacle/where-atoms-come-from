"""Fresh definitions extracted from v_core.py (verifier code, independent of t24_lib)."""
import numpy as np
from fractions import Fraction

# ---------------------------------------------------------------- groups (fresh)
def group_D4():
    # element = (a, b) meaning r^a s^b ; s r s = r^-1 ; encode a*2 + b  (DIFFERENT from lane: lane used a + 4b)
    n = 8
    enc = lambda a, b: (a % 4) * 2 + (b % 2)
    dec = lambda g: (g // 2, g % 2)
    MUL = np.zeros((n, n), dtype=np.int64)
    for g1 in range(n):
        a1, b1 = dec(g1)
        for g2 in range(n):
            a2, b2 = dec(g2)
            MUL[g1, g2] = enc(a1 + (a2 if b1 == 0 else -a2), b1 + b2)
    INV = np.array([[h for h in range(n) if MUL[g, h] == enc(0, 0)][0] for g in range(n)], dtype=np.int64)
    E = enc(0, 0)
    # associativity check (full)
    ok_assoc = all(MUL[MUL[a, b], c] == MUL[a, MUL[b, c]] for a in range(n) for b in range(n) for c in range(n))
    assert ok_assoc
    # conjugacy classes
    classes = []
    seen = set()
    for g in range(n):
        if g in seen: continue
        cl = sorted({int(MUL[k, MUL[g, INV[k]]]) for k in range(n)})
        seen |= set(cl); classes.append(cl)
    # characters (5 irreps): 1-dim eps^a del^b, and 2-dim (2 on e, -2 on r^2, else 0)
    chars = []
    for eps in (1, -1):
        for dl in (1, -1):
            chars.append(np.array([eps ** dec(g)[0] * dl ** dec(g)[1] for g in range(n)], dtype=np.int64))
    two = np.zeros(n, dtype=np.int64); two[enc(0, 0)] = 2; two[enc(2, 0)] = -2
    chars.append(two)
    r2 = enc(2, 0); s = enc(0, 1)
    return dict(n=n, MUL=MUL, INV=INV, E=E, classes=classes, chars=chars, r2=r2, s=s, enc=enc, dec=dec)

def group_Z2():
    MUL = np.array([[0, 1], [1, 0]], dtype=np.int64)
    return dict(n=2, MUL=MUL, INV=np.array([0, 1], dtype=np.int64), E=0,
                classes=[[0], [1]], chars=[np.array([1, 1], np.int64), np.array([1, -1], np.int64)],
                r2=None, s=None)

# ---------------------------------------------------------------- carrier (fresh layout)
class Car:
    """config index = ((u1*n + u0)*n + h1)*n + h0  -- REVERSED vs the lane."""
    def __init__(self, G):
        self.G = G; n = G["n"]; self.n = n; self.N = n ** 4
        idx = np.arange(self.N)
        self.h0 = idx % n; r = idx // n
        self.h1 = r % n; r = r // n
        self.u0 = r % n; self.u1 = r // n
        self.M = G["MUL"]; self.I = G["INV"]
    def pack(self, h0, h1, u0, u1):
        n = self.n
        return ((u1 * n + u0) * n + h1) * n + h0
    def A0(self, k):
        M, I = self.M, self.I
        return self.pack(M[k, self.h0], M[self.h1, I[k]], M[M[k, self.u0], I[k]], self.u1)
    def A1(self, k):
        M, I = self.M, self.I
        return self.pack(M[self.h0, I[k]], M[k, self.h1], self.u0, M[M[k, self.u1], I[k]])
    def dB0(self):
        M, I = self.M, self.I
        return (self.u0 == M[M[self.h0, self.u1], I[self.h0]]).astype(np.int64)
    def dB1(self):
        M, I = self.M, self.I
        return (self.u1 == M[M[self.h1, self.u0], I[self.h1]]).astype(np.int64)
    def comp(self, e):
        return {"h0": self.h0, "h1": self.h1, "u0": self.u0, "u1": self.u1}[e]
    def with_comp(self, e, new):
        c = {"h0": self.h0, "h1": self.h1, "u0": self.u0, "u1": self.u1}
        c = dict(c); c[e] = new
        return self.pack(c["h0"], c["h1"], c["u0"], c["u1"])

