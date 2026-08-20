"""Exact fixed-point-counting trace engine extracted from v_core.py (verifier code)."""
import itertools
import numpy as np
from fractions import Fraction

# ---------------------------------------------------------------- exact trace engine
def sector_trace(car, sector, Mspec):
    """Tr(P_s M) exactly.  P_s = Q0^{a0} Q1^{a1} D0^{b0} D1^{b1}, Q^1 = A-avg, Q^0 = I - A-avg,
    D^1 = diagB, D^0 = 1 - diagB.  Mspec: ("perm", pi) or ("diag", d) or ("permdiag", pi, d) meaning diag then perm? --
    we only ever need pure perm or pure diag."""
    a0, a1, b0, b1 = sector
    n, N = car.n, car.N
    d = (car.dB0() if b0 else 1 - car.dB0()) * (car.dB1() if b1 else 1 - car.dB1())
    idx = np.arange(N)
    total = Fraction(0)
    terms0 = [(1, True)] if a0 == 1 else [(1, False), (-1, True)]   # (sign, use A0 average?)
    terms1 = [(1, True)] if a1 == 1 else [(1, False), (-1, True)]
    for s0, useA0 in terms0:
        for s1, useA1 in terms1:
            ks0 = range(n) if useA0 else [None]
            ks1 = range(n) if useA1 else [None]
            w = Fraction(s0 * s1)
            if useA0: w /= n
            if useA1: w /= n
            acc = 0
            for k0 in ks0:
                p0 = car.A0(k0) if k0 is not None else idx
                for k1 in ks1:
                    p01 = p0[car.A1(k1)] if k1 is not None else p0
                    if Mspec[0] == "perm":
                        pim = Mspec[1]
                        comp = p01[pim]           # A0 A1 M as index map
                        fixed = comp == idx
                        acc += int(d[pim[fixed]].sum()) if fixed.any() else 0
                    else:                          # diag
                        m = Mspec[1]
                        fixed = p01 == idx
                        acc += int((d[fixed] * m[fixed]).sum()) if fixed.any() else 0
            total += w * acc
    return total

def eig_traces(car, Mspec):
    """{eig k (H eigenvalue -k): exact Tr(P_E M)} via the 16 sectors."""
    out = {}
    for s in itertools.product((0, 1), repeat=4):
        k = sum(s)
        out[k] = out.get(k, Fraction(0)) + sector_trace(car, s, Mspec)
    return out

