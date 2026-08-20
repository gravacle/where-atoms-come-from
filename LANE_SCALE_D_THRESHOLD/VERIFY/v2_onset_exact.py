"""V2 -- THE ONSET SWEEP, DONE EXACTLY, AND AN ORDINARY EXPLANATION FOR IT.

The lane's only arm that separates "number of records" from "density of records" is S6/S6b:
random stabiliser carriers at FIXED n, k swept, k*(n,r) = smallest k where more than half the
ensemble has reach_r > 0.  Its D-17 kill ("k* grows with n at every radius") rests on that
table.  Two things are checked here.

(1) reach_r > 0 IS EXACTLY "the code has distance <= r".
    reach_r > 0  <=>  some region of r qubits supports a Pauli in N(S) that is not in S
                 <=>  min weight of N(S)\S  <=  r  <=>  d <= r.
    So the whole onset arm is a distance computation, and it can be done EXACTLY by
    enumerating the 3^w C(n,w) Paulis of weight w <= r -- no region sampling at all.
    That removes the nsamp = 120 estimator whose coverage falls 15x across the lane's n range.

(2) FIRST-MOMENT COUNTING for a random stabiliser code predicts the onset with no free
    parameters.  A uniformly random Pauli lies in N(S) with probability 2^(k-n); there are
    3^r C(n,r) Paulis of weight exactly r, so the expected number of weight-<=r elements of
    N(S) crosses 1 at
         k*_pred(n,r) = n - log2( sum_{w<=r} 3^w C(n,w) ).
    If that reproduces the measured k*, the onset is textbook random-code weight enumeration
    (Gilbert-Varshamov), not a record-count phenomenon -- and it explains WHY k* must move
    with n: the predicted k* is n minus a term that grows only like r*log2(n).

The random-code generator is the lane's own (transvection transport), copied verbatim so this
tests the MEASUREMENT, not the ensemble.
"""
import sys, time
from math import comb, log2
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD")
import numpy as np
import carriers as C
from f2 import rank

OUT = []
def P(s=""):
    print(s); OUT.append(s)

# ---------------------------------------------------------------- lane's generator, verbatim
def random_code(n, k, rng, ntrans=None):
    need = n - k
    S = []
    for i in range(need):
        v = [0] * (2 * n); v[n + i] = 1; S.append(v)
    if ntrans is None: ntrans = 3 * n
    for _ in range(ntrans):
        v = [int(x) for x in rng.integers(0, 2, size=2 * n)]
        if not any(v): continue
        S = [[(w[j] + C.sp(w, v, n) * v[j]) % 2 for j in range(2 * n)] for w in S]
    if any(C.sp(S[i], S[j], n) for i in range(need) for j in range(i + 1, need)): return None
    if rank(S, 2 * n) != need: return None
    return S

# ---------------------------------------------------------------- exact distance <= r test
def bitpack(S, n):
    """each stabiliser -> (x_mask, z_mask) as ints"""
    return [(sum(1 << q for q in range(n) if s[q]),
             sum(1 << q for q in range(n) if s[n + q])) for s in S]

def make_reducer(Sb):
    """row-reduced basis of S as 2n-bit ints, for the 'is P in S' test"""
    basis = []
    for (x, z) in Sb:
        v = x | (z << 64)
        for b in basis:
            if (v ^ b) < v: v ^= b
        if v: basis.append(v); basis.sort(reverse=True)
    return basis

def in_S(x, z, basis):
    v = x | (z << 64)
    for b in basis:
        if (v ^ b) < v: v ^= b
    return v == 0

def syn_tables(Sb, n):
    """syndrome (as an int bitmask over the stabilisers) of X_q, Z_q, Y_q for every qubit"""
    m = len(Sb)
    sx = [0] * n; sz = [0] * n
    for q in range(n):
        bx = bz = 0
        for i, (x, z) in enumerate(Sb):
            if (z >> q) & 1: bx |= 1 << i        # X_q anticommutes with s iff s has Z at q
            if (x >> q) & 1: bz |= 1 << i
        sx[q] = bx; sz[q] = bz
    return sx, sz

def distance_le(S, n, r):
    """True iff min weight of N(S)\\S is <= r.  Exact enumeration of weight <= r Paulis."""
    Sb = bitpack(S, n); basis = make_reducer(Sb); sx, sz = syn_tables(Sb, n)
    typ = ((1, 0), (0, 1), (1, 1))            # X, Z, Y on a qubit
    def synt(q, t): return (sx[q] if t[0] else 0) ^ (sz[q] if t[1] else 0)
    from itertools import combinations, product
    for w in range(1, r + 1):
        for qs in combinations(range(n), w):
            for ts in product(typ, repeat=w):
                s = 0
                for q, t in zip(qs, ts): s ^= synt(q, t)
                if s: continue
                x = z = 0
                for q, t in zip(qs, ts):
                    if t[0]: x |= 1 << q
                    if t[1]: z |= 1 << q
                if not in_S(x, z, basis): return True
    return False

t0 = time.time()
P("=" * 120)
P("V2  THE ONSET SWEEP DONE EXACTLY (no region sampling), AND THE COUNTING LAW")
P("=" * 120)

# ---------------------------------------------------------------- self-check on known codes
P()
P("SELF-CHECK: exact distance test on codes whose distance is KNOWN.  A wrong answer here")
P("            invalidates everything below, so it is printed before any sweep.")
P("  %-16s %-6s %-8s %-8s %-8s %-10s" % ("carrier", "n", "d<=1", "d<=2", "d<=3", "verdict"))
P("  " + "-" * 62)
ok = True
for car, dtrue in ((C.family_A(6), 2), (C.family_B(2), 2), (C.family_C(1), 3), (C.family_C(2), 3)):
    n = car["n"]; S = car["stabs"]
    v = [distance_le(S, n, r) for r in (1, 2, 3)]
    good = (v == [r >= dtrue for r in (1, 2, 3)])
    ok = ok and good
    P("  %-16s %-6d %-8s %-8s %-8s %-10s" % (car["label"], n, v[0], v[1], v[2],
                                             "OK d=%d" % dtrue if good else "WRONG"))
# a distance-1 control: one bare qubit tensored on, so d = 1 must be detected
n = 5; S = [[1,1,1,1,0] + [0]*5, [0]*5 + [1,1,1,1,0]]
v1 = distance_le(S, n, 1)
P("  %-16s %-6d %-8s %-8s %-8s %-10s" % ("[[4,2,2]]+bare", 5, v1, True, True,
                                         "OK d=1" if v1 else "WRONG"))
ok = ok and v1
P("  self-check %s" % ("PASSED" if ok else "FAILED -- NO CONCLUSION"))
if not ok:
    open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD/VERIFY/v2_onset_exact.txt",
         "w").write("\n".join(OUT) + "\n")
    sys.exit(0)

# ---------------------------------------------------------------- the sweep
P()
P("EXACT ONSET k*(n,r) = smallest k at which > half of the ensemble has d <= r.")
P("  ensemble = 32 random stabiliser carriers per (n,k), the lane's own generator, seed 2026.")
P("  LANE  = the number s6b reported (nsamp=120 region sampling).")
P("  PRED  = n - log2( sum_{w<=r} 3^w C(n,w) ), the first-moment counting law, NO free params.")
P()
P("  %-5s %-4s %-9s %-9s %-9s %-9s" % ("n", "r", "k*_EXACT", "k*_LANE", "k*_PRED", "exact-pred"))
P("  " + "-" * 52)
LANE = {(10,1):5,(10,2):1,(10,3):1,(12,1):6,(12,2):3,(12,3):1,(16,1):10,(16,2):5,(16,3):4,
        (20,1):14,(20,2):9,(20,3):8,(24,1):17,(24,2):13,(24,3):12,(28,1):21,(28,2):17,(28,3):15}
NENS = 32
rows = []
for n in (10, 12, 16, 20, 24):
    rng = np.random.default_rng(2026)
    codes = {}
    for k in range(1, n):
        cs = []
        for e in range(NENS):
            S = random_code(n, k, rng)
            if S is not None: cs.append(S)
        codes[k] = cs
    for r in (1, 2, 3):
        kstar = None
        for k in range(1, n):
            cs = codes[k]
            if not cs: continue
            frac = sum(1 for S in cs if distance_le(S, n, r)) / len(cs)
            if frac > 0.5:
                kstar = k; break
        pred = n - log2(sum(3 ** w * comb(n, w) for w in range(1, r + 1)))
        lane = LANE.get((n, r))
        P("  %-5d %-4d %-9s %-9s %-9.2f %-9s" %
          (n, r, kstar, lane, pred, ("%+.2f" % (kstar - pred)) if kstar else "-"))
        rows.append((n, r, kstar, lane, pred))
    P("    (t=%.0fs)" % (time.time() - t0))

P()
P("READ, from the numbers above and nowhere else:")
ex1 = [(n, r, ks, ln, pr) for (n, r, ks, ln, pr) in rows if r in (1, 2)]
d1 = [abs(ks - pr) for (n, r, ks, ln, pr) in ex1 if ks]
d3 = [abs(ks - pr) for (n, r, ks, ln, pr) in rows if r == 3 and ks]
P("  * |k*_exact - k*_pred| at r=1,2 : mean %.2f  max %.2f" % (np.mean(d1), np.max(d1)))
P("  * |k*_exact - k*_pred| at r=3   : mean %.2f  max %.2f" % (np.mean(d3), np.max(d3)))
agree = sum(1 for (n, r, ks, ln, pr) in rows if ln is not None and ks == ln)
P("  * lane k* reproduced exactly by the exact computation: %d of %d rows" % (agree, len(rows)))

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD/VERIFY/v2_onset_exact.txt",
     "w").write("\n".join(OUT) + "\n")
P("total %.1fs" % (time.time() - t0))
