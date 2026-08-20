"""V3 -- THREE TARGETED ATTACKS, EACH RUN RATHER THAN ARGUED.

A.  D-16 DENOMINATOR AUDIT.  The finding states "chi_SPREAD (=q1 value) is the D-16 denominator,
    printed on every row" and "Because the couplings commute, SPREAD equals the one-record value
    exactly".  In S9 the spread layout is  assign = [i % 3]  with a THREE-qubit bath, so at k = 4
    and k = 6 the spread layout already puts TWO records on site 0.  Recomputed here.

B.  IS THE CLAUSE-(v) ZERO A FUNCTION OF k AT ALL?  Two carriers with the SAME n and the SAME k
    but different code distance are put side by side.  If reach1 differs between them at equal k,
    then reach1/cf1_max is a function of the CARRIER, not of the record count, and "cf1_max = 0
    at every k up to 126" is a statement about the families chosen, not about k.
    (The lane's positive control for this zero varies the REGION RADIUS, not k.)

C.  IS reach_r CONSTANT IN k, OR CONSTANT IN k AND n AND EVERYTHING ELSE?  reach_r is computed
    across the whole of family A and family B; if it depends only on r, the k-sweep could not
    have registered a k-threshold in it under any circumstances.
"""
import sys, time
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD")
import numpy as np
from math import comb
import carriers as C
import battery as BAT

OUT = []
def P(s=""):
    print(s); OUT.append(s)

Z = np.array([[1, 0], [0, -1]], dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
TIMES = np.linspace(1.0, 13.0, 25)

def vN(r):
    e = np.linalg.eigvalsh(r); e = e[e > 1e-13]
    return float(-(e * np.log2(e)).sum())

def site_rho(e, lam, c, beta, t):
    p = np.exp(-beta * np.array([e, -e])); p = p / p.sum()
    r0 = np.diag(p).astype(complex)
    w, V = np.linalg.eigh(e * Z + lam * c * X)
    U = (V * np.exp(-1j * w * t)) @ V.conj().T
    return U @ r0 @ U.conj().T

def chi_site(q, e, lam, beta, t):
    cs = {}
    for j in range(q + 1): cs[q - 2 * j] = cs.get(q - 2 * j, 0) + comb(q, j)
    states = {c: site_rho(e, lam, c, beta, t) for c in cs}
    halves = []
    for s in (1, -1):
        cnt = {}
        for j in range(q):
            c = s + (q - 1 - 2 * j); cnt[c] = cnt.get(c, 0) + comb(q - 1, j)
        h = float(2 ** (q - 1))
        for c in cnt:
            if c not in states: states[c] = site_rho(e, lam, c, beta, t)
        halves.append(sum((m / h) * states[c] for c, m in cnt.items()))
    av = 0.5 * (halves[0] + halves[1])
    return max(vN(av) - 0.5 * (vN(halves[0]) + vN(halves[1])), 0.0)

def chi_avg(q, e=1.0, lam=0.8, beta=2.0):
    return float(np.mean([chi_site(q, e, lam, beta, t) for t in TIMES]))

t0 = time.time()
P("=" * 120)
P("V3  DENOMINATOR AUDIT, k-FREEDOM OF THE CLAUSE-(v) ZERO, AND THE reach_r FORMULA")
P("=" * 120)

# ------------------------------------------------------------------ A
P()
P("A.  D-16 DENOMINATOR AUDIT OF S9 (the lane's only non-analytic dynamical arm)")
P("    S9 spread layout is assign=[i%3] on a 3-qubit bath.  q0 = records landing on site 0.")
P("    %-12s %-4s %-16s %-5s %-14s %-14s %-16s" %
  ("carrier", "k", "assign", "q0", "chi_avg(q0)", "chi_avg(1)", "S9 printed SPREAD"))
P("    " + "-" * 86)
S9_SPREAD = {2: 0.52152730, 4: 0.13640869, 6: 0.13640869}
c1 = chi_avg(1)
for k in (2, 4, 6):
    assign = [i % 3 for i in range(k)]
    q0 = sum(1 for a in assign if a == assign[0])
    P("    %-12s %-4d %-16s %-5d %-14.8f %-14.8f %-16.8f"
      % ("[[%d,%d,2]]" % (k + 2, k), k, str(assign), q0, chi_avg(q0), c1, S9_SPREAD[k]))
P("    -> the S9 'SPREAD' column is the q0-record value, and q0 = 2 at k = 4 and k = 6.")
P("       It is NOT the one-record value at those k, so the printed commut/spread ratios")
P("       0.2616 / 0.3311 / 0.1782 have a DENOMINATOR THAT CHANGES VENUE between rows.")
P()
P("    Same suppression re-normalised against a TRUE one-record spread (D-16 as stated):")
P("    %-12s %-4s %-14s %-16s %-16s" % ("carrier", "k", "chi_COMMUTING", "lane commut/spread",
                                        "commut/chi_avg(1)"))
P("    " + "-" * 66)
LANE = {2: (0.13640869, 0.26155618), 4: (0.04516257, 0.33108283), 6: (0.02430704, 0.17819278)}
for k in (2, 4, 6):
    cc, lr = LANE[k]
    P("    %-12s %-4d %-14.8f %-16.6f %-16.6f" % ("[[%d,%d,2]]" % (k + 2, k), k, cc, lr, cc / c1))
P("    -> corrected, the suppression is MONOTONE in k (0.2616, 0.0866, 0.0466); the")
P("       'the ratio does not trend' reading came from the moving denominator.")

# ------------------------------------------------------------------ B
P()
P("B.  IS THE CLAUSE-(v) ZERO A FUNCTION OF k?  SAME n, SAME k, DIFFERENT CARRIER.")
def bare_padded(m, pad):
    """m blocks of [[4,2,2]] plus `pad` BARE qubits: n = 4m+pad, k = 2m+pad, distance 1."""
    n = 4 * m + pad
    S = []
    for b in range(m):
        x = [0] * (2 * n); z = [0] * (2 * n)
        for j in range(4 * b, 4 * b + 4):
            x[j] = 1; z[n + j] = 1
        S.append(x); S.append(z)
    return dict(name="P", label="[[4,2,2]]^%d+%dbare" % (m, pad), n=n, dim=2 ** n, stabs=S)

P("    %-22s %-5s %-5s %-9s %-9s %-9s %-9s" %
  ("carrier", "n", "k", "reach1", "reach2", "reach3", "prot_r"))
P("    " + "-" * 70)
def prot(car):
    for r in range(1, car["n"] + 1):
        if max(BAT.local_logical_dim(car, reg) for reg in BAT.regions(car, r)) > 0: return r
    return None
def c513_padded(pad):
    """one [[5,1,3]] block plus `pad` BARE qubits: n = 5+pad, k = 1+pad, distance 1."""
    n = 5 + pad
    pat = ["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"]
    S = []
    for p in pat:
        v = [0] * (2 * n)
        for j, ch in enumerate(p):
            if ch in "XY": v[j] = 1
            if ch in "ZY": v[n + j] = 1
        S.append(v)
    return dict(name="P", label="[[5,1,3]]+%dbare" % pad, n=n, dim=2 ** n, stabs=S)

pairs = [(C.family_A(6), bare_padded(1, 2)),          # n=6  k=4  both
         (C.family_A(10), bare_padded(1, 6)),         # n=10 k=8  both
         (C.family_B(2), c513_padded(3))]             # n=8  k=4  both
for a, b in pairs:
    for car in (a, b):
        n = car["n"]; k = n - len(car["stabs"])
        rs = [max(BAT.local_logical_dim(car, reg) for reg in BAT.regions(car, r)) for r in (1, 2, 3)]
        P("    %-22s %-5d %-5d %-9d %-9d %-9d %-9s" % (car["label"], n, k, rs[0], rs[1], rs[2], prot(car)))
    P("    " + "." * 70)
P("    -> at IDENTICAL (n, k) the clause-(v) radius-1 quantity is 0 on one carrier and > 0 on")
P("       the other.  reach1/cf1_max is a function of the CARRIER'S DISTANCE, not of k.")

# ------------------------------------------------------------------ C
P()
P("C.  reach_r ACROSS FAMILY A AND B: does it depend on k at all, or only on r?")
P("    %-16s %-5s %-5s %-8s %-8s %-8s %-8s" % ("carrier", "n", "k", "reach1", "reach2", "reach3", "reach4"))
P("    " + "-" * 62)
vals = {}
for n in (4, 6, 8, 10, 12, 14, 16, 20):
    car = C.family_A(n); k = n - 2
    rs = [max(BAT.local_logical_dim(car, reg) for reg in BAT.regions(car, r)) for r in (1, 2, 3, 4)]
    vals.setdefault("A", []).append(tuple(rs))
    P("    %-16s %-5d %-5d %-8d %-8d %-8d %-8d" % (car["label"], n, k, *rs))
for m in (1, 2, 3, 4, 6):
    car = C.family_B(m); k = 2 * m
    rs = [max(BAT.local_logical_dim(car, reg) for reg in BAT.regions(car, r)) for r in (1, 2, 3, 4)]
    vals.setdefault("B", []).append(tuple(rs))
    P("    %-16s %-5d %-5d %-8d %-8d %-8d %-8d" % (car["label"], car["n"], k, *rs))
for fam in ("A", "B"):
    same = len(set(vals[fam][1:])) == 1
    P("    family %s: reach vector identical for every k >= 4 ?  %s   (%s)"
      % (fam, same, sorted(set(vals[fam]))))
P("    -> for family A, reach_r = 2(r-1) with NO k and NO n in it: the quantity the lane swept")
P("       to k = 126 is a closed form in the region radius alone.")

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD/VERIFY/v3_spread_and_kfree.txt",
     "w").write("\n".join(OUT) + "\n")
P("total %.1fs" % (time.time() - t0))
