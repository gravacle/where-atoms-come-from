"""STEP 6 -- push the BASIS-FREE object (all logical classes at minimum weight) past n = 8.

The exact class set has 4^k - 1 members, so the distance matrix is 4^k x 4^k: exact to n = 8
(4095 classes), 65535 classes at n = 10 which is a 34 GB matrix.  Beyond n = 8 a uniform random
SUBSAMPLE of classes is used.  SC-13 validates the subsample against the exact value at n = 8
before any subsampled number is reported.

CONTROL IN THE SAME TABLE: the FREE carrier (k unentangled qubits, H = 0) treated identically --
its 4^k - 1 non-trivial Pauli classes, same subsample size, same statistics.
"""
import sys
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_B_METRIC")
from lib_scaleb import *
import numpy as np, itertools

OUT = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); OUT.append(s)

rng = np.random.default_rng(777)

def fast_dc(D):
    D2 = D ** 2
    return -0.5 * (D2 - D2.mean(1, keepdims=True) - D2.mean(0, keepdims=True) + D2.mean())

def dimsD(D):
    B = fast_dc(D)
    w = np.linalg.eigvalsh((B + B.T) / 2)[::-1]
    pos = w[w > 1e-10]
    if not len(pos): return 0, 0.0, 0.0
    c = np.cumsum(pos) / pos.sum()
    neg = float(-w[w < -1e-10].sum())
    return int(np.searchsorted(c, 0.90) + 1), float(pos.sum() ** 2 / (pos ** 2).sum()), \
           neg / (neg + pos.sum())

def coeff_to_pauli(bits, gens, n):
    v = [0] * (2 * n)
    for c, g in zip(bits, gens):
        if c: v = pauli_mul(v, g, n)
    return v

def canonical_min_weight(v, n):
    best = None
    for g in stab_group(n):
        w = pauli_mul(v, g, n)
        key = (len(support(w, n)), tuple(w))
        if best is None or key < best[0]: best = (key, w)
    return best[1]

P("=" * 110)
P("LANE_SCALE_B_METRIC  STEP 6 -- BASIS-FREE LOGICAL-CLASS GEOMETRY, PUSHED AS FAR AS IT GOES")
P("=" * 110)
P("")
P("SC-13  does a random subsample of classes reproduce the exact value?  n = 8, 4095 classes.")
P("  %-10s %-10s %-10s %-10s" % ("sample", "d90 CODE", "dPR CODE", "cdim CODE"))
P("  " + "-" * 44)
n = 8
stab, pairs = carrier(n); vs, lab = record_vectors(pairs, n); k = len(pairs)
allcls = []
for bits in itertools.product((0, 1), repeat=2 * k):
    if not any(bits): continue
    allcls.append(canonical_min_weight(coeff_to_pauli(bits, vs, n), n))
for sz in [500, 1000, 2000, len(allcls)]:
    idx = rng.choice(len(allcls), size=min(sz, len(allcls)), replace=False)
    cc = [allcls[i] for i in idx]
    D = 1.0 - M_support(cc, n); np.fill_diagonal(D, 0.0)
    d90, dpr, negf = dimsD(D)
    P("  %-10d %-10d %-10.2f %-10.2f" % (len(cc), d90, dpr, corr_dim(D)))

P("")
P("6B.  CLASS-LEVEL GEOMETRY vs n.  1500 random classes per carrier (exact set when smaller).")
P("     CODE = [[n,n-2,2]] logical classes at minimum weight.  FREE = k free qubits, all")
P("     non-trivial single-carrier Pauli classes.  Both use the SUPPORT (Jaccard) relation")
P("     and the LETTER (Hamming) relation, d = 1 - J and d = Hamming/n.")
P("")
P("  %-4s %-4s %-9s | %-26s | %-26s" %
  ("n", "k", "#classes", "CODE  d90J dPRJ cdimJ d90H", "FREE  d90J dPRJ cdimJ d90H"))
P("  " + "-" * 78)
SAMP = 1500
for n in [4, 6, 8, 10, 12, 14]:
    stab, pairs = carrier(n); vs, lab = record_vectors(pairs, n); k = len(pairs)
    total = 4 ** k - 1
    cells = []
    for tag in ("CODE", "FREE"):
        cls = []
        if total <= SAMP:
            for bits in itertools.product((0, 1), repeat=2 * k):
                if not any(bits): continue
                if tag == "CODE":
                    cls.append(canonical_min_weight(coeff_to_pauli(bits, vs, n), n))
                else:
                    cls.append([int(b) for b in bits])
        else:
            seen = set()
            while len(cls) < SAMP:
                bits = [int(b) for b in rng.integers(0, 2, 2 * k)]
                if not any(bits): continue
                t = tuple(bits)
                if t in seen: continue
                seen.add(t)
                cls.append(canonical_min_weight(coeff_to_pauli(bits, vs, n), n) if tag == "CODE" else bits)
        nn = n if tag == "CODE" else k
        DJ = 1.0 - M_support(cls, nn); np.fill_diagonal(DJ, 0.0)
        DH = M_hamming(cls, nn); np.fill_diagonal(DH, 0.0)
        d90J, dprJ, _ = dimsD(DJ)
        d90H, _, _ = dimsD(DH)
        cells.append("%5d %5.2f %5.2f %5d" % (d90J, dprJ, corr_dim(DJ), d90H))
    P("  %-4d %-4d %-9d | %-26s | %-26s" % (n, k, min(total, SAMP), cells[0], cells[1]))

P("")
P("  READ: if a geometry emerged from the mutual relations of many records, d90 would approach")
P("  a CONSTANT as n grows.  Fill this from the numbers above, never in advance.")

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_B_METRIC/06_classes_large.txt", "w").write("\n".join(OUT) + "\n")
