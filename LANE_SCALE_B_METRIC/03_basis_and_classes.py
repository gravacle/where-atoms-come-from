"""STEP 3 -- IS THE LOW-DIMENSIONAL SUPPORT GEOMETRY OF STEP 2 REAL, OR AN ARTIFACT OF WHICH
2k RECORDS WERE CHOSEN AS A BASIS?

Step 2 found: under the support/letter relations, the CODE carrier's intrinsic dimension
SATURATES near 3-4 while the FREE control's grows like 2k.  Two ways that could be an artifact:

  ARTIFACT 1 -- THE BASIS.  The records of the carrier are the WHOLE logical Pauli group;
    the 2k that symplectic_logicals happens to return are one symplectic basis out of
    |Sp(2k,F2)| of them.  CONTROL: apply a random element of Sp(2k,F2) (a product of random
    symplectic transvections) to the basis.  The result is still 2k independent records with a
    non-degenerate pairing -- SC-1 is re-checked on every draw.  If the dimension moves, the
    step-2 number described the Gram-Schmidt order, not the carrier.

  ARTIFACT 2 -- THE REPRESENTATIVE.  A record is a class modulo the stabiliser group.
    CONTROL: drop the basis entirely and use ALL 4^k - 1 non-trivial logical classes, each
    with its canonical MINIMUM-WEIGHT representative.  That object has no basis in it at all.

CONTROL IN THE SAME TABLE: FREE = k unentangled qubits, H = 0, same construction throughout.
"""
import sys
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_B_METRIC")
from lib_scaleb import *
import numpy as np
import itertools

OUT = []
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); OUT.append(s)

rng = np.random.default_rng(31337)

def fast_double_centre(D):
    D2 = D ** 2
    rm = D2.mean(axis=1, keepdims=True); cm = D2.mean(axis=0, keepdims=True); gm = D2.mean()
    return -0.5 * (D2 - rm - cm + gm)

def dims(D):
    B = fast_double_centre(D)
    w = np.linalg.eigvalsh((B + B.T) / 2)[::-1]
    pos = w[w > 1e-10]
    if len(pos) == 0: return 0, 0.0, 0.0
    c = np.cumsum(pos) / pos.sum()
    d90 = int(np.searchsorted(c, 0.90) + 1)
    dpr = float(pos.sum() ** 2 / (pos ** 2).sum())
    neg = float(-w[w < -1e-10].sum()); negf = neg / (neg + pos.sum())
    return d90, dpr, negf

def sp2k(a, b, k):
    return sum(a[i] * b[k + i] + a[k + i] * b[i] for i in range(2 * k) if i < k) % 2

def sp_coeff(a, b, k):
    """symplectic form on F_2^{2k} coefficient space, convention (x_1..x_k | z_1..z_k)"""
    return sum(a[i] * b[k + i] + a[k + i] * b[i] for i in range(k)) % 2

def random_symplectic_basis(k, rng, ntrans=40):
    """A random element of Sp(2k,F2) as a product of transvections, applied to the standard
       symplectic basis of the COEFFICIENT space.  Returns 2k coefficient vectors."""
    basis = [[1 if j == i else 0 for j in range(2 * k)] for i in range(2 * k)]
    for _ in range(ntrans):
        v = [int(x) for x in rng.integers(0, 2, 2 * k)]
        if not any(v): continue
        basis = [[(b[i] + sp_coeff(b, v, k) * v[i]) % 2 for i in range(2 * k)] for b in basis]
    return basis

def coeff_to_pauli(coefs, gens, n):
    """map an F2 coefficient vector over the 2k record generators to an (x|z) Pauli vector"""
    v = [0] * (2 * n)
    for c, g in zip(coefs, gens):
        if c:
            v = pauli_mul(v, g, n)
    return v

def canonical_min_weight(v, n):
    """the minimum-weight representative of the logical class of v -- basis-free and
       representative-free (ties broken lexicographically)."""
    best = None
    for g in stab_group(n):
        w = pauli_mul(v, g, n)
        key = (len(support(w, n)), tuple(w))
        if best is None or key < best[0]:
            best = (key, w)
    return best[1]

P("=" * 118)
P("LANE_SCALE_B_METRIC  STEP 3 -- BASIS ARTIFACT AND REPRESENTATIVE ARTIFACT CONTROLS")
P("=" * 118)

# ============================================================ 3A random symplectic basis
P("")
P("3A.  RANDOM SYMPLECTIC BASIS OF RECORDS.  Same carrier, same 2k-dimensional record algebra,")
P("     a DIFFERENT choice of which 2k records form the basis.  SC-1 re-checked every draw.")
P("     8 random draws per n.  Reported: d90 for the support relation (d) and the letter")
P("     relation (e), as GS = the Gram-Schmidt basis of step 2, RS = random symplectic basis.")
P("")
P("  %-4s %-4s | %-28s | %-28s | %-16s" %
  ("n", "2k", "(d) SUPPORT  d90  GS / RS(min-max)", "(e) LETTERS  d90  GS / RS(min-max)",
   "FREE control d90"))
P("  " + "-" * 96)

def f2_nondeg(vs, n):
    m = len(vs)
    A = [[sp_form(vs[i], vs[j], n) for j in range(m)] for i in range(m)]
    for c in range(m):
        p = next((r for r in range(c, m) if A[r][c]), None)
        if p is None: return False
        A[c], A[p] = A[p], A[c]
        for r in range(m):
            if r != c and A[r][c]:
                A[r] = [(x + y) % 2 for x, y in zip(A[r], A[c])]
    return True

NS = [4, 6, 8, 10, 12, 14, 16, 18, 20]
sc1_all = True
summary_rs = {}
for n in NS:
    stab, pairs = carrier(n); vs, lab = record_vectors(pairs, n); k = len(pairs); m = 2 * k
    gs_d = dims(1.0 - M_support(vs, n))[0]
    gs_e = dims(M_hamming(vs, n))[0]
    rs_d, rs_e = [], []
    for _ in range(8):
        cb = random_symplectic_basis(k, rng)
        nvs = [coeff_to_pauli(c, vs, n) for c in cb]
        if not f2_nondeg(nvs, n):
            sc1_all = False
            continue
        rs_d.append(dims(1.0 - M_support(nvs, n))[0])
        rs_e.append(dims(M_hamming(nvs, n))[0])
    # FREE control gets the identical treatment: its records span the same algebra
    fvs = []
    for q in range(k):
        v = [0] * (2 * k); v[q] = 1; fvs.append(v)
    for q in range(k):
        v = [0] * (2 * k); v[k + q] = 1; fvs.append(v)
    free_d = dims(1.0 - M_support(fvs, k))[0]
    summary_rs[n] = (gs_d, min(rs_d), max(rs_d), gs_e, min(rs_e), max(rs_e), free_d)
    P("  %-4d %-4d | %10d  /  %2d-%-2d %8s | %10d  /  %2d-%-2d %8s | %-16d" %
      (n, m, gs_d, min(rs_d), max(rs_d), "", gs_e, min(rs_e), max(rs_e), "", free_d))
P("")
P("  SC-1 (non-degenerate pairing) held on EVERY random symplectic draw: %s" % sc1_all)

# ============================================================ 3B all logical classes
P("")
P("=" * 118)
P("3B.  NO BASIS AT ALL -- every one of the 4^k - 1 non-trivial logical classes, each taken at")
P("     its canonical MINIMUM-WEIGHT representative.  This object contains no basis choice and")
P("     no representative choice.  CONTROL: the FREE carrier's 4^k - 1 non-trivial Pauli classes.")
P("")
P("  %-4s %-4s %-8s | %-30s | %-30s" % ("n", "k", "#classes",
                                        "CODE   d90  dPR   neg   cdim", "FREE   d90  dPR   neg   cdim"))
P("  " + "-" * 92)
CLASS_N = [4, 6, 8]
for n in CLASS_N:
    stab, pairs = carrier(n); vs, lab = record_vectors(pairs, n); k = len(pairs)
    cls, fcls = [], []
    for bits in itertools.product((0, 1), repeat=2 * k):
        if not any(bits): continue
        cls.append(canonical_min_weight(coeff_to_pauli(bits, vs, n), n))
        fv = [0] * (2 * k)
        for i, b in enumerate(bits):
            if b: fv[i] = 1
        fcls.append(fv)
    cells = []
    for tag, cc, nn in [("CODE", cls, n), ("FREE", fcls, k)]:
        D = 1.0 - M_support(cc, nn)
        np.fill_diagonal(D, 0.0)
        d90, dpr, negf = dims(D)
        cd = corr_dim(D)
        cells.append("%4d %5.2f %6.3f %6s" % (d90, dpr, negf,
                                              " nan" if np.isnan(cd) else "%.2f" % cd))
    P("  %-4d %-4d %-8d | %-30s | %-30s" % (n, k, len(cls), cells[0], cells[1]))

# ============================================================ 3C is it a line?
P("")
P("=" * 118)
P("3C.  WHAT SHAPE IS THE GRAM-SCHMIDT SUPPORT GEOMETRY?  MDS coordinate 1 against record")
P("     weight, n = 20.  If the 'geometry' is just the weight staircase that Gram-Schmidt")
P("     produced, coordinate 1 is a monotone function of weight and nothing else is there.")
P("")
n = 20
stab, pairs = carrier(n); vs, lab = record_vectors(pairs, n)
D = 1.0 - M_support(vs, n); np.fill_diagonal(D, 0.0)
B = fast_double_centre(D)
w, V = np.linalg.eigh((B + B.T) / 2)
idx = np.argsort(w)[::-1]
c1 = V[:, idx[0]] * np.sqrt(max(w[idx[0]], 0)); c2 = V[:, idx[1]] * np.sqrt(max(w[idx[1]], 0))
wts = [len(support(v, n)) for v in vs]
P("  %-5s %-8s %-10s %-10s" % ("rec", "weight", "MDS c1", "MDS c2"))
order = np.argsort(c1)
for i in order:
    P("  %-5s %-8d %-10.4f %-10.4f" % (lab[i], wts[i], c1[i], c2[i]))
cc = np.corrcoef(np.array(wts, float), c1)[0, 1]
P("")
P("  |corr(weight, MDS coordinate 1)| = %.4f" % abs(cc))
P("  variance carried by coordinate 1 = %.4f of the positive MDS variance"
  % (w[idx[0]] / w[w > 1e-10].sum()))

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_B_METRIC/03_basis_and_classes.txt", "w").write("\n".join(OUT) + "\n")
