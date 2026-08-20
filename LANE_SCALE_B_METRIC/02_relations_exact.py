"""STEP 2 -- EXACT relation matrices between records, their spectra, and whether an
intrinsic dimension STABILISES as n grows.

Four relation definitions, all computed EXACTLY (F2 / set algebra, no sampling):
  (a) SYMPLECTIC   S_ij = <c_i,c_j>            1 iff R_i and R_j anticommute   GAUGE-INVARIANT
  (c) CODE OVERLAP C_ij = |Tr(P_g R_iR_j)|/Tr(P_g)                              GAUGE-INVARIANT
  (d) SUPPORT      J_ij = Jaccard of the physical supports                      gauge-DEPENDENT
  (e) LETTERS      Hm_ij = Pauli-letter Hamming distance / n                    gauge-DEPENDENT

STATED MONOTONE TRANSFORMS TO A DISTANCE (fixed in advance, never fitted):
  (a) d = 1 - S            (b) d = sqrt(1 - C)       (d) d = 1 - J        (e) d = Hm

CONTROLS CARRIED IN THE SAME TABLE (D-15):
  RAND   a random symmetric matrix of the SAME SIZE and SAME off-diagonal density and value set
  FREE   a STRUCTURELESS carrier: k independent unentangled qubits, H = 0, records = the
         2k single-qubit Paulis X_q, Z_q.  Same number of records, no code, no entanglement.
D-17 (vary the venue's own scale): n is swept 4..20, i.e. 2k = 4..36 records.
"""
import sys
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_B_METRIC")
from lib_scaleb import *
import numpy as np

OUT = []
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); OUT.append(s)

NS = [4, 6, 8, 10, 12, 14, 16, 18, 20]
rng = np.random.default_rng(20260819)

# ------------------------------------------------------------------ the free-qubit control
def free_carrier(k):
    """k independent qubits, H = 0.  records = X_q and Z_q.  Returned as (x|z) over k qubits."""
    vs, lab = [], []
    for q in range(k):
        v = [0] * (2 * k); v[q] = 1; vs.append(v); lab.append("x%d" % (q + 1))
    for q in range(k):
        v = [0] * (2 * k); v[k + q] = 1; vs.append(v); lab.append("z%d" % (q + 1))
    return vs, lab

def free_codespace_overlap(vs, k):
    """H = 0 -> the whole space is the ground space, P_g = I.  Tr(R_iR_j)/dim = delta_ij."""
    m = len(vs)
    return np.eye(m)

# ------------------------------------------------------------------ the analysis of one matrix
def analyse(D, name):
    sym = is_symmetric(D)
    m = D.shape[0]
    nviol, worst = triangle_violations(D) if m <= 40 else (-1, -1.0)
    B = double_centre(D)
    st = dim_stats(B)
    cd = corr_dim(D)
    return dict(name=name, m=m, sym=sym, nviol=nviol, worst=worst,
                d90=st["d_frac"], dpr=st["d_pr"], negfrac=st.get("neg_frac", 0.0),
                corrdim=cd, evals=st["evals"], spread=float(D.max() - D[~np.eye(m, dtype=bool)].min()))

def spec_str(w, k=8):
    return "[" + " ".join("%+.3f" % x for x in w[:k]) + (" ...]" if len(w) > k else "]")

P("=" * 118)
P("LANE_SCALE_B_METRIC  STEP 2 -- EXACT RELATION MATRICES AND THEIR GEOMETRY")
P("=" * 118)

# =================================================================== 2A raw relation matrices
P("")
P("2A.  THE RELATION MATRICES THEMSELVES (n = 6 shown in full; all n computed)")
n = 6
stab, pairs = carrier(n); vs, lab = record_vectors(pairs, n)
mats = dict(a=M_symplectic(vs, n), c=M_codespace_overlap(vs, n),
            d=M_support(vs, n), e=M_hamming(vs, n))
for key, title in [("a", "(a) SYMPLECTIC  <c_i,c_j>"), ("c", "(c) CODE-SPACE OVERLAP"),
                   ("d", "(d) SUPPORT JACCARD"), ("e", "(e) LETTER HAMMING / n")]:
    P("")
    P("   " + title + "     rows/cols = " + " ".join(lab))
    M = mats[key]
    for i in range(len(lab)):
        P("   %-4s " % lab[i] + " ".join("%6.3f" % M[i, j] for j in range(len(lab))))

# =================================================================== 2B spectra
P("")
P("=" * 118)
P("2B.  EIGENVALUE SPECTRUM OF EACH RELATION MATRIX AT EACH n   (control columns in the same table)")
P("     RAND = random symmetric matrix, same size, same off-diagonal density and value set")
P("     FREE = k unentangled qubits, H = 0, records = the 2k single-qubit Paulis")
P("")

rows = []
for n in NS:
    stab, pairs = carrier(n); vs, lab = record_vectors(pairs, n); k = len(pairs); m = 2 * k
    fvs, flab = free_carrier(k)
    Mn = dict(a=M_symplectic(vs, n), c=M_codespace_overlap(vs, n),
              d=M_support(vs, n), e=M_hamming(vs, n))
    Mf = dict(a=M_symplectic(fvs, k), c=free_codespace_overlap(fvs, k),
              d=M_support(fvs, k), e=M_hamming(fvs, k))
    dens = sparsity(Mn["a"])
    Mr = random_control(m, dens, rng, vals=[1.0])
    rows.append((n, k, m, Mn, Mf, Mr, dens, vs, lab))

for key, title in [("a", "(a) SYMPLECTIC"), ("c", "(c) CODE OVERLAP"),
                   ("d", "(d) SUPPORT JACCARD"), ("e", "(e) LETTER HAMMING")]:
    P("")
    P("  %s" % title)
    P("  %-4s %-4s %-9s %-46s %-46s" % ("n", "2k", "offdiagdens", "spectrum, CODE carrier (top 8)",
                                        "spectrum, FREE control (top 8)"))
    P("  " + "-" * 112)
    for (n, k, m, Mn, Mf, Mr, dens, vs, lab) in rows:
        wn = np.linalg.eigvalsh(Mn[key])[::-1]
        wf = np.linalg.eigvalsh(Mf[key])[::-1]
        P("  %-4d %-4d %-9.4f %-46s %-46s" % (n, m, sparsity(Mn[key]), spec_str(wn), spec_str(wf)))
    P("  RAND control spectrum at each n (same size and density as (a)):")
    for (n, k, m, Mn, Mf, Mr, dens, vs, lab) in rows:
        wr = np.linalg.eigvalsh(Mr)[::-1]
        P("  %-4d %-4d %-9.4f %-46s" % (n, m, sparsity(Mr), spec_str(wr)))
    break_after = True

# =================================================================== 2C geometry test
P("")
P("=" * 118)
P("2C.  IS IT A GEOMETRY?  metric test + intrinsic dimension, CODE vs FREE vs RAND, same table")
P("     d90  = # MDS eigenvalues carrying 90% of the positive variance")
P("     dPR  = participation ratio of the positive MDS eigenvalues (a smooth dimension)")
P("     neg  = fraction of |MDS spectrum| that is NEGATIVE (0 => exactly Euclidean-embeddable)")
P("     tri  = # of ordered triples violating the triangle inequality (0 => a metric)")
P("     cdim = correlation dimension from the pair-distance distribution (nan if degenerate)")
P("")

def to_dist(key, M):
    m = M.shape[0]
    if key == "a":  D = 1.0 - M
    elif key == "c": D = np.sqrt(np.maximum(1.0 - np.abs(M), 0.0))
    elif key == "d": D = 1.0 - M
    else:            D = M.copy()
    np.fill_diagonal(D, 0.0)
    return D

for key, title in [("a", "(a) SYMPLECTIC     d = 1 - S"),
                   ("c", "(c) CODE OVERLAP   d = sqrt(1-|C|)"),
                   ("d", "(d) SUPPORT        d = 1 - J"),
                   ("e", "(e) LETTERS        d = Hamming/n")]:
    P("")
    P("  %s" % title)
    P("  %-4s %-4s | %-24s | %-24s | %-24s" % ("n", "2k", "CODE  d90 dPR   neg  tri cdim",
                                               "FREE  d90 dPR   neg  tri cdim",
                                               "RAND  d90 dPR   neg  tri cdim"))
    P("  " + "-" * 112)
    for (n, k, m, Mn, Mf, Mr, dens, vs, lab) in rows:
        cells = []
        for tag, M in [("CODE", Mn[key]), ("FREE", Mf[key]),
                       ("RAND", random_control(m, sparsity(Mn[key]), rng,
                                               vals=sorted(set(np.round(Mn[key][~np.eye(m, dtype=bool)], 6))
                                                           - {0.0}) or [1.0]))]:
            D = to_dist(key, M)
            r = analyse(D, tag)
            cd = "  nan" if np.isnan(r["corrdim"]) else "%5.2f" % r["corrdim"]
            cells.append("%3d %5.2f %5.3f %4d %s" % (r["d90"], r["dpr"], r["negfrac"], r["nviol"], cd))
        P("  %-4d %-4d | %-24s | %-24s | %-24s" % (n, m, cells[0], cells[1], cells[2]))

# =================================================================== 2D gauge dependence
P("")
P("=" * 118)
P("2D.  GAUGE TEST -- a record is only defined MODULO THE STABILISER GROUP.")
P("     Multiply every logical by a randomly drawn stabiliser element and recompute.")
P("     A relation that CHANGES is bookkeeping about the representative, not a relation")
P("     between records.  10 random gauge draws per n; max change and the d90 range reported.")
P("")
P("  %-4s %-4s | %-26s | %-26s | %-26s | %-26s" %
  ("n", "2k", "(a) max-change  d90 range", "(c) max-change  d90 range",
   "(d) max-change  d90 range", "(e) max-change  d90 range"))
P("  " + "-" * 140)
grng = np.random.default_rng(4242)
for (n, k, m, Mn, Mf, Mr, dens, vs, lab) in rows:
    S = stab_group(n)
    chg = {key: 0.0 for key in "acde"}
    d90s = {key: set() for key in "acde"}
    for key in "acde":
        d90s[key].add(analyse(to_dist(key, Mn[key]), "base")["d90"])
    for _ in range(10):
        gvs = [pauli_mul(v, S[grng.integers(0, 4)], n) for v in vs]
        Mg = dict(a=M_symplectic(gvs, n), c=M_codespace_overlap(gvs, n),
                  d=M_support(gvs, n), e=M_hamming(gvs, n))
        for key in "acde":
            chg[key] = max(chg[key], float(np.max(np.abs(Mg[key] - Mn[key]))))
            d90s[key].add(analyse(to_dist(key, Mg[key]), "g")["d90"])
    cells = []
    for key in "acde":
        cells.append("%10.4f      %d-%d" % (chg[key], min(d90s[key]), max(d90s[key])))
    P("  %-4d %-4d | %-26s | %-26s | %-26s | %-26s" % (n, m, cells[0], cells[1], cells[2], cells[3]))

# =================================================================== 2E drift
P("")
P("=" * 118)
P("2E.  DOES THE INTRINSIC DIMENSION STABILISE OR DRIFT WITH n?")
P("     A geometry means d90 approaches a CONSTANT.  Bookkeeping means d90 tracks 2k.")
P("")
P("  %-6s %-6s | %-30s | %-30s" % ("n", "2k", "CODE  d90 (a/c/d/e)", "d90 / (2k-1)  (a/c/d/e)"))
P("  " + "-" * 84)
drift = {key: [] for key in "acde"}
for (n, k, m, Mn, Mf, Mr, dens, vs, lab) in rows:
    ds = {}
    for key in "acde":
        ds[key] = analyse(to_dist(key, Mn[key]), key)["d90"]
        drift[key].append((m, ds[key]))
    P("  %-6d %-6d | %-30s | %-30s" %
      (n, m, "  ".join("%2d" % ds[key] for key in "acde"),
       "  ".join("%.3f" % (ds[key] / (m - 1)) for key in "acde")))
P("")
for key, title in [("a", "(a) SYMPLECTIC"), ("c", "(c) CODE OVERLAP"),
                   ("d", "(d) SUPPORT"), ("e", "(e) LETTERS")]:
    xs = np.array([x for x, _ in drift[key]], float); ys = np.array([y for _, y in drift[key]], float)
    sl = np.polyfit(xs, ys, 1)[0]
    P("  %-18s d90 vs 2k slope = %+.4f   %s" %
      (title, sl, "SATURATES (slope ~ 0)" if abs(sl) < 0.05 else "GROWS WITH n -- no fixed dimension"))

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_B_METRIC/02_relations_exact.txt", "w").write("\n".join(OUT) + "\n")
