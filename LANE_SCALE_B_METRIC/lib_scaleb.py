"""LANE_SCALE_B_METRIC -- shared construction for the [[n,n-2,2]] scaling family.

Nothing here nominates a logical operator.  Every logical comes from
record_model.symplectic_logicals; every clause is checked with record_model's own
clause_iii / clause_iv.
"""
import sys, itertools
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
import numpy as np
from record_model import (RecordModel, Environment, symplectic_logicals, xz_to_matrix,
                          eigenspaces, clause_iii, clause_iv)

I2 = np.eye(2, dtype=complex)
Xm = np.array([[0, 1], [1, 0]], dtype=complex)
Zm = np.array([[1, 0], [0, -1]], dtype=complex)
Ym = np.array([[0, -1j], [1j, 0]], dtype=complex)

# ------------------------------------------------------------------ F2 / Pauli symbolics
def sp_form(a, b, n):
    """symplectic form over F2: 1 iff the two Paulis ANTICOMMUTE"""
    return sum(a[i] * b[n + i] + a[n + i] * b[i] for i in range(n)) % 2

def pauli_mul(a, b, n):
    """(x|z) of the product, ignoring phase"""
    return [(a[i] + b[i]) % 2 for i in range(2 * n)]

def support(v, n):
    return frozenset(q for q in range(n) if v[q] or v[n + q])

def letters(v, n):
    """the Pauli letter at each physical qubit"""
    return ''.join('IXZY'[v[q] + 2 * v[n + q]] for q in range(n))

def stab_group(n):
    """the four elements of the stabiliser group of [[n,n-2,2]], as (x|z) vectors"""
    s1 = [1] * n + [0] * n
    s2 = [0] * n + [1] * n
    return [[0] * (2 * n), s1, s2, pauli_mul(s1, s2, n)]

# ------------------------------------------------------------------ the carrier
def carrier(n):
    """H = -(X^n + Z^n) for the [[n,n-2,2]] code, plus the logical pairs."""
    stab = [[1] * n + [0] * n, [0] * n + [1] * n]
    pairs = symplectic_logicals(stab, n)
    return stab, pairs

def record_vectors(pairs, n):
    """the 2k records: every X_i then every Z_i, as (x|z) vectors, with labels"""
    vs, lab = [], []
    for i, (x, z) in enumerate(pairs):
        vs.append(list(x)); lab.append("X%d" % (i + 1))
    for i, (x, z) in enumerate(pairs):
        vs.append(list(z)); lab.append("Z%d" % (i + 1))
    return vs, lab

def hamiltonian(n):
    Xn = xz_to_matrix([1] * n + [0] * n, n)
    Zn = xz_to_matrix([0] * n + [1] * n, n)
    return -(Xn + Zn)

# ------------------------------------------------------------------ exact relation matrices
def M_symplectic(vs, n):
    m = len(vs)
    return np.array([[sp_form(vs[i], vs[j], n) for j in range(m)] for i in range(m)], dtype=float)

def M_codespace_overlap(vs, n):
    """Tr(P_g R_i R_j)/Tr(P_g), EXACTLY.

       P_g = (1/4) sum_{g in S} g, and Tr(g P) = 0 for Paulis unless gP ~ I, i.e. P in S.
       So the ratio is +-1 when R_iR_j lies in the stabiliser group and 0 otherwise.  The sign
       is recovered from the actual matrices at small n and checked against this."""
    S = [tuple(g) for g in stab_group(n)]
    m = len(vs)
    out = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            p = tuple(pauli_mul(vs[i], vs[j], n))
            out[i, j] = 1.0 if p in S else 0.0
    return out

def M_support(vs, n):
    """Jaccard overlap of the physical supports -- GAUGE-DEPENDENT by construction."""
    m = len(vs); sup = [support(v, n) for v in vs]
    out = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            u = len(sup[i] | sup[j])
            out[i, j] = (len(sup[i] & sup[j]) / u) if u else 0.0
    return out

def M_hamming(vs, n):
    """Pauli-letter Hamming distance / n -- a genuine metric on strings, gauge-DEPENDENT."""
    m = len(vs); L = [letters(v, n) for v in vs]
    out = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            out[i, j] = sum(1 for a, b in zip(L[i], L[j]) if a != b) / n
    return out

# ------------------------------------------------------------------ geometry statistics
def is_symmetric(M, tol=1e-9):
    return float(np.max(np.abs(M - M.T)))

def triangle_violations(D, tol=1e-9):
    """exact count over all ordered triples of d(i,k) > d(i,j)+d(j,k) + tol"""
    m = D.shape[0]; bad = 0; worst = 0.0
    for i in range(m):
        for j in range(m):
            for k in range(m):
                v = D[i, k] - D[i, j] - D[j, k]
                if v > tol:
                    bad += 1; worst = max(worst, v)
    return bad, worst

def double_centre(D):
    """B = -1/2 J D^2 J  (classical MDS)"""
    m = D.shape[0]
    J = np.eye(m) - np.ones((m, m)) / m
    return -0.5 * J @ (D ** 2) @ J

def dim_stats(B, frac=0.90):
    """intrinsic-dimension read-outs from the MDS Gram matrix"""
    w = np.linalg.eigvalsh((B + B.T) / 2)[::-1]
    pos = w[w > 1e-10]
    tot = pos.sum() if len(pos) else 0.0
    if tot <= 0:
        return dict(d_frac=0, d_pr=0.0, neg_mass=float(-w[w < -1e-10].sum()), evals=w)
    c = np.cumsum(pos) / tot
    d_frac = int(np.searchsorted(c, frac) + 1)
    d_pr = float(pos.sum() ** 2 / (pos ** 2).sum())          # participation ratio
    neg = float(-w[w < -1e-10].sum())
    return dict(d_frac=d_frac, d_pr=d_pr, neg_mass=neg, evals=w,
                neg_frac=float(neg / (neg + tot)))

def corr_dim(D):
    """correlation dimension from the pair-distance distribution: slope of
       log C(r) vs log r over the middle of the range.  Returns nan if degenerate."""
    m = D.shape[0]
    d = np.array([D[i, j] for i in range(m) for j in range(i + 1, m)])
    d = d[d > 1e-12]
    if len(d) < 6 or d.max() / d.min() < 1.2:
        return float('nan')
    rs = np.exp(np.linspace(np.log(np.percentile(d, 10)), np.log(np.percentile(d, 90)), 12))
    C = np.array([(d <= r).mean() for r in rs])
    ok = (C > 0) & (C < 1)
    if ok.sum() < 3:
        return float('nan')
    return float(np.polyfit(np.log(rs[ok]), np.log(C[ok]), 1)[0])

def sparsity(M):
    m = M.shape[0]
    off = M[~np.eye(m, dtype=bool)]
    return float((np.abs(off) > 1e-9).mean())

def random_control(m, dens, rng, symmetric=True, vals=None):
    """a random matrix of the SAME SIZE and SAME OFF-DIAGONAL DENSITY"""
    A = np.zeros((m, m))
    for i in range(m):
        for j in range(i + 1, m):
            if rng.random() < dens:
                v = rng.choice(vals) if vals is not None else rng.random()
                A[i, j] = v
                A[j, i] = v if symmetric else (rng.choice(vals) if vals is not None else rng.random())
    return A
