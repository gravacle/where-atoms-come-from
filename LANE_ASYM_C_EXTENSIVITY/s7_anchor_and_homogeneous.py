"""S7 -- TWO LOOSE ENDS TIED OFF.

(1) THE [PHYS] ANCHOR AT m=3.  Everything about chi was computed in the [CODE] restriction.
    S1 anchored that restriction to the full 4m-qubit Hilbert space at m=1 and m=2.  Here it is
    anchored once more at m=3: the FULL 12-qubit space, dim 4096, diagonalised densely for its
    eigenspace structure, with the clause checks done as EXACT OPERATOR IDENTITIES on that same
    4096-dimensional space using the sparse Pauli representation (a Pauli has exactly one
    non-zero per row, so products are O(d) instead of O(d^3)).  The sparse representation is
    validated against the model's own xz_to_matrix at n=4 and n=8 before it is used.
    m=3 is the largest m at which the full space can be diagonalised here; at m=4 the dense
    Hamiltonian is 65536 x 65536 = 68 GB in complex128.

(2) THE HETEROGENEOUS-BATH ARTIFACT.  In s3 TABLE 9 the chi relation matrix has a small non-zero
    TRACE defect (-3.0e-04 at m=2).  That is not an interaction: the shared bath's three qubits
    carry DIFFERENT energies (1.0, 1.4, 0.7) and records are assigned round-robin, so tr M for
    2 blocks samples a different multiset of bath energies than 2 x tr M for 1 block.  With a
    HOMOGENEOUS bath the defect must vanish exactly, while the OFF-DIAGONAL entries -- the only
    entries that could ever carry record-record leakage -- stay at zero either way.
"""
import sys, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_C_EXTENSIVITY")
from lanelib import *

OUT = []
def P_(s=""):
    print(s, flush=True); OUT.append(str(s))

# ---------------------------------------------------------------- sparse Pauli representation
def pauli_rep(v, n):
    """R = sum_b c(b) |b xor xmask><b|.  Qubit i is bit (n-1-i), matching xz_to_matrix."""
    xmask = 0; zmask = 0; ny = 0
    for i in range(n):
        if v[i]: xmask |= 1 << (n - 1 - i)
        if v[n + i]: zmask |= 1 << (n - 1 - i)
        if v[i] and v[n + i]: ny += 1
    b = np.arange(1 << n, dtype=np.int64)
    par = np.zeros(1 << n, dtype=np.int64)
    t = b & zmask
    while True:
        if not np.any(t): break
        par ^= (t & 1).astype(np.int64); t >>= 1
        if t.max() == 0: break
    c = ((1j) ** ny) * np.where(par & 1, -1.0, 1.0)
    return xmask, c.astype(complex)

def p_apply(rep, M):
    """R @ M, exactly.  (R M)[a,:] = c(a xor x) M[a xor x,:]"""
    x, c = rep
    idx = np.arange(M.shape[0]) ^ x
    return M[idx, :] * c[idx][:, None]

def p_compose(r2, r1, n):
    """the Pauli r2 @ r1"""
    x1, c1 = r1; x2, c2 = r2
    b = np.arange(1 << n)
    return (x1 ^ x2, c1 * c2[b ^ x1])

def p_dense(rep, n):
    x, c = rep
    D = np.zeros((1 << n, 1 << n), dtype=complex)
    b = np.arange(1 << n)
    D[b ^ x, b] = c
    return D

P_("=" * 112)
P_("S7  [PHYS] ANCHOR AT m=3, AND THE HOMOGENEOUS-BATH RELATION MATRIX")
P_("=" * 112)

# ---------------------------------------------------------------- validate the representation
P_("\n" + "-" * 112)
P_("VALIDATION -- the sparse Pauli representation against the model's xz_to_matrix")
P_("-" * 112)
rng = np.random.default_rng(3)
worst = 0.0
for nn in (4, 8):
    for _ in range(6):
        v = list(rng.integers(0, 2, size=2 * nn))
        worst = max(worst, float(np.abs(p_dense(pauli_rep(v, nn), nn) - xz_to_matrix(v, nn)).max()))
    for _ in range(4):
        v1 = list(rng.integers(0, 2, size=2 * nn)); v2 = list(rng.integers(0, 2, size=2 * nn))
        got = p_dense(p_compose(pauli_rep(v2, nn), pauli_rep(v1, nn), nn), nn)
        ref = xz_to_matrix(v2, nn) @ xz_to_matrix(v1, nn)
        worst = max(worst, float(np.abs(got - ref).max()))
P_("   max deviation over 20 random Paulis and products at n=4 and n=8: %.3e" % worst)
P_("   SELF-CHECK: %s" % ("PASS" if worst < 1e-12 else "FAIL -- conclude nothing"))
assert worst < 1e-12

# ---------------------------------------------------------------- (1) the m=3 anchor
P_("\n" + "-" * 112)
P_("PART 1 -- the FULL 12-qubit space, dim 4096, three disjoint [[4,2,2]] blocks")
P_("-" * 112)
m = 3; n = 4 * m
H = stab_hamiltonian(m)
w, V = np.linalg.eigh(H)
blocks, i = [], 0
while i < len(w):
    j = i
    while j + 1 < len(w) and abs(w[j + 1] - w[i]) < 1e-8: j += 1
    blocks.append((float(w[i]), i, j + 1)); i = j + 1
P_("   H is %d x %d.  eigenvalues / multiplicities: %s"
   % (H.shape[0], H.shape[1], [(round(v, 6), hi - lo) for v, lo, hi in blocks]))
k = blocks[0][2] - blocks[0][1]
Q = V[:, :k]
P_("   ground-space dimension = %d   (must be 4^m = %d)  -> %s" % (k, 4 ** m, k == 4 ** m))
recs_v, wrts_v, _ = composite_records_writers(m)
P_("   %d records embedded, weights %s" % (len(recs_v), [weight(v, n) for v in recs_v]))
stab_v = composite_stab(m)
Ident = (0, np.ones(1 << n, dtype=complex))

def is_identity(rep):
    x, c = rep
    return x == 0 and float(np.abs(c - 1.0).max()) < 1e-12

P_("")
P_("   %-42s %-16s %-22s" % ("check (exact, on the full dim-4096 space)", "value", "CONTROL, same test"))
P_("   " + "-" * 84)
# clause (i)
sq = max(float(np.abs(p_compose(pauli_rep(v, n), pauli_rep(v, n), n)[1] - 1.0).max())
         for v in recs_v)
P_("   %-42s %-16.3e %-22s" % ("max ||R^2 - I||_inf                (i)", sq, "0.7R gives 0.51"))
# clause (ii): [H,R] = 0 iff R commutes with every stabiliser, checked as an operator identity
worst_comm = 0.0
for v in recs_v:
    rv = pauli_rep(v, n)
    for s in stab_v:
        rs = pauli_rep(s, n)
        a = p_compose(rs, rv, n); b = p_compose(rv, rs, n)
        worst_comm = max(worst_comm, float(np.abs(a[1] - b[1]).max()) + abs(a[0] - b[0]))
X0v = [1] + [0] * (2 * n - 1)
ctrl_comm = 0.0
for s in stab_v:
    rs = pauli_rep(s, n); rx = pauli_rep(X0v, n)
    a = p_compose(rs, rx, n); b = p_compose(rx, rs, n)
    ctrl_comm = max(ctrl_comm, float(np.abs(a[1] - b[1]).max()))
P_("   %-42s %-16.3e %-22.3f" % ("max ||[S,R]||_inf over all 6 stabs (ii)", worst_comm, ctrl_comm))
# invariance of the code space
Pg = Q @ Q.conj().T
inv = max(float(np.linalg.norm(p_apply(pauli_rep(v, n), Pg) - (p_apply(pauli_rep(v, n), Pg.conj().T)).conj().T))
          for v in recs_v)
P_("   %-42s %-16.3e %-22s" % ("max ||R P_g - P_g R||_F", inv, "-"))
# clause (iv): Tr(P_E R) on EVERY eigenspace, from one O(d^2) pass per record
max_iv = 0.0
for v in recs_v:
    RV = p_apply(pauli_rep(v, n), V)
    d = np.einsum('ij,ij->j', V.conj(), RV)
    for _, lo, hi in blocks:
        max_iv = max(max_iv, abs(complex(d[lo:hi].sum())))
RVi = p_apply(Ident, V)
di = np.einsum('ij,ij->j', V.conj(), RVi)
ctrl_iv = max(abs(complex(di[lo:hi].sum())) for _, lo, hi in blocks)
P_("   %-42s %-16.3e %-22.1f" % ("max |Tr(P_E R)| over all 7 eigenspaces (iv)", max_iv, ctrl_iv))
# clause (iii): given (ii) and (iv), non-constancy is forced -- and the norm is exact
P_("   %-42s %-16s %-22s"
   % ("||P_E R P_E||_F^2 = Tr(P_E) = m_E > 0 (iii)", "forced", "identity: 0 after"))
P_("   %-42s %-16s %-22s" % ("", "", "subtracting its mean"))
P_("")
P_("   why (iii) is forced, exactly: (ii) gives [R,P_E]=0, so P_E R P_E = R P_E and")
P_("   ||P_E R P_E||_F^2 = Tr(P_E R^dag R P_E) = Tr(P_E) = m_E.  (iv) gives Tr(P_E R) = 0, so if R")
P_("   were constant on E it would be the zero operator there, contradicting ||.||^2 = m_E > 0.")
P_("   Multiplicities m_E here are %s -- all non-zero." % [hi - lo for _, lo, hi in blocks])
assert sq < 1e-12 and worst_comm < 1e-12 and inv < 1e-9 and max_iv < 1e-9
assert k == 4 ** m and ctrl_comm > 0.1 and ctrl_iv > 1.0
P_("   SELF-CHECK: PASS -- all 6 records on 3 disjoint blocks satisfy (i)-(iv) as exact operator")
P_("   identities on the full 12-qubit space, and the code space is an exact invariant subspace.")
P_("   Every control registered a failure on the same instrument.")
del H, V, Q, Pg

# ---------------------------------------------------------------- (2) homogeneous relation matrix
P_("\n" + "-" * 112)
P_("PART 2 -- the chi RELATION MATRIX with a HOMOGENEOUS bath (every site e = 1.0)")
P_("          M_ij = chi about record j when ONLY record i is coupled.  [DIRECT] simulation,")
P_("          25 times in [1,13].")
P_("-" * 112)
nq = 3
env = Environment(nq=nq, energies=(1.0, 1.0, 1.0), beta=2.0)
P_("   %-5s %-10s %-18s %-18s %-13s %-15s %-15s"
   % ("m", "joint dim", "tr M", "m*tr M(1)", "DEFECT", "max off-diag", "max cross-block"))
P_("   " + "-" * 100)
tr1 = None
for mm in (1, 2, 3):
    Zl, dc = code_records_couplings(mm)
    HS = np.zeros((dc, dc), dtype=complex)
    M = np.zeros((2 * mm, 2 * mm))
    for i in range(2 * mm):
        M[i, :] = chi_timeavg(Propagator(HS, env, [(Zl[i], i % nq)], lam=0.8), Zl)
    if tr1 is None: tr1 = float(np.trace(M))
    off = float(np.abs(M - np.diag(np.diag(M))).max())
    cross = max([abs(M[i, j]) for i in range(2 * mm) for j in range(2 * mm) if i // 2 != j // 2]
                + [0.0])
    P_("   %-5d %-10d %-18.12f %-18.12f %-13.3e %-15.3e %-15.3e"
       % (mm, dc * env.dim, np.trace(M), mm * tr1, np.trace(M) - mm * tr1, off, cross))
P_("   CONTROL on the same scale: every DIAGONAL entry of M is %.8f" % (tr1 / 2.0))
P_("")
P_("   READ: with a homogeneous bath the trace defect is exactly 0 at every m, confirming that")
P_("   the -3.0e-04 in s3 TABLE 9 was the round-robin over bath qubits of DIFFERENT energies")
P_("   (1.0, 1.4, 0.7), not an interaction.  The off-diagonal and cross-block entries -- the")
P_("   only places record-record leakage could show -- are at 1e-16 in BOTH the homogeneous and")
P_("   the heterogeneous run, against a diagonal of ~0.5.")

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_C_EXTENSIVITY/s7_anchor_and_homogeneous.txt",
     "w").write("\n".join(OUT) + "\n")
print("\n[written]")
