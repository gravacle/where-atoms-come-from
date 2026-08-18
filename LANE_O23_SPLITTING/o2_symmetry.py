"""LANE O-2/O-3, PART 1.  DOES AN ORDINARY (0-FORM) SYMMETRY DEGENERACY SPLIT AT FIRST ORDER?

W-61 MEASURED 2.0e-06 splitting at eps=1e-06, linear, on an 8x8 random Hermitian H0 with a two-fold
degeneracy and a fully generic Hermitian perturbation.  That is a measurement of GENERIC HERMITIAN
perturbation theory, NOT of a symmetry-protected degeneracy under a symmetry-breaking perturbation.
Row O-2 asks for the proof.  This script tests the hypotheses of the proof.

THE STATEMENT UNDER TEST (the honest form; the unqualified form is FALSE and this script shows why).

  Fix a degenerate eigenspace E0 of H0, projector P, n = dim E0 >= 2.  Fix a REAL VECTOR SPACE V of
  admissible Hermitian perturbations.  Define the FIRST-ORDER MAP

        Phi : V -> {traceless Hermitian operators on E0},   Phi(V) = PVP - (tr PVP / n) P.

  T1  (Rellich / degenerate perturbation theory) the n eigenvalues emerging from E0 under
      H(eps)=H0+eps V are E + eps*lam_i + O(eps^2) with lam_i = eig(PVP|E0).  Hence
             E0 SPLITS AT FIRST ORDER  <=>  Phi(V) != 0.
  T2  (Schur)  if E0 is G-isotypic carrying an IRREP rho and [V,U(g)]=0 for all g, then Phi(V)=0,
      and in fact the degeneracy is exact to ALL orders.
  T3  (SELECTION RULE)  decompose V into irreducible tensor components under the adjoint action of G.
      Phi(V) depends ONLY on the components lying in irreps lam contained in rho (x) rho*.
      Components in any other irrep are in ker Phi EXACTLY.  This is where the first-order matrix
      element VANISHES BY A FURTHER SYMMETRY -- and the "further symmetry" is G itself.
  T4  (GENERICITY, made precise)  ker Phi is a linear SUBSPACE of V of codimension r = rank Phi.
      So "generic" can only ever mean "outside a codimension-r subspace", and the claim
      HAS CONTENT ONLY RELATIVE TO V.  first-order splitting is generic in V  <=>  r >= 1.
      If r = 0 NO perturbation in V splits at first order, however generic.
  T5  (LOCALITY -- the half that makes C-2 true)  if G is a compact CONNECTED group whose generators
      are sums of operators on contractible regions, Q^a = sum_x q^a_x, and rho is nontrivial, then
      r >= 1 already for V = {contractible-support operators}.  A LOCAL perturbation splits at first
      order.  PROOF: P Q^a P = Q^a|E0 = d.rho(X^a), which is traceless and (for rho nontrivial on the
      identity component) nonzero for some a.  If every P q^a_x P were a scalar their sum would be a
      scalar; a nonzero traceless operator is not.  Hence some single term q^a_x has Phi != 0.  QED

TESTS RUN HERE
  SELF-CHECK 0   a case with a hand-computable answer: H0 = 0 on C^2, V = sigma_z -> splitting = 2 eps
                 exactly, slope 1.000.
  SELF-CHECK 1   the adjoint-Casimir tensor decomposition reproduces the rep-theoretic multiplicities
                 of End(H) for two independently known cases.
  TEST A  (T2,T3) spin-1 (x) spin-1/2, E0 = the j=1/2 DOUBLET.  rho (x) rho* = 0 + 1, so only RANK-1
                 components can split it.  rank-1 must split; rank-2 (a crystal field) must NOT,
                 exactly, at any strength.  POSITIVE CONTROL: the rank-1 run.
  TEST B  (T3)   spin-1 (x) spin-1, E0 = the j=1 TRIPLET, INTEGER spin so there is no Kramers
                 theorem in play.  rho (x) rho* = 0+1+2, so a RANK-3 perturbation must give
                 Phi = 0 exactly and NO first-order splitting.  POSITIVE CONTROL: rank-2 -> slope 1.
  TEST C  (T4)   measure rank Phi for V = all Hermitians:  must equal n^2 - 1 exactly.
  TEST D  (T5)   measure rank Phi for V = single-site operators on a genuine many-body 0-form
                 symmetry (SU(2) Heisenberg ring) and on an SSB Ising ring:  must be >= 1.
"""
import numpy as np
np.set_printoptions(precision=4, suppress=True)
rng = np.random.default_rng(2026)

# ---------------------------------------------------------------- basic tools
I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)


def spin_ops(twoj):
    """spin matrices for spin j = twoj/2, dimension twoj+1."""
    j = twoj / 2.0
    d = twoj + 1
    m = np.array([j - k for k in range(d)])
    Jz = np.diag(m).astype(complex)
    Jp = np.zeros((d, d), complex)
    for k in range(1, d):
        Jp[k - 1, k] = np.sqrt(j * (j + 1) - m[k] * (m[k] + 1))
    Jm = Jp.conj().T
    return (Jp + Jm) / 2, (Jp - Jm) / (2j), Jz


def kron(*ms):
    out = np.array([[1.0 + 0j]])
    for m in ms:
        out = np.kron(out, m)
    return out


def herm(M):
    return (M + M.conj().T) / 2


def eigencluster(H, P0, n):
    """the n eigenvalues of H whose eigenvectors have the largest overlap with range(P0)."""
    ev, U = np.linalg.eigh(H)
    w = np.einsum('ij,jk,ki->i', U.conj().T, P0, U).real
    order = np.argsort(-w)[:n]
    return np.sort(ev[order])


def splitting(H, P0, n):
    e = eigencluster(H, P0, n)
    return e[-1] - e[0]


def slope(H0, V, P0, n, eps_list):
    """log-log slope of splitting vs eps."""
    xs, ys = [], []
    for eps in eps_list:
        s = splitting(H0 + eps * V, P0, n)
        if s > 1e-13:
            xs.append(np.log(eps))
            ys.append(np.log(s))
    if len(xs) < 2:
        return float('nan'), 0.0
    A = np.polyfit(xs, ys, 1)
    resid = np.max(np.abs(np.polyval(A, xs) - ys))
    return A[0], resid


def phi(V, Q):
    """traceless part of PVP written in the E0 basis Q (columns = orthonormal basis of E0)."""
    B = Q.conj().T @ V @ Q
    n = B.shape[0]
    return B - np.trace(B) / n * np.eye(n)


def first_order_rank(basisV, Q, tol=1e-9):
    """rank of Phi restricted to the real span of basisV (a list of Hermitian matrices)."""
    rows = []
    n = Q.shape[1]
    for V in basisV:
        A = phi(V, Q)
        # real coordinates on traceless Hermitian n x n  (dimension n^2-1)
        v = []
        for a in range(n):
            for b in range(a, n):
                if a == b:
                    v.append(A[a, a].real)
                else:
                    v.append(A[a, b].real)
                    v.append(A[a, b].imag)
        rows.append(v)
    M = np.array(rows)
    if M.size == 0:
        return 0, 0.0
    s = np.linalg.svd(M, compute_uv=False)
    smax = s.max() if s.size else 0.0
    r = int(np.sum(s > tol * max(1.0, smax)))
    return r, smax


def degenerate_space(H, which=0, tol=1e-9):
    """return (energy, orthonormal basis Q, projector P) of the `which`-th distinct eigenvalue."""
    ev, U = np.linalg.eigh(H)
    groups = []
    i = 0
    while i < len(ev):
        j = i
        while j + 1 < len(ev) and abs(ev[j + 1] - ev[i]) < tol * max(1.0, abs(ev[i])):
            j += 1
        groups.append((ev[i], U[:, i:j + 1]))
        i = j + 1
    E, Q = groups[which]
    return E, Q, Q @ Q.conj().T, len(groups)


# ------------------------------------------- adjoint-Casimir tensor decomposition
def adjoint_rank_projector(Sops, dim, k):
    """projector (as a list of basis operators) onto the rank-k irreducible tensor operators.

    The adjoint action of SU(2) on operators has Casimir  J^2_ad(X) = sum_a [S^a,[S^a,X]] with
    eigenvalue k(k+1) on a rank-k tensor operator.  Diagonalise the superoperator.
    """
    N = dim
    # superoperator on the N^2-dimensional operator space, vec(X) column-major consistent with reshape
    Sup = np.zeros((N * N, N * N), complex)
    Id = np.eye(N, dtype=complex)
    for S in Sops:
        # [S,[S,X]] = S S X - 2 S X S + X S S
        Sup += np.kron(S @ S, Id) - 2 * np.kron(S, S.T) + np.kron(Id, (S @ S).T)
    Sup = herm(Sup)
    w, Uv = np.linalg.eigh(Sup)
    target = k * (k + 1)
    sel = np.abs(w - target) < 1e-7 * max(1.0, target)
    cols = Uv[:, sel]
    ops = [cols[:, c].reshape(N, N) for c in range(cols.shape[1])]
    return ops, w


def rank_multiplicities(Sops, dim):
    _, w = adjoint_rank_projector(Sops, dim, 0)
    out = {}
    for val in w:
        k = int(round((-1 + np.sqrt(1 + 4 * max(val, 0))) / 2))
        out[k] = out.get(k, 0) + 1
    return out


def random_rank_k_herm(Sops, dim, k, rng):
    ops, _ = adjoint_rank_projector(Sops, dim, k)
    if not ops:
        return None
    M = sum(rng.normal() * o for o in ops)
    M = herm(M)
    # herm() can project out of the rank-k space if the space is not closed under dagger;
    # for SU(2) the rank-k space IS closed under dagger up to a sign convention, verified below.
    return M


print("=" * 100)
print("LANE O-2 PART 1.  FIRST-ORDER SPLITTING OF AN ORDINARY (0-FORM) SYMMETRY DEGENERACY")
print("=" * 100)

# ============================================================= SELF-CHECK 0
print("\nSELF-CHECK 0 -- a case whose answer is known by hand.")
H0 = np.zeros((2, 2), complex)
V = SZ.copy()
_, Q0, P0, _ = degenerate_space(H0)
eps_list = np.array([1e-6, 1e-5, 1e-4, 1e-3])
sl, res = slope(H0, V, P0, 2, eps_list)
exact_ok = all(abs(splitting(H0 + e * V, P0, 2) - 2 * e) < 1e-12 * max(1.0, 2 * e) for e in eps_list)
print(f"    H0 = 0 on C^2, V = sigma_z.  exact answer: splitting = 2*eps, slope 1.")
print(f"    measured slope = {sl:.6f}   max log-log residual = {res:.2e}")
print(f"    splitting == 2*eps to machine precision at every eps: {exact_ok}")
print(f"    SELF-CHECK 0: {'PASS' if abs(sl - 1) < 1e-6 and exact_ok else 'FAIL'}")

# ============================================================= SELF-CHECK 1
print("\nSELF-CHECK 1 -- the adjoint-Casimir tensor decomposition against known multiplicities.")
# case (a): a single spin-1.  End(C^3) = rank 0 + rank 1 + rank 2, dims 1+3+5 = 9.
S1 = spin_ops(2)
m_a = rank_multiplicities(list(S1), 3)
exp_a = {0: 1, 1: 3, 2: 5}
# case (b): spin-1 (x) spin-1/2 = j=3/2 + j=1/2.  End = (3/2 + 1/2)(x)(3/2 + 1/2)*
#   = (3/2 x 3/2) + (3/2 x 1/2) + (1/2 x 3/2) + (1/2 x 1/2)
#   = (0+1+2+3) + (1+2) + (1+2) + (0+1)
#   ranks: 0 twice(dim2), 1 four times(dim12), 2 three times(dim15), 3 once(dim7) -> 2+12+15+7 = 36
SA = [kron(s, I2) for s in spin_ops(2)]
SB = [kron(np.eye(3, dtype=complex), s) for s in spin_ops(1)]
Stot6 = [SA[i] + SB[i] for i in range(3)]
m_b = rank_multiplicities(Stot6, 6)
exp_b = {0: 2, 1: 12, 2: 15, 3: 7}
print(f"    single spin-1, End(C^3):        measured {dict(sorted(m_a.items()))}   expected {exp_a}")
print(f"    spin-1 (x) spin-1/2, End(C^6):  measured {dict(sorted(m_b.items()))}   expected {exp_b}")
ok1 = (m_a == exp_a) and (m_b == exp_b)
print(f"    SELF-CHECK 1: {'PASS' if ok1 else 'FAIL'}")

# ============================================================= TEST A
print("\n" + "-" * 100)
print("TEST A (T2, T3).  E0 = the j=1/2 DOUBLET of spin-1 (x) spin-1/2.   n = 2, rho = spin-1/2.")
print("   rho (x) rho* = rank 0 + rank 1.  SELECTION RULE: only a RANK-1 component can split it.")
print("-" * 100)
Ssq6 = sum(S @ S for S in Stot6)
Hsym = herm(Ssq6)
Eg, Q6, P6, ngroups = degenerate_space(Hsym)
print(f"    H0 = (S_tot)^2.  distinct levels = {ngroups}.  lowest level E = {Eg:.4f}, degeneracy n = {Q6.shape[1]}")
assert Q6.shape[1] == 2, "expected the j=1/2 doublet at the bottom"
n = 2

# verify the rank-k spaces are dagger-closed so that herm() does not leak between ranks
for k in (1, 2):
    ops, _ = adjoint_rank_projector(Stot6, 6, k)
    Pk = np.zeros((36, 36), complex)
    for o in ops:
        v = o.reshape(-1, 1)
        Pk += v @ v.conj().T
    leak = 0.0
    for o in ops:
        d = o.conj().T
        vd = d.reshape(-1, 1)
        leak = max(leak, float(np.linalg.norm(vd - Pk @ vd)))
    print(f"    rank-{k} operator space is closed under dagger: leak = {leak:.2e}  "
          f"{'PASS' if leak < 1e-9 else 'FAIL'}")

cases = []
V1 = random_rank_k_herm(Stot6, 6, 1, rng)
V2 = random_rank_k_herm(Stot6, 6, 2, rng)
V3 = random_rank_k_herm(Stot6, 6, 3, rng)
Vinv = herm(sum(SA[i] @ SB[i] for i in range(3)))          # G-invariant (rank 0)
Vgen = herm(rng.normal(size=(6, 6)) + 1j * rng.normal(size=(6, 6)))
for name, V in [("rank-1 (magnetic, T-odd)", V1),
                ("rank-2 (crystal field)", V2),
                ("rank-3", V3),
                ("rank-0 = G-INVARIANT S_A.S_B", Vinv),
                ("generic Hermitian (breaks everything)", Vgen)]:
    if V is None:
        continue
    V = V / np.linalg.norm(V) * np.linalg.norm(Hsym)
    ph = np.linalg.norm(phi(V, Q6))
    eps_l = np.array([1e-5, 1e-4, 1e-3, 1e-2])
    sl, res = slope(Hsym, V, P6, n, eps_l)
    s_at = splitting(Hsym + 1e-3 * V, P6, n)
    cases.append((name, ph, s_at, sl))
    print(f"    {name:<40s} ||Phi(V)|| = {ph:10.3e}   splitting(eps=1e-3) = {s_at:10.3e}   slope = {sl:.4f}")
print("    EXPECTED: rank-1 and generic -> Phi != 0, slope 1.  rank-2, rank-3, invariant -> Phi = 0 EXACTLY.")
a_ok = (cases[0][1] > 1e-6 and cases[0][3] > 0.99 and                    # rank 1 splits, POSITIVE CONTROL
        cases[1][1] < 1e-12 and cases[2][1] < 1e-12 and cases[3][1] < 1e-12)
print(f"    TEST A: {'PASS' if a_ok else 'FAIL'}")

print("\n    NOTE, AND IT MATTERS FOR THE CONTRAST O-2 IS TRYING TO DRAW:")
print("    the rank-2 perturbation has Phi = 0 to machine precision yet the doublet splits with")
print(f"    SLOPE {cases[1][3]:.4f}, i.e. at SECOND order.  A 0-form symmetry degeneracy can therefore")
print("    have splitting exponent 2.  The exponent alone does NOT distinguish 0-form from 1-form.")

print("\n    KRAMERS SUB-TEST.  the system has ONE half-integer spin, so T^2 = -1 on the whole space")
print("    and the j=1/2 level is a KRAMERS doublet.  Kramers' theorem predicts EXACT degeneracy at")
print("    every strength for any T-INVARIANT V.  Note SU(2)-rank does NOT fix T-parity here (the")
print("    rank-2 space mixes the j-blocks), so T-parity must be imposed explicitly.")
def expm_herm(A, coeff):
    """exp(coeff * A) for Hermitian A, via its spectral decomposition (no scipy)."""
    w, U = np.linalg.eigh(herm(A))
    return U @ np.diag(np.exp(coeff * w)) @ U.conj().T


Uy = expm_herm(Stot6[1], -1j * np.pi)       # T = Uy * (complex conjugation)


def T_conj(V):
    return Uy @ V.conj() @ Uy.conj().T


print(f"      check T^2 = -1 on the whole space: ||Uy conj(Uy) + I|| = "
      f"{np.linalg.norm(Uy @ Uy.conj() + np.eye(6)):.2e}")
print(f"      check H0 is T-invariant:           ||T H T^-1 - H||    = "
      f"{np.linalg.norm(T_conj(Hsym) - Hsym):.2e}")
Vg = herm(rng.normal(size=(6, 6)) + 1j * rng.normal(size=(6, 6)))
Veven = herm(Vg + T_conj(Vg))
Vodd = herm(Vg - T_conj(Vg))
Veven = Veven / np.linalg.norm(Veven) * np.linalg.norm(Hsym)
Vodd = Vodd / np.linalg.norm(Vodd) * np.linalg.norm(Hsym)
print(f"      ||Phi(V_T-even)|| = {np.linalg.norm(phi(Veven, Q6)):.3e}    "
      f"||Phi(V_T-odd)|| = {np.linalg.norm(phi(Vodd, Q6)):.3e}  <- POSITIVE CONTROL")
kr_ok = True
for eps in (1e-3, 1e-2, 1e-1, 0.5, 1.0):
    se = splitting(Hsym + eps * Veven, P6, n)
    so = splitting(Hsym + eps * Vodd, P6, n)
    kr_ok &= se < 1e-10
    print(f"      eps = {eps:5.3f}:  T-EVEN splitting = {se:.3e}     "
          f"T-ODD splitting = {so:.3e}   <- POSITIVE CONTROL")
print(f"    KRAMERS SUB-TEST: {'PASS' if kr_ok else 'FAIL'}  -- a T-invariant perturbation, however")
print("    large and however badly it breaks SU(2), does not split the doublet at ANY order.")

# ============================================================= TEST B
print("\n" + "-" * 100)
print("TEST B (T3).  E0 = the j=1 TRIPLET of spin-1 (x) spin-1.  n = 3, rho = spin-1, INTEGER spin,")
print("   so NO Kramers theorem.  rho (x) rho* = 0+1+2.  A RANK-3 perturbation must give Phi = 0.")
print("-" * 100)
SA9 = [kron(s, np.eye(3, dtype=complex)) for s in spin_ops(2)]
SB9 = [kron(np.eye(3, dtype=complex), s) for s in spin_ops(2)]
St9 = [SA9[i] + SB9[i] for i in range(3)]
mult9 = rank_multiplicities(St9, 9)
# 9 = j2 + j1 + j0.  End = sum over pairs (ja,jb) of ja(x)jb.  Counting DIMENSIONS by rank k:
#   k=0: from 2x2,1x1,0x0                      = 1+1+1                 = 3
#   k=1: from 2x2,2x1,1x2,1x1,1x0,0x1          = 3*6                   = 18
#   k=2: from 2x2,2x1,1x2,2x0,0x2,1x1          = 5*6                   = 30
#   k=3: from 2x2,2x1,1x2                      = 7*3                   = 21
#   k=4: from 2x2                              = 9                     = 9      total 81
exp9 = {0: 3, 1: 18, 2: 30, 3: 21, 4: 9}
print(f"    End(C^9) tensor multiplicities: {dict(sorted(mult9.items()))}   expected {exp9}   "
      f"{'PASS' if mult9 == exp9 else 'FAIL'}   (sum = {sum(mult9.values())}, must be 81)")
Ssq9 = herm(sum(S @ S for S in St9))
# pick out the j=1 triplet: eigenvalue j(j+1) = 2
ev9, U9 = np.linalg.eigh(Ssq9)
sel = np.abs(ev9 - 2.0) < 1e-8
Q9 = U9[:, sel]
P9 = Q9 @ Q9.conj().T
print(f"    j=1 triplet found: dimension {Q9.shape[1]} at S^2 = {ev9[sel][0]:.4f}")
H9 = Ssq9
for k in (1, 2, 3, 4):
    Vk = random_rank_k_herm(St9, 9, k, rng)
    if Vk is None:
        continue
    Vk = Vk / np.linalg.norm(Vk) * np.linalg.norm(H9)
    ph = np.linalg.norm(phi(Vk, Q9))
    sl, res = slope(H9, Vk, P9, 3, np.array([1e-5, 1e-4, 1e-3, 1e-2]))
    print(f"    rank-{k} perturbation:  ||Phi(V)|| = {ph:10.3e}   slope = {sl:.4f}   "
          f"splitting(1e-3) = {splitting(H9 + 1e-3 * Vk, P9, 3):.3e}")
# explicit rank-3 statistics over many draws -- a zero must not be a lucky draw
ops3, _ = adjoint_rank_projector(St9, 9, 3)
worst = 0.0
for _ in range(200):
    Vk = herm(sum(rng.normal() * o for o in ops3))
    Vk = Vk / max(np.linalg.norm(Vk), 1e-300)
    worst = max(worst, float(np.linalg.norm(phi(Vk, Q9))))
print(f"    200 random RANK-3 perturbations: max ||Phi(V)|| over all draws = {worst:.3e}")
ops2, _ = adjoint_rank_projector(St9, 9, 2)
best2 = 0.0
for _ in range(200):
    Vk = herm(sum(rng.normal() * o for o in ops2))
    Vk = Vk / max(np.linalg.norm(Vk), 1e-300)
    best2 = max(best2, float(np.linalg.norm(phi(Vk, Q9))))
print(f"    200 random RANK-2 perturbations: max ||Phi(V)|| over all draws = {best2:.3e}  <- POSITIVE CONTROL")
b_ok = worst < 1e-10 and best2 > 1e-3
print(f"    TEST B: {'PASS' if b_ok else 'FAIL'}   (rank-3 vanishes identically; rank-2 does not)")

# ============================================================= TEST C
print("\n" + "-" * 100)
print("TEST C (T4).  rank of Phi on V = ALL Hermitian operators.  Theory: exactly n^2 - 1.")
print("-" * 100)
for label, Q in [("doublet n=2 (dim-6 space)", Q6), ("triplet n=3 (dim-9 space)", Q9)]:
    N = Q.shape[0]
    basis = []
    for a in range(N):
        Eaa = np.zeros((N, N), complex); Eaa[a, a] = 1; basis.append(Eaa)
        for b in range(a + 1, N):
            M = np.zeros((N, N), complex); M[a, b] = 1; M[b, a] = 1; basis.append(M)
            M = np.zeros((N, N), complex); M[a, b] = -1j; M[b, a] = 1j; basis.append(M)
    r, _ = first_order_rank(basis, Q)
    nq = Q.shape[1]
    print(f"    {label:<28s}  dim V = {len(basis):4d}   rank Phi = {r:3d}   n^2-1 = {nq**2 - 1}   "
          f"{'PASS' if r == nq**2 - 1 else 'FAIL'}")
print("    => the non-splitting set is a SUBSPACE of codimension n^2-1 >= 3: closed, measure zero,")
print("       nowhere dense.  THAT is the precise content of the word 'generic' in row O-2.")

# ============================================================= TEST D
print("\n" + "-" * 100)
print("TEST D (T5).  MANY-BODY 0-FORM SYMMETRY: is there a CONTRACTIBLE (single-site) perturbation")
print("   with Phi != 0?  Theory says yes for connected G whose generators are sums of local terms.")
print("-" * 100)


def site_op(op, k, L):
    ms = [I2] * L
    ms[k] = op
    return kron(*ms)


def heisenberg_ring(L, J=1.0):
    H = np.zeros((2 ** L, 2 ** L), complex)
    for k in range(L):
        for s in (SX, SY, SZ):
            H += J * site_op(s, k, L) @ site_op(s, (k + 1) % L, L) / 4
    return herm(H)


d1_report = {}
for L, J, tag in [(4, -1.0, "FERROMAGNETIC (ground multiplet = spin L/2 irrep)"),
                  (6, -1.0, "FERROMAGNETIC (ground multiplet = spin L/2 irrep)"),
                  (5, +1.0, "ANTIFERROMAGNETIC odd ring (ground multiplet has half-integer spin)")]:
    H = heisenberg_ring(L, J)
    E, Q, P, ng = degenerate_space(H)
    nq = Q.shape[1]
    onesite = [site_op(s, k, L) for k in range(L) for s in (SX, SY, SZ)]
    r_local, smax = first_order_rank(onesite, Q)
    print(f"    Heisenberg ring L={L}, J={J:+.0f}  {tag}")
    print(f"      ground level E={E:.6f}, degeneracy n={nq}   (expected multiplet dim "
          f"{L + 1 if J < 0 else '2 or 4'})")
    print(f"      rank Phi on SINGLE-SITE operators = {r_local}   n^2-1 = {nq**2 - 1}")
    print(f"      (when E0 is a MULTIPLICITY-FREE irrep, single-site ops are rank-1 tensors and by")
    print(f"       Wigner-Eckart reach only the rank-1 part of End(E0), so rank Phi = 3 exactly;")
    print(f"       with multiplicity > 1 they also reach the multiplicity indices and rank Phi is larger.)")
    if nq > 1:
        best = max(np.linalg.norm(phi(V, Q)) for V in onesite)
        Vloc = max(onesite, key=lambda V: np.linalg.norm(phi(V, Q)))
        sl, res = slope(H, Vloc, P, nq, np.array([1e-6, 1e-5, 1e-4, 1e-3]))
        d1_report[f"Heis L={L} J={J:+.0f}"] = (best, sl)
        print(f"      FIRST-ORDER LOCAL DEFECT  D1 = max over single-site V of ||Phi(V)|| = {best:.4f}")
        print(f"      splitting slope for that single-site V = {sl:.4f}   (log-log residual {res:.1e})")
        print(f"      A SINGLE-SITE OPERATOR SPLITS THE MULTIPLET AT FIRST ORDER: "
              f"{'CONFIRMED' if best > 1e-6 and abs(sl - 1) < 1e-3 else 'NOT CONFIRMED'}")
        # the explicit T5 construction: the generator itself, broken into its local terms
        Sz_tot = sum(site_op(SZ, k, L) for k in range(L)) / 2
        pg = phi(Sz_tot, Q)
        loc_terms = [phi(site_op(SZ, k, L) / 2, Q) for k in range(L)]
        print(f"      T5 CONSTRUCTION: ||Phi(Q^z_total)|| = {np.linalg.norm(pg):.4f} (nonzero, traceless: "
              f"|tr| = {abs(np.trace(pg)):.1e});  its L local terms have ||Phi|| = "
              f"{[round(float(np.linalg.norm(t)), 4) for t in loc_terms]}")

print("\n    SSB CONTROL -- Ising ring (finite group Z2, spontaneously broken).")
for L in (4, 6):
    H = np.zeros((2 ** L, 2 ** L), complex)
    for k in range(L):
        H += -site_op(SZ, k, L) @ site_op(SZ, (k + 1) % L, L)
    H = herm(H)
    E, Q, P, ng = degenerate_space(H)
    nq = Q.shape[1]
    onesite = [site_op(s, k, L) for k in range(L) for s in (SX, SY, SZ)]
    r_local, smax = first_order_rank(onesite, Q)
    best = max(np.linalg.norm(phi(V, Q)) for V in onesite)
    Vloc = max(onesite, key=lambda V: np.linalg.norm(phi(V, Q)))
    sl, _ = slope(H, Vloc, P, nq, np.array([1e-6, 1e-5, 1e-4, 1e-3]))
    print(f"    Ising ring L={L}: n={nq}, rank Phi on single-site ops = {r_local}, "
          f"best ||Phi|| = {best:.4f}, slope = {sl:.4f}")

# ============================================================= TEST E
print("\n" + "-" * 100)
print("TEST E.  WHERE T5's PROOF BREAKS: FINITE GROUPS.  T5 used the LIE GENERATORS Q^a = sum_x q^a_x.")
print("   A FINITE group has no generators, so the proof does not apply and the conclusion is not")
print("   established for finite G.  This is a PROBE, not a proof: sweep random G-invariant")
print("   Hamiltonians for G = S_3 acting ON-SITE in its 2-dim standard irrep, and ask whether any")
print("   of them has a degenerate ground multiplet that NO single-site operator splits at first")
print("   order (i.e. first-order local defect D1 = 0 with n >= 2).")
print("-" * 100)
th = 2 * np.pi / 3
r2 = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]], dtype=complex)
s2 = np.array([[1, 0], [0, -1]], dtype=complex)
S3 = [np.eye(2, dtype=complex), r2, r2 @ r2, s2, r2 @ s2, r2 @ r2 @ s2]
# confirm it is a group of order 6 and non-abelian
prods = []
for a in S3:
    for b in S3:
        prods.append(a @ b)
closed = all(min(np.linalg.norm(p - g) for g in S3) < 1e-10 for p in prods)
nonab = max(np.linalg.norm(a @ b - b @ a) for a in S3 for b in S3) > 1e-6
print(f"    S_3 standard rep: closed under multiplication = {closed}, non-abelian = {nonab}   "
      f"{'PASS' if closed and nonab else 'FAIL'}")
for L in (3, 4):
    dim = 2 ** L
    Ug = [kron(*([g] * L)) for g in S3]
    onesite = [site_op(s, k, L) for k in range(L) for s in (SX, SY, SZ)]
    found_zero = 0
    tried = 0
    worst_D1 = None
    dist = {}
    for _ in range(300):
        A = herm(rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim)))
        Hinv = sum(U @ A @ U.conj().T for U in Ug) / len(Ug)
        Hinv = herm(Hinv)
        chk = max(np.linalg.norm(U @ Hinv - Hinv @ U) for U in Ug)
        if chk > 1e-9:
            continue
        E, Qg, Pg, ng = degenerate_space(Hinv)
        nq = Qg.shape[1]
        dist[nq] = dist.get(nq, 0) + 1
        if nq < 2:
            continue
        tried += 1
        D1 = max(float(np.linalg.norm(phi(V, Qg))) for V in onesite)
        worst_D1 = D1 if worst_D1 is None else min(worst_D1, D1)
        if D1 < 1e-9:
            found_zero += 1
    print(f"    L = {L} sites (dim {dim}): ground-multiplet dimensions seen = {dict(sorted(dist.items()))}")
    print(f"      degenerate draws examined = {tried};  draws with D1 = 0 (no single-site splitter) "
          f"= {found_zero}")
    print(f"      SMALLEST first-order local defect D1 over all degenerate draws = "
          f"{('%.4f' % worst_D1) if worst_D1 is not None else 'n/a'}")
print("    READ THIS CORRECTLY: a sweep finding no counterexample is NOT a proof.  T5 is PROVED for")
print("    connected compact G and is a CONJECTURE for finite G.  What would settle it: show that a")
print("    finite on-site 0-form symmetry cannot produce a ground multiplet with no local order")
print("    parameter -- equivalently that d_R = 1 always.  That is the missing step, and it is the")
print("    ONLY thing standing between this lane and a complete proof of C-2.")

print("\n" + "=" * 100)
print("WHAT PART 1 ESTABLISHES")
print("=" * 100)
print("""  (1) O-2's claim is TRUE ONLY IN THE FORM 'rank Phi >= 1 on the admissible space V'.  Stated as
      'a symmetry-breaking perturbation splits the multiplet at first order' it is FALSE: TEST A and
      TEST B exhibit perturbations that break the symmetry maximally and give Phi = 0 IDENTICALLY,
      by the selection rule lam subset rho (x) rho*.  These are not measure-zero accidents -- they
      are whole subspaces of perturbations (200/200 random draws gave zero).
  (2) The genericity statement that IS true: the non-splitting perturbations form a linear subspace
      ker Phi of codimension rank Phi.  Generic = outside that subspace.  The claim is empty unless
      V is named.
  (3) The half that makes C-2 true is T5, LOCALITY: for a 0-form symmetry the generators are SUMS OF
      LOCAL TERMS, so some single local term is already non-scalar on the multiplet.  TEST D measures
      rank Phi >= 1 on single-site operators for a genuine many-body 0-form symmetry.  This -- not
      genericity -- is what kills 0-form symmetry as a record carrier.""")
