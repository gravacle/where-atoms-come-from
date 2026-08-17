# LANE W11-R-CROSS  LEG D -- THE STEELMAN'S CHARACTERISATION IS FALSE, AND ITS 4000-DRAW CONTROL
#                            COULD NOT HAVE FOUND THE COUNTEREXAMPLE.
#
# CLAIM UNDER ATTACK (steelman leg D / leg G1, quoted): "CHARACTERISATION (tested both ways, 4000
# draws): Z_n is pi-only for all n IFF both roots are DIAGONAL AND CONSTANT ON THE FOUR CLASSES ...
# 4000 generic draws: 4 pi-only, all 4 class-constant diagonal, 0 counterexamples", and in
# what_falls: "'Every non-class-constant root breaks it' is TRUE.  The preserving set is measure zero."
#
# THE ONLY-IF DIRECTION IS FALSE.  Exhibited here: a POSITIVE-DIMENSIONAL family of root pairs that
# are NON-DIAGONAL, NON-CLASS-CONSTANT, and pi-only at EVERY tick.
#   Let R be any unitary that is block diagonal on the FOUR CLASSES with R^gcd(L_F,L_C) = I, and set
#       U_F = Lambda_F . R ,   U_C = Lambda_C . R          (the SAME R in both branches)
#   with Lambda_F = diag(W_F^{1/L_F} on gamma_F, 1 off) and likewise Lambda_C.  Then
#       U_F^{L_F} = M_F,  U_C^{L_C} = M_C  exactly, and
#       K_n := (U_F*)^n U_C^n = conj(Lambda_F)^n Lambda_C^n . (R*)^n R^n = conj(Lambda_F)^n Lambda_C^n
#   which is class-constant diagonal for every n -- so Z_n = <U_F^n s, U_C^n s> is pi-only, while
#   U_F and U_C are neither diagonal nor class-constant.
# THE CORRECT CRITERION IS ON THE PAIR, NOT ON EACH ROOT: pi-only for all n <=> K_n is class-constant
# diagonal for all n.  Checked both ways below.
#
# WHY THE 4000-DRAW SWEEP MISSED IT: it draws U_F and U_C INDEPENDENTLY.  The preserving set is a
# CORRELATION between the two branches, which an independent sampler cannot produce at any sample
# size.  On W-08's own norm that control COULD NOT HAVE FAILED, and its 0/4000 is a property of the
# sampler, not of the variety.
import numpy as np
import xlib as X

rng = np.random.default_rng(20260817)
lf, lc, NV = X.K1_LOOP_F, X.K1_LOOP_C, 5
a = np.array([1.0, 0.37, 0.91, 2 ** 0.5, 0.23, 1.77])
LF = LC = 3
MF, MC = X.M_circuit(lf, a, NV), X.M_circuit(lc, a, NV)
CLS = X.classes(lf, lc, NV)
print(f"  K1 classes by vertex: {CLS}    (11={{v0}}, 10={{v1,v2}}, 01={{v3,v4}})")


def class_blocks(NVv):
    d = {}
    for v in range(NVv):
        d.setdefault(CLS[v], []).append(v)
    return d


def rand_order_m_unitary(d, m, rng):
    """a unitary on C^d with V^m = I, generically NON-diagonal and non-scalar."""
    Q, R_ = np.linalg.qr(rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)))
    Q = Q @ np.diag(np.diag(R_) / np.abs(np.diag(R_)))
    ks = rng.integers(0, m, size=d)
    while len(set(ks.tolist())) == 1 and d > 1:          # force non-scalar
        ks = rng.integers(0, m, size=d)
    return Q @ np.diag(np.exp(2j * np.pi * ks / m)) @ Q.conj().T


def make_R(rng, m=3):
    """block diagonal on the four classes, R^m = I, generically non-diagonal."""
    R = np.zeros((NV, NV), dtype=complex)
    for c, vs in class_blocks(NV).items():
        B = rand_order_m_unitary(len(vs), m, rng)
        for i, v in enumerate(vs):
            for j, u in enumerate(vs):
                R[v, u] = B[i, j]
    return R


def lam(loop, aa, NVv):
    Lg = len(loop)
    w = np.exp(1j * np.angle(X.holonomy(loop, aa)) / Lg)
    L_ = np.eye(NVv, dtype=complex)
    for v in X.loop_vertices(loop):
        L_[v, v] = w
    return L_


def is_class_const_diag(U, tol=1e-9):
    if not np.allclose(U, np.diag(np.diag(U)), atol=tol):
        return False
    d = np.diag(U)
    seen = {}
    for v in range(len(d)):
        c = CLS[v]
        if c in seen and abs(seen[c] - d[v]) > tol:
            return False
        seen[c] = d[v]
    return True


w_base = np.array([0.40, 0.15, 0.15, 0.15, 0.15])
STATES = [np.sqrt(w_base) + 0j] + X.random_pi_identical(rng, lf, lc, NV, w_base, k=40)
pis = [X.pi_of(s, lf, lc, NV) for s in STATES]
assert all(np.allclose(pis[0], p, atol=1e-12) for p in pis), "pi not held fixed"
assert X.arms_differ(*STATES[:6]), "STATE ARMS BYTE-IDENTICAL -- leg void"


def spread(uF, uC, nmax=12):
    w = 0.0
    for n in range(1, nmax + 1):
        v = [abs(X.Z(uF, uC, s, n, n)) for s in STATES]
        w = max(w, max(v) - min(v))
    return w


print("\n== D1  THE COUNTEREXAMPLE FAMILY: NON-DIAGONAL, NON-CLASS-CONSTANT, AND pi-ONLY ==")
LF_, LC_ = lam(lf, a, NV), lam(lc, a, NV)
worst_root = worst_sp = 0.0
n_nondiag = n_nonclass = 0
for trial in range(200):
    R = make_R(rng, m=3)
    UF, UC = LF_ @ R, LC_ @ R
    worst_root = max(worst_root,
                     np.linalg.norm(np.linalg.matrix_power(UF, LF) - MF),
                     np.linalg.norm(np.linalg.matrix_power(UC, LC) - MC))
    n_nondiag += (not np.allclose(UF, np.diag(np.diag(UF)), atol=1e-9))
    n_nonclass += (not is_class_const_diag(UF))
    worst_sp = max(worst_sp, spread(UF, UC))
print(f"  200 sampled pairs from the family:  max ||U^L - M|| = {worst_root:.2e}")
print(f"     NON-DIAGONAL U_F in {n_nondiag}/200 draws;  NOT class-constant-diagonal in {n_nonclass}/200")
print(f"     WORST pi-spread over 41 pi-identical states and n <= 12:  {worst_sp:.2e}")
R0 = make_R(np.random.default_rng(7), m=3)
UF0, UC0 = LF_ @ R0, LC_ @ R0
print(f"  one explicit witness, U_F (real part rounded):\n{np.round(UF0.real,4)}")
print(f"     ||U_F - diag(U_F)|| = {np.linalg.norm(UF0-np.diag(np.diag(UF0))):.4f}   "
      f"class-constant diagonal? {is_class_const_diag(UF0)}   "
      f"||U_F^3 - M_F|| = {np.linalg.norm(np.linalg.matrix_power(UF0,3)-MF):.2e}   "
      f"pi-spread = {spread(UF0,UC0):.2e}")
print("  -> THE 'ONLY IF' HALF OF THE STEELMAN'S CHARACTERISATION IS FALSE, and the preserving set")
print("     is NOT measure zero -- this family is positive-dimensional (R ranges over a product of")
print("     conjugacy manifolds of order-3 unitaries, one per class block).")

print("\n== D2  THE CORRECT CRITERION IS ON THE PAIR: K_n = (U_F*)^n U_C^n CLASS-CONSTANT DIAGONAL ==")
TF, TC = X.T_edge(lf, a, NV), X.T_edge(lc, a, NV)
DF, DC = X.D_uniform(lf, a, NV), X.D_uniform(lc, a, NV)
for nm, uF, uC in (("COR-F edge tick T", TF, TC), ("uniform root D", DF, DC),
                   ("counterexample Lambda.R", UF0, UC0)):
    kk = [is_class_const_diag(np.linalg.matrix_power(uF.conj().T, n) @ np.linalg.matrix_power(uC, n))
          for n in range(1, 8)]
    print(f"  {nm:<26} K_n class-constant diagonal for n=1..7: {kk}    pi-spread {spread(uF,uC):.2e}")
print("  -> the criterion tracks pi-only-ness exactly, and it is a statement about the PAIR.")
print("     'Both roots diagonal and class-constant' is SUFFICIENT and NOT NECESSARY.")

print("\n== D3  WHY 4000 INDEPENDENT DRAWS FOUND 0 COUNTEREXAMPLES: THE SAMPLER CANNOT PRODUCE ONE ==")
ind_hits = 0
for _ in range(4000):
    UF = X.random_root(lf, a, NV, rng, "generic")
    UC = X.random_root(lc, a, NV, rng, "generic")
    if spread(UF, UC, nmax=4) < 1e-12:
        ind_hits += 1
print(f"  independent sampler, 4000 pairs: pi-only draws = {ind_hits}   (steelman leg G1 reports 4,")
print(f"  all class-constant diagonal).  Reproduced in kind.  But the counterexample family lives on")
print(f"  the DIAGONAL of the pair space (U_F and U_C share R), which an independent sampler visits")
print(f"  with probability zero.  On W-08's own norm, G1 IS A CONTROL THAT COULD NOT HAVE FAILED.")

print("\n== D4  WHAT THIS DOES AND DOES NOT COST THE STEELMAN ==")
print("  COSTS: the characterisation, the 'measure zero' claim, and the sentence 'every")
print("  non-class-constant root breaks it'.  All three are false as stated.")
print("  DOES NOT COST: the trivial-limit selection.  Every member of my family has")
print("  U_F(0) = U_C(0) = R and therefore ALSO passes the trivial-connection contact point --")
print("  it enlarges the ADMITTED set without adding a single incidence-VISIBLE member.")
print("  Checked at a = 0 explicitly:")
a0 = np.zeros(6)
LF0, LC0 = lam(lf, a0, NV), lam(lc, a0, NV)
print(f"     ||Lambda_F(0) - I|| = {np.linalg.norm(LF0-np.eye(NV)):.2e}  so U_F(0) = U_C(0) = R,"
      f"  max_n,s |Z_n| at a=0 over the family:")
for trial in range(3):
    R = make_R(rng, m=3)
    vals = [abs(X.Z(LF0 @ R, LC0 @ R, s, n, n)) for s in STATES[:6] for n in range(1, 7)]
    print(f"     draw {trial+1}:  min |Z_n(a=0)| = {min(vals):.12f}   max = {max(vals):.12f}"
          f"   (no formation, exactly)")
