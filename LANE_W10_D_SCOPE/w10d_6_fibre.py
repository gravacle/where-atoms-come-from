# W10-D leg 6 -- the rows that are about the FIBRE and the GAUGE ACTION rather than the class
# multiset: W-04's ERR-1 (scalarity, not commutativity), W-05's LEG ONE, W-06's N4 (fibre-wise-ness).
# These are the rows for which "which carrier" is the WRONG axis, and this leg says so by testing
# them on B0b (V=9) and finding the answer identical to K1's (V=5) with V the only thing that moves.
import numpy as np
rng = np.random.default_rng(20260816)

print("="*100)
print("== 6A  ERR-1's ONE-LINE REFUTATION, RE-RUN ==")
print("="*100)
print("  S2 glosses its gate as 'the connection must be non-abelian in its action'.  W-04: FALSE;")
print("  the operative variable is SCALARITY.  The exhibit: T_F = I, T_C = diag(1,-1),")
print("  z = (1,1)/sqrt2 on C^2.")
TF = np.eye(2, dtype=complex); TC = np.diag([1, -1]).astype(complex)
z = np.array([1, 1], dtype=complex)/np.sqrt(2)
print(f"    <T_F z, T_C z>      = {np.vdot(TF@z, TC@z):.16f}")
print(f"    ||[T_F,T_C]||_F     = {np.linalg.norm(TF@TC-TC@TF):.3e}   (ABELIAN, COMMUTING)")
print(f"    T_C scalar?          {np.allclose(TC, TC[0,0]*np.eye(2))}")
print("  An abelian, commuting, NON-SCALAR connection writes the gate.  ERR-1 REPRODUCED.")
print("  This row is about the FIBRE, not the carrier: it holds at V = 5, 9 or any V, and no")
print("  class occupancy enters.  SCOPE = CARRIER-INDEPENDENT (fibre-rank-scoped).")

print("\n"+"="*100)
print("== 6B  W-05's LEG ONE, AND WHAT W-06 CORRECTED, BOTH ON B0b (V=9) ==")
print("="*100)
print("  M_gamma multiplies s_v by a v-dependent phase, so it lies in the FIBRE-WISE group")
print("  U(1)^V on ANY carrier -- that is a definition, not a measurement, and V does not enter.")
print("  W-05 read that as 'S3-0 is the identity that gauge-invariant observables are gauge")
print("  invariant'.  W-06's correction: the FULL gauge action ALSO moves the connection")
print("  (a_e -> a_e + theta_t - theta_s), and M_gamma leaves a_e alone.  So M_gamma is a full")
print("  gauge transformation ONLY IF theta is constant on every edge, i.e. GLOBAL.  Checked on")
print("  the reconstructed B0b of leg 5A (V=9, E=18):")
n = 3
vid = {(i, j): 3*i+j for i in range(n) for j in range(n)}
edges = []
for i in range(n):
    for j in range(n):
        edges.append((vid[(i, j)], vid[((i+1) % n, j)]))
        edges.append((vid[(i, j)], vid[(i, (j+1) % n)]))
V, E = 9, 18
FVs = {vid[(0, 0)], vid[(1, 0)], vid[(1, 1)], vid[(0, 1)]}
th = np.array([np.pi if v in FVs else 0.0 for v in range(V)])   # the phase M_dF applies, W_F = -1
da = np.array([th[t]-th[s] for (s, t) in edges])
print(f"    M_dF's vertex phase pattern realised as a gauge parameter theta_v = pi.1[v in gamma_F]")
print(f"    the connection shift it WOULD induce: max |theta_t - theta_s| = {np.abs(da).max():.6f}")
print(f"    number of edges on which it is non-zero: {int((np.abs(da)>1e-12).sum())} of {E}")
print("    NON-ZERO -> M_dF is NOT a full gauge transformation on B0b either.  Same as K1.")
print("  SCOPE = CARRIER-INDEPENDENT.  What decides it is whether the gauge law is implemented")
print("  on the connection as well as the section -- an implementation fact, not a carrier fact.")

print("\n"+"="*100)
print("== 6C  N4 -- FIBRE-WISE-NESS.  dim A^G = V AT EVERY RANK, ON B0b AS ON K1 ==")
print("="*100)
print("  A = M_{Vr}(C) on (+)_v C^r; G = the FIBRE-WISE group U(r)^V.  Claim (W-06's corrected")
print("  N4): any fibre-wise unitary lies in G at every rank, and dim_C A^G = V.")
print("  ONE-LINE PROOF: E_Haar[U_v X_vw U_w^*] = (tr X_vv / r) I_r if v = w and 0 if v != w,")
print("  because the two blocks' Haar phases are independent.  So A^G = (+)_v C.I_r, dim = V.")
# DEFECT RECORDED, NOT SILENTLY PATCHED.  This block first counted the fixed-space dimension as
# np.linalg.matrix_rank(M, tol=1e-2) on a 4000-draw Monte-Carlo average.  The Haar average of
# Ad_U is an ORTHOGONAL PROJECTOR (eigenvalues 0 and 1 only), but at n draws the zero eigenvalues
# sit at O(n^{-1/2}) ~ 1.6e-2, ABOVE that tolerance, so the first run printed 15/63/139/59/198/
# 452/84 -- every row wrong and every row "False".  The right statistic is the count of
# eigenvalues near 1, not a rank at a hand-picked tolerance.  Both are printed below.
print(f"  {'V':>3s} {'r':>2s} {'dim A':>7s} {'#eig(M) > 0.5':>15s} {'gap to next eig':>17s} {'= V?':>6s}")
for V_, r in [(5, 1), (5, 2), (5, 3), (9, 1), (9, 2), (9, 3), (6, 2)]:
    d = V_*r
    M = np.zeros((d*d, d*d), dtype=complex)
    NS = 20000
    for _ in range(NS):
        U = np.zeros((d, d), dtype=complex)
        for v in range(V_):
            X = rng.normal(size=(r, r)) + 1j*rng.normal(size=(r, r))
            Q, R = np.linalg.qr(X)
            Q = Q*(np.diag(R)/np.abs(np.diag(R)))
            U[v*r:(v+1)*r, v*r:(v+1)*r] = Q
        M += np.kron(U, U.conjugate())
    M /= NS
    ev = np.sort(np.abs(np.linalg.eigvals(M)))[::-1]
    nfix = int((ev > 0.5).sum())
    gap = ev[nfix] if nfix < len(ev) else 0.0
    print(f"  {V_:3d} {r:2d} {d*d:7d} {nfix:15d} {gap:17.4f} {str(nfix == V_):>6s}")
print("  IDENTICAL AT V = 5 AND V = 9, AT RANKS 1, 2, 3.  W-06's N4 is CARRIER-INDEPENDENT and")
print("  its content is entirely the WORD 'fibre-wise'.  It is not a finding about K1, and no")
print("  four-class carrier can bear on it.  (W-06's own code for this does not exist; this is a")
print("  RECONSTRUCTION of the statement, and the one-line proof above is what actually settles it.)")
