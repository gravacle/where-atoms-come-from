# W10-D REFUTER 1 -- LEG R6.  TWO MORE THINGS INSIDE ROWS MARKED CARRIER_INDEPENDENT.
#
# R6-A  LEG 4C's BRANCH OPERATORS ARE NOT THE BRANCH OPERATORS.  D-26 exhibits W-05's LEG TWO
#       ("the record slot is already inside the carrier, dim_C = 4, overlap agreeing with Z").
#       W-01/S3/W-05 define  (M_gamma s)(v) = W(gamma) s(v) for v on the loop -- multiplication
#       by the LOOP HOLONOMY.  Lane D's leg 4C multiplies by the vertex's full CHARACTER
#       chi = u^a v^b instead:
#           a = amp*np.array([ch[i] if EXP[CLS[i]][0] else 1 for i in range(4)])
#       For class 11 that is e^{i(c-f)}, not W_F = e^{if}.  The consequence is visible IN LANE
#       D'S OWN PRINTED TABLE and is not commented on: the |<a,b>| column and the Z_1 column
#       disagree (0.749367 vs 0.628335 on B0b), where W-05's registered claim is that they AGREE
#       to 6.2e-17.  Implemented correctly they agree exactly.  THE ROW'S CONCLUSION SURVIVES --
#       rank 2 either way -- but its exhibit does not exhibit the object it names.
#
# R6-B  LEG 6C's Haar-projector statistic, replaced by an EXACT one.  Lane D had to correct this
#       block once already (its D-16 self-confound).  dim A^G is a two-line linear-algebra fact
#       and needs no Monte Carlo at all; here it is settled exactly, at the same (V,r).
import numpy as np
from itertools import product

CLS = ('00', '10', '01', '11')
EXP = {'00': (0, 0), '10': (1, 0), '01': (0, 1), '11': (1, 1)}

print("="*100)
print("== R6-A  THE BRANCH OPERATORS IN LEG 4C ==")
print("="*100)
print("  CONVENTIONS (lane D's own PUBLISHED_CONVENTIONS): u = conj(W_F) = e^{-if}, v = W_C = e^{ic},")
print("  chi_ab = u^a v^b.  M_dF multiplies s_v by W_F = e^{+if} for v in gamma_F; M_c by W_C.")
print("  Then <M_dF s, M_c s> = sum_v |s_v|^2 conj(W_F)^{[F]} W_C^{[C]} = sum p_ab u^a v^b = Z_1.")
print("  IT IS AN IDENTITY.  Any implementation for which |<a,b>| != |Z_1| has the wrong operator.")
f0, c0 = 1.3, 2.0
CARR = [("B1  K1  (3cl)", np.array([0, 2/5, 2/5, 1/5])),
        ("B1q spec (3cl)", np.array([1/7, 3/7, 3/7, 0])),
        ("B1p brdg (2cl)", np.array([0, 1/2, 1/2, 0])),
        ("B0b torus(4cl)", np.array([4/9, 2/9, 1/9, 2/9])),
        ("B4  spin (4cl)", np.array([1/6, 1/6, 1/6, 3/6]))]
ch = np.array([np.exp(1j*(-EXP[cl][0]*f0 + EXP[cl][1]*c0)) for cl in CLS])
WF, WC = np.exp(1j*f0), np.exp(1j*c0)
print(f"\n  {'carrier':15s} {'|Z_1|':>10s} {'CORRECT |<a,b>|':>16s} {'dev':>9s} "
      f"{'LANE D |<a,b>|':>15s} {'dev from Z_1':>13s} {'rank both':>10s}")
for nm, p in CARR:
    amp = np.sqrt(p)
    Z1 = sum(p[i]*ch[i] for i in range(4))
    a_ok = amp*np.array([WF if EXP[CLS[i]][0] else 1 for i in range(4)])
    b_ok = amp*np.array([WC if EXP[CLS[i]][1] else 1 for i in range(4)])
    a_D = amp*np.array([ch[i] if EXP[CLS[i]][0] else 1 for i in range(4)])
    b_D = amp*np.array([ch[i] if EXP[CLS[i]][1] else 1 for i in range(4)])
    ov_ok, ov_D = abs(np.vdot(a_ok, b_ok)), abs(np.vdot(a_D, b_D))
    r_ok = np.linalg.matrix_rank(np.stack([a_ok, b_ok]), tol=1e-12)
    r_D = np.linalg.matrix_rank(np.stack([a_D, b_D]), tol=1e-12)
    print(f"  {nm:15s} {abs(Z1):10.6f} {ov_ok:16.6f} {abs(ov_ok-abs(Z1)):9.1e} "
          f"{ov_D:15.6f} {abs(ov_D-abs(Z1)):13.6f} {str((int(r_ok), int(r_D))):>10s}")
print("\n  The CORRECT operator reproduces W-05's registered 'overlap agreeing with Z' to 1e-16 on")
print("  every carrier including both four-class ones.  Lane D's differs from its own adjacent")
print("  Z_1 column by up to 0.72 and the page does not remark on it.")
print("  dim span = 2 UNDER BOTH implementations, so D-26's CONCLUSION (dim_C = 4, carrier-")
print("  independent, conditioned on G != {1}) SURVIVES -- but on the corrected exhibit, which")
print("  additionally reproduces the registered overlap the lane's own version contradicts.")

print("\n"+"="*100)
print("== R6-B  dim A^G FOR THE FIBRE-WISE GROUP, EXACTLY -- NO MONTE CARLO ==")
print("="*100)
print("  A = M_{Vr}(C) acting on (+)_v C^r;  G = U(r)^V acting by conjugation.")
print("  Block condition: U_v X_vw U_w^* = X_vw for all U.")
print("    v = w : U X U^* = X for all U in U(r)  =>  X = scalar          (dim 1 per vertex)")
print("    v != w: take U_v = I, U_w = e^{i th} I  =>  e^{-i th} X = X    =>  X = 0")
print("  so A^G = (+)_v C.I_r and dim_C A^G = V, for every V and every r.  Settled exactly below")
print("  by solving the fixed-point equations for a FINITE generating set that already forces the")
print("  answer (per-vertex diagonal 4th-root phases + a per-vertex 2-cycle), via exact ranks.")
print(f"  {'V':>3s} {'r':>2s} {'dim A':>7s} {'dim A^G (exact solve)':>22s} {'= V?':>6s} "
      f"{'smallest sing. val. gap':>24s}")
for V, r in [(5, 1), (5, 2), (5, 3), (9, 1), (9, 2), (9, 3), (6, 2), (6, 3), (4, 4)]:
    d = V*r
    gens = []
    for v in range(V):
        for gk in range(3):
            U = np.eye(d, dtype=complex)
            if gk == 0:                                   # scalar phase i on block v
                blk = 1j*np.eye(r, dtype=complex)         #   kills all off-diagonal blocks
            elif gk == 1:                                 # diagonal 4th-root phases
                blk = np.diag([1j**(j % 4) for j in range(r)]).astype(complex)
            else:                                         # cyclic permutation of the fibre
                blk = np.zeros((r, r), dtype=complex)
                for j in range(r):
                    blk[j, (j+1) % r] = 1.0
            U[v*r:(v+1)*r, v*r:(v+1)*r] = blk
            gens.append(U)
    # X in A^G  <=>  U X - X U = 0 for every generator ; stack the commutator maps
    rows = []
    for U in gens:
        rows.append(np.kron(U, np.eye(d)) - np.kron(np.eye(d), U.T))
    M = np.vstack(rows)
    sv = np.linalg.svd(M, compute_uv=False)
    tol = max(M.shape)*np.finfo(float).eps*sv[0]
    nz = int((sv > tol).sum())
    dimfix = d*d - nz
    gap = sv[nz-1]/sv[nz] if nz < len(sv) and sv[nz] > 0 else float('inf')
    print(f"  {V:3d} {r:2d} {d*d:7d} {dimfix:22d} {str(dimfix == V):>6s} {gap:24.3e}")
print("  dim A^G = V AT EVERY (V, r) TESTED, INCLUDING TWO LANE D DID NOT RUN (6,3) AND (4,4),")
print("  with a singular-value gap of 1e+15 or better rather than the 0.01 gap of a 20000-draw")
print("  Monte-Carlo projector.  D-16 and D-18 SURVIVE, on a statistic that cannot be tuned by a")
print("  tolerance -- which is what went wrong in lane D's first run of this block.")
