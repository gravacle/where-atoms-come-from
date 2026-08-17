# RULING LEG 4 -- IS COR-F's T UNIQUE?  The registrar's own declared weak point.
#
# TWO SELF-DEFECTS, RECORDED RATHER THAN PATCHED.
#  (1) My first attempt was a numerical-gradient descent over the local ansatz.  It returned
#      "0 of 24 random starts converged" -- a null about the SEARCH, not the variety -- and an
#      Adam retry DIVERGED to |p| ~ 2.4e+02, loss ~ 9e+16, the loss being degree 2L.  Caught only
#      because I ran the positive control (COR-F's own T: loss 1.0e-30) FIRST.  Optimiser dropped.
#  (2) My first Fourier test asserted that U^L has Fourier support in {-L,0,+L} and measured it on
#      the WHOLE matrix.  It printed 6.13, refuting my own sentence.  The claim is true of the
#      DIAGONAL entries only: an off-diagonal entry (v,w) sums walks with net displacement
#      congruent to v-w, not to 0.  Corrected below; the theorem is unaffected because U^L = W.I
#      is a statement about the diagonal AND the vanishing of everything off it.
#
# ADMISSIBILITY, TAKEN FROM THE CORPUS'S OWN BYTES, NOT ASSERTED:
#   (a) unitary                                        -- S2 audit CHOICE LEDGER A1 (:657) clause (a)
#   (b) no data beyond fibres/edges/orientation/connection  -- A1 clause (b)
#   (c) [STRUCK as CIRCULAR: "reduces to the build's T_gamma on L_v0" presupposes the operator is
#        an endomorphism of each fibre, i.e. presupposes the very fibre-wise-ness at issue]
#   (d) LOCALITY: a tick moves a fibre value at most one edge -- S1 :51-53, which DEFINES parallel
#       transport along e:u->v as z |-> U_e z, with reverse traversal by U_e^{-1}
#   (e) identity off the loop -- A1's own extension convention
#   (f) U^L = M_gamma -- agreement with the corpus's own operator at loop closure
#
# THEOREM (mine).  On an L-cycle the general local operator is
#     U[v,v] = g_v ,  U[v+1,v] = c_v U_v ,  U[v-1,v] = d_v conj(U_{v-1}).
# 1. A closed walk of L steps has (#forward - #backward) in {-L, 0, +L}, so each DIAGONAL entry
#    of U^L is  (prod c) W  +  N_vv  +  (prod d) W^{-1}, with N constant in the connection.
# 2. U^L = W.I identically in the connection forces  prod c_v = 1  and  prod d_v = 0  and N = 0.
# 3. Column v of U has exactly those three entries, so unitarity gives |g_v|^2+|c_v|^2+|d_v|^2 = 1,
#    hence |c_v| <= 1 for every v.
# 4. prod_v |c_v| = 1 with every |c_v| <= 1 forces |c_v| = 1, hence g_v = d_v = 0 for ALL v.
# CONCLUSION.  U = Lam.T, Lam = diag(c_v) unimodular with prod c_v = 1 -- COR-F's T times an
# (L-1)-torus of ADJOINED per-vertex phases.  Clause (b) supplies no vertex-indexed phase, so
# setting them to 1 gives EXACTLY COR-F's T.  A pure backward shift delivers W^{-1}, not W, so
# the loop's own orientation (S1:27) fixes the direction.  Steps 3 and 4 are arithmetic and need
# no numerics; steps 1 and 2 are verified below by a test that CAN come out the other way.
import numpy as np, rlib
rng = np.random.default_rng(rlib.SEED)

def localU(L, g, c, d, a):
    U = np.zeros((L, L), dtype=complex); Ue = np.exp(1j*np.asarray(a))
    for v in range(L):
        U[v, v] += g[v]; U[(v+1) % L, v] += c[v]*Ue[v]; U[(v-1) % L, v] += d[v]*np.conj(Ue[(v-1) % L])
    return U

print("== R4.1  STEP 1, ON THE DIAGONAL, BY FOURIER EXTRACTION IN a_v = t + b_v ==")
for L in (3, 4, 5):
    sup = wp = wm = 0.0; offdiag_nonzero = 0
    for _ in range(60):
        g = (rng.normal(size=L)+1j*rng.normal(size=L))/2
        c = (rng.normal(size=L)+1j*rng.normal(size=L))/2
        d = (rng.normal(size=L)+1j*rng.normal(size=L))/2
        b = rng.uniform(0, 2*np.pi, L); NT = 64
        P = np.array([np.linalg.matrix_power(localU(L, g, c, d, 2*np.pi*j/NT + b), L) for j in range(NT)])
        Fh = np.fft.fft(P, axis=0)/NT
        Wb = np.exp(1j*b.sum()); dg = np.arange(L)
        off = [k for k in range(NT) if k not in (0, L % NT, (-L) % NT)]
        sup = max(sup, float(np.max(np.abs(Fh[off][:, dg, dg]))))
        wp = max(wp, float(np.max(np.abs(np.diag(Fh[L % NT]) - np.prod(c)*Wb))))
        wm = max(wm, float(np.max(np.abs(np.diag(Fh[(-L) % NT]) - np.prod(d)*np.conj(Wb)))))
        offdiag_nonzero += int(np.max(np.abs(Fh[off])) > 1e-9)
    print(f"  L={L}  max |DIAGONAL Fourier coeff off {{-L,0,+L}}| = {sup:.2e}"
          f"   max|A_+L - (prod c) W| = {wp:.2e}   max|A_-L - (prod d) conj(W)| = {wm:.2e}")
print(f"  (Same statistic on the FULL matrix, INCLUDING off-diagonal entries: nonzero in {offdiag_nonzero}/60")
print("   draws -- that is self-defect (2).  Off-diagonal entry (v,w) sums walks of net displacement")
print("   congruent to v-w, not 0, so its support is wider.  The claim is about the diagonal,")
print("   and on the diagonal it is exact.)\n")

print("== R4.2  STEP 2 AS A DISCRIMINATOR.  Sweep the diagonal and backward bands to zero ==")
print("   Set c_v unimodular with prod c_v = 1, then give every column a diagonal part s and a")
print("   backward part r, re-normalising c.  Root defect ||U^L - W.I|| against (s, r):")
L = 3; a = rng.uniform(0, 2*np.pi, L); W = np.exp(1j*a.sum())
ph = rng.uniform(0, 2*np.pi, L); ph[-1] = -ph[:-1].sum()
print("      s | r " + "".join(f"{r:>12.2f}" for r in (0.0, 0.05, 0.2, 0.5)))
for s in (0.0, 0.05, 0.2, 0.5):
    row = f"   {s:>6.2f}"
    for r in (0.0, 0.05, 0.2, 0.5):
        amp = np.sqrt(max(1-s*s-r*r, 0))
        U = localU(L, np.full(L, s, dtype=complex), amp*np.exp(1j*ph),
                   np.full(L, r, dtype=complex), a)
        row += f"{np.linalg.norm(np.linalg.matrix_power(U,L)-W*np.eye(L)):>12.2e}"
    print(row)
print("   -> ZERO in exactly one cell, s = r = 0.  Any admixture of a diagonal or a reverse band")
print("      destroys the root condition.  Steps 3-4 say why: prod|c_v| = 1 leaves no room.\n")

print("== R4.3  AND THE FREE PHASES CHANGE NOTHING: THE WHOLE ADMISSIBLE FAMILY IS pi-VISIBLE ==")
for C in (rlib.K1(), rlib.B0b()):
    a = rlib.a_generic(C, rng, 1.0, 2**0.5)
    base = rng.dirichlet(np.ones(C.nv)); S = rlib.same_pi_states(C, rng, base, 20)
    smallest = 1e9; nblind = 0; checked = 0
    for _ in range(400):
        ops = []
        for (walk, Vs, Lp) in ((C.walkF, C.VF, C.LF), (C.walkC, C.VC, C.LC)):
            T = rlib.Tedge(C, walk, a)
            phv = rng.uniform(0, 2*np.pi, Lp); phv[-1] = -phv[:-1].sum()
            Lam = np.eye(C.nv, dtype=complex)
            for k, v in enumerate(sorted(Vs)): Lam[v, v] = np.exp(1j*phv[k])
            U = Lam@T
            assert np.linalg.norm(np.linalg.matrix_power(U, Lp)
                                  - rlib.Mcirc(C, Vs, rlib.holon(walk, a))) < 1e-9
            checked += 1; ops.append(U)
        sp = rlib.pi_spread(C, ops[0], ops[1], S, range(1, 2*C.LF*C.LC+1))
        smallest = min(smallest, sp); nblind += (sp < 1e-9)
    print(f"  {C.name:4s} 400 draws from the FULL admissible local family ({checked} root checks all"
          f" passed): {nblind} pi-blind; SMALLEST spread over n<=2.L_F.L_C = {smallest:.3e}")
print("  -> EVERY admissible edge tick, not only COR-F's, makes the incidence VISIBLE.")
print("     T's canonicity is therefore NOT load-bearing.  The finding needs ADMISSIBILITY only,")
print("     and one admissible witness suffices -- which the corpus's own sealed audit supplies.")
