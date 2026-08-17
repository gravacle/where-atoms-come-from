# verify_core.py -- LANE W20_R_LEDGER.  PROVE THE MACHINERY BEFORE ANY ARM RUNS.
# Every check prints a NUMBER.  A silent pass is not a pass.
import numpy as np, itertools, sys
from core_w20r import *

LOG = []
def P(*a):
    s = " ".join(str(x) for x in a); LOG.append(s); print(s, flush=True)
def rule(t=""):
    P("\n" + "=" * 100); P(t); P("=" * 100)

rule("V1 -- THE CARRIER AND THE TWO FORCED SUBSPACES")
P("V = %d  L = %d  edges = %s" % (V_N, L, EDGES))
deg = [0]*V_N
for (a,b) in EDGES: deg[a]+=1; deg[b]+=1
P("degrees %s  cubic=%s" % (deg, all(d==3 for d in deg)))
P("dim Z_1 (cycles)   = %d   |Z_1| = %d" % (len(gf2_rref(CYC)), len(CYC)))
P("dim B^1 (cocycles) = %d   |B^1| = %d" % (len(COCYC_BASIS), len(COCYC)))
odd = sum(1 for z in CYC for c in COCYC if pop(z & c) % 2)
P("odd pairings <z,c> over Z_1 x B^1 : %d of %d  (0 => mutual annihilators)" % (odd, len(CYC)*len(COCYC)))

rule("V2 -- THE PHYSICAL SECTOR IS BUILT CORRECTLY:  G_v |psi> = eta_v |psi>, PRINTED")
for nm, eta in [("vacuum", [1]*8), ("charges at {0,4}", [-1,1,1,1,-1,1,1,1])]:
    be = PhysSector(eta)
    P("%-18s dim = %d" % (nm, be.dim))
    worst = 0.0
    for v in range(V_N):
        for i in range(be.dim):
            e = np.zeros(be.dim, dtype=complex); e[i] = 1
            r = be.apply(e, STAR[v], 0) - eta[v]*e
            worst = max(worst, float(np.linalg.norm(r)))
    P("   max || G_v |e_i> - eta_v |e_i> ||  over all v,i = %.12f" % worst)
    assert worst < 1e-12

rule("V3 -- THE EFFECTIVE OPERATOR SPACE IS EXACTLY B(H_phys):  2^10 = 1024 = 32^2")
be = PhysSector([1]*8)
allgens = [mk(1<<l, 0) for l in range(L)] + [mk(0, z) for z in gf2_rref(CYC).values()]
FULLALG = Algebra(be, allgens, N_PHYS, "B(H)")
P("label-space dim n = %d   -> algebra dim = %d   (32^2 = %d)" % (FULLALG.n, FULLALG.dim, 32*32))
P("radical dim r = %d   symplectic pairs s = %d   max entropy = %d bits" % (FULLALG.r, FULLALG.s, FULLALG.maxent))
assert FULLALG.n == 10 and FULLALG.r == 0 and FULLALG.s == 5

rule("V4 -- ENTROPY CALIBRATION.  A pure state must give S(B(H)) = 0 EXACTLY.")
rng = np.random.default_rng(11)
for t in range(3):
    psi = rng.normal(size=32) + 1j*rng.normal(size=32); psi /= np.linalg.norm(psi)
    S = FULLALG.entropy(psi)
    P("   Haar pure state %d :  S(rho | B(H)) = %.12f bits   (must be 0)" % (t, S))
    assert abs(S) < 1e-8
# maximally mixed check via a 2-qubit-like subalgebra with a random state: S <= maxent
rule("V5 -- THE FOUR CHANNELS OF A_S, AND THE PRE-REGISTERED DIMENSIONS REPRODUCED")
S_LINKS = (1<<1)|(1<<2)|(1<<3)
W_S = (1<<1)|(1<<2)|(1<<3)
A_S   = alg_links(be, S_LINKS, N_PHYS, "A_S FULL")
A_CEN = Algebra(be, [mk((1<<1)|(1<<2),0), mk((1<<1)|(1<<3),0)], N_PHYS, "A_S CENTRE")
A_BLK = Algebra(be, [mk(1<<1,0), mk(0,W_S)], N_PHYS, "A_S BLOCK")
A_MAG = Algebra(be, [mk(0,W_S)], N_PHYS, "A_S MAG")
for A, want_dim, want_max in [(A_S,16,3),(A_CEN,4,2),(A_BLK,4,1),(A_MAG,2,1)]:
    P("%-12s n=%d dim=%-4d r=%d s=%d  maxent=%d bits   (pre-registered dim %d, max %d)"
      % (A.label, A.n, A.dim, A.r, A.s, A.maxent, want_dim, want_max))
    assert A.dim == want_dim and A.maxent == want_max

rule("V6 -- Z(A_S) IS THE CENTRE OF A_S, COMPUTED NOT ASSERTED")
cen = [v for v in gf2_span(list(A_S.W.values())) if all(omega(v,w)==0 for w in A_S.W.values())]
P("centre label set size = %d  -> centre dim = %d  (pre-registration says 4)" % (len(cen), len(cen)))
P("centre labels as (a,b) with a reduced mod B^1: %s" % [(bits(sa(v)), bits(sb(v))) for v in cen])
assert len(cen) == 4

rule("V7 -- Z(A_S) *IS* THE SURFACE ALGEBRA.  OPERATOR IDENTITY, NORM PRINTED.")
SIG = (1<<0)|(1<<4)|(1<<5)
for nm, eta in [("vacuum", [1]*8), ("charges {0,4}", [-1,1,1,1,-1,1,1,1]), ("charges {4,5}", [1,1,1,1,-1,-1,1,1])]:
    b2 = PhysSector(eta)
    def M(a,bb):
        out = np.zeros((32,32), dtype=complex)
        for i in range(32):
            e = np.zeros(32, dtype=complex); e[i]=1
            out[:,i] = b2.apply(e,a,bb)
        return out
    n1 = np.linalg.norm(M((1<<1)|(1<<2),0) - eta[0]*M(1<<0,0))
    n2 = np.linalg.norm(M((1<<1)|(1<<3),0) - eta[1]*M(1<<4,0))
    n3 = np.linalg.norm(M((1<<2)|(1<<3),0) - eta[2]*M(1<<5,0))
    n4 = np.linalg.norm(M(SIG,0) - eta[0]*eta[1]*eta[2]*np.eye(32))
    P("%-16s ||X1X2 - eta0 X0||=%.10f  ||X1X3 - eta1 X4||=%.10f  ||X2X3 - eta2 X5||=%.10f  ||X^Sigma - flux*I||=%.10f"
      % (nm, n1, n2, n3, n4))
    assert max(n1,n2,n3,n4) < 1e-10

rule("V8 -- THE VACUITY THEOREM V1, MEASURED:  I(Z(A_S) : A_Sigma) = H(Z(A_S)) FOR EVERY STATE")
A_SIG = alg_links(be, SIG, N_PHYS, "A_Sigma")
P("A_Sigma: n=%d dim=%d r=%d s=%d maxent=%d bits   D(Sigma)=%d" % (A_SIG.n, A_SIG.dim, A_SIG.r, A_SIG.s, A_SIG.maxent, D_of(SIG)))
rng = np.random.default_rng(5)
for t in range(4):
    psi = rng.normal(size=32) + 1j*rng.normal(size=32); psi /= np.linalg.norm(psi)
    I_, _ = MI(psi, A_CEN, A_SIG)
    P("   random state %d :  H(Z(A_S)) = %.9f   I(Z(A_S):A_Sigma) = %.9f   gap = %.2e"
      % (t, A_CEN.entropy(psi), I_, abs(I_ - A_CEN.entropy(psi))))
    assert abs(I_ - A_CEN.entropy(psi)) < 1e-7
P(">>> CONFIRMED AS A THEOREM, NOT A FINDING.  This is FIT-1 and it is now quarantined.")

rule("V9 -- FRAGMENTS: D(F) REPRODUCED FROM THE PRE-REGISTRATION")
def LM(ls):
    m = 0
    for l in ls: m |= 1 << l
    return m
E_ENV = LM([0,4,5,6,7,8,9,10,11])
prim = {"F1":[0,4,5], "F2":[7,8,9], "F3":[11], "F4":[6,10]}
seco = {"G1":[0,4,5], "G2":[6,9], "G3":[7,11], "G4":[8,10]}
for nm, d in [("PRIMARY", prim), ("SECONDARY", seco)]:
    P("%s partition:" % nm)
    u = 0
    for k, ls in d.items():
        m = LM(ls); u |= m
        A = alg_links(be, m, N_PHYS, k)
        P("   %s links %-12s |F|=%d  D(F)=%d  cycles-inside=%d  n=%d dim=%d r=%d s=%d maxent=%d"
          % (k, ls, len(ls), D_of(m), len(gf2_rref(cycles_inside(m))), A.n, A.dim, A.r, A.s, A.maxent))
    P("   union == E_env : %s" % (u == E_ENV))
    assert u == E_ENV
A_ENV = alg_links(be, E_ENV, N_PHYS, "A_env")
P("A_env: |F|=9 D=%d cycles-inside-dim=%d  n=%d dim=%d r=%d s=%d maxent=%d bits"
  % (D_of(E_ENV), len(gf2_rref(cycles_inside(E_ENV))), A_ENV.n, A_ENV.dim, A_ENV.r, A_ENV.s, A_ENV.maxent))

rule("V10 -- JOIN != UNION, MEASURED (the G3 v G4 magnetic gap the pre-registration predicted)")
G3 = alg_links(be, LM([7,11]), N_PHYS, "G3"); G4 = alg_links(be, LM([8,10]), N_PHYS, "G4")
J = join(G3, G4); U = alg_links(be, LM([7,8,10,11]), N_PHYS, "A_{7,8,10,11}")
P("dim(G3 v G4) = %d (n=%d)    dim A_{7,8,10,11} = %d (n=%d)   gap = %d bit(s) of label space"
  % (J.dim, J.n, U.dim, U.n, U.n - J.n))
P("the extra generator is the cycle {7,8,10,11}: is it inside the JOIN label space? %s"
  % (gf2_red(mk(0, LM([7,8,10,11])), J.W) == 0))

rule("V11 -- MONOTONICITY GUARD: A_Sigma is a SUBALGEBRA of A_env (so Delta_surf >= 0 is forced)")
inside = all(gf2_red(v, A_ENV.W) == 0 for v in A_SIG.W.values())
P("every A_Sigma label lies in A_env's label space : %s" % inside)
assert inside

rule("V12 -- LANCZOS AGREES WITH EXACT DIAGONALISATION ON THE 32-DIM SECTOR")
for g2 in [0.20, 0.645497, 3.00]:
    H = H_matrix(be, g2)
    ev = np.linalg.eigvalsh(H)
    E0, psi = lanczos_ground(be, g2, m=30, seed=3, restarts=3)
    P("   g2=%.6f   eigh E0 = %.12f   lanczos E0 = %.12f   diff = %.2e   gap E1-E0 = %.9f"
      % (g2, ev[0], E0, abs(ev[0]-E0), ev[1]-ev[0]))
    assert abs(ev[0]-E0) < 1e-8

rule("V13 -- THE FREE (UNGAUGED) BACKEND: 4096 DIMS, LANCZOS CONVERGENCE PRINTED")
fb = FreeSpace()
for g2 in [0.20, 3.00]:
    E0, psi = lanczos_ground(fb, g2, m=120, seed=3, restarts=4)
    res = float(np.linalg.norm(H_matvec(fb, g2, psi) - E0*psi))
    P("   g2=%.4f   E0 = %.12f   || (H-E0)psi || = %.3e   norm = %.12f" % (g2, E0, res, float(np.linalg.norm(psi))))
    assert res < 1e-7

rule("V14 -- W-19's PERRON-FROBENIUS CLAIM, RE-MEASURED ON THIS CARRIER")
P("|<free ground state | projected into vacuum sector>| -- the overlap W-19 reported as 1.000000000")
for g2 in [0.20, 0.645497, 1.00, 3.00]:
    E0f, psif = lanczos_ground(fb, g2, m=150, seed=3, restarts=4)
    # project onto the vacuum physical sector and measure the retained norm
    # |psi_c> coefficient = sqrt(128) * <c|psif> summed with signs; vacuum eta => all signs +1
    coeff = np.zeros(32, dtype=complex)
    for c in range(1 << L):
        coeff[be.orb[c]] += psif[c] * be.sgn[c]
    coeff = coeff / np.sqrt(128.0)
    nrm = float(np.linalg.norm(coeff))
    Hs = H_matrix(be, g2); evs = np.linalg.eigvalsh(Hs)
    P("   g2=%.6f  || P_vac psi_free || = %.9f     E0(free) = %.9f   E0(sector) = %.9f  equal=%s"
      % (g2, nrm, E0f, evs[0], abs(E0f-evs[0]) < 1e-7))

rule("V15 -- SEED-STABILITY OF THE ENTROPY MACHINERY (same state, permuted generator order)")
rng = np.random.default_rng(99)
psi = rng.normal(size=32) + 1j*rng.normal(size=32); psi /= np.linalg.norm(psi)
base = A_S.entropy(psi)
import random
vals = []
for t in range(5):
    g = [mk(1<<1,0), mk(1<<2,0), mk(1<<3,0), mk(0,W_S)]
    random.Random(t).shuffle(g)
    vals.append(Algebra(be, g, N_PHYS, "perm").entropy(psi))
P("S(rho|A_S) under 5 generator orderings: %s" % ["%.12f" % v for v in vals])
P("max spread = %.2e   (the block decomposition must not depend on generator order)" % (max(vals)-min(vals)))
assert max(vals)-min(vals) < 1e-9

P("\nALL CORE CHECKS PASSED.")
open("OUT_verify_core.txt","w").write("\n".join(LOG) + "\n")
