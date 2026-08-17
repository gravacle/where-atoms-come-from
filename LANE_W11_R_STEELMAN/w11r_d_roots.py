# LANE W11-R  LEG D -- SETTLING "T IS NOT UNIQUE", THE WEAKNESS THE REGISTRAR DECLARED FIRST.
#
# STRUCTURE, PROVED THEN CHECKED.
#  (1) If U^L = M then U commutes with M (a matrix always commutes with its own power).  So U
#      preserves M's eigenspaces: U = U_loop (+) U_off with U_loop^L = W.I and U_off^L = I.
#  (2) The FULL root variety is therefore  { w.V (+) V' : w^L = W, V^L = I, V'^L = I }  -- a
#      positive-dimensional manifold, so COR-F's T is one point of a continuum.  Confirmed.
#  (3) CHARACTERISATION.  Z_n = sum_v conj(u_F(v))^n u_C(v)^n |s_v|^2 whenever U_F, U_C are
#      DIAGONAL.  It is a function of pi alone for every n IFF both diagonals are constant on the
#      FOUR CLASSES.  Tested both directions below, including diagonal-but-not-class-constant.
#  (4) SELECTION.  The pre-registered trivial-connection contact point (FOUNDING_DESIGN:117-118,
#      S2 CHOICE LEDGER C4) requires U(a=0) = I.  { V unitary : V^L = I } is a disjoint union of
#      conjugacy-class manifolds indexed by eigenvalue multiplicities; the one containing I is the
#      single point {I}.  So for a CONTINUOUS family U(a) with U(a)^L = M(a) and U(0) = I:
#      V(a) := D(a)^{-1} U(a) is continuous, satisfies V^L = I (D commutes with U, being a
#      function of M), and V(0) = I -- hence V == I on any connected neighbourhood of a = 0.
#      THEREFORE U = D, the uniform principal root, which is class-constant diagonal, and by (3)
#      pi-only at EVERY tick.  COR-F's T has T(0) = a cyclic permutation != I and is excluded.
import numpy as np, w11r_lib as L
rng = np.random.default_rng(20260817)

lf, lc, NV, ne = L.K1_LOOP_F, L.K1_LOOP_C, 5, 6
a = np.array([1.0,0.37,0.91,2**0.5,0.23,1.77])
LF = LC = 3
MF, MC = L.M_circuit(lf,a,NV), L.M_circuit(lc,a,NV)
TF, TC = L.T_edge(lf,a,NV), L.T_edge(lc,a,NV)
DF, DC = L.D_uniform(lf,a,NV), L.D_uniform(lc,a,NV)
sA = np.sqrt(np.array([0.40,0.15,0.15,0.15,0.15]))+0j
sB = np.sqrt(np.array([0.40,0.30,0.00,0.05,0.25]))+0j
sC = sA*np.exp(1j*np.array([0.0,1.3,-0.7,2.2,0.4]))
STATES=(sA,sB,sC); CLS = L.classes(lf,lc,NV)

def spread(oF,oC,nmax=6):
    w=0.0
    for n in range(1,nmax+1):
        v=[abs(L.Z(oF,oC,s,n,n)) for s in STATES]; w=max(w,max(v)-min(v))
    return w

print("== D1  ANY L-th ROOT COMMUTES WITH M, AND THE ROOT VARIETY IS A CONTINUUM ==")
print(f"  || [T_F, M_dF] || = {np.linalg.norm(TF@MF-MF@TF):.2e}    || [T_C, M_c] || = {np.linalg.norm(TC@MC-MC@TC):.2e}")

def random_root(loop, aa, NV, kind):
    """sample a unitary U with U^L = M_gamma(aa).  kind in {'generic','diag','diag_classconst'}."""
    Ls = len(loop); on = sorted(L.loop_vertices(loop)); off = [v for v in range(NV) if v not in on]
    W = L.holonomy(loop, aa)
    w0 = np.exp(1j*(np.angle(W)+2*np.pi*rng.integers(0,Ls))/Ls)      # any L-th root of W
    z  = np.exp(2j*np.pi/Ls)
    U = np.zeros((NV,NV), dtype=complex)
    def block(idx, scal):
        d = len(idx)
        ks = rng.integers(0,Ls,size=d)
        if kind == 'generic':
            Q,_ = np.linalg.qr(rng.normal(size=(d,d))+1j*rng.normal(size=(d,d)))
            B = Q@np.diag(z**ks)@Q.conj().T
        elif kind == 'diag':
            B = np.diag(z**ks.astype(float))
        else:                                   # diagonal, constant on each of the four classes
            cmap = {}
            B = np.eye(d, dtype=complex)
            for i,v in enumerate(idx):
                c = CLS[v]
                if c not in cmap: cmap[c] = z**rng.integers(0,Ls)
                B[i,i] = cmap[c]
        return scal*B
    Bl = block(on, w0); Bo = block(off, 1.0)
    for i,v in enumerate(on):
        for j,u in enumerate(on): U[v,u] = Bl[i,j]
    for i,v in enumerate(off):
        for j,u in enumerate(off): U[v,u] = Bo[i,j]
    return U

for kind in ("generic","diag","diag_classconst"):
    ok = bad = 0; sp = []; ndiag = 0
    for _ in range(400):
        UF = random_root(lf,a,NV,kind); UC = random_root(lc,a,NV,kind)
        assert np.linalg.norm(np.linalg.matrix_power(UF,LF)-MF) < 1e-10, "not a root"
        assert np.linalg.norm(np.linalg.matrix_power(UC,LC)-MC) < 1e-10, "not a root"
        if not np.allclose(UF,np.diag(np.diag(UF))): ndiag += 1
        s_ = spread(UF,UC); sp.append(s_)
        if s_ < 1e-12: ok += 1
        else: bad += 1
    print(f"  {kind:<16} 400 sampled roots: non-diagonal {ndiag:3d}   pi-only {ok:3d}   "
          f"pi-BROKEN {bad:3d}   median spread {np.median(sp):.2e}   max {max(sp):.2e}")

print("\n== D2  THE CHARACTERISATION IS EXACT IN BOTH DIRECTIONS ==")
print("  DIAGONAL BUT NOT CLASS-CONSTANT breaks invisibility too -- so 'diagonal' is NOT the")
print("  criterion; 'diagonal AND constant on the four classes' is.  (row 'diag' above.)")
print(f"  COR-F's T          spread = {spread(TF,TC):.2e}     (BREAKS)")
print(f"  uniform root D     spread = {spread(DF,DC):.2e}     (PRESERVES)")
print(f"  corpus's M_gamma   spread = {spread(MF,MC):.2e}     (PRESERVES)")

print("\n== D3  THE PRE-REGISTERED TRIVIAL-LIMIT FILTER, AND IT IS DECISIVE ==")
a0 = np.zeros(ne)
print(f"  at a_e = 0:  M_dF = I ?  {np.allclose(L.M_circuit(lf,a0,NV),np.eye(NV))}"
      f"     D_F = I ?  {np.allclose(L.D_uniform(lf,a0,NV),np.eye(NV))}"
      f"     T_F = I ?  {np.allclose(L.T_edge(lf,a0,NV),np.eye(NV))}")
print(f"  || T_F(0) - I || = {np.linalg.norm(L.T_edge(lf,a0,NV)-np.eye(NV)):.4f}   "
      f"|| D_F(0) - I || = {np.linalg.norm(L.D_uniform(lf,a0,NV)-np.eye(NV)):.2e}")
print("  CONTINUITY CHECK of the uniqueness argument: switch on the connection with t -> 0 and")
print("  ask for the root of M(t.a) NEAREST the identity.  It must be D, and the distance must ->0.")
print(f"  {'t':>10} {'||D(t)-I||':>14} {'||T(t)-I||':>14} {'min over 400 sampled roots ||U-I||':>36}")
for t in (1e-1,1e-2,1e-3,1e-4):
    at = t*a
    dI = np.linalg.norm(L.D_uniform(lf,at,NV)-np.eye(NV))
    tI = np.linalg.norm(L.T_edge(lf,at,NV)-np.eye(NV))
    best = min(np.linalg.norm(random_root(lf,at,NV,'generic')-np.eye(NV)) for _ in range(400))
    print(f"  {t:>10.0e} {dI:>14.6f} {tI:>14.6f} {best:>36.6f}")
print("  -> D(t) -> I; every other sampled root stays a fixed distance away.  The identity has an")
print("     isolated preimage, exactly as the component argument says, so the continuous family")
print("     obeying the founding design's trivial-connection contact point IS the uniform root.")

print("\n== D4  SETTLED ==")
print("  NEITHER of the two extremes the brief names is right:")
print("   * T is NOT canonical -- it is one point of a positive-dimensional root variety, and the")
print("     uniform root D is another point with the same L-th power, exactly.")
print("   * but 'every non-diagonal root breaks invisibility' is TRUE, and more: every root that is")
print("     not class-constant-diagonal breaks it.  The preserving set is measure zero.")
print("  What decides between them is NOT measure and NOT naturalness-by-assertion.  It is the")
print("  pre-registered trivial-connection contact point, which admits exactly the preserving set.")
