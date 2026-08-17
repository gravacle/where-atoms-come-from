# RULING LEG 3 -- T's CANONICITY.  The registrar's own declared weak point, and the brief's
# stated decision rule.  I settle it on the corpus's OWN admissibility clause, S2 audit
# CHOICE LEDGER A1 (:657): "(a) unitary, (b) uses no data beyond fibres/edges/orientation/
# connection, (c) reduces to the build's T_gamma on L_v0."  Clause (c) is CIRCULAR (it
# presupposes the operator is an endomorphism of each fibre) and I strike it; (a) and (b) stand.
import numpy as np, rlib, itertools
rng = np.random.default_rng(rlib.SEED)
def haar(n, rng):
    z=(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))/np.sqrt(2)
    q,r=np.linalg.qr(z); return q*(np.diag(r)/np.abs(np.diag(r)))

print("== R3.1  THE STRUCTURE OF THE ROOT VARIETY.  U^L = M on the loop  <=>  U = w.V, V^L = I ==")
print("   Hence {n : U^n is diagonal} is closed under addition and contains L: it is d.Z, d | L.")
for C in (rlib.K1(), rlib.B0b()):
    for (walk, Vs, L, tag) in ((C.walkF,C.VF,C.LF,'F'), (C.walkC,C.VC,C.LC,'C')):
        obs, worst = {}, 0.0
        vs = sorted(Vs)
        for _ in range(1200):
            a = rng.uniform(0,2*np.pi,len(C.edges)); W = rlib.holon(walk,a)
            G = haar(L,rng); ex = rng.integers(0,L,L)
            V = G@np.diag(np.exp(2j*np.pi*ex/L))@G.conj().T
            w = np.exp(1j*(np.angle(W)+2*np.pi*rng.integers(0,L))/L)
            U = np.eye(C.nv,dtype=complex); U[np.ix_(vs,vs)] = w*V
            worst = max(worst, np.linalg.norm(np.linalg.matrix_power(U,L)-rlib.Mcirc(C,Vs,W)))
            ds = [n for n in range(1,L+1) if np.linalg.norm(np.linalg.matrix_power(U,n)
                   - np.diag(np.diag(np.linalg.matrix_power(U,n)))) < 1e-9]
            d = min(ds) if ds else L
            obs[d] = obs.get(d,0)+1
        print(f"  {C.name:4s} gamma_{tag} (L={L}):  max||U^L - M|| = {worst:.2e}   observed d: "
              f"{dict(sorted(obs.items()))}   every d divides L: {all(L%d==0 for d in obs)}")
print()

print("== R3.2  THE TWO CRITERIA THAT ACTUALLY DISCRIMINATE, AND ONE THAT DOES NOT ==")
L = 3
print("  (i) CONTINUITY IN THE CONNECTION.  Degree argument: a continuous h:U(1)->U(1) with")
print("      h(W)^L = W would give L.deg(h) = 1, impossible for L>1.  So EVERY branch rule has a cut.")
C = rlib.K1(); vs = sorted(C.VC)
prevT = prevD = None; stepT = stepD = 0.0
a0 = rng.uniform(0,2*np.pi,6)
for t in np.linspace(0, 2*np.pi, 4001):
    a = a0.copy(); a[3] = a0[3] + t                    # wind W_C once around U(1)
    W = rlib.holon(C.walkC,a)
    T = rlib.Tedge(C, C.walkC, a); D = rlib.Droot(C, C.VC, W, 3)
    if prevT is not None:
        stepT = max(stepT, np.linalg.norm(T-prevT)); stepD = max(stepD, np.linalg.norm(D-prevD))
    prevT, prevD = T, D
print(f"      max step over a 4001-point sweep of one full winding:  T = {stepT:.3e}   D = {stepD:.3e}")
print(f"      -> COR-F's T is CONTINUOUS in the connection.  The fibre-wise root D is NOT.")
print("  (ii) CARRIER DATA ONLY (A1 clause (b)).  Every entry must be built from the U_e.")
a1 = rng.uniform(0,2*np.pi,6)
tgt = np.angle(np.exp(1j*(a1[3]+a1[4]+a1[5])/3))
hits = [m for m in itertools.product(range(-4,5),repeat=3)
        if abs(np.angle(np.exp(1j*(m[0]*a1[3]+m[1]*a1[4]+m[2]*a1[5])))-tgt) < 1e-9]
print(f"      integer m in [-4,4]^3 with prod U_e^m_e = W^(1/3):  {hits}   (exhaustive, 729 candidates)")
print(f"      COR-F's T's entries: T[v,u] = U_e -- a SINGLE edge transport, m = one-hot.")
print("      -> D needs a LIFT of the holonomy to R.  S1 sec4 publishes W itself as the complete")
print("         gauge-invariant content; it publishes no lift.  D FAILS A1(b).  T PASSES.")
print("  (iii) GAUGE COVARIANCE (COR-J's unledgered premise) DOES *NOT* DISCRIMINATE:")
wT=wD=wM=0.0
for _ in range(600):
    a=rng.uniform(0,2*np.pi,6); th=rng.uniform(0,2*np.pi,5)
    ag=np.array([a[j]+th[t]-th[s] for j,(s,t) in enumerate(C.edges)])
    s=rng.normal(size=5)+1j*rng.normal(size=5); g=np.exp(1j*th)
    W,Wg=rlib.holon(C.walkC,a),rlib.holon(C.walkC,ag)
    wT=max(wT,np.linalg.norm(rlib.Tedge(C,C.walkC,ag)@(g*s)-g*(rlib.Tedge(C,C.walkC,a)@s)))
    wD=max(wD,np.linalg.norm(rlib.Droot(C,C.VC,Wg,3)@(g*s)-g*(rlib.Droot(C,C.VC,W,3)@s)))
    wM=max(wM,np.linalg.norm(rlib.Mcirc(C,C.VC,Wg)@(g*s)-g*(rlib.Mcirc(C,C.VC,W)@s)))
print(f"      covariance defect over 600 gauge transforms:  T = {wT:.2e}   D = {wD:.2e}   M = {wM:.2e}")
print("      All three pass.  COR-J excludes none of them.\n")

print("== R3.3  T IS *NOT* UNIQUE.  Three admissible-looking rivals, all with U^L = M_gamma ==")
for C in (rlib.K1(), rlib.B0b()):
    a = rlib.a_generic(C, rng, 1.0, 2**0.5)
    base = rng.dirichlet(np.ones(C.nv)); S = rlib.same_pi_states(C, rng, base, 24)
    WF,WC = rlib.holon(C.walkF,a), rlib.holon(C.walkC,a)
    MF,MC = rlib.Mcirc(C,C.VF,WF), rlib.Mcirc(C,C.VC,WC)
    TF,TC = rlib.Tedge(C,C.walkF,a), rlib.Tedge(C,C.walkC,a)
    DF,DC = rlib.Droot(C,C.VF,WF,C.LF), rlib.Droot(C,C.VC,WC,C.LC)
    # rival 1: the (L-1)-torus  A = Lam.T  with Lam diagonal on the loop, loop product 1
    def torus(walk, Vs, L, T):
        ph = rng.uniform(0,2*np.pi,L); ph[-1] = -ph[:-1].sum()
        Lam = np.eye(C.nv,dtype=complex)
        for k,v in enumerate(sorted(Vs)): Lam[v,v] = np.exp(1j*ph[k])
        return Lam@T
    AF,AC = torus(C.walkF,C.VF,C.LF,TF), torus(C.walkC,C.VC,C.LC,TC)
    # rival 2: a Haar-rotated root  U = w G diag(zeta^j) G* on the loop
    def haarroot(Vs, W, L):
        G=haar(L,rng); ex=rng.permutation(L)
        V=G@np.diag(np.exp(2j*np.pi*ex/L))@G.conj().T
        U=np.eye(C.nv,dtype=complex); U[np.ix_(sorted(Vs),sorted(Vs))]=np.exp(1j*np.angle(W)/L)*V
        return U
    HF,HC = haarroot(C.VF,WF,C.LF), haarroot(C.VC,WC,C.LC)
    rows = [("M_gamma  (the corpus's convention)", MF, MC),
            ("T        (COR-F's edge tick)      ", TF, TC),
            ("D        (fibre-wise L-th root)   ", DF, DC),
            ("Lam.T    ((L-1)-torus rival)      ", AF, AC),
            ("Haar root (generic member)        ", HF, HC)]
    print(f"  --- {C.name} ---")
    for nm,UF,UC in rows:
        rL = np.linalg.norm(np.linalg.matrix_power(UF,C.LF)-MF)
        rLs = "n/a (M is the operator, not a root of itself)" if nm.startswith("M_gamma") else f"{rL:.1e}"
        fw = np.allclose(UF, np.diag(np.diag(UF)))
        okQ,_ = rlib.is_class_constant_diag(C, rlib.Qrel(UF,UC,1))
        sp = rlib.pi_spread(C,UF,UC,S,range(1,13))
        print(f"    {nm}  ||U^L-M||={rLs:>4s}  diag={str(fw):5s}  Q_1 class-const={str(okQ):5s}  pi-spread(n<=12)={sp:.2e}")
print()
print("== R3.4  THE BRIEF'S DISJUNCTION IS NOT EXHAUSTIVE, AND BOTH HORNS ARE FALSE AS WRITTEN ==")
print("  HORN 1 'a different EQUALLY NATURAL edge tick restores invisibility, so Reading B falls':")
print("     D restores it -- but D moves NO fibre value along any edge, is DISCONTINUOUS in the")
print("     connection, and fails A1(b).  It is not an edge tick.  It is the corpus's own")
print("     fibre-wise class-constant convention at a finer clock: diag(D)_v = w^(1[v in gamma]).")
for C in (rlib.K1(), rlib.B0b()):
    a=rlib.a_generic(C,rng,1.0,2**0.5); W=rlib.holon(C.walkC,a)
    D=rlib.Droot(C,C.VC,W,C.LC); T=rlib.Tedge(C,C.walkC,a)
    e=np.zeros(C.nv,dtype=complex); e[sorted(C.VC)[0]]=1.0
    print(f"     {C.name:4s} |<e0, T e0>| = {abs(np.vdot(e,T@e)):.3f}  (T MOVES the excitation)   "
          f"|<e0, D e0>| = {abs(np.vdot(e,D@e)):.3f}  (D moves NOTHING)")
print("  HORN 2 'every unitary with U^L = M except the DIAGONAL ones breaks invisibility':")
print("     FALSE.  Correlated non-diagonal pairs are pi-blind (leg R2.2 mode 2).  The census")
print("     below shows why an INDEPENDENT sampler can never find them -- a control that")
print("     could not have failed, in either direction.")
for C in (rlib.K1(), rlib.B0b()):
    a=rlib.a_generic(C,rng,1.0,2**0.5)
    WF,WC=rlib.holon(C.walkF,a),rlib.holon(C.walkC,a)
    def hroot(Vs,W,L):
        G=haar(L,rng); ex=rng.integers(0,L,L)
        V=G@np.diag(np.exp(2j*np.pi*ex/L))@G.conj().T
        U=np.eye(C.nv,dtype=complex); U[np.ix_(sorted(Vs),sorted(Vs))]=np.exp(1j*np.angle(W)/L)*V
        return U
    ind=sum(rlib.is_class_constant_diag(C, rlib.Qrel(hroot(C.VF,WF,C.LF),hroot(C.VC,WC,C.LC),1))[0]
            for _ in range(2000))
    print(f"     {C.name:4s} INDEPENDENT sampling of the root variety, 2000 pairs: {ind} pi-blind at n=1")
print("     0 or near-0 is a fact about the SAMPLER, not about the variety: the pi-blind set is")
print("     the CORRELATED locus and has measure zero for an independent sampler AT ANY SIZE.")
