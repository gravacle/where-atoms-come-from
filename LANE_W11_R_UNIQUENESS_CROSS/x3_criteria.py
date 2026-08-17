# LANE W-11-R-T-CROSS, leg X3 — AUDIT THE CANONICITY ARGUMENT.  My lens is STEELMAN COR-F, so I
# discount the uniqueness lane's KILLS as hard as its confirmations, and I check whether its
# rivals smuggle in structure the carrier does not supply.  Both directions produce findings.
import numpy as np, itertools
rng = np.random.default_rng(11072026)
NV=5; EDG=[(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)]
LF=[(0,1,0),(1,2,1),(2,0,2)]; LC=[(0,3,3),(3,4,4),(4,0,5)]; VF=[0,1,2]; VC=[0,3,4]
def T_edge(loop,a):
    U=np.exp(1j*a); on={s for s,_,_ in loop}; T=np.zeros((NV,NV),dtype=complex)
    for v in range(NV):
        if v not in on: T[v,v]=1.0
    for (s_,d_,e) in loop: T[d_,s_]=U[e]
    return T
def M_of(vs,W):
    M=np.eye(NV,dtype=complex)
    for v in vs: M[v,v]=W
    return M
def holF(a): return np.exp(1j*(a[0]+a[1]+a[2]))
def root_op(W,vs,L,branch=0):
    w=np.exp(1j*(np.angle(W)+2*np.pi*branch)/L); X=np.eye(NV,dtype=complex)
    for v in vs: X[v,v]=w
    return X,w

a = rng.uniform(0,2*np.pi,6)
WF = holF(a); T = T_edge(LF,a); M = M_of(VF,WF); R,w0 = root_op(WF,VF,3)

print("== X3.1  CRITERION (b) AS THE UNIQUENESS LANE ACTUALLY TESTED IT IS VOID ==")
print("   Its leg T3.1 declares (b) 'DECISIVE' and its returned verdict says: \"COR-F's T is")
print("   unchanged (4.97e-16) while the rival's L-th root multiplies by a primitive root of")
print("   unity.\"  ITS OWN SEALED OUTPUT ON THE SAME LINE READS  |w(a')-w(a)| = 0.000000  and")
print("   w(a')/w(a) = 1.000000-0.000000j.  The labels 'root: CHANGED' and '(= a cube root of")
print("   unity)' are hard-coded strings sitting next to numbers that refute them.  Reproduced:")
worstT=worstR=0.0
for _ in range(2000):
    aa=rng.uniform(0,2*np.pi,6); j=rng.integers(0,6); k=rng.integers(-3,4)
    ab=aa.copy(); ab[j]+=2*np.pi*k
    worstT=max(worstT,np.linalg.norm(T_edge(LF,ab)-T_edge(LF,aa)))
    _,w1=root_op(holF(aa),VF,3); _,w2=root_op(holF(ab),VF,3)
    worstR=max(worstR,abs(w2-w1))
print(f"      2000 draws, a_e -> a_e + 2pi.k for a random edge and k in [-3,3]:")
print(f"         max || T(a') - T(a) ||   = {worstT:.2e}      COR-F: unchanged")
print(f"         max |  w(a') - w(a) |    = {worstR:.2e}      ROOT: ALSO UNCHANGED")
print("   REASON: the lane's own root_op takes the principal root OF THE HOLONOMY W, i.e.")
print("   exp(i.angle(W)/L).  W is a point of U(1) and is unchanged by a 2pi shift of any a_e, so")
print("   the root is unchanged too.  The 2pi test cannot discriminate BY CONSTRUCTION, and the")
print("   verdict's claim that it does is contradicted by the lane's own sealed .OUT file.")
print("   CONSEQUENCE: 'R is a function of a LIFT of the connection to R, which S1 does not carry'")
print("   IS FALSE.  R is a function of W alone, and W is exactly what S1 sec4 publishes as the")
print("   complete gauge-invariant content of the connection.  CRITERION (b) FALLS.\n")

print("== X3.2  WHAT SURVIVES OF (b), STATED CORRECTLY -- AND IT IS NOT NOTHING ==")
print("   The lane's OTHER argument in the same leg is sound and is a different argument:")
print("   no CONTINUOUS h : U(1) -> U(1) has h(W)^L = W (deg(h).L = 1 is unsolvable in Z), so the")
print("   root rule has a cut.  Independent exhibition, with COR-F's T as the control arm:")
th=np.linspace(0,2*np.pi,4001)[:-1]
wr=np.exp(1j*np.angle(np.exp(1j*th))/3)
jr=np.max(np.abs(np.diff(np.append(wr,wr[0]))))
tj=max(np.linalg.norm(T_edge(LF,np.array([t,0,0,0,0,0]))-T_edge(LF,np.array([t+2*np.pi/4000,0,0,0,0,0])))
       for t in th[::40])
print(f"      W swept once around U(1), 4000 points: max step in the root  = {jr:.4f}  (a JUMP)")
print(f"                                             max step in COR-F's T = {tj:.2e}  (continuous)")
print("   AND THAT REPAIRS CRITERION (a), WHICH THE LANE DISMISSED AS '(b) IN DISGUISE'.  (a) says")
print("   the tick is a graph automorphism composed with EDGE PHASES, i.e. its entries are")
print("   monomials prod U_e^{m_e}, m_e in Z.  Monomials are exactly the CONTINUOUS characters of")
print("   the connection torus; w is discontinuous, so w is no monomial.  Exhaustive check:")
found=[]
aset=[rng.uniform(0,2*np.pi,6) for _ in range(12)]
for m in itertools.product(range(-3,4),repeat=3):
    ok=all(abs(np.exp(1j*(m[0]*aa[0]+m[1]*aa[1]+m[2]*aa[2])) -
               np.exp(1j*np.angle(holF(aa))/3))<1e-9 for aa in aset)
    if ok: found.append(m)
print(f"      integer exponent vectors m in [-3,3]^3 with prod U_e^m = W^(1/3) on 12 connections: {found}")
print("      (none exists for any range: 3m_e = 1 has no integer solution).  So (a) DISCRIMINATES")
print("      and is NOT (b) in disguise.  THE LANE UNDERSOLD THE ONE CRITERION THAT SURVIVES ITS")
print("      OWN ARITHMETIC AND OVERSOLD THE ONE THAT DOES NOT.\n")

print("== X3.3  A CRITERION NEITHER LANE PROPOSED: IS THE TICK LOCAL *IN THE CONNECTION*? ==")
print("   S1 sec3 gives the connection edge by edge.  A tick that moves a value across edge e")
print("   should depend on a_e and on nothing else.  Jacobian support, computed numerically:")
def entry_deps(build, a, eps=1e-6):
    """for each non-zero, non-trivial entry of the tick, how many edge phases does it depend on?"""
    base=build(a); out={}
    for u in range(NV):
        for v in range(NV):
            if abs(base[u,v])<1e-12 or (u==v and abs(base[u,v]-1)<1e-12): continue
            k=0
            for e in range(6):
                ap=a.copy(); ap[e]+=eps
                if abs(build(ap)[u,v]-base[u,v])/eps>1e-6: k+=1
            out[(u,v)]=k
    return out
dT=entry_deps(lambda x: T_edge(LF,x), a)
dR=entry_deps(lambda x: root_op(holF(x),VF,3)[0], a)
dM=entry_deps(lambda x: M_of(VF,holF(x)), a)
print(f"      COR-F tick T : entries {sorted(dT)} depend on {sorted(dT.values())} edge phases each")
print(f"      root tick R  : entries {sorted(dR)} depend on {sorted(dR.values())} edge phases each")
print(f"      circuit M    : entries {sorted(dM)} depend on {sorted(dM.values())} edge phases each")
print(f"      max edges per entry:  T = {max(dT.values())}   R = {max(dR.values())}   M = {max(dM.values())}")
print("   Each ENTRY of T is a single U_e; R's every entry needs all three.  So a sub-circuit tick")
print("   built as an L-th root of the holonomy is NON-LOCAL IN THE CONNECTION at sub-circuit")
print("   resolution: it must consult edges the fibre has not reached yet.  M_gamma is non-local")
print("   too, but M_gamma is a WHOLE-CIRCUIT operator, so for it that is not a defect.")
print("   THIS IS A DISCRIMINATOR THAT IS NEITHER (d) NOR CONTINUITY, AND IT IS THE CARRIER'S OWN")
print("   DATA LAYOUT (S1:16-22, one phase per edge).  It repairs part of what X3.1 removed.\n")

print("== X3.4  CRITERION (d), AND ITS CIRCULARITY, TESTED RATHER THAN CONCEDED ==")
ok_T=ok_R=0.0
for k in range(1,7):
    Tk=np.linalg.matrix_power(T,k); Rk=np.linalg.matrix_power(R,k)
    for i,(s_,d_,e) in enumerate(LF):
        tgt=LF[(i+k-1)%3][1]
        P=np.prod([np.exp(1j*a[LF[(i+j)%3][2]]) for j in range(k)])
        ok_T=max(ok_T,abs(Tk[tgt,s_]-P)); ok_R=max(ok_R,abs(Rk[tgt,s_]-P))
print(f"      max |T^k[v_{{i+k}},v_i] - transport along that k-path|, k<=6 = {ok_T:.2e}  (reproduced)")
print(f"      max |R^k[v_{{i+k}},v_i] - transport along that k-path|, k<=6 = {ok_R:.2e}  (reproduced)")
print("   THE CIRCULARITY IS REAL AND I DO NOT LET THE LANE OFF IT: (d) says 'a tick must be one")
print("   edge of parallel transport', which is 'T is canonical because it is parallel transport'.")
print("   A reader who holds the object is the HOLONOMY -- a closed-loop invariant -- rejects (d)")
print("   outright.  BUT THE CORPUS ITSELF DOES NOT GET THAT READING FOR FREE: its own W-01 row")
print("   extends the transport off L_v0 to Gamma(L) = C^5 precisely so it can act on SECTIONS,")
print("   and S3 needs Omega_N = prod Z_k, a PROCESS in time.  (d) is not neutral, and neither is")
print("   its denial.  X3.3 is the discriminator that does NOT presuppose the process reading.\n")

print("== X3.5  DISCOUNTING THE LANE'S KILL OF CHOICE LEDGER A1 -- ONE WITNESS SMUGGLES, ONE DOES NOT ==")
print("   A1 (S2 audit :657) claims M_gamma is the UNIQUE extension that is (a) unitary, (b) uses")
print("   no data beyond fibres/edges/orientation/connection, (c) reduces to W. on L_v0.")
print("   The lane's X1 = M_gamma . (phased swap of v1,v2 ALONG e2) SMUGGLES: it singles out ONE")
print("   edge of the loop.  The carrier supplies the edge LIST, not a marked edge, and A1's own")
print("   clause (b) is exactly about what data may be used.  I DISCOUNT X1.  Its X2 = diag(W, W^2,")
print("   conj(W), 1, 1) needs an ORDERING of the loop's vertices -- which orientation does supply")
print("   (S1:27), so X2 stands.  AND HERE IS A THIRD THAT NEEDS NOTHING AT ALL:")
X3 = np.diag([WF, np.conj(WF), np.conj(WF), 1, 1]).astype(complex)   # W at the root, conj(W) elsewhere
th=rng.uniform(0,2*np.pi,5); ag=np.array([a[j]+th[t]-th[s] for j,(s,t) in enumerate(EDG)])
Xg = np.diag([holF(ag), np.conj(holF(ag)), np.conj(holF(ag)), 1, 1]).astype(complex)
Dm = np.diag(np.exp(1j*th))
print(f"      X3 = diag(W, conj(W), conj(W), 1, 1)")
print(f"         unitary  ||X*X - I||                = {np.linalg.norm(X3.conj().T@X3-np.eye(5)):.2e}")
print(f"         gauge-covariant ||X(a^th) - D X D*|| = {np.linalg.norm(Xg-Dm@X3@np.linalg.inv(Dm)):.2e}")
print(f"         (c) |X[v0,v0] - W|                  = {abs(X3[0,0]-WF):.2e}")
print(f"         acts on every loop vertex, fixes v3,v4: {np.allclose(X3[3:,3:],np.eye(2))}")
print(f"         uses NO ordering, NO marked edge, NO root choice -- only W and the loop's")
print(f"         vertex SET, and the root is distinguished by the incidence itself (S1:14).")
print(f"         || X3 - M_gamma || = {np.linalg.norm(X3-M):.3f}")
print("      -> A1's UNIQUENESS IS FALSE, AND IT IS FALSE WITHOUT SMUGGLING ANYTHING.  The lane's")
print("         finding SURVIVES my discount, on a cleaner witness than either of its own.")
print("         The missing axiom on THIS witness is HOMOGENEITY (the same scalar, and specifically")
print("         W itself, at every loop vertex).  FIBRE-WISE-NESS is needed only to kill the lane's")
print("         discounted X1, so its 'A1 NEEDS TWO AXIOMS' is one axiom too many on the witnesses")
print("         that survive scrutiny.  THE CONCLUSION IS UNHARMED: W-10 0.4.1's phrase is")
print("         'FIBRE-WISE SCALAR MULTIPLICATION BY THE WHOLE-CIRCUIT HOLONOMY', and homogeneity")
print("         is its second half.  The lane's finding lands -- on the other half of the phrase")
print("         than it named.")
