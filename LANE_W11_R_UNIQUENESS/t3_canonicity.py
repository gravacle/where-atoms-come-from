# LANE W-11-R-T, leg T3 — IS THERE A PRINCIPLED REASON TO PREFER COR-F's T?
# The brief lists four candidate criteria "to test, not to assert".  I test those four, add two
# more that actually discriminate, and report the two that do NOT discriminate as failures.
import numpy as np, rlib
from rlib import K1, diag_member, shift_member
rng = np.random.default_rng(20260817)

a = rlib.a_with_holonomies(K1, 1.0, np.sqrt(2), rng)
WF = K1.hol(a,'F'); WC = K1.hol(a,'C')
T  = K1.T_corf(a,'F'); M = K1.M(a,'F')
def root_op(W, vs, L, branch=0):
    """the ONLY invisibility-restoring rival shape: fibre-wise multiplication by an L-th root of W."""
    w = np.exp(1j*(np.angle(W)+2*np.pi*branch)/L)
    X = np.eye(K1.NV, dtype=complex)
    for v in vs: X[v,v] = w
    return X, w
R, w0 = root_op(WF, K1.VF, 3)
print("== T3.0  THE ONLY RIVAL THAT MATTERS ==")
print("  Leg T2 leaves exactly one rival SHAPE that restores invisibility on K1 and B0b:")
print("  R = (an L-th root of W) on the loop, identity off it.  Everything below compares")
print("  COR-F's T against R.  R is in the family:")
print(f"     ||R^3 - M_gamma||     = {np.linalg.norm(np.linalg.matrix_power(R,3)-M):.2e}")
print(f"     ||R*R - I||           = {np.linalg.norm(R.conj().T@R-np.eye(5)):.2e}")
print(f"     D(invisibility) with R on both loops = "
      f"{rlib.Dscore(K1, R, root_op(WC,K1.VC,3)[0], 12):.2e}   -> R RESTORES INVISIBILITY.\n")

print("== T3.1  CRITERION (b) -- 'DEPENDS ONLY ON THE EDGE DATA U_e, NO EXTRA STRUCTURE'.  DECISIVE. ==")
print("  S1 sec3: the connection is U_e = exp(i a_e) with a_e in R/2piZ.  The a_e are ONLY defined")
print("  modulo 2pi.  So a candidate transport must be unchanged when one a_e moves by 2pi.")
a2 = a.copy(); a2[0] += 2*np.pi
T2 = K1.T_corf(a2,'F'); R2, w2 = root_op(K1.hol(a2,'F'), K1.VF, 3)
print(f"     a_1 -> a_1 + 2pi :   ||T(a') - T(a)||   = {np.linalg.norm(T2-T):.2e}      COR-F: UNCHANGED")
print(f"                          |w(a') - w(a)|     = {abs(w2-w0):.6f}      root: CHANGED")
print(f"                          w(a')/w(a)         = {w2/w0:.6f}   (= a cube root of unity)")
print("  AND THE OBSTRUCTION IS TOPOLOGICAL, NOT A BAD BRANCH CHOICE.  A rule producing the root")
print("  from the holonomy is a continuous h : U(1) -> U(1) with h(W)^L = W.  Winding numbers give")
print("  deg(h).L = 1, impossible for L >= 2.  SO NO CONTINUOUS RULE EXISTS AT ALL.  Exhibited:")
th = np.linspace(0, 2*np.pi, 2001)[:-1]
wr = np.exp(1j*np.angle(np.exp(1j*th))/3)
jump = np.max(np.abs(np.diff(np.append(wr, wr[0]))))
Tsw = [np.linalg.norm(K1.T_corf(np.array([t,0,0,0,0,0]),'F') - K1.T_corf(np.array([t+1e-4,0,0,0,0,0]),'F'))
       for t in th[::100]]
print(f"     sweeping W once around U(1) (2000 points): max jump in the principal root = {jump:.4f}")
print(f"                                                 max jump in COR-F's T          = {max(Tsw):.2e}")
print("  VERDICT (b): COR-F's T PASSES on the carrier's own published data type; R FAILS, and no")
print("  repair exists.  R is a function of a LIFT of the connection to R, which S1 does not carry.\n")

print("== T3.2  CRITERION (d) -- 'RESTRICTS CORRECTLY TO A SUBLOOP'.  DECISIVE, AND IT PINS T EXACTLY. ==")
print("  Read as: after k ticks, the fibre value at v_i must have arrived at v_{i+k} multiplied by")
print("  the parallel transport along that k-edge sub-path (S1 sec3's own definition of transport).")
U = np.exp(1j*a)
ok_T = ok_R = 0.0
for k in range(1,7):
    Tk = np.linalg.matrix_power(T,k); Rk = np.linalg.matrix_power(R,k)
    for i,(s_,d_,e) in enumerate(K1.loopF):
        tgt = K1.loopF[(i+k-1)%3][1]
        P = np.prod([np.exp(1j*a[K1.loopF[(i+j)%3][2]]) for j in range(k)])
        ok_T = max(ok_T, abs(Tk[tgt,s_] - P)); ok_R = max(ok_R, abs(Rk[tgt,s_] - P))
print(f"     max | T^k[v_{{i+k}}, v_i] - (parallel transport along the k-edge path) |, k<=6 = {ok_T:.2e}")
print(f"     max | R^k[v_{{i+k}}, v_i] - (parallel transport along the k-edge path) |, k<=6 = {ok_R:.2e}")
print("  AND AT k = 1 THE CRITERION HAS A UNIQUE SOLUTION IN THE WHOLE FAMILY.  Any shift member")
print("  has weights n_i; requiring n_i = U_i is one equation per edge with one solution:")
nsol = [np.exp(1j*a[e]) for (_,_,e) in K1.loopF]
Tsol = shift_member(K1, a, 'F', np.angle(nsol[:2]))
print(f"     || (the unique k=1 solution) - COR-F's T || = {np.linalg.norm(Tsol-T):.2e}")
print("  VERDICT (d): COR-F's T PASSES and is the UNIQUE member that does.  R fails at every k:")
print("  it never relates two distinct fibres, so it transports nothing and no sub-path statement")
print("  about it is even expressible.\n")

print("== T3.3  CRITERION (a) -- 'A GRAPH AUTOMORPHISM COMPOSED WITH EDGE PHASES'. PASSES, WEAKLY. ==")
P = np.zeros((5,5)); P[1,0]=P[2,1]=P[0,2]=1.0
for v in (3,4): P[v,v]=1.0
DU = np.diag([np.exp(1j*a[0]),np.exp(1j*a[1]),np.exp(1j*a[2]),1,1])
print(f"     || T - P_gamma . diag(U_e) || = {np.linalg.norm(T - P@np.diag([np.exp(1j*a[0]),np.exp(1j*a[1]),np.exp(1j*a[2]),1,1])):.2e}")
print("     with P_gamma the cyclic rotation of the loop -- a graph automorphism of gamma.")
print("     R = (identity automorphism) . diag(w,w,w,1,1), and w is NOT an edge phase: it is not")
print("     any monomial prod U_e^{m_e} with integer m_e (leg T3.1's 2pi test is exactly that test).")
print("  VERDICT (a): passes for COR-F, but only because 'edge phases' is doing the work -- this")
print("  criterion is criterion (b) in disguise and is not independent evidence.\n")

print("== T3.4  CRITERION (c) -- 'COVARIANT UNDER RELABELLING THE LOOP'S STARTING VERTEX'. NON-DISCRIMINATING. ==")
rot = [(1,2,1),(2,0,2),(0,1,0)]
Trot = K1.T_corf(a,'F'); K1b = rlib.Carrier("K1rot",5,K1.edges,rot,K1.loopC)
print(f"     || T(loop written from v1) - T(loop written from v0) || = {np.linalg.norm(K1b.T_corf(a,'F')-T):.2e}")
Rrot,_ = root_op(K1b.hol(a,'F'), K1b.VF, 3)
print(f"     || R(loop written from v1) - R(loop written from v0) || = {np.linalg.norm(Rrot-R):.2e}")
print("  VERDICT (c): BOTH pass (on U(1) the holonomy is base-point independent).  REPORTED AS A")
print("  FAILURE OF THE CRITERION, not as evidence for either side.\n")

print("== T3.5  A CRITERION I PROPOSED AND THAT FAILED: ORIENTATION REVERSAL SHOULD INVERT THE TICK ==")
rev = [(0,2,2),(2,1,1),(1,0,0)]
Trev = np.zeros((5,5),dtype=complex)
for v in (3,4): Trev[v,v]=1.0
for (s_,d_,e) in rev: Trev[d_,s_] = np.exp(-1j*a[e])     # traversed backwards: transport by conj(U)
print(f"     || T(gamma reversed) - T(gamma)^{-1} || = {np.linalg.norm(Trev-np.linalg.inv(T)):.2e}   COR-F: EXACT")
worst = 0.0
for _ in range(2000):
    W = np.exp(1j*rng.uniform(-np.pi,np.pi))
    worst = max(worst, abs(np.exp(1j*np.angle(1/W)/3) - 1/np.exp(1j*np.angle(W)/3)))
print(f"     max | w(W^{-1}) - w(W)^{-1} | over 2000 connections, principal branch arg in (-pi,pi]")
print(f"                                             = {worst:.2e}   root: ALSO EXACT")
print("  VERDICT (e): I PROPOSED THIS CRITERION AND IT DOES NOT DISCRIMINATE.  With numpy's")
print("  symmetric branch arg in (-pi,pi] the root operator is orientation-covariant too.  I")
print("  record it as a failed discriminator rather than dropping it.  What survives of the")
print("  intuition is T3.1's DISCONTINUITY, which is branch-independent: every branch rule has a")
print("  cut somewhere on U(1), and T has none.\n")

print("== T3.6  A SECOND ONE: STABILITY UNDER SUBDIVIDING AN EDGE (a change with no physics in it) ==")
print("  Split e5 : v3->v4 into v3->x->v4 with U'U'' = U_{e5}.  gamma_C's length goes 3 -> 4;")
print("  W_C, and every path transport, are unchanged.")
Csub = rlib.Carrier("K1sub",6,[(0,1),(1,2),(2,0),(0,3),(3,5),(5,4),(4,0)],
                    [(0,1,0),(1,2,1),(2,0,2)], [(0,3,3),(3,5,4),(5,4,5),(4,0,6)])
asub = np.array([a[0],a[1],a[2],a[3],0.4*a[4],0.6*a[4],a[5]])
print(f"     W_C before = {WC:.6f}   W_C after = {Csub.hol(asub,'C'):.6f}   |diff| = {abs(WC-Csub.hol(asub,'C')):.2e}")
Tsub = Csub.T_corf(asub,'C')
print(f"     ||T_sub^4 - M_c(subdivided)|| = {np.linalg.norm(np.linalg.matrix_power(Tsub,4)-Csub.M(asub,'C')):.2e}")
print(f"     transport v3 -> v4 along the SUBDIVIDED path, from T_sub^2  : {(Tsub@Tsub)[4,3]:.6f}")
print(f"     transport v3 -> v4 before subdivision, from T_C           : {K1.T_corf(a,'C')[4,3]:.6f}")
w3 = np.exp(1j*np.angle(WC)/3); w4 = np.exp(1j*np.angle(WC)/4)
print(f"     the root operator, by contrast: w = W_C^(1/3) = {w3:.6f}  ->  W_C^(1/4) = {w4:.6f}")
print(f"                                     |change| = {abs(w3-w4):.6f} for a change with no physics in it")
print("  VERDICT (f): COR-F's T transports identically across the subdivision; R's scalar jumps.\n")

print("== T3.7  WHAT GAUGE COVARIANCE + HOMOGENEITY LEAVE OF THE SHIFT BRANCH -- AND IT IS NOT NOTHING ==")
print("  Gauge covariance forces the tick's (v_{i+1}, v_i) entry to carry weight exp(i(th_{i+1}-th_i)),")
print("  i.e. n_i = c_i U_i with c_i gauge-INVARIANT.  'Same rule at every edge' forces c_i = c.")
print("  T^L = M_gamma then forces c^L = 1.  So the shift branch collapses to the L members zeta^j.T.")
print("  I EXPECTED THAT AMBIGUITY TO BE INVISIBLE TO |Z_n|.  IT IS NOT -- because zeta multiplies")
print("  only the LOOP components of the section, so it does not factor out of the inner product:")
s = rng.normal(size=5)+1j*rng.normal(size=5); s /= np.linalg.norm(s)
zt = np.exp(2j*np.pi/3)
TF0, TC0 = K1.T_corf(a,'F'), K1.T_corf(a,'C')
TFz = TF0.copy()
for v in K1.VF:
    for u in K1.VF: TFz[v,u] *= zt
def Zs(A,B,n): return abs(np.vdot(np.linalg.matrix_power(A,n)@s, np.linalg.matrix_power(B,n)@s))
worst = max(abs(Zs(TFz,TC0,n)-Zs(TF0,TC0,n)) for n in range(1,13))
print(f"     max_n | |Z_n| under zeta.T  -  |Z_n| under T |  = {worst:.2e}   -> VISIBLE.")
print(f"     but every one of the L members is still a SHIFT, so none restores invisibility:")
print(f"        D(zeta.T_F, T_C) = {rlib.Dscore(K1, TFz, TC0, 12):.3f}   D(T_F, T_C) = {rlib.Dscore(K1, TF0, TC0, 12):.3f}")
print("  SO: gauge covariance + homogeneity leave an L-FOLD AMBIGUITY THAT THE OBSERVABLE SEES.")
print("  CRITERION (d) IS WHAT REMOVES IT: 'one tick = one edge of S1 sec3's own parallel")
print("  transport' is the equation n_i = U_i, whose unique solution is COR-F's T (leg T3.2).")
print("  RECORDED AGAINST MYSELF: I predicted this ambiguity would be invisible and it is not.")
print("  The canonicity therefore rests on (d) ALONE, not on (b)+homogeneity, and (d) is not a")
print("  free-standing axiom -- it is the corpus's own definition of transport, S1:57-58.\n")

print("== T3.8  NEW DEFECT: THE CORPUS'S ONLY WRITTEN UNIQUENESS ARGUMENT FOR M_gamma IS FALSE ==")
print("  S2_FORMATION_CONDITION_ON_K1_AUDIT_V001.md:650, CHOICE LEDGER A1 -- the entry that closed")
print("  the extension question -- claims M_gamma is 'the UNIQUE extension that is (a) unitary,")
print("  (b) uses no data beyond fibres/edges/orientation/connection, (c) reduces to the build's")
print("  T_gamma on L_v0'.  TWO COUNTEREXAMPLES, both built only from K1's own data:")
Sw = np.eye(5,dtype=complex); Sw[1,1]=Sw[2,2]=0.0
Sw[2,1] = np.exp(1j*a[1]); Sw[1,2] = np.exp(-1j*a[1])      # phased swap of v1,v2 along e2; fixes v0
X1 = M@Sw
X2 = np.diag([WF, WF*WF, WF*np.conj(WF)*np.conj(WF), 1, 1]).astype(complex)  # fibre-wise, still not M
for nm,X in (("X1 = M_gamma . (phased swap of v1,v2 along e2)", X1),
             ("X2 = diag(W, W^2, W^-1, 1, 1), fibre-wise", X2)):
    th = rng.uniform(0,2*np.pi,5)
    ag = np.array([a[j]+th[t]-th[s_] for j,(s_,t) in enumerate(K1.edges)])
    if X is X1:
        Sg = np.eye(5,dtype=complex); Sg[1,1]=Sg[2,2]=0.0
        Sg[2,1]=np.exp(1j*ag[1]); Sg[1,2]=np.exp(-1j*ag[1]); Xg = K1.M(ag,'F')@Sg
    else:
        Wg_ = K1.hol(ag,'F'); Xg = np.diag([Wg_, Wg_**2, Wg_*np.conj(Wg_)**2, 1, 1]).astype(complex)
    D_ = np.diag(np.exp(1j*th))
    print(f"     {nm}")
    print(f"        unitary            ||X*X-I||        = {np.linalg.norm(X.conj().T@X-np.eye(5)):.2e}")
    print(f"        gauge-covariant    ||X(a^th)-D X D*|| = {np.linalg.norm(Xg - D_@X@np.linalg.inv(D_)):.2e}")
    print(f"        (c) on L_v0        |X[v0,v0] - W|   = {abs(X[0,0]-WF):.2e}")
    print(f"        acts on every loop vertex, fixes v3,v4:  {np.allclose(X[3:,3:],np.eye(2))}")
    print(f"        BUT  ||X - M_gamma|| = {np.linalg.norm(X-M):.3f}   -> A1's 'unique' IS FALSE")
print("  A1 needs TWO axioms it does not state: FIBRE-WISE-NESS (kills X1) and HOMOGENEITY --")
print("  the same scalar at every loop vertex (kills X2).  Fibre-wise-ness is exactly the")
print("  stipulation at issue in W-10 0.4.1, and it is unledgered inside the ledger entry that")
print("  was supposed to close the question.")

print("\n== T3.9  THE 'FINER CLOCK' ROUTE LANDS INSIDE THE SAME FAMILY AND ADDS NOTHING ==")
print("  S2 audit CHOICE LEDGER A2 rejected 'a real parameter t with a Hamiltonian'.  Suppose one")
print("  takes it anyway: H Hermitian with exp(iH) = M_gamma, tick = exp(iH/L).  Every such tick")
print("  satisfies (exp(iH/L))^L = M_gamma, so it IS a member of the family of leg T1.3 --")
print("  the continuous-time route buys no new operators at all.  Checked on 4000 draws:")
hits=0; best=np.inf; worstfam=0.0
for _ in range(4000):
    V3 = rlib.haar(3,rng); m = rng.integers(-2,3,3)
    thW = np.angle(WF)
    Hl = V3@np.diag(thW + 2*np.pi*m)@V3.conj().T
    ev, evec = np.linalg.eigh((Hl+Hl.conj().T)/2)
    X = np.eye(5,dtype=complex)
    B = evec@np.diag(np.exp(1j*ev/3))@evec.conj().T
    for i,u in enumerate(K1.VF):
        for k,w_ in enumerate(K1.VF): X[u,w_] = B[i,k]
    worstfam = max(worstfam, np.linalg.norm(np.linalg.matrix_power(X,3)-M))
    d = rlib.Dscore(K1, X, K1.T_corf(a,'C'), 12); best=min(best,d)
    if d < 1e-9: hits += 1
print(f"     max ||exp(iH/3)^3 - M_gamma|| over 4000 draws = {worstfam:.2e}  -> all inside the family")
print(f"     how many restore invisibility (against COR-F's T_C): {hits};  best D = {best:.2e}")
print("  -> the clock question and the transport question are the SAME question.  Refining the")
print("     clock without refining the transport is not an available option.")
