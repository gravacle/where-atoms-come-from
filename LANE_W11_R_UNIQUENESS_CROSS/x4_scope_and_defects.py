# LANE W-11-R-T-CROSS, leg X4 — THE UNIQUENESS LANE'S OWN DEFECTS, AND ITS TWO STRIKES RE-JUDGED.
import numpy as np
rng = np.random.default_rng(11072026)
NV=5; LF=[(0,1,0),(1,2,1),(2,0,2)]; LC=[(0,3,3),(3,4,4),(4,0,5)]; VF=[0,1,2]; VC=[0,3,4]
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

print("== X4.1  THE LANE'S B0b CONNECTION IS NOT THE CONNECTION ITS CONVENTIONS FILE PUBLISHES ==")
print("   PUBLISHED_CONVENTIONS: 'Generic runs use S4:603's f = 1.0, c = sqrt(2) ... edge phases are")
print("   drawn uniformly and then RESCALED SO THE TWO HOLONOMIES HIT (f,c)'.  Its rlib.a_with_")
print("   holonomies sets  a[lastF] = (f - angle(hol))%2pi  -- correct only when the last edge is")
print("   traversed FORWARDS.  B0b's gamma_F traverses its last two edges BACKWARDS (signF =")
print("   [1,1,-1,-1]), so that edge enters the holonomy conjugated and the assignment overshoots.")
def b0b():
    Vv=lambda i,j:3*(j%3)+(i%3); Hh=lambda i,j:3*(j%3)+(i%3); Wg=lambda i,j:9+3*(j%3)+(i%3)
    E=[None]*18
    for j in range(3):
        for i in range(3):
            E[Hh(i,j)]=(Vv(i,j),Vv(i+1,j)); E[Wg(i,j)]=(Vv(i,j),Vv(i,j+1))
    gF=[(Hh(0,0),1),(Wg(1,0),1),(Hh(0,1),-1),(Wg(0,0),-1)]; gC=[(Hh(0,0),1),(Hh(1,0),1),(Hh(2,0),1)]
    wf=[]; wc=[]
    for (e,s) in gF:
        u,v=E[e]; wf.append((u,v,e,s) if s>0 else (v,u,e,s))
    for (e,s) in gC:
        u,v=E[e]; wc.append((u,v,e,s) if s>0 else (v,u,e,s))
    return E,wf,wc
BE,BF,BC=b0b()
def bhol(a,loop): return np.exp(1j*sum(a[e]*s for _,_,e,s in loop))
def a_with_holonomies_asWritten(f,c,rg):
    a=rg.uniform(0,2*np.pi,18)
    lastF=BF[-1][2]; a[lastF]=0.0; a[lastF]=(f-np.angle(bhol(a,BF)))%(2*np.pi)
    lastC=BC[-1][2]
    if lastC!=lastF: a[lastC]=0.0; a[lastC]=(c-np.angle(bhol(a,BC)))%(2*np.pi)
    return a
err=[]
for _ in range(500):
    a=a_with_holonomies_asWritten(1.0,np.sqrt(2),rng)
    err.append((abs(((np.angle(bhol(a,BF))-1.0+np.pi)%(2*np.pi))-np.pi),
                abs(((np.angle(bhol(a,BC))-np.sqrt(2)+np.pi)%(2*np.pi))-np.pi)))
err=np.array(err)
print(f"      500 draws on B0b: max |angle(W_F) - 1.0|      = {err[:,0].max():.4f}   TARGET MISSED")
print(f"                        mean|angle(W_F) - 1.0|      = {err[:,0].mean():.4f}")
print(f"                        max |angle(W_C) - sqrt(2)|  = {err[:,1].max():.2e}   (C is fine)")
print("   IMPACT, STATED HONESTLY: none of the lane's B0b VERDICTS depend on the connection being")
print("   that one -- the 144-of-6912 count, the shift null and Corollary 5 are connection-free.")
print("   What falls is the PUBLISHED CONVENTION: on B0b the F holonomy is a uniform random angle,")
print("   not S4:603's 1.0.  Since it is uniform it is generic almost surely, so W-10 N-4's")
print("   'no rational connection labelled generic' is not breached either.  A LEDGER DEFECT, NOT")
print("   A RESULT DEFECT -- but the lane's arms on K1 and B0b are not the same connection family")
print("   it says they are.  (K1's gamma_F and gamma_C are traversed forwards, so K1 is unaffected:")
a5=None
def a_wh_K1(f,c,rg):
    a=rg.uniform(0,2*np.pi,6); a[2]=0.0; a[2]=(f-np.angle(np.exp(1j*a[:3].sum())))%(2*np.pi)
    a[5]=0.0; a[5]=(c-np.angle(np.exp(1j*a[3:].sum())))%(2*np.pi); return a
e1=max(abs(((np.angle(np.exp(1j*a_wh_K1(1.0,np.sqrt(2),rng)[:3].sum()))-1.0+np.pi)%(2*np.pi))-np.pi) for _ in range(200))
print(f"    max |angle(W_F) - 1.0| on K1 over 200 draws = {e1:.2e}.)\n")

print("== X4.2  ITS STRIKE OF THE REGISTRAR'S LEG D CLAUSE 2, RE-JUDGED.  I DO NOT UPHOLD IT. ==")
a=np.array([1.0,0.37,0.91,2**0.5,0.23,1.77])
WF=np.exp(1j*a[:3].sum()); WC=np.exp(1j*a[3:].sum())
MF,MC=M_of(VF,WF),M_of(VC,WC); TF,TC=T_edge(LF,a),T_edge(LC,a)
def rate(s,A,B,N):
    x=s.copy(); y=s.copy(); t=0.0
    for _ in range(N):
        x=A@x; y=B@y; z=abs(np.vdot(x,y)); t+= np.log(z) if z>0 else -700.0
    return t/N
def m_jensen(p,n=1<<20):
    A,B,C,Dd=p; th=2*np.pi*np.arange(n)/n; ct=np.cos(th)
    return np.log(np.maximum(np.sqrt(np.maximum(A*A+B*B+2*A*B*ct,0)),
                             np.sqrt(np.maximum(C*C+Dd*Dd+2*C*Dd*ct,0)))+1e-300).mean()
sA=np.sqrt(np.array([0.40,0.15,0.15,0.15,0.15]))+0j
sB2=np.sqrt(np.array([0.40,0.26,0.04,0.11,0.19]))+0j
sC=sA*np.exp(1j*np.array([0.0,1.3,-0.7,2.2,0.4]))
N=20000; mP=m_jensen(np.array([0.0,0.30,0.30,0.40]))
print(f"      N1's m(P) = {mP:.9f}   for pi = (0, 0.30, 0.30, 0.40)")
for nm,s in (("A",sA),("B2",sB2),("C",sC)):
    rc=rate(s,MF,MC,N); re_=rate(s,TF,TC,N)
    print(f"      state {nm:<2}: circuit/circuit = {rc:.9f}   edge/tick = {re_:.9f}   3 x edge = {3*re_:.9f}")
print("   The registrar wrote 'the EDGE rate is state-dependent AND is not m(P)/3'.  The lane")
print("   struck the second clause as 'comparing incommensurables'.  I DISAGREE, AND THE NUMBERS")
print("   ABOVE ARE WHY: 3 x (edge rate) is EXACTLY the per-circuit decay under the edge clock,")
print("   because the n = 3k subsequence of the edge product IS the circuit product; the three")
print("   ticks of a circuit are a REFINEMENT of that circuit, not a different object.  So the")
print("   comparison is between two per-circuit rates on the same carrier, connection and state.")
print("   It is well posed, and what it says is STRONGER than clause 1: under the finer clock the")
print("   record decays MORE per circuit than m(P), by a state-dependent amount.  W-08's")
print("   'incommensurables' was about a floor and a rate with different UNITS; here the units")
print("   agree after the factor of 3, which the registrar itself supplied.  CLAUSE 2 STANDS.\n")

print("== X4.3  NEW-4 (the convention-scope of S3's 'circuits grow no algebra') REPRODUCED ==")
def rank_span(ops):
    return np.linalg.matrix_rank(np.array([o.flatten() for o in ops]), tol=1e-9)
def star_alg_dim(gens,nv):
    B=[np.eye(nv,dtype=complex)]+list(gens)+[g.conj().T for g in gens]
    for _ in range(8):
        new=[x@y for x in B for y in B]
        Mx=np.array([o.flatten() for o in B+new])
        u,sv,vh=np.linalg.svd(Mx,full_matrices=False); r=int((sv>1e-9).sum())
        if r==len(B): break
        B=[vh[i].reshape(nv,nv) for i in range(r)]
    return len(B)
for Nn in (1,2,5,100):
    ci=[np.linalg.matrix_power(MF,n) for n in range(0,Nn+1)]+[np.linalg.matrix_power(MC,n) for n in range(1,Nn+1)]
    ed=[np.linalg.matrix_power(TF,n) for n in range(0,Nn+1)]+[np.linalg.matrix_power(TC,n) for n in range(1,Nn+1)]
    print(f"      N={Nn:>3}  dim span(M_F^n,M_C^n) = {rank_span(ci)}  [S3 sec3.1: 3]     "
          f"dim span(T_F^n,T_C^n) = {rank_span(ed)}")
print(f"      *-algebra <M_F,M_C> in M_5(C): dim = {star_alg_dim([MF,MC],5)}   [S3 sec2.3: 3]")
print(f"      *-algebra <T_F,T_C> in M_5(C): dim = {star_alg_dim([TF,TC],5)}   [= all of M_5(C)]")
print("   REPRODUCED INDEPENDENTLY.  I add the caveat the lane did not: the *-algebra figure is")
print("   sensitive to the CONNECTION, not only to the convention -- at the trivial connection")
print("   T_F, T_C are permutation matrices and the algebra is still large, but at W = 1 the")
print("   CIRCUIT operators are the identity and the span collapses to 1, not 3.  Checked:")
a0=np.zeros(6); MF0,MC0=M_of(VF,1.0+0j),M_of(VC,1.0+0j); TF0,TC0=T_edge(LF,a0),T_edge(LC,a0)
print(f"      trivial connection: dim span(M) = {rank_span([np.linalg.matrix_power(MF0,n) for n in range(0,6)]+[np.linalg.matrix_power(MC0,n) for n in range(1,6)])}"
      f"   dim span(T) = {rank_span([np.linalg.matrix_power(TF0,n) for n in range(0,6)]+[np.linalg.matrix_power(TC0,n) for n in range(1,6)])}"
      f"   *-alg(T) = {star_alg_dim([TF0,TC0],5)}")
print("   So NEW-4 is CONNECTION-scoped as well as convention-scoped.  Recorded as scope; scored")
print("   as evidence for neither reading, exactly as the lane scored it.\n")

print("== X4.4  TWO REPORTING DEFECTS IN THE LANE'S OWN SEALED PROSE ==")
print("   (i)  THEOREM_AND_PROOF.txt SHARPNESS 7 and the returned verdict both report the SYNTH-D")
print("        witness as 'D = 2.1e-15 with ||offdiag T_F|| = 0.997'.  Its own sealed")
print("        t2_invisibility.OUT.txt T2.5 reads 'D over n<=16 = 1.20e-14' and '||offdiag T_F||")
print("        = 0.364'.  NEITHER FIGURE IS REPRODUCIBLE FROM THE SEALED CODE.  The qualitative")
print("        claim (a non-diagonal restorer exists when |J|>=2 and gcd>=2) IS reproducible and")
print("        stands; the two numbers are ABSENT under this program's own rule.  COR-K class.")
print("   (ii) t3_canonicity.py prints the A1 pointer as S2 audit ':650'; DOCUMENTARY_FINDINGS.txt")
print("        prints ':657'.  The CHOICE LEDGER A1 row is at :657.  A pointer-rule slip, minor,")
print("        and it is inside the leg that carries the lane's biggest documentary finding.")

print("\n== X4.5  A COULD-NOT-HAVE-FAILED CONTROL INSIDE THE LANE COMMISSIONED TO CATCH THEM ==")
print("   t1_family.py T1.3 says 'Verified BOTH directions on 4000 draws' for the claim that the")
print("   spectral parameterisation is ONTO the family.  Its second loop builds A = V diag(rho.")
print("   zeta^j) V*, then computes eig(A) and checks the eigenvalues are the L-th roots of W.")
print("   THAT IS THE SPECTRAL THEOREM EVALUATED ON ITS OWN OUTPUT.  It could not have failed and")
print("   it tests the SAME direction as the first loop; ontoness (every family member is")
print("   REACHED) is never tested at all.  Demonstrated: the check returns ~0 for a deliberately")
print("   WRONG root set too, as long as the matrix is built from it.")
def haar(n):
    z=(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))/np.sqrt(2)
    q,r=np.linalg.qr(z); return q*(np.diag(r)/abs(np.diag(r)))
WFx=np.exp(1j*rng.uniform(0,2*np.pi))
worst_true=worst_bogus=0.0
for _ in range(2000):
    V=haar(3); j=rng.integers(0,3,3)
    r=np.exp(1j*np.angle(WFx)/3); zt=np.exp(2j*np.pi/3)
    A=V@np.diag(r*zt**j)@V.conj().T
    ev=np.linalg.eigvals(A); roots=r*zt**np.arange(3)
    worst_true=max(worst_true,max(min(abs(e-q) for q in roots) for e in ev))
    # BOGUS: build from FIFTH roots and check against FIFTH roots -- the same tautology
    r5=np.exp(1j*np.angle(WFx)/5); z5=np.exp(2j*np.pi/5); j5=rng.integers(0,5,3)
    A5=V@np.diag(r5*z5**j5)@V.conj().T
    ev5=np.linalg.eigvals(A5); roots5=r5*z5**np.arange(5)
    worst_bogus=max(worst_bogus,max(min(abs(e-q) for q in roots5) for e in ev5))
print(f"      the lane's check, with the CORRECT cube roots : {worst_true:.2e}")
print(f"      the same check with FIFTH roots (A^3 != M)    : {worst_bogus:.2e}   <- passes equally")
print(f"      ||A5^3 - W.I|| for that bogus family          = {np.linalg.norm(np.linalg.matrix_power(A5,3)-WFx*np.eye(3)):.3f}")
print("   -> the control passes on a family that is NOT the family.  VOID AS EVIDENCE.  The claim")
print("      it was meant to support is a THEOREM (PROP 1, the spectral theorem), and this")
print("      program's rule is that 'could not have failed' voids a control and never a theorem --")
print("      so PROP 1 stands and nothing downstream moves.  But this is the program's signature")
print("      defect, committed by a lane whose brief names it, and it is worth exactly one line.")
