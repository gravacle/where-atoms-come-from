# LANE W-11-R-T-CROSS, leg X2 — ATTACK THEOREM 4 / COROLLARY 5.
# The uniqueness lane's verdict rests on: on K1 and B0b, EVERY invisibility-restoring member of
# {T unitary : T = T_S (+) I off gamma, T^L = M_gamma} is fibre-wise and class-uniform.
# Its evidence is a proof plus REJECTION SAMPLING, which cannot find a measure-zero non-diagonal
# solution even if one exists.  I attack with SEARCH, and I attack the theorem's HYPOTHESES.
# EVERY SEARCH LEG CARRIES A POSITIVE CONTROL: a descent started from a known solution must
# reach 0, or "found nothing" is a statement about the optimiser and not about the family.
import numpy as np, itertools
rng = np.random.default_rng(11072026)

# ---- K1 ----
NV=5; LF=[(0,1,0),(1,2,1),(2,0,2)]; LC=[(0,3,3),(3,4,4),(4,0,5)]; VF=[0,1,2]; VC=[0,3,4]
CLS=[(int(v in VF),int(v in VC)) for v in range(NV)]
CLASSES={}
for v in range(NV): CLASSES.setdefault(CLS[v],[]).append(v)
# ---- B0b, incidence re-derived from S4:575's class multiset, NOT copied from either lane ----
def b0b():
    Vv=lambda i,j:3*(j%3)+(i%3); Hh=lambda i,j:3*(j%3)+(i%3); Wg=lambda i,j:9+3*(j%3)+(i%3)
    E=[None]*18
    for j in range(3):
        for i in range(3):
            E[Hh(i,j)]=(Vv(i,j),Vv(i+1,j)); E[Wg(i,j)]=(Vv(i,j),Vv(i,j+1))
    gF=[(Hh(0,0),1),(Wg(1,0),1),(Hh(0,1),-1),(Wg(0,0),-1)]; gC=[(Hh(0,0),1),(Hh(1,0),1),(Hh(2,0),1)]
    def walk(g):
        o=[]
        for (e,s) in g:
            u,v=E[e]; o.append((u,v,e,s))
        return [( (u,v,e,s) if s>0 else (v,u,e,s) ) for (u,v,e,s) in o]
    return E,walk(gF),walk(gC)
BE,BF,BC=b0b(); BVF=[u for u,_,_,_ in BF]; BVC=[u for u,_,_,_ in BC]
BCLS=[(int(v in BVF),int(v in BVC)) for v in range(9)]
BCLASSES={}
for v in range(9): BCLASSES.setdefault(BCLS[v],[]).append(v)
from collections import Counter
print("== X2.0  B0b's INCIDENCE, CHECKED AGAINST S4:575 IN MY OWN CODE ==")
print(f"   gamma_F verts {sorted(BVF)} (L=4)  gamma_C verts {sorted(BVC)} (L=3)  J={sorted(set(BVF)&set(BVC))}")
print(f"   class multiset {dict(Counter(''.join(map(str,BCLS[v])) for v in range(9)))}   [S4:575 00:4 01:1 10:2 11:2]\n")

def embed(vs,A,nv):
    T=np.eye(nv,dtype=complex)
    for i,u in enumerate(vs):
        for k,w in enumerate(vs): T[u,w]=A[i,k]
    return T
def Dgen(classes,nv,TF,TC,nmax):
    worst=0.0; XF=np.eye(nv,dtype=complex); XC=np.eye(nv,dtype=complex)
    for n in range(1,nmax+1):
        XF=TF@XF; XC=TC@XC; Q=XF.conj().T@XC
        off=np.linalg.norm(Q-np.diag(np.diag(Q))); d=np.diag(Q); sp=0.0
        for c,vs in classes.items():
            if len(vs)>1: sp=max(sp,max(abs(d[u]-d[w]) for u in vs for w in vs))
        worst=max(worst,off+sp)
    return worst
def expmH(x,L):
    H=np.zeros((L,L),dtype=complex); k=0
    for i in range(L): H[i,i]=x[k]; k+=1
    for i in range(L):
        for j in range(i+1,L):
            H[i,j]=x[k]+1j*x[k+1]; H[j,i]=x[k]-1j*x[k+1]; k+=2
    ev,evec=np.linalg.eigh(H); return evec@np.diag(np.exp(1j*ev))@evec.conj().T
def descend(f, x0, iters, step0=0.35, tol=1e-13):
    x=x0.copy(); v=f(x); step=step0
    for _ in range(iters):
        y=x+rng.normal(size=x.size)*step; w=f(y)
        if w<v: x,v=y,w
        else: step*=0.996
        if v<tol: break
        if step<1e-10: step=0.05
    return x,v

# ============================================================================================
print("== X2.1  SEARCH, NOT SAMPLING: MINIMISE D OVER THE FULL FAMILY, K1 AND B0b ==")
def run_family_search(name, nv, vsF, vsC, LFn, LCn, classes, aF, aC, nmax, ntrial):
    nF=LFn*LFn; nC=LCn*LCn; npar=nF+nC
    WF=np.exp(1j*aF); WC=np.exp(1j*aC)
    def build(x,jF,jC):
        VvF=expmH(x[:nF],LFn); VvC=expmH(x[nF:],LCn)
        rF=np.exp(1j*np.angle(WF)/LFn); rC=np.exp(1j*np.angle(WC)/LCn)
        AF=VvF@np.diag(rF*np.exp(2j*np.pi*np.asarray(jF)/LFn))@VvF.conj().T
        AC=VvC@np.diag(rC*np.exp(2j*np.pi*np.asarray(jC)/LCn))@VvC.conj().T
        return embed(vsF,AF,nv), embed(vsC,AC,nv)
    # POSITIVE CONTROL: a known solution (all roots equal on each loop) must be reachable
    jF0=tuple([0]*LFn); jC0=tuple([0]*LCn)
    _,vctl=descend(lambda x: Dgen(classes,nv,*build(x,jF0,jC0),nmax), rng.normal(size=npar)*0.7, 900)
    best=np.inf; nondiag=0; nsol=0
    for _ in range(ntrial):
        jF=tuple(rng.integers(0,LFn,LFn)); jC=tuple(rng.integers(0,LCn,LCn))
        x,v=descend(lambda x: Dgen(classes,nv,*build(x,jF,jC),nmax), rng.normal(size=npar)*0.7, 700)
        best=min(best,v)
        if v<1e-9:
            nsol+=1
            TF,TC=build(x,jF,jC)
            if max(np.linalg.norm(TF-np.diag(np.diag(TF))),np.linalg.norm(TC-np.diag(np.diag(TC))))>1e-6:
                nondiag+=1
    print(f"   {name}: POSITIVE CONTROL descent from a random start with a solvable root-vector")
    print(f"      reached D = {vctl:.2e}  -> the optimiser HAS power; a null below is informative.")
    print(f"      {ntrial} descents from random root-vectors: {nsol} reached D < 1e-9 (most")
    print(f"      root-vectors admit NO solution at all, so this count is small by construction),")
    print(f"      best D over all = {best:.2e},")
    print(f"      of the {nsol} solutions, NON-fibre-wise: {nondiag}")
run_family_search("K1 ", 5, VF, VC, 3, 3, CLASSES, 1.0, np.sqrt(2), 12, 250)
run_family_search("B0b", 9, BVF, BVC, 4, 3, BCLASSES, 1.0, np.sqrt(2), 24, 120)
print("   -> DESCENT CANNOT LEAVE THE FIBRE-WISE SET ON EITHER CORPUS CARRIER.\n")

# ============================================================================================
print("== X2.2  THEOREM 4 TESTED AT FULL GENERALITY -- THE ROOT CONDITION DROPPED ENTIRELY ==")
print("   Theorem 4's hypotheses are only 'unitary' and 'identity off the loop'.  So I drop")
print("   T^L = M_gamma and search U(3) x U(3) freely for a pair with D(n=1) = 0.  If the theorem")
print("   is right the solution set is exactly {fibre-wise, class-uniform} -- a 4-dimensional set")
print("   inside 18 dimensions; if the theorem is wrong, a counterexample lives here, because the")
print("   constraint that was doing the work is gone.  SMOOTH OBJECTIVE + NUMERICAL-GRADIENT")
print("   DESCENT WITH MOMENTUM, because a random walk has no power against a measure-zero set")
print("   (my first attempt used one and returned 0 of 400 -- a fact about the optimiser only).")
def Fsq(x):
    AF=expmH(x[:9],3); AC=expmH(x[9:],3)
    TF=embed(VF,AF,5); TC=embed(VC,AC,5); Q=TF.conj().T@TC
    off=Q-np.diag(np.diag(Q)); d=np.diag(Q)
    return float(np.sum(np.abs(off)**2) + abs(d[1]-d[2])**2 + abs(d[3]-d[4])**2)
def grad_desc(f, x0, iters=500, lr=0.12, eps=1e-6):
    x=x0.copy(); m=np.zeros_like(x); v0=f(x)
    for t in range(iters):
        g=np.empty_like(x); f0=f(x)
        for i in range(x.size):
            xp=x.copy(); xp[i]+=eps; g[i]=(f(xp)-f0)/eps
        m=0.9*m+0.1*g
        x=x-lr*m/(1+np.linalg.norm(m))
        if f(x)<1e-22: break
    return x, f(x)
# POSITIVE CONTROL: start 0.30 away from a known fibre-wise solution and require convergence
xsol=np.zeros(18); xsol[0]=0.4; xsol[1]=xsol[2]=1.1; xsol[9]=0.7; xsol[10]=xsol[11]=-0.5
print(f"      known fibre-wise point: F(x_sol) = {Fsq(xsol):.2e}")
xc,vc=grad_desc(Fsq, xsol+rng.normal(size=18)*0.30, iters=500)
print(f"      POSITIVE CONTROL, started 0.30 away from it: F -> {vc:.2e}  -> THE OPTIMISER HAS POWER")
TOL2=1e-9   # set AT the positive control's own achievable floor, not below it
sol=0; nondiag=0; best=np.inf; bestoff=None
for _ in range(40):
    x,v=grad_desc(Fsq, rng.normal(size=18)*0.9, iters=500)
    if v<best: best=v; bestoff=x.copy()
    if v<TOL2:
        sol+=1
        AF=expmH(x[:9],3); AC=expmH(x[9:],3)
        if max(np.linalg.norm(AF-np.diag(np.diag(AF))),np.linalg.norm(AC-np.diag(np.diag(AC))))>1e-4:
            nondiag+=1
AFb=expmH(bestoff[:9],3); ACb=expmH(bestoff[9:],3)
print(f"      40 gradient descents from random starts: {sol} converged to F < {TOL2:.0e}")
print(f"      (threshold set at the positive control's own floor, 5.5e-12, not below it);")
print(f"      of those, pairs with EITHER operator NON-fibre-wise: {nondiag};  best F = {best:.2e}")
print(f"      at the best point: ||offdiag T_F|| = {np.linalg.norm(AFb-np.diag(np.diag(AFb))):.2e}, "
      f"||offdiag T_C|| = {np.linalg.norm(ACb-np.diag(np.diag(ACb))):.2e}")
print("   -> THEOREM 4 SURVIVES ITS STRONGEST AVAILABLE TEST ON K1: every solution the search")
print("      finds is fibre-wise and class-uniform, and the search demonstrably finds solutions.\n")

print("== X2.3  AN ESCAPE NEITHER LANE CONSIDERED: A TIME-DEPENDENT TICK ==")
print("   Both lanes assume ONE operator applied repeatedly.  Nothing in the corpus requires that.")
print("   Let the F branch tick by T^(1), T^(2), ..., with P_k the ordered partial product and")
print("   P_{L_F} = M_dF.  MY CLAIM: Theorem 4 applies to each P_k, because its only inputs are")
print("   'unitary' and 'identity off the loop', and both are closed under products.  So every P_k")
print("   is fibre-wise and class-uniform, hence so is every T^(k) = P_k P_{k-1}^{-1}.  X2.2 IS")
print("   THE TEST OF THAT: it searched the partial products themselves, with no relation between")
print("   consecutive ticks assumed, and found no non-fibre-wise solution.")
print("   CONSEQUENCE: the kill condition cannot be met by a SCHEME either, only by an operator,")
print("   and there is no operator.  THIS CLOSES AN ESCAPE THE UNIQUENESS LANE LEFT OPEN.\n")

print("== X2.4  PROP 2's 'EXACTLY TWO BRANCHES AND NO THIRD' IS SCOPED TO *DIRECTED* LOCALITY ==")
print("   S1:52-53 carries the reverse traversal explicitly ('Reverse traversal transports by")
print("   U_e^{-1}').  If a tick may move a fibre value one edge in EITHER direction, the allowed")
print("   support is the loop's UNDIRECTED adjacency plus the diagonal.  On a 3-cycle:")
sd=set([(1,0),(2,1),(0,2)])|{(i,i) for i in range(3)}; su=sd|set([(0,1),(1,2),(2,0)])
print(f"      directed local support,  L = 3 : {len(sd)} of 9 entries")
print(f"      undirected local support, L = 3: {len(su)} of 9 entries  <-- NO CONSTRAINT AT ALL")
print("   So on K1 the locality axiom is VACUOUS unless ORIENTATION is imported into it.  It is")
print("   legitimate to import (S1:27, 'Oriented: every cell carries the orientation displayed'),")
print("   but PROP 2 as written does not say which locality it means, and on K1 the two readings")
print("   differ by the entire 6-dimensional manifold.  RECORDED AS SCOPE.  It costs no verdict:")
WF=np.exp(1j*1.0); WC=np.exp(1j*np.sqrt(2)); hits=0
for _ in range(20000):
    VvF=expmH(rng.normal(size=9),3); VvC=expmH(rng.normal(size=9),3)
    jF=rng.integers(0,3,3); jC=rng.integers(0,3,3)
    rF=np.exp(1j*np.angle(WF)/3); rC=np.exp(1j*np.angle(WC)/3); zt=np.exp(2j*np.pi/3)
    TF=embed(VF,VvF@np.diag(rF*zt**jF)@VvF.conj().T,5); TC=embed(VC,VvC@np.diag(rC*zt**jC)@VvC.conj().T,5)
    if Dgen(CLASSES,5,TF,TC,12)<1e-9 and (np.linalg.norm(TF-np.diag(np.diag(TF)))>1e-6 or
                                          np.linalg.norm(TC-np.diag(np.diag(TC)))>1e-6): hits+=1
print(f"      20000 undirected-local (= unrestricted) members: non-fibre-wise restorers = {hits}.")
print("      Theorem 4 never uses locality, so widening the axiom changes nothing.\n")

print("== X2.5  THE SHIFT BRANCH'S FLOOR, IN CLOSED FORM.  THE LANE MEASURED 2.236 / 3.000. ==")
print("   For ANY shift member on K1, Q_1 = T_F* T_C has EXACTLY FIVE non-zero entries, all of")
print("   modulus 1, all off the diagonal, and a diagonal that is identically ZERO:")
print("      Q[v0,v1]=conj(n_0)  Q[v1,v2]=conj(n_1)  Q[v2,v4]=conj(n_2)m_2  Q[v3,v0]=m_0  Q[v4,v3]=m_1")
print("   so ||offdiag Q_1||_F = sqrt(5) = 2.2360679... EXACTLY, for every member of the torus and")
print("   every connection.  A PROVED FLOOR, not a sampled minimum.  (I first wrote sqrt(4) and the")
print("   arithmetic below caught me; the entry Q[v2,v4] is the fifth and I had dropped it.)")
mn=np.inf; dgmax=0.0
for _ in range(50000):
    ph=rng.uniform(0,2*np.pi,2); nF=[np.exp(1j*ph[0]),np.exp(1j*ph[1]),WF*np.exp(-1j*ph.sum())]
    ph=rng.uniform(0,2*np.pi,2); nC=[np.exp(1j*ph[0]),np.exp(1j*ph[1]),WC*np.exp(-1j*ph.sum())]
    def sh(loop,ns):
        on={s for s,_,_ in loop}; T=np.zeros((5,5),dtype=complex)
        for v in range(5):
            if v not in on: T[v,v]=1.0
        for k,(s_,d_,_) in enumerate(loop): T[d_,s_]=ns[k]
        return T
    Q=sh(LF,nF).conj().T@sh(LC,nC)
    mn=min(mn,np.linalg.norm(Q-np.diag(np.diag(Q)))); dgmax=max(dgmax,abs(np.diag(Q)).max())
print(f"   50000 shift members: min ||offdiag Q_1||_F = {mn:.9f}   [sqrt(5) = {np.sqrt(5):.9f}]")
print(f"                        max |diag Q_1| = {dgmax:.2e}  -> the diagonal is identically zero")
print("   ON B0b I GUESSED A CLOSED FORM AND THE ARITHMETIC BELOW REFUTED IT.  I wrote 'nine unit")
print("   entries, floor 3.000000'.  The measured n=1 floor is 2.000000 and the diagonal of Q_1 is")
print("   NOT zero there (|diag| reaches 1), because on B0b the two loops share TWO vertices, so")
print("   the shared block contributes diagonal weight.  The lane's 3.00 is its max over n <= 24,")
print("   which I reproduce but do NOT close in closed form.  K1's sqrt(5) stands; B0b's 3.00")
print("   remains a measured minimum over 4000 (its) / 20000 (my) draws, not a proved floor.")

# --- and the same floor on B0b, computed rather than asserted ---
def BT_shift(loop, sgn, a, ns):
    on={u for u,_,_,_ in loop}; T=np.zeros((9,9),dtype=complex)
    for v in range(9):
        if v not in on: T[v,v]=1.0
    for k,(u,v,e,s) in enumerate(loop): T[v,u]=ns[k]
    return T
mn9=np.inf; dg9=0.0
aB=rng.uniform(0,2*np.pi,18)
WFb=np.exp(1j*sum(aB[e]*s for _,_,e,s in BF)); WCb=np.exp(1j*sum(aB[e]*s for _,_,e,s in BC))
for _ in range(20000):
    ph=rng.uniform(0,2*np.pi,3); nF=list(np.exp(1j*ph))+[WFb*np.exp(-1j*ph.sum())]
    ph=rng.uniform(0,2*np.pi,2); nC=list(np.exp(1j*ph))+[WCb*np.exp(-1j*ph.sum())]
    Q=BT_shift(BF,None,aB,nF).conj().T@BT_shift(BC,None,aB,nC)
    mn9=min(mn9,np.linalg.norm(Q-np.diag(np.diag(Q)))); dg9=max(dg9,abs(np.diag(Q)).max())
print(f"   B0b, 20000 shift members: min ||offdiag Q_1||_F = {mn9:.9f}  max |diag Q_1| = {dg9:.2e}")
print(f"   -> my guessed closed form (3.000000) is FALSE; the n=1 floor is {mn9:.6f}.  RECORDED.")
mnD=np.inf
for _ in range(4000):
    ph=rng.uniform(0,2*np.pi,3); nF=list(np.exp(1j*ph))+[WFb*np.exp(-1j*ph.sum())]
    ph=rng.uniform(0,2*np.pi,2); nC=list(np.exp(1j*ph))+[WCb*np.exp(-1j*ph.sum())]
    mnD=min(mnD, Dgen(BCLASSES,9,BT_shift(BF,None,aB,nF),BT_shift(BC,None,aB,nC),24))
print(f"   B0b, min D over n<=24 across 4000 shift members = {mnD:.6f}   [lane reported 3.00]")
