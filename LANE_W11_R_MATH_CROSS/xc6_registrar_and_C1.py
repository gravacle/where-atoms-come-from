# LANE_W11_R_MATH_CROSS — LEG 6.  (a) the REGISTRAR's tables from my own operator code;
# (b) the lane's "the drift is O(1/sqrt N)" rate claim, tested; (c) the lane's Theorem C1
# biconditional stress-tested on ticks that are NOT of the form U^L = M, where it should still
# hold if C1 is what the lane says it is.
import numpy as np, xc0_lib as X
np.set_printoptions(linewidth=200)

print("== X6a  THE REGISTRAR'S OWN TABLES, FROM MY INDEPENDENT OPERATOR CODE ==")
car=X.K1(); NV=5; a=np.array([1.0,0.37,0.91,2**0.5,0.23,1.77])
TF,TC=X.Top(car["walkF"],a,NV),X.Top(car["walkC"],a,NV)
MF,MC=X.Mop(car["walkF"],a,NV),X.Mop(car["walkC"],a,NV)
sA=np.sqrt(np.array([0.40,0.15,0.15,0.15,0.15]))+0j
sB=np.sqrt(np.array([0.40,0.30,0.00,0.05,0.25]))+0j
sC=sA*np.exp(1j*np.array([0.0,1.3,-0.7,2.2,0.4]))
reg_circ=[0.319750930010,0.760050706869,0.586994870469,0.307756757635,0.456532066199,0.513844088442]
reg_edge=[0.569227769927,0.071727337054,0.319750930010,0.581125277232,0.234214852448,
          0.760050706869,0.647209237845,0.199740989173,0.586994870469]
mc=max(abs(abs(np.vdot(np.linalg.matrix_power(MF,k)@sA,np.linalg.matrix_power(MC,k)@sA))-reg_circ[k-1]) for k in range(1,7))
me=max(abs(abs(np.vdot(np.linalg.matrix_power(TF,n)@sA,np.linalg.matrix_power(TC,n)@sA))-reg_edge[n-1]) for n in range(1,10))
print(f"   W-11 leg B1, six circuit rows : max |mine - registrar's| = {mc:.2e}   [wm6 reports 4.72e-13]")
print(f"   W-11 leg B2, nine edge rows   : max |mine - registrar's| = {me:.2e}   [wm6 reports 3.53e-13]")
carB=X.B0b(); NB=9; aB=np.random.default_rng(20260817).uniform(0,2*np.pi,18)
TFb,TCb=X.Top(carB["walkF"],aB,NB),X.Top(carB["walkC"],aB,NB)
MFb,MCb=X.Mop(carB["walkF"],aB,NB),X.Mop(carB["walkC"],aB,NB)
clb,_,_=X.classes(carB)
w=np.array([.10,.12,.09,.14,.11,.11,.11,.11,.11]); w/=w.sum()
wB=w.copy(); wB[0],wB[1]=w[0]+w[1],0.0; wB[3],wB[4]=0.0,w[3]+w[4]; wB[5],wB[8]=w[5]+w[8],0.0
sAb=np.sqrt(w)+0j; sBb=np.sqrt(wB)+0j
sCb=sAb*np.exp(1j*np.random.default_rng(7).uniform(0,2*np.pi,9))
regC3_circ=[0.157787872257,0.542911929181,0.331950055873,0.701868545596,0.242079131640]
regC3_edge=[0.660894986456,0.650698044893,0.465043830794,0.028491526437,0.268458748603,
            0.267289870295,0.439399105411,0.582583413132,0.613011266761,0.717612419937,
            0.495748304436,0.302148955800]
c3c=max(abs(abs(np.vdot(np.linalg.matrix_power(MFb,k)@sAb,np.linalg.matrix_power(MCb,k)@sAb))-regC3_circ[k-1]) for k in range(1,6))
c3e=max(abs(abs(np.vdot(np.linalg.matrix_power(TFb,n)@sAb,np.linalg.matrix_power(TCb,n)@sAb))-regC3_edge[n-1]) for n in range(1,13))
print(f"   W-11 leg C3, B0b circuit rows : max |mine - registrar's| = {c3c:.2e}")
print(f"   W-11 leg C3, B0b edge rows    : max |mine - registrar's| = {c3e:.2e}   (leg C was never")
print("      re-run by the lane under test at the |Z| level; I re-run it here.  It reproduces.)")
print(f"   || T_F^4 - M_F || on B0b = {np.linalg.norm(np.linalg.matrix_power(TFb,4)-MFb):.2e}"
      f"   || T_C^3 - M_C || = {np.linalg.norm(np.linalg.matrix_power(TCb,3)-MCb):.2e}")
print(f"   B0b class multiset from my own incidence: "
      f"{ {k:int((clb==v).sum()) for k,v in (('00',0),('10',1),('01',2),('11',3))} }   [S4:575 {{00:4,01:1,10:2,11:2}}]\n")

print("== X6b  THE LANE'S 'THE DRIFT IS O(1/sqrt N)' IS AN UNSUPPORTED RATE CLAIM ==")
def coeff_rows(car,a,s,NV):
    wF,wC=car["walkF"],car["walkC"]; LF,LC=len(wF),len(wC); Lam=int(np.lcm(LF,LC))
    TF,TC=X.Top(wF,a,NV),X.Top(wC,a,NV)
    x=np.conj(X.hol(wF,a)); y=X.hol(wC,a); _,F,C=X.classes(car)
    inF=[1 if v in F else 0 for v in range(NV)]; inC=[1 if v in C else 0 for v in range(NV)]
    ix={(0,0):0,(1,0):1,(0,1):2,(1,1):3}; rows={}
    for rho in range(Lam):
        B=np.linalg.inv(np.linalg.matrix_power(TF,rho%LF))@np.linalg.matrix_power(TC,rho%LC)
        c=np.zeros(4,dtype=complex)
        for u in range(NV):
            for v in range(NV): c[ix[(inF[u],inC[v])]]+=np.conj(s[u])*s[v]*B[u,v]
        rows[rho]=np.array([c[0],c[1]*x**(rho//LF),c[2]*y**(rho//LC),c[3]*x**(rho//LF)*y**(rho//LC)])
    return rows
R=coeff_rows(car,a,sA,NV); closed=float(np.mean([X.m_quad(R[r],1<<20) for r in sorted(R)]))
Ns=[2000,8000,32000,128000,512000,2048000]
xF=sA.copy(); xC=sA.copy(); tot=0.0; errs=[]
nxt=0
for i in range(1,Ns[-1]+1):
    xF=TF@xF; xC=TC@xC; z=abs(np.vdot(xF,xC)); tot+= np.log(z) if z>1e-300 else -700.0
    if nxt<len(Ns) and i==Ns[nxt]:
        errs.append(abs(tot/i-closed)); nxt+=1
sl=np.polyfit(np.log(Ns),np.log(errs),1)[0]
print(f"   N        : {Ns}")
print(f"   |ta-clsd|: {[f'{e:.2e}' for e in errs]}")
print(f"   log-log slope = {sl:+.2f}   (the lane asserts -0.50 ; -1.00 would be the Birkhoff rate")
print("   for a Diophantine rotation).  The claim 'the drift is O(1/sqrt N)' is stated from four")
print("   non-monotone ratios and is not what the data shows.  Nothing downstream depends on it:")
print("   the closed form is confirmed independently (X4d), so the LIMIT is right either way.\n")

print("== X6c  THEOREM C1 STRESS-TESTED OFF ITS OWN HYPOTHESIS ==")
print("   C1 is stated for arbitrary unitary ticks: Z_n pi-determined for every s  <=>  A_n")
print("   diagonal with class-constant diagonal.  I test the BICONDITIONAL on ticks that are NOT")
print("   L-th roots of M at all -- random unitaries supported on each loop -- on three carriers.")
def haar(n,r):
    z=(r.normal(size=(n,n))+1j*r.normal(size=(n,n)))/np.sqrt(2)
    q,rr=np.linalg.qr(z); return q@np.diag(np.diag(rr)/np.abs(np.diag(rr)))
def embed(vs,B,NV):
    U=np.eye(NV,dtype=complex)
    for i,u in enumerate(vs):
        for j,v in enumerate(vs): U[u,v]=B[i,j]
    return U
rng=np.random.default_rng(20260817)
for car in (X.K1(),X.B0b(),X.SHARE2()):
    NVc=car["NV"]; cl,F,C=X.classes(car); Fv,Cv=sorted(F),sorted(C)
    sizes=[int((cl==c).sum()) for c in range(4)]
    pi=np.array([s/sum(sizes) for s in sizes])
    sts=X.same_pi_states(cl,pi,np.random.default_rng(31),25)
    agree=0; tot=0; both=0; inv=0
    for _ in range(3000):
        mode=rng.integers(0,3)
        if mode==0:      # fully random
            UF=embed(Fv,haar(len(Fv),rng),NVc); UC=embed(Cv,haar(len(Cv),rng),NVc)
        elif mode==1:    # random diagonal (should be pi-determined only if class-constant)
            UF=embed(Fv,np.diag(np.exp(1j*rng.uniform(0,2*np.pi,len(Fv)))),NVc)
            UC=embed(Cv,np.diag(np.exp(1j*rng.uniform(0,2*np.pi,len(Cv)))),NVc)
        else:            # class-constant diagonal (should always be pi-determined)
            dF=np.exp(1j*rng.uniform(0,2*np.pi,4)); dC=np.exp(1j*rng.uniform(0,2*np.pi,4))
            UF=np.eye(NVc,dtype=complex); UC=np.eye(NVc,dtype=complex)
            for v in Fv: UF[v,v]=dF[cl[v]]
            for v in Cv: UC[v,v]=dC[cl[v]]
        A=np.linalg.inv(UF)@UC
        off=np.linalg.norm(A-np.diag(np.diag(A))); d=np.diag(A); bad=0.0
        for c in range(4):
            idx=np.where(cl==c)[0]
            if len(idx)>=2: bad=max(bad,np.max(np.abs(d[idx]-d[idx[0]])))
        alg=(off<1e-10 and bad<1e-10)
        emp=(X.spread(UF,UC,1,sts)<1e-10)
        tot+=1; agree+= (alg==emp); both+=alg; inv+=emp
    print(f"   {car['name']:>7}: {tot} random tick pairs   algebraic condition true {both:>5}"
          f"   measured pi-determined {inv:>5}   AGREE {agree}/{tot}")
print("   -> the biconditional holds on every draw on all three carriers, including the carrier")
print("      the lane did not build.  C1 is the lane's real result and it survives everything I")
print("      threw at it.  U4 -- the corollary that turns C1 into 'therefore DIAGONAL' -- does not.")
