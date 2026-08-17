# REFUTER 1 of LANE W-08 / M2 — LENS: MATHEMATICS.  Leg R2: THE SCALING CLAIM F(K) ~ K^{-2/d_eff}.
#   R2.1  reproduce leg F1's table with M2's OWN naive residues, and again with an exact-residue
#         protocol, and refit.  (Leg F1 does NOT use the 3-limb protocol its own conventions page
#         publishes; leg B does.  Establish whether that changes the fitted exponents.)
#   R2.2  THE ROW LEG F1 LEAVES OUT.  F_VWA is dim H = 2 by M2's own B0 classification.  Put it in
#         the same table with the same fit and read the exponent.
#   R2.3  IS THE CLUSTER SEPARATION REAL?  40 fresh d_eff=2 pairs and 40 fresh d_eff=1 pairs,
#         same 5-point fit, same K-grid.  Do the two clusters overlap?
# PRECISION.  Residues by a 2-limb split of a 60-digit mpmath value (k*a_hi and k*a_mid exact for
# k < 2^27); every reported argmin re-checked in EXACT big-integer arithmetic against a 2^-120
# rational truncation of the constant.  Double precision otherwise, and said so.
import numpy as np, mpmath as mp
mp.mp.dps = 60
KS=[10**3,10**4,10**5,10**6,10**7]

def limbs(a):
    a=mp.frac(mp.mpf(a))
    a1=float(mp.nint(a*2**26))/2**26
    r=a-mp.mpf(a1)
    a2=float(mp.nint(r*2**52))/2**52
    a3=float(r-mp.mpf(a2))
    return a1,a2,a3
def frac_exact(k,L):
    a1,a2,a3=L
    t1=k*a1; t1-=np.floor(t1)
    t2=k*a2; t2-=np.floor(t2)
    t=t1+t2+k*a3
    return t-np.floor(t)
def frac_naive(k,a):                      # exactly what leg F1's floors() does
    return (k*float(a))%1.0

def gap(u,v,W):
    w11,w10,w01=W
    du=u-np.round(u); dv=v-np.round(v); duv=u+v; duv-=np.round(duv)
    S=np.minimum(4.0*(w11*w10*np.sin(np.pi*dv)**2+w11*w01*np.sin(np.pi*du)**2
                      +w10*w01*np.sin(np.pi*duv)**2),1.0)
    return S/(1.0+np.sqrt(np.maximum(0.0,1.0-S)))

def floors(al,be,W,mode="exact",KS=KS,chunk=10**6,want_arg=False):
    """running min of 1-|Z_k| at the checkpoints.  mode='exact' -> limb split; 'naive' -> leg F1."""
    La,Lb=limbs(al),limbs(be); fa,fb=float(mp.frac(mp.mpf(al))),float(mp.frac(mp.mpf(be)))
    out=[]; run=np.inf; arg=-1; lo=0
    for K in KS:
        while lo<K:
            hi=min(lo+chunk,K); k=np.arange(lo+1,hi+1,dtype=np.float64)
            if mode=="exact": u,v=frac_exact(k,La),frac_exact(k,Lb)
            else:             u,v=frac_naive(k,fa),frac_naive(k,fb)
            g=gap(u,v,W); j=int(np.argmin(g))
            if g[j]<run: run=float(g[j]); arg=lo+1+j
            lo=hi
        out.append(run)
    return (out,arg) if want_arg else out

def exact_check(al,be,W,k):
    """1-|Z_k| in 60-digit arithmetic with EXACT big-integer reduction of k*alpha mod 1."""
    w11,w10,w01=W; res=[]
    for a in (al,be):
        P=int(mp.floor(mp.frac(mp.mpf(a))*mp.mpf(2)**120))       # alpha ~ P/2^120, exact rational
        r=(k*P)%(2**120); res.append(mp.mpf(r)/mp.mpf(2)**120)
    u,v=res
    du=u-mp.nint(u); dv=v-mp.nint(v); duv=u+v-mp.nint(u+v)
    S=4*(w11*w10*mp.sin(mp.pi*dv)**2+w11*w01*mp.sin(mp.pi*du)**2+w10*w01*mp.sin(mp.pi*duv)**2)
    return S/(1+mp.sqrt(max(mp.mpf(0),1-S)))

phi=(1+mp.sqrt(5))/2; tt=2*mp.cos(2*mp.pi/7)
Q=2999999
vwa_a=mp.mpf(1234567)/Q+mp.mpf('1e-18')*(mp.sqrt(2)-1)     # leg B's OWN F_VWA pair, verbatim
vwa_b=mp.mpf( 765431)/Q+mp.mpf('1e-18')*(mp.sqrt(3)-1)
rng=np.random.default_rng(20260816); RND=[(mp.mpf(float(rng.random())),mp.mpf(float(rng.random()))) for _ in range(3)]
CONN=[("B_S3RES  f=2.0,c=11/10",1/mp.pi,mp.mpf(11)/(20*mp.pi),1),
      ("E_W07GEN 2pi.phi,2pi.phi^2",mp.frac(phi),mp.frac(phi**2),1),
      ("C_BADAPP cubic pair",mp.frac(tt),mp.frac(tt**2),2),
      ("D_RAND1",RND[0][0],RND[0][1],2),
      ("D_RAND2",RND[1][0],RND[1][1],2),
      ("D_RAND3",RND[2][0],RND[2][1],2),
      ("F_VWA  <- LEG F1 OMITS THIS ROW",vwa_a,vwa_b,2)]
RS=[("RS-G",(0.4,0.3,0.3)),("RS-P",(0.5,0.0,0.5))]

print("== R2.1 / R2.2  LEG F1's TABLE, REPRODUCED AND THEN COMPLETED ==")
print("   'naive' = leg F1's own (k*alpha)%1.0 in float64.  'exact' = the 3-limb protocol that")
print("   this lane's PUBLISHED_CONVENTIONS promises and leg B uses.  d_eff: RS-G -> dim H;")
print("   RS-P -> 1 (p10 = 0 kills one character).")
print(f"   {'connection':<32}{'ready':<6}{'d_eff':>5}" + "".join(f"{f'F(1e{e})':>12}" for e in range(3,8))
      + f"{'fit(exact)':>11}{'fit(naive)':>11}{'theory':>8}")
FIT={}
for lab,al,be,dim in CONN:
    for rt,W in RS:
        deff=dim if rt=="RS-G" else 1
        fe,arg=floors(al,be,W,"exact",want_arg=True); fn=floors(al,be,W,"naive")
        se=np.polyfit(np.log10(KS),np.log10(np.maximum(fe,1e-300)),1)[0]
        sn=np.polyfit(np.log10(KS),np.log10(np.maximum(fn,1e-300)),1)[0]
        FIT[(lab[:8],rt)]=(se,sn,fe,arg,deff)
        print(f"   {lab:<32}{rt:<6}{deff:>5}"+"".join(f"{x:>12.3e}" for x in fe)
              +f"{se:>11.3f}{sn:>11.3f}{-2.0/deff:>8.1f}")
print()
print("   EXACT-ARITHMETIC RE-CHECK of every reported argmin (big-int k*alpha mod 2^120):")
for (lab,rt),(se,sn,fe,arg,deff) in FIT.items():
    al,be=[(a,b) for l,a,b,_ in CONN if l[:8]==lab][0]
    W=dict(RS)[rt]
    ex=exact_check(al,be,W,arg)
    rel=abs(mp.mpf(fe[-1])-ex)/ex if ex>0 else mp.mpf(0)
    print(f"     {lab:<9}{rt:<6} k*={arg:>8}  float64 F(1e7)={fe[-1]:.6e}  exact={mp.nstr(ex,7):>12}"
          f"  rel.err={float(rel):.2e}")
print()
print("   THE OMITTED ROW, READ OFF:  F_VWA has dim H = 2 (leg B's own B0 block says so, no")
print("   relation to height 60) and its 5-point fit is printed above.  Compare the d_eff = 2")
print("   cluster leg F1 reports: -0.970, -1.027, -1.080, -0.845.\n")

print("== R2.3  IS THE 'TWO CLUSTERS DO NOT OVERLAP' EVIDENCE ROBUST?  FRESH SAMPLES. ==")
print("   40 rank-2 pairs (alpha,beta ~ U[0,1)) and 40 rank-1 pairs (beta = (n1*alpha - m)/n2,")
print("   n1,n2 coprime <= 12, so 'n1 alpha - n2 beta = m' and dim H = 1), IDENTICAL fit, K-grid,")
print("   ready state RS-G, code path.  Seed 424242.  ONE VARIABLE MOVES: dim H.")
r=np.random.default_rng(424242); W=(0.4,0.3,0.3)
S2=[];S1=[]
for i in range(40):
    a=mp.mpf(float(r.random())); b=mp.mpf(float(r.random()))
    S2.append(np.polyfit(np.log10(KS),np.log10(np.maximum(floors(a,b,W),1e-300)),1)[0])
for i in range(40):
    a=mp.mpf(float(r.random()))
    while True:
        n1=int(r.integers(1,13)); n2=int(r.integers(1,13))
        if np.gcd(n1,n2)==1: break
    m=int(r.integers(0,4))
    b=mp.frac((n1*a-m)/n2)
    S1.append(np.polyfit(np.log10(KS),np.log10(np.maximum(floors(a,b,W),1e-300)),1)[0])
S1=np.array(S1);S2=np.array(S2)
print(f"   dim H = 2 (theory -1.0): n=40  mean {S2.mean():+.3f}  sd {S2.std():.3f}  "
      f"min {S2.min():+.3f}  max {S2.max():+.3f}")
print(f"   dim H = 1 (theory -2.0): n=40  mean {S1.mean():+.3f}  sd {S1.std():.3f}  "
      f"min {S1.min():+.3f}  max {S1.max():+.3f}")
ov = S1.max() > S2.min()
print(f"   OVERLAP?  max(dimH=1 fits) = {S1.max():+.3f} vs min(dimH=2 fits) = {S2.min():+.3f}"
      f"   ->  {'CLUSTERS OVERLAP' if ov else 'clusters separate'}")
print(f"   dim H=1 fits shallower than leg F1's shallowest d_eff=1 row (-1.379): "
      f"{int((S1>-1.379).sum())}/40")
print(f"   dim H=1 fits shallower than leg F1's deepest d_eff=2 row (-1.080): "
      f"{int((S1>-1.080).sum())}/40")
print(f"   dim H=2 fits deeper than -1.379: {int((S2<-1.379).sum())}/40")
mis=int((S1>-1.5).sum())+int((S2<-1.5).sum())
print(f"   misclassification rate of a -1.5 threshold applied to the fitted exponent: {mis}/80")
np.save("r2_cluster_fits.npy",np.vstack([S2,S1]))
print("   (fits saved to r2_cluster_fits.npy)")
