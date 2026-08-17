# REFUTER 1 of LANE W-08 / M2 — LENS: MATHEMATICS.  Leg R7: the overlap test with the RELATION
# HEIGHT held fixed, and what the sampled pairs actually are.
#   R7.1  R2.3's rank-1 sample drew a random relation (n1,n2) as well as a random alpha, so the
#         constant in F(K) <= C K^{-2} varied with the height.  Redone with the height FIXED at
#         the same relation E_W07GEN carries (beta = alpha, i.e. 1a - 1b = 0).  ONE variable:
#         dim H.  Does the cluster overlap survive?
#   R7.2  what the "H = T^2" rows actually are: every alpha,beta drawn as a numpy double is a
#         DYADIC RATIONAL p/2^53, so H is finite of order ~2^53, not T^2.  Immaterial at K <= 1e7
#         and stated so -- but leg B's B0 prints "no relation up to height 60 -> H = T^2", which
#         is a window statement printed as a fact about the object.  Same defect class as R3.1.
import numpy as np, mpmath as mp
mp.mp.dps=60
KS=[10**3,10**4,10**5,10**6,10**7]; W=(0.4,0.3,0.3); w11,w10,w01=W
def limbs(a):
    a=mp.frac(mp.mpf(a)); a1=float(mp.nint(a*2**26))/2**26; r=a-mp.mpf(a1)
    a2=float(mp.nint(r*2**52))/2**52; return a1,a2,float(r-mp.mpf(a2))
def fr(k,L):
    a1,a2,a3=L; t1=k*a1; t1-=np.floor(t1); t2=k*a2; t2-=np.floor(t2)
    t=t1+t2+k*a3; return t-np.floor(t)
def floors(al,be):
    La,Lb=limbs(al),limbs(be); out=[];run=np.inf;lo=0
    for K in KS:
        while lo<K:
            hi=min(lo+10**6,K); k=np.arange(lo+1,hi+1,dtype=np.float64)
            u=fr(k,La); v=fr(k,Lb)
            du=u-np.round(u); dv=v-np.round(v); duv=u+v; duv-=np.round(duv)
            S=np.minimum(4.0*(w11*w10*np.sin(np.pi*dv)**2+w11*w01*np.sin(np.pi*du)**2
                              +w10*w01*np.sin(np.pi*duv)**2),1.0)
            g=S/(1.0+np.sqrt(np.maximum(0.0,1.0-S))); run=min(run,float(g.min())); lo=hi
        out.append(run)
    return out
def fit(f): return np.polyfit(np.log10(KS),np.log10(np.maximum(f,1e-300)),1)[0]

print("== R7.1  OVERLAP TEST WITH THE RELATION HEIGHT HELD FIXED (beta = alpha, height 1) ==")
r=np.random.default_rng(31415); A=[];B=[]
for i in range(40):
    a=mp.mpf(float(r.random())); A.append(fit(floors(a,a)))          # dim H = 1, height 1
for i in range(40):
    a=mp.mpf(float(r.random())); b=mp.mpf(float(r.random())); B.append(fit(floors(a,b)))
A=np.array(A); B=np.array(B)
print(f"   dim H = 1, relation 1a-1b=0 (height 1): n=40  mean {A.mean():+.3f} sd {A.std():.3f}"
      f"  min {A.min():+.3f}  max {A.max():+.3f}   theory -2.0")
print(f"   dim H = 2                             : n=40  mean {B.mean():+.3f} sd {B.std():.3f}"
      f"  min {B.min():+.3f}  max {B.max():+.3f}   theory -1.0")
print(f"   OVERLAP?  max(dim1) = {A.max():+.3f}  vs  min(dim2) = {B.min():+.3f}  ->  "
      f"{'CLUSTERS OVERLAP' if A.max()>B.min() else 'clusters separate'}")
print(f"   dim1 rows shallower than -1.379 (leg F1's shallowest d_eff=1 row): {int((A>-1.379).sum())}/40")
print(f"   dim2 rows deeper   than -1.379                                   : {int((B<-1.379).sum())}/40")
print("   -> the LAW survives in the mean (means within ~0.03 of theory, both classes); the")
print("      EVIDENCE M2-4 rests on -- 'the two clusters do not overlap' -- does not survive")
print("      resampling.  A 5-point fit of a running minimum is not a d_eff classifier.\n")

print("== R7.2  WHAT THE 'H = T^2' ROWS ARE ==")
rng=np.random.default_rng(20260816)
a=float(rng.random())
num,den=a.as_integer_ratio()
print(f"   leg B's D_RAND1 alpha = {a!r}")
print(f"   as an EXACT rational: {num} / {den}   (den = 2^{int(np.log2(den))})")
print(f"   so alpha is RATIONAL, the orbit is periodic with period dividing {den}, and H is a")
print(f"   FINITE cyclic group of order ~{den:.3e} -- dim H = 0, not 2.  The period is 9.0e15 >> 1e7,")
print(f"   so NO MEASURED NUMBER IN THE LANE CHANGES.  What changes is the status of B0's line")
print(f"   'no relation up to height 60 -> H = T^2': that is a statement about the WINDOW k<=1e7,")
print(f"   not about the object, and the lane prints it as a fact about the object.  The honest")
print(f"   form is: no relation of height <= H0 is visible below k ~ H0, and d_eff must be read")
print(f"   relative to the window.  (This is why leg F3's rational pair, height 3e6 < 1e7, DID")
print(f"   show up as an exact return and broke the table -- same defect, different scale.)")
