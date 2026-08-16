# W-07 leg C — (1) the control S3's test provides; (2) how deep the near-return goes off the
# published connection; (3) the 52-partition sweep W-06 offers as its uniqueness result.
import numpy as np, itertools
rng=np.random.default_rng(20260816)
EDGES=[(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)]; FACE_V={0,1,2}; CYC_V={0,3,4}
TREE={1:(0,),2:(0,1),3:(3,),4:(3,4)}
def dress(s,a):
    u=np.exp(1j*np.asarray(a)); t=np.array(s,dtype=complex)
    for v,p in TREE.items():
        w=1.0+0j
        for e in p: w*=u[e]
        t[v]=s[v]/w
    return t

print("== C1  THE CONTROL: what S3's own (undressed diagonal) test returns on the same branches ==")
a=np.array([np.pi/3]*3+[np.pi/2]*3); s=rng.normal(size=5)+1j*rng.normal(size=5); s/=np.linalg.norm(s)
WF=np.exp(1j*sum(a[:3])); WC=np.exp(1j*sum(a[3:]))
x=np.array([ (WF if v in FACE_V else 1)*s[v] for v in range(5)])
y=np.array([ (WC if v in CYC_V  else 1)*s[v] for v in range(5)])
print(f"  max_v | |x_v|^2 - |y_v|^2 |  (every diagonal / S3-gauge-invariant observable) = {np.abs(np.abs(x)**2-np.abs(y)**2).max():.3e}")
t_x,t_y=dress(x,a),dress(y,a)
print(f"  | conj(t_2) t_3 [x] - conj(t_2) t_3 [y] |  (the DRESSED observable)            = {abs(np.conj(t_x[2])*t_x[3]-np.conj(t_y[2])*t_y[3]):.12f}")
print("  -> W-06's core claim stands independently: the dressed algebra SEES what S3's test cannot.\n")

print("== C2  HOW DEEP DOES THE NEAR-RETURN GO?  amp*|rho^k - 1|, rho = conj(W_F)/W_C ==")
print(f"  {'K':>9} {'PUBLISHED (rho=-i)':>26} {'GENERIC (rho=e^{2pi i/phi^2})':>32}")
phi=(1+5**0.5)/2
for K in [10**3,10**4,10**5,10**6,10**7]:
    k=np.arange(1,K+1)
    dpub=np.abs((-1j)**k-1); dgen=np.abs(np.exp(2j*np.pi*(1/phi**2)*k)-1)
    print(f"  {K:>9} {dpub.min():>26.3e} {dgen.min():>32.3e}   (zeros: {int((dpub<1e-12).sum())} vs {int((dgen<1e-12).sum())})")
print("  -> published: EXACTLY zero on K/4 cells at every K.  generic: never zero, floor falls like ~1/K.\n")

print("== C3  THE 52-PARTITION SWEEP — re-run, and its degrees of freedom counted ==")
def partitions(c):
    if len(c)==1: yield [c]; return
    first,rest=c[0],c[1:]
    for p in partitions(rest):
        for i in range(len(p)): yield p[:i]+[[first]+p[i]]+p[i+1:]
        yield [[first]]+p
d1=np.zeros((5,6))
for j,(s_,t_) in enumerate(EDGES): d1[t_,j]+=1; d1[s_,j]-=1
from collections import Counter
cnt=Counter(); winners=[]
for P in partitions(list(range(5))):
    B=np.zeros((5,len(P)))
    for i,blk in enumerate(P):
        for v in blk: B[v,i]=1
    inv=6-np.linalg.matrix_rank((d1.T@B))          # invariants = E - rank(action of the subgroup)
    cnt[inv]+=1
    if inv==2: winners.append([sorted(b) for b in P])
print(f"  partitions swept: {sum(cnt.values())}   (Bell(5) = 52)")
print(f"  invariant-count distribution: {dict(sorted(cnt.items()))}")
print(f"  partitions giving 2 = b1 + #faces : {len(winners)}  ->  {winners}")
print(f"  one-line reason: a k-block subgroup acts through <= k-1 parameters (constants act trivially),")
print(f"  so rank <= k-1; rank 4 forces k = 5. THE SWEEP COULD NOT HAVE RETURNED ANY OTHER ANSWER.")
