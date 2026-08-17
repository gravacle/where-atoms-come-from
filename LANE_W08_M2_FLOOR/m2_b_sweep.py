# LANE W-08 / M2 leg B — THE SWEEP.  F(K) = min_{k<=K} (1-|Z_k|), the full distribution of
# (1-|Z_k|), the partial sums of (1-|Z_k|), and sum log|Z_k|, on every connection class.
# READY STATE HELD FIXED WITHIN EACH TABLE.  ONE VARIABLE MOVES: the connection's (alpha,beta).
#
# PRECISION PROTOCOL (double is the default; this is where it would have bitten):
#  * 1-|Z_k| is NEVER computed as 1 - |Z_k|.  It is computed from the EXACT identity of leg A,
#      1-|Z|^2 = 4[ w11 w10 sin^2(pi v) + w11 w01 sin^2(pi u) + w10 w01 sin^2(pi(u+v)) ],
#    each sin^2 evaluated at the NEAREST-INTEGER residue, so relative accuracy survives to the
#    smallest near-return.  Then 1-|Z| = (1-|Z|^2)/(1+|Z|).
#  * k*alpha mod 1 uses a 3-limb split of a 60-digit mpmath value of alpha: alpha = a1+a2+a3 with
#    a1 a multiple of 2^-26 and a2 a multiple of 2^-52.  k*a1 and k*a2 are EXACT products for
#    k < 2^27, their fractional parts are exact, and their sum is exactly representable.  Residual
#    error in u is ~1.1e-16 ABSOLUTE, set by the last addition of k*a3.  Near-returns are therefore
#    resolvable to |u| ~ 1e-15, i.e. 1-|Z| ~ 1e-29.  Everything reported here is >= 1e-16.
#  * Rational classes are done in EXACT INTEGER arithmetic instead (k*num mod den).
import numpy as np, mpmath as mp, pickle
mp.mp.dps = 60
KMAX = 10**7
BOUNDS = [0,10**3,10**4,10**5,10**6] + [i*10**6 for i in range(2,11)]
CHECKPOINTS = [10**3,10**4,10**5,10**6,10**7]
EPS = [10.0**-e for e in range(1,15)]

def limbs(a):
    a1 = float(mp.nint(a*2**26))/2**26
    r  = a - mp.mpf(a1)
    a2 = float(mp.nint(r*2**52))/2**52
    a3 = float(r - mp.mpf(a2))
    return a1,a2,a3

def frac_irr(k, L):
    a1,a2,a3 = L
    t1 = k*a1; t1 -= np.floor(t1)          # exact
    t2 = k*a2; t2 -= np.floor(t2)          # exact
    t  = t1 + t2 + k*a3                    # only this line rounds (~1.1e-16)
    return t - np.floor(t)

phi = (1+mp.sqrt(5))/2
tt  = 2*mp.cos(2*mp.pi/7)                  # t^3 + t^2 - 2t - 1 = 0; Q(t) totally real cubic
rng = np.random.default_rng(20260816)
RND = [(mp.mpf(float(rng.random())), mp.mpf(float(rng.random()))) for _ in range(3)]
Q = 2999999                                # very-well-approximable rank-2 pair: rational core of
vwa_a = mp.mpf(1234567)/Q + mp.mpf('1e-18')*(mp.sqrt(2)-1)   # denominator Q, perturbed 1e-18,
vwa_b = mp.mpf( 765431)/Q + mp.mpf('1e-18')*(mp.sqrt(3)-1)   # i.e. BELOW the double ulp of alpha.

CLASSES = [
 ("A_S1PUB   f=pi, c=3pi/2         (S1 published)",     'rat', (1,2),          (3,4)),
 ("B_S3RES   f=2.0, c=11/10        (exactly resonant)", 'irr', 1/mp.pi,        mp.mpf(11)/(20*mp.pi)),
 ("C_BADAPP  (t,t^2), t=2cos(2pi/7)  cubic field",      'irr', mp.frac(tt),    mp.frac(tt**2)),
 ("D_RAND1   uniform, master seed 20260816 stream 1",   'irr', RND[0][0],      RND[0][1]),
 ("D_RAND2   uniform, master seed 20260816 stream 2",   'irr', RND[1][0],      RND[1][1]),
 ("D_RAND3   uniform, master seed 20260816 stream 3",   'irr', RND[2][0],      RND[2][1]),
 ("E_W07GEN  f=2pi*phi, c=2pi*phi^2 (W-07 'GENERIC')",  'irr', mp.frac(phi),   mp.frac(phi**2)),
 ("F_VWA     very well approximable, Q=2999999",        'irr', vwa_a,          vwa_b),
]
READY = [("RS-G", "(p11,p10,p01)=(0.4,0.3,0.3)  all three classes live", (0.4,0.3,0.3)),
         ("RS-P", "S1 published p=(1/2,0,0,1/4,1/4) -> (0.5,0,0.5)",     (0.5,0.0,0.5))]

print("== B0  THE PAIRS AND THEIR RATIONAL RELATIONS — dim H IS THE OPERATIVE VARIABLE ==")
print("   H := closure{ (x^k,y^k) } <= T^2.  dim H = 2 - (number of independent relations")
print("   n1*alpha + n2*beta in Z).  Search |n1|,|n2| <= 60, tolerance 1e-25.")
DIM = {}
for lab,kind,A,B in CLASSES:
    a = mp.mpf(A[0])/A[1] if kind=='rat' else mp.frac(A)
    b = mp.mpf(B[0])/B[1] if kind=='rat' else mp.frac(B)
    isr_a = min((q for q in range(1,5000) if abs(a*q-mp.nint(a*q))<mp.mpf('1e-25')), default=None)
    isr_b = min((q for q in range(1,5000) if abs(b*q-mp.nint(b*q))<mp.mpf('1e-25')), default=None)
    rels=[]
    for n1 in range(0,61):
        for n2 in range(-60,61):
            if n1==0 and n2<=0: continue
            v = n1*a+n2*b
            if abs(v-mp.nint(v))<mp.mpf('1e-25') and np.gcd(n1,abs(n2))==1:
                rels.append((n1,n2,int(mp.nint(v))))
    dim = 0 if (isr_a and isr_b) else (1 if rels else 2)
    DIM[lab[:9]] = dim
    ordr = (int(np.lcm(isr_a,isr_b)) if dim==0 else None)
    print(f"   {lab:<52} dim H = {dim}"
          + (f"   |H| = {ordr}  (ATTAINED at every k = 0 mod {ordr})" if dim==0 else "")
          + (f"   relation: {rels[0][0]}a {rels[0][1]:+d}b = {rels[0][2]}" if dim==1 else "")
          + ("   no relation up to height 60  -> H = T^2" if dim==2 else ""))
print()

results = {}
for rtag, rlab, W in READY:
    w11,w10,w01 = W
    c1,c2,c3 = w11*w10, w11*w01, w10*w01
    print(f"================ READY STATE {rtag}: {rlab} ================")
    print(f"{'connection':<52}{'K':>9}{'F(K)=min(1-|Z_k|)':>21}{'argmin':>10}"
          f"{'SUM(1-|Z_k|)':>15}{'SUM/K':>10}{'-(1/K)SUMlog|Z_k|':>19}")
    for lab,kind,A,B in CLASSES:
        La = Lb = None
        if kind=='irr': La,Lb = limbs(mp.frac(A)), limbs(mp.frac(B))
        run_min=np.inf; run_arg=-1; run_sum=0.0; run_log=0.0
        counts=np.zeros(len(EPS),dtype=np.int64); hist=np.zeros(60,dtype=np.int64); records=[]
        first=True
        for bi in range(len(BOUNDS)-1):
            lo,hi = BOUNDS[bi],BOUNDS[bi+1]
            k = np.arange(lo+1,hi+1,dtype=np.float64)
            if kind=='rat':
                ki=np.arange(lo+1,hi+1,dtype=np.int64)
                u=(ki*A[0] % A[1])/A[1]; v=(ki*B[0] % B[1])/B[1]
            else:
                u=frac_irr(k,La); v=frac_irr(k,Lb)
            duv=u+v; duv-=np.round(duv); du=u-np.round(u); dv=v-np.round(v)
            S = 4.0*(c1*np.sin(np.pi*dv)**2 + c2*np.sin(np.pi*du)**2 + c3*np.sin(np.pi*duv)**2)
            S = np.minimum(S,1.0)
            absZ = np.sqrt(np.maximum(0.0,1.0-S))
            g = S/(1.0+absZ)                                   # = 1-|Z_k|
            run_sum += float(g.sum())
            run_log += float(np.log(absZ).sum())               # -inf if any |Z_k| = 0 exactly
            for i,e in enumerate(EPS): counts[i]+=int((g<e).sum())
            hist += np.histogram(np.log10(np.maximum(g,1e-30)),bins=60,range=(-30,0))[0]
            acc=np.minimum(np.minimum.accumulate(g),run_min)
            nr=np.flatnonzero(np.r_[acc[0]<run_min, acc[1:]<acc[:-1]])
            for j in nr: records.append((int(lo+1+j),float(g[j])))
            if acc[-1]<run_min:
                jj=int(np.argmin(g)); run_min=float(g[jj]); run_arg=lo+1+jj
            if hi in CHECKPOINTS:
                print(f"{(lab if first else ''):<52}{hi:>9}{run_min:>21.6e}{run_arg:>10}"
                      f"{run_sum:>15.4f}{run_sum/hi:>10.6f}"
                      f"{(-run_log/hi if np.isfinite(run_log) else np.inf):>19.9f}")
                first=False
        results[(rtag,lab[:9])]=dict(counts=counts.tolist(),hist=hist.tolist(),
                                     records=records[:2000],dim=DIM[lab[:9]],
                                     Fmin=run_min,argmin=run_arg,S=run_sum,L=run_log)
        print(f"{'':<52}{'N(eps):':>9} " + " ".join(f"{e:>8.0e}" for e in EPS))
        print(f"{'':<52}{'':>9} " + " ".join(f"{c:>8d}" for c in counts))
        print()
with open("m2_b_sweep.pkl","wb") as fh: pickle.dump(results,fh)
print("wrote m2_b_sweep.pkl  (N(eps), log10 histogram, record-setters, totals per class)")
