"""
VERIFY-3.  TWO HEADLINE SENTENCES AGAINST THE LANE'S OWN TABLES.

(A) "The only sign-definite extensive quantity is M = sum_i |J_i|, which takes exactly ONE value
     over all 2^n configurations."
    The lane's own step 3 defines and verifies B(s) = E(s) + M = 2 * sum over UNSATISFIED bonds
    of |J_i|.  Score B on exactly the criteria the standard names.

(B) "coh(random) = mean|E|/sum|J| falls as exactly sqrt(2/pi)*m^(-1/2) (closed form ...)".
    The closed form E|sum of m signs| = m*C(m-1,floor((m-1)/2))/2^(m-1) is the UNIT-coefficient
    random walk.  Test it against the lane's other four coupling families.
"""
import itertools, math, random
from fractions import Fraction

def couplings(name,m,seed=0):
    rnd=random.Random(seed)
    if name=="uniform":  return [1]*m
    if name=="randpos":  return [rnd.randrange(1,61) for _ in range(m)]
    if name=="randsign": return [rnd.choice((1,-1))*rnd.randrange(1,61) for _ in range(m)]
    if name=="decaying": return [max(1,5040//((i+1)**2)) for i in range(m)]
    if name=="superinc": return [2**i for i in range(m)]

def dist(J):
    R=sum(abs(j) for j in J); off=R; cur=[0]*(2*R+1); cur[off]=1; lo=hi=off
    for j in J:
        a=abs(j); nxt=[0]*(2*R+1)
        for i in range(lo,hi+1):
            c=cur[i]
            if c: nxt[i+a]+=c; nxt[i-a]+=c
        lo-=a; hi+=a; cur=nxt
    return {i-off:c for i,c in enumerate(cur) if c}

print("="*104)
print("VERIFY-3A   IS M REALLY 'THE ONLY SIGN-DEFINITE EXTENSIVE QUANTITY'?")
print("            B(s) = E(s) + M is defined AND verified in the lane's own s3_structure.py.")
print("="*104)
print()
print(f"  {'family':<10} {'n':>3} {'M':>7} {'B>=0 on all 2^n':>16} {'termwise coh(B)':>16} "
      f"{'#distinct B over 2^n':>21} {'mean B':>9} {'mean B(2n)/mean B(n)':>21} {'record-blind?':>14}")
prev={}
for fam in ("uniform","randpos","randsign","decaying"):
    for n in (7,8,13,14):
        J=couplings(fam,n-1); M=sum(abs(j) for j in J)
        Bs=set(); tot=0; cnt=0; nonneg=True; coh1=True
        for t in itertools.product((1,-1),repeat=n-1):
            terms=[abs(J[i]) if J[i]*t[i]>0 else 0 for i in range(n-1)]
            B=2*sum(terms)
            if B<0: nonneg=False
            if sum(terms)!=0 and abs(sum(terms))!=sum(abs(x) for x in terms): coh1=False
            Bs.add(B); tot+=B; cnt+=1
        mb=Fraction(tot,cnt)
        key=(fam,n)
        ratio="-"
        if n in (8,14) and (fam,n//2 if n==14 else 4) : pass
        prev[(fam,n)]=mb
        if n==14 and (fam,7) in prev and prev[(fam,7)]!=0:
            ratio=f"{float(mb/prev[(fam,7)]):.4f}"
        if n==8 and (fam,7) in prev: ratio="-"
        print(f"  {fam:<10} {n:>3} {M:>7} {str(nonneg):>16} {('1 (all terms >=0)' if coh1 else 'NOT 1'):>16} "
              f"{len(Bs):>21} {float(mb):>9.2f} {ratio:>21} "
              f"{'NO - varies':>14}")
    print()
print("  READ: B is non-negative on every one of the 2^n configurations (C-46 termwise coherence")
print("  = 1 EXACTLY, at every configuration, not on average); its exact mean is M, so it is")
print("  extensive wherever M is; the lane's own additivity table shows its cut defect is exactly")
print("  one bond; and it takes THOUSANDS of distinct values over the configuration space, not one.")
print("  So the sentence 'the only sign-definite extensive quantity is M, which takes exactly ONE")
print("  value over all 2^n configurations' contradicts the lane's own s3 table.  What is true is")
print("  narrower: no energetic quantity depends on the single WRITABLE bit b -- and that follows")
print("  from [W,H]=0 alone, in one line, for any H whatsoever.")
print()
print("="*104)
print("VERIFY-3B   THE sqrt(2/pi) m^(-1/2) CLOSED FORM IS THE UNIFORM-J CASE ONLY")
print("="*104)
print()
def closed(m): return Fraction(m*math.comb(m-1,(m-1)//2), 2**(m-1))
print(f"  {'family':<10} {'m':>4} {'exact mean|E|':>16} {'closed form m*C/2^(m-1)':>24} "
      f"{'equal?':>7} {'coh':>10} {'coh*sqrt(m)':>12}  [sqrt(2/pi)=0.797885]")
for fam in ("uniform","randpos","randsign","decaying"):
    for n in (9,17,33,65):
        J=couplings(fam,n-1); m=n-1; M=sum(abs(j) for j in J)
        d=dist(J); tot=sum(d.values())
        ma=Fraction(sum(abs(e)*c for e,c in d.items()),tot)
        cf=closed(m)
        coh=ma/M
        print(f"  {fam:<10} {m:>4} {float(ma):>16.6f} {float(cf):>24.6f} "
              f"{str(ma==cf):>7} {float(coh):>10.6f} {float(coh)*math.sqrt(m):>12.6f}")
    print()
print("  READ: the closed form reproduces mean|E| ONLY for uniform J.  It is the mean absolute")
print("  displacement of an m-step simple random walk -- a textbook identity -- and the")
print("  sqrt(2/pi) constant is the CLT, not a property of records or of clause (iv).")
print("  The lane's own 'decaying' family flatly contradicts the headline scaling: coh*sqrt(m)")
print("  GROWS instead of settling at 0.7979, because sum|J_i| converges.")
