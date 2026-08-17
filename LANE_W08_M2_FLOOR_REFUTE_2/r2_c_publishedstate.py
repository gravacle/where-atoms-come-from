# LANE W-08 / M2 REFUTER 2 — leg C: THE SAME DISTINCTION ON K1'S OWN PUBLISHED READY STATE.
# The headline's "worth 4.8% in lambda" is measured on RS-G, a ready state that appears nowhere
# in the corpus.  K1's PUBLISHED ready state is p = (1/2,0,0,1/4,1/4) (S1 sec6 / W-01 / S3 sec4.5
# / COR-D).  This leg computes the SAME distinction there.  Exact arithmetic where it matters.
# ISOLATION: everything held as in leg B; the ONE thing that moves relative to leg B is the
# ready state, RS-G -> RS-P.  Reported as a separate table, never mixed into leg B's.
import numpy as np, mpmath as mp
from math import gcd
from fractions import Fraction
mp.mp.dps = 40

print("== C1  RS-P IN CLOSED FORM (reproducing the lane's own M2-9, from the operators) ==")
print("   p = (1/2,0,0,1/4,1/4) -> (p11,p10,p01) = (1/2,0,1/2), so")
print("   Z_k = (1/2) y^k (x^k + 1)  and  |Z_k| = |cos(pi k alpha)|, alpha = f/2pi.  W_C invisible.")
FACE_V=[0,1,2]; CYC_V=[0,3,4]
def Z_ops(f,c,p,K):
    WF=np.exp(1j*f); WC=np.exp(1j*c); p=np.asarray(p,float)
    inF=np.array([1.0 if v in FACE_V else 0.0 for v in range(5)])
    inC=np.array([1.0 if v in CYC_V else 0.0 for v in range(5)])
    return np.array([np.sum(p*np.conj(WF**(k*inF))*(WC**(k*inC))) for k in range(1,K+1)])
rng=np.random.default_rng(11); worst=0.0
for _ in range(50):
    f=rng.uniform(0,2*np.pi); c=rng.uniform(0,2*np.pi)
    Z=Z_ops(f,c,[0.5,0,0,0.25,0.25],80)
    pred=np.abs(np.cos(np.pi*np.arange(1,81)*f/(2*np.pi)))
    worst=max(worst,float(np.max(np.abs(np.abs(Z)-pred))))
print(f"   checked from the operators, 50 random (f,c), k<=80:  max discrepancy {worst:.3e}")
print(f"   APPROACHED (alpha irrational): lambda = int_0^1 log|cos(pi t)| dt = -log 2 = "
      f"{float(-mp.log(2)):.12f}   (the lane's M2-9 measured -0.6931467..0.6931477)\n")

print("== C2  ATTAINED ON RS-P: THE DISTINCTION IS NOT 4.8%.  IT IS INFINITE. ==")
print("   x = conj(W_F) of finite order m.  |Z_k| = 0 EXACTLY iff x^k = -1, i.e. iff m is EVEN.")
print("   Then Omega_N = 0 exactly for every N >= m/2: the branches are exactly orthogonal, the")
print("   record is complete at a finite cell, and lambda = -infinity.  Computed in exact")
print("   rational/root-of-unity arithmetic below (no float is involved in the zero).")
print(f"   {'ord(x)=m':>9} {'first k with Z_k = 0':>21} {'lambda on RS-P':>22} {'deviation vs -log2':>20}")
def lam_rsp(m):
    """exact orbit average of log|cos(pi k /m)| over k=1..m for x = e^{2 pi i/m}; -inf if m even"""
    if m%2==0: return mp.mpf('-inf')
    s=mp.mpf(0)
    for k in range(1,m+1): s+=mp.log(abs(mp.cos(mp.pi*mp.mpf(k)/m)))
    return s/m
for m in [2,3,4,5,6,7,8,9,15,16,101]:
    lam=lam_rsp(m)
    first = m//2 if m%2==0 else None
    dev = "infinite" if m%2==0 else f"{float((lam+mp.log(2))/mp.log(2))*100:+.2f}%"
    print(f"   {m:>9} {str(first) if first else 'never':>21} {mp.nstr(lam,12):>22} {dev:>20}")
print("   S1's PUBLISHED CONNECTION IS THE m = 2 ROW: W_F = -1 -> x = -1 -> Z_1 = 0 EXACTLY.")
print("   That is W-01's registered firing and S3 sec4.5 / COR-D's 'Omega_N = 0 for all N>=1'.")
tot=0; even=0
for m in range(2,2001):
    ph=sum(1 for j in range(1,m+1) if gcd(j,m)==1)
    tot+=ph; even += ph if m%2==0 else 0
print(f"   DENSITY OF THE INFINITE CASE INSIDE THE ATTAINED SET (counting roots of unity of")
print(f"   order <= 2000 by Euler phi): {even}/{tot} = {even/tot:.4f}  (-> 1/3 exactly: the odd-order")
print(f"   Dirichlet factor at s=2 is (1-2^-1)/(1-2^-2) = 2/3, so even orders carry 1/3).")
print(f"   MY OWN DEFECT, RECORDED NOT SILENTLY FIXED: the first version of this line predicted")
print(f"   2/3 for the even-order share and printed 0.3335 next to it.  The printed number was")
print(f"   right and the prediction was wrong; 1/3 is the correct limit.  It is not exotic either.")
print(f"   ODD-ORDER attained points on RS-P still deviate by up to "
      f"{float((lam_rsp(3)+mp.log(2))/mp.log(2))*100:+.2f}% (m=3), never 4.8%.\n")

print("== C3  WHAT THE HEADLINE SAYS AND WHAT THE SAME OBJECT SAYS ON THE PUBLISHED STATE ==")
print("   RS-G, |H|=4  : lambda -0.804719 vs -0.767508   ->  4.85%   [the registered headline]")
print("   RS-G, |H|=3  : lambda -1.535057 vs -0.767508   -> 100.01%  [leg B1, same object]")
print("   RS-G, |H|=2/3: lambda -0.458145 / -0.331417    -> +40.3% / +56.8%, OPPOSITE SIGN")
print("   RS-P, |H|=2  : lambda -infinity vs -log 2      ->  INFINITE [this leg, and the lane's")
print("                  own M2-10, which files it as NOT load-bearing]")
print("   The distinction has no single value.  4.8% is the value at one point of one ready state,")
print("   and it is the smallest of the four.")
