# LANE W-08 / M2 REFUTER 2 — leg D: THE SCHEDULE AXIS, AND WHETHER "ABSOLUTE" MEANS ANYTHING.
# The lane: "nearly cosmetic on the honest schedule (4.8%) ... ABSOLUTE on the schedule axis".
# The two halves are measured with DIFFERENT YARDSTICKS: a relative % on one axis and an
# exact-zero test on the other.  This leg applies each yardstick to BOTH axes.
# ISOLATION LEDGER (this leg): HELD - carrier K1, RS-G, K, the schedule RULE (write the delta*K
# cells of smallest 1-|Z_k|), the kernel, double precision with a two-limb phase split.
#   D1 MOVES: the order n of the attained point AND, in lockstep, the density delta = 1/n at
#             which the approached connection is charged.  Both columns move together BY DESIGN
#             and the row says so: this is a MATCHED-DENSITY comparison, not an isolation of n.
#   D2 MOVES: K only, with delta_K = K^{-1/2}; connection and ready state fixed.
import numpy as np, mpmath as mp
mp.mp.dps=40
W11,W10,W01=0.4,0.3,0.3
LAM_T2=-0.767507880357776

def frac_split(alpha,k):
    """two-limb frac(k*alpha): a_hi exact in double for k<2^27, residual carries ~1e-17."""
    a_hi=np.round(alpha*2**26)/2**26; a_lo=alpha-a_hi
    return (k*a_hi)%1.0 + k*a_lo
def gvals(alpha,beta,K,chunk=2*10**6):
    """1-|Z_k| for k=1..K, from |Z|^2 = 1 - sum_{j<l} w_j w_l |chi_j^k-chi_l^k|^2,
       chi_0=xy, chi_F=x, chi_C=y with x=e^{-2pi i alpha}, y=e^{+2pi i beta} (verified leg A1)."""
    out=np.empty(K); done=0
    while done<K:
        c=min(chunk,K-done); k=np.arange(done+1,done+c+1,dtype=np.float64)
        u=frac_split(alpha,k); v=frac_split(beta,k)
        du=u-np.round(u); dv=v-np.round(v); duv=u+v; duv=duv-np.round(duv)
        S=4*(W11*W10*np.sin(np.pi*dv)**2+W11*W01*np.sin(np.pi*du)**2+W10*W01*np.sin(np.pi*duv)**2)
        S=np.minimum(S,1.0)
        out[done:done+c]=S/(1.0+np.sqrt(np.maximum(0.0,1.0-S))); done+=c
    return out
def gvals_rat(a,b,n,K):
    k=np.arange(1,K+1,dtype=np.int64); u=(k*a%n)/n; v=(k*b%n)/n
    du=u-np.round(u); dv=v-np.round(v); duv=u+v; duv=duv-np.round(duv)
    S=np.minimum(4*(W11*W10*np.sin(np.pi*dv)**2+W11*W01*np.sin(np.pi*du)**2
                    +W10*W01*np.sin(np.pi*duv)**2),1.0)
    return S/(1.0+np.sqrt(np.maximum(0.0,1.0-S)))

astar=float(np.mod(2*np.cos(2*np.pi/7),1)); bstar=float(np.mod((2*np.cos(2*np.pi/7))**2,1))
K=10**7
print("== D1  THE MATCHED-DENSITY COMPARISON THE LANE DID NOT RUN ==")
print("   An ATTAINED connection of order n has EXACTLY ONE blank positive-density schedule and")
print("   its density is EXACTLY 1/n (|Z_k| = 1 iff n | k).  So 'positive density' is not a free")
print("   parameter: it is 1/n.  Charge the APPROACHED connection at the SAME density and the")
print("   'absolute' gap becomes a number that goes to zero with n -- exactly like the honest-")
print("   schedule gap the lane calls nearly cosmetic.   K = 1e7, RS-G, cubic pair as approached.")
g_app=gvals(astar,bstar,K); g_app.sort()
print(f"   {'n = |H|':>8} {'delta=1/n':>11} {'ATTAINED total':>16} {'APPROACHED total':>18} "
      f"{'per written cell':>18} {'|Omega| approached':>19} {'honest-axis gap':>17}")
for n in [2,4,5,16,256,1024,4096,65536,1000003]:
    a=int(round(n*astar))%n; b=int(round(n*bstar))%n
    from math import gcd
    order=n//gcd(gcd(a,b),n)
    lam_n=float(np.mean(np.log(np.maximum(1e-300,1.0-gvals_rat(a,b,n,n)))))
    J=max(1,int(K/n)); tot=float(g_app[:J].sum())
    print(f"   {order:>8} {1.0/n:>11.3e} {0.0:>16.6e} {tot:>18.6e} {tot/J:>18.6e} "
          f"{np.exp(-tot):>19.6e} {abs(lam_n-LAM_T2)/abs(LAM_T2)*100:>16.2f}%")
print("   READ THE LAST TWO COLUMNS TOGETHER.  Both gaps -- the 'absolute' one and the 'cosmetic'")
print("   one -- are functions of the SAME variable n and both go to zero with it.  At n = 65536")
print("   the adversary on an APPROACHED connection, writing at the same density the attained")
print("   blank schedule uses, leaves |Omega| = 0.999+ : operationally as blank as exactly blank.")
print("   The dichotomy 'cosmetic here, absolute there' is an artefact of using a RELATIVE")
print("   yardstick on one axis and an EXACT-ZERO yardstick on the other.\n")

print("== D2  THE STIPULATION DOING ALL THE WORK: 'POSITIVE DENSITY' IS NOT THE CORPUS'S TEST ==")
print("   W-02's registered mechanism sentence: 'All the content sits in the divergence")
print("   Sum (1 - z_n) = inf.'  DIVERGENCE is the corpus's durability test, not density.")
print("   On an APPROACHED connection take delta_K = K^{-1/2}: the number of cells written grows")
print("   without bound (J = sqrt(K)) and the record stays BOUNDED.  Measured:")
print(f"   {'K':>10} {'J = K^1/2':>10} {'SUM over the J best cells':>27} {'|Omega|':>12}")
for e in [4,5,6,7]:
    Kx=10**e; gx=gvals(astar,bstar,Kx); gx.sort(); J=int(Kx**0.5)
    t=float(gx[:J].sum()); print(f"   {Kx:>10} {J:>10} {t:>27.6f} {np.exp(-t):>12.6f}")
print("   -> unboundedly many writes, a record that never accumulates: durability fails on an")
print("      APPROACHED connection too, by the corpus's own criterion.  The lane's separation")
print("      needs the extra stipulation 'at FIXED positive density', which appears nowhere in")
print("      the register and which the lane itself concedes it cannot prove (self-flag 3).\n")

print("== D3  THE OTHER YARDSTICK, APPLIED TO THE AXIS THE LANE CALLED COSMETIC ==")
print("   On the ATTAINED connection the honest schedule ALSO contains exactly-blank cells at")
print("   positive density -- a fraction 1/|H| of the record's own cells contribute EXACTLY 0.")
for n,lab in [(4,"A_S1PUB |H|=4"),(2,"|H|=2"),(3,"|H|=3")]:
    a,b=(2,3) if n==4 else ((1,1) if n==2 else (1,2))
    g=gvals_rat(a,b,n,10**6)
    print(f"   {lab:<16} cells with 1-|Z_k| = 0 exactly, k <= 1e6: {int((g<=0).sum()):>8} "
          f"= {(g<=0).mean():.4f} of the record")
print("   So 'a positive-density set of exactly blank cells' is a property of the ATTAINED")
print("   CONNECTION, visible on the honest schedule, not a discovery about schedules.  And it")
print("   follows in one line from attainment: if |Z_{k0}| = 1 then all three characters agree")
print("   at k0, hence at every multiple of k0, hence |Z_{m k0}| = 1 for all m.  No computation,")
print("   no Diophantine input, no sorting of 1e6 cells is needed to know it.\n")

print("== D4  AND THE SAME ONE-LINE ARGUMENT BOUNDS THE ADVERSARY'S REACH ON BOTH SIDES ==")
print("   The lane's E3 delta=0.25 row is the attained blank set exactly (density 1/4 = 1/|H|).")
print("   For delta > 1/|H| the attained adversary is forced off the blank set and its total")
print("   jumps to a POSITIVE, linear-in-K number -- the lane's own delta=0.5 row, 1.5e5 at")
print("   K = 1e6, which is LARGER THAN ALL FOUR approached rows at the same delta")
print("   (lane E3, delta=0.5: attained 1.500e5 vs 1.330e5, 1.330e5, 1.330e5, 1.090e5).")
g4=gvals_rat(2,3,4,10**6); g4.sort()
for d in [0.25,0.26,0.3,0.5]:
    print(f"      attained |H|=4, delta={d:<5}: total over the best delta*K cells = "
          f"{g4[:int(d*10**6)].sum():.6e}")
print("   -> the attained connection is blank ONLY at density <= 1/|H| and is the WORST carrier")
print("      for the adversary above it.  'ATTAINED => the adversary wins' is true on a single")
print("      density and false immediately above it.")
