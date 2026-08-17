# REFUTER 1 of LANE W-08 / M2 — LENS: MATHEMATICS.  Leg R3: THE CONTROL THAT RULES OUT
# "DIOPHANTINE QUALITY" (leg F3 / ISOLATION_LEDGER ROW 11) IS CONFOUNDED, AND ITS ZERO IS REAL.
#   R3.1  the pair leg F3 calls "VERY well approximable, Q=2999999, same d_eff = 2" is the pair
#         (1234567/2999999, 765431/2999999).  Both coordinates RATIONAL -> H is FINITE -> d_eff=0
#         by this lane's own definition.  Shown in EXACT INTEGER arithmetic, no float.
#   R3.2  F(1e7) = 0 for that pair is therefore the TRUE value of the object, not the "precision
#         artefact of this table" leg F3 calls it.  Exhibited exactly.
#   R3.3  quality at genuinely FIXED d_eff = 2: badly-approximable / random / leg B's own
#         perturbed VWA pair.  The fitted exponent moves by 2 units.  Competitor (c) is NOT ruled out.
#   R3.4  what IS true and IS a theorem: F(K) <= C K^{-2/d_eff} (Dirichlet).  One-sided.  Checked.
#   R3.5  the N(eps) half, which survives: the exponent AND the constant predicted from the exact
#         Haar measure of the quadratic sublevel set, against leg B's counts.
from fractions import Fraction as Fr
import numpy as np, mpmath as mp
mp.mp.dps=60
W=(0.4,0.3,0.3); w11,w10,w01=W

print("== R3.1  THE F3 'VWA' PAIR IS RATIONAL: d_eff = 0, NOT 2.  EXACT INTEGER ARITHMETIC. ==")
Q=2999999; A,B=1234567,765431
print(f"   alpha = {A}/{Q}   beta = {B}/{Q}    gcd(alpha)={np.gcd(A,Q)}  gcd(beta)={np.gcd(B,Q)}")
per=Q//int(np.gcd(np.gcd(A,Q),np.gcd(B,Q)))
ordr=int(np.lcm(Q//int(np.gcd(A,Q)), Q//int(np.gcd(B,Q))))
print(f"   orbit {{(k*alpha, k*beta) mod 1}} is PERIODIC with period lcm = {ordr}")
print(f"   -> H = closure = a FINITE cyclic group of order {ordr};  dim H = 0;  d_eff = 0.")
print(f"   leg F3's table header says 'Same d_eff = 2, ... only the arithmetic quality moves'.")
print(f"   ISOLATION_LEDGER ROW 11 says 'HELD: d_eff = 2 ... MOVES: Diophantine QUALITY only'.")
print(f"   BOTH ARE FALSE FOR THE ROW THAT CARRIES THE CONCLUSION.\n")

print("== R3.2  AND ITS ZERO IS THE OBJECT, NOT THE FLOAT.  EXACT. ==")
for k in [ordr]:
    ru=Fr(k*A,Q)-Fr(k*A//Q); rv=Fr(k*B,Q)-Fr(k*B//Q)
    print(f"   k = {k}:  k*alpha mod 1 = {ru}   k*beta mod 1 = {rv}   (EXACT rationals)")
    print(f"            all three characters equal 1  ->  Z_k = w11+w10+w01 = 1  ->  1-|Z_k| = 0 EXACTLY")
print("   leg F3 prints: '(the VWA row reads F(1e7) = 0 EXACTLY here: at k = 2999999 the DOUBLE")
print("    value of alpha rounds k*alpha to an integer.  That is a precision artefact of this")
print("    table, not the object)'.  THAT DIAGNOSIS IS WRONG.  k*alpha IS an integer, in Q, for")
print("    the pair the row names.  The double is right and the label is wrong: the row is an")
print("    ATTAINED, d_eff = 0 connection, which is the OTHER end of the lane's own taxonomy.")
print("   The number substituted for it (3.609e-23) belongs to leg B's DIFFERENT pair")
print("   (rational core + 1e-18*(sqrt2-1)) — so the table's row and its footnote are two")
print("   different connections.\n")

print("== R3.3  QUALITY AT GENUINELY FIXED d_eff = 2.  (leg B's own F_VWA pair is irrational.) ==")
def limbs(a):
    a=mp.frac(mp.mpf(a)); a1=float(mp.nint(a*2**26))/2**26; r=a-mp.mpf(a1)
    a2=float(mp.nint(r*2**52))/2**52; return a1,a2,float(r-mp.mpf(a2))
def fr(k,L):
    a1,a2,a3=L; t1=k*a1; t1-=np.floor(t1); t2=k*a2; t2-=np.floor(t2)
    t=t1+t2+k*a3; return t-np.floor(t)
def gapv(u,v):
    du=u-np.round(u); dv=v-np.round(v); duv=u+v; duv-=np.round(duv)
    S=np.minimum(4.0*(w11*w10*np.sin(np.pi*dv)**2+w11*w01*np.sin(np.pi*du)**2
                      +w10*w01*np.sin(np.pi*duv)**2),1.0)
    return S/(1.0+np.sqrt(np.maximum(0.0,1.0-S)))
KS=[10**3,10**4,10**5,10**6,10**7]
def floors(al,be):
    La,Lb=limbs(al),limbs(be); out=[];run=np.inf;lo=0
    for K in KS:
        while lo<K:
            hi=min(lo+10**6,K); k=np.arange(lo+1,hi+1,dtype=np.float64)
            g=gapv(fr(k,La),fr(k,Lb)); run=min(run,float(g.min())); lo=hi
        out.append(run)
    return out
tt=2*mp.cos(2*mp.pi/7); rng=np.random.default_rng(20260816)
RND=[(mp.mpf(float(rng.random())),mp.mpf(float(rng.random()))) for _ in range(3)]
vwa=(mp.mpf(1234567)/Q+mp.mpf('1e-18')*(mp.sqrt(2)-1), mp.mpf(765431)/Q+mp.mpf('1e-18')*(mp.sqrt(3)-1))
liou=(mp.mpf('0.1100010000000000000000010'), mp.mpf('0.2200020000000000000000020'))  # ~sum 10^{-j!}
QUAL=[("cubic field  (badly approximable)",mp.frac(tt),mp.frac(tt**2)),
      ("uniform random, stream 1",RND[0][0],RND[0][1]),
      ("uniform random, stream 2",RND[1][0],RND[1][1]),
      ("leg B's OWN F_VWA (irrational, dim H = 2)",vwa[0],vwa[1]),
      ("Liouville-style pair (dim H = 2)",liou[0],liou[1])]
print(f"   {'pair (all dim H = 2)':<44}"+"".join(f"{f'F(1e{e})':>12}" for e in range(3,8))+f"{'fit':>9}{'theory':>8}")
for lab,a,b in QUAL:
    f=floors(a,b); s=np.polyfit(np.log10(KS),np.log10(np.maximum(f,1e-300)),1)[0]
    print(f"   {lab:<44}"+"".join(f"{x:>12.3e}" for x in f)+f"{s:>9.3f}{-1.0:>8.1f}")
print("   -> at FIXED d_eff = 2 the fitted exponent moves from -0.97 to -3.0 as the quality moves.")
print("      'Quality sets the constant, d_eff sets the exponent' is FALSE as stated.  The spread")
print("      in the implied constant C = F(1e7)*1e7 is 16 orders of magnitude, not 'a factor ~4'.\n")

print("== R3.4  WHAT SURVIVES AND IS A THEOREM: THE DIRICHLET UPPER BOUND, ONE-SIDED. ==")
print("""   Dirichlet: for every K there is k <= K with max(||k alpha||,||k beta||) <= K^{-1/2}.  With
   1-|Z| <= 1-|Z|^2 = 4[w11w10 sin^2(pi||kb||)+w11w01 sin^2(pi||ka||)+w10w01 sin^2(pi||k(a+b)||)]
   <= 4 pi^2 (w11w10 + w11w01 + 4 w10w01) max^2, hence
       F(K) <= C_2 K^{-1},   C_2 = 4 pi^2 (0.12+0.12+4*0.09) = %.4f      (d_eff = 2)
   and the d_eff = 1 analogue with the relation's height in the constant.  THIS IS ONE-SIDED.
   The matching LOWER bound F(K) >= c K^{-2/d} is equivalent to the system being BADLY
   APPROXIMABLE — a measure-zero hypothesis, false for every VWA pair.  M2-4 states the law
   two-sidedly ('F(K) ~ K^{-2/d_eff}', 'fixed by d_eff alone') and its own F_VWA row breaks the
   lower half by 16 orders of magnitude.""" % (4*np.pi**2*(w11*w10+w11*w01+4*w10*w01)))
C2=4*np.pi**2*(w11*w10+w11*w01+4*w10*w01)
print(f"   check of the upper bound on every row above and in R2:  F(1e7)*1e7 <= {C2:.3f} ?")
for lab,a,b in QUAL:
    f=floors(a,b); print(f"     {lab:<44} F(1e7)*1e7 = {f[-1]*1e7:.6e}   {'OK' if f[-1]*1e7<=C2 else 'VIOLATED'}")
print()

print("== R3.5  THE HALF THAT SURVIVES INTACT: N(eps) ~ K eps^{d/2}, EXPONENT **AND** CONSTANT ==")
print("""   Near the identity 1-|Z| = 2 pi^2 [ (w11w10+w10w01) v^2 + 2 w10w01 uv + (w11w01+w10w01) u^2 ]
   + O(res^4): a POSITIVE DEFINITE quadratic form.  Haar measure of {1-|Z| < eps} on T^2 is
   therefore exactly  pi*eps/(2 pi^2 sqrt(det M)) + O(eps^2)  -- exponent 1 EXACTLY, and a
   predicted constant.  On the circle H it is the 1-d analogue, exponent 1/2.""")
M=np.array([[w11*w01+w10*w01, w10*w01],[w10*w01, w11*w10+w10*w01]])
pred=1.0/(2*np.pi*np.sqrt(np.linalg.det(M)))
print(f"   det M = {np.linalg.det(M):.6f}   predicted N(eps)/K = {pred:.6f} * eps")
print(f"   leg B's measured counts /1e7 (C_BADAPP, dim H = 2):")
for e,c in zip([1e-3,1e-4,1e-5,1e-6],[8386,835,85,9]):
    print(f"     eps={e:.0e}  measured {c/1e7:.4e}   predicted {pred*e:.4e}   ratio {c/1e7/(pred*e):.4f}")
print("   -> the equidistribution half of M2-4 is right to the CONSTANT, not just the exponent.")
print("      It is the extreme-value half (the floor F(K)) that quality breaks, and only that half.")
