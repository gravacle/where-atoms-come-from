#!/usr/bin/env python3
"""REFUTER LANE (charge axis) -- part 3.  Number-level probes."""
import numpy as np
from itertools import product
from math import log, pi, gcd
from fractions import Fraction

a_v = np.array([1,1,1,0,0]); b_v = np.array([1,0,0,1,1]); V=5
def lat_rank(vs):
    vs=[v for v in vs if v!=(0,0)]
    if not vs: return 0
    for i in range(len(vs)):
        for j in range(i+1,len(vs)):
            if vs[i][0]*vs[j][1]-vs[i][1]*vs[j][0]!=0: return 2
    return 1
def drank(P):
    P=sorted(set(P))
    if len(P)<=1: return 0
    b=P[0]; return lat_rank([(p[0]-b[0],p[1]-b[1]) for p in P])

print("="*96); print("P1.  WHERE COULD '24 of 343' HAVE COME FROM?  every reading enumerated.")
print("="*96)
box=range(-3,4)
def cnt(pred): return sum(1 for t in product(box,repeat=3) if pred(t))
P=lambda t:[(t[0],t[0]),(t[1],0),(0,t[2])]
print("  |S|=3 distinct AND rank<=1                     :", cnt(lambda t: len(set(P(t)))==3 and drank(P(t))<=1))
print("  all three charges nonzero AND rank<=1          :", cnt(lambda t: all(x!=0 for x in t) and drank(P(t))<=1))
print("  rank<=1, any support size                      :", cnt(lambda t: drank(P(t))<=1))
print("  rank<=1 and |S|>=2                             :", cnt(lambda t: len(set(P(t)))>=2 and drank(P(t))<=1))
print("  |S|=3 distinct AND rank==1 (excl rank 0)       :", cnt(lambda t: len(set(P(t)))==3 and drank(P(t))==1))
print("  per-vertex {-3..3}^5, |S|>=3 and rank<=1       :",
      sum(1 for q in product(box,repeat=5)
          if len(set([(q[v]*a_v[v],q[v]*b_v[v]) for v in range(V)]))>=3
          and drank([(q[v]*a_v[v],q[v]*b_v[v]) for v in range(V)])<=1))
print("  the 6 solutions of s t = r(s+t) in {-3..3}^3   :",
      [t for t in product(box,repeat=3) if all(x!=0 for x in t) and t[1]*t[2]==t[0]*(t[1]+t[2])])
print("  NO reading of the claim's own box returns 24.")

print("\n"+"="*96); print("P2.  THE RECORD'S CHARGE NUMBER -1.200555 IS WRONG.  EXACT VALUE = log(3/10).")
print("="*96)
print("  q=(1,2,2,2,2) gives E = {(1,1),(2,0),(0,2)}.  Relative to (1,1) the exponents are")
print("  (0,0), (1,-1), (-1,1): Delta = <(1,-1)>, rank 1.  So |Z| = |0.4 + 0.3 w + 0.3 w^-1|")
print("  with w = chi_(1,-1), and lambda_B = m(0.3 + 0.4 w + 0.3 w^2) over the circle.")
r=np.roots([0.3,0.4,0.3]); print("  roots of 0.3 z^2 + 0.4 z + 0.3 :", r, " |roots| =", np.abs(r))
print("  both roots on the unit circle  =>  m = log(0.3) =", log(0.3))
for f,c,N in ((1.0,np.sqrt(2),4_000_000),(2.0,1.1,4_000_000),(0.7,2.9,4_000_000)):
    ks=np.arange(1,N+1)
    z=0.4*np.exp(1j*ks*(-f+c))+0.3*np.exp(1j*ks*(-2*f))+0.3*np.exp(1j*ks*(2*c))
    print(f"  direct schedule-B at (f,c)=({f},{c}), N={N}: {np.mean(np.log(np.abs(z))):.9f}")
print("  RECORD (W-03, established): -1.200555.   TRUTH: -1.203972804325936 = log(3/10).")
print("  discrepancy 3.4e-03 -- far outside the 1e-6 agreement the corpus works to.")

print("\n"+"="*96); print("P3.  HOMOGENEOUS CHARGE KILLS FORMATION AT EVERY q-TORSION CONNECTION (exact).")
print("="*96)
print("  G is computed EXACTLY as a subgroup of Q/Z: connection (f,c) = 2pi(alpha,beta),")
print("  alpha,beta rational; chi_(m,n) = exp(2 pi i(-m alpha + n beta)).")
def G_order(alpha,beta,pts,q):
    """order of G = <chi_x/chi_y> for exponents q*pts; exact, via Q/Z."""
    E=[(Fraction(q)*Fraction(m),Fraction(q)*Fraction(n)) for m,n in pts]
    vals=set()
    for i in range(len(E)):
        for j in range(len(E)):
            d=(E[i][0]-E[j][0], E[i][1]-E[j][1])
            vals.add((-d[0]*alpha + d[1]*beta) % 1)
    N=1
    for v in vals: N=N*v.denominator//gcd(N,v.denominator)
    ks=[int(v*N) for v in vals]; g=0
    for k in ks: g=gcd(g,k)
    g=gcd(g,N)
    return N//g if g else 1
pts=[(1,1),(1,0),(0,1)]
print("   (f,c)/2pi      q=1   q=2   q=3   q=4   q=5   q=6      |G| (1 means NO FORMATION EVER)")
for al,be in [(Fraction(1,2),Fraction(1,2)),(Fraction(1,3),Fraction(1,3)),
              (Fraction(1,4),Fraction(1,4)),(Fraction(1,5),Fraction(2,5)),
              (Fraction(1,6),Fraction(1,6))]:
    row=[G_order(al,be,pts,q) for q in (1,2,3,4,5,6)]
    print(f"   ({al},{be})".ljust(16), "  ".join(f"{x:4d}" for x in row))
print("  Every row: |S| = 3 and rank Delta = 2 throughout, yet |G| collapses to 1 for the")
print("  homogeneous charges q that are multiples of the connection's order.  S4-1 says")
print("  'rank 2' for all of these.  The formation criterion G != {1} FLIPS under a")
print("  homogeneous charge -- the exact case C2a declares safe.")
