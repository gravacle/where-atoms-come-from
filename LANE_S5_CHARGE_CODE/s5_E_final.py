"""S5-CHARGE  part E : the ten carriers under charge, the restricted multiplicity decision,
and the numerical cross-checks."""
import numpy as np
from itertools import product, permutations
from fractions import Fraction
from s5lib import *
np.set_printoptions(linewidth=210)
def hdr(s): print("\n"+"="*100+"\n"+s+"\n"+"="*100)

# ---------------- carriers (B3 corrected) ----------------
C={}
C['B1 ']=(CW("B1",5,[(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)],[{0:1,1:1,2:1}]),{0:1,1:1,2:1},{3:1,4:1,5:1})
C['B2 ']=(CW("B2",5,[(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)],[{0:1,1:1,2:1},{3:1,4:1,5:1}]),{0:1,1:1,2:1},{3:1,4:1,5:1})
C['B1p']=(CW("B1p",6,[(0,1),(1,2),(2,0),(3,4),(4,5),(5,3),(0,3)],[{0:1,1:1,2:1}]),{0:1,1:1,2:1},{3:1,4:1,5:1})
C['B1q']=(CW("B1q",7,[(0,1),(1,2),(2,0),(3,4),(4,5),(5,3),(0,6),(6,3)],[{0:1,1:1,2:1}]),{0:1,1:1,2:1},{3:1,4:1,5:1})
es=[(0,5),(5,1),(1,6),(6,2),(2,7),(7,0),(0,8),(8,3),(3,9),(9,4),(4,10),(10,0)]
C['B1s']=(CW("B1s",11,es,[{i:1 for i in range(6)}]),{i:1 for i in range(6)},{i:1 for i in range(6,12)})
ge=[];ix={}
for i in range(3):
    for j in range(3): ix[('h',i,j)]=len(ge); ge.append((3*i+j,3*((i+1)%3)+j))
for i in range(3):
    for j in range(3): ix[('v',i,j)]=len(ge); ge.append((3*i+j,3*i+(j+1)%3))
gf=[{ix[('h',i,j)]:1, ix[('v',(i+1)%3,j)]:1, ix[('h',i,(j+1)%3)]:-1, ix[('v',i,j)]:-1} for i in range(3) for j in range(3)]
T=CW("T3x3",9,ge,gf); row0={ix[('h',0,0)]:1,ix[('h',1,0)]:1,ix[('h',2,0)]:1}
C['B0a']=(T,gf[1],row0); C['B0b']=(T,gf[0],row0)
oe=[(0,1),(0,2),(0,3),(0,4),(0,1),(0,2),(0,3),(0,4),(1,2),(2,3),(3,4),(4,1)]
of=[{0:1,8:1,1:-1},{1:1,9:1,2:-1},{2:1,10:1,3:-1},{3:1,11:1,0:-1},
    {4:1,8:1,5:-1},{5:1,9:1,6:-1},{6:1,10:1,7:-1},{7:1,11:1,4:-1}]
C['B3 ']=(CW("B3",5,oe,of),{0:1,8:1,1:-1},{6:1,10:1,3:-1})
C['B4 ']=(CW("B4",6,[(0,2),(2,1),(1,3),(3,0),(0,4),(4,1),(1,5),(5,0)],
             [{0:1,1:1,2:1,3:1},{0:1,1:1,2:1,3:1},{4:1,5:1,6:1,7:1},{4:1,5:1,6:1,7:1}]),
          {0:1,1:1,2:1,3:1},{0:1,1:1,5:-1,4:-1})
C['B5 ']=(CW("B5",4,[(0,1),(1,2),(2,3),(3,0)],[{0:1,1:1,2:1,3:1},{0:1,1:1,2:1,3:1}]),
          {0:1,1:1,2:1,3:1},None)

hdr("E1.  S4-B REDONE UNDER CHARGE.  TEN CARRIERS x THREE CHARGE REGIMES, SENSE U (p_v = 1/V).\n"
    "     REGIME 1 unit charge  |  REGIME 2 homogeneous q = 2  |  REGIME 3 q_v = 1 + (v mod 2)")
print(f"  {'carrier':6}{'chi':>4}{'b1':>3}{'b2':>3} {'classes (unit charge)':>34} "
      f"{'#pts':>5}{'rk':>3} {'lam_B^gen R1':>14} {'#pts':>5}{'rk':>3} {'lam_B^gen R2':>14} {'#pts':>5}{'rk':>3} {'lam_B^gen R3':>14}")
for nm in ['B0a','B0b','B3 ','B1 ','B4 ','B5 ','B2 ','B1p','B1q','B1s']:
    K,gF,gC = C[nm]; B=K.betti()
    if gC is None:
        print(f"  {nm:6}{B['chi']:>4}{B['b1']:>3}{B['b2']:>3} {'-- no gamma_C (b1 = 0) --':>34}   formation datum does not exist at ANY charge")
        continue
    a,b = K.classes(gF,gC); p=np.ones(K.nV)/K.nV
    cls={}
    for v in range(K.nV): cls[(int(a[v]),int(b[v]))]=cls.get((int(a[v]),int(b[v])),0)+1
    out=[]
    for q in [np.ones(K.nV,dtype=int), np.full(K.nV,2), np.array([1+(v%2) for v in range(K.nV)])]:
        E=exponents_from_charge(a,b,q); Es,ps=support_exponents(E,p)
        out.append((len(set(map(tuple,Es))), len(difference_lattice(Es)), mahler_generic(Es,ps,Nx=8192)))
    print(f"  {nm:6}{B['chi']:>4}{B['b1']:>3}{B['b2']:>3} {str(cls):>34} "
          + " ".join(f"{o[0]:>5}{o[1]:>3} {o[2]:>14.9f}" for o in out))
print("""  READ IT.  Under unit charge S4 found SIX distinct SENSE-U values across the family.
  Homogeneous charge q=2 changes NONE of them (Theorem C-3 where |S|=3; and the 4-class
  carriers B0a/B0b/B4/B1q also do not move, because q*corners is the corners rescaled and the
  reduction to the Delta basis is the same polynomial).  The alternating charge, which is
  inhomogeneous WITHIN classes, moves seven of the nine.""")

hdr("E2.  DIRECT SCHEDULE-B SIMULATION UNDER CHARGE -- the exact values cross-checked.")
print(f"  (f,c) = (1.0, sqrt2);  N = 2e6;  no closed form used in the 'direct' column.")
print(f"  {'configuration':38} {'exact / quadrature':>20} {'direct N=2e6':>16} {'dev':>10}")
pc=np.array([0.4,0.3,0.3])
cases=[("K1 unit charge (0.4,0.3,0.3)",np.array([[1,1],[1,0],[0,1]]),pc),
       ("K1 q=(1,2,2,2,2)",np.array([[1,1],[2,0],[0,2]]),pc),
       ("K1 q=(2,1,1,1,1)",np.array([[2,2],[1,0],[0,1]]),pc),
       ("K1 q=(1,1,2,1,1) 4 pts",np.array([[1,1],[1,0],[2,0],[0,1]]),np.array([.4,.15,.15,.3])),
       ("K1 q=(1,1,2,1,2) 5 pts",np.array([[1,1],[1,0],[2,0],[0,1],[0,2]]),np.array([.4,.15,.15,.15,.15])),
       ("4 corners charge (0,1,1,2)",np.array([[0,0],[1,0],[0,1],[2,2]]),np.array([.25]*4))]
for lab,E,p in cases:
    ex = mahler_generic(E,p,Nx=32768)
    di = lambda_B_direct(E,p,1.0,np.sqrt(2),N=2_000_000)
    print(f"  {lab:38} {ex:>20.9f} {di:>16.9f} {abs(ex-di):>10.1e}")

hdr("E3.  INTEGER PERIODICITY IN THE WINDING (the design's named S4 external contact) UNDER CHARGE.")
rng=np.random.default_rng(5150050)
w=0.0
for _ in range(500):
    q=rng.integers(-3,4,5); f,c=rng.uniform(0,2*np.pi,2); k1,k2=rng.integers(-5,6,2)
    K,gF,gC=C['B1 ']; a,b=K.classes(gF,gC); E=exponents_from_charge(a,b,q)
    w=max(w,abs(abs(Z_of_k(1,E,np.array([.4,.15,.15,.15,.15]),f,c)[0])
                -abs(Z_of_k(1,E,np.array([.4,.15,.15,.15,.15]),f+2*np.pi*k1,c+2*np.pi*k2)[0])))
print(f"  500 random (charge, connection, winding) samples, seed 5150050 : max deviation {w:.3e}")
print("  Integer periodicity is EXACT at every integer charge, because E_v in Z^2 always.")
print("  It would FAIL at fractional charge -- which is why the charge must be an integer, and")
print("  that is the only quantization statement this construction supports.")

hdr("E4.  CHARGE vs MULTIPLICITY -- restricted to the loop re-designations the carrier ALLOWS.")
print("""  S2 section 2: the based loop group of K1's 1-skeleton is FREE OF RANK 2 on (gamma_F, gamma_C),
  so a change of designated pair is an element of GL2(Z).  But S4's CHOICE LEDGER C4 pins the
  pair: gamma_F must BOUND and gamma_C must NOT.  gamma_F^m gamma_C^n bounds iff n = 0, so the
  admissible changes are gamma_F -> gamma_F^{+-1}, gamma_C -> gamma_C^{+-1} gamma_F^k -- a
  TRIANGULAR subgroup of GL2(Z), not all of it.  The equivalence test is run under BOTH.""")
def GLfull(m=3):
    out=[]
    for A in product(range(-m,m+1),repeat=4):
        M=np.array(A).reshape(2,2)
        if abs(int(round(np.linalg.det(M.astype(float)))))==1: out.append(M)
    return out
def GLtri(m=6):
    out=[]
    for s1 in (1,-1):
        for s2 in (1,-1):
            for k in range(-m,m+1):
                out.append(np.array([[s1,k],[0,s2]]))
                out.append(np.array([[s1,0],[k,s2]]))
    return out
def equiv(E1,p1,E2,p2,GL):
    E1=np.asarray(E1);E2=np.asarray(E2);n=len(E1)
    if len(E2)!=n: return None
    for M in GL:
        for s in (1,-1):
            F=s*(M@E1.T).T
            for perm in permutations(range(n)):
                t=E2[perm[0]]-F[0]
                if all(tuple(F[i]+t)==tuple(E2[perm[i]]) for i in range(n)) and \
                   all(abs(p1[i]-p2[perm[i]])<1e-12 for i in range(n)):
                    return (M,s,tuple(t),perm)
    return None
Emult=np.array([[2,1],[1,0],[0,1]]); pm=np.array([0.4,0.3,0.3])
GF=GLfull(3); GT=GLtri(6)
full=[];tri=[]
for q in product(range(-6,7),repeat=3):
    E=np.array([[q[0],q[0]],[q[1],0],[0,q[2]]])
    if len(set(map(tuple,E)))<3: continue
    if equiv(E,pm,Emult,pm,GF): full.append(q)
    if equiv(E,pm,Emult,pm,GT): tri.append(q)
print(f"\n  target multiplicity configuration E = {[tuple(x) for x in Emult]} , p = (0.4,0.3,0.3)")
print(f"  charge assignments (q11,q10,q01) in [-6,6]^3 matching it")
print(f"     under the FULL GL2(Z) (all loop bases, C4 ABANDONED)  : {len(full)}   e.g. {full[:6]}")
print(f"     under the TRIANGULAR subgroup (C4 RESPECTED)          : {len(tri)}   {tri[:6]}")
print("""  DECISION.  If one is willing to re-designate which loop is the curvature -- i.e. to give up
  S4's C4 and K1's whole reason for existing (S1 section 5) -- then charge can reproduce a
  multiplicity configuration.  If C4 is kept, IT CANNOT: charge and multiplicity are
  INEQUIVALENT modalities.  C11 was closed on an unexamined premise, and the answer depends on
  a choice made three stages earlier.""")

hdr("E5.  WHAT HAPPENS TO S3's THEOREM S3-2 UNDER CHARGE.")
print("""  S3-2: 'Take k_n = n and a ready state with p_0, p_1+p_2, p_3+p_4 > 0.  If (W_F,W_C) != (1,1)
  then sum (1-|Z_n|) = infinity and lambda < 0.'
  UNDER CHARGE THIS IS FALSE AS STATED.  Its proof step is: '|Z| = 1 requires the three unit
  vectors uv, u, v to coincide on the support; with p fully supported that forces u = v = 1.'
  Under charge the three characters are u^{q0}v^{q0}, u^{q1}, v^{q2}, and they coincide on a
  POSITIVE-DIMENSIONAL subgroup of T^2 whenever the exponent vectors are collinear.
  CORRECTED STATEMENT (Theorem C-1(iii)):  the crossing escapes recurrence iff Delta(S) is NOT
  contained in L -- full support is neither necessary nor sufficient.
  EXHIBITED: q = (1,2,2,2,2), full support, W_F = e^{2i} != 1, W_C = e^{-2i} != 1 :""")
E=np.array([[1,1],[2,0],[0,2]]); p=np.array([0.4,0.3,0.3])
zz=np.abs(Z_of_k(np.arange(1,2001),E,p,2.0,-2.0))
print(f"     min |Z_k| over k <= 2000 = {zz.min():.15f}   sum(1-|Z_k|) over k<=2000 = {(1-zz).sum():.3e}")
print(f"     |Omega_2000| = prod |Z_k| = {np.exp(np.log(zz).sum()):.15f}   THE RECORD NEVER FORMS.")
print(f"     Every one of the nine requirements that depends on |Omega_N| -> 0 fails here.")
print("     The escape is CONDITIONAL IN A THIRD WAY, flagged nowhere in the corpus: THE CHARGE.")
