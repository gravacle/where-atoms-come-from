"""S5-CHARGE  part A : the charge axis on K1 (S4-A recomputed)."""
import numpy as np
from itertools import product
from fractions import Fraction
from s5lib import *

np.set_printoptions(linewidth=200)
K1 = CW("K1", 5, [(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)], [{0:1,1:1,2:1}])
gF = {0:1,1:1,2:1}; gC = {3:1,4:1,5:1}
a,b = K1.classes(gF,gC)
pv = np.array([0.4,0.15,0.15,0.15,0.15])
pc = np.array([0.4,0.3,0.3]); Ecorner = np.array([[1,1],[1,0],[0,1]])

def hdr(s): print("\n"+"="*94+"\n"+s+"\n"+"="*94)

# ----------------------------------------------------------------------------------------
hdr("A1.  WHAT 'CHARGE' CAN LEGITIMATELY MEAN -- THREE CANDIDATES, TWO OF THEM VACUOUS")
# ----------------------------------------------------------------------------------------
print("""CANDIDATE 1  PER-EDGE charge q_e, charged holonomy exp(i sum_e q_e a_e).
   Gauge:  a_e -> a_e + theta_tgt - theta_src.  Around gamma_F = e1.e2.e3 (v0->v1->v2->v0):
     delta = q1(t1-t0) + q2(t2-t1) + q3(t0-t2) = t0(q3-q1) + t1(q1-q2) + t2(q2-q3)
   which vanishes for all theta IFF q1 = q2 = q3.   => PER-EDGE CHARGE IS NOT GAUGE INVARIANT
   unless it is constant around each loop, in which case it is a PER-LOOP charge.""")
rng = np.random.default_rng(5150001)   # SEED
worst_bad = 0.0; worst_ok = 0.0
for _ in range(2000):
    aa = rng.uniform(0,2*np.pi,6); th = rng.uniform(0,2*np.pi,5)
    src=[0,1,2,0,3,4]; tgt=[1,2,0,3,4,0]
    aa2 = aa + th[tgt]-th[src]
    qbad = np.array([1,2,3,1,1,1]); qok = np.array([2,2,2,3,3,3])
    hF = lambda A,q: np.exp(1j*(q[0]*A[0]+q[1]*A[1]+q[2]*A[2]))
    worst_bad = max(worst_bad, abs(hF(aa2,qbad)-hF(aa,qbad)))
    worst_ok  = max(worst_ok,  abs(hF(aa2,qok)-hF(aa,qok)))
print(f"  2000 random (connection,gauge) pairs, seed 5150001:")
print(f"     q=(1,2,3) on the face loop : max |W_F(gauged) - W_F| = {worst_bad:.6f}   NOT INVARIANT")
print(f"     q=(2,2,2) on the face loop : max |W_F(gauged) - W_F| = {worst_ok:.3e}   INVARIANT")
print("""
CANDIDATE 2  PER-LOOP charge (Q_F,Q_C):  W_F -> W_F^{Q_F}, W_C -> W_C^{Q_C}, i.e. (f,c) -> (Q_F f, Q_C c).
   This is a SURJECTIVE endomorphism of the connection torus T^2 whenever Q_F,Q_C != 0.
   It relabels points of T^2 and changes NO function's range.  VACUOUS as a modality.

CANDIDATE 3  PER-VERTEX FIBRE charge q_v:  U(1) acts on L_v by z -> z^{q_v}, so W(gamma)
   acts on L_v as W(gamma)^{q_v} and the exponent vector is  E_v = q_v (a_v, b_v).
   Gauge acts fibrewise as g_v^{q_v}; |s_v|^2 = p_v is untouched; Z_k is unchanged.
   THIS IS THE ONLY NON-VACUOUS READING, and it is the one S2:173-175 names.
   NOTE, and it is load-bearing: q_v may vary with v ONLY because W-01's operator
   (M_gamma s)(v) = W(gamma) s(v) is a DIAGONAL multiplication, not an edge transport.
   Edge transport L_u -> L_v would force q_u = q_v.  Charge inhomogeneity is definable
   here exactly to the extent that W-01 replaced transport by multiplication.""")

# gauge invariance of the charged functional, checked at the matrix level
rng = np.random.default_rng(5150002)
worst=0.0
for _ in range(200):
    q = rng.integers(-3,4,5)
    f,c = rng.uniform(0,2*np.pi,2)
    ph = rng.uniform(0,2*np.pi,5); s = np.sqrt(pv)*np.exp(1j*ph)
    g = np.exp(1j*rng.uniform(0,2*np.pi,5))
    E = exponents_from_charge(a,b,q)
    for k in (1,2,7):
        z1 = Z_of_k_matrix(k,a,b,q,s,f,c); z2 = Z_of_k_matrix(k,a,b,q,g*s,f,c)
        z3 = Z_of_k(k,E,pv,f,c)[0]
        worst = max(worst, abs(z1-z2), abs(z1-z3))
print(f"\n  CHARGED functional: 200 random (q in [-3,3]^5, connection, section, gauge), seed 5150002")
print(f"     max | Z_k(gauged) - Z_k |  and  | matrix - closed form |   =  {worst:.3e}")

# ----------------------------------------------------------------------------------------
hdr("A2.  HOMOGENEOUS CHARGE q: k -> qk.  IT DISCRIMINATES THE TWO SCHEDULES.")
# ----------------------------------------------------------------------------------------
print("  Z_k(q homogeneous) = Z_{qk}(unit charge)  -- identity, checked:")
w=0.0
rng = np.random.default_rng(5150003)
for _ in range(500):
    f,c = rng.uniform(0,2*np.pi,2); q=int(rng.integers(1,6)); k=int(rng.integers(1,50))
    E = exponents_from_charge(a,b,np.full(5,q))
    w=max(w, abs(Z_of_k(k,E,pv,f,c)[0]-Z_of_k(k*q,exponents_from_charge(a,b,np.ones(5,dtype=int)),pv,f,c)[0]))
print(f"     500 random samples, seed 5150003 : max deviation {w:.3e}")
print(f"\n  {'q':>3} {'lambda_A = log|Z_1|':>22} {'lambda_B (direct N=2e6)':>26} {'lambda_B generic exact':>24}")
f0,c0 = 1.0, np.sqrt(2)      # generic pair (1, sqrt2): rationally independent with 2pi
for q in (1,2,3,5,7):
    E = exponents_from_charge(a,b,np.full(5,q))
    lA = lambda_A(E,pv,f0,c0)
    lBd = lambda_B_direct(E,pv,f0,c0,N=2_000_000)
    lBe = mahler_generic(*support_exponents(E,pv), Nx=8192)
    print(f"  {q:>3} {lA:>22.9f} {lBd:>26.9f} {lBe:>24.9f}")
print("""  READ IT: lambda_A moves with q (it is log|Z_1| at the rescaled point (qf,qc));
  lambda_B does NOT (the closure of <(u^q,v^q)> has the same identity component as that of
  <(u,v)>, so on the generic set both are T^2).  HOMOGENEOUS CHARGE SEPARATES THE SCHEDULES.
  This bears on S4 FLAG F1 / CHOICE LEDGER C1, which S4 declared undecidable.""")

# ----------------------------------------------------------------------------------------
hdr("A3.  THEOREM C-1 : THE FORMATION CRITERION AT ARBITRARY EXPONENT DATA")
# ----------------------------------------------------------------------------------------
print("""  Setting: weights p_v>0 on a support S, exponent vectors E_v in Z^2, connection (u,v),
     Z_k = sum_v p_v u^{k m_v} v^{k n_v},
     L     = { (m,n) : u^m v^n = 1 }                      (relation lattice)
     Delta = < E_x - E_y : x,y in S >  <= Z^2             (difference lattice)
     G     = < chi_x/chi_y : x,y in S > <= U(1)           (W-02's formation group)

  (i)   G = phi(Delta) where phi(m,n) = u^m v^n;  hence  G  =  Delta / (Delta ^ L).
  (ii)  |Z_k| = 1  <=>  k*Delta  <=  L.
        [equality in |sum p_v z_v| <= 1 for unit z_v and a strictly positive convex weight
         holds iff all z_v coincide, i.e. iff phi(k(E_x-E_y))=1 for all x,y in S.]
  (iii) FORMATION NEVER OCCURS  <=>  Delta <= L  <=>  G = {1}.
  (iv)  sup_k |Z_k| = 1 is ATTAINED  <=>  G is FINITE (take k = exp(G)); if G is infinite,
        |Z_k| < 1 for every k but sup = 1 by equidistribution (W-02's recurrence).
  (v)   rank G = rank Delta - rank(Delta ^ L).
  (vi)  lambda_B = int_{L^perp} log|Z| dHaar, and since |Z| = |sum p_v chi_{E_v - E_{v0}}|
        depends only on the images of the DIFFERENCES in Z^2/L -- i.e. only on G --
             lambda_B  =  int_{Ghat} log| sum_v p_v psi(g_v) | dHaar_Ghat ,  g_v = [E_v - E_{v0}].
        [restriction L^perp -> Ghat is a surjection of compact groups, so it pushes Haar to Haar.]

  W-02's criterion  FORMATION <=> G != {1}  SURVIVES CHARGE VERBATIM.
  What fails under charge is S4-1, which COMPUTES rank G from the four corners of a square.""")

# --- (ii)/(iii) verified by enumeration over exponent configurations and connections ---
rng = np.random.default_rng(5150004)
tests=0; bad=0
for trial in range(4000):
    nS = int(rng.integers(2,5))
    E = rng.integers(-3,4,(nS,2))
    p = rng.uniform(0.1,1,nS); p/=p.sum()
    # exact rational connection so that L is exact
    den = int(rng.integers(1,13))
    pf = Fraction(int(rng.integers(0,den)),den); pc_ = Fraction(int(rng.integers(0,den)),den)
    f = 2*np.pi*float(pf); c = 2*np.pi*float(pc_)
    L = relation_lattice(f,c,exact_pair=(pf,pc_))
    D = difference_lattice(E)
    # G trivial  <=>  Delta <= L
    def inlat(g, B):
        if len(B)==0: return g[0]==0 and g[1]==0
        if len(B)==1:
            (m,n)=B[0]
            if m!=0:
                if g[0]%m: return False
                t=g[0]//m; return t*n==g[1]
            return g[0]==0 and n!=0 and g[1]%n==0
        B=np.array(B);
        sol=np.linalg.solve(B.T.astype(float), np.array(g,dtype=float))
        return np.allclose(sol, np.round(sol), atol=1e-9)
    Gtriv = all(inlat(g,L) for g in (D if D else []))
    zz = np.abs(Z_of_k(np.arange(1,61),E,p,f,c))
    never = bool(np.all(zz>1-1e-12))
    tests+=1
    if Gtriv != never: bad+=1
print(f"  (iii) verified: 4000 random (configuration, rational connection) pairs, seed 5150004")
print(f"        'Delta <= L' vs '|Z_k| = 1 for k<=60' :  {tests-bad} agree, {bad} disagree")

# --- (vi): lambda_B depends only on G-data.  Two different (E,L) with the same quotient data ---
print("\n  (vi) verified on a matched pair -- SAME quotient data, DIFFERENT (charge, lattice):")
# unit charge on the (1,1)-resonance  vs  charge (1,2,2,2,2) at a GENERIC connection
Echg = np.array([[1,1],[2,0],[0,2]])
v1 = lambda_B_exact(Ecorner,pc,((1,1),))
v2 = mahler_generic(Echg,pc,Nx=16384)
print(f"     unit charge, L = <(1,1)>            lambda_B = {v1:.12f}")
print(f"     charge q=(1,2,2,2,2), L = 0 (generic) lambda_B = {v2:.12f}")
print(f"     |difference| = {abs(v1-v2):.3e}    <-- CHARGE AND RESONANCE ARE INTERCHANGEABLE")
lBd = lambda_B_direct(Echg,pc,1.0,np.sqrt(2),N=4_000_000)
print(f"     direct schedule-B under charge, N=4e6, (f,c)=(1,sqrt2) : {lBd:.9f}")
print(f"     W-03's critic reported -1.200555 for q=(1,2,2,2,2); the converged value is")
print(f"     {v1:.9f} = m(0.3 z^2 + 0.4 z + 0.3) = log(0.3) (both roots on |z|=1).")
print(f"     log(0.3) = {np.log(0.3):.12f}")
