"""S5-CHARGE  part B : the systematic charge sweep on K1, and everything S4-A claimed."""
import numpy as np
from itertools import product, permutations
from fractions import Fraction
from s5lib import *
np.set_printoptions(linewidth=200)

K1 = CW("K1", 5, [(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)], [{0:1,1:1,2:1}])
gF={0:1,1:1,2:1}; gC={3:1,4:1,5:1}
a,b = K1.classes(gF,gC)
pv = np.array([0.4,0.15,0.15,0.15,0.15])
pc = np.array([0.4,0.3,0.3]); Ecorner = np.array([[1,1],[1,0],[0,1]])
def hdr(s): print("\n"+"="*96+"\n"+s+"\n"+"="*96)

def combo_of(D):
    """rank-1 Delta -> the single real combination of (f,c) that |Z| sees."""
    d,e = primitive(D[0])
    d1,d2 = d
    return f"{-d1:+d}*f {d2:+d}*c   (period 2pi/{e})" if e>1 else f"{-d1:+d}*f {d2:+d}*c"

# --------------------------------------------------------------------------------------
hdr("B1.  CLASS-HOMOGENEOUS CHARGE SWEEP ON K1.  q = (q11, q10, q01) in {-3..3}^3, 343 runs.\n"
    "     weights (p11,p10,p01) = (0.4,0.3,0.3)  [S3's own state].  FULL SUPPORT THROUGHOUT.")
# --------------------------------------------------------------------------------------
rows=[]; rk_count={0:0,1:0,2:0}
for q11,q10,q01 in product(range(-3,4),repeat=3):
    E = np.array([[q11,q11],[q10,0],[0,q01]])
    Es,ps = E, pc
    D = difference_lattice(Es)
    rk = len(D); rk_count[rk]+=1
    rows.append((q11,q10,q01,rk,D))
print(f"  rank Delta over the 343 assignments:  rank0 = {rk_count[0]}, rank1 = {rk_count[1]}, rank2 = {rk_count[2]}")
print(f"  UNIT CHARGE (1,1,1) is one of the {rk_count[2]} rank-2 points.")
print("\n  THE RANK-1 AND RANK-0 ASSIGNMENTS AT FULL SUPPORT (S4-1 says these cannot exist):")
n1=0
for (q11,q10,q01,rk,D) in rows:
    if rk<=1 and (q11,q10,q01)!=(0,0,0):
        n1+=1
        if n1<=18:
            if rk==0: print(f"     q={(q11,q10,q01)!s:>12}  rank Delta = 0  -> NO FORMATION AT ANY CONNECTION")
            else:     print(f"     q={(q11,q10,q01)!s:>12}  rank Delta = 1  -> lambda sees only  {combo_of(D)}")
print(f"     ... {n1} such assignments in total (of 343), all with |S| = 3.")

print("\n  EXACT GENERIC RATES FOR A REPRESENTATIVE SET (Nx = 16384 where 2-variable):")
reps = [(1,1,1),(1,2,2),(2,1,1),(1,3,3),(2,3,3),(1,-1,-1),(3,2,2),(2,2,2),(1,1,2),(1,2,3),(0,1,1),(1,0,1),(2,1,3)]
print(f"     {'q=(q11,q10,q01)':>18} {'rkD':>4} {'index':>6}  {'lambda_B generic':>18}  {'lambda_A(1,sqrt2)':>18}  sees")
for q in reps:
    E = np.array([[q[0],q[0]],[q[1],0],[0,q[2]]])
    D = difference_lattice(E); rk=len(D)
    lg = mahler_generic(E,pc,Nx=16384)
    lA = lambda_A(E,pc,1.0,np.sqrt(2))
    idx = lattice_index(list(np.array([E[i]-E[0] for i in (1,2)])))
    sees = "both f and c" if rk==2 else (combo_of(D) if rk==1 else "nothing (never forms)")
    print(f"     {str(q):>18} {rk:>4} {str(idx):>6}  {lg:>18.9f}  {lA:>18.9f}  {sees}")

# --------------------------------------------------------------------------------------
hdr("B2.  ARMED FALSIFIER F-A, RE-RUN UNDER CHARGE.\n"
    "     S4: 'F-A does not fire on full support'.  UNDER CHARGE IT FIRES ON FULL SUPPORT.")
# --------------------------------------------------------------------------------------
qFA = np.array([1,2,2,2,2])
E = exponents_from_charge(a,b,qFA)
print(f"  q = {qFA} on K1  ->  E_v = {[tuple(x) for x in E]}")
D = difference_lattice(np.array([[1,1],[2,0],[0,2]]))
print(f"  Delta = {D}  (rank {len(D)})    => |Z_k| is a function of  {combo_of(D)}  ALONE.")
print("  Algebraically:  Z_1 = 0.4 u v + 0.3 u^2 + 0.3 v^2 = u v ( 0.4 + 0.3 (u/v) + 0.3 (v/u) )")
print("                       = u v ( 0.4 + 0.6 cos(f+c) )        [u = e^{-if}, v = e^{ic}]")
rng=np.random.default_rng(5150010)
w=0.0
for _ in range(4000):
    f,c = rng.uniform(0,2*np.pi,2)
    lhs = abs(Z_of_k(1,np.array([[1,1],[2,0],[0,2]]),pc,f,c)[0])
    rhs = abs(0.4+0.6*np.cos(f+c))
    w=max(w,abs(lhs-rhs))
print(f"  |Z_1| vs |0.4 + 0.6 cos(f+c)| over 4000 random connections, seed 5150010 : max dev {w:.3e}")
print("\n  DIRECT WITNESSES -- pairs with the SAME product f+c but wildly different (f,c):")
for (f,c) in [(0.4,2.6),(1.5,1.5),(2.9,0.1),(-1.0,4.0)]:
    print(f"     f={f:+.3f} c={c:+.3f}  f+c={f+c:.3f}   lambda_A = {lambda_A(np.array([[1,1],[2,0],[0,2]]),pc,f,c):.12f}")
print("  S4's four witnesses at unit charge separated these to 5e-2.  Under charge they agree to 1e-12.")
print("\n  FOURIER SUPPORT of |Z_1|^2 on a 64x64 grid in (f,c)  [grid: f_j = 2pi j/64, c_k = 2pi k/64]:")
for (label,Ex) in [("unit charge q=(1,1,1,1,1)",np.array([[1,1],[1,0],[0,1]])),
                   ("charge   q=(1,2,2,2,2)",np.array([[1,1],[2,0],[0,2]])),
                   ("charge   q=(1,1,1,2,2)",np.array([[1,1],[1,0],[0,2]]))]:
    G=64
    fs = 2*np.pi*np.arange(G)/G
    FF,CC = np.meshgrid(fs,fs,indexing='ij')
    Zg = np.zeros_like(FF,dtype=complex)
    for i in range(len(Ex)):
        Zg += pc[i]*np.exp(1j*(-Ex[i,0]*FF + Ex[i,1]*CC))
    A = np.fft.fft2(np.abs(Zg)**2)/G**2
    sup = [( (m if m<=G//2 else m-G), (n if n<=G//2 else n-G) ) for m in range(G) for n in range(G) if abs(A[m,n])>1e-9]
    print(f"     {label:28s} support (m,n) = {sorted(sup)}   lattice rank = {lattice_rank(sup)}")
print("  PROPOSITION (charge form).  The Fourier support of |Z_1|^2 on the dual of T^2 is exactly")
print("  the difference SET {E_x - E_y}.  Hence lambda_A is a function of a single combination")
print("  iff rank Delta <= 1.  At unit charge with |S|>=3 the corners force rank 2 (S4's proof).")
print("  Under charge the corners are gone and rank 1 is reachable with full support.  QED")

# --------------------------------------------------------------------------------------
hdr("B3.  THE NON-FORMATION LOCUS AND THE EXACT-FIRING LOCUS GAIN A DIMENSION.")
# --------------------------------------------------------------------------------------
print("  Unit charge, full support:  no formation only at (f,c) = (0,0)  -- ONE POINT of T^2.")
print("  Charge q=(1,2,2,2,2):  Delta = <(1,-1)>, so Delta <= L  <=>  u v^{-1} = 1  <=>  f + c = 0 (2pi)")
print("  which is a CIRCLE.  Checked: |Z_k| = 1 for all k on that circle, for non-trivial holonomies:")
Ec2 = np.array([[1,1],[2,0],[0,2]])
for f in (0.3, 1.0, 2.0, np.pi, 5.0):
    c = -f
    zz = np.abs(Z_of_k(np.arange(1,301),Ec2,pc,f,c))
    print(f"     f={f:+.5f} c={c:+.5f}  (W_F,W_C) = ({np.exp(1j*f):.4f},{np.exp(1j*c):.4f})   min|Z_k| k<=300 = {zz.min():.15f}")
print("  DIMENSION OF THE NON-FORMATION LOCUS:  0 at unit charge, 1 under charge.")
print("\n  EXACT FIRING (lambda = -infinity).  Unit charge: exactly 2 points on T^2 (S4 section 3.1).")
print("  Under q=(1,2,2,2,2):  Z_1 = 0 <=> 0.4 + 0.6 cos(f+c) = 0 <=> f+c = +-arccos(-2/3):")
r = np.arccos(-2/3)
for f in (0.0, 1.0, 2.5):
    for sgn in (+1,-1):
        c = sgn*r - f
        print(f"     f={f:+.4f}  c={c:+.4f}   |Z_1| = {abs(Z_of_k(1,Ec2,pc,f,c)[0]):.3e}")
print(f"  arccos(-2/3) = {r:.9f}.  TWO CIRCLES, not two points: the firing locus is codimension 1.")

# --------------------------------------------------------------------------------------
hdr("B4.  S4-1 REPLACED.  THEOREM C-2 (the affine rank theorem) AND ITS ENUMERATION.")
# --------------------------------------------------------------------------------------
print("""  THEOREM C-2.  With E_v = q_v (a_v,b_v) and support S,
       rank G  =  rank Delta(S) - rank(Delta(S) ^ L),
       rank Delta(S) = dim of the AFFINE span over Q of { E_v : v in S } in Q^2.
  Corollaries:
    (C2a) If the charge is HOMOGENEOUS (q_v = q != 0 for all v in S), the configuration is
          q * (corners of the unit square) and S4-1 holds VERBATIM: rank 2 iff |S| >= 3.
          S4-1 is therefore a theorem about CHARGE HOMOGENEITY, not about charge 1.
    (C2b) With INHOMOGENEOUS charge, |S| >= 3 no longer implies rank 2: any number of
          exponent vectors on one affine line gives rank 1, and on one point gives rank 0.
    (C2c) The four-case |S|=2 table of S4-1 is replaced by: for |S| = 2, G = <phi(E_1-E_2)>,
          a single character, whose (m,n) is now an ARBITRARY element of Z^2 rather than one
          of the four differences of square corners.""")
print("\n  C2a checked: homogeneous charge q on all vertices, all 15 support subsets, q = 1..4:")
corners = [(0,0),(1,0),(0,1),(1,1)]
ok=True
for q in range(1,5):
    cnt={0:0,1:0,2:0}
    for mask in range(1,16):
        S=[corners[i] for i in range(4) if mask>>i & 1]
        E=np.array([[q*x,q*y] for (x,y) in S])
        cnt[len(difference_lattice(E))]+=1
    print(f"     q={q}:  rank2 = {cnt[2]}, rank1 = {cnt[1]}, rank0 = {cnt[0]}   (S4-1 predicts 5 / 6 / 4)")
    ok &= (cnt=={2:5,1:6,0:4})
print(f"     S4-1's 5/6/4 holds for every homogeneous charge tested: {ok}")
print("\n  C2b checked: the critic's counterexample and a systematic search on K1-like supports:")
print(f"     E = (1,0),(2,0),(3,0)  [|S|=3]  ->  Delta = {difference_lattice(np.array([[1,0],[2,0],[3,0]]))}  rank "
      f"{lattice_rank([(1,0),(2,0)])}   AGAINST S4-1")
cnt_bad=0; ex=[]
for q in product(range(0,4),repeat=3):
    E=np.array([[q[0],q[0]],[q[1],0],[0,q[2]]])
    if len(difference_lattice(E))<2 and len(set(map(tuple,E)))==3:
        cnt_bad+=1; ex.append(q)
print(f"     charges in {{0..3}}^3 on K1's three classes with |S|=3 (3 DISTINCT exponent vectors)")
print(f"     but rank Delta < 2 :  {cnt_bad} assignments, e.g. {ex[:8]}")

# --------------------------------------------------------------------------------------
hdr("B5.  MULTISET INVARIANCE.  PROVED AT UNIT CHARGE; BROKEN BY CHARGE.")
# --------------------------------------------------------------------------------------
print("""  PROOF at unit charge (rank L = 0, four corners).  By Jensen in y,
      lambda_B = m(p00 + p10 x + p01 y + p11 x y) = int_0^1 log max(|p00+p10 x|, |p01+p11 x|) dt.
  Now |A + B e^{i t}|^2 = A^2 + B^2 + 2 A B cos t is SYMMETRIC IN (A,B) POINTWISE IN t.
  Hence the integrand is invariant under (p00 p10) and under (p01 p11) separately; it is
  invariant under (p00 p01)(p10 p11) (swap the two arguments of max); and the substitution
  x <-> y is measure preserving and realises (p10 p01).  The transpositions (p00 p10),
  (p10 p01), (p01 p11) generate the full symmetric group S_4.  QED  -- this proves W-03's
  24/24 permutation invariance rather than sampling it.""")
rng=np.random.default_rng(5150020)
Ecor4 = np.array([[0,0],[1,0],[0,1],[1,1]])
w=0.0
base = np.array([0.4,0.3,0.2,0.1])
v0 = mahler_generic(Ecor4,base,Nx=16384)
vals=[]
for pi_ in permutations(range(4)):
    vals.append(mahler_generic(Ecor4,base[list(pi_)],Nx=16384))
print(f"  CHECKED: p = (0.4,0.3,0.2,0.1) on the four corners, all 24 permutations, Nx=16384")
print(f"     spread = {max(vals)-min(vals):.3e}   value = {v0:.12f}")
print("\n  UNDER CHARGE THE INVARIANCE DIES.  Exponents (1,1),(2,0),(0,2) [q=(1,2,2,2,2)],")
print("  weights permuted among the three positions -- these are the SAME MULTISET of weights:")
E3 = np.array([[1,1],[2,0],[0,2]])
seen=set()
for perm in permutations([0.4,0.3,0.3]):
    seen.add(round(mahler_generic(E3,np.array(perm),Nx=16384),12))
print(f"     weights (0.4,0.3,0.3) permuted : distinct lambda values = {sorted(seen)}")
for perm in [(0.4,0.3,0.3),(0.3,0.4,0.3),(0.3,0.3,0.4)]:
    print(f"       p = {perm}  ->  lambda = {mahler_generic(E3,np.array(perm),Nx=16384):.12f}")
print("  At unit charge the same three weights on classes (1,1),(1,0),(0,1) give one value:")
for perm in [(0.4,0.3,0.3),(0.3,0.4,0.3),(0.3,0.3,0.4)]:
    print(f"       p = {perm}  ->  lambda = {mahler_generic(Ecorner,np.array(perm),Nx=16384):.12f}")

# --------------------------------------------------------------------------------------
hdr("B6.  THE PINCH/SPECTATOR SYMMETRY UNDER CHARGE.")
# --------------------------------------------------------------------------------------
print("""  W-03: multiplying Z_k by conj(u)^k conj(v)^k fixes |Z_k| and maps class (a,b)->(1-a,1-b),
  so 00<->11 and 10<->01 exactly, at every connection.
  GENERAL FORM:  |Z| is invariant under  E_v -> t - E_v  for ANY t in Z^2
     [ E -> E + t  is multiplication by a unimodular monomial;  E -> -E is conjugation ].
  W-03's involution is the case t = (1,1).  SO THE SYMMETRY SURVIVES CHARGE -- as a symmetry
  of the FUNCTIONAL.  What does NOT survive is its reading as a carrier statement:""")
rng=np.random.default_rng(5150030); w=0.0
for _ in range(2000):
    nS=int(rng.integers(2,5)); E=rng.integers(-3,4,(nS,2)); p=rng.uniform(.1,1,nS); p/=p.sum()
    t=rng.integers(-3,4,2); f,c=rng.uniform(0,2*np.pi,2); k=int(rng.integers(1,30))
    w=max(w,abs(abs(Z_of_k(k,E,p,f,c)[0])-abs(Z_of_k(k,t-E,p,f,c)[0])))
print(f"  E -> t - E  leaves |Z_k| fixed : 2000 random samples, seed 5150030, max dev {w:.3e}")
rays = set()
for q in range(-6,7):
    for (x,y) in [(0,0),(1,0),(0,1),(1,1)]:
        rays.add((q*x,q*y))
print("\n  REALIZABILITY.  A charge configuration lives on the FOUR RAYS Z*(0,0), Z*(1,0),")
print("  Z*(0,1), Z*(1,1).  For the involution E -> t - E to map charge configurations to")
print("  charge configurations we need t - R contained in R.  Testing t = (T,T):")
for T in range(0,5):
    bad=[e for e in [(q*x,q*y) for q in range(0,4) for (x,y) in [(0,0),(1,0),(0,1),(1,1)]]
         if (T-e[0],T-e[1]) not in rays]
    print(f"     t=({T},{T}): exponents with charge <= 3 whose image leaves the rays: {len(bad)}  e.g. {bad[:4]}")
print("""  At unit charge every exponent has entries in {0,1} and t=(1,1) maps the four corners to
  themselves -- the pinch (1,1) to the spectator (0,0).  AT CHARGE >= 2 THE INVOLUTION LEAVES
  THE REALIZABLE SET: (1,1)-(2,0) = (-1,1) is on no ray.  The identity 'pinch = spectator'
  is therefore a UNIT-CHARGE COINCIDENCE, not a structural fact about the two vertex roles.""")
