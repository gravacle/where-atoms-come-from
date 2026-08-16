"""R1 -- the attack on:  'lambda is a formation rate of the carrier / of anything the
cell structure represents.  STOPS AT: the collapse of a spanning tree of the 1-skeleton.'"""
import numpy as np, sys, itertools
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_R_MAPS_REFUTER")
from rmlib import *

def hr(t): print("\n" + "=" * 78 + "\n" + t + "\n" + "=" * 78)
def row(tag, K):
    b, r = K.betti(); pi = K.pi(); nS, rG, S = rank_G(pi)
    lA = lambda_A(pi, 2.0, 1.1)
    lB = cassaigne_maillot(pi[3], pi[1], pi[2]) if pi[0] == 0 else float('nan')
    print(f"{tag:34s} V={len(K.verts)} E={len(K.edges)} F={len(K.faces)} chi={K.chi():+d} "
          f"b={b} d1d2={K.check_d2_zero()} pi={np.round(pi,6)} |S|={nS} rkG={rG} "
          f"FORMS={forms(pi)!s:5s} lA={lA:+.12f} lB={lB:+.12f}")
    return pi, lA, lB

# ============================================================ REPRODUCE THE CLAIM
hr("R1.0  THE CLAIM, REPRODUCED EXACTLY (SENSE U = uniform vertex weights)")
KU = K1()
piU, lAU, lBU = row("K1  (as handed, SENSE U)", KU)
# spanning tree of the 1-skeleton: e1,e2,e4,e5  (4 edges, spans all 5 vertices)
T = ["e1", "e2", "e4", "e5"]
KT = collapse_tree(KU, T, "K1 / spanning tree")
piT, lAT, lBT = row("K1/T  (tree collapsed)", KT)
print("\n  remaining edges after the collapse:", [(e[0], e[1] + '->' + e[2]) for e in KT.edges])
print("  loopF ->", KT.loopF, "   loopC ->", KT.loopC)
print("  |Z_k| on K1/T at 12 random connections, k=1..40:")
rng = np.random.default_rng(31415926)          # SEED PUBLISHED
worst = 0.0
for _ in range(12):
    WF, WC = np.exp(1j * rng.uniform(0, 2*np.pi)), np.exp(1j * rng.uniform(0, 2*np.pi))
    for k in range(1, 41):
        worst = max(worst, abs(abs(KT.Z_matrix(WF, WC, k)) - 1.0))
print("    max | |Z_k| - 1 | = %.3e   (so lambda = 0 under EVERY schedule)" % worst)

# ============================================================ ATTACK 1
hr("R1.1  ATTACK 1 -- THE CONTROL COULD NOT HAVE FAILED (vacuous by identity)")
print("""THEOREM R-1 (one-vertex vacuity).  Let K' be ANY CW complex with exactly one
vertex w, rank-one fibre, U(1) connection, and two designated loops.  Then
  (i) the incidence map v |-> (a_v,b_v) has a one-point image, so |S| = 1 and
      G = <chi_a/chi_b : a,b in S> = <chi_w/chi_w> = {1}: FORMS = False;
  (ii) any ready state has p_w = 1, so Z_k = chi_w^k * 1 with |chi_w| = 1,
       hence |Z_k| = 1 for every k, every connection, every charge, every schedule;
  (iii) lambda = lim (1/N) sum log|Z_{k_n}| = 0 EXACTLY, for every schedule.
No property of any MAP, of homotopy type, of Betti numbers, of a chain map, or of
the tree is used.  The values are properties of the TARGET's vertex count alone. QED

Corollary: the claim's entire evidence block is the statement 'the target has one
vertex', restated.  This is the corpus's SECOND vacuous control (W-03 convicted S4's
Control 1 of the same defect).  Demonstrated on three targets that share nothing:""")
# (a) a carrier that never formed anyway -- nothing to kill
Kn = K1(p=dict(v0=1.0, v1=0.0, v2=0.0, v3=0.0, v4=0.0), name="K1, p on the root only")
row("  (a) K1, all weight on v0", Kn)
row("  (a) collapsed", collapse_tree(Kn, T))
# (b) a big unrelated carrier: an n-gon wedge n-gon at one vertex, n=100
def wedge(n, m, nameo="wedge"):
    verts = ["w"] + [f"a{i}" for i in range(1, n)] + [f"b{i}" for i in range(1, m)]
    A = ["w"] + [f"a{i}" for i in range(1, n)] + ["w"]
    B = ["w"] + [f"b{i}" for i in range(1, m)] + ["w"]
    edges, lF, lC = [], [], []
    for i in range(n):
        edges.append((f"A{i}", A[i], A[i+1])); lF.append((f"A{i}", +1))
    for i in range(m):
        edges.append((f"B{i}", B[i], B[i+1])); lC.append((f"B{i}", +1))
    return Carrier(nameo, verts, edges, [], lF, lC, {v: 1.0 for v in verts})
W = wedge(100, 100, "100-gon wedge 100-gon (199 verts)")
row("  (b) 199-vertex wedge", W)
Wt = W
for e in [f"A{i}" for i in range(99)] + [f"B{i}" for i in range(99)]:
    Wt = contract(Wt, e)
Wt.name = "collapsed"
row("  (b) collapsed to 1 vertex", Wt)
print("""
  -> identical output (FORMS False, rkG 0, lA 0, lB 0) on a carrier that NEVER
     formed and on a 199-vertex carrier with no relation to K1.  A control whose
     output is fixed by the target's vertex count tests nothing about the map.""")

# ============================================================ ATTACK 2
hr("R1.2  ATTACK 2 -- THE WALL IS IN THE WRONG PLACE.  HOMOTOPY EQUIVALENCES SIT ON BOTH SIDES")
print("Every row below is a contraction of ONE edge: a quotient by a contractible")
print("subcomplex, hence a homotopy equivalence, with the claim's own certificate")
print("(chain map defect d1.d2 = 0, Betti preserved).  SENSE U weights.\n")
base = row("K1  (baseline, SENSE U)", KU)
for e in ["e1", "e2", "e3", "e4", "e5", "e6"]:
    Kc = contract(KU, e)
    row(f"  contract {e}", Kc)
print("""
READ THE TABLE.
  * contract e2 (merges v1,v2 -- BOTH in class (1,0)) and contract e5 (merges v3,v4
    -- both in class (0,1)) are HOMOTOPY EQUIVALENCES under which pi, lambda_A and
    lambda_B are EXACTLY unchanged, bit for bit.
  * contract e1/e3/e4/e6 are HOMOTOPY EQUIVALENCES under which lambda_A MOVES, while
    every diagnostic the claim offers is unchanged: Betti (1,1,0), chain-map defect 0,
    |S| = 3, rank G = 2, FORMS = True.  Formation is fully alive and lambda has moved.
  => 'homotopy equivalence' is not the wall.  Homotopy equivalences fall on both sides.""")

hr("R1.3  THE MAXIMAL CLASS-PRESERVING COLLAPSE -- the wall is at 3 vertices, not 1")
K32 = contract(contract(KU, "e2"), "e5", "K1 / {e2,e5}")
row("K1 / {e2,e5}  (3 vertices)", K32)
print("  pi identical to K1's to the bit:", np.array_equal(K32.pi(), KU.pi()))
print("  lambda_A identical to the bit  :",
      lambda_A(K32.pi(), 2.0, 1.1) == lambda_A(KU.pi(), 2.0, 1.1))
print("""  This is the LARGEST subcomplex of the spanning tree whose collapse preserves
  everything.  The tree T = {e1,e2,e4,e5} contains {e2,e5}; collapsing the OTHER two
  tree edges is what kills formation.  The claim's wall ('the spanning tree') is two
  edges too late: the destructive content of the tree collapse is e1 and e4, and each
  of those, taken alone, leaves FORMS = True.""")

# ============================================================ ATTACK 3
hr("R1.4  ATTACK 3 -- THE NORMALISATION IS DOING THE WORK (SENSE U vs SENSE C)")
KC = K1(p=dict(v0=0.4, v1=0.15, v2=0.15, v3=0.15, v4=0.15), name="K1 (S3's own p)")
print("Same six one-edge contractions, but with the CORPUS'S OWN published ready state")
print("p = (0.4,0.15,0.15,0.15,0.15) (S3 sec.3.3 / S4 sec.3) instead of uniform:\n")
row("K1  (baseline, SENSE C)", KC)
for e in ["e1", "e2", "e3", "e4", "e5", "e6"]:
    row(f"  contract {e}", contract(KC, e))
print("""
  Under the corpus's OWN weights, contracting a single edge moves lambda_B as well as
  lambda_A -- e.g. e1: -0.767507880 -> %.12f, a shift of %.3e -- with FORMS = True,
  |S| = 3, rank G = 2, Betti (1,1,0), chain-map defect 0.  The claim reports lambda_B
  as surviving every move until the full collapse ONLY because uniform weights make
  K1's four class weights {0,2/5,2/5,1/5} and the contracted {0,1/5,2/5,2/5} the SAME
  MULTISET, and lambda_B is a multiset function (W-03, of record).  The claim's silent
  choice of uniform weights is exactly what hides the earlier wall.""" % (
      cassaigne_maillot(*[contract(KC,'e1').pi()[i] for i in (3,1,2)]),
      abs(cassaigne_maillot(*[contract(KC,'e1').pi()[i] for i in (3,1,2)]) - (-0.767507880358))))

# ============================================================ ATTACK 4
hr("R1.5  ATTACK 4 -- TRAFFIC THE OTHER WAY: DESTROY THE HOMOTOPY TYPE, KEEP lambda EXACTLY")
KU2 = Carrier("K1 + second face (both triangles filled)",
              KU.verts, KU.edges,
              [("F", [("e1", +1), ("e2", +1), ("e3", +1)]),
               ("F2", [("e4", +1), ("e5", +1), ("e6", +1)])],
              KU.loopF, KU.loopC, {v: 1.0 for v in KU.verts})
row("K1  (b1=1, homotopy type S^1)", KU)
row("K1 + F2  (b1=0, CONTRACTIBLE)", KU2)
print("""  chi 0 -> +1, b1 1 -> 0: the homotopy type is destroyed outright (S^1 -> point),
  the flat holonomy is no longer flat -- and lambda_A and lambda_B do not move by one
  bit.  So NON-homotopy-equivalences preserve lambda exactly while homotopy
  equivalences (R1.2) change it.  A wall with traffic in both directions is not a wall.
  (This is S4's registered CONTROL 1, already ruled VACUOUS BY IDENTITY at W-03.)""")

hr("R1.6  ATTACK 5 -- lambda = 0 AND G = {1} ARE NOT SIGNATURES OF THE COLLAPSE")
print("Uncollapsed K1, connection W_F = W_C = 1 (f = c = 0):")
print("  |Z_k| for k=1..10 :", [round(abs(KU.Z_matrix(1+0j, 1+0j, k)), 15) for k in range(1, 11)])
print("  lambda_A = %.12f   lambda_B = %.12f" % (lambda_A(piU, 0.0, 0.0),
                                                 lambda_B_direct(piU, 0.0, 0.0, 200000)))
print("  characters chi_(1,0)=chi_(0,1)=chi_(1,1)=1  ->  G = {1}  ->  FORMS = False")
print("""  The claim's ENTIRE signature -- |Z_k| = 1 for every k, lambda = 0 exactly,
  G = {1}, FORMS False -- is produced on the UNCOLLAPSED five-vertex carrier by
  moving the connection to one point of T^2.  The signature reports degeneracy of
  the CHARACTER SET, not anything about the complex or the map.""")

hr("R1.7  ATTACK 6 -- THERE IS NO MORPHISM OF CARRIERS, SO THERE IS NO 'STOPS AT'")
print("""A carrier is (complex, rank-one Hermitian bundle, U(1) connection, ready state).
The tree collapse psi: K1 -> K1/T is a map of COMPLEXES.  For it to be a map of
CARRIERS it must carry the state space Gamma(L) = C^V:
   push  C^5 -> C^1 : no linear isometry exists (rank).                 max rank 1 < 5
   pull  C^1 -> C^5 : psi^*s has norm sqrt(5)|s(w)|, not |s(w)|.        NOT an isometry
So the ready state does not transport; the lane must CHOOSE a pushforward of p, and it
chose summation.  That choice is an IMPORT the corpus does not define.  W-03 already
records of record: 'No cellular collapse, quotient, or map connects any two of the ten
complexes.  The family is a LIST.'  The claim's wall is drawn across a morphism that,
by the corpus's own registered finding, does not exist.""")

hr("R1.8  THE ACTUAL WALL, STATED AND VERIFIED")
print("""THEOREM R-2.  Let psi: K -> K' be any cellular map with psi(gamma_F) = gamma'_F,
psi(gamma_C) = gamma'_C and psi_* p = p'.  Then pi' = pi implies Z'_k = Z_k for every k
and every connection with matched holonomies, hence lambda' = lambda EXACTLY under
EVERY schedule.  Conversely lambda moves only when pi moves.  Formation dies iff |S|
drops to 1, i.e. iff psi merges ALL occupied classes.
The wall is: 'psi fails to be injective on the occupied loop-incidence classes.'
It is not homotopy, not homeomorphism, not Betti, not the chain map, not the tree.""")
ok = True
for e in ["e2", "e5"]:
    Kc = contract(KU, e)
    for _ in range(200):
        WF, WC = np.exp(1j*rng.uniform(0, 2*np.pi)), np.exp(1j*rng.uniform(0, 2*np.pi))
        k = int(rng.integers(1, 50))
        ok &= abs(Kc.Z_matrix(WF, WC, k) - KU.Z_matrix(WF, WC, k)) < 1e-12
print("  Z_k(K1) == Z_k(K1/e2) and Z_k(K1/e5) at 400 random (W_F,W_C,k):", ok)
print("  (seed 31415926, k drawn uniform on 1..49)")

hr("R1.9  SUBDIVISION -- THE WALL IS EARLIER STILL (a HOMEOMORPHISM moves lambda)")
Ks = Carrier("K1, e2 subdivided (same space)",
             KU.verts + ["w"],
             [("e1", "v0", "v1"), ("e2a", "v1", "w"), ("e2b", "w", "v2"),
              ("e3", "v2", "v0"), ("e4", "v0", "v3"), ("e5", "v3", "v4"), ("e6", "v4", "v0")],
             [("F", [("e1", +1), ("e2a", +1), ("e2b", +1), ("e3", +1)])],
             [("e1", +1), ("e2a", +1), ("e2b", +1), ("e3", +1)],
             [("e4", +1), ("e5", +1), ("e6", +1)],
             {v: 1.0 for v in KU.verts + ["w"]})
row("K1            (SENSE U)", KU)
row("K1 subdivided (SENSE U)", Ks)
print("""  Same topological space, same homotopy type, same everything -- and under the
  claim's own SENSE U normalisation lambda_A and lambda_B both move.  The earliest
  wall is not even at a homotopy equivalence: it is at a REFINEMENT OF THE CELL
  STRUCTURE that changes nothing at all.  (S4 CONTROL 3 already recorded this.)""")
