import numpy as np, math, sys, itertools
sys.path.insert(0,'/private/tmp/claude-501/-Users-bgm-MB-Work/549fc5c6-445d-490b-bf95-ac9313727b33/scratchpad')
from rx import *
rng = np.random.default_rng(RNG_SEED)

def show(w):
    return "(%.6f,%.6f,%.6f,%.6f)"%(w[(0,0)],w[(1,0)],w[(0,1)],w[(1,1)])

print("="*78); print("A1  IS lambda_B A FUNCTION OF THE UNORDERED MULTISET?  (S4/W-03 says yes)")
print("="*78)
for wt in [(0.10,0.20,0.30,0.40), (0.05,0.11,0.37,0.47), (0.0,0.3,0.3,0.4), (0.25,)*4]:
    vals = [mahler2(*perm) for perm in set(itertools.permutations(wt))]
    print(f"  weights {wt}  distinct perms {len(vals):2d}  spread = {max(vals)-min(vals):.3e}   lambda = {vals[0]:.12f}")
print("""
  PROOF (mine, not the corpus's).  m(a1+a2X+a3Y+a4XY) = (1/2pi) int log max(|a1+a2 e^{it}|,
  |a3+a4 e^{it}|) dt  (Jensen in Y).  For REAL a_i, |a1+a2 e^{it}| = |a2+a1 e^{it}| POINTWISE
  (both = sqrt(a1^2+a2^2+2 a1 a2 cos t)), so the transpositions (a1 a2) and (a3 a4) are
  pointwise identities.  X<->Y is a monomial substitution of the torus, giving (a2 a3).
  <(12),(23),(34)> = S_4.   Hence lambda_B is a SYMMETRIC function of the four class weights. QED
  => the loop-incidence LABELS are invisible to lambda: only the MULTISET survives.""")

print()
print("="*78)
print("A2  A CLASS-MERGING ELEMENTARY EDGE COLLAPSE THAT PRESERVES lambda EXACTLY")
print("="*78)
K, L = K1(), K1_collapse_e4()
print("  phi : K1 -> K1/e4  collapses e4 (v0->v3), identifying v0 [class (1,1)] with v3 [class (0,1)].")
print("  It is an elementary edge collapse and a homotopy equivalence.  CLASS-MERGING.")
print("  L4 incidence d1 =\n", L.d1())
print("  L4 betti (b0,b1,b2,chi) =", L.betti(), " classes:", L.classes())

def push(K, L, phi, p):
    q = {v: 0.0 for v in L.V}
    for v in K.V: q[phi[v]] += p[v]
    return q

phi4 = {'v0':'w','v3':'w','v1':'v1','v2':'v2','v4':'v4'}

cases = [("S3/S4 state  p=(0.4,.15,.15,.15,.15)", {'v0':0.4,'v1':0.15,'v2':0.15,'v3':0.15,'v4':0.15}),
         ("SENSE U      p=(0.2,0.2,0.2,0.2,0.2)", {v:0.2 for v in K.V}),
         ("MY STATE     p=(0.25,.15,.15,.20,.25)", {'v0':0.25,'v1':0.15,'v2':0.15,'v3':0.20,'v4':0.25})]
for lbl, p in cases:
    wK = class_weights(K,p); wL = class_weights(L, push(K,L,phi4,p))
    lK, lL = lambdaB_exact(wK), lambdaB_exact(wL)
    fK, rK, SK = formation(wK); fL, rL, SL = formation(wL)
    print(f"\n  {lbl}")
    print(f"    source class wts {show(wK)}   ->  target {show(wL)}")
    print(f"    lambda_B  source {lK:.12f}   target {lL:.12f}   GAP = {abs(lK-lL):.3e}")
    print(f"    formation {fK} rank G {rK}    ->  {fL} rank G {rL}")
print("""
  READ IT.  The SAME class-merging elementary edge collapse has gap 1.70e-01 on the corpus's
  own S3/S4 state and EXACTLY 0.000e+00 on SENSE U -- the uniform-on-vertices state that produced
  the claim's own source number -0.756573585640.  The gap is a property of the (map, STATE) PAIR,
  not of the map.  The exactness is not numerical luck: the
  collapse moves mass p(v3) from class (0,1) to class (1,1), and p(v0)=p(v4) makes the new pair
  (p01',p11') = (p01,p11) TRANSPOSED, which A1 proves lambda cannot see.""")

# a whole codim-1 family, not one point
print("  A ONE-PARAMETER FAMILY of exactly-preserving states (all weights strictly positive):")
for t in (0.05,0.10,0.18,0.30):
    z = (1-2*t)/3.0
    p = {'v0':z,'v1':t,'v2':t,'v3':1-2*t-2*z,'v4':z}
    wK = class_weights(K,p); wL = class_weights(L,push(K,L,phi4,p))
    print(f"    p={tuple(round(p[v],6) for v in K.V)}  gap = {abs(lambdaB_exact(wK)-lambdaB_exact(wL)):.3e}")

print()
print("="*78)
print("A3  A CLASS-COMPATIBLE MAP THAT DOES *NOT* PRESERVE lambda  (SENSE U)")
print("="*78)
Ks, Kt = K1_subdivided(), K1()
# collapse the forest {e1b, e2b, e3a, e4b, e5b, e6a}: m1->v1, m2->v2, m3->v2, m4->v3, m5->v4, m6->v4
phis = {'v0':'v0','v1':'v1','v2':'v2','v3':'v3','v4':'v4',
        'm1':'v1','m2':'v2','m3':'v2','m4':'v3','m5':'v4','m6':'v4'}
cs, ct = Ks.classes(), Kt.classes()
ok = all(cs[v]==ct[phis[v]] for v in Ks.V)
print("  phi : B1s -> K1  collapses the six subdivision half-edges.  CLASS-COMPATIBLE:", ok)
print("  It is a homotopy equivalence AND the two carriers are HOMEOMORPHIC spaces.")
pu_s = {v:1.0/len(Ks.V) for v in Ks.V}; pu_t = {v:1.0/len(Kt.V) for v in Kt.V}
ws, wt_ = class_weights(Ks,pu_s), class_weights(Kt,pu_t)
print(f"    SENSE U  source {show(ws)}  lambda = {lambdaB_exact(ws):.12f}")
print(f"    SENSE U  target {show(wt_)}  lambda = {lambdaB_exact(wt_):.12f}")
print(f"    GAP = {abs(lambdaB_exact(ws)-lambdaB_exact(wt_)):.3e}      (S4 Control 3: 3.181e-02)")
wpf = class_weights(Kt, push(Ks,Kt,phis,pu_s))
print(f"    pushforward of the source state: {show(wpf)}  lambda = {lambdaB_exact(wpf):.12f}  GAP = 0 BY IDENTITY")
print("""
  So 'class-compatible => lambda preserved EXACTLY' holds only under the PUSHFORWARD (=SENSE C)
  normalisation, where the class weights are preserved BY CONSTRUCTION and Z_k^K = Z_k^L is an
  IDENTITY, not a measurement.  Under SENSE U -- S4's other sense, recorded in its CHOICE LEDGER
  as C2 'OPEN and load-bearing' -- the same class-compatible map moves lambda by 3.18e-02.
  And the claim's own source number -0.756573585640 IS the SENSE U value of K1.""")
