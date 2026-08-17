# W10-D leg 4 -- the claims that are about the CARRIER'S ARITHMETIC AND ITS COMPLEX, not about the
# class multiset: W-02's minimality/cost, W-05's LEG TWO and LEG THREE, W-06's dressed restoration
# and its wedge route, W-03's "topology is inert", W-08's schedule exponent.
import numpy as np
from itertools import permutations

rng = np.random.default_rng(20260816)
CLS = ('00', '10', '01', '11')
EXP = {'00': (0, 0), '10': (1, 0), '01': (0, 1), '11': (1, 1)}

def Zk(p, f, c, k):
    return sum(p[i]*np.exp(1j*k*(-EXP[CLS[i]][0]*f + EXP[CLS[i]][1]*c)) for i in range(4))

print("="*100)
print("== 4A  FIX A DEGENERACY IN LEG 3D BEFORE USING IT ==")
print("="*100)
print("  Leg 3D counted 4 permutations preserving |Z_k| pointwise on B0b.  B0b has p10 = p11 =")
print("  2/9, so two of those four are the SAME weight vector as the other two -- an accidental")
print("  degeneracy of B0b's multiset, not a symmetry.  Recorded, not silently patched.  Re-run")
print("  on a GENERIC four-class weight vector with no repeated entry:")
pg = np.array([0.31, 0.17, 0.29, 0.23]); fv, cv = 1.3, 2.0
base = np.abs(Zk(pg, fv, cv, np.arange(1, 9)))
surv = [q for q in permutations(range(4))
        if np.abs(np.abs(Zk(pg[list(q)], fv, cv, np.arange(1, 9)))-base).max() < 1e-12]
print(f"    generic p = {pg}:  |Z_k|-preserving permutations = {len(surv)} of 24 -> {surv}")
print("    (0,1,2,3) = identity;  (3,2,1,0) = W-03's involution 00<->11, 10<->01.")
print("    TWO of 24 at the |Z_k| level; TWENTY-FOUR of 24 at the lambda level.  Leg 3D's")
print("    conclusion is unchanged and its count is corrected from 4 to 2.")

print("\n"+"="*100)
print("== 4B  THE ARITHMETIC CLAIMS.  V = 5 IS PRIME.  V = 9 AND V = 6 ARE NOT. ==")
print("="*100)
print("  S2 Theorem 2 (via W-01): no tensor factorisation of the state space because dim = 5 is")
print("  prime.  S3 CHOICE LEDGER C1 rejects the adjunction-free alternatives (b),(c) on 5 | 9")
print("  and 5 | 7 failing.  W-05 LEG THREE: 'an accident of five vertices being prime'.")
print("  W-06's rebuild route is the wedge sequence V = 4k+1.  ALL FOUR are statements about 5.")
print(f"  {'carrier':16s} {'V':>3s} {'V prime?':>9s} {'factorisations of C^V':>26s} {'V = 4k+1?':>10s}")
for nm, V in [("B1  K1", 5), ("B1p bridged", 6), ("B4  spindle", 6), ("B1q spectator", 7),
              ("B0b torus", 9), ("B1s subdivided", 11)]:
    fac = [f"C^{d} (x) C^{V//d}" for d in range(2, V) if V % d == 0 and d <= V//d]
    isp = all(V % d for d in range(2, int(V**.5)+1))
    print(f"  {nm:16s} {V:3d} {str(isp):>9s} {(', '.join(fac) if fac else 'NONE (prime)'):>26s} "
          f"{str((V-1) % 4 == 0):>10s}")
print("\n  ON BOTH FOUR-CLASS CARRIERS THE PRIME ARGUMENT IS UNAVAILABLE:")
print("    B0b: C^9 = C^3 (x) C^3 exists, and 9 | 81 = 9^2, so M_9^(x)N is a directed system with")
print("         UNITAL embeddings and NO adjunction at all -- S3's C1(b)/(c) rejection ground")
print("         ('5 does not divide 9, 5 does not divide 7') has no analogue on B0b.")
print("    B4 : C^6 = C^2 (x) C^3.  A genuine tensor factorisation of the WHOLE space exists.")
print("  S2 THEOREM 2, S3's C1 REJECTIONS, W-02's 'cost: one qubit per cell, PROVED MINIMAL',")
print("  and W-06's V = 4k+1 WEDGE ROUTE ARE ALL K1-SCOPED, AND THEY FAIL ON BOTH FOUR-CLASS")
print("  CARRIERS THE CORPUS OWNS -- not because the theorems are wrong but because their")
print("  hypothesis (V prime, V = 4k+1) is a property of the number 5.")
print("  What is NOT K1-scoped: the FLOOR dim_C >= 4 from 'two non-parallel branch vectors'.")
print("  That needs only two occupied classes with distinct characters -- checked next.")

print("\n"+"="*100)
print("== 4C  W-05's LEG TWO ON A FOUR-CLASS CARRIER: IS THE SLOT ALREADY INSIDE? ==")
print("="*100)
print("  W-05: the compression of M_V(C) to span{M_dF s, M_c s} has dim_C = 4, both branches pure,")
print("  overlap = Z.  Ingredients: two branch vectors that are NOT parallel.  Off K1:")
CARR = {'B1  K1  (3cl)': (5, np.array([0, 2/5, 2/5, 1/5])),
        'B1q spec (3cl)': (7, np.array([1/7, 3/7, 3/7, 0])),
        'B1p brdg (2cl)': (6, np.array([0, 1/2, 1/2, 0])),
        'B0b torus(4cl)': (9, np.array([4/9, 2/9, 1/9, 2/9])),
        'B4  spin (4cl)': (6, np.array([1/6, 1/6, 1/6, 3/6]))}
print(f"  {'carrier':16s} {'|<a,b>|':>10s} {'||a||':>8s} {'||b||':>8s} {'dim span':>9s} {'dim_C':>6s} {'Z_1':>10s}")
for nm, (V, p) in CARR.items():
    f0, c0 = 1.3, 2.0
    ch = np.array([np.exp(1j*(-EXP[cl][0]*f0 + EXP[cl][1]*c0)) for cl in CLS])
    amp = np.sqrt(p)
    a = amp*np.array([ch[i] if EXP[CLS[i]][0] else 1 for i in range(4)])   # M_dF s
    b = amp*np.array([ch[i] if EXP[CLS[i]][1] else 1 for i in range(4)])   # M_c s
    # collapse to the class-level vectors; the compression's rank is the rank of [a;b]
    ov = np.vdot(a, b)
    r = np.linalg.matrix_rank(np.stack([a, b]), tol=1e-12)
    print(f"  {nm:16s} {abs(ov):10.6f} {np.linalg.norm(a):8.6f} {np.linalg.norm(b):8.6f} "
          f"{r:9d} {r*r:6d} {abs(Zk(p,f0,c0,1)):10.6f}")
print("  dim_C = (dim span)^2 = 4 wherever the two branch vectors are independent, which happens")
print("  exactly when two occupied classes have distinct characters -- i.e. exactly when G != {1}.")
print("  W-05's LEG TWO IS CARRIER-INDEPENDENT, conditioned on G != {1}.  On B1p (2 classes,")
print("  {10,01}) it still holds; the class that kills it is |S| = 1, on any carrier.")

print("\n"+"="*100)
print("== 4D  W-03's 'THE CARRIER'S TOPOLOGY IS INERT', EXHIBITED ON THE FOUR-CLASS PAIR ==")
print("="*100)
print("  B0b: chi = 0, b1 = 2, b2 = 1, V = 9, F = 9.   B4: chi = 2, b1 = 1, b2 = 2, V = 6, F = 4.")
print("  Every topological invariant differs.  At SENSE C (class weights 1/4 each) both give")
print("  P = (1+x)(1+y)/4 and lambda = log(1/4) EXACTLY:")
pc = np.array([.25, .25, .25, .25])
for nm in ('B0b', 'B4'):
    k = np.arange(1, 400001)
    Z = Zk(pc, 1.3, 2.0, k)
    print(f"    {nm}  SENSE-C direct (1/N)log|Omega_N| at N=4e5: {np.log(np.abs(Z)).mean():+.9f}"
          f"   log(1/4) = {np.log(0.25):+.9f}")
print("  Identical to the last printed digit BECAUSE THE INPUT IS IDENTICAL -- both rows feed the")
print("  SAME weight 4-vector to the SAME functional.  THIS IS A ZERO-VARIABLE CONTROL AND IS")
print("  REPORTED AS ONE.  It confirms nothing about topology; it re-exhibits W-03's own ruling")
print("  that d2 enters the functional NOWHERE, which is an analytic fact, not an experiment.")
print("  The only non-vacuous version: B0b and B4 at their OWN SENSE-U weights give DIFFERENT")
print("  rates (log(4/9) vs log(1/2)) -- and that difference is entirely the class multiset,")
print("  which is incidence, not topology.  'TOPOLOGY IS INERT' IS CARRIER-INDEPENDENT AND")
print("  UNTESTABLE BY CONSTRUCTION: no experiment on any carrier can bear on it.")

print("\n"+"="*100)
print("== 4E  W-08's SCHEDULE EXPONENT ON FOUR CLASSES ==")
print("="*100)
print("  The adversary writes only the sqrt(K) cells of SMALLEST 1-|Z_k| and accumulates O(1).")
print("  Same code path, same K sweep; the ONE thing that moves is the weight 4-vector.")
print(f"  {'carrier':16s} {'K=1e4':>10s} {'K=1e5':>10s} {'K=1e6':>10s} {'K=1e7':>10s}  {'honest k_n=n at 1e6':>20s}")
for nm, (V, p) in CARR.items():
    row = []
    for K in (10**4, 10**5, 10**6, 10**7):
        kk = np.arange(1, K+1)
        d = 1.0 - np.abs(Zk(p, 1.3, 2.0, kk))
        m = int(round(np.sqrt(K)))
        idx = np.argpartition(d, m)[:m]
        row.append(-np.log(np.abs(Zk(p, 1.3, 2.0, kk[idx]))).sum())
    kk = np.arange(1, 10**6+1)
    hon = -np.log(np.abs(Zk(p, 1.3, 2.0, kk))).sum()
    print(f"  {nm:16s} " + " ".join(f"{r:10.4f}" for r in row) + f"  {hon:20.1f}")
print("  The adversarial accumulation is O(1) and K-INDEPENDENT on every carrier including both")
print("  four-class ones; the honest schedule diverges linearly.  W-08's schedule statement is")
print("  CARRIER-INDEPENDENT.  Its EXPONENT (K^{-1/2} at d_eff = 2, K^{-1/3} at d_eff = 1)")
print("  depends on d_eff = rank of the relation lattice, which is 2 for every occupied set of")
print("  size >= 3 (leg 3A) -- so it is the same on three- and four-class carriers.")
