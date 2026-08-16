"""R-CHARGE REFUTER -- the attack.  Every section names the clause of the claim it hits."""
import itertools
import numpy as np
import rc_lib as R

Ecls = R.class_vectors(5, R.K1_EDGES, R.K1_GAMMA_F, R.K1_GAMMA_C)
p_S3 = np.array([0.4, 0.15, 0.15, 0.15, 0.15])          # S3's ready state, pi = (.4,.3,.3)
PI = "pi = (p_00,p_10,p_01,p_11) = (0, 0.3, 0.3, 0.4)"


def Echarge(q):
    return (Ecls.T * np.asarray(q, dtype=np.int64)).T


print("=" * 100)
print("A1.  CLAUSE 'CHARGE IS THE FIRST MODALITY THAT VARIES SOMETHING OTHER THAN pi'")
print("     Two objects already in the corpus vary lambda at pi HELD EXACTLY FIXED.")
print("=" * 100)
print("     ", PI, "   unit charge, carrier K1, held fixed in every row below.")
pts = [(1.0, np.sqrt(2)), (2.0, 1.1), (np.pi, 3 * np.pi / 2), (0.4, 2.6), (1.5, 1.5)]
print("      (f,c)                 lambda_A = log|Z_1|      lambda_B (generic/exact)")
for (f, c) in pts:
    z = abs(R.Zk(Ecls, p_S3, f, c, 1))
    la = np.log(z) if z > 0 else float("-inf")
    print(f"      ({f:.6f},{c:.6f})   {la:>18.9f}      {R.lambda_B_generic(Ecls, p_S3, Nx=8192):>14.9f}")
print("      THE CONNECTION moves lambda_A over the whole of [-inf, 0] at fixed pi.")
print("      THE SCHEDULE moves lambda at fixed pi AND fixed connection:")
f, c = 1.0, np.sqrt(2)
laA = np.log(abs(R.Zk(Ecls, p_S3, f, c, 1)))
laB = R.lambda_B_generic(Ecls, p_S3, Nx=16384)
print(f"        at (1,sqrt2):  schedule A {laA:.9f}   schedule B {laB:.9f}   gap {abs(laA-laB):.9f}")
print("      Neither the connection nor the schedule is pi.  Both are of record BEFORE charge:")
print("      S4 section 2.2 defines both schedules; W-03's corrected headline names TWO")
print("      determinants, 'the pushforward pi AND (under the canonical clock) the relation")
print("      lattice'.  W-03 also records that the critic RAN winding, which varies the")
print("      circuit index, not pi.  CHARGE IS AT BEST THIRD.")

print()
print("=" * 100)
print("A2.  CLAUSE 'IT DESTROYS THE FOUR-CLASS TAXONOMY OUTRIGHT'")
print("     (a) homogeneous charge is a SCHEDULE REPARAMETRISATION, not a new object")
print("=" * 100)
rng = np.random.default_rng(90210002)
worst = 0.0
for _ in range(500):
    q = int(rng.integers(1, 7))
    f, c = rng.uniform(0, 2 * np.pi, 2)
    k = int(rng.integers(1, 30))
    worst = max(worst, abs(R.Zk(Echarge([q] * 5), p_S3, f, c, k) - R.Zk(Ecls, p_S3, f, c, q * k)))
print(f"      Z_k(homogeneous charge q) == Z_{{qk}}(unit charge):  max dev {worst:.3e}"
      f"   (500 samples, seed 90210002)")
print("      So homogeneous charge q IS the clock k -> qk.  It cannot be a modality that")
print("      'varies something other than pi' unless the SCHEDULE already was one.")
print()
print("     (b) the taxonomy survives EVERY class-homogeneous charge, exactly")
base = R.lambda_B_generic(Ecls, p_S3, Nx=16384)
rows = []
for (q0, q1, q3) in [(1, 1, 1), (2, 2, 2), (3, 3, 3), (2, 1, 1), (1, 2, 3), (3, 1, 2), (-1, 2, 5)]:
    q = [q0, q1, q1, q3, q3]
    lam = R.lambda_B_generic(Echarge(q), p_S3, Nx=16384)
    rows.append((q, lam, abs(lam - base)))
for q, lam, d in rows:
    print(f"      q = {str(q):<20} lambda_B^gen = {lam:.12f}   |dev from unit charge| = {d:.2e}")
print("      Class-homogeneous charge NEVER moves the generic rate on K1.  The 'four-class'")
print("      object is exactly what survives; what breaks is charge INHOMOGENEOUS WITHIN a")
print("      class -- which is not a property of charge, it is a refusal to be a class function.")
print()
print("     (c) how much of the taxonomy actually dies: enumeration over q in {0,1,2}^5")
cnt_rank = {0: 0, 1: 0, 2: 0}
same = 0
vals = {}
for q in itertools.product(range(3), repeat=5):
    Eq = Echarge(q)
    bas = R.delta_lattice(Eq, p_S3)
    cnt_rank[len(bas)] += 1
    lam = R.lambda_B_generic(Eq, p_S3, Nx=1024)
    if abs(lam - base) < 1e-6:
        same += 1
    vals[round(lam, 6)] = vals.get(round(lam, 6), 0) + 1
print(f"      243 charge assignments: rank Delta = 2 in {cnt_rank[2]}, = 1 in {cnt_rank[1]},"
      f" = 0 in {cnt_rank[0]}")
print(f"      lambda_B^gen equal to the unit-charge value in {same} of 243;"
      f" {len(vals)} distinct values")
print("      S4-1 has TWO implications.  Only ONE of them fails under charge:")
print("        SURVIVES at every exponent map E:  |S| = 1  =>  Delta = 0  =>  no formation, ever")
print("        SURVIVES at every exponent map E:  |S| = 2  =>  rank Delta <= 1  (two points have")
print("                                           one difference; nothing about corners is used)")
print("        FAILS under inhomogeneous charge:  |S| >= 3  =>  rank 2")
enum = {0: 0, 1: 0, 2: 0}
rng2 = np.random.default_rng(90210003)
bad2 = 0
for _ in range(4000):
    k = int(rng2.integers(1, 4))
    Erand = rng2.integers(-4, 5, size=(k, 2))
    pr = np.ones(k) / k
    b = R.delta_lattice(Erand, pr)
    enum[len(b)] += 1
    if k <= 2 and len(b) > 1:
        bad2 += 1
print(f"      4000 random exponent sets with |S| <= 3 (seed 90210003): rank counts {enum};"
      f" |S|<=2 with rank 2: {bad2}")
print("      'OUTRIGHT' IS FALSE.  One implication of S4-1 dies; two live; and the criterion")
print("      itself (W-02's G != 1) is untouched by charge -- as the claim's own next clause says.")

print()
print("=" * 100)
print("A3.  CLAUSE 'THE CORRECT CRITERION AT ARBITRARY CHARGE IS Delta(S) not-contained-in L'")
print("     (a) it is not a criterion for FORMATION in W-01's sense -- the sense the")
print("         project's own first ruling fixed.  It is necessary, never sufficient.")
print("=" * 100)
p_heavy = np.array([0.8, 0.05, 0.05, 0.05, 0.05])
mn = min(abs(R.Zk(Ecls, p_heavy, 1.0, np.sqrt(2), k)) for k in range(1, 20001))
bas = R.delta_lattice(Ecls, p_heavy)
print(f"      p = (0.8,0.05,0.05,0.05,0.05), unit charge, (f,c)=(1,sqrt2)")
print(f"        Delta basis = {bas}  -> Delta = Z^2, L = 0, so Delta not-contained-in L: "
      f"criterion says FORMATION")
print(f"        min_k |Z_k| over k <= 20000 = {mn:.12f}   (the overlap NEVER vanishes)")
print(f"        lambda_B^gen = {R.lambda_B_generic(Ecls, p_heavy, Nx=16384):.12f}  (< 0: it DECAYS)")
print("      W-01's ruling of record is a CONVEX-HULL condition -- 0 in the hull of the")
print("      characters -- and it is a condition on the WEIGHTS (it needs max class weight")
print("      <= 1/2).  No containment of lattices can carry a weight condition.  The lattice")
print("      statement decides lambda < 0; it does not decide firing.  The claim states an iff.")
print()
print("     (b) and the containment is CLOCK-NORMALISED -- unstated in the claim.")
A_, B_, M_ = 2, 3, 4                      # f = 2pi*2/4 = pi, c = 2pi*3/4 = 3pi/2  (S1 section 6)
f1, c1 = 2 * np.pi * A_ / M_, 2 * np.pi * B_ / M_
bas = R.delta_lattice(Ecls, p_S3)
print(f"      S1 section 6's OWN published connection: W_F = -1, W_C = -i, i.e."
      f" (f,c) = (2pi*{A_}/{M_}, 2pi*{B_}/{M_})")
print(f"        L = {{(m,n) : -{A_}m + {B_}n = 0 mod {M_}}};  Delta = Z^2 basis {bas}")
print(f"        Delta contained in L ? {R.delta_subset_L(bas, A_, B_, M_)}  -> criterion says FORMATION")
for M in (1, 2, 4):
    zs = [abs(R.Zk(Ecls, p_S3, f1, c1, M * n)) for n in range(1, 501)]
    om = float(np.sum(np.log(np.maximum(zs, 1e-300))))
    print(f"        schedule k_n = {M}n :  min|Z| = {min(zs):.12f}   sum(1-|Z|) = "
          f"{sum(1-z for z in zs):.3e}   log|Omega_500| = {om:.6f}")
print("      Under the clock k_n = 4n the record NEVER FORMS on this connection, while")
print("      Delta not-contained-in L holds.  The correct clock-free statement is")
print("      'formation <=> k Delta not-contained-in L for infinitely many k of the schedule'.")
print("      W-03's corrected headline says 'under the canonical clock' and the claim drops it.")
print("      S4 FLAG F1 / CHOICE LEDGER C1: which clock is physical is OPEN.")

print()
print("=" * 100)
print("A4.  CLAUSE 'WHICH UNIFIES S4-1 AND THE RELATION-LATTICE THEOREM INTO ONE STATEMENT'")
print("     G = Delta/(Delta ^ L) DOES NOT DETERMINE lambda_B.  Same Delta, same L, same")
print("     weights, two values of lambda 0.29 apart -- both reachable by per-vertex charge.")
print("=" * 100)
cases = [
    ("q = (1,2,2,2,2)", [1, 2, 2, 2, 2]),
    ("q = (2,3,3,6,6)", [2, 3, 3, 6, 6]),
    ("q = (2,6,6,3,3)", [2, 6, 6, 3, 3]),
    ("q = (3,4,4,12,12)", [3, 4, 4, 12, 12]),
]
for name, q in cases:
    Eq = Echarge(q)
    bas = R.delta_lattice(Eq, p_S3)
    lam = R.lambda_B_generic(Eq, p_S3, Nx=16384)
    print(f"      {name:<18} E = {[tuple(int(t) for t in x) for x in Eq[[0,1,3]]]}"
          f"  Delta basis {bas}  rank {len(bas)}  lambda_B^gen = {lam:.12f}")
print(f"      log(0.3) = {np.log(0.3):.12f}    log(0.4) = {np.log(0.4):.12f}")
print("      Every row above has rank Delta = 1, generic L = 0, hence G = Delta ~ Z in every")
print("      row -- the SAME group, the SAME quotient, the SAME weights (0.4,0.3,0.3).")
print("      lambda_B is not constant across them.  A statement whose content is the quotient")
print("      G CANNOT contain the relation-lattice theorem, which is a statement about the")
print("      VALUE of lambda_B.  The containment is binary; the rate needs the embedded")
print("      exponent multiset, not the abstract quotient.  'One statement' is two statements.")

print()
print("=" * 100)
print("A5.  CLAUSE 'd2 STILL ENTERS NOWHERE, SO CONTROL 1 STAYS VACUOUS BY IDENTITY AT EVERY")
print("     CHARGE'  --  A CONTROL THAT COULD NOT HAVE FAILED, RE-RUN 1024 TIMES.")
print("=" * 100)
print("      Arguments of the functional Z_k: (p_v), (E_v) = (q_v a_v, q_v b_v), (f,c), k.")
print("      Control 1 (fill the second triangle) changes: the face list, hence d2, hence chi,")
print("      b1, b2.  It changes NO vertex, NO edge, NO loop, hence no a_v, no b_v, no p_v,")
print("      no q_v, no (f,c).  Every argument is fixed, so Z_k is fixed pointwise in k.")
print("      Verified as an identity rather than a statistic:")
rng3 = np.random.default_rng(90210005)
worst = 0.0
for _ in range(1024):
    q = rng3.integers(-4, 5, size=5)
    f, c = rng3.uniform(0, 2 * np.pi, 2)
    k = int(rng3.integers(1, 50))
    # 'filled' and 'unfilled' differ ONLY in the face list; the exponent data is built from
    # the 1-skeleton and the designated loops, which are identical.
    E_unfilled = Echarge(q)
    E_filled = Echarge(q)
    worst = max(worst, abs(R.Zk(E_filled, p_S3, f, c, k) - R.Zk(E_unfilled, p_S3, f, c, k)))
print(f"        max |Z_k(filled) - Z_k(unfilled)| over 1024 charges = {worst:.1e}"
      f"   (seed 90210005)")
print("      This is 0 because the two sides are the SAME EXPRESSION, not because charge was")
print("      tested.  A control whose two arms are literally the same function of the same")
print("      arguments cannot fail at any charge, any connection, any schedule or any carrier.")
print("      Reporting its survival as a result is the exact defect W-03 convicted S4 of.")
print()
print("      AND THE UNIVERSAL IS OVER THE WRONG INDEX SET.  Nothing indexed by charge can")
print("      make d2 enter, because d2 has no slot in Z_k.  The object that gives d2 a slot is")
print("      an ACTION on the 2-cells -- W-04's RECOMMENDATION OF RECORD, test (1).")
print("      ILLUSTRATION (not actual surface; the action is an IMPORT, see IMPORT AUDIT):")
print("      Gibbs weight exp(beta * sum over FILLED faces of Re W), K1 (one filled) vs B2")
print("      (both filled), same pi, same charge, ensemble average of lambda_A:")
NG = 2048
th = (np.arange(NG) + 0.5) * 2 * np.pi / NG
FF, CC = np.meshgrid(th, th, indexing="ij")
for q in ([1, 1, 1, 1, 1], [1, 2, 2, 2, 2], [1, 1, 2, 1, 1]):
    Eq = Echarge(q)
    lam = np.zeros_like(FF)
    for j in range(5):
        pass
    Zg = sum(p_S3[j] * np.exp(1j * (-Eq[j, 0] * FF + Eq[j, 1] * CC)) for j in range(5))
    lamA = np.log(np.abs(Zg) + 1e-300)
    for beta in (0.0, 1.0, 2.0):
        w1 = np.exp(beta * np.cos(FF))                       # K1: face F only
        w2 = np.exp(beta * (np.cos(FF) + np.cos(CC)))        # B2: both triangles filled
        e1 = float(np.sum(w1 * lamA) / np.sum(w1))
        e2 = float(np.sum(w2 * lamA) / np.sum(w2))
        print(f"        q={str(q):<16} beta={beta:<4} E[lambda_A] K1 = {e1:.9f}   B2 = {e2:.9f}"
              f"   |diff| = {abs(e1-e2):.3e}")
print("      At beta = 0 the fill is invisible; at beta > 0 it is not, and the size of the")
print("      gap depends on the charge.  'd2 enters nowhere' is a statement about a")
print("      ZERO-ACTION construction, not a statement that survives 'at every charge'.")
