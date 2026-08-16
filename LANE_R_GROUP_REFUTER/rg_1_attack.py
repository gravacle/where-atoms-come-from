#!/usr/bin/env python3
"""
rg_1_attack.py -- THE ATTACK ON "S4-1 STOPS AT d = 2".

THREE INDEPENDENT KILLS.

  A1.  d = 2 IS NOT SUFFICIENT.  A rank-TWO fibre with a SCALAR representation
       (Lambda = {1}, i.e. chi (+) chi) reproduces S4-1 verbatim, bit for bit.
       So the wall is not at fibre rank.

  A2.  d = 2 IS NOT NECESSARY -- S4-1 ALREADY FAILS AT d = 1.  Non-uniform
       charge on rank-ONE fibres with a SCALAR abelian U(1) connection breaks
       BOTH halves of S4-1: |S| = 3 with rank G = 1, and (worse) |S| = 1 with
       formation occurring, against S4-1's "|S| = 1 -> never".

  A3.  THE PROPOSED REPLACEMENT IS FALSE IN THE CASE THE CLAIM NAMES.
       "A single occupied class already generates rank 2" holds only on the
       NON-COMMUTING locus.  For abelian G at any d -- and for SU(2) itself
       with commuting holonomies -- a single occupied class gives rank 1.

GRIDS / SEEDS: Haar M = 4096 midpoint; direct N = 400000; rng seed 20260816.
"""
import numpy as np, math, itertools
from rg_lib import *

CLS = k1_classes()

def report(label, md, note=""):
    sup = sorted(md.keys())
    rk = lattice_rank(sup)
    print(f"  {label}")
    print(f"     support ({len(sup)}): {sup}")
    print(f"     rank G = {rk}   {note}")
    return sup, rk

print("=" * 78)
print("A1.  d = 2 IS NOT SUFFICIENT.  RANK-TWO FIBRE, SCALAR REPRESENTATION.")
print("     Structure group U(1); fibre C^2; representation chi (+) chi,")
print("     weight set Lambda = {1,1}, i.e. rho(W) = e^{i theta} I_2.")
print("=" * 78)

A_d1, B_d1 = op_u1([1], sign=-1), op_u1([1], sign=+1)
A_sc, B_sc = op_u1([1, 1], sign=-1), op_u1([1, 1], sign=+1)

pvC = {(1, 1): 0.4, (1, 0): 0.3, (0, 1): 0.3}
md_d1 = modes(CLS, sections_from_class_weights(pvC, d=1), A_d1, B_d1)
rng = np.random.default_rng(20260816)
# random fibre DIRECTIONS at d=2 -- the thing that breaks section 2's cancellation
md_sc = modes(CLS, sections_from_class_weights(pvC, d=2, rng=rng), A_sc, B_sc)

s1, r1 = report("d = 1, U(1) unit charge, SENSE C", md_d1)
s2, r2 = report("d = 2, SCALAR rep chi(+)chi, SENSE C, RANDOM fibre directions", md_sc)
print(f"     supports identical: {s1 == s2}")
l1 = lam_haar(md_d1, M=4096); l2 = lam_haar(md_sc, M=4096)
print(f"     lambda(d=1) = {l1:.15f}")
print(f"     lambda(d=2) = {l2:.15f}      |difference| = {abs(l1-l2):.3e}")

print()
print("  THEOREM A1 (proved, not sampled).  If rho(W_F) = e^{if} I_d and")
print("  rho(W_C) = e^{ic} I_d then")
print("     <s_v, rho(W_F)^{-ka} rho(W_C)^{kb} s_v> = |s_v|^2 e^{ik(-af+bc)}")
print("  identically, so Z_k is the d=1 functional with p_v = |s_v|^2 and EVERY")
print("  conclusion of S4-1 holds verbatim at every d.  The fibre DIRECTION")
print("  cancels exactly -- S4 section 2's cancellation survives d=2 unharmed.")
print("  Exhaustive check of the S4-1 enumeration at d = 2, scalar rep:")
cnt = {0: 0, 1: 0, 2: 0}
corners = [(0, 0), (1, 0), (0, 1), (1, 1)]
for r in range(1, 5):
    for S in itertools.combinations(corners, r):
        # build a 5-vertex toy carrier realising exactly the classes in S
        cls = list(S)
        secs = []
        for _ in cls:
            v = rng.normal(size=2) + 1j * rng.normal(size=2)
            secs.append(v / np.linalg.norm(v) / math.sqrt(len(cls)))
        m = modes(cls, secs, A_sc, B_sc)
        cnt[lattice_rank(sorted(m.keys()))] += 1
print(f"     rank2 = {cnt[2]}  rank1 = {cnt[1]}  rank0 = {cnt[0]}"
      f"   [S4-1 predicts 5 / 6 / 4]   MATCH: {(cnt[2],cnt[1],cnt[0])==(5,6,4)}")

print()
print("=" * 78)
print("A2.  d = 2 IS NOT NECESSARY.  S4-1 FAILS AT d = 1 UNDER CHARGE.")
print("     Structure group U(1), rank-ONE fibres, SCALAR transport, ABELIAN.")
print("     Per-vertex charge q_v: the fibre at v carries the rep z -> z^{q_v}.")
print("=" * 78)

def charge_modes(pv_per_vertex, q):
    """d=1, per-vertex charge.  Returns modes dict.  Exact and closed form."""
    ops = [(op_u1([qq], sign=-1), op_u1([qq], sign=+1)) for qq in q]
    secs = [np.array([math.sqrt(p)], dtype=complex) for p in pv_per_vertex]
    return modes_general(CLS, secs, ops)

print("\n  A2a.  THE CHARGE VECTOR ALREADY OF RECORD, q = (1,2,2,2,2), SENSE C.")
p_senseC = [0.4, 0.15, 0.15, 0.15, 0.15]        # (0.4,0.3,0.3) over the classes
md_q = charge_modes(p_senseC, [1, 2, 2, 2, 2])
sup, rk = report("q = (1,2,2,2,2), classes {(1,1),(1,0),(0,1)} all occupied",
                 md_q, note="  <-- |S| = 3 OCCUPIED CLASSES")
print(f"     S4-1 predicts rank 2 for |S| = 3.  COMPUTED rank = {rk}.")
print(f"     S4-1 FALSIFIED AT d = 1: {rk != 2}")
lq = lam_haar(md_q, M=4096)
lq_d = lam_direct(md_q, 1.0, math.sqrt(2.0), N=400000)
print(f"     lambda (Haar T^2, 4096^2) = {lq:.12f}")
print(f"     lambda (direct N=4e5)     = {lq_d:.12f}")
print(f"     EXACT: Z_k = e^{{ik(c-f)}}(0.4 + 0.6 cos k(f+c)), so")
print(f"       lambda = m(0.3 + 0.4 z + 0.3 z^2) = log(0.3) = {math.log(0.3):.15f}")
print(f"       (both roots of 0.3z^2+0.4z+0.3 lie ON |z|=1: product = 1,")
print(f"        discriminant = {0.4**2-4*0.09:.3f} < 0, so Jensen gives log 0.3)")
print(f"     dev(Haar, exact)   = {abs(lq-math.log(0.3)):.3e}")
print(f"     dev(direct, exact) = {abs(lq_d-math.log(0.3)):.3e}")

print("\n  A2b.  THE STRONGER HALF: |S| = 1 AND FORMATION OCCURS ANYWAY.")
print("        All weight on ONE class, (1,0) = vertices v1,v2, charges 1 and 2.")
p_one = [0.0, 0.5, 0.5, 0.0, 0.0]
md_1 = charge_modes(p_one, [1, 1, 2, 1, 1])
sup1, rk1 = report("|S| = 1 (only class (1,0) occupied), charges q_v1=1, q_v2=2",
                   md_1, note="  <-- S4-1 says G = {1}, NO FORMATION EVER")
print(f"     G != trivial: {len(sup1) > 1}   -> FORMATION OCCURS.")
l_1 = lam_haar(md_1, M=4096)
print(f"     lambda (Haar) = {l_1:.12f}   EXACT Jensen: log max(0.5,0.5) = {math.log(0.5):.15f}")
print(f"     dev = {abs(l_1-math.log(0.5)):.3e}")
print("     S4-1's |S| = 1 clause -- 'there is no formation, ever' -- IS FALSE")
print("     at d = 1, rank-one fibre, U(1), scalar transport, abelian group.")

print()
print("=" * 78)
print("A3.  THE PROPOSED REPLACEMENT IS FALSE WHERE THE CLAIM EXTENDS IT.")
print("     CLAIM: 'A single occupied class already generates rank 2', and")
print("     'the character support is a(-Lambda_F) + b(Lambda_C)' (a product).")
print("=" * 78)

def single_class_test(label, A, B, d, direction=None, seed=7):
    r = np.random.default_rng(seed)
    if direction is None:
        v = r.normal(size=d) + 1j * r.normal(size=d)
    else:
        v = np.asarray(direction, dtype=complex)
    v = v / np.linalg.norm(v)
    md = modes([(1, 1)], [v], A, B)
    return report(label, md)

print("\n  A3a.  ABELIAN G = U(1), d = 2, weight set Lambda = {1,2}.")
Aab, Bab = op_u1([1, 2], sign=-1), op_u1([1, 2], sign=+1)
sup, rk = single_class_test("single occupied class (1,1)", Aab, Bab, 2)
print("     CLAIM's product set a(-Lambda)+b(Lambda) = "
      f"{sorted(itertools.product([-1,-2],[1,2]))}  -> rank 2")
print(f"     TRUE support = {sup}  -> rank {rk}.  CLAIM OVER-COUNTS BY 2 MODES.")
print(f"     'A single occupied class already generates rank 2' is FALSE here: {rk != 2}")

print("\n  A3b.  SU(2) -- THE CLAIM'S OWN GROUP -- WITH COMMUTING HOLONOMIES.")
Ac, Bc = op_su2([0, 0, 1], sign=-1), op_su2([0, 0, 1], sign=+1)
sup, rk = single_class_test("single occupied class (1,1), both axes = z", Ac, Bc, 2)
print("     CLAIM's product set {-1,+1} x {-1,+1} = 4 points -> rank 2")
print(f"     TRUE support = {sup}  -> rank {rk}.")

print("\n  A3c.  SU(2) NON-COMMUTING: the claim's headline IS true here.")
An, Bn = op_su2([0, 0, 1], sign=-1), op_su2([1, 0, 0], sign=+1)
sup, rk = single_class_test("single occupied class (1,1), axes z and x", An, Bn, 2)
print("     -> the headline is a statement about the NON-COMMUTING LOCUS only.")

print("\n  A3d.  AND EVEN THERE IT DEPENDS ON THE READY STATE, NOT THE GROUP.")
print("        Same non-commuting SU(2), single class, fibre direction = an")
print("        eigenvector of A (i.e. of rho(W_F)).")
ev = An.V[:, 0].copy()
sup, rk = single_class_test("s_v = eigenvector of rho(W_F)", An, Bn, 2, direction=ev)
print(f"     rank {rk}, not 2, at the SAME group, SAME representation, SAME class.")
print("     The CLAIM's replacement is indexed by (occupied classes, Lambda).")
print("     The true support is also a function of the READY STATE's fibre")
print("     DIRECTIONS -- a variable that does not exist at d = 1 and which")
print("     W-03 ALREADY RECORDED as moving the observable ('six ready states")
print("     with identical class weights give |Z_1| spread 0.4247').")
