# X4 — THREE OF LANE C's LEGS ARE IDENTITIES, AND ONLY ONE OF THE THREE IS MARKED AS ONE.
# Lane C's own rule (its PUBLISHED_CONVENTIONS, "WHAT IS A THEOREM AND WHAT IS EVIDENCE") marks
# legs 1, 5 and 7A as [T].  X4 shows legs 4b, 7B and the leg-3A "identity check" belong on that
# list too -- and that leg 4b's identity is about the GAUGE GROUP, not about COR-F's T at all.
import numpy as np
from x_lib import *

rng = np.random.default_rng(4242)

print("== X4a  LEG 4b: 'T^m is gauge-INVARIANT iff L | m' IS THE COMMUTANT OF THE DIAGONAL GROUP ==")
print("   Lane C's leg 4 verifies, with 200 random gauge transforms at each of 14 ticks on two")
print("   carriers, that {m : T^m gauge-invariant} = {m : L | m}, and calls the coincidence with")
print("   leg 1's invisibility sublattice 'the stipulation, correctly identified -- NEW'.")
print("   But T(g.a) = G T(a) G* (leg 4a, its own number), so")
print("        T^m(g.a) = G T^m(a) G*  for EVERY m,")
print("   and X is fixed by conjugation by ALL diagonal unitaries  <=>  X IS DIAGONAL.")
print("   That is a fact about U(1)^V, true of EVERY matrix, and it needs no T.  Exhibited:")
worst_dia, worst_off = 0.0, np.inf
for trial in range(400):
    n = rng.integers(3, 8)
    X = rng.normal(size=(n,n)) + 1j*rng.normal(size=(n,n))
    Xd = np.diag(np.diag(X))
    dev_d = dev_o = 0.0
    for _ in range(40):
        th = rng.uniform(0, 2*np.pi, n); G = np.diag(np.exp(1j*th))
        dev_d = max(dev_d, np.linalg.norm(G@Xd@G.conj().T - Xd))
        dev_o = max(dev_o, np.linalg.norm(G@X @G.conj().T - X ))
    worst_dia = max(worst_dia, dev_d); worst_off = min(worst_off, dev_o)
print(f"   400 random matrices: max ||G X_diag G* - X_diag|| = {worst_dia:.2e}   (invariant)")
print(f"                        min over non-diagonal X of max ||G X G* - X|| = {worst_off:.3e}")
print("   => 'gauge-INVARIANT OPERATOR' and 'DIAGONAL OPERATOR' ARE THE SAME PREDICATE.")
print("   Leg 4b therefore re-derives leg 1a in different words; the 'same set, two descriptions'")
print("   line is literally right and is not a second, independent identification.")
print("   AND THE REGISTER ALREADY HOLDS THE NAME: W-05 leg one -- 'M_gamma is LITERALLY AN ELEMENT")
print("   OF THE GAUGE GROUP U(1)^V' -- and W-06's corrected N4 mechanism -- 'not scalar")
print("   multiplication: FIBRE-WISE-NESS.  Any fibre-wise unitary lies in the local gauge group.'")

print("\n== X4b  LEG 3A's 'IDENTITY CHECK' IS TRUE OF ANY SEQUENCE OF NUMBERS WHATEVER ==")
print("   Lane C files as a NEW DEFECT: \"the registrar's 'EDGE rescaled x3' ROW IS A UNIT ERROR,")
print("   exhibited as an identity: 3 x (mean over all n) = sum of the three residue means,")
print("   0.00e+00 at N=3e5\".  That is the definition of an arithmetic mean split by residue.")
for trial in range(3):
    x = rng.normal(size=300000)                       # ARBITRARY numbers, no physics at all
    m = [x[r::3].mean() for r in range(3)]
    print(f"   random sequence #{trial+1}:  3*mean(all) - sum(residue means) = "
          f"{abs(3*x.mean() - sum(m)):.2e}")
print("   The identity carries NO information about the registrar's row.  What the row actually")
print("   claims -- 'the EDGE per-tick rate is not m(P)/3' -- is ARITHMETICALLY TRUE and lane C")
print("   reproduces it (its own -0.864991741 vs m(P)/3 = -0.255836).  The charge that survives is")
print("   a charge of PRESUPPOSITION ('x3 presupposes one factor per edge tick'), which is a")
print("   framing objection, NOT a unit error.  Filing it under new_defects as 'A UNIT ERROR'")
print("   overstates it: no unit is misconverted anywhere in w11_d_rate.py.")
mP = m_jensen(np.array([0.,.30,.30,.40]))
print(f"   for the record: m(P) = {mP:.9f},  m(P)/3 = {mP/3:.9f},  registrar's edge/tick = -0.864939")

print("\n== X4c  LEG 3B's HEADER SAYS 'ONE VARIABLE: THE CARRIER'.  DIFF THE ARMS. ==")
K, B = K1(), B0b()
wB = np.array([.10,.12,.09,.14,.11,.11,.11,.11,.11]); wB /= wB.sum()
armK = dict(carrier="K1", nv=K.nv, LF=K.LF, LC=K.LC, pi=np.array([0.,.30,.30,.40]),
            nedges=len(K.edges))
armB = dict(carrier="B0b", nv=B.nv, LF=B.LF, LC=B.LC, pi=pi_of(B, np.sqrt(wB)+0j),
            nedges=len(B.edges))
aK = generic_conn(K, np.random.default_rng(7+K.nv)); aB = generic_conn(B, np.random.default_rng(7+B.nv))
print(f"   arm 1: {armK}")
print(f"   arm 2: {armB}")
print(f"   connection arm 1 (6 phases)  = {np.round(aK,6)}")
print(f"   connection arm 2 (18 phases) = {np.round(aB,6)}")
print(f"   m(P) arm 1 = {m_jensen(armK['pi']):.9f}     m(P) arm 2 = {m_jensen(armB['pi']):.9f}")
moved = [k for k in ("carrier","nv","LF","LC","nedges") if armK[k]!=armB[k]]
print(f"   THINGS THAT MOVE BETWEEN THE TWO ARMS: {moved} + pi + the connection + m(P)")
print("   SIX+ variables move.  The header 'ONE VARIABLE: the carrier (K1 vs B0b)' is false as")
print("   written; pi is a FREE choice, not forced by the carrier (only its zero pattern is).")
print("   Lane C's self_flag reports the 30x spread difference as unscored and 'not chased' --")
print("   correct -- but the leg is labelled as an isolated one-variable comparison and is not.")
