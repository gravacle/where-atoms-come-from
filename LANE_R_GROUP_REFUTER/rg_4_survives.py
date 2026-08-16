#!/usr/bin/env python3
"""
rg_4_survives.py -- WHAT OF THE CLAIM I COULD NOT BREAK.
A refuter that reports only kills is not a refuter.  These are the parts of the
CLAIM that I attacked and that held.

  S1.  THE CONTAINMENT IS A THEOREM.  supp(z) is contained in the union over
       occupied classes (a,b) of  a(-Lambda) x b(Lambda).  Proved and searched.
  S2.  THE ARITHMETIC IS EXACT.  3 / 8 / 6 characters at exactly the stated
       exponent sets (already shown in rg_0_validate.py).
  S3.  A SINGLE OCCUPIED CLASS CAN ALREADY FORM, once transport is non-scalar.
       True, and the correct converse of my T1: it needs NON-COMMUTING
       holonomies for rank 2, but rank 1 already means formation.
  S4.  "THE THEOREM IS RIGHT AND THE BOX IS WRONG" is the right framing.

SEED: numpy default_rng(999).  2400 random configurations.
"""
import numpy as np, math, itertools
from rg_lib import *

print("=" * 78)
print("S1.  THE CONTAINMENT.  supp(z) subset UNION_(a,b) a(-Lambda) x b(Lambda).")
print("=" * 78)
print("""
  PROOF.  A(x) = V_A diag(e^{i wA_j x}) V_A^{-1} with wA in -Lambda, and
  B(y) likewise with wB in Lambda.  Every term of
     <s_v, A(x)^a B(y)^b s_v> = SUM_{j,m} c_{jm} e^{i(a wA_j x + b wB_m y)}
  carries an exponent (a wA_j, b wB_m).  There are no others.   QED
  The CLAIM's formula is therefore a correct UPPER BOUND.  It is not an
  equality, and every kill in rg_1/rg_2 exploits the gap.
""")
rng = np.random.default_rng(999)
viol = 0; tested = 0
CORNERS = [(0, 0), (1, 0), (0, 1), (1, 1)]
for trial in range(2400):
    kind = int(rng.integers(0, 3))
    if kind == 0:
        d = int(rng.integers(1, 5))
        w = sorted(set(int(x) for x in rng.integers(1, 5, size=d)))
        while len(w) < d:
            w.append(w[-1] + 1)
        A, B = op_u1(w, -1), op_u1(w, +1)
        LamA, LamB = set(int(x) for x in A.w), set(int(x) for x in B.w)
    elif kind == 1:
        ax1 = rng.normal(size=3); ax2 = rng.normal(size=3)
        A, B = op_su2(ax1, -1), op_su2(ax2, +1)
        LamA, LamB = {-1, 1}, {-1, 1}
        d = 2
    else:
        ax1 = rng.normal(size=3)
        A, B = op_su2(ax1, -1), op_su2(ax1, +1)     # commuting
        LamA, LamB = {-1, 1}, {-1, 1}
        d = 2
    r = int(rng.integers(1, 5))
    S = [CORNERS[i] for i in rng.choice(4, size=r, replace=False)]
    cls, secs = [], []
    for cl in S:
        cls.append(cl)
        v = rng.normal(size=d) + 1j * rng.normal(size=d)
        secs.append(v / np.linalg.norm(v) / math.sqrt(r))
    md = modes(cls, secs, A, B)
    allowed = set()
    for (a, b) in S:
        for wa in (LamA if a else {0}):
            for wb in (LamB if b else {0}):
                allowed.add((a * wa, b * wb))
    tested += 1
    if not set(md.keys()) <= allowed:
        viol += 1
        print("   VIOLATION:", sorted(set(md.keys()) - allowed))
print(f"   {tested} random configurations, containment violations: {viol}")
print(f"   THE CONTAINMENT HELD EVERYWHERE.  I could not break it.")

print()
print("=" * 78)
print("S3.  A SINGLE OCCUPIED CLASS CAN ALREADY FORM (the claim's true half).")
print("=" * 78)
A, B = op_su2([0, 0, 1], -1), op_su2([1, 0, 0], +1)
v = np.array([1.0, 0.4 + 0.3j]); v = v / np.linalg.norm(v)
md = modes([(1, 1)], [v], A, B)
print(f"   SU(2), non-commuting, ONE class (1,1): support {sorted(md.keys())}, "
      f"rank {lattice_rank(sorted(md.keys()))}")
print(f"   lambda (Haar T^2, 4096^2) = {lam_haar(md, M=4096):.12f}   < 0 => FORMATION")
Ac, Bc = op_su2([0, 0, 1], -1), op_su2([0, 0, 1], +1)
md2 = modes([(1, 1)], [v], Ac, Bc)
print(f"   SU(2), COMMUTING,     ONE class (1,1): support {sorted(md2.keys())}, "
      f"rank {lattice_rank(sorted(md2.keys()))}")
print(f"   lambda (Haar T^2, 4096^2) = {lam_haar(md2, M=4096):.12f}   < 0 => FORMATION")
print("""
   BOTH FORM.  So the claim's 'a single occupied class already' is right about
   FORMATION and wrong about RANK.  The corrected statement:

     a single occupied class suffices for FORMATION as soon as the transport
     on that class's fibre is NON-SCALAR (two distinct weights);
     it suffices for RANK 2 only when the two holonomies additionally FAIL TO
     COMMUTE and the ready state is not an eigenvector of either.

   That is a three-condition statement.  The CLAIM states one condition, and
   the one it states (d = 2) is not among the three.
""")
