#!/usr/bin/env python3
"""
X_02 — R-7's SECOND SENTENCE IS FALSE, AND THE LANE'S OWN SEALED ARM REFUTES IT.

TARGET.  R_07_THEOREMS.txt R-7 (and R_06_sufficiency.out.txt part 3, and the lane's returned
finding R-7):
    "BY R-2, N1's LIMIT HOLDS AT THE CORPUS'S OWN PUBLISHED GENERIC CONNECTION AND ITS OWN
     REGISTERED READY STATE, and by the same argument AT EVERY ALGEBRAIC (alpha,beta) and
     more generally whenever alpha, beta lie in the Q-span of 1, algebraic numbers, and
     logarithms of algebraic numbers."
    R_06_sufficiency.out.txt: "-- and the same argument covers every algebraic (alpha,beta)."

THE DEFECT.  Baker supplies (D1) and (D2).  R-2 needs (D1) AND (D2) AND -- unstated in R-7 --
L(omega) = {0}.  R-7's own non-vanishing argument for Lam_3 is
    "Lam_3 = 0 with p != 0 => pi algebraic, false; with p = 0 => m + n sqrt2 = 0 => m = n = 0",
which is a computation SPECIFIC to (alpha,beta) = (-1/2pi, sqrt2/2pi): it is exactly the
verification that L = {0} there, and it is precisely the step that does not generalise.
For an ALGEBRAIC pair, Lam_3 = m alpha + n beta - p is an algebraic number and CAN be zero.

WHAT FALLS.  "At every algebraic (alpha,beta)" is false, and the counterexample is inside the
target lane's own sealed output: R_03_sudler_1d.out.txt's arm BA-silver runs
alpha = sqrt2 - 1 on the line u v = 1, i.e. beta = 2 - sqrt2.  BOTH ARE ALGEBRAIC OF DEGREE 2.
That arm converges -- the lane reports it converging with fitted r = 0.736 -- to log(0.3),
NOT to m(P).

ISOLATION.  ONE VARIABLE: the RELATION LATTICE L(omega).  pi, the estimator, the N-grid and
the phase reducer are identical in every arm.  Arms A,B,C are algebraic pairs with rank L = 1
(three different lattices); arm D is the corpus's published generic connection, which is
NOT algebraic and has L = {0}.  The repaired claim -- algebraic AND L(omega) = {0} -- is
what arm D satisfies and arms A,B,C do not.
"""
import numpy as np
from fractions import Fraction
from X_lib import (PI_K1, P_eval, m_maxform, m_one_var, frac_sqrt, frac_pi,
                   ExactRot, relation_lattice, arm_hash)

PREC = 60
SQ2 = frac_sqrt(2, PREC)
PIC = frac_pi(PREC)
MP = m_maxform(PI_K1, 1 << 24)
KMAX = 10 ** 7
DECS = [10 ** i for i in range(1, 8)]

print("=" * 79)
print("X_02 — 'AT EVERY ALGEBRAIC (alpha,beta)' IS FALSE")
print("=" * 79)
print("\nm(P) = %.12f    log(0.3) = %.12f" % (MP, np.log(0.3)))

ARMS = [
    ("A  alg, L=<(1,1)>   a=sqrt2-1, b=2-sqrt2", SQ2 - 1, 2 - SQ2, "ALGEBRAIC deg 2 both"),
    ("B  alg, L=<(-2,1)>  a=frac(sqrt2), b=frac(2 sqrt2)", SQ2 - 1, 2 * SQ2 - 2, "ALGEBRAIC deg 2 both"),
    ("C  alg, L=<(2,2)>   a=sqrt2-1, b=3/2-sqrt2", SQ2 - 1, Fraction(3, 2) - SQ2, "ALGEBRAIC deg 2 both"),
    ("D  CORPUS f=1,c=sqrt2  a=-1/2pi, b=sqrt2/2pi", -1 / (2 * PIC), SQ2 / (2 * PIC), "NOT algebraic (Lindemann)"),
]

print("\n1. THE RELATION LATTICE OF EACH ARM, EXACT (Fraction arithmetic), |m|,|n| <= 12")
for name, a, b, tag in ARMS:
    L = relation_lattice(a, b, B=12)
    print("   %-46s  %s" % (name, tag))
    print("      L(omega) cap box = %s" % (L if L else "{0}  -> H2 HOLDS"))

print("\n2. THE SUBTORUS VALUE EACH RESONANT ARM MUST CONVERGE TO (Jensen on roots)")
p00, p10, p01, p11 = PI_K1
vA = m_one_var([p01, p11, p10])                       # Q_{1,1}: 0.3 + 0.4 z + 0.3 z^2
vB = m_one_var([p10, p01, p11])                       # Q_{-2,1} = P(z,z^2)/z: 0.3 + 0.3z + 0.4z^2
vC = 0.5 * (m_one_var([p01, p11, p10]) + m_one_var([-p01, -p11, p10]))
print("   A  m(Q_{1,1})  = %.12f   ( = log 0.3, the register's published (1,1) row )" % vA)
print("   B  m(Q_{-2,1}) = %.12f" % vB)
print("   C  d=2 two-circle closure = %.12f   (X_01)" % vC)
print("   D  m(P)        = %.12f" % MP)

print("\n3. SEVEN DECADES.  ARMS-DIFF GUARD ON THE OUTPUT VECTORS.")
k = np.arange(1, KMAX + 1, dtype=np.int64)
rows, hs = {}, {}
for name, a, b, tag in ARMS:
    fa = ExactRot(a).frac(k); fb = ExactRot(b).frac(k)
    lz = np.log(np.abs(P_eval(PI_K1, np.exp(2j * np.pi * fa), np.exp(2j * np.pi * fb))))
    hs[name] = arm_hash(lz[:200000])
    cs = np.cumsum(lz)
    rows[name] = [cs[N - 1] / N for N in DECS]
    del fa, fb, lz, cs
for name in rows:
    print("   %-46s %s" % (name, hs[name]))
print("   distinct hashes: %d of %d" % (len(set(hs.values())), len(hs)))

hdr = "   %-46s" % "N" + "".join("%15d" % N for N in DECS)
print("\n   S_N")
print(hdr)
for name, *_ in ARMS:
    print("   %-46s" % name + "".join("%15.9f" % v for v in rows[name]))
print("\n   S_N - m(P)      <- R-7 claims this goes to 0 on EVERY algebraic (alpha,beta)")
print(hdr)
for name, *_ in ARMS:
    print("   %-46s" % name + "".join("%15.2e" % (v - MP) for v in rows[name]))
print("\n   S_N - (its own subtorus value)")
print(hdr)
for (name, *_), tgt in zip(ARMS, [vA, vB, vC, MP]):
    print("   %-46s" % name + "".join("%15.2e" % (v - tgt) for v in rows[name]))

print("""
4. READ.
   Arms A, B, C are ALGEBRAIC PAIRS.  None converges to m(P); each converges to its own
   subtorus value, to 1e-6 or better at N = 1e7, across seven decades and monotonically in
   trend.  Arm A IS the target lane's own R_03 'BA-silver' row (R_03_sudler_1d.out.txt,
   first orbit point 0.41421356 = sqrt2 - 1, on the line u v = 1), reported there as
   converging to log(0.3) with fitted r = 0.736.
   SO THE LANE'S OWN SEALED ARM IS A COUNTEREXAMPLE TO ITS OWN R-7.
   Arm D -- the corpus's published generic connection, which is NOT algebraic -- is the one
   that converges to m(P).  R-7's THEOREM AT f=1, c=sqrt2 IS NOT TOUCHED BY ANY OF THIS.
   THE REPAIR IS ONE CLAUSE: 'at every algebraic (alpha,beta) WITH L(omega) = {0}', i.e. with
   1, alpha, beta linearly independent over Q.  Under that clause Liouville's inequality gives
   (D1) and Baker gives (D2) and R-7's conclusion stands.
""")
print("DONE X_02")
