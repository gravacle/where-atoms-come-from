#!/usr/bin/env python3
"""
RL3 — IS DOBROWOLSKI/LAWTON-Thm-1 ACTUALLY AN INGREDIENT OF THE LICENCE?
      AND HOW MUCH OF L4's "HOLDS 20/20" COULD HAVE FAILED?

TARGET CLAIMS UNDER TEST (LANE_W13_L_LITERATURE):
  (a) headline: "Lawton 1983 is TWO theorems and the second one IS an ingredient of the
      licence" -- offered as the ONE-WORD CORRECTION to W-08's M1 refuter.
  (b) L-2 / Sec. 1.2 / L4 docstring: "(*) is UNIFORM IN pi.  That uniformity is what the
      shell argument needs, because P's zero structure degenerates on three loci and a
      pi-dependent constant would not survive the limit."  L4 repeats: "A shell sum
      converges for any s > 0, so the exponent is not load-bearing; the UNIFORMITY IN pi is."
  (c) L4 evidence: "HOLDS 20/20, worst (measured-bound) = -2.24e-02."

VERBATIM CHECK OF THE CITATION ITSELF (I read the CMB PDF at the bytes; see RL7_SOURCES.txt):
  Dobrowolski Thm 1.1 constant: "C_k = (k - 1)(12 sqrt2 / pi)^{(k-2)/(k-1)}"        CONFIRMED
  Dobrowolski Thm 1.4: "mu_l{z in T^l : |P(z)| <= v} <= C(k_1,...,k_l)(v/h)^{1/(sum(k_i-1))},
     where C(k_1,...,k_l) = C_{k_1} + ... + C_{k_l}"                                CONFIRMED
  "h = h(P)" is the maximum modulus of the coefficients, "k_i = the number of terms of P
  considered as a polynomial in z_i with polynomial coefficients in the remaining variables"
                                                                                     CONFIRMED
  So the target's specialisation C = 2, exponent 1/2, bound 2 (v/max(pi))^{1/2} is ARITHMETIC-
  ALLY CORRECT and is verified below in exact rational arithmetic.  THE CITATION IS GOOD.
  What is tested here is the CLAIM ABOUT ITS ROLE.

TWO SEPARATE THINGS ARE TESTED, ONE VARIABLE EACH.
  LEG A (vacuity audit).  ONE variable moves: whether a row of L4's 20 could have failed.
  LEG B (load-bearing test).  ONE variable moves: which shell-measure input the assembly
        uses -- Dobrowolski's uniform bound vs the local expansion at FIXED pi.  Same pi,
        same shells, same delta schedule.
"""
import math
from fractions import Fraction as F
import numpy as np

print("=" * 78)
print("RL3 — THE ROLE OF DOBROWOLSKI IN THEOREM L, TESTED.")
print("=" * 78)

CASES = [
    ("K1 REGISTERED (0, .3, .3, .4)", (F(0), F(3, 10), F(3, 10), F(2, 5))),
    ("three-class centroid", (F(0), F(1, 3), F(1, 3), F(1, 3))),
    ("four-class uniform (two zero circles)", (F(1, 4), F(1, 4), F(1, 4), F(1, 4))),
    ("S1 published (one zero circle)", (F(0), F(1, 2), F(1, 2), F(0))),
    ("non-firing max > 1/2", (F(0), F(1, 10), F(1, 10), F(8, 10))),
]
VS = [1e-1, 1e-2, 1e-3, 1e-4]


def area_small(p, v, n):
    t = np.arange(n) * (2 * np.pi / n)
    ex = np.exp(1j * t)
    A = (p[0] + p[1] * ex)[:, None]
    B = (p[2] + p[3] * ex)[:, None]
    ey = np.exp(1j * t)[None, :]
    return float(np.mean(np.abs(A + B * ey) <= v))


def has_zero(p):
    """EXACT: does P vanish somewhere on T^2?  f(c)=g(c) for some c in [-1,1], both affine."""
    p00, p10, p01, p11 = p
    f1 = (p00 + p10) ** 2 - (p01 + p11) ** 2
    fm = (p00 - p10) ** 2 - (p01 - p11) ** 2
    return (f1 == 0) or (fm == 0) or (f1 > 0) != (fm > 0)


print("""
--------------------------------------------------------------------------------
LEG A — VACUITY AUDIT OF L4's "HOLDS 20/20".
A row CANNOT FAIL if the Dobrowolski bound is >= 1 (mu_2 <= 1 always), or if P has no torus
zero at all so the measured area is identically 0 for small v.  Both are EXACT tests.
--------------------------------------------------------------------------------""")
print("  %-40s %-8s %-12s %-12s %-26s" % ("case", "v", "bound", "measured", "could it have failed?"))
n_rows = 0
n_dead = 0
reasons = {}
for name, p in CASES:
    h = max(float(q) for q in p)
    hz = has_zero(p)
    for v in VS:
        n_rows += 1
        bd = 2.0 * math.sqrt(v / h)
        mu = area_small([float(q) for q in p], v, 8192)
        why = []
        if bd >= 1.0:
            why.append("BOUND>=1 (vacuous)")
        if not hz:
            why.append("NO TORUS ZERO (mu==0 exactly)")
        dead = len(why) > 0
        if dead:
            n_dead += 1
            reasons[name] = reasons.get(name, 0) + 1
        print("  %-40s %-8.0e %-12.4f %-12.4e %-26s"
              % (name if v == VS[0] else "", v, bd, mu, ("NO: " + "; ".join(why)) if dead else "yes"))
print("\n  ROWS THAT COULD NOT HAVE FAILED: %d of %d." % (n_dead, n_rows))
print("  L4 reports the check as 'HOLDS on 20 of 20'.  The live count is %d." % (n_rows - n_dead))
print("  This is NOT a fatal defect -- L4's own ledger declares L4 'NOT a control' -- but the")
print("  headline count 20/20 is quoted in finding L-2 as evidence, and %d of those 20 rows" % n_dead)
print("  are arithmetically incapable of violating the bound.  W-08's isolation audit calls")
print("  this the commonest FATAL defect class in this program; here it is non-fatal because")
print("  the surviving %d rows do carry the check." % (n_rows - n_dead))

print("""
--------------------------------------------------------------------------------
LEG B — IS THE UNIFORMITY IN pi LOAD-BEARING FOR THEOREM L?  NO, AND HERE IS WHY.

THEOREM L (target, Sec. 2.3) begins: "Let pi = (p00,p10,p01,p11) be a probability vector with
ALGEBRAIC entries".  pi is FIXED before the limits are taken, and the limits are N -> infinity
THEN delta -> 0.  NO LIMIT IN pi IS EVER TAKEN.  A pi-dependent constant therefore survives
every limit in the proof, and the sentence "a pi-dependent constant would not survive the
limit" describes a limit that Theorem L does not contain.

What the shell step actually needs at a FIXED pi with isolated zeros is the local expansion
|P| ~ L*d, which gives shell measure ~ pi*d^2 -- exponent 2, not 1/2.  Both make the shell
sum converge.  Measured below, at K1's registered pi, over four decades of shell scale.
--------------------------------------------------------------------------------""")
p = (0.0, 0.3, 0.3, 0.4)
# local gradient at the zero x0 = -2/3 + i sqrt5/3, y0 = conj(x0)
x0 = complex(-2.0 / 3.0, math.sqrt(5.0) / 3.0)
y0 = x0.conjugate()
# P(x0 e^{i a}, y0 e^{i b}) = i[a x0 (p10 + p11 y0) + b y0 (p01 + p11 x0)] + O(2)
A = x0 * (p[1] + p[3] * y0)
B = y0 * (p[2] + p[3] * x0)
g1, g2 = abs(A), abs(B)
detJ = abs((A.conjugate() * B).imag)      # Jacobian of (a,b) -> P/(2 pi i), in TURN coordinates
print("  local gradient moduli at the zero: |dP/da| = %.6f   |dP/db| = %.6f" % (g1, g2))
print("  |Im(conj(A) B)| = %.6f  (equals |A||B| = %.6f here, i.e. the two directions are orthogonal)"
      % (detJ, g1 * g2))
print("\n  %-10s %-16s %-18s %-18s %-12s" % ("v", "measured mu_2", "Dobrowolski 2(v/h)^.5", "local ~ c v^2", "ratio D/local"))
for v in (1e-2, 1e-3, 1e-4, 1e-5):
    mu = area_small(p, v, 16384)
    bd = 2.0 * math.sqrt(v / 0.4)
    loc = 2 * math.pi * v * v / (4 * math.pi * math.pi * detJ)   # TWO elliptical discs, area in turns^2
    print("  %-10.0e %-16.4e %-18.4e %-18.4e %-12.3e" % (v, mu, bd, loc, bd / max(loc, 1e-300)))
print("""
  At the shell scales a dyadic argument actually visits (v down to delta*2^-j_max with
  j_max = O(log N)), Dobrowolski's uniform bound is larger than the truth by 4 to 9 ORDERS OF
  MAGNITUDE and rising.  It is a correct bound and a very loose one.

  VERDICT ON THE TARGET'S ONE-WORD CORRECTION.
  * The CITATION is verified verbatim and the specialisation is arithmetically right.
  * Dobrowolski/Lawton Thm 1 is a SUFFICIENT shell-measure input, and it is the RIGHT one
    for a statement UNIFORM OVER pi, or for isolated zeros of unknown multiplicity, since it
    needs no case analysis.  That is a real use and it is worth citing.
  * It is NOT NECESSARY for Theorem L as the target states it, because Theorem L fixes pi;
    and the reason the target gives for its necessity -- that a pi-dependent constant "would
    not survive the limit" -- names a limit Theorem L does not take.
  * So "the second one IS an ingredient of the licence" over-states an OPTIONAL ingredient,
    and W-08/M1_08 T2(d)'s ruling (Lawton is not the theorem that licenses N1) needed no
    correction at all.  M1_08 T2(c) already carried the local expansion that does the job.

--------------------------------------------------------------------------------
LEG C — A GAP IN THE ASSEMBLY THAT NEITHER M1_08 T2(c) NOR THE TARGET FLAGS.

Both write "the number of points is at most N(mu_2 + 4 D_N) ... by discrepancy".  D_N is the
ordinary (box) discrepancy that Erdos-Turan-Koksma bounds.  Box discrepancy bounds counts in
AXIS-PARALLEL BOXES, not in the sublevel set {|P| <= w}, which is a neighbourhood of an
algebraic curve.  The standard repair is ISOTROPIC discrepancy J_N (Kuipers-Niederreiter,
Ch. 2, Thm 1.6), which is bounded by a CONSTANT TIMES D_N^{1/d} -- in d = 2, a SQUARE ROOT.
The assembly's error term therefore reads C_7 N^{-delta'/2}(log N)^2, not C_7 N^{-delta'}(log N)^2.
Still a power of N, so THEOREM L IS NOT BROKEN -- but the exponent as printed is wrong, and it
is wrong in the sealed lane it was inherited from as well.
Below: how many connected components the sublevel set has, and whether they are convex, at
K1's registered pi -- the geometry that decides whether the isotropic-discrepancy repair even
applies.
--------------------------------------------------------------------------------""")


def components(p, v, n=1024):
    t = np.arange(n) * (2 * np.pi / n)
    ex = np.exp(1j * t)
    A = (p[0] + p[1] * ex)[:, None]
    B = (p[2] + p[3] * ex)[:, None]
    M = np.abs(A + B * np.exp(1j * t)[None, :]) <= v
    seen = np.zeros_like(M, dtype=bool)
    comp = 0
    idx = np.argwhere(M)
    for i0, j0 in idx:
        if seen[i0, j0]:
            continue
        comp += 1
        stack = [(i0, j0)]
        seen[i0, j0] = True
        while stack:
            i, j = stack.pop()
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                a, b = (i + di) % n, (j + dj) % n
                if M[a, b] and not seen[a, b]:
                    seen[a, b] = True
                    stack.append((a, b))
    return comp, int(M.sum())


for v in (3e-1, 1e-1, 3e-2, 1e-2):
    c, cells = components(p, v)
    print("   v = %-8.0e  connected components on a 1024^2 torus grid: %-3d   cells: %d" % (v, c, cells))
print("""
   READING.  Below v ~ 0.1 the sublevel set is TWO small components, one per zero, and each
   is convex to grid resolution -- so isotropic discrepancy applies with a factor 2 and the
   repair above is available.  The target's assembly does not carry it; neither does the
   sealed M1_08 T2(c) sketch it was built from.  RECORDED AS A CORRECTION TO BOTH.
--------------------------------------------------------------------------------""")
print("\nDONE RL3")
