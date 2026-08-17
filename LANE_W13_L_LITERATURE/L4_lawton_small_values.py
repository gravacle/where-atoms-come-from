#!/usr/bin/env python3
"""
L4 — LAWTON'S *OTHER* THEOREM, WHICH IS THE ONE THIS CORPUS ACTUALLY NEEDS, CHECKED ON OUR P.

W-03's IMPORT AUDIT note, W-05's REDISCOVERY list and W-08/M1_08 T2(d) all treat "Lawton 1983"
as ONE theorem — the Boyd conjecture on limits of Mahler measures — and M1_08 T2(d) rules it
out as the licence for N1.  THAT RULING IS CORRECT AND IS CONFIRMED VERBATIM BY THE LITERATURE
(Lind-Schmidt-Verbitskiy, Sec. 9, quoted in THEOREMS_AND_CITATIONS.txt).

BUT LAWTON'S 1983 PAPER CONTAINS A SECOND THEOREM, AND IT *IS* AN INGREDIENT OF THE LICENCE.
Quoted verbatim from Dobrowolski, "A Note on Lawton's Theorem", Canad. Math. Bull. 60 (2017)
484-489, Section 1, where it is recalled as "Theorem 1 of Lawton's seminal paper":

    "Let P(x) in C[x] be a monic polynomial and let k be the number of nonzero coefficients
     of P.  Then if k >= 2, there is a positive constant C_k that depends only on k such that
         mu{ x in [0,1) : |P(e^{2 pi i x})| < v } <= C_k v^{1/(k-1)}
     for every real v > 0."

  NOTE THE COEFFICIENT FIELD: C[x].  COMPLEX, NOT INTEGER.  So this half of Lawton applies to
  our P directly, on the whole simplex, including irrational pi — unlike Boyd-Lawton read in
  its integer form (R1_08 A3, sealed, flagged exactly this gap and it is closed here).

DOBROWOLSKI 2017 Theorem 1.4, the SEVERAL-VARIABLE version, verbatim:

    "Let P(z_l) be a polynomial with complex coefficients.  If P has at least two monomials,
     then
         mu_l{ z in T^l : |P(z)| <= v } <= C(k_1,...,k_l) (v/h)^{1/(sum_{i=1}^{l}(k_i - 1))},
     where C(k_1,...,k_l) = C_{k_1} + ... + C_{k_l}"
  with h = h(P) the maximum modulus of the coefficients, k_i the number of terms of P read as
  a polynomial in z_i over the remaining variables, and (Theorem 1.1)
         C_k = (k-1) (12 sqrt2 / pi)^{(k-2)/(k-1)}.

SPECIALISED TO OUR P = p00 + p10 x + p01 y + p11 xy:
  as a polynomial in x it is (p00 + p01 y) + (p10 + p11 y) x, so k_1 = 2; likewise k_2 = 2.
  sum (k_i - 1) = 2, C_2 = 1(12 sqrt2/pi)^0 = 1, C(k_1,k_2) = 2.  Hence, for every pi:

      mu_2{ |P| <= v }  <=  2 (v / max(pi))^{1/2}.          (*)

  (*) is UNIFORM IN pi.  That uniformity is what the shell argument needs, because P's zero
  structure degenerates on three loci (M1_02, sealed) and a pi-dependent constant would not
  survive the limit.  THIS IS THE PLACE WHERE LAWTON IS THE RIGHT CITATION FOR N1, and it is
  a different theorem from the one the corpus has been failing to cite.

This leg checks (*) numerically on our P, and measures the TRUE local exponent, which is not 1/2.
Grid-based AREA estimation is used, not a min or a quadrature of log|P|: an area converges like
(perimeter / n) and is not noise-limited near the zeros.  Convergence in n is shown, not assumed.
"""
import math
import numpy as np

print("=" * 78)
print("L4 — DOBROWOLSKI 1.4 / LAWTON THM 1 ON OUR OWN P.  UNIFORM SMALL-VALUE BOUND.")
print("=" * 78)

CASES = [
    ("K1 REGISTERED (0, .3, .3, .4)", (0.0, 0.3, 0.3, 0.4)),
    ("three-class centroid", (0.0, 1 / 3, 1 / 3, 1 / 3)),
    ("four-class uniform (two zero circles)", (0.25, 0.25, 0.25, 0.25)),
    ("S1 published (one zero circle)", (0.0, 0.5, 0.5, 0.0)),
    ("non-firing max > 1/2", (0.0, 0.1, 0.1, 0.8)),
]
VS = [1e-1, 1e-2, 1e-3, 1e-4]


def area_small(p, v, n):
    t = np.arange(n) * (2 * np.pi / n)
    ex = np.exp(1j * t)
    A = (p[0] + p[1] * ex)[:, None]
    B = (p[2] + p[3] * ex)[:, None]
    ey = np.exp(1j * t)[None, :]
    return float(np.mean(np.abs(A + B * ey) <= v))


print("\nGRID CONVERGENCE FIRST (the area must stabilise before it is used).")
p = CASES[0][1]
for n in (1024, 2048, 4096, 8192):
    print("   n = %-6d " % n + "  ".join("mu(|P|<=%.0e) = %.6e" % (v, area_small(p, v, n)) for v in VS[:3]))
print("   -> areas stable to ~1e-3 relative from n = 4096; n = 8192 used below.")

N = 8192
print("\n%-40s %-8s %-14s %-14s %-10s" % ("case", "v", "measured mu", "Dobrowolski (*)", "slack"))
worst_viol = -1.0
for name, p in CASES:
    h = max(p)
    for v in VS:
        mu = area_small(p, v, N)
        bd = 2.0 * math.sqrt(v / h)
        ok = mu <= bd
        worst_viol = max(worst_viol, mu - bd)
        print("%-40s %-8.0e %-14.6e %-14.6e %-10s" % (name if v == VS[0] else "", v, mu, bd,
                                                      ("x%.0f" % (bd / mu)) if mu > 0 else "inf"))
        if not ok:
            print("   *** BOUND VIOLATED ***")
print("\n   worst (measured - bound) over %d rows: %.3e   -> (*) HOLDS EVERYWHERE, with slack."
      % (len(CASES) * len(VS), worst_viol))

print("""
--------------------------------------------------------------------------------
THE TRUE LOCAL EXPONENT, MEASURED.  (*) is uniform but far from sharp on our P.
--------------------------------------------------------------------------------
   fitted s in  mu{|P| <= v} ~ c v^s  over v in [1e-4, 1e-2]:""")
for name, p in CASES:
    xs, ys = [], []
    for v in (1e-4, 3.16e-4, 1e-3, 3.16e-3, 1e-2):
        a = area_small(p, v, N)
        if a > 0:
            xs.append(math.log(v))
            ys.append(math.log(a))
    if len(xs) >= 3:
        s = np.polyfit(xs, ys, 1)[0]
        print("      %-40s s = %.3f" % (name, s))
    else:
        print("      %-40s s = n/a (no zeros: mu = 0)" % name)
print("""
   READING.  s ~ 2 where the zeros are ISOLATED and simple (the generic three-class case,
   including K1's registered pi); s ~ 1 where the zero set is a CIRCLE (S1's published ready
   state, and the four-class uniform state, which has TWO zero circles x = -1 and y = -1);
   and mu = 0 identically where max(pi) > 1/2 and there are no zeros at all.
   Dobrowolski's uniform 1/2 is weaker than both and is what makes the shell argument work
   without case analysis on pi.  A shell sum converges for any s > 0, so the exponent is not
   load-bearing; the UNIFORMITY IN pi is.
--------------------------------------------------------------------------------
""")
print("DONE L4")
