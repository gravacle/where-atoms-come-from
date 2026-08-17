#!/usr/bin/env python3
"""
R_02 — THE ONE NUMBER IN LANE C's C-4 THAT SEPARATES TWO HYPOTHESES, COMPUTED IN EXACT
       RATIONAL ARITHMETIC INSTEAD OF QUADRATURE.

C-4 claims the RESONANT arm of B0b converges to the subtorus Mahler measure -0.810930204535,
NOT to m(P) = log(4/9) = -0.810930216216, a separation of 1.168e-08.  Lane C computes that
number two ways, BOTH float64: numpy.roots on a degree-31 polynomial, and a 2^22-node
trapezoid.  Under this lens both are suspect until an exact computation agrees.

THE EXACT COMPUTATION.  On B0b, |B| = |A|/2 POINTWISE (lane C's reversal identity, which is
true).  Write log|P| = log|A| + log|1 + tau y| with tau = B/A, |tau| = 1/2 everywhere.
On the rank-one orbit closure H = {(x,y) = (z^n, z^-m)} with gcd(m,n) = 1:

     tau y  =  psi(z^n) z^-m / 2 ,      psi(w) = (1 + 2w)/(2 + w) ,
and psi is a Blaschke factor: analytic on |w| <= 1 (pole at w = -2) with |psi| = 1 on |w| = 1,
so every Taylor coefficient of psi^k is bounded by 1 in modulus.

     INT_H log|A|  =  log(2/9) + INT log|2 + z^n| dz  =  log(2/9) + log 2  =  log(4/9) = m(P),
     INT_H log|1 + tau y|  =  Re SUM_{k>=1} (-1)^{k+1} <(tau y)^k> / k ,
     <(tau y)^k>  =  2^-k [z^{mk}] psi(z^n)^k  =  2^-k c^{(k)}_{mk/n}   nonzero only if n | k.

     THEOREM R2.   lambda_(m,n) - m(P)  =  SUM_{r>=1} (-1)^{nr+1} c^{(nr)}_{mr} / (n r 2^{nr}),
     where c^{(K)}_j = [w^j] psi(w)^K is an EXACT RATIONAL and |c| <= 1.

This is a convergent series of exact rationals whose r-th term is < 2^-nr, so two terms give
the value to 2^-2n.  NO QUADRATURE, NO ROOT-FINDER.  It also PROVES the Boyd-Lawton decay law
on B0b is GEOMETRIC, ~ 2^-n / n, not a power of n -- which is what lane C fitted (slope -7.7)
and correctly declined to defend, but for the wrong reason (it blamed the noise floor).
"""
from fractions import Fraction as Fr
import numpy as np
import mpmath as mp
import sys

mp.mp.dps = 60

def psi_pow_coeffs(K, J):
    """Taylor coefficients [w^0..w^J] of ((1+2w)/(2+w))^K, EXACT Fractions.
       ((1+2w)/(2+w))^K = 2^-K (1+2w)^K (1+w/2)^-K ."""
    # (1+2w)^K
    num = [Fr(0)] * (J + 1)
    ck = Fr(1)
    for i in range(0, min(K, J) + 1):
        # binom(K,i) 2^i
        num[i] = Fr(_binom(K, i)) * Fr(2) ** i
    # (1+w/2)^-K  = sum_j binom(-K, j) (w/2)^j = sum_j (-1)^j binom(K+j-1, j) w^j / 2^j
    den = [Fr(0)] * (J + 1)
    for j in range(J + 1):
        den[j] = Fr((-1) ** j * _binom(K + j - 1, j), 2 ** j)
    out = [Fr(0)] * (J + 1)
    for i in range(J + 1):
        if num[i] == 0:
            continue
        for j in range(J + 1 - i):
            out[i + j] += num[i] * den[j]
    inv = Fr(1, 2 ** K)
    return [inv * o for o in out]

_bc = {}
def _binom(n, k):
    if k < 0 or k > n:
        return 0
    key = (n, k)
    if key in _bc:
        return _bc[key]
    r = 1
    for i in range(k):
        r = r * (n - i) // (i + 1)
    _bc[key] = r
    return r

def gap_exact(m, n, rmax=3):
    """THEOREM R2: lambda_(m,n) - m(P) for B0b, as an exact Fraction sum plus a tail bound."""
    tot = Fr(0)
    for r in range(1, rmax + 1):
        K = n * r; j = m * r
        c = psi_pow_coeffs(K, j)[j]
        tot += Fr((-1) ** (K + 1)) * c / Fr(K * 2 ** K)
    # tail bound: |c| <= 1, so |tail| <= sum_{r>rmax} 1/(n r 2^{n r}) <= 2 / (n (rmax+1) 2^{n(rmax+1)})
    tail = Fr(2, n * (rmax + 1) * 2 ** (n * (rmax + 1)))
    return tot, tail

def m_subtorus_roots_mp(p, m, n, dps=40):
    """Independent high-precision check: mpmath.polyroots on p01 + p00 z^m + p11 z^n + p10 z^{m+n}."""
    old = mp.mp.dps; mp.mp.dps = dps
    deg = m + n
    coef = [mp.mpf(0)] * (deg + 1)
    coef[0] += mp.mpf(p[2].numerator) / p[2].denominator
    coef[m] += mp.mpf(p[0].numerator) / p[0].denominator
    coef[n] += mp.mpf(p[3].numerator) / p[3].denominator
    coef[deg] += mp.mpf(p[1].numerator) / p[1].denominator
    r = mp.polyroots(coef[::-1], maxsteps=200, extraprec=400)
    v = mp.log(coef[deg]) + sum(mp.log(abs(z)) for z in r if abs(z) > 1)
    mp.mp.dps = old
    return v

B0b = (Fr(4, 9), Fr(2, 9), Fr(1, 9), Fr(2, 9))

if __name__ == "__main__":
    mP = mp.log(mp.mpf(4) / 9)
    print("=" * 104)
    print("R_02 SECTION 1 — THE RESONANT SUBTORUS LIMIT AT THE CORPUS'S OWN RESONANCE (11,20),")
    print("                 BY THREE METHODS, ONE OF WHICH IS EXACT RATIONAL ARITHMETIC.")
    print("=" * 104)
    m_, n_ = 11, 20
    g, tail = gap_exact(m_, n_, rmax=3)
    gmp = mp.mpf(g.numerator) / g.denominator
    print(f"  EXACT (Theorem R2, 3 terms):  lambda_(11,20) - m(P) = {g}")
    print(f"                                                     = {mp.nstr(gmp, 20)}")
    print(f"     tail bound beyond r = 3:  |tail| <= {mp.nstr(mp.mpf(tail.numerator)/tail.denominator, 4)}")
    print(f"     leading coefficient c^(20)_11 = [w^11] psi(w)^20 = "
          f"{psi_pow_coeffs(20, 11)[11]}  = {mp.nstr(mp.mpf(psi_pow_coeffs(20,11)[11].numerator)/psi_pow_coeffs(20,11)[11].denominator, 12)}")
    lam_exact = mP + gmp
    print(f"     => lambda_(11,20) = m(P) + gap = {mp.nstr(lam_exact, 22)}")
    lam_mp = m_subtorus_roots_mp(B0b, m_, n_, dps=40)
    print(f"  mpmath polyroots dps=40:       lambda_(11,20) = {mp.nstr(lam_mp, 22)}")
    print(f"  lane C (numpy.roots, float64): lambda_(11,20) = -0.810930204535")
    print(f"  lane C (2^22 trapezoid):       lambda_(11,20) = -0.810930204535")
    print(f"\n  |exact - mpmath polyroots| = {mp.nstr(abs(lam_exact - lam_mp), 6)}")
    print(f"  |exact - lane C's printed 12 places| = "
          f"{mp.nstr(abs(lam_exact - mp.mpf('-0.810930204535')), 6)}   "
          f"(half-ulp of the 12th place = 5e-13)")
    print(f"  m(P) = log(4/9) = {mp.nstr(mP, 22)}")
    print(f"  SEPARATION lambda_(11,20) - m(P) = {mp.nstr(gmp, 10)}  "
          f"(lane C reports 1.168e-08)")

    print("\n" + "=" * 104)
    print("R_02 SECTION 2 — SO THE SEPARATION IS REAL.  BUT CAN LANE C's MEASUREMENT SEE IT?")
    print("=" * 104)
    avg = mp.mpf("-0.810930225")     # C_04's resonant B0b average at K = 1e7, as printed
    print(f"  C_04 prints, for the RESONANT arm on B0b at K = 1e7:  average = {avg}")
    print(f"     |average - subtorus limit| = {mp.nstr(abs(avg - lam_exact), 4)}")
    print(f"     |average - m(P)|           = {mp.nstr(abs(avg - mP), 4)}")
    print(f"     ratio = {mp.nstr(abs(avg-lam_exact)/abs(avg-mP), 4)}")
    print("""
  THE MEASUREMENT IS 2.3x CLOSER TO m(P) THAN TO THE LIMIT LANE C SAYS IT CONVERGES TO.
  At K = 1e7 the Birkhoff average's own finite-K wobble (5.4e-08 on the SAME carrier's
  DIOPHANTINE arm, where the target is known to be m(P)) is FOUR TIMES the separation being
  claimed (1.168e-08).  So C_04's resonant row DOES NOT MEASURE the distinction it is offered
  for.  The distinction is real -- Theorem R2 settles it in exact rationals -- but it is
  carried entirely by the closed-form subtorus computation, NOT by the K = 1e7 run.
  C-4's evidence line ("resonant -0.810930225 -> subtorus -0.810930204535") reads as though
  the run resolved it.  It did not, and on B4 the same row is worse: |avg - subtorus| =
  7.26e-09 against |avg - m(P)| = 7.75e-09, a 6% difference on a 5.6e-08 noise floor.""")

    print("\n" + "=" * 104)
    print("R_02 SECTION 3 — THE BOYD-LAWTON DECAY LAW ON B0b IS GEOMETRIC, NOT A POWER OF n.")
    print("            Theorem R2's exact gaps vs lane C's float64 ladder (C_06 LEG B).")
    print("=" * 104)
    LANE = {(1,1):-0.466321295399, (1,2):-0.916658991929, (2,3):-0.787515156935,
            (3,5):-0.808347305186, (5,8):-0.810871539969, (8,13):-0.810933046885,
            (13,21):-0.810930210766, (21,34):-0.810930216216, (34,55):-0.810930216216,
            (55,89):-0.810930216216, (89,144):-0.810930216216, (144,233):-0.810930216216,
            (233,377):-0.810930216215, (377,610):-0.810930216215}
    print(f"  {'(m,n)':>10s} {'EXACT gap (Thm R2)':>24s} {'2^-n / n':>12s} {'lane C |lam-m(P)|':>19s} "
          f"{'lane C gap SIGNED':>19s}  agree?")
    for (m_, n_) in [(1,1),(1,2),(2,3),(3,5),(5,8),(8,13),(13,21),(21,34),(34,55),(55,89),
                     (89,144),(144,233),(233,377),(377,610)]:
        rmax = 6 if n_ <= 4 else (3 if n_ <= 12 else 1)
        g, tail = gap_exact(m_, n_, rmax=rmax)
        gv = float(mp.mpf(g.numerator) / g.denominator)
        lane_signed = LANE[(m_, n_)] - float(mP)
        ok = abs(gv - lane_signed) < max(3e-13, 2e-3 * abs(gv))
        print(f"  {str((m_,n_)):>10s} {gv:24.3e} {2.0**(-n_)/n_:12.3e} {abs(lane_signed):19.3e} "
              f"{lane_signed:19.3e}  {ok}")
    print("""
  EXACT AND MEASURED AGREE IN SIGN AND MAGNITUDE ON EVERY ROW WHERE THE MEASUREMENT IS ABOVE
  ITS OWN NOISE FLOOR.  The exact column falls like 2^-n, so:

    * lane C's fitted "slope = -7.7 over 4-5 relations" for the zero-free rows is not a slow
      power law measured through a noise floor -- IT IS NOT A POWER LAW AT ALL.  Fitting
      log(err) vs log(n) on a geometric sequence returns a slope that GROWS with the window,
      so the number is a window artefact in exactly COR-E's defect class.  Lane C declines to
      defend it (D-4) for the right reason but names the wrong cause: the noise floor is not
      what makes the exponent meaningless; the functional form is.
    * B0b's Boyd-Lawton convergence needs NO root-finder and has NO noise floor.  At (377,610)
      the exact gap is ~2^-610 = 1e-184; lane C's numpy.roots returns 1.36e-12, which is 172
      orders of magnitude of pure float64 noise reported as a convergence measurement.""")
    sys.exit(0)
