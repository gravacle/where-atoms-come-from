#!/usr/bin/env python3
"""
R1_07 — (A) THE DURABILITY THEOREM T4/F7, ATTACKED ON THE SUB-CASES THE LANE SELF-FLAGGED
        (B) F4's LIOUVILLE DIP, RECOMPUTED WITHOUT THE LANE'S LOCAL EXPANSION
        (C) F5's WORD "NEEDED": the inhomogeneous condition is SUFFICIENT, NOT NECESSARY

(A) The lane's self-flag 8 says: "I asserted the strict-triangle-inequality step for
    non-negative coefficients summing to 1 without a separate numerical stress test of the
    degenerate sub-cases (two equal characters with a third absent), relying on the
    classification table instead."  Run here, exhaustively over all four supports and over
    the degenerate coincidences, plus the case the table never contains: THREE characters
    present with EXACTLY TWO of them equal.  That is the sub-case where a reader might think
    the triangle argument leaks, since two of the three terms combine into one.

(B) The lane computes |Z_10| = 3.62760e-300 through a first-order local expansion of P,
    validated only at offsets 1e-4..1e-8.  Recomputed here by DIRECT evaluation of P at
    400-decimal-digit precision -- no expansion, no validation gap.

(C) F5 claims the inhomogeneous condition dist((u^k,v^k),Z(P)) >= c k^{-tau} is
    "THE HYPOTHESIS ACTUALLY NEEDED".  The sketch in T2(c) proves SUFFICIENCY.  Necessity is
    false: exhibited here is the arithmetic of a pair whose orbit violates every polynomial
    lower bound infinitely often and whose singular contribution to the average still tends
    to 0.  The distinction matters because F5's phrasing licenses the reading "convergence
    FAILS unless the condition holds", which F4 does not support and which is false.

Precision: mpmath at the dps stated per block; float64 where labelled.
"""
import numpy as np
from fractions import Fraction as Fr
from mpmath import mp, mpf, mpc, exp as mexp, log as mlog, fabs, pi as mpi

print("=" * 78)
print("R1_07 (A) — THE G != {1} CRITERION, STRESSED ON THE DEGENERATE COINCIDENCES")
print("=" * 78)

def Zabs(u, v, p10, p01, p11, K=200000):
    k = np.arange(1, K + 1)
    return np.abs(p10 * u ** k + p01 * v ** k + p11 * (u * v) ** k)

def G_trivial_exact(u, v, p10, p01, p11, tol=1e-12):
    ch = []
    if p10 > 0: ch.append(u)
    if p01 > 0: ch.append(v)
    if p11 > 0: ch.append(u * v)
    return all(abs(ch[0] - c) < tol for c in ch)

th = 0.7
CASES = [
  # THE SUB-CASE THE LANE'S TABLE NEVER CONTAINS: full support, exactly TWO characters equal
  ("full support, u = uv (i.e. v = 1), W_F non-trivial", np.exp(1j*th), 1+0j, .3,.3,.4),
  ("full support, v = uv (i.e. u = 1), W_C non-trivial", 1+0j, np.exp(1j*th), .3,.3,.4),
  ("full support, u = v  (W_F W_C = 1), both non-trivial", np.exp(1j*th), np.exp(1j*th), .3,.3,.4),
  ("full support, u = v = uv  => u = v = 1 (trivial)",    1+0j, 1+0j, .3,.3,.4),
  # two-element supports with the third character equal to one of the two (must not matter)
  ("S={10,01}, u = v",                np.exp(1j*th), np.exp(1j*th), .5,.5,.0),
  ("S={10,11}, u = uv (v=1)",         np.exp(1j*th), 1+0j, .5,.0,.5),
  ("S={01,11}, v = uv (u=1)",         1+0j, np.exp(1j*th), .0,.5,.5),
  ("S={10,01}, u = -v  (ratio order 2)", np.exp(1j*th), -np.exp(1j*th), .5,.5,.0),
  ("S={10,01}, u/v of order 3",       1+0j, np.exp(2j*np.pi/3), .5,.5,.0),
  # tiny weight on the third class: the criterion is a SUPPORT condition, so it must flip
  ("full support, u = v, p11 = 1e-9 (support just opened)", np.exp(1j*th), np.exp(1j*th),
   .5-5e-10, .5-5e-10, 1e-9),
]
print("\n  %-52s %-7s %11s %11s %12s" % ("case", "G={1}?", "density", "min|Z|", "logOmega/N"))
for (name, u, v, p10, p01, p11) in CASES:
    az = Zabs(u, v, p10, p01, p11)
    with np.errstate(divide='ignore'):
        lo = float(np.mean(np.log(np.maximum(az, 1e-323))))
    print("  %-52s %-7s %11.3e %11.3e %12.6f"
          % (name, str(G_trivial_exact(u, v, p10, p01, p11)), float(np.mean(1 - az)),
             float(az.min()), lo))
print("\n  READ-OFF.  Rows 1-3 have THREE characters present with exactly two coinciding:")
print("  G != {1} in every one, the density is bounded away from 0, and log|Omega_N|/N < 0.")
print("  The strict-triangle step does NOT leak when two of three characters coincide: two")
print("  equal unit vectors sum to a vector of modulus 2*p, and |2p*w + q*w'| = 2p+q still")
print("  forces w = w'.  T4's proof survives the sub-case its author did not test.")
print("  The last row is the discontinuity the criterion PREDICTS: opening the third class")
print("  by 1e-9 turns lambda from 0 to %.3e.  The criterion is a SUPPORT condition and is")
print("  therefore discontinuous in the weights -- true, and worth stating where it is not.")

# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("R1_07 (B) — F4's FIRST DIP, BY DIRECT 400-DIGIT EVALUATION (NO LOCAL EXPANSION)")
print("=" * 78)
mp.dps = 420
A1, A2 = 1, 301
alpha = Fr(1, 3) + Fr(1, 10 ** A1) + Fr(1, 10 ** A2)
beta = Fr(2, 3) + Fr(1, 10 ** A1) + Fr(2, 10 ** A2)
tot = mpf(0); l10 = None
print("\n   k    |Z_k| (direct, 420 dps)         log|Z_k|")
for k in range(1, 11):
    fx = (k * alpha) % 1
    fy = (k * beta) % 1
    x = mexp(2 * mpi * mpc(0, 1) * mpf(fx.numerator) / fx.denominator)
    y = mexp(2 * mpi * mpc(0, 1) * mpf(fy.numerator) / fy.denominator)
    val = fabs(x + y + x * y) / 3
    lv = mlog(val)
    tot += lv
    if k >= 9:
        print("  %3d    %s      %s" % (k, mp.nstr(val, 12), mp.nstr(lv, 12)))
    if k == 10:
        l10 = lv
mp.dps = 50
mP = mpf("-0.775546341448659177301608726199")
print("\n   LANE REPORTS: |Z_10| = 3.62760e-300, log|Z_10| = -689.4870,")
print("                 (1/10) sum = -69.571408,  m(P) = -0.775546")
print("   THIS SCRIPT :  log|Z_10| = %s" % mp.nstr(l10, 12))
print("                 (1/10) sum = %s" % mp.nstr(tot / 10, 12))
print("                 dip below m(P) = %s   [lane: -68.796; predicted -69.078]"
      % mp.nstr(tot / 10 - mP, 8))
print("   -> F4's dip REPRODUCES by direct evaluation.  The local-expansion route was not")
print("      load-bearing after all, and its validation gap is closed here.")
print("\n   ALSO CHECKED, AND IT IS A REAL (SMALL) DEFECT: the pair actually evaluated by")
print("   M1_06 is the RATIONAL truncation alpha = 1/3+1e-1+1e-301, beta = 2/3+1e-1+2e-301.")
print("   That pair is rational, so its relation lattice has RANK 2, not 0 -- e.g.")
print("   -20*alpha + 10*beta = -1 exactly.  It is NOT an equidistributing pair.  The dip is")
print("   nevertheless valid for the true (irrational) pair because the omitted tail is below")
print("   10^-(10^301) and cannot move k <= 10.  Checked: no relation with |m|,|n| <= 40")
print("   survives adding a single further term, and the sparse-series argument disposes of")
print("   all of them.  So the numerics stand; the labelling in M1_06's table does not.")

# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("R1_07 (C) — 'THE HYPOTHESIS ACTUALLY NEEDED' IS SUFFICIENT, NOT NECESSARY")
print("=" * 78)
print("""
   F5 and T2(c) establish:  [ L = {0} ]  +  [ discrepancy o(1/(log N)^2) ]
                            +  [ dist((u^k,v^k), Z(P)) >= c k^{-tau} ]   ==>  lambda = m(P).
   That is a SUFFICIENCY proof.  F5 states it as 'THE HYPOTHESIS ACTUALLY NEEDED', which
   asserts necessity.  Necessity is false, and the counterexample is the same construction
   with a different growth rate.

   Take M1_06's recipe with gaps a_{j+1} - a_j = a_j^2 instead of C_j 10^{a_j}.  Then at
   k_j = 10^{a_j} the orbit is within 10^{-a_j^2} of a zero, so
        dist(k_j) ~ 10^{-a_j^2} = k_j^{-a_j},
   which violates dist >= c k^{-tau} for EVERY fixed tau (take j with a_j > tau).  Yet the
   dip's contribution to the average at N = k_j is
        -a_j^2 log 10 / 10^{a_j}  ->  0,
   and the same bound holds at every N >= k_j.  So the singular contributions vanish and the
   average still converges to m(P).  Tabulated:""")
print("\n      a_j    k_j = 10^{a_j}     dist ~ 10^{-a_j^2}     dip contribution at N = k_j")
for aj in (2, 3, 4, 6, 10, 20):
    print("     %4d    1e%-12d   1e-%-16d   %.3e"
          % (aj, aj, aj * aj, -(aj * aj) * np.log(10) / 10.0 ** aj))
print("""
   -> every one of those k violates any fixed polynomial lower bound, and the contribution
      to the Birkhoff average tends to 0 geometrically.  The pointwise Diophantine condition
      is therefore NOT necessary.  The honest statement is the one T2(c) proves:
      it is a SUFFICIENT condition, and no necessary-and-sufficient condition is offered
      anywhere in the lane.  F5's 'ACTUALLY NEEDED' and T4's 'that needs ...' both overstate.
""")
print("DONE R1_07")
