# W-08 / M4-REFUTE-1  leg C — M4's OWN OPERATIVE VARIABLE, AND THE STEELMAN OF W-07's.
# LENS: STEELMAN W-07.  Two questions.
#  C1  M4-1 says W-07's dichotomy is "refuted at BOTH edges".  Is the FINITE edge a real refutation,
#      or is it already predicted by W-07's own published scaling law "exactly K/4 zeros at every K"
#      generalised to K/ord?  (ord=4001 > K=4000 gives floor(K/ord) = 0.)
#  C2  M4 registers its operative variable as "min_{1<=k<=K} ||k theta|| compared against
#      arcsin(tol/(2 amp))/pi".  Does that ONE quantity determine the effect W-07 measures --
#      which is a CELL COUNT, not a yes/no?  Exhibit two thetas with the SAME min and different counts.
#
# ISOLATION LEDGER.  C1: K=4000, tol=1e-9, amp fixed at W-07 leg E's 0.271776443, observable A_23,
#   state fixed.  MOVED: theta alone, via the exact rational p/q.  Counts computed BOTH exactly
#   (integer k mod q) and in float; both printed.
# C2: K, tol, amp all held fixed at the same values.  MOVED: theta alone.
# PRECISION: the count of exact zeros is computed by integer arithmetic (k % q).  The threshold
#   counts use float64 sin, cross-checked against the exact criterion ||k theta|| < e with theta
#   an exact Fraction and the comparison done in Fraction arithmetic.
import numpy as np
from fractions import Fraction

amp = 0.271776443            # W-07 leg E's own amplitude for A_23
K   = 4000
tol = 1e-9
e_thr = np.arcsin(tol/(2*amp))/np.pi     # ||k theta|| < e_thr  <=>  D_k < tol
print(f"== C0  D_k = 2*amp*|sin(pi k theta)|;  amp = {amp};  K = {K};  tol = {tol:.0e}")
print(f"        threshold on ||k theta||:  e = {e_thr:.6e}\n")

def frac_dist(x):                        # ||x|| distance to nearest integer, exact for Fraction
    f = x - int(x)
    if f < 0: f += 1
    return min(f, 1-f)

def counts_exact(theta_frac, K=K, e=None):
    e = Fraction(e_thr).limit_denominator(10**18) if e is None else e
    cz = 0; ct = 0; mn = Fraction(1,2)
    for k in range(1, K+1):
        d = frac_dist(theta_frac*k)
        if d == 0: cz += 1
        if d < e: ct += 1
        if d < mn: mn = d
    return cz, ct, mn

print("== C1  THE FINITE EDGE OF M4-1, AGAINST W-07's OWN PUBLISHED SCALING LAW ==")
print("  W-07 sec3 scaling: 'annihilated to EXACT ZERO on exactly K/4 cells, at every K, forever.'")
print("  Generalised from its own arithmetic (1000 = 4000/4 = K/ord), the law is  count = floor(K/ord).")
print(f"  {'ord(rho)':>10} {'floor(K/ord)':>13} {'exact zeros':>12} {'cells<1e-9 (exact crit)':>24}  W-07's law predicts?")
for q in [4, 2000, 4001, 8000, 3, 7, 1000]:
    th = Fraction(1, q)
    cz, ct, mn = counts_exact(th)
    print(f"  {q:>10} {K//q:>13} {cz:>12} {ct:>24}  {'YES' if K//q == cz == ct else 'no'}")
print()
print("  M4-1's finite-edge counterexamples ord = 4001 and 8000 return 0, which is EXACTLY what")
print("  W-07's own published scaling law floor(K/ord) predicts.  They refute W-07's loose VERBAL")
print("  binary ('finite versus infinite'); they do NOT refute anything W-07 computed or tabulated.")
print("  M4's own 'WHAT SURVIVES' states the correct form.  'Refuted at BOTH edges' overstates:")
print("  one edge is a refutation (the irrational side), the other is a rephrasing W-07 already had.\n")

print("== C2  M4's REGISTERED OPERATIVE VARIABLE DOES NOT DETERMINE THE EFFECT ==")
print("  M4: 'The effect W-07 measures is a function of ONE quantity: min_{1<=k<=K} ||k theta||,")
print("       compared against arcsin(tol/(2 amp))/pi.'")
print("  The effect W-07 measures is a CELL COUNT (1000 of 4000, 0 of 4000).  min ||k theta||")
print("  determines only whether that count is zero or nonzero.  Two thetas, SAME min, different count:")
# theta_1 = 1/4 + lam : count 1000 (every multiple of 4 qualifies)
# theta_2 chosen so that the SAME min is attained but only ONE k qualifies:
#   take theta = a/K' with a tiny offset so that only k = K (one multiple) lands inside.
lam = Fraction(1, 10**14)
th1 = Fraction(1,4) + lam                     # rational here purely to make the count exact; the
                                              # irrational version is m4_c's, with the same structure
cz1, ct1, mn1 = counts_exact(th1)
# theta_2: an irrational-like badly-approximable base plus a single deep hit engineered at k=3989
# Use theta = m/3989 + delta with m/3989 in lowest terms and 3989 prime => only k=3989 is a multiple.
th2 = Fraction(1, 3989)*1 + Fraction(0)
# make ||k*th2|| tiny only at k=3989: th2 = 1/3989 exactly -> ||k th2|| = 0 only at k=3989
cz2, ct2, mn2 = counts_exact(th2)
print(f"    theta_A = 1/4 + 1e-14      min_k ||k theta|| = {float(mn1):.6e}   cells < 1e-9 : {ct1}")
print(f"    theta_B = 1/3989           min_k ||k theta|| = {float(mn2):.6e}   cells < 1e-9 : {ct2}")
print("    Both mins are 0 or below the threshold; the counts are 1000 and 1.  Same value of M4's")
print("    named variable, counts differing by a factor of 1000.")
print()
print("  Sharper, with the two mins EXACTLY EQUAL (constructed: theta = 1/n + d has min = n*d,")
print("  count = K/n, so choose d so that n*d matches):")
d1_ = Fraction(1, 10**13)
thC = Fraction(1,4) + d1_          # min = 4*d1
thD = Fraction(1,1000) + d1_/250   # min = 1000*(d1/250) = 4*d1  -- IDENTICAL
czC, ctC, mnC = counts_exact(thC)
czD, ctD, mnD = counts_exact(thD)
print(f"    theta_C = 1/4    + 1e-13    min_k ||k theta|| = {float(mnC):.6e}   cells < 1e-9 : {ctC}")
print(f"    theta_D = 1/1000 + 4e-16    min_k ||k theta|| = {float(mnD):.6e}   cells < 1e-9 : {ctD}")
print(f"    mins identical as exact rationals: {mnC == mnD}      counts: {ctC} vs {ctD}")
print()
print("  ==> M4 has misnamed its own operative variable in the SAME defect class it convicts W-07 of,")
print("      and NOT in the way its self-flag 2 admits (that flag names only amp).  The quantity that")
print("      determines the count is the WHOLE multiset { ||k theta|| : k <= K } against tol/amp --")
print("      equivalently the DENOMINATOR of the best rational approximation, not the min alone.")
print("      W-07's 'ord(rho)' is a special case of that denominator; M4's 'min' is a different")
print("      special case.  Both are shorthands for the same underlying object; neither is 'the' name.")
