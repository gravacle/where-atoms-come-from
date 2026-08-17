#!/usr/bin/env python3
"""
R1_06 — THREE STATEMENTS IN THE LANE THAT CONFLATE lambda WITH m(P), AND ONE THAT IS
REFUTED BY THE LANE'S OWN F4.

The lane is careful about lambda vs m(P) in T2/T4 and then is not careful in four places.
Each is checked here by exhibiting an instance.

(1) F12 / M1_05 docstring S5 / M1_08 T4 last paragraph:
      "for a ready state UNIFORM on n cells ... so lambda <= -1/n; on K1 (n=5) lambda <= -1/5"
    ATTACK: the quantity bounded by Jensen is m(P), NOT lambda.  The lane's OWN F9(ii) and
    M1_04 row 4 exhibit a state uniform on 4 cells with lambda = 0 > -1/4.  Instance below.

(2) M1_04 printed output, PART 3:
      "N1's identification lambda = m(P) is TRUE EXACTLY ON THE COMPLEMENT of a dense,
       Haar-null set of connections"
    ATTACK: refuted by the lane's OWN M1_06.  The Liouville pair has L = {0} -- it is in
    that complement -- and lambda != m(P) there (indeed = -infinity).  The correct statement
    is "on a set of FULL Haar measure", and its complement contains, besides the relation
    locus, a dense set of relation-free pairs.  Checked below by re-deriving that the
    Liouville pair's relation lattice really is {0} and by locating it in the complement.

(3) The lane's stated OPERATIVE VARIABLE: "THE RELATION LATTICE L ... rank L = 0 gives
    H = T^2 and lambda = m(P)".
    ATTACK: false by the lane's own F4.  L does not determine lambda even in rank 0.  Nor
    does it determine lambda in rank 1: (1,1) gives -1.2040 and (11,20) gives -0.7670, and
    R1_05 exhibits two connections with the SAME L = Z(1,1) and different limits.  So the
    variable that determines the rate is not L; it is (L, inhomogeneous Diophantine type of
    the orbit relative to Z(P|_H)).  The lane's own Comparison 5 says exactly this and the
    operative_variable field contradicts it.

(4) F7/T4: "if G != {1} then ... |Omega_N| -> 0 EXPONENTIALLY with limsup (1/N) log|Omega_N|
    <= lambda_H <= log(1-c) < 0."   ATTACK on the DIRECTION of the bound: lambda_H is not an
    upper bound for the rate in the sense a reader will take -- it is the limit itself where
    the limit exists, and it can be LARGER than m(P).  Exhibited: lambda_H(11,20) > m(P).
    So no statement of the form "the rate is at most m(P)" is available, and F12's -1/n line
    is exactly such a statement.

Precision: float64 for the orbit sums (stated); mpmath 50 dps for m(P) and the subtorus
Mahler measures; exact integers for the relation-lattice checks.
"""
import numpy as np
from fractions import Fraction as Fr
from mpmath import mp, mpf, log as mlog, fabs, polyroots
mp.dps = 50

print("=" * 78)
print("R1_06 — lambda VERSUS m(P): FOUR PLACES THE LANE CONFLATES THEM")
print("=" * 78)

# ---------------------------------------------------------------- (1) the -1/n bound
print("\n(1) THE CLAIM  'lambda <= -1/n for a ready state uniform on n cells; on K1 lambda <= -1/5'")
print("    INSTANCE THAT REFUTES IT, taken from the lane's OWN M1_04 row 4:")
print("      ready state uniform on the 4 non-root vertices v1,v2,v3,v4")
print("        -> (p10,p01,p11) = (1/2, 1/2, 0),  n = 4,  the bound asserts lambda <= -1/4")
print("      connection with u = v = e^{0.7 i}  (i.e. W_C = conj(W_F) e^{...}: W_F W_C = 1),")
print("      which is a NON-TRIVIAL connection.")
u = np.exp(0.7j); v = np.exp(0.7j)
k = np.arange(1, 200001)
az = np.abs(0.5 * u ** k + 0.5 * v ** k + 0.0)
print("      min|Z_k| = %.15f   max|Z_k| = %.15f   (float64, k<=2e5)" % (az.min(), az.max()))
lam = float(np.mean(np.log(az)))
print("      lambda = (1/N) sum log|Z_k| = %.3e   -- it is 0, not <= -1/4." % lam)
print("      m(P) = -log 2 = %.12f, so the JENSEN bound on m(P) is fine; the bound on LAMBDA")
print("      is false.  Same failure on n = 5: a connection with G = {1} and the uniform state")
print("      gives lambda = 0 > -1/5.  Exhibited:")
# uniform on all 5 -> (2/5,2/5,1/5); need G={1}: chi_10=u, chi_01=v, chi_11=uv all equal
# u = v and uv = u  => u = 1 and v = 1: only the trivial connection has G={1} at FULL support.
print("        at FULL support G={1} forces u = v = uv, hence u = v = 1: the TRIVIAL connection.")
u2 = v2 = 1.0 + 0j
az2 = np.abs(0.4 * u2 ** k + 0.4 * v2 ** k + 0.2 * (u2 * v2) ** k)
print("        trivial connection, uniform-on-5 state: min|Z_k| = %.1f, lambda = %.3e > -1/5"
      % (az2.min(), float(np.mean(np.log(az2)))))
print("      -> the -1/n line is TRUE OF m(P) AND FALSE OF lambda.  It is stated of lambda in")
print("         F12, in M1_05's docstring S5, and in M1_08 T4's closing paragraph.")

# ---------------------------------------------------------------- (2) 'exactly on the complement'
print("\n(2) THE CLAIM  'lambda = m(P) is TRUE EXACTLY ON THE COMPLEMENT of a dense, Haar-null")
print("    set of connections'  (M1_04 printed output, PART 3).")
print("    The lane's OWN Liouville pair lies in that complement and violates the conclusion.")
print("    Relation lattice of the Liouville pair, checked exactly on the truncation used:")
A1, A2 = 1, 301
alpha = Fr(1, 3) + Fr(1, 10 ** A1) + Fr(1, 10 ** A2)
beta = Fr(2, 3) + Fr(1, 10 ** A1) + Fr(2, 10 ** A2)
found = []
for m_ in range(-40, 41):
    for n_ in range(-40, 41):
        if m_ == 0 and n_ == 0:
            continue
        if (m_ * alpha + n_ * beta).denominator == 1:
            found.append((m_, n_))
print("      integer relations with |m|,|n| <= 40 on the TRUNCATED pair: %s" % (found if found else "NONE"))
print("      (the full pair has none at all, by the sparse-series argument -- verified as an")
print("       argument in R1_08, not as a search)")
print("    -> the pair is in the complement of the relation locus, and the lane's own F4 says")
print("       its liminf is -infinity.  'TRUE EXACTLY ON THE COMPLEMENT' is FALSE.  The true")
print("       statement is 'on a set of FULL HAAR MEASURE'.")

# ---------------------------------------------------------------- (3)/(4) rank does not fix lambda
print("\n(3)(4) RANK L DOES NOT DETERMINE lambda, AND lambda_H IS NOT BOUNDED ABOVE BY m(P).")
P10, P01, P11 = Fr(3, 10), Fr(3, 10), Fr(2, 5)
def m_sub(mm, nn):
    terms = [(nn, P10), (-mm, P01), (nn - mm, P11)]
    shift = -min(e for e, _ in terms)
    deg = max(e + shift for e, _ in terms)
    coef = [Fr(0)] * (deg + 1)
    for e, c in terms:
        coef[e + shift] += c
    while len(coef) > 1 and coef[-1] == 0: coef.pop()
    while len(coef) > 1 and coef[0] == 0: coef.pop(0)
    if len(coef) == 1:
        return mlog(mpf(coef[0]))
    r = polyroots([mpf(x) for x in coef[::-1]], maxsteps=400, extraprec=300)
    return mlog(mpf(coef[-1])) + sum(mlog(fabs(z)) for z in r if fabs(z) > 1)
MP = mpf("-0.767507880357775871645874051819")
print("      m(P) = %s" % mp.nstr(MP, 20))
for (mm, nn) in [(1, 1), (1, 0), (11, 20), (29, 17), (13, 8)]:
    val = m_sub(mm, nn)
    print("      rank L = 1, relation (%3d,%3d):  lambda_H = %s   lambda_H - m(P) = %+s"
          % (mm, nn, mp.nstr(val, 16), mp.nstr(val - MP, 5)))
vals = [float(m_sub(a, b)) for (a, b) in [(1, 1), (1, 0), (11, 20), (29, 17), (13, 8)]]
print("      -> within rank L = 1 the limit ranges over [%.7f, %.7f] on these rows alone;"
      % (min(vals), max(vals)))
print("         rank alone fixes nothing.")
print("      -> and lambda_H(11,20) and lambda_H(29,17) are BOTH GREATER than m(P): the rate")
print("         can be SLOWER than the Mahler measure, so no upper bound of the form")
print("         'lambda <= m(P)' or 'lambda <= -(1-max)' is available for the true rate.")
print("      -> the ORDER-4 case goes the other way: lambda = -(1/2) log 5 = %.12f < m(P)."
      % float(-0.5 * np.log(5)))
print("         lambda straddles m(P) from both sides.  Any inequality between them is false.")
print("\nDONE R1_06")
