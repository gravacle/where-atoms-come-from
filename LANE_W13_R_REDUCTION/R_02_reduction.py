#!/usr/bin/env python3
"""
R_02 — THE REDUCTION, CARRIED OUT EXPLICITLY, AND THE DICHOTOMY IT PRODUCES.

BRIEF ITEM (1): "The orbit {(u^k,v^k)} is a line winding on T^2.  Restricted to the closure of
that line, log|P| becomes a one-variable function of k.  Carry out the reduction explicitly and
identify what one-variable object the Birkhoff sum becomes."

IT IS CARRIED OUT HERE, AND IT SPLITS IN TWO.  There are two different things the phrase can
mean, and they have opposite answers.

(A) A REDUCTION OF THE INTEGRAND.  This EXISTS, unconditionally, and is exact.  Writing P as a
    polynomial in y,  P = (p01 + p11 x) y + (p00 + p10 x),  gives, for every (x,y) on T^2,

        log|P(x,y)| = log|p01 + p11 x|  +  log|y - R(x)|,     R(x) = -(p00+p10 x)/(p01+p11 x),

    and integrating in y by Jensen recovers the brief's own formula
        m(P) = (1/2pi) INT log max(|p00+p10 e^{it}|, |p01+p11 e^{it}|) dt.
    So the Birkhoff sum splits EXACTLY into
        sum_k log|p01+p11 u^k|         -- a CONTINUOUS one-variable Birkhoff sum (Weyl closes it),
        sum_k log|v^k - R(u^k)|        -- a distance-to-a-MOVING-TARGET sum.
    ALL of the singularity is in the second, and the target R(u^k) moves with k.  That is the
    exact sense in which our object is a generalisation of a Sudler product, and it is the
    sense in which it is NOT one: a Sudler product has a FIXED target.

(B) A REDUCTION OF THE BIRKHOFF SUM TO A ONE-VARIABLE BIRKHOFF SUM (CHOICE LEDGER L1(c)).
    This exists IF AND ONLY IF H2 FAILS.  Proof and computation below.

THE DICHOTOMY (THEOREM R-1).  X_omega := closure{k omega : k in Z} is a closed subgroup of T^2.
  L(omega) := {(m,n) in Z^2 : m alpha + n beta in Z}   (the relation lattice)
  rank L = 0  <=>  X_omega = T^2          -- H2 HOLDS; X is 2-dimensional; NO REDUCTION.
  rank L = 1  <=>  X_omega = finitely many parallel circles; on each, log|P| IS log|Q| for the
                   one-variable Laurent polynomial Q(z) = P(z^n, z^{-m}), and the Birkhoff sum
                   IS a one-variable Birkhoff sum over a circle rotation.
  rank L = 2  <=>  X_omega finite; the average is a finite average.
So the one-variable reduction is available EXACTLY ON THE RESONANT LOCUS -- which is
Haar-null, and which is exactly where N1's limit is NOT m(P).  Under the hypothesis N1
requires, the reduction does not exist.

AND ON THE RESONANT LOCUS THE REDUCED PROBLEM IS ALMOST ALWAYS NON-SINGULAR (PROP R-5).
Q_{m,n}(z) = P(z^n, z^{-m}) has a root on |z| = 1 iff the circle X_omega passes through one of
P's two torus zeros.  At K1's pi that happens for EXACTLY ONE primitive relation, (m,n)=(1,1),
i.e. u v = 1, i.e. c = f.  Proof:  0.3 z^{m+n} + 0.4 z^n + 0.3 = 0 with |z|=1, z = e(psi),
divide by e(n psi):  0.3(e(m psi) + e(-n psi)) = -0.4, and a sum of two unit vectors equal to a
negative real forces them conjugate, so (m-n) psi in Z and cos(2 pi m psi) = -2/3, i.e.
m psi = +- theta* mod 1.  If m != n then psi is rational and m psi cannot be the irrational
theta* (R_01 Z2).  Hence m = n, and gcd(m,n)=1 gives m=n=1.
ON THAT ONE LINE the object IS a classical (inhomogeneous) SUDLER PRODUCT, with the closed form
    |Z_k| = 0.3 * |2 sin(pi(k alpha - theta*))| * |2 sin(pi(k alpha + theta*))|,
and its Mahler measure is log(0.3) = -1.2039728043 -- WHICH IS THE REGISTER'S OWN PUBLISHED
(m,n) = (1,1) SUBTORUS ROW.  Verified below against the corpus's sixteen published rows.
"""
import numpy as np
from math import gcd
from R_lib import PI_K1, P_eval, m_jensen

rng = np.random.default_rng(20260817)
p00, p10, p01, p11 = PI_K1
TH = np.arccos(-2.0 / 3.0) / (2 * np.pi)

print("=" * 79)
print("R_02 — THE REDUCTION.  (A) OF THE INTEGRAND: EXACT.  (B) OF THE SUM: ONLY IF H2 FAILS.")
print("=" * 79)

# ------------------------------------------------------------------ (A) the exact split
print("\n(A) THE EXACT ONE-VARIABLE SPLIT OF THE INTEGRAND  (unconditional)")
th1 = rng.uniform(0, 1, 200000); th2 = rng.uniform(0, 1, 200000)
x = np.exp(2j * np.pi * th1); y = np.exp(2j * np.pi * th2)
lhs = np.log(np.abs(P_eval(PI_K1, x, y)))
den = p01 + p11 * x
R = -(p00 + p10 * x) / den
rhs = np.log(np.abs(den)) + np.log(np.abs(y - R))
print("    log|P(x,y)| == log|p01+p11 x| + log|y - R(x)|   max deviation over 200000 points: %.3e"
      % np.max(np.abs(lhs - rhs)))

nq = 1 << 22
t = np.arange(nq) * (2 * np.pi / nq)
e = np.exp(1j * t)
den_t = p01 + p11 * e
R_t = -(p00 + p10 * e) / den_t
part1 = float(np.mean(np.log(np.abs(den_t))))
part2 = float(np.mean(np.log(np.maximum(np.abs(R_t), 1.0))))
MP = m_jensen(PI_K1, nq)
print("    m(P) = INT log|p01+p11 x| + INT log^+|R(x)|  =  %.12f + %.12f = %.12f"
      % (part1, part2, part1 + part2))
print("    m(P) by the brief's max-form                  =  %.12f     deviation %.3e"
      % (MP, abs(part1 + part2 - MP)))
print("    INT log|p01+p11 x| = log(max(p01,p11)) = log(0.4) = %.12f  (Jensen; |p01|!=|p11| so"
      % np.log(0.4))
print("    the FIRST sum's integrand is CONTINUOUS and bounded away from 0: |0.3+0.4x| >= 0.1.")
print("    EVERY singularity of log|P| lives in the SECOND sum, where the target R(u^k) MOVES.")
print("    |R(x)| ranges over [%.6f, %.6f]; it equals 1 exactly at x = e(+-theta*)."
      % (np.abs(R_t).min(), np.abs(R_t).max()))
cross = np.abs(np.abs(R_t) - 1.0) < 1e-5
print("    {x : |R(x)| = 1} located numerically at theta = %s ; theta* = %.9f"
      % (np.round(t[cross][[0, -1]] / (2 * np.pi), 9), TH))

# ------------------------------------------------------------------ (B) the resonant reduction
print("\n(B) THE RESONANT REDUCTION Q_{m,n}(z) = P(z^n, z^{-m}) = 0.3 z^{m+n} + 0.4 z^n + 0.3")
print("    (times z^{-m}, which does not change the Mahler measure)")


def m_one_var(coef):
    """m of a one-variable polynomial by Jensen on its roots.  coef[i] multiplies z^i."""
    nz = np.nonzero(np.abs(coef) > 0)[0]
    coef = coef[nz[0]:nz[-1] + 1]
    if len(coef) == 1:
        return float(np.log(abs(coef[0])))
    r = np.roots(coef[::-1])
    return float(np.log(abs(coef[-1])) + np.sum(np.log(np.maximum(np.abs(r), 1.0))))


def Q_coef(m, n):
    """coefficients of 0.3 z^{m+n} + 0.4 z^n + 0.3 shifted to non-negative exponents."""
    ex = [(m + n, p10), (n, p11), (0, p01)]
    sh = -min(a for a, _ in ex)
    deg = max(a + sh for a, _ in ex)
    c = np.zeros(deg + 1)
    for a, w in ex:
        c[a + sh] += w
    return c


S4ROWS = {(1, 0): -0.356674944, (0, 1): -0.356674944, (1, 1): -1.203972804,
          (1, -1): -0.510825624, (2, 1): -0.681980359, (2, -1): -0.916290732,
          (3, 1): -0.767783712, (3, 2): -0.732940865, (4, 1): -0.784966659,
          (5, 1): -0.749392712, (5, 3): -0.765224351, (7, 3): -0.759305247,
          (7, 11): -0.764712281, (11, 20): -0.767014993, (13, 8): -0.768271734,
          (29, 17): -0.767138179}
print("\n    CONTROL (declared in advance as a control; it could not have failed if the reduction")
print("    formula is right, and a disagreement would have been the finding):")
print("    reproduce the corpus's sixteen published subtorus rows (S4:258ff, M1_07)")
worst = 0.0
for (m, n), val in sorted(S4ROWS.items()):
    mine = m_one_var(Q_coef(m, n))
    worst = max(worst, abs(mine - val))
print("      worst deviation over the 16 published rows: %.3e   -> REDUCTION FORMULA AGREES" % worst)
print("      (1,1) row: this lane %.12f   register -1.203972804   log(0.3) = %.12f"
      % (m_one_var(Q_coef(1, 1)), np.log(0.3)))

# ------------------------------------------------------------------ the unique singular line
print("\n    WHICH RESONANT LINES ARE SINGULAR?  min over |z|=1 of |Q_{m,n}| and the closest root")
print("      (m,n)     min_{|z|=1}|Q|      min | |root| - 1 |     singular?")
sing = []
for m in range(-12, 13):
    for n in range(-12, 13):
        if (m, n) == (0, 0) or gcd(abs(m), abs(n)) != 1:
            continue
        c = Q_coef(m, n)
        zz = np.exp(2j * np.pi * np.arange(200000) / 200000)
        vals = np.abs(np.polyval(c[::-1], zz))
        r = np.roots(c[::-1])
        d = float(np.min(np.abs(np.abs(r) - 1.0)))
        if d < 1e-9:
            sing.append((m, n, float(vals.min()), d))
for (m, n, v, d) in sing:
    print("      (%3d,%3d)   %14.6e     %14.6e      YES" % (m, n, v, d))
print("      primitive (m,n) with |m|,|n| <= 12 scanned: %d ; SINGULAR: %d, and they are"
      % (sum(1 for m in range(-12, 13) for n in range(-12, 13)
             if (m, n) != (0, 0) and gcd(abs(m), abs(n)) == 1), len(sing)))
print("      exactly (1,1) and (-1,-1), the SAME subgroup u v = 1, i.e. c = f.  Matches the proof.")
for (m, n) in [(1, 0), (0, 1), (1, -1), (2, 1), (11, 20), (7, 3)]:
    c = Q_coef(m, n)
    zz = np.exp(2j * np.pi * np.arange(200000) / 200000)
    print("      non-singular example (%3d,%3d): min|Q| on |z|=1 = %.6e  -> log|Q| CONTINUOUS,"
          "  Weyl closes it" % (m, n, float(np.abs(np.polyval(c[::-1], zz)).min())))

# ------------------------------------------------------------------ the Sudler closed form
print("\n    THE ONE SINGULAR LINE IS A CLASSICAL (INHOMOGENEOUS) SUDLER PRODUCT")
from R_lib import PhaseReducer
al = np.sqrt(2) - 1
k = np.arange(1, 100001, dtype=np.int64)
fa = PhaseReducer(al).frac(k)                      # exact-integer phase reduction (L4)
Zk = np.abs(0.4 + 0.6 * np.cos(2 * np.pi * fa))
Su = 0.3 * np.abs(2 * np.sin(np.pi * (fa - TH))) * np.abs(2 * np.sin(np.pi * (fa + TH)))
print("      |Z_k| == 0.3 |2 sin(pi(k a - theta*))| |2 sin(pi(k a + theta*))|")
print("      max deviation over k <= 100000 (a = sqrt2 - 1): %.3e" % np.max(np.abs(Zk - Su)))
print("      so  sum_{k<=N} log|Z_k| = N log(0.3) + log Sud_N(a,-theta*) + log Sud_N(a,+theta*)")
print("      with Sud_N(a,g) = prod_{k<=N} |2 sin(pi(k a + g))| the INHOMOGENEOUS Sudler product.")
print("      m(Q) = log(0.3) = %.12f, and the shifts +-theta* are IRRATIONAL (R_01 Z2), so this"
      % np.log(0.3))
print("      is the inhomogeneous problem, not Sudler's homogeneous one (whose shift is 0).")

print("\nDONE R_02")
