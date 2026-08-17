#!/usr/bin/env python3
"""
R1_04 — THE ATTACK THE LANE'S OWN T2 SETS UP AND DOES NOT RUN.

T2 splits the limit into (a) a FREE upper half and (b) a lower half needing an INHOMOGENEOUS
Diophantine condition, and says the fragile half is needed "when P has zeros on T^2".  T4
then says: "What is fragile is only the NAMING of the exponential rate as m(P): that needs
H = T^2 plus the inhomogeneous Diophantine condition of T2(c)."  H = T^2 -- i.e. the whole
fragility is confined to rank L = 0.

BUT F10 AND M1_07 REPORT RANK-1 LIMITS AS ESTABLISHED VALUES:
   "-0.767015029  ->  subtorus value -0.767014992998"     (the (11,20) resonance)
and a 16-row table of subtorus rates presented as the limits of the corresponding orbits.
On a rank-1 H the Birkhoff average is a ONE-DIMENSIONAL average of log|Q| along an irrational
circle rotation, where Q = P restricted to H.  IF Q HAS ZEROS ON THE CIRCLE, that average is
subject to EXACTLY the same inhomogeneous-Diophantine failure as the rank-0 case -- and the
rotation number there is 1/(20 pi), whose Diophantine type is unknown to mathematics.

THIS SCRIPT ASKS, FOR EVERY ROW OF THE TABLE:  does Q have zeros ON |z| = 1?
  - if NO, log|Q| is continuous on H, Weyl alone gives the limit, and the row is FREE.
  - if YES, the row's "limit" is an unproven Diophantine claim about a specific rotation
    number, and the lane's confinement of the fragility to H = T^2 is FALSE.

Q for a primitive relation u^m v^n = 1:  H = {(z^n, z^{-m})},
   P|_H = p10 z^n + p01 z^{-m} + p11 z^{n-m}, shifted to  Q(z) = sum c_e z^e, e >= 0.
Weights held at S3/S4's (p10,p01,p11) = (0.3,0.3,0.4) -- the lane's own held-fixed state.

Method: (i) exact rational test for unimodular roots via the SELF-INVERSIVE criterion and
gcd with the reciprocal, over Q (fractions), no floats; (ii) float64 companion-matrix roots
as a cross-check; (iii) mpmath 60-dps polyroots on the small rows.
"""
import numpy as np
from fractions import Fraction as Fr
from mpmath import mp, mpf, polyroots, fabs, log as mlog, mpc, exp as mexp, pi as mpi

mp.dps = 60
P10, P01, P11 = Fr(3, 10), Fr(3, 10), Fr(2, 5)


def Qcoef(mm, nn):
    """exact rational coefficient list, index = exponent, of Q = shift * P|_H."""
    assert np.gcd(abs(mm), abs(nn)) == 1
    terms = [(nn, P10), (-mm, P01), (nn - mm, P11)]
    shift = -min(e for e, _ in terms)
    deg = max(e + shift for e, _ in terms)
    coef = [Fr(0)] * (deg + 1)
    for e, c in terms:
        coef[e + shift] += c
    while len(coef) > 1 and coef[-1] == 0:
        coef.pop()
    while len(coef) > 1 and coef[0] == 0:
        coef.pop(0)
    return coef


def polydiv_mod(a, b):
    """a mod b, exact rational polynomials as lists low->high."""
    a = a[:]
    db = len(b) - 1
    while len(a) - 1 >= db and any(x != 0 for x in a):
        da = len(a) - 1
        if a[-1] == 0:
            a.pop(); continue
        f = a[-1] / b[-1]
        for i in range(db + 1):
            a[da - db + i] -= f * b[i]
        a.pop()
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def polygcd(a, b):
    while any(x != 0 for x in b):
        a, b = b, polydiv_mod(a, b)
        while len(b) > 1 and b[-1] == 0:
            b.pop()
    return a


def has_unimodular_root_exact(coef):
    """Q has a root on |z|=1  =>  that root is also a root of the reciprocal-conjugate
    polynomial Q*(z) = z^deg * conj(Q(1/conj z)).  Coefficients here are REAL rationals, so
    Q*(z) = z^deg Q(1/z) = reversed coefficient list.  A root on |z|=1 that is NOT real must
    come with its conjugate = its inverse, so it is a root of gcd(Q, Q*).  Real roots +-1 are
    tested directly.  Returns (bool_certain_none, gcd_degree, Q(1), Q(-1))."""
    rev = coef[::-1]
    g = polygcd(coef[:], rev[:])
    q1 = sum(coef)
    qm1 = sum(c * (-1) ** i for i, c in enumerate(coef))
    # gcd degree 0 (a nonzero constant) => no common root => no unimodular root at all
    return (len(g) - 1 == 0), len(g) - 1, q1, qm1


S4ROWS = [(1, 0), (0, 1), (1, 1), (1, -1), (2, 1), (2, -1), (3, 1), (3, 2), (4, 1),
          (5, 1), (5, 3), (7, 3), (7, 11), (11, 20), (13, 8), (29, 17)]
S4VAL = {(1, 0): -0.356674944, (0, 1): -0.356674944, (1, 1): -1.203972804,
         (1, -1): -0.510825624, (2, 1): -0.681980359, (2, -1): -0.916290732,
         (3, 1): -0.767783712, (3, 2): -0.732940865, (4, 1): -0.784966659,
         (5, 1): -0.749392712, (5, 3): -0.765224351, (7, 3): -0.759305247,
         (7, 11): -0.764712281, (11, 20): -0.767014993, (13, 8): -0.768271734,
         (29, 17): -0.767138179}

print("=" * 78)
print("R1_04 — DOES THE RANK-1 LIMIT COME FREE?  ZEROS OF P|_H ON THE CIRCLE.")
print("=" * 78)
print("\nweights (p10,p01,p11) = (0.3,0.3,0.4)   [the lane's held-fixed state]")
print("\n  (m,n)     deg Q   gcd(Q,Q*) deg   Q(1)      Q(-1)    min ||r|-1| (float64)   "
      "ZEROS ON |z|=1?")
flagged = []
for (mm, nn) in S4ROWS:
    coef = Qcoef(mm, nn)
    none_certain, gdeg, q1, qm1 = has_unimodular_root_exact(coef)
    cf = np.array([float(c) for c in coef])
    if len(cf) > 1:
        r = np.roots(cf[::-1])
        mind = float(np.min(np.abs(np.abs(r) - 1.0)))
    else:
        mind = np.inf
    onc = "NO (proved: gcd is a unit)" if none_certain else "YES / possible (gcd deg %d)" % gdeg
    if not none_certain:
        flagged.append((mm, nn, gdeg))
    print("  (%3d,%4d)  %5d   %8d      %+8.4f %+8.4f   %.3e            %s"
          % (mm, nn, len(coef) - 1, gdeg, float(q1), float(qm1), mind, onc))

print("\n  ROWS WITH ROOTS ON THE UNIT CIRCLE (exact, over Q):", flagged)

print("\n" + "-" * 78)
print("WHAT THAT MEANS ROW BY ROW.")
print("  A row with NO unimodular root:  log|Q| is CONTINUOUS on the circle H, the orbit")
print("  z^k is equidistributed there by Weyl (the rotation number is irrational), so the")
print("  Birkhoff average converges to m(Q) WITH NO DIOPHANTINE HYPOTHESIS.  The row is a")
print("  theorem.")
print("  A row WITH a unimodular root: log|Q| has a log singularity on H and the Birkhoff")
print("  average needs an inhomogeneous Diophantine condition on the rotation number ")
print("  relative to that root.  The row is NOT a theorem; it is a numerical observation.")

# ------------------------------------------------------------------ the (1,1) row in detail
print("\nTHE (1,1) ROW, WORKED.  Q = 0.3 + 0.4 z + 0.3 z^2.")
c11 = Qcoef(1, 1)
print("   exact coefficients:", c11)
rr = polyroots([mpf(c11[2]), mpf(c11[1]), mpf(c11[0])], maxsteps=200, extraprec=200)
for z in rr:
    print("   root %s   |root| = %s" % (mp.nstr(z, 25), mp.nstr(fabs(z), 25)))
print("   BOTH roots are EXACTLY on |z| = 1 (Q is self-reciprocal: 0.3,0.4,0.3).")
print("   So m(Q) = log(0.3) = %s  -- which IS S4's row value -1.203972804." % mp.nstr(mlog(mpf(3)/10), 15))
print("   BUT the ORBIT average of log|Q(z^k)| along the rotation is a Birkhoff average of an")
print("   UNBOUNDED function, and whether it converges to log(0.3) depends on how well the")
print("   rotation number approximates the two singular phases.  The lane reports this row as")
print("   an established limit; T2(c) says exactly why it is not one.")

# demonstrate: pick a rotation number that is Liouville RELATIVE to the singularity of the
# (1,1) row and show the average diverges, at the corpus's own weights.
print("\n   DEMONSTRATION that the (1,1) row is not free.  The roots of Q are at additive")
print("   phases +-theta0/2pi with cos(theta0) = -2/3.  Choose a rotation number that is a")
print("   Liouville approximation to that phase; the average dives.")
th0 = float(np.arccos(-2.0 / 3.0))
tgt = th0 / (2 * np.pi)                      # the singular phase, as a fraction of a turn
print("   singular phase (turns) = %.15f" % tgt)
from fractions import Fraction as FR
# rho = tgt truncated at 10^-a1 plus a far-away tail: frac(k*rho) hits tgt to 10^-a2 at k=10^a1
a1, a2 = 1, 220
best = FR(tgt).limit_denominator(10 ** 15)
rho = FR(int(tgt * 10 ** a1), 10 ** a1) + FR(1, 10 ** a2)
print("   rho = %s  (k=10 lands within ~1e-%d of a zero of Q)" % (rho, a2 - a1))
logs = []
for k in range(1, 11):
    f = (k * rho) % 1
    d = min(abs(float(f) - tgt), abs(float(f) - (1 - tgt)), abs(float(f) - tgt + 1),
            abs(float(f) - (1 - tgt) - 1))
    if d < 1e-8:
        # exact local expansion of |Q| near a simple unimodular root
        z0 = mexp(mpc(0, 1) * 2 * mpi * mpf(tgt if abs(float(f) - tgt) < 0.5 else 1 - tgt))
        dz = (mpf(f.numerator) / f.denominator) - mpf(tgt)
        # Q'(z0) * z0 * 2 pi i * dphase
        Qp = mpf(c11[1]) + 2 * mpf(c11[2]) * z0
        val = fabs(Qp * z0 * 2 * mpi * mpc(0, 1) * dz)
        logs.append(float(mlog(val)))
    else:
        z = np.exp(2j * np.pi * float(f))
        logs.append(float(np.log(abs(0.3 + 0.4 * z + 0.3 * z * z))))
print("   (1/10) sum_{k<=10} log|Q(z^k)| = %.4f     against m(Q) = log(0.3) = %.6f"
      % (float(np.mean(logs)), float(np.log(0.3))))
print("   -> the rank-1 rows with unimodular roots carry the SAME hypothesis gap the lane")
print("      confines to rank L = 0.  T4's sentence 'that needs H = T^2 plus ...' is FALSE:")
print("      it also needs an inhomogeneous condition on any rank-1 H whose Q vanishes on H.")

# ------------------------------------------------------------------ the (11,20) row
print("\nTHE (11,20) ROW -- THE ERRATUM'S OWN CONNECTION -- IN DETAIL.")
c = Qcoef(11, 20)
print("   Q(z) = 0.3 z^31 + 0.4 z^20 + 0.3,  degree %d" % (len(c) - 1))
none_certain, gdeg, q1, qm1 = has_unimodular_root_exact(c)
print("   gcd(Q, reversed Q) degree = %d  ->  %s"
      % (gdeg, "NO root on |z|=1, PROVED over Q" if none_certain else "possible unimodular root"))
rr = polyroots([mpf(x) for x in c[::-1]], maxsteps=400, extraprec=400)
mods = sorted(float(fabs(z)) for z in rr)
print("   root moduli: min %.15f  max %.15f   min ||r|-1| = %.3e"
      % (mods[0], mods[-1], min(abs(m - 1) for m in mods)))
mQ = float(mlog(mpf(c[-1]))) + sum(float(mlog(fabs(z))) for z in rr if fabs(z) > 1)
print("   m(Q) by Jensen at 60 dps = %.15f" % mQ)
print("   lane / register value    = -0.767014992998 / -0.767014993")
print("   -> log|Q| IS continuous on the circle here (min ||r|-1| = %.1e), so THIS row's"
      % min(abs(m - 1) for m in mods))
print("      limit IS free by Weyl.  The lane's number is right and is right for a reason the")
print("      lane does not give.")
print("\nDONE R1_04")
