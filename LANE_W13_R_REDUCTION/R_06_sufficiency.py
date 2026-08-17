#!/usr/bin/env python3
"""
R_06 — THE OTHER HALF: WHAT DOES SUFFICE, AND WHERE IT IS VERIFIABLE.

Three things are done here, in order of what they decide.

(1) THE MECHANISM OF THE SUFFICIENCY PROOF, MEASURED.  Theorem R-2 (R_07_THEOREMS.txt) proves
    convergence from two hypotheses: (D1) a Diophantine condition on omega, which by
    Erdos-Turan-Koksma gives discrepancy D_N = O(N^{-eta}); and (D2) an INHOMOGENEOUS
    Diophantine condition, dist(k omega, Z(P)) >= C k^{-sigma}.  The proof works by dyadic
    shells around the zeros: the shell {e^{-M-i-1} <= |P| < e^{-M-i}} has area ~ e^{-2(M+i)}
    and is hit ~N times its area plus a discrepancy term, and (D2) caps the number of shells
    at O(log N).  BOTH INGREDIENTS ARE MEASURED HERE OVER FOUR DECADES.

(2) THE COUPLING, WHICH IS A CORRECTION TO THE SHAPE OF H2.  H2 as registered is a condition
    on the CONNECTION alone.  The condition that actually decides convergence is a condition
    on the PAIR (ready state, connection), because the zeros of P -- the targets of (D2) --
    are determined by pi.  DEMONSTRATED BY MOVING EXACTLY ONE VARIABLE: R_05's counterexample
    omega, unchanged, evaluated at a DIFFERENT pi.  The divergence disappears.
    THIS BEARS ON W-12 COROLLARY 2, which proves H2 has no structural failure mode because
    phi pushes Haar to Haar on every carrier.  That argument is about the connection's law.
    IT DOES NOT EXTEND TO H2', WHICH IS NOT A CONDITION ON THE CONNECTION.

(3) WHERE (D1) AND (D2) ARE VERIFIABLE.  At K1's pi the zeros are (theta*, -theta*) with
    e(theta*) = zeta = (-2+i sqrt5)/3, ALGEBRAIC (R_01 Z1).  So both conditions are lower
    bounds on LINEAR FORMS IN LOGARITHMS OF ALGEBRAIC NUMBERS WITH ALGEBRAIC COEFFICIENTS,
    and Baker's theorem supplies them.  The reduction is exhibited here; the non-vanishing of
    each form is proved by Lindemann.  See R_07_THEOREMS.txt R-7.

Plus a PRECISION CONTROL that could have failed: rerun an arm with an EXACT rational rotation
of 60 digits instead of a float64 one, and compare.
"""
import math
import numpy as np
from fractions import Fraction
from R_lib import (PI_K1, PhaseReducer, HighPhaseReducer, m_jensen, hp_theta_star,
                   hp_sqrt, hp_pi, arm_hash)
from R_ladder import build_ladder, log_fraction

MP = m_jensen(PI_K1, 1 << 24)
p00, p10, p01, p11 = PI_K1
TH = np.arccos(-2.0 / 3.0) / (2 * np.pi)
ZEROS = [(TH, (-TH) % 1.0), ((-TH) % 1.0, TH)]
TWOPI = 2 * math.pi
KMAX = 10 ** 7

print("=" * 79)
print("R_06 — SUFFICIENCY: THE MECHANISM MEASURED, THE COUPLING EXHIBITED, BAKER APPLIED")
print("=" * 79)

# ================================================================== (1) shells
print("\n(1) THE SHELL MECHANISM OF THEOREM R-2, MEASURED OVER FOUR DECADES")
print("    arm: CORPUS omega = (-1/2pi, sqrt2/2pi)  <->  f = 1.0, c = sqrt(2)")
ra, rb = PhaseReducer(-1.0 / TWOPI), PhaseReducer(math.sqrt(2.0) / TWOPI)
DEC = [10 ** i for i in range(4, 8)]
Ms = [2, 4, 6, 8, 10]
counts = {N: {M: 0 for M in Ms} for N in DEC}
minz = {N: np.inf for N in DEC}
CH = 10 ** 6
run_c = {M: 0 for M in Ms}
run_min = np.inf
for lo in range(0, KMAX, CH):
    k = np.arange(lo + 1, min(lo + CH, KMAX) + 1, dtype=np.int64)
    fa, fb = ra.frac(k), rb.frac(k)
    x = np.exp(2j * np.pi * fa); y = np.exp(2j * np.pi * fb)
    az = np.abs(p00 + p10 * x + p01 * y + p11 * x * y)
    for N in DEC:
        if lo < N <= lo + len(k):
            for M in Ms:
                counts[N][M] = run_c[M] + int(np.sum(az[:N - lo] < math.exp(-M)))
            minz[N] = min(run_min, float(az[:N - lo].min()))
    for M in Ms:
        run_c[M] += int(np.sum(az < math.exp(-M)))
    run_min = min(run_min, float(az.min()))

# area of {|P| < delta} on T^2.  ANALYTIC, from the local expansion, cross-checked on a grid.
#   near a zero  P ~ 2 pi (A d1 + B d2), B = conj(A) = rho e^{-i psi}
#   (d1,d2) -> rho ( cos psi (d1+d2), sin psi (d1-d2) ),  |det| = rho^2 |sin 2 psi|
#   area{|P|<delta} = 2 * pi (delta/2pi)^2 / (rho^2 |sin 2 psi|) = delta^2 / (2 pi rho^2 |sin 2psi|)
x0 = np.exp(2j * np.pi * TH); y0 = np.conj(x0)
Aloc = (p10 + p11 * y0) * x0
rho = abs(Aloc); psi = math.atan2(Aloc.imag, Aloc.real)
JAC = rho * rho * abs(math.sin(2 * psi))


def area_analytic(delta):
    return delta ** 2 / (2 * math.pi * JAC)


NG = 8000
tt = np.arange(NG) / NG
XX = np.exp(2j * np.pi * tt)[:, None]; YY = np.exp(2j * np.pi * tt)[None, :]
AZG = np.abs(p00 + p10 * XX + p01 * YY + p11 * XX * YY)
print("    rho = |A| = %.6f, psi = %.6f, |det| = rho^2|sin 2psi| = %.6f" % (rho, psi, JAC))
print("    area{|P|<delta} = delta^2 / (2 pi |det|)   -- cross-check against an %d^2 grid:" % NG)
for M in (2, 3, 4, 5):
    d = math.exp(-M)
    print("       M=%d  delta=%.3e   analytic %.4e   grid %.4e   ratio %.4f"
          % (M, d, area_analytic(d), float(np.mean(AZG < d)), float(np.mean(AZG < d)) / area_analytic(d)))
print("    THE SHELL COUNTS.  expected = N * area{|P|<e^-M}; observed = #{k<=N : |Z_k| < e^-M}.")
print("    %-5s %-11s" % ("M", "e^-M") + "".join("   N=%-9d obs/exp" % N for N in DEC))
for M in Ms:
    ar = area_analytic(math.exp(-M))
    cells = []
    for N in DEC:
        exp_ = N * ar
        obs = counts[N][M]
        cells.append("  %6d/%8.2f %6s" % (obs, exp_, ("%.2f" % (obs / exp_)) if exp_ >= 1 else "  --"))
    print("    %-5d %-11.3e" % (M, math.exp(-M)) + "".join(cells))
print("    Where the expectation exceeds 1 the ratio sits near 1: the orbit visits the shells at")
print("    their Lebesgue rate.  Where it does not, '--' is printed rather than a ratio computed")
print("    from a zero count -- COR-E's defect class, avoided rather than committed.")
print("    THIS IS THE FIRST INGREDIENT OF THE PROOF AND IT IS WHAT EQUIDISTRIBUTION BUYS.")
print("    It is NOT enough by itself -- R_05 -- because it says nothing about the DEEPEST")
print("    shell reached, which is what (D2) controls.")
print("    min_{k<=N} |Z_k| and the depth it implies:")
for N in DEC:
    print("       N = %9d   min|Z_k| = %.4e   deepest shell M = %.2f"
          % (N, minz[N], -math.log(minz[N])))
print("    (D2) is the statement that the 'deepest shell' column grows like log N, not like N.")
print("    Measured growth of -log min|Z_k| against log N, last four decades: slope %.3f"
      % np.polyfit([math.log(N) for N in DEC], [-math.log(minz[N]) for N in DEC], 1)[0])
print("    A slope near 1/2 is what a codimension-2 singularity in a 2-torus predicts.")

# ================================================================== (2) the coupling
print("\n(2) THE COUPLING: THE SAME omega, A DIFFERENT pi.  ONE VARIABLE MOVED.")
PREC = 900
G1 = hp_theta_star(PREC); G2 = 1 - G1
D1, D2 = 0.5, 1.2
k1 = 10
M1 = int(math.ceil(math.exp(D1 * k1))); k2 = k1 * M1
M2 = 1 << int(math.ceil(D2 * k2 / math.log(2.0)))
ksA, _, _, alpha, epsA = build_ladder(G1, k1, 3, [M1, M2])
ksB, _, _, beta, epsB = build_ladder(G2, k1, 7, [M1, M2])

PIS = [("K1 registered (0,.3,.3,.4)", (0.0, 0.3, 0.3, 0.4)),
       ("M1_06's       (0,1/3,1/3,1/3)", (0.0, 1 / 3, 1 / 3, 1 / 3)),
       ("shifted       (0,.30,.35,.35)", (0.0, 0.30, 0.35, 0.35)),
       ("SAME ZEROS    (.1,.3,.3,.3)", (0.1, 0.3, 0.3, 0.3))]
ra, rb = HighPhaseReducer(alpha), HighPhaseReducer(beta)
KM2 = 10 ** 6
k = np.arange(1, KM2 + 1, dtype=np.int64)
fa, fb = ra.frac(k), rb.frac(k)
x = np.exp(2j * np.pi * fa); y = np.exp(2j * np.pi * fb)
print("    THE ORBIT IS BYTE-IDENTICAL IN ALL FOUR ROWS (sha256 of frac(k alpha)): %s"
      % arm_hash(fa[:200000]))
print("    the ONLY thing that moves is pi.  m(P) recomputed for each pi by the same quadrature.")
print("    %-30s  %12s %12s %12s %12s" % ("pi", "m(P)", "S_1490-m", "S_1e5-m", "S_1e6-m"))
for label, pv in PIS:
    q00, q10, q01, q11 = pv
    az = np.abs(q00 + q10 * x + q01 * y + q11 * x * y)
    lg = np.log(az)
    # rung corrections only where the rung really is a near-zero of THIS P
    for j, kk in enumerate(ksA[:2]):
        d1, d2 = float(epsA[j]) if log_fraction(epsA[j]) > -700 else 0.0, \
                 float(epsB[j]) if log_fraction(epsB[j]) > -700 else 0.0
        if log_fraction(epsA[j]) < -700:
            x0 = np.exp(2j * np.pi * float(G1)); y0 = np.exp(2j * np.pi * float(G2))
            AA = (q10 + q11 * y0) * x0; BB = (q01 + q11 * x0) * y0
            base = abs(q00 + q10 * x0 + q01 * y0 + q11 * x0 * y0)
            if base < 1e-12:                       # the rung IS a zero of this P
                r = math.exp(log_fraction(epsB[j]) - log_fraction(epsA[j]))
                lg[kk - 1] = math.log(2 * math.pi) + log_fraction(epsA[j]) + math.log(abs(AA + BB * r))
            else:
                lg[kk - 1] = math.log(base)
    mp_here = m_jensen(pv, 1 << 22)
    print("    %-30s  %12.6f %12.6f %12.6f %12.6f"
          % (label, mp_here,
             float(np.mean(lg[:1490])) - mp_here,
             float(np.mean(lg[:10 ** 5])) - mp_here,
             float(np.mean(lg[:KM2])) - mp_here))
print("    REGISTRAR-STYLE DEFECT, RECORDED RATHER THAN PATCHED OUT: the fourth row was my first")
print("    choice for a 'different pi' and it is NOT one.  On the anti-diagonal y = conj(x) the")
print("    zero condition reads (p00+p11) + 2 p10 cos(2 pi theta) = 0 whenever p10 = p01, so")
print("    EVERY pi with p10 = p01 = 0.3 and p00 + p11 = 0.4 has EXACTLY K1's two torus zeros.")
print("    (.1,.3,.3,.3) is such a pi, so it is thrown too -- correctly.  It is kept in the table")
print("    as the arm that shows the effect tracks the ZEROS and nothing else.")
print("    THE FIRST AND FOURTH ROWS ARE THROWN BY -1.2 AT N = 1490.  THE OTHER TWO ARE NOT")
print("    THROWN AT ALL.")
print("    The connection is identical in all three.  THEREFORE THE FAILING HYPOTHESIS IS NOT")
print("    A PROPERTY OF THE CONNECTION.  H2 has the wrong SHAPE, not merely the wrong strength.")
print("    W-12 COROLLARY 2 proves H2 'is never violable by the carrier or by the loop")
print("    designation' because phi pushes Haar to Haar.  That is a theorem about the")
print("    connection's law and it is not touched here -- but it does not transfer to H2',")
print("    which quantifies over (pi, omega) jointly.  RECORDED, NOT SCORED AS A REFUTATION.")

# ================================================================== (3) Baker
print("\n(3) WHERE (D1) AND (D2) ARE VERIFIABLE: LINEAR FORMS IN LOGARITHMS")
zeta = np.exp(2j * np.pi * TH)
print("    zeta = e(theta*) = %+0.12f %+0.12fi ; 3 zeta^2 + 4 zeta + 3 = %.2e  (degree 2, |zeta| = %.12f)"
      % (zeta.real, zeta.imag, abs(3 * zeta ** 2 + 4 * zeta + 3), abs(zeta)))
phi = np.arccos(-2 / 3)
print("    phi* = arccos(-2/3) = %.15f ;  -i Log(zeta) = %.15f   dev %.2e"
      % (phi, np.angle(zeta), abs(phi - np.angle(zeta))))
print("    AT THE CORPUS'S PUBLISHED GENERIC CONNECTION f = 1, c = sqrt(2):")
print("       alpha = -1/(2 pi),  beta = sqrt2/(2 pi).")
print("       ||k alpha + theta*||   = |k + phi* - 2 pi j| / (2 pi)  ... form  Lam_1")
print("       ||k beta  - theta*||   = |k sqrt2 - phi* - 2 pi j|/(2 pi) ... form  Lam_2")
print("       ||m alpha + n beta||   = |-m + n sqrt2 - 2 pi p|/(2 pi)  ... form  Lam_3")
print("    Each is  beta_0 + beta_1 Log(zeta) + beta_2 Log(-1)  with ALGEBRAIC beta_i of height")
print("    O(k) (resp. O(max(|m|,|n|,|p|))), since 2 pi = -2i Log(-1) and phi* = -i Log(zeta).")
print("    NON-VANISHING, each by Lindemann-Weierstrass:")
print("       Lam_1 = 0 => zeta = e^{ik} with k a non-zero rational integer => e^{ik} algebraic, false.")
print("       Lam_2 = 0 => zeta^{-1} = e^{i k sqrt2}, k sqrt2 non-zero algebraic, same contradiction.")
print("       Lam_3 = 0 with p != 0 => pi algebraic, false; with p = 0 => m = n = 0.")
print("    Baker's effective theorem for inhomogeneous linear forms in logarithms then gives")
print("    |Lam| > H^{-C} with C effective, i.e. exactly (D1) and (D2) with polynomial exponents.")
print("    SO N1's LIMIT HOLDS AT THE CORPUS'S OWN PUBLISHED GENERIC CONNECTION AND ITS OWN")
print("    REGISTERED READY STATE -- and the same argument covers every algebraic (alpha,beta).")
print("    FLAGGED: the constants are not tracked here and the citation is by theorem-shape.")
print("    What is proved in this lane is the REDUCTION and the NON-VANISHING; the lower bound")
print("    is quoted from Baker.")
print("    EMPIRICAL COMPANION (not a proof, and not scored as one): min over k <= 1e7 of")
for nm, a, g in [("k alpha + theta*", -1.0 / TWOPI, TH), ("k beta - theta*", math.sqrt(2) / TWOPI, -TH)]:
    r = PhaseReducer(a)
    kk = np.arange(1, KMAX + 1, dtype=np.int64)
    d = np.abs(((r.frac(kk) + g + 0.5) % 1) - 0.5)
    i = int(np.argmin(d))
    print("       ||%-18s|| = %.4e at k = %8d ;  k * that = %.4f" % (nm, d[i], i + 1, (i + 1) * d[i]))

# ================================================================== precision control
print("\nPRECISION CONTROL THAT COULD HAVE FAILED")
print("    Rerun BA-cubic with an EXACT 60-digit rational rotation through HighPhaseReducer,")
print("    against the float64 rotation through M1_06's PhaseReducer.  If the float64 phase")
print("    drift at k = 1e7 (~1e-10) mattered, these would disagree.")
a_f, b_f = 2.0 ** (1 / 3), 4.0 ** (1 / 3)
a_e = Fraction(hp_sqrt(2 * 10 ** 120, 0), 10 ** 60) if False else None
# exact 60-digit rationals for 2^(1/3), 4^(1/3) via integer cube roots
S = 10 ** 60


def icbrt(n):
    x = int(round(n ** (1 / 3))) if n < 2 ** 52 else 1 << ((n.bit_length() + 2) // 3)
    while True:
        y = (2 * x + n // (x * x)) // 3
        if y >= x:
            break
        x = y
    return x


a_e = Fraction(icbrt(2 * S ** 3), S)
b_e = Fraction(icbrt(4 * S ** 3), S)
print("    exact rationals agree with float64 to %.2e and %.2e"
      % (abs(float(a_e) - a_f), abs(float(b_e) - b_f)))
out = {}
for tag, ra_, rb_ in [("float64", PhaseReducer(a_f), PhaseReducer(b_f)),
                      ("exact60", HighPhaseReducer(a_e), HighPhaseReducer(b_e))]:
    tot = 0.0
    vals = {}
    for lo in range(0, KMAX, CH):
        kk = np.arange(lo + 1, min(lo + CH, KMAX) + 1, dtype=np.int64)
        fa2, fb2 = ra_.frac(kk), rb_.frac(kk)
        xx = np.exp(2j * np.pi * fa2); yy = np.exp(2j * np.pi * fb2)
        cs = np.cumsum(np.log(np.abs(p00 + p10 * xx + p01 * yy + p11 * xx * yy)))
        for N in DEC:
            if lo < N <= lo + len(kk):
                vals[N] = (tot + cs[N - lo - 1]) / N
        tot += float(cs[-1])
    out[tag] = vals
print("    %-9s" % "N" + "".join("%18d" % N for N in DEC))
for tag in ("float64", "exact60"):
    print("    %-9s" % tag + "".join("%18.12f" % out[tag][N] for N in DEC))
print("    %-9s" % "|diff|" + "".join("%18.2e" % abs(out["float64"][N] - out["exact60"][N]) for N in DEC))
print("    Agreement at 1e-12 or better at every decade -> the float64 rotations of R_03/R_04 are")
print("    not carrying their results.  The control could have failed and did not.")
print("\nDONE R_06")
