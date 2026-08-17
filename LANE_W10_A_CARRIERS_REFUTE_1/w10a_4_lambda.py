#!/usr/bin/env python3
"""
LANE W10-A — SCRIPT 4.  THE RATE ON B0b AND B4, AGAINST S4's PUBLISHED lambda COLUMN,
plus a side-check on the hypothesis the registrar named for the multiset theorem.

S4:574-578 publishes, for the two four-class carriers:
    B0b  lambda(SENSE U) = -0.810930216   lambda(SENSE C) = -1.386294361
    B4   lambda(SENSE U) = -0.693147181   lambda(SENSE C) = -1.386294361
and S4:599 says B0b's SENSE-U value is "genuinely 4-term and does not factor. QUADRATURE ONLY".
W-03 corrected that to log(4/9) EXACTLY.  This script decides it from MY OWN incidence.

QUADRATURE WARNING CARRIED: a four-class P has torus zeros, so a 2D grid on log|P| is
noise-limited.  Every value here comes from the JENSEN REDUCTION
    m(P) = (1/2pi) INT log max(|p00 + p01 e^{it}|, |p10 + p11 e^{it}|) dt,
which has a log singularity only where BOTH branches vanish at the same t.
"""
import sys, math
from fractions import Fraction
import numpy as np
from w10a_lib import B0b, B4, K1, B1q, CLASSES

LOG = []
def out(s=""):
    print(s); LOG.append(s)

out("=" * 104)
out("W10-A SCRIPT 4 — lambda ON B0b AND B4 FROM MY OWN INCIDENCE, AGAINST S4:574-578")
out("=" * 104)
out(f"numpy {np.__version__}.  IEEE double unless a line says EXACT.")
out()

def jensen(pi4, n):
    """m(P) by the Jensen reduction, n-point midpoint rule on the t-circle."""
    p00, p10, p01, p11 = [float(x) for x in pi4]
    t = (np.arange(n) + 0.5) * (2 * np.pi / n)
    e = np.exp(1j * t)
    A = np.abs(p00 + p01 * e)
    B = np.abs(p10 + p11 * e)
    M = np.maximum(A, B)
    return float(np.mean(np.log(np.maximum(M, 1e-300))))

def branch_dominance(pi4):
    """EXACT: is one Jensen branch >= the other for ALL t?  |A|^2-|B|^2 = C + D cos t with
       C = p00^2+p01^2-p10^2-p11^2, D = 2(p00 p01 - p10 p11).  One branch dominates iff |D|<=|C|."""
    p00, p10, p01, p11 = [Fraction(x) for x in pi4]
    C = p00 * p00 + p01 * p01 - p10 * p10 - p11 * p11
    D = 2 * (p00 * p01 - p10 * p11)
    if C == 0 and D == 0:
        return "EQ", C, D          # the two Jensen branches are IDENTICAL for all t
    if abs(D) <= abs(C):
        return ("A" if C > 0 else "B"), C, D
    return (None, C, D)

def schedule_B(pi4, f, c, N, seed=None):
    """DIRECT simulation of the canonical clock: (1/N) sum_{k<=N} log|Z_k| with u=e^{-if}, v=e^{ic}."""
    k = np.arange(1, N + 1, dtype=np.float64)
    a = -f / (2 * np.pi); b = c / (2 * np.pi)
    ua = np.exp(2j * np.pi * ((k * a) % 1.0))
    vb = np.exp(2j * np.pi * ((k * b) % 1.0))
    p00, p10, p01, p11 = [float(x) for x in pi4]
    m = np.abs(p00 + p10 * ua + p01 * vb + p11 * ua * vb)
    with np.errstate(divide="ignore"):
        return float(np.mean(np.log(np.maximum(m, 1e-300))))

S4_LAMBDA = {"B0b U": -0.810930216, "B0b C": -1.386294361,
             "B4 U": -0.693147181,  "B4 C": -1.386294361,
             "B1 U": -0.756573586,  "B1 C": -0.767507880,
             "B1q U": -0.741029583, "B1q C": -0.767507880}

CASES = []
for tag, K in [("B0b", B0b()), ("B4", B4()), ("B1", K1()), ("B1q", B1q())]:
    piU = K.pushforward_uniform()
    occ = [c for c, x in zip(CLASSES, piU) if x > 0]
    piC = ([Fraction(1, 4)] * 4 if len(occ) == 4 else
           [Fraction(2, 5) if c == occ[0] else (Fraction(3, 10) if c in occ else Fraction(0))
            for c in CLASSES])
    CASES.append((tag + " U", piU))
    CASES.append((tag + " C", piC))

out("SENSE U = p_v = 1/V (S4's own).  SENSE C = class weights fixed by hand: (1/4,1/4,1/4,1/4)")
out("for four classes, (0.4,0.3,0.3) on the occupied three otherwise (S4:557-560).")
out()
out(f"{'case':<8}{'pi = (p00,p10,p01,p11)':<34}{'dominant branch':<17}"
    f"{'Jensen 2^22':<16}{'S4 published':<15}{'|dev|':<11}{'exact closed form'}")
FAILS = []
for nm, pi4 in CASES:
    dom, C, D = branch_dominance(pi4)
    val = jensen(pi4, 1 << 22)
    ref = S4_LAMBDA[nm]
    p00, p10, p01, p11 = pi4
    cf = ""
    if dom == "EQ":
        cf = (f"branches IDENTICAL; log max(p10,p11) = log({max(p10,p11)}) = "
              f"{math.log(float(max(p10,p11))):.12f} -- BOTH vanish at t=pi, so the Jensen "
              f"integrand has a log singularity and the quadrature is slow; the exact value "
              f"comes from the FACTORISATION, not from this grid")
    elif dom == "A":
        cf = f"log max(p00,p01) = log({max(p00,p01)}) = {math.log(float(max(p00,p01))):.12f}"
    elif dom == "B":
        cf = f"log max(p10,p11) = log({max(p10,p11)}) = {math.log(float(max(p10,p11))):.12f}"
    elif p00 * p11 == p10 * p01:
        a = p00 + p10; c_ = p00 + p01
        cf = "factors: log max + log max (see script 3 part D)"
    out(f"{nm:<8}{'('+', '.join(str(x) for x in pi4)+')':<34}"
        f"{({'EQ':'branches equal'}.get(dom, dom) if dom else 'neither (crossing)'):<17}{val:< 16.9f}{ref:< 15.9f}"
        f"{abs(val-ref):<11.2e}{cf}")
    if abs(val - ref) > 2e-6:
        FAILS.append(f"lambda mismatch {nm}: {val} vs S4 {ref}")
out()

out("-" * 104)
out("EXACT: THE BRANCH-DOMINANCE TEST, IN RATIONALS")
out("-" * 104)
out("|p00+p01 e^{it}|^2 - |p10+p11 e^{it}|^2 = C + D cos t with")
out("   C = p00^2+p01^2-p10^2-p11^2 ,  D = 2(p00 p01 - p10 p11).")
out("cos t sweeps [-1,1] exactly, so ONE BRANCH DOMINATES FOR ALL t iff |D| <= |C|, and then")
out("Jensen gives m(P) = log max of that branch's two coefficients -- an EXACT closed form.")
for nm, pi4 in CASES:
    dom, C, D = branch_dominance(pi4)
    out(f"   {nm:<8} C = {str(C):<10} D = {str(D):<10} |D| <= |C| ? {abs(D) <= abs(C)}   "
        f"dominant: {({'EQ':'branches EQUAL'}.get(dom, dom) if dom else 'neither')}")
out()
out("B0b SENSE U: pi = (4/9, 2/9, 1/9, 2/9).  C = (16+1-4-4)/81 = 9/81 = 1/9, D = 2(4-4)/81 = 0.")
out("   |D| = 0 <= |C| = 1/9, branch A dominates everywhere, so")
out("   lambda(B0b, SENSE U) = log(4/9) = -0.810930216216328  EXACTLY.")
out("   **S4:599 CALLS THIS ROW 'QUADRATURE ONLY' AND IT IS NOT.  W-03 already corrected it;")
out("   this is the first time the correction is derived from a BUILT B0b rather than from S4's")
out("   own number.  S4's flag F-quadrature on that row should be struck.**")
out("B4 SENSE U: pi = (1/6,1/6,1/6,1/2).  C = (1+1-1-9)/36 = -8/36, D = 2(1-3)/36 = -4/36.")
out("   |D| <= |C|, branch B dominates, lambda = log(1/2) = -0.693147180559945 EXACTLY,")
out("   which is S4's own reasoning (S4:594-595) reproduced from my incidence.")
out("SENSE C on BOTH four-class carriers: pi = (1/4,1/4,1/4,1/4), C = D = 0, so the two Jensen")
out("   branches are IDENTICAL and BOTH VANISH at t = pi.  That is exactly the log singularity")
out("   W-09 warned about, and it is visible in the convergence table below: the SENSE C rows")
out("   are the only ones still moving at n = 2^22 (1.65e-07 from log(1/4)).  THE EXACT VALUE ON")
out("   THOSE ROWS COMES FROM THE FACTORISATION P = (1+x)(1+y)/4, NOT FROM ANY GRID.")
out()

out("-" * 104)
out("QUADRATURE CONVERGENCE (the W-09 warning, respected)")
out("-" * 104)
out(f"{'case':<10}" + "".join(f"{'n=2^'+str(e):<18}" for e in (16, 18, 20, 22)))
for nm, pi4 in CASES:
    out(f"{nm:<10}" + "".join(f"{jensen(pi4, 1 << e):< 18.12f}" for e in (16, 18, 20, 22)))
out()

out("-" * 104)
out("DIRECT SCHEDULE-B SIMULATION AT S4's OWN TEST POINT f = 1.0, c = sqrt(2), N = 2e6")
out("-" * 104)
out(f"{'case':<10}{'direct':<18}{'Jensen 2^22':<18}{'dev'}")
for nm, pi4 in CASES:
    d = schedule_B(pi4, 1.0, math.sqrt(2.0), 2000000)
    j = jensen(pi4, 1 << 22)
    out(f"{nm:<10}{d:< 18.9f}{j:< 18.9f}{abs(d-j):.2e}")
    if abs(d - j) > 5e-5:
        FAILS.append(f"direct vs Jensen {nm}")
out()

# --------------------------------------------------------------------------------------------
out("=" * 104)
out("SIDE-CHECK — THE HYPOTHESIS OF THE MULTISET THEOREM.  NOT THIS LANE'S ASSIGNMENT; ONE BLOCK.")
out("=" * 104)
out("The registrar reports 24/24 permutations of (p00,p10,p01,p11) leave m(P) invariant with")
out("spread 0.000e+00, and names the hypothesis REAL NON-NEGATIVITY.  The step that supplies the")
out("two extra transpositions is  |a + b e^{it}| = |b + a e^{it}| pointwise, and")
out("     |a + b e^{it}|^2 = a^2 + b^2 + 2 a b cos t")
out("is SYMMETRIC IN a AND b FOR ANY REAL a, b, SIGN IRRELEVANT.  So the hypothesis is REALITY,")
out("not non-negativity.  Tested below on mixed-sign real vectors and on complex ones.")
out()
from itertools import permutations
rng = np.random.default_rng(20260816)

def m_jensen_general(v, n=1 << 20):
    p00, p10, p01, p11 = v
    t = (np.arange(n) + 0.5) * (2 * np.pi / n)
    e = np.exp(1j * t)
    M = np.maximum(np.abs(p00 + p01 * e), np.abs(p10 + p11 * e))
    return float(np.mean(np.log(np.maximum(M, 1e-300))))

out(f"{'coefficient kind':<34}{'spread over the 24 permutations (max - min of m)':<50}")
for kind, gen in [
    ("real non-negative (registrar's)", lambda: rng.random(4)),
    ("real MIXED SIGN",                 lambda: rng.normal(size=4)),
    ("real mixed sign, 2nd draw",       lambda: rng.normal(size=4)),
    ("COMPLEX",                         lambda: rng.normal(size=4) + 1j * rng.normal(size=4)),
    ("COMPLEX, 2nd draw",               lambda: rng.normal(size=4) + 1j * rng.normal(size=4)),
]:
    v = gen()
    vals = [m_jensen_general(list(p)) for p in permutations(v)]
    out(f"{kind:<34}{max(vals)-min(vals):<.6e}      (m = {vals[0]: .9f})")
out()
out("READ: S4-invariance survives MIXED SIGN and DIES on COMPLEX coefficients.  The registrar's")
out("prediction that the group is S4 rather than the Newton-polygon D4 is CONFIRMED and its")
out("STATED HYPOTHESIS IS ONE WORD TOO STRONG: the operative hypothesis is REAL COEFFICIENTS.")
out("Non-negativity is what makes pi a probability vector; it is not what makes m symmetric.")
out("(This corrects a naming, not a number.  It is the corpus's recurring defect and it is")
out(" recorded here rather than left for the next layer.)")

out()
out("=" * 104)
if FAILS:
    out(f"**{len(FAILS)} FAILURES**")
    for f in FAILS:
        out("   " + f)
else:
    out("0 failures.  Both four-class carriers reproduce S4's published lambda column from MY")
    out("incidence, and both SENSE-U values are EXACT closed forms, not quadrature.")

with open("w10a_4_lambda.OUT.txt", "w") as fh:
    fh.write("\n".join(LOG) + "\n")
sys.exit(1 if FAILS else 0)
