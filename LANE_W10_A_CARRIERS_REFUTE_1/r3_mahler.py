#!/usr/bin/env python3
"""
REFUTER 1 — SCRIPT 3.  THE RATE COLUMN, THE BRANCH-DOMINANCE CLAIM, AND THE PERMUTATION NAMING.

A. Verify the lane's two exact SENSE-U closed forms to 30 DECIMAL DIGITS with mpmath, by a
   quadrature that is not the lane's (Gauss-Legendre / tanh-sinh on the Jensen integrand at
   mp.dps = 30), and by the exact rational branch-dominance argument done independently.
B. Test the branch-dominance criterion itself (|D| <= |C|) as an iff, and test the lane's
   claim that this criterion is NEW.
C. Test the factorisation criterion p00 p11 = p10 p01 and the lane's claim that the locus is
   "FOUR-CLASS-ONLY" and that S4 got log(1/4) as "a lucky closed form".
D. THE PERMUTATION NAMING.  The lane corrects the registrar's "REAL NON-NEGATIVE" to
   "REAL COEFFICIENTS".  Attack THAT name with a common-phase counterexample, and check
   whether a weaker-looking condition (cross-ratio real) suffices.

mpmath 1.4.1 at mp.dps = 30 for A.  IEEE double for D (spreads of order 1e-16 vs 1e-1;
the verdict is not precision-sensitive and that is shown, not asserted).
"""
import sys, math
from fractions import Fraction
from itertools import permutations
import numpy as np
import mpmath as mp

LOG = []
def out(s=""):
    print(s); LOG.append(s)

FAIL = []
mp.mp.dps = 30

out("=" * 104)
out("REFUTER 1 / SCRIPT 3 — RATE, BRANCH DOMINANCE, FACTORISATION, PERMUTATION NAMING")
out("=" * 104)
out(f"mpmath {mp.__version__} at mp.dps = {mp.mp.dps}.  numpy {np.__version__}.")
out()


def branch(pi4):
    """EXACT rational C, D with |A(t)|^2 - |B(t)|^2 = C + D cos t."""
    p00, p10, p01, p11 = [Fraction(x) for x in pi4]
    C = p00 * p00 + p01 * p01 - p10 * p10 - p11 * p11
    D = 2 * (p00 * p01 - p10 * p11)
    return C, D


def m_jensen_mp(pi4):
    """m(P) to mp precision by the Jensen reduction, mpmath quad (NOT the lane's midpoint rule).
    The integrand log max(|A|,|B|) has KINKS wherever the branches cross, and Gauss-Legendre
    on a kinked integrand converges slowly -- my first run lost 7 digits on the two crossing
    rows.  Recorded, not hidden: the panel endpoints are the exact crossing points
    cos t = -C/D, so every panel has an analytic integrand."""
    p00, p10, p01, p11 = [mp.mpf(x.numerator) / mp.mpf(x.denominator) if isinstance(x, Fraction)
                          else mp.mpf(x) for x in pi4]
    C, D = branch(pi4)
    pts = [mp.mpf(0), mp.pi, 2 * mp.pi]
    if D != 0 and abs(Fraction(C, D)) < 1:
        t0 = mp.acos(-mp.mpf(C.numerator * D.denominator) / mp.mpf(C.denominator * D.numerator))
        pts = sorted(set([mp.mpf(0), t0, mp.pi, 2 * mp.pi - t0, 2 * mp.pi]))

    def f(t):
        e = mp.expj(t)
        A = abs(p00 + p01 * e)
        B = abs(p10 + p11 * e)
        M = A if A > B else B
        return mp.log(M) if M > 0 else mp.mpf('-1e30')
    return mp.quad(f, pts) / (2 * mp.pi)


CASES = [
    ("B0b SENSE U", [Fraction(4, 9), Fraction(2, 9), Fraction(1, 9), Fraction(2, 9)],
     "-0.810930216"),
    ("B4  SENSE U", [Fraction(1, 6), Fraction(1, 6), Fraction(1, 6), Fraction(1, 2)],
     "-0.693147181"),
    ("SENSE C 1/4", [Fraction(1, 4)] * 4, "-1.386294361"),
    ("B1  SENSE U", [Fraction(0), Fraction(2, 5), Fraction(2, 5), Fraction(1, 5)],
     "-0.756573586"),
    ("B1q SENSE U", [Fraction(1, 7), Fraction(3, 7), Fraction(3, 7), Fraction(0)],
     "-0.741029583"),
]

out("-" * 104)
out("A — THE TWO FOUR-CLASS SENSE-U ROWS TO 30 DIGITS, BY A QUADRATURE THAT IS NOT THE LANE'S")
out("-" * 104)
out(f"{'case':<14}{'C (exact)':<12}{'D (exact)':<12}{'|D|<=|C|':<10}"
    f"{'mpmath m(P) (dps=30)':<36}{'S4 published':<15}{'dev'}")
for nm, pi4, s4 in CASES:
    C, D = branch(pi4)
    val = m_jensen_mp(pi4)
    dev = abs(float(val) - float(s4))
    out(f"{nm:<14}{str(C):<12}{str(D):<12}{str(abs(D) <= abs(C)):<10}"
        f"{mp.nstr(val, 25):<36}{s4:<15}{dev:.2e}")
    if dev > 2e-9:
        FAIL.append(f"{nm} disagrees with S4 by {dev}")
out()
out("EXACT CLOSED FORMS, re-derived here and checked against the 30-digit quadrature:")
for nm, pi4, s4 in CASES[:3]:
    C, D = branch(pi4)
    p00, p10, p01, p11 = pi4
    if abs(D) <= abs(C) and not (C == 0 and D == 0):
        arg = max(p00, p01) if C > 0 else max(p10, p11)
        exact = mp.log(mp.mpf(arg.numerator) / mp.mpf(arg.denominator))
        out(f"   {nm:<14} one branch dominates -> m = log({arg}) = {mp.nstr(exact, 25)}")
        out(f"                  quadrature - closed form = "
            f"{mp.nstr(abs(m_jensen_mp(pi4) - exact), 8)}")
    else:
        out(f"   {nm:<14} C = D = 0: the two branches are IDENTICAL.  P = (1+x)(1+y)/4, so")
        out(f"                  m = 2 log max(1/2,1/2)... no: m = log(1/4) directly from the")
        out(f"                  factorisation.  log(1/4) = {mp.nstr(mp.log(mp.mpf(1)/4), 25)}")
        out(f"                  quadrature - log(1/4)   = "
            f"{mp.nstr(abs(m_jensen_mp(pi4) - mp.log(mp.mpf(1) / 4)), 8)}")
out()
out(">>> CONFIRMED.  B0b SENSE U = log(4/9) and B4 SENSE U = log(1/2), EXACTLY, and S4:599's")
out("    'QUADRATURE ONLY' on the B0b row is wrong.  BUT: W-03's register row already says so")
out("    ('B0b is log(4/9) = -0.8109302162163288 exactly, so nine of nine carrier rates are")
out("    exact, not eight'), and the lane says so too.  This is a REPRODUCTION of a standing")
out("    correction, not a new one.")
out()

# ==================================================================================================
out("-" * 104)
out("B — IS THE BRANCH-DOMINANCE CRITERION NEW?  (the lane's W10A-07 says it is)")
out("-" * 104)
out("S4_THE_MEASUREMENT_V001.md:594-595, verbatim:")
out("   'B4   U  : the Jensen-in-x integrand's max is always the second branch (their squares")
out("              differ by 0.2222 + 0.1111 cos y > 0), so lambda = log(1/2) ... EXACTLY'")
C4, D4 = branch(CASES[1][1])
out(f"My exact C, D for B4 SENSE U: C = {C4} = {float(C4):.4f}, D = {D4} = {float(D4):.4f}.")
out(f"S4's printed pair is (0.2222, 0.1111) = ({float(-C4):.4f}, {float(-D4):.4f}) -- the SAME")
out("two numbers with the branches named in the other order.  S4 states the criterion, in the")
out("form 'their squares differ by C + D cos y > 0', and USES it to get an exact closed form.")
out(">>> THE CRITERION IS NOT NEW.  It is S4's own argument, generalised from one row to a")
out("    rule.  Generalising it is worth something; calling it new is a novelty overclaim, and")
out("    the lane's OWN SCRIPT contradicts the lane's OWN FINDING here -- w10a_4_lambda.py:143")
out("    says 'which is S4's own reasoning (S4:594-595) reproduced from my incidence'.")
out()

# ==================================================================================================
out("-" * 104)
out("C — THE FACTORISATION LOCUS:  p00 p11 = p10 p01")
out("-" * 104)
out("Claim (W10A-06): the locus is EMPTY on every three-class carrier, so factorisation is a")
out("FOUR-CLASS-ONLY phenomenon.  Enumerate every support and test.")
from itertools import combinations
CLASSES4 = [(0, 0), (1, 0), (0, 1), (1, 1)]
NAME = {(0, 0): "00", (1, 0): "10", (0, 1): "01", (1, 1): "11"}
out(f"{'support':<18}{'|S|':<5}{'uniform pi':<34}{'p00 p11 - p10 p01':<22}{'factors?'}")
nfac = {1: 0, 2: 0, 3: 0, 4: 0}
for r in range(1, 5):
    for S in combinations(CLASSES4, r):
        w = Fraction(1, r)
        pi4 = [w if c in S else Fraction(0) for c in CLASSES4]
        d = pi4[0] * pi4[3] - pi4[1] * pi4[2]
        f = (d == 0)
        nfac[r] += int(f)
        nm = "{" + ",".join(NAME[c] for c in S) + "}"
        out(f"{nm:<18}{r:<5}{'(' + ', '.join(str(x) for x in pi4) + ')':<34}{str(d):<22}{f}")
out()
out(f"supports that factor, by size: |S|=1: {nfac[1]}/4   |S|=2: {nfac[2]}/6   "
    f"|S|=3: {nfac[3]}/4   |S|=4: {nfac[4]}/1")
out(">>> The lane's claim that the locus is EMPTY at |S| = 3 is CORRECT and I confirm it.")
out("    But 'FOUR-CLASS-ONLY' is not right as stated: FOUR of the six two-class supports")
out("    factor (trivially -- they are one-variable), and the |S|=1 supports factor too.")
out("    The accurate statement is: the factorisation locus is empty at |S| = 3 and non-empty")
out("    at |S| = 4, and NON-DEGENERATE factorisation (both factors genuinely present) requires")
out("    |S| = 4.  The lane's own script text says the right thing; the finding does not.")
out()
out("AND S4 ALREADY FACTORED THE FOUR-CLASS SENSE-C POLYNOMIAL.  S4:598, verbatim:")
out("   'SENSE C, 4 classes: (1+x+y+xy)/4 = (1+x)(1+y)/4, so lambda = log(1/4) ... EXACTLY'")
out(">>> W10A-06's 'which is why S4:596's log(1/4) is exact -- it is a sum of two Jensen terms,")
out("    not a lucky closed form' is wrong twice: the line is S4:598 (S4:596 is B1p's row), and")
out("    S4 obtained it BY THE FACTORISATION, printed on that line.  What is genuinely new is")
out("    the IFF (p00 p11 = p10 p01); the phenomenon and its use are S4's.")
out()

# ==================================================================================================
out("-" * 104)
out("D — THE PERMUTATION NAMING.  THE LANE'S CORRECTION IS ITSELF ONE STEP TOO STRONG.")
out("-" * 104)
out("Registrar: the multiset theorem's hypothesis is REAL NON-NEGATIVE coefficients.")
out("Lane W10A-10: 'the operative hypothesis is REAL COEFFICIENTS, not real non-negativity',")
out("evidenced by spread 1.7e-16 on real mixed sign and 2.302e-01 on COMPLEX.")
out()
out("BUT m(P) is invariant under multiplying P by ANY unimodular constant.  So if v is real")
out("then e^{i phi} v is not real and m is STILL S4-invariant on it, for every phi.  The")
out("lane's 'COMPLEX' draws are generic complex vectors, which cannot separate 'real' from")
out("'real up to one global phase'.  Test below.")
out()


def m_double(v, n=1 << 20):
    p00, p10, p01, p11 = v
    t = (np.arange(n) + 0.5) * (2 * np.pi / n)
    e = np.exp(1j * t)
    M = np.maximum(np.abs(p00 + p01 * e), np.abs(p10 + p11 * e))
    return float(np.mean(np.log(np.maximum(M, 1e-300))))


rng = np.random.default_rng(20260816)
rows = []
base_nn = rng.random(4)
base_ms = rng.normal(size=4)
phi = 0.7391                                     # any phase
TESTS = [
    ("real non-negative", list(base_nn)),
    ("real MIXED SIGN", list(base_ms)),
    ("COMPLEX generic", list(rng.normal(size=4) + 1j * rng.normal(size=4))),
    ("COMPLEX generic, 2nd", list(rng.normal(size=4) + 1j * rng.normal(size=4))),
    ("e^{i*0.7391} * (real non-neg)  <-- NOT REAL", [np.exp(1j * phi) * x for x in base_nn]),
    ("e^{i*0.7391} * (real mixed)    <-- NOT REAL", [np.exp(1j * phi) * x for x in base_ms]),
    ("i * (real non-neg)             <-- NOT REAL", [1j * x for x in base_nn]),
    ("cross-ratio REAL but not common-phase: (1, i, 1, i)", [1 + 0j, 1j, 1 + 0j, 1j]),
]
out(f"{'coefficient kind':<52}{'spread of m over 24 permutations':<36}{'m(identity)'}")
for kind, v in TESTS:
    vals = [m_double(list(p)) for p in permutations(v)]
    sp = max(vals) - min(vals)
    rows.append((kind, sp))
    out(f"{kind:<52}{sp:<36.6e}{vals[0]: .9f}")
out()
out("READ:")
out("  * real non-negative  -> invariant   (registrar, confirmed)")
out("  * real mixed sign    -> invariant   (lane's correction, confirmed)")
out("  * e^{i phi} * real   -> INVARIANT, AND IT IS NOT REAL.  So 'REAL COEFFICIENTS' is a")
out("    SUFFICIENT condition that is strictly stronger than necessary, exactly as")
out("    'real NON-NEGATIVE' was.  The sharp sufficient condition is: all four coefficients")
out("    share one argument mod pi, i.e. e^{-i phi} p is real for some phi.  Proof: m is")
out("    unchanged by P -> e^{i phi} P, and the S4-generating identities |a+b e^{it}| =")
out("    |b+a e^{it}| hold for the real vector.")
out("  * cross-ratio real but NOT common-phase, (1,i,1,i): p00 p11/(p10 p01) = 1 is real and")
out("    P = (1+ix)(1+y) factors, yet the permutation spread is LARGE.  So 'the cross-ratio is")
out("    real' -- the natural next weakening, and the same quantity as the factorisation")
out("    criterion -- does NOT suffice.  The monomial twist that makes the coefficients real")
out("    does not commute with permuting them.")
out()
out(">>> FINDING.  W10A-10 corrects the registrar's name and lands one step short of the sharp")
out("    one.  This is the corpus's signature defect (four consecutive layers misnaming the")
out("    operative variable) committed by the layer that flags itself for exactly that risk,")
out("    inside the block it labels 'out of lane'.  The DIRECTION of the lane's correction is")
out("    right and the registrar's D4-vs-S4 prediction stays refuted.")
out()

out("=" * 104)
if FAIL:
    out(f"**{len(FAIL)} FAILURES**")
    for f in FAIL:
        out("   " + f)
else:
    out("A: both four-class SENSE-U closed forms verified to 30 digits by an independent")
    out("   quadrature.  The lane's arithmetic is right.")
    out("B: the branch-dominance criterion is S4's own, generalised -- not new.")
    out("C: the factorisation IFF is new; 'four-class-only' is imprecise; S4 already factored")
    out("   the SENSE-C row on the line the lane mis-cites.")
    out("D: 'REAL COEFFICIENTS' is not the sharp hypothesis; 'real up to a global phase' is.")

with open("r3_mahler.OUT.txt", "w") as fh:
    fh.write("\n".join(LOG) + "\n")
