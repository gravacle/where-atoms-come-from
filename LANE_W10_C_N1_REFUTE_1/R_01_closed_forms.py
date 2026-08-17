#!/usr/bin/env python3
"""
R_01 — EVERY CLOSED FORM AND EVERY QUOTED DIGIT IN LANE C, RE-DERIVED INDEPENDENTLY.

Lens: assume each digit is quadrature-limited until shown otherwise.  This script never uses
quadrature.  It uses exact rational arithmetic and mpmath at dps=60, and it re-derives each
closed form from the ready state rather than checking the lane's stated one.

THE GENERAL THEOREM THIS SCRIPT ESTABLISHES, WHICH IS STRICTLY STRONGER THAN LANE C's C-2:

   THEOREM R1 (DOMINANCE).  For non-negative weights with hi*lo >= 0 (lane C's own zero-set
   criterion, NEGATED), one Jensen branch dominates the other at every point of T, hence
        m(P) = log max(p00,p10)    if hi >= 0 and lo >= 0,
        m(P) = log max(p01,p11)    if hi <= 0 and lo <= 0,
   a CLOSED FORM IN ONE LINE, with no reciprocity identity and no reversal locus.
   Conversely if hi*lo < 0 the two branches cross transversally, log max has a corner, and
   m(P) is NOT of that form.

   PROOF.  |A|^2 - |B|^2 is AFFINE in c = cos t with values hi*(p00+p10+p01+p11) at c=+1 and
   lo*(|p00-p10|+|p01-p11|) at c=-1; both bracket factors are >= 0, so the sign of |A|^2-|B|^2
   at the two endpoints is the sign of hi and of lo.  Same sign at both ends + affine =>
   one sign throughout.  Then max(|A|,|B|) is that branch and Jensen in x finishes.  []

   CONSEQUENCE.  m(P) is elementary on the COMPLEMENT of the zero-having set, i.e. on 3/4 of
   the ready-state simplex, and needs quadrature on exactly the 1/4 where P has a torus zero.
   S4:599's "genuinely 4-term and does not factor, QUADRATURE ONLY" is refuted by this one
   line; lane C's reversal-locus mechanism is a MEASURE-ZERO SUB-CASE of it.
"""
from fractions import Fraction as Fr
import mpmath as mp
import numpy as np
import sys

mp.mp.dps = 60
rng = np.random.default_rng(20260816 + 101)

def hi_lo(p):
    return (p[0] + p[1]) - (p[2] + p[3]), abs(p[0] - p[1]) - abs(p[2] - p[3])

def closed_form_or_none(p):
    """THEOREM R1.  Returns (value_as_mpf, which_branch) or (None, None) if branches cross."""
    hi, lo = hi_lo(p)
    if hi >= 0 and lo >= 0:
        q = max(p[0], p[1]); return mp.log(mp.mpf(q.numerator) / q.denominator), "A = p00+p10 x"
    if hi <= 0 and lo <= 0:
        q = max(p[2], p[3]); return mp.log(mp.mpf(q.numerator) / q.denominator), "B = p01+p11 x"
    return None, None

def m_jensen_mp(p, n):
    """Midpoint trapezoid in mpmath.  Reference only; never a reported value on its own."""
    p00, p10, p01, p11 = [mp.mpf(q.numerator) / q.denominator for q in p]
    s = mp.mpf(0)
    for k in range(n):
        c = mp.cos(2 * mp.pi * (k + mp.mpf(1) / 2) / n)
        A2 = p00**2 + p10**2 + 2*p00*p10*c
        B2 = p01**2 + p11**2 + 2*p01*p11*c
        s += mp.log(A2 if A2 > B2 else B2) / 2
    return s / n

CARRIERS = [
    ("B0b  SENSE U", (Fr(4,9), Fr(2,9), Fr(1,9), Fr(2,9)), "-0.810930216216", "S4:599 / C-2"),
    ("B4   SENSE U", (Fr(1,6), Fr(1,6), Fr(1,6), Fr(3,6)), "-0.693147181",    "S4:582 / C-2"),
    ("SENSE C 4cls", (Fr(1,4), Fr(1,4), Fr(1,4), Fr(1,4)), "-1.386294361120", "S4:582 / C-6"),
    ("K1   SENSE U", (Fr(0),   Fr(2,5), Fr(2,5), Fr(1,5)), "-0.756573586",    "S4:582 control"),
    ("B0b* SENSE U", (Fr(2,9), Fr(1,9), Fr(3,9), Fr(3,9)), "-0.987918288038", "C-5"),
    ("pX  (4,2,3,1)/10", (Fr(4,10), Fr(2,10), Fr(3,10), Fr(1,10)), "-0.916291", "C-10"),
]

if __name__ == "__main__":
    print("=" * 100)
    print("R_01 SECTION 1 — THEOREM R1 APPLIED.  WHICH LANE-C VALUES ARE CLOSED FORMS AT ALL?")
    print("=" * 100)
    print(f"  {'ready state':16s} {'hi':>8s} {'lo':>8s} {'hi*lo':>8s}  {'closed form?':>34s}  lane C's quoted digits")
    for name, p, quoted, src in CARRIERS:
        hi, lo = hi_lo(p)
        v, br = closed_form_or_none(p)
        s = f"{mp.nstr(v, 18)}  [{br}]" if v is not None else "NO — branches cross, QUADRATURE"
        print(f"  {name:16s} {str(hi):>8s} {str(lo):>8s} {str(hi*lo):>8s}  {s:>34s}   {quoted}  ({src})")
    print("""
  READ THE TABLE.  FOUR of the six are closed forms by THEOREM R1, including B4, pX and
  SENSE C, NONE of which is on lane C's reversal locus.  (K1 crosses, but K1 is THREE-class:
  with p00 = 0 the polynomial is a monomial times a 3-term one and Cassaigne-Maillot applies,
  which is how S4 got it.  B0b* -- the carrier lane C had to BUILD to make the hypothesis bite
  -- is genuinely four-term AND crossing, and is the one value in the whole lane that needs
  quadrature.)  Lane C's C-2 attributes B0b's exactness to the reversal identity |a+bx|=|b+ax|;
  that identity is TRUE and is NOT the mechanism -- B0b is exact because ONE BRANCH DOMINATES,
  which is C-3's own criterion negated.  The reversal locus has measure ZERO in the simplex;
  the dominance set has measure 3/4.  SO S4:599's sentence "genuinely 4-term and does not
  factor, QUADRATURE ONLY" is refuted for B0b by a condition that holds on 3/4 of ready states,
  and would have been CORRECT on lane C's own B0b*.""")

    print("\n" + "=" * 100)
    print("R_01 SECTION 2 — THE DIGIT AUDIT.  mp.dps = %d.  Exact values to 40 places." % mp.mp.dps)
    print("=" * 100)
    exact = {
        "log(4/9)      B0b":       mp.log(mp.mpf(4)/9),
        "-log 2        B4":        -mp.log(2),
        "-log 4        SENSE C":   -mp.log(4),
        "log(2/5)      pX":        mp.log(mp.mpf(2)/5),
        "(1/4)log(5/243) B0b fin": mp.log(mp.mpf(5)/243)/4,
        "-(3/4) log 3  B4 finite": -3*mp.log(3)/4,
        "-log(9/4)/... K1 fin":    None,
    }
    for k, v in exact.items():
        if v is not None:
            print(f"  {k:26s} = {mp.nstr(v, 40)}")
    print("""
  AGAINST WHAT LANE C AND S4 PRINT:""")
    checks = [
        ("S4:599  B0b  quotes  -0.810930216216 (12 dp)", mp.log(mp.mpf(4)/9), "-0.810930216216"),
        ("S4:582  B4   quotes  -0.693147181    ( 9 dp)", -mp.log(2), "-0.693147181"),
        ("S4:582  SENSE C quotes -1.386294361120 (12dp)", -mp.log(4), "-1.386294361120"),
        ("C-4     B0b finite order -0.970905883", mp.log(mp.mpf(5)/243)/4, "-0.970905883"),
        ("C-4     B4  finite order -0.823959217", -3*mp.log(3)/4, "-0.823959217"),
        ("C-10    pX  m(P)         -0.916291",    mp.log(mp.mpf(2)/5), "-0.916291"),
    ]
    worst = mp.mpf(0)
    for lbl, v, quoted in checks:
        dp = len(quoted.split(".")[1])
        rounded = mp.nstr(v, 25)
        # round the exact value to the same number of decimal places and compare strings
        r = mp.mpf(quoted)
        d = abs(v - r)
        tol = mp.mpf(10) ** (-dp) / 2
        ok = d <= tol
        worst = max(worst, d / tol)
        print(f"    {lbl:48s}  |exact - quoted| = {mp.nstr(d,4):>12s}   "
              f"half-ulp of last printed place = {mp.nstr(tol,3):>10s}   CORRECTLY ROUNDED: {ok}")
    print(f"    WORST (|exact-quoted| / half-ulp) over all six = {mp.nstr(worst,6)}   (must be <= 1)")

    print("\n" + "=" * 100)
    print("R_01 SECTION 3 — THEOREM R1 TESTED AGAINST QUADRATURE ON 200000 RANDOM READY STATES.")
    print("            float64 Jensen at n = 2^16 vs the one-line closed form.  Seed 20260816+101.")
    print("=" * 100)
    NT = 200000
    w = rng.dirichlet(np.ones(4), size=NT)
    HI = (w[:,0]+w[:,1]) - (w[:,2]+w[:,3])
    LO = np.abs(w[:,0]-w[:,1]) - np.abs(w[:,2]-w[:,3])
    dom = HI*LO >= 0
    cf = np.where(HI >= 0, np.log(np.maximum(w[:,0], w[:,1])), np.log(np.maximum(w[:,2], w[:,3])))
    n = 1 << 16
    t = 2*np.pi*(np.arange(n)+0.5)/n; c = np.cos(t)
    def jensen_batch(ww):
        out = np.empty(len(ww))
        for i in range(0, len(ww), 2000):
            b = ww[i:i+2000]
            A2 = (b[:,0]**2+b[:,1]**2)[:,None] + 2*(b[:,0]*b[:,1])[:,None]*c[None,:]
            B2 = (b[:,2]**2+b[:,3]**2)[:,None] + 2*(b[:,2]*b[:,3])[:,None]*c[None,:]
            out[i:i+2000] = np.mean(0.5*np.log(np.maximum(A2, B2)), axis=1)
        return out
    idx_dom = np.where(dom)[0][:20000]
    idx_cross = np.where(~dom)[0][:20000]
    jd = jensen_batch(w[idx_dom]); jc = jensen_batch(w[idx_cross])
    ed = np.abs(jd - cf[idx_dom]); ec = np.abs(jc - cf[idx_cross])
    print(f"  DOMINANCE ARM   (hi*lo >= 0): {len(idx_dom)} draws   max|closed - Jensen(2^16)| = {ed.max():.3e}"
          f"   median = {np.median(ed):.3e}")
    print(f"  CROSSING ARM    (hi*lo <  0): {len(idx_cross)} draws   max|WRONG formula - Jensen| = {ec.max():.3e}"
          f"   median = {np.median(ec):.3e}")
    print(f"  fraction of the {NT} draws in the crossing (quadrature-needed) arm = {1-dom.mean():.6f}"
          f"   [Theorem: EXACTLY 1/4]")
    print("""  The dominance arm agrees to quadrature accuracy on every draw; the crossing arm does
  not, by a median of ~1e-2.  THE ARMS ARE NOT BYTE-IDENTICAL AND THE CONTROL IS NOT VACUOUS:
  the SAME formula is applied in both and only the criterion's sign moves.""")
    assert ed.max() < 1e-12, "THEOREM R1 FAILS ON THE DOMINANCE ARM"
    assert np.median(ec) > 1e-4, "CROSSING ARM INDISTINGUISHABLE — CONTROL WOULD BE VACUOUS"

    print("\n" + "=" * 100)
    print("R_01 SECTION 4 — HOW BAD IS THE TRAPEZOID WHERE QUADRATURE IS ACTUALLY NEEDED?")
    print("            The corner (crossing) gives n^-2; the log singularity gives n^-1.")
    print("=" * 100)
    B0bs = (Fr(2,9), Fr(1,9), Fr(3,9), Fr(3,9))
    SC   = (Fr(1,4), Fr(1,4), Fr(1,4), Fr(1,4))
    # high-accuracy reference for B0b*: split the integral at the corner and use Gauss-Legendre
    p00, p10, p01, p11 = [mp.mpf(q.numerator)/q.denominator for q in B0bs]
    A0 = p00**2 + p10**2 - p01**2 - p11**2
    B0 = 2*(p00*p10 - p01*p11)
    ccross = -A0/B0
    tc = mp.acos(ccross)
    def f(t):
        c = mp.cos(t)
        A2 = p00**2 + p10**2 + 2*p00*p10*c
        B2 = p01**2 + p11**2 + 2*p01*p11*c
        return mp.log(A2 if A2 > B2 else B2)/2
    ref = (mp.quad(f, [0, tc]) + mp.quad(f, [tc, mp.pi]))/mp.pi
    print(f"  B0b* corner at cos t = {mp.nstr(ccross,20)}  (t = {mp.nstr(tc,20)})")
    print(f"  B0b* m(P) by CORNER-SPLIT Gauss-Legendre, dps={mp.mp.dps}: {mp.nstr(ref, 30)}")
    print(f"       lane C prints -0.987918288038 ; difference = {mp.nstr(abs(ref+mp.mpf('0.987918288038')),4)}")
    print(f"\n  {'n':>10s}   {'B0b* corner: |trap - ref|':>28s}   {'SENSE C log-sing: |trap - (-log4)|':>36s}")
    for e in (8, 10, 12, 14, 16, 18, 20):
        n = 2**e
        t = 2*np.pi*(np.arange(n)+0.5)/n; cc = np.cos(t)
        def trap(p):
            q = [float(z) for z in p]
            A2 = q[0]**2+q[1]**2+2*q[0]*q[1]*cc
            B2 = q[2]**2+q[3]**2+2*q[2]*q[3]*cc
            return float(np.mean(0.5*np.log(np.maximum(A2, B2))))
        e1 = abs(trap(B0bs) - float(ref)); e2 = abs(trap(SC) - float(-mp.log(4)))
        print(f"  2^{e:<2d}={n:>7d}   {e1:28.3e}   {e2:36.3e}")
    print("""  Ratios halve per doubling for SENSE C (n^-1: log singularity) and quarter for B0b*
  (n^-2: Lipschitz corner).  THREE regimes, and lane C's self-flag calls the taxonomy 'fitted,
  not proved'.  It is not fitted here: the corner is EXHIBITED at cos t above, and the
  singularity is the whole circle x = -1.""")
    sys.exit(0)
