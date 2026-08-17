"""LANE W-13 / Z_REFUTE -- r3: THE CUSTODY AND SCOPE ATTACK.
Lane Z's Z-13 convicts the program of UNDER-READING: a sealed correction sat unread for eleven
lane-directories.  This script asks whether lane Z under-read, and whether its own corpus rows
say what it says they say.  Three findings, each checked at the bytes:
  (1) B0b's class multiset is misquoted from S4:575, and the misquote converts the corpus's
      flagship NON-FACTORING four-term state into one that FACTORS.
  (2) Theorem Z4 (the anti-diagonal) and the closed form cos s0 = -2/3 depend on a LABELLING
      that S4's own definition of SENSE C does not fix.
  (3) THE HEADLINE'S SIGN IS INVERTED, and the inversion is written out in the very lane Z
      cites for priority.  Lane Z's own citation, M3-2(d), says the opposite in one line."""
import sys, math
import numpy as np
from fractions import Fraction as F
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from r0_lib import (strat_sorted, zero_angles, zero_points, mahler_1d, jensen_branches,
                    has_torus_zero_invariant, fr, S4_575, K1_REG, S1_PUB, SENSEC4,
                    B0b_LANE, B0b_S4, CENTROID, Pabs)

W = 96
def hdr(s):
    print("=" * W); print(s); print("=" * W)

hdr("r3  WHAT LANE Z's CORPUS ROWS ACTUALLY SAY, AND WHICH WAY ITS HEADLINE POINTS")
print("numpy", np.__version__, "\n")

# =========================================================================================
print("-" * W)
print("(a) DEFECT 1 -- B0b's MULTISET IS MISQUOTED FROM S4:575, AND THE MISQUOTE FACTORS")
print("    THE POLYNOMIAL THAT S4 AND W-10 N-3 BOTH CERTIFY AS NON-FACTORING.\n")
print("    S4_THE_MEASUREMENT_V001.md:575 row B0b reads   {00:4, 01:1, 10:2, 11:2}")
print("      -> pi = (p00,p10,p01,p11) = (4/9, 2/9, 1/9, 2/9)")
print("    LANE_W13_Z_ZEROSET/z0_lib.py NAMED entry reads {00:4,10:2,01:2,11:1}")
print("      -> pi = (4/9, 2/9, 2/9, 1/9)        <- p01 and p11 TRANSPOSED\n")
for lbl, p in (("S4:575 as written ", B0b_S4), ("LANE Z as coded   ", B0b_LANE)):
    p00, p10, p01, p11 = p
    print(f"    {lbl} pi = {tuple(str(q) for q in p)}   p00*p11 = {p00*p11}   "
          f"p10*p01 = {p10*p01}   FACTORS: {p00*p11 == p10*p01}")
print()
print("    W-10 defect N-3, REGISTER_V001.md (W-10 row), verbatim: \"The non-factoring half is")
print("    exact (p00*p11 = 8/81 != 2/81 = p10*p01)\".  That arithmetic is TRUE of S4's state")
print("    and FALSE of lane Z's, whose two cross-products are BOTH 4/81.  Lane Z's z3(e)")
print("    quotes N-3 as \"the counterexample the corpus already owns\" while running a state")
print("    for which N-3's own certificate of non-factoring fails.")
print()
print("    LANE Z's STATE FACTORS EXPLICITLY:  9*P = 4 + 2x + 2y + xy = (2+x)(2+y).")
xs = np.exp(1j * np.linspace(0, 2 * np.pi, 997, endpoint=False))
ys = np.exp(1j * np.linspace(0, 2 * np.pi, 991, endpoint=False))
X, Y = np.meshgrid(xs, ys)
lhs = np.abs(4 + 2 * X + 2 * Y + X * Y)
rhs = np.abs((2 + X) * (2 + Y))
print(f"      max | 4+2x+2y+xy - (2+x)(2+y) | over a 997x991 torus grid = "
      f"{float(np.max(np.abs(lhs-rhs))):.3e}")
print()
print("    WHY IT SURVIVED: the transposition is WITHIN the pair {p01,p11}, and the Jensen")
print("    branch |B(t)|^2 = p01^2+p11^2+2 p01 p11 cos t is SYMMETRIC in that pair.  So")
print("    |A|, |B|, min|P|, m(P) and the stratum are IDENTICAL for the two states, and every")
print("    figure lane Z computes from B0b is unaffected.  Checked, not asserted:")
t = np.linspace(0, 2 * np.pi, 200001)
a1, b1 = jensen_branches(B0b_S4, t); a2, b2 = jensen_branches(B0b_LANE, t)
print(f"      max | |A|_S4 - |A|_laneZ | = {float(np.max(np.abs(a1-a2))):.3e}")
print(f"      max | |B|_S4 - |B|_laneZ | = {float(np.max(np.abs(b1-b2))):.3e}")
print(f"      min_T^2 |P| : S4 {float(np.min(np.abs(a1-b1))):.12f}   "
      f"lane Z {float(np.min(np.abs(a2-b2))):.12f}   (both 1/9 = {1/9:.12f})")
print(f"      m(P)        : S4 {mahler_1d(B0b_S4):.12f}   lane Z {mahler_1d(B0b_LANE):.12f}   "
      f"log(4/9) = {math.log(4/9):.12f}")
print(f"      stratum     : S4 {strat_sorted(B0b_S4)}   lane Z {strat_sorted(B0b_LANE)}")
print("""
    SO THE FINDING SURVIVES AND THE CLAIM DOES NOT.  Z-6's existential witness is still a
    witness (the multiset is the same, so max p_a = 4/9 <= 1/2 and the zero set is still
    empty), and Z-14's stratum for B0b is still EMPTY.  What falls is Z-14's own header --
    "THE GEOMETRIC TYPE OF EVERY CARRIER ROW OF S4:575" -- because one of the ten rows is not
    the row S4 wrote; and what falls with it is the implicit claim that a Z_k SEQUENCE
    computed from that entry (z4(b) and the EMPTY row of z4(c)) is B0b's.  IT IS NOT: the two
    states have identical Jensen branches but DIFFERENT Z_k, because Z_k depends on (u^k,v^k)
    separately and not only on the pair-symmetric moduli.  Lane Z's only EMPTY-stratum
    convergence exhibit therefore runs on a polynomial that splits into two independent
    one-variable factors -- the easiest possible case, and not the corpus's own.
""")
rng = np.random.default_rng(20260817)
mx = 0.0
for _ in range(20000):
    s, tt = rng.random(2) * 2 * np.pi
    v1 = float(Pabs(B0b_S4, s, tt)); v2 = float(Pabs(B0b_LANE, s, tt))
    mx = max(mx, abs(v1 - v2))
print(f"    max | |P|_S4(x,y) - |P|_laneZ(x,y) | over 20000 random torus points = {mx:.6f}")
print("      (0 would mean the two states are the same function; they are not)\n")

# =========================================================================================
print("-" * W)
print("(b) DEFECT 2 -- THEOREM Z4 AND THE CLOSED FORM cos s0 = -2/3 REST ON A LABELLING")
print("    S4 DOES NOT FIX.  S4:566 defines SENSE C as \"(0.4,0.3,0.3) for 3 classes\" and")
print("    says nothing about WHICH class carries the 0.4.  There are three assignments.")
print("    The MULTISET quantities are the same in all three, exactly as W-03/N2 requires.")
print("    THE ZERO LOCATION IS NOT, AND THEOREM Z4's HYPOTHESIS p10 = p01 HOLDS IN ONLY TWO")
print("    OF THE THREE.\n")
print(f"    {'assignment (p00,p10,p01,p11)':<36s} {'stratum':<7s} {'cos s0':>9s} {'m(P)':>14s} "
      f"{'x0*y0':>26s} {'ANTI-DIAG?':>11s}")
for lbl, p in (("0.4 on p11  <- lane Z / the brief", fr(F(0), F(3,10), F(3,10), F(2,5))),
               ("0.4 on p01", fr(F(0), F(3,10), F(2,5), F(3,10))),
               ("0.4 on p10", fr(F(0), F(2,5), F(3,10), F(3,10)))):
    st = strat_sorted(p)
    ang = zero_angles(p); pts = zero_points(p)
    p00, p10, p01, p11 = p
    C = p00**2 + p10**2 - p01**2 - p11**2
    D = 2 * (p00 * p10 - p01 * p11)
    cs = F(-C, 1) / D
    x0, y0 = pts[0]
    prod = x0 * y0
    print(f"    {lbl:<36s} {st:<7s} {str(cs):>9s} {mahler_1d(p):14.12f} "
          f"{f'{prod.real:+.9f}{prod.imag:+.9f}i':>26s} {str(abs(prod-1) < 1e-12):>11s}")
print("""
    ==> stratum and m(P) are identical in all three (multiset invariance, as N2 says), so
    lane Z's Z-1/Z-4/Z-8 verdicts are labelling-independent and stand.  BUT Z-2's headline
    closed form and Z-10 (THEOREM Z4, "no connection satisfying H2 ever lands exactly on a
    zero") are stated as facts about "K1's REGISTERED pi", and they are facts about ONE OF
    THREE READINGS OF S4's SENSE C.  Under the third reading x0 y0 != 1, the two zeros are
    NOT on the anti-diagonal subtorus, and the exact-hit obstruction Theorem Z4 asserts does
    not exist.  THE OPERATIVE VARIABLE FOR Z-10 IS p10 = p01, NOT "K1's REGISTERED pi", AND
    LANE Z NAMES IT CORRECTLY IN THE THEOREM'S BODY ("indeed whenever p10 = p01") AND
    INCORRECTLY IN ITS HEADLINE AND IN FINDING Z-10's CLAIM FIELD.
""")

# =========================================================================================
print("-" * W)
print("(c) DEFECT 3 -- THE HEADLINE'S SIGN.  Lane Z's headline calls the EMPTY 3/4 \"THE")
print("    DECISIVE STRUCTURAL FACT\" because there H2 settles N1 and m(P) is an elementary")
print("    logarithm.  BOTH HALVES ARE TRUE AND TOGETHER THEY SAY THE OPPOSITE.\n")
print("""    On the EMPTY stratum m(P) = log(p_max): N1's content there is 'the rate is the log
    of the largest class weight', and NO MAHLER MEASURE IS INVOLVED.  The 2-variable Mahler
    measure -- the object N1 is FOR, the thing that "inherits the entropy theory of algebraic
    Z^d-actions wholesale" -- is only non-elementary on the crossing quarter.  So:

        THE 3/4 WHERE LANE Z SETTLES THE CONVERGENCE IS EXACTLY THE 3/4 WHERE N1 HAS NO
        MAHLER-MEASURE CONTENT, AND THE 1/4 WHERE N1 HAS CONTENT IS EXACTLY THE 1/4 LANE Z
        LEAVES OPEN.  The two facts are the same fact with opposite signs.

    AND THIS IS NOT NEW.  LANE_W08_M3_ZEROSET/m3_2_fourclass.OUT.txt (d), sealed 2026-08-16,
    the lane Z itself convicts the register of not reading, says it in one line:

      "(so the Mahler measure is a LOGARITHM OF A WEIGHT off the firing region, and only
       inside it does the dilogarithm/Cassaigne-Maillot regime begin.  The firing region is
       exactly where the analytically interesting rate lives.)"

    LANE Z CITES THAT EXACT BLOCK (Z-7's evidence field cites "M3-2(d)") AND INVERTS ITS SIGN.
    Z-13 is lane Z's own finding that this program fails by UNDER-READING.  Z-7 and the
    headline are that failure, one directory later, in the lane that found it.

    AND THE IMPORT IS MISSING TOO.  The three-class half of lane Z's THEOREM Z2 -- 'no
    crossing => m = log max' -- is the elementary branch of CASSAIGNE-MAILLOT (1997), which
    this corpus cites EIGHTEEN TIMES (REGISTER_V001.md:313 counts them) and which S4:590-597
    uses for five of its own nine carrier rates.  `grep -ril 'cassaigne\\|maillot\\|dilog'`
    over LANE_W13_Z_ZEROSET returns NOTHING.  Z-7 presents the statement as one 'nobody
    generalised'; what nobody did was cite the 1997 theorem that is half of it.
""")
print("    CHECKED: is m(P) the log of a weight on the crossing quarter?  If it ever were,")
print("    the inversion above would be wrong.  200 random TWO-stratum states:")
rng = np.random.default_rng(20260818)
worst_gap = 9.9; ncheck = 0
for _ in range(4000):
    w = rng.dirichlet(np.ones(4))
    p = tuple(F(float(v)).limit_denominator(10**6) for v in w)
    ssum = sum(p)
    p = tuple(q / ssum for q in p)
    if strat_sorted(p) != 'TWO':
        continue
    ncheck += 1
    m = mahler_1d(p, 1 << 18)
    gap = min(abs(m - math.log(float(q))) for q in p if q > 0)
    worst_gap = min(worst_gap, gap)
    if ncheck >= 200:
        break
print(f"      {ncheck} TWO-stratum states: SMALLEST |m(P) - log(any weight)| = {worst_gap:.6e}")
print(f"      (on the EMPTY stratum the same quantity is 0 to 2.2e-16 -- M3-2(d), 3029 states)")
print("      ==> on the crossing quarter m(P) is NEVER the log of a class weight.\n")

# =========================================================================================
print("-" * W)
print("(d) AND THE MEASURE THE 3/4 IS TAKEN IN.  Lane Z's 3/4 is LEBESGUE measure on the")
print("    simplex.  The corpus has no measure on ready states -- W-14 says every ready state")
print("    in the corpus was CHOSEN, with ONE exception, SENSE U, where pi is fixed by the")
print("    carrier's own class sizes and is 'the one state rule that is not a stipulation'.")
print("    SO THE ONLY NON-STIPULATED SAMPLE THE CORPUS HAS IS S4:575's CARRIER COLUMN.")
print("    HERE IS ITS STRATUM CENSUS -- computed from S4's rows AS WRITTEN, not from lane Z's")
print("    transcription.\n")
print(f"    {'carrier row (S4:575, verbatim)':<48s} {'pi (p00,p10,p01,p11)':<30s} {'stratum':<9s}")
cens = {}
for name, p in S4_575:
    st = strat_sorted(p)
    cens[st] = cens.get(st, 0) + 1
    print(f"    {name:<48s} {str(tuple(str(q) for q in p)):<30s} {st:<9s}")
extra = [("SENSE C, three classes  = N1 AS REGISTERED", K1_REG),
         ("SENSE C, four classes   (S4:597)", SENSEC4),
         ("S1 sec6 published ready state (M1_08 T1)", S1_PUB)]
for name, p in extra:
    st = strat_sorted(p)
    cens[st] = cens.get(st, 0) + 1
    print(f"    {name:<48s} {str(tuple(str(q) for q in p)):<30s} {st:<9s}")
n = sum(cens.values())
nsing = n - cens.get('EMPTY', 0)
distinct = {}
for _n, _p in S4_575 + extra:
    distinct[_p] = strat_sorted(_p)
dn = len(distinct)
dsing = sum(1 for v in distinct.values() if v != 'EMPTY')
print(f"\n    CENSUS over the corpus's own {n} published (carrier, state) rows: {cens}")
print(f"    SINGULAR (zero set non-empty): {nsing} of {n} = {nsing/n:.3f}")
print(f"    DEDUPED (B1 = B2 = B3 are one state): {dsing} of {dn} DISTINCT states = "
      f"{dsing/dn:.3f}   -- reported both ways so the count cannot be read as inflated")
print(f"    LANE Z's LEBESGUE FIGURE:      EMPTY on 3/4 = 0.750 of the simplex")
print("""
    ==> THE CORPUS'S OWN SAMPLE POINTS THE OTHER WAY, BY A FACTOR OF THREE.  Lane Z's Z-9
    ("on 3/4 of all FORMING states there is no singularity at all... the corpus does not live
    in the hard regime by necessity; it lives there at its two registered states") is true of
    Lebesgue measure and false of every state the corpus has ever computed with but two.
    THIS IS NOT A CORRECTION TO THE ARITHMETIC -- the 1/4 and 3/4 are exact and this refuter
    re-derived them -- IT IS A CORRECTION TO WHAT THEY WEIGH.  A measure nobody adopted is
    doing the work in the headline.
    AND A SECOND, SMALLER SCOPE SLIP IN THE SAME FINDING: Z-9 says m(P) = 0 "except at the
    three corners max p_a = 1".  On the FOUR-class simplex, where the 1/4 and 3/4 live, there
    are FOUR such corners.  Three is M1_08 T4's count and M1_08 T4 is K1-scoped (p00 = 0).
    Z-9 mixes a three-class corner count with a four-class volume.
""")
print("DONE r3")
