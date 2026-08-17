# W-08 / M4-REFUTE-2 leg B — LENS 2 (MATHEMATICS).
# DOES THE DECISIVE TEST DISCRIMINATE?  M4-1 calls m4_c C1 a refutation of W-07's name "at BOTH
# edges".  A test refutes a NAME only on rows where the two rival names PREDICT DIFFERENT THINGS.
# This leg computes, row by row, what each rival name predicts, and how far the counterexample
# reaches in K.
#
# THE TWO RIVAL NAMES, stated so they can be evaluated:
#   W-07 (sec 3, its own sharp form):  cells = floor(K / ord(rho));  "1000 = 4000/4 is ord(rho)".
#   M4   (operative_variable):         cells = #{k<=K : ||k theta|| < arcsin(tol/(2 amp))/pi}.
#
# ISOLATION LEDGER.
#  B1 HELD FIXED: K=4000, tol=1e-9, amp = |s_2||s_3| for W-07 leg E's seed-20260816 state,
#     observable A_23, exact rational arithmetic.  MOVED: ord(rho) = q alone, q = 1..8000.
#  B2 HELD FIXED: everything in B1.  MOVED: nothing -- B2 is a proof-by-exhaustion over B1's output.
#  B3 HELD FIXED: the C4 statement.  MOVED: the tolerance d alone, into the range the statement
#     quantifies over ("any d > 0") but its proof does not.
#  B4 HELD FIXED: S1's published connection (ord(rho)=4), observable A_23, K, tol.
#     MOVED: the ready state alone -- generic (seed 20260816) versus S1's PUBLISHED p.
#  B5 HELD FIXED: the counterexample family theta = 1/4 + lam, tol = 1e-9, amp, observable.
#     MOVED: K alone, over W-07's own published grid K = 4e3, 1e4, 1e5, 1e6, 1e7.
from fractions import Fraction
import math, numpy as np

rng = np.random.default_rng(20260816)
s = rng.normal(size=5) + 1j*rng.normal(size=5); s /= np.linalg.norm(s)
amp = float(abs(s[2])*abs(s[3])); K = 4000; TOL = 1e-9
eps = math.asin(TOL/(2*amp))/math.pi          # ||k theta|| < eps  <=>  D_k < TOL
print(f"  amp = {amp!r}   eps = arcsin(tol/2amp)/pi = {eps:.6e}\n")

print("== B1  THE RATIONAL SIDE, EXACTLY.  q = ord(rho) from 1 to 8000, theta = 1/q. ==")
print("  For theta = p/q in lowest terms, min_{k<=K}||k theta|| = 0 if q<=K, and >= 1/q if q>K.")
def cells_rational(q):
    # exact: ||k/q|| = min(k mod q, q - k mod q)/q ; count those < eps
    c = 0
    for k in range(1, K+1):
        r = k % q; d = min(r, q-r)
        if Fraction(d, q) < Fraction(eps).limit_denominator(10**15): c += 1
    return c
bad = []
for q in list(range(1, 60)) + [100, 999, 1000, 1999, 2000, 2001, 3999, 4000, 4001, 5000, 8000]:
    c = cells_rational(q); pred_w07 = K//q
    if c != pred_w07: bad.append((q, c, pred_w07))
print(f"  q values checked: 1..59 plus 100,999,1000,1999,2000,2001,3999,4000,4001,5000,8000")
print(f"  rows where observed cell count != W-07's floor(K/ord) : {len(bad)}   {bad}")
print("  W-07 sec3's OWN formula ('1000 = 4000/4 is ord(rho)') therefore PREDICTS, exactly:")
for q in (4, 2000, 4001, 8000):
    print(f"      ord(rho) = {q:<5} ->  floor(4000/{q}) = {K//q:<5}  observed {cells_rational(q)}")
print("  m4_c C1's four rational rows return 1000 / 0 / 0 / 2.  ALL FOUR ARE W-07's OWN FORMULA.")
print("  On rationals of SMALL denominator the two names agree: ord(rho)<=K <=> min||k theta||=0.")
print("  (They part company on rationals of denominator > 1/eps = 1.7e9; see leg D block D4, which")
print("   exhibits a FINITE-order counterexample and needs no irrationality at all.)")
print("  So rows 4001 and 8000 -- M4's entire 'FINITE edge' -- discriminate NOTHING between the two")
print("  names.  They refute the four-word paraphrase 'finite versus infinite', not the measurement.\n")

print("== B2  WHICH ROWS OF THE DECISIVE TEST ACTUALLY DISCRIMINATE? ==")
LAM13 = (2**0.5-1)*1e-13; LAM16 = (2**0.5-1)*1e-16
def cells_irr(th):
    c = 0
    for k in range(1, K+1):
        f = (k*th) % 1.0; d = min(f, 1-f)
        if d < eps: c += 1
    return c
print(f"  {'row':<32} {'W-07 name predicts':>20} {'M4 name predicts':>18} {'discriminates?':>15}")
for tag, w07, m4 in [
    ("ord=4    (S1 published)", "1000 (=K/4)", "1000"),
    ("ord=4001                ", "0 (=floor(K/4001))", "0"),
    ("ord=8000                ", "0 (=floor(K/8000))", "0"),
    ("ord=2000                ", "2 (=floor(K/2000))", "2"),
    ("irrational 1/4+lam e-13 ", "0 (ord infinite)", f"{cells_irr(0.25+LAM13)}"),
    ("irrational 1/4+lam e-16 ", "0 (ord infinite)", f"{cells_irr(0.25+LAM16)}")]:
    disc = "YES" if w07.split()[0] != m4 else "no"
    print(f"  {tag:<32} {w07:>20} {m4:>18} {disc:>15}")
print("  TWO of six rows carry any discriminating power, and leg A showed one of those two (e-16)")
print("  is realised at 46% error by its own float constructor.  The decisive test is ONE ROW.\n")

print("== B3  C4's QUANTIFIER IS WRONG.  'for any K, any n <= K and any d > 0' ==")
print("  C4 claims: there is an irrational theta with #{k<=K : |exp(2pi i k theta)-1| < d} = floor(K/n).")
print("  Its proof needs the non-qualifying terms bounded below by 2 sin(pi/n) AND d <= that bound.")
print("  Counterexamples inside C4's own quantifier, at theta = 1/n itself and hence on a whole")
print("  neighbourhood of irrationals around it:")
for n, d, KK in [(2, 2.5, 10), (4, 1.5, 4000), (3, 1.8, 100), (4, 1.4143, 4000)]:
    cnt = sum(1 for k in range(1, KK+1) if abs(complex(math.cos(2*math.pi*k/n), math.sin(2*math.pi*k/n)) - 1) < d)
    print(f"    n={n:<3} d={d:<7} K={KK:<6}  2 sin(pi/n) = {2*math.sin(math.pi/n):.6f}   count = {cnt:<6} floor(K/n) = {KK//n}")
print("  The statement is TRUE for d <= 2 sin(pi/n) and FALSE above it.  As written it is false;")
print("  the correction is one clause and it does not touch M4-1's application (d/amp = 3.68e-9).\n")

print("== B4  M4's 'ONLY SURVIVING FORM' IS FALSE AT S1's OWN PUBLISHED READY STATE ==")
print("  M4: 'D_k = 0 EXACTLY for some k <= K  iff  ord(rho) is finite AND ord(rho) <= K.'")
print("  D_k = amp |rho^k - 1| with amp = |s_u| |s_v|.  The 'iff' silently assumes amp != 0.")
p_pub = [Fraction(1,2), 0, 0, Fraction(1,4), Fraction(1,4)]     # S1 sec6 / W-01 / S3, published
print(f"    S1's PUBLISHED ready state p = (1/2, 0, 0, 1/4, 1/4)  ->  |s_1| = |s_2| = 0 EXACTLY.")
print(f"    Observable A_23 = conj(t_2) t_3 :  amp = |s_2||s_3| = 0 EXACTLY.")
print(f"    Hence D_k = 0 for EVERY k and EVERY rho, including every irrational one.")
print(f"    The right-to-left direction of M4's iff therefore fails on the corpus's own published")
print(f"    state -- the same state M4's own leg G4 computes with.  On K1 at that state the ordered")
print(f"    pairs with amp = 0 are those meeting {{v1,v2}}:")
zero_pairs = [(u,v) for u in range(5) for v in range(5) if u!=v and (p_pub[u]==0 or p_pub[v]==0)]
print(f"      pairs with amp = 0 : {len(zero_pairs)} of 20  -> {zero_pairs}")
print("  Correct form:  amp != 0  =>  ( D_k = 0 for some k <= K  <=>  ord(rho) finite and <= K ).")
print("  M4 also drops the companion hypothesis tol < 2 amp; m4_g's own code carries min(1.0, tol/(2amp)),")
print("  so the lane knew, and the finding text does not say it.\n")

print("== B5  HOW FAR IN K DOES THE COUNTEREXAMPLE REACH?  W-07's OWN GRID. ==")
print("  W-07 sec3 publishes a K-scaling table at K = 1e3,1e4,1e5,1e6,1e7 and states the claim")
print("  there: 'exact zero on exactly K/4 cells, at EVERY K, FOREVER' vs 'never zero at any k'.")
print("  M4's counterexample is built at K=4000 with lam = (sqrt2-1)e-13.  Run the SAME lam over")
print("  W-07's grid, changing K ALONE:")
lam = (2**0.5-1)*1e-13
print(f"  {'K':>10} {'cells<1e-9 (irrational)':>24} {'K/4 (ord=4)':>14} {'ratio to K/4':>14}")
for KK in [4000, 10**4, 10**5, 10**6, 10**7]:
    # ||k theta|| = ||k/4 + k lam||; qualifying k are the multiples of 4 with 4m*lam < eps
    mmax = int(eps/(4*lam))
    c = min(KK//4, mmax)
    print(f"  {KK:>10} {c:>24} {KK//4:>14} {c/(KK//4):>14.4f}")
print(f"  The count SATURATES at {int(eps/(4*lam))} and the density collapses like 1/K.")
print("  At K = 1e6 the irrational row gives 0.35% of the cells the ord=4 row gives; at 1e7, 0.035%.")
print("  To hold '1000 of 4000'-style agreement out to K = 1e7 the counterexample needs")
print(f"    lam < eps/K = {eps/1e7:.3e}, i.e. below the ~2e-17 float noise floor leg A measured in")
print("  m4_c's own constructor.  The counterexample is a FIXED-(K,tol) object.  M4-1 says it holds")
print("  'on W-07's own ... k-range'; W-07's k-range is 1e3..1e7, not 4000.")
