# W-08 / M4-REFUTE-2 leg D — LENS 2 (MATHEMATICS).  THE NAMED SCALAR DOES NOT DETERMINE THE
# REPORTED EFFECT.
#
# M4's operative_variable field, verbatim: "The effect W-07 measures is a function of ONE
# quantity: min_{1<=k<=K} ||k*theta|| ... compared against arcsin(tol/(2*amp))/pi."
#
# But the effect M4 REPORTS in every row of its decisive test, and the effect W-07 registered,
# is a CELL COUNT (#{k<=K : D_k < tol}), not the min.  A count is not a function of a min.
# This leg exhibits two irrational theta with the SAME min_{k<=K}||k theta|| -- to the last bit --
# and cell counts 1000 and 4.  Both are verified in EXACT arithmetic (Fraction + a 60-digit pi),
# not float, because the perturbations are ~1e-16 and float64 cannot carry them (leg A).
#
# ISOLATION LEDGER.  HELD FIXED: K = 4000, tol = 1e-9, amp = |s_2||s_3| at W-07 leg E's own
# seed-20260816 state, observable A_23, the counting function, and the value of
# min_{k<=K}||k theta|| itself -- held equal BY CONSTRUCTION to 22 significant digits.
# MOVED: the DENOMINATOR q of the nearby rational alone (q = 4 versus q = 1000).
import numpy as np, math
from fractions import Fraction
from decimal import Decimal, getcontext, ROUND_FLOOR
getcontext().prec = 60
PI = Decimal("3.14159265358979323846264338327950288419716939937510582097494")
def frac1(x): return x - x.to_integral_value(rounding=ROUND_FLOOR)
def nrm(x):
    f = frac1(x); return f if f <= Decimal("0.5") else 1-f
def dsin(x):
    t=x; ssum=x; n=1
    while abs(t) > Decimal("1e-58"):
        t = -t*x*x/Decimal((2*n)*(2*n+1)); ssum += t; n += 1
    return ssum

rng = np.random.default_rng(20260816)
s = rng.normal(size=5)+1j*rng.normal(size=5); s/=np.linalg.norm(s)
amp = Decimal(repr(float(abs(s[2])*abs(s[3])))); K = 4000; TOL = Decimal("1e-9")
EPS = Decimal(repr(math.asin(float(TOL)/(2*float(amp)))/math.pi))
print(f"  K = {K}   tol = 1e-9   amp = {amp}   eps = {EPS:.6e}")
print("  (D_k < tol  <=>  ||k theta|| < eps, exactly; see m4_c C3 and leg A.)\n")

SQ2 = Decimal(2).sqrt()
lam = (SQ2-1)*Decimal("1e-13")          # irrational
mu  = lam/250                           # irrational; chosen so 1000*mu == 4*lam exactly
thA = Decimal("0.25")   + lam           # nearby rational 1/4,    denominator q = 4
thB = Decimal("0.001")  + mu            # nearby rational 1/1000, denominator q = 1000

def profile(th):
    mn=None; cnt=0; kmn=None
    for k in range(1, K+1):
        d = nrm(Decimal(k)*th)
        if mn is None or d < mn: mn, kmn = d, k
        if d < EPS: cnt += 1
    return mn, cnt, kmn

print("== D1  TWO IRRATIONAL theta WITH IDENTICAL min_{k<=4000} ||k theta|| ==")
mA, cA, kA = profile(thA); mB, cB, kB = profile(thB)
print(f"  theta_A = 1/4    + (sqrt2-1)e-13      irrational, ord(rho) = INFINITY")
print(f"  theta_B = 1/1000 + (sqrt2-1)e-13/250  irrational, ord(rho) = INFINITY")
print()
print(f"  {'':<10} {'min_(k<=K) ||k theta||':>30} {'at k':>7} {'cells D_k < 1e-9':>18}")
print(f"  {'theta_A':<10} {mA:>30.24f} {kA:>7} {cA:>18}")
print(f"  {'theta_B':<10} {mB:>30.24f} {kB:>7} {cB:>18}")
print(f"  min_A - min_B = {mA-mB:.1e}   (EXACTLY zero: 4*lam = 1000*mu by construction)")
print(f"  min D_A = 2 amp sin(pi min_A) = {2*amp*dsin(PI*mA):.6e}")
print(f"  min D_B = 2 amp sin(pi min_B) = {2*amp*dsin(PI*mB):.6e}   -- identical")
print()
print("  SAME NAMED SCALAR.  COUNTS 1000 AND 4.  M4's 'the effect is a function of ONE quantity,")
print("  min_{k<=K}||k theta||' is FALSE for the effect its own decisive test reports.\n")

print("== D2  WHAT DOES DETERMINE THE COUNT?  THE DENOMINATOR OF THE NEARBY RATIONAL. ==")
print("  For theta = p/q + delta with q*K*|delta| < eps and gcd(p,q)=1, the qualifying k are")
print("  exactly the multiples of q, so   count = floor(K/q)   -- and q is the ORDER of the")
print("  root of unity exp(2 pi i p/q).  Checked over a family, one variable moved (q):")
print(f"  {'q':>7} {'theta = 1/q + lam/(q/4)':>28} {'count':>8} {'floor(K/q)':>12} {'min ||k theta||':>20}")
for q in (2, 4, 5, 8, 10, 100, 1000, 2000):
    d = lam*4/q
    th = Decimal(1)/Decimal(q) + d
    m, c, kk = profile(th)
    print(f"  {q:>7} {'1/'+str(q)+' + delta':>28} {c:>8} {K//q:>12} {m:>20.18f}")
print("  count = floor(K/q) in every row while min ||k theta|| is HELD FIXED at 1.657e-13.")
print("  THE COUNT IS A FUNCTION OF q AND THE MIN IS A FUNCTION OF delta.  They are independent")
print("  coordinates of the same point, and M4 named only the second while reporting only the first.\n")

print("== D3  WHAT THIS DOES TO THE TWO NAMES ==")
print("  W-07 registered:  count = K / ord(rho).            Wrong only in that ord(rho) must be")
print("                    read as the order q of the NEAREST root of unity within tolerance,")
print("                    and the count is floor(K/q), zero when q > K.")
print("  M4 registered:    the effect is a function of min_{k<=K}||k theta|| alone.  Refuted by D1.")
print("  CORRECT NAME, and it contains both:  the effect at (K, tol) is a function of the PAIR")
print("      ( q , delta )  =  ( denominator of the best rational approximation to theta with")
print("                          q <= K , the residual |theta - p/q| ),")
print("  with, whenever q K |delta| stays below 1/2,")
print("      count = min( floor(K/q) , ceil(eps/(q|delta|)) - 1 )      and      min = q|delta|.")
print("  W-07's name is the q coordinate; M4's is the delta coordinate.  W-07 named the one that")
print("  carries its own headline number 1000; M4 named the one that carries the min, then")
print("  convicted W-07 of misnaming.  Formula checked against brute force on the D2 family:")
def count_formula(q, delta):
    import math as _m
    return min(K//q, _m.ceil(float(EPS)/(q*float(delta))) - 1)
ok = True
for q in (2,4,5,8,10,100,1000,2000):
    d = lam*4/q; th = Decimal(1)/Decimal(q) + d
    _, c, _ = profile(th)
    f = count_formula(q, d)
    if c != f: ok = False
    print(f"      q={q:<6} brute force {c:<6} formula {f:<6} {'OK' if c==f else '***MISMATCH***'}")
print(f"  formula reproduces brute force on all rows: {ok}\n")

print("== D4  AND THE COUNTEREXAMPLE DID NOT NEED TO BE IRRATIONAL AT ALL ==")
print("  M4 built its infinite-order side out of an irrational theta and defended it with a")
print("  closed-form argument about irrationality (m4_g G3b).  None of that is needed: a RATIONAL")
print("  theta of large finite order does it, entirely inside exact arithmetic.")
th_rat = Fraction(1,4) + Fraction(1, 10**13)
print(f"    theta = 1/4 + 1e-13 = {th_rat}   ->  ord(rho) = {th_rat.denominator}  (FINITE, and > K)")
cnt = 0; mn = None
for k in range(1, K+1):
    r = (k*th_rat) % 1
    d = min(r, 1-r)
    if mn is None or d < mn: mn = d
    if d < Fraction(EPS): cnt += 1
print(f"    min_(k<=4000) ||k theta|| = {float(mn):.6e}    cells D_k < 1e-9 : {cnt}")
print(f"    W-07's own formula floor(K/ord) predicts {K//th_rat.denominator}.  Observed {cnt}.")
print("    So W-07's 'count = K/ord(rho)' is refuted by a FINITE-ORDER ratio, with no appeal to")
print("    irrationality and no float anywhere -- a cleaner refutation than M4's, and one that")
print("    also refutes M4's own framing of the question as 'finite versus infinite order'.")
print("    It is the (q, delta) pair that governs, and NEITHER lane named it.")
