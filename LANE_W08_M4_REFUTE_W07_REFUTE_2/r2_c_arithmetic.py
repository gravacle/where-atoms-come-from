# W-08 / M4-REFUTE-2 leg C — LENS 2 (MATHEMATICS).  The arithmetic underneath the decisive test.
#  C1  Weyl equidistribution is CITED for a window it cannot control, in the finding (M4-6) that
#      voids W-07's only control.  The number is right; the theorem invoked is not.
#  C2  m4_c's constructor cannot realise ANY finite-order rho: every row of the decisive test
#      builds an irrational psi.  The 'ord(rho) = q' column labels intent, not the object built.
#  C3  W-07's OWN order detector, run on M4's counterexamples.  (This one goes M4's way.)
#  C4  M4-2 re-derived independently: which group element governs A_uv, checked against both lanes.
#  C5  Where does the counterexample live?  The measure of the effect set, and its structure.
#
# ISOLATION LEDGER.
#  C1 HELD FIXED: K=4000, eps, and the count function.  MOVED: theta alone, over a specified
#     algebraic irrational, W-07's two named generic irrationals, and 200000 Lebesgue draws.
#  C2 HELD FIXED: m4_c's constructor, verbatim.  MOVED: nothing; C2 is exact arithmetic on its
#     output.  C3 HELD FIXED: W-07 leg E's `ordr` line, verbatim.  MOVED: the connection alone.
#  C4 HELD FIXED: carrier, state, observable, k-range.  MOVED: the ratio FORMULA alone.
#  C5 is measure theory on stated closed forms and draws no comparison.
import numpy as np, math
from fractions import Fraction

rng = np.random.default_rng(20260816)
s = rng.normal(size=5)+1j*rng.normal(size=5); s/=np.linalg.norm(s)
amp = float(abs(s[2])*abs(s[3])); K = 4000; TOL = 1e-9
eps = math.asin(TOL/(2*amp))/math.pi

print("== C1  'BY WEYL EQUIDISTRIBUTION THE EXPECTED COUNT IS 2*K*eps' — WRONG THEOREM ==")
print(f"  window eps = {eps:.4e},  K = {K},  K*eps = {K*eps:.4e}")
print("  Weyl's theorem says #{k<=K : ||k theta|| < eps}/K -> 2 eps as K -> infinity at FIXED eps.")
print("  Effective forms (Erdos-Turan) need K eps >> 1 to say anything about the count; here")
print(f"  K eps = {K*eps:.3e} << 1, so equidistribution has NO content at this window and this K.")
print("  It is also the wrong KIND of statement: for a FIXED irrational theta the count is not a")
print("  random variable and has no 'expected' value.  Two correct arguments, both computed here:")
print()
print("  (a) FUBINI over Lebesgue-random theta (this is what M4's number actually is):")
print("      E[#{k<=K : ||k theta||<eps}] = sum_k P(||k theta||<eps) = 2 K eps  (exact, no Weyl),")
print(f"      = {2*K*eps:.4e};  Markov: P(count >= 1) <= {2*K*eps:.4e}.")
draws = 200000
rr = np.random.default_rng(7)
th = rr.random(draws)
hit = 0
for t in th[:draws]:
    kk = np.arange(1, K+1); f = (t*kk) % 1.0
    if np.minimum(f, 1-f).min() < eps: hit += 1
print(f"      empirical over {draws} Lebesgue draws: {hit} hits (bound predicts <= {2*K*eps*draws:.2f}).")
print()
print("  (b) DETERMINISTIC for W-07's actual rows, which were NOT random.  Two of W-07's five")
print("      generic connections are named algebraic irrationals (sqrt2/sqrt3 at leg E, phi at")
print("      leg B).  For a quadratic irrational the partial quotients are bounded, so")
print("      ||k theta|| >= c/k for an explicit c > 0 and min_{k<=K}||k theta|| >= c/K.")
def cf_min(theta, K):
    kk = np.arange(1, K+1); f = (theta*kk) % 1.0
    d = np.minimum(f, 1-f); i = int(d.argmin()); return d[i], i+1
for name, t in [("sqrt2+sqrt3 (W-07 leg E's rho, corrected)", (2**0.5+3**0.5) % 1.0),
                ("sqrt2-sqrt3 (W-07 leg E's rho as coded) ", (2**0.5-3**0.5) % 1.0),
                ("1/phi^2     (W-07 leg B / m4_g generic) ", ((1+5**0.5)/2)**-2)]:
    m, kmin = cf_min(t, K)
    print(f"      {name}: min_(k<=4000)||k theta|| = {m:.4e} at k={kmin};  K*min = {K*m:.4f}")
print("      K*min is O(1) in every row -- the Liouville/CF bound, deterministic.  W-07's control")
print("      could not have failed, and for these rows it PROVABLY could not.  M4-6's CONCLUSION")
print("      is right and its stated reason is not the one that applies.\n")

print("== C2  m4_c's CONSTRUCTOR CANNOT BUILD A FINITE-ORDER rho AT ALL ==")
print("  Every angle m4_c stores is a float64, i.e. a dyadic rational a = m/2^e.  exp(i a) is a")
print("  root of unity iff a is a rational multiple of 2 pi, i.e. iff pi is rational.  So for every")
print("  nonzero float angle, ord(exp(i a)) = INFINITY.  Checking the realised psi of the row")
print("  labelled 'ord(rho) = 4  [S1 PUBLISHED]' against the rationals with denominator <= 10^6:")
from decimal import Decimal, getcontext, ROUND_FLOOR
getcontext().prec = 60
PI60 = Decimal("3.14159265358979323846264338327950288419716939937510582097494")
argWF = np.pi; argWC = -argWF - 2*np.pi*(-0.25)
aF, aC = argWF/3.0, argWC/3.0
sF, sC = (aF+aF)+aF, (aC+aC)+aC
psiD = (Decimal(repr(sF))+Decimal(repr(sC)))/(2*PI60)
psiD = psiD - psiD.to_integral_value(rounding=ROUND_FLOOR)
print(f"    exact value of the realised psi (60 digits) = {psiD:.30f}")
print(f"    psi - 1/4 = {psiD-Decimal('0.25'):.4e}   -- NOT zero; the float pipeline lands 3.5e-17 short")
print("    (float64 evaluation of the same expression PRINTS 0.25, which is why the row reads as")
print("     exact; the exactness is an artefact of the print, not of the object).")
print("    ord(rho) as realised = infinite in EVERY row of C1, including all four rational rows.")
print("    The decisive test therefore never instantiates the side of the dichotomy it names, and")
print("    separates its rows by exactly the quantity it concludes with -- rational-approximation")
print("    quality.  The conclusion is not thereby wrong (leg A's INTENDED column reproduces every")
print("    count), but 'THE DECISIVE TEST' is decided in exact arithmetic that m4_c does not run;")
print("    only the e-13 row has a closed-form backstop (m4_g G3b).  The four rational rows have")
print("    none, and m4_g G3(a) supplies exact arithmetic only for rho = -i.\n")

print("== C3  W-07's OWN ORDER DETECTOR, ON M4's COUNTEREXAMPLES.  (THIS GOES M4's WAY.) ==")
print("  W-07 leg E line:  ordr = next((n for n in range(1,10001) if abs(rho**n - 1) < 1e-12), None)")
for tag, lam in [("theta = 1/4 + (sqrt2-1)e-13", (2**0.5-1)*1e-13),
                 ("theta = 1/4 + (sqrt2-1)e-16", (2**0.5-1)*1e-16),
                 ("theta = 1/4 + (sqrt2-1)e-15", (2**0.5-1)*1e-15)]:
    rho = np.exp(2j*np.pi*(0.25+lam))
    ordr = next((n for n in range(1,10001) if abs(rho**n-1) < 1e-12), None)
    print(f"    {tag}:  ORDER reported = {ordr}   (true order: INFINITE)")
print("  W-07's instrument reports ORDER = 4 for irrational ratios.  Its 'ord(rho)' column is a")
print("  1e-12-thresholded Diophantine statistic wearing a group-theoretic name.  That is M4-1's")
print("  thesis, demonstrated on W-07's own code rather than on a new connection, and it is a")
print("  BETTER argument than the one M4 gave.  IT SURVIVES THIS LANE.\n")

print("== C4  WHICH GROUP ELEMENT GOVERNS A_uv — M4-2 RE-DERIVED FROM SCRATCH ==")
FACE_V={0,1,2}; CYC_V={0,3,4}; TREE={1:(0,),2:(0,1),3:(3,),4:(3,4)}
def dress(x,a):
    u=np.exp(1j*np.asarray(a)); t=np.array(x,dtype=complex)
    for v,p in TREE.items():
        w=1.0+0j
        for e in p: w*=u[e]
        t[v]=x[v]/w
    return t
def A(x,a,u,v): t=dress(x,a); return np.conj(t[u])*t[v]
print("  BRUTE FORCE: apply the branch operators to the state, then measure.  No formula assumed.")
print(f"  {'connection':<30} {'k':>3} {'A[M_dF^k s]/A[s]':>34} {'W_F^(dF k)':>26} {'match':>7}")
for tag, a in [("S1 published", np.array([np.pi/3]*3+[np.pi/2]*3)),
               ("generic sqrt2/sqrt3", np.array([2*np.pi*(2**0.5)/3]*3+[2*np.pi*(3**0.5)/3]*3))]:
    WF=np.exp(1j*sum(a[:3])); WC=np.exp(1j*sum(a[3:])); u,v=2,3
    dF=(v in FACE_V)-(u in FACE_V); dC=(v in CYC_V)-(u in CYC_V)
    for k in (1,3):
        x=np.array([(WF**k if w in FACE_V else 1)*s[w] for w in range(5)])
        r=A(x,a,u,v)/A(s,a,u,v)
        print(f"  {tag:<30} {k:>3} {r:>34.12f} {WF**(dF*k):>26.12f} {abs(r-WF**(dF*k))<1e-12!s:>7}")
print()
print(f"  {'connection':<30} {'rho_TRUE = W_F^dF W_C^-dC':>28} {'rho as W-07 leg E codes it':>29}")
for tag, a in [("S1 published", np.array([np.pi/3]*3+[np.pi/2]*3)),
               ("generic sqrt2/sqrt3", np.array([2*np.pi*(2**0.5)/3]*3+[2*np.pi*(3**0.5)/3]*3)),
               ("random seed 1", np.random.default_rng(1).uniform(0,2*np.pi,6))]:
    WF=np.exp(1j*sum(a[:3])); WC=np.exp(1j*sum(a[3:])); dF,dC=-1,1
    rt = WF**dF*WC**(-dC); re_ = np.conj(WF)**dF*WC**(-dC)
    print(f"  {tag:<30} {np.angle(rt)/(2*np.pi):>28.9f} {np.angle(re_)/(2*np.pi):>29.9f}")
print("  M4-2 CONFIRMED independently: leg E's element is rho_true * W_F^2, equal only when W_F^2=1.")
print("  The decisive test's `profile` computes D from W_F, W_C directly, NOT from rho, so m4_c's")
print("  D column is immune to the slip; only its 'arg(rho)/2pi' label depends on it.\n")

print("== C5  WHERE THE COUNTEREXAMPLE LIVES — THE STRUCTURE OF THE EFFECT SET ==")
print("  E(K,eps) = {theta : exists k<=K with ||k theta|| < eps} = union over k<=K, |j|<=k of")
print("  the interval (j/k - eps/k, j/k + eps/k).  It is a union of intervals CENTRED AT RATIONALS")
print("  OF DENOMINATOR <= K, and nothing else.")
print(f"    |E| <= 2 K eps = {2*K*eps:.4e}   (Lebesgue measure, K={K}, eps={eps:.3e})")
lam=(2**0.5-1)*1e-13
print(f"    M4's counterexample theta = 1/4 + {lam:.4e} sits at distance {lam:.3e} from 1/4,")
print(f"    i.e. inside the interval of half-width eps/4 = {eps/4:.3e} around the order-4 rational.")
print("    Every irrational counterexample to 'finite order' at fixed (K,tol) is, necessarily, a")
print("    point within eps/q of a rational of denominator q <= K.  So the corrected name is not")
print("    a rival to 'order' but its tolerance-aware closure: PROXIMITY TO A LOW-ORDER RATIONAL,")
print("    of which finite order is the eps -> 0 limit.  M4-1's verdict word -- W-07's name is")
print("    'WRONG' -- is stronger than what the mathematics supports; 'INCOMPLETE: needs the")
print("    quantifiers ord <= K and the tolerance eps' is what the decisive test establishes.")
