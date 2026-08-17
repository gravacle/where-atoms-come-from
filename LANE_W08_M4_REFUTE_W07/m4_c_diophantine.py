# W-08 / M4 leg C — THE DECISIVE TEST.  IS ord(rho) THE OPERATIVE VARIABLE, OR IS THE NAME WRONG?
#
# W-07: "THE OPERATIVE VARIABLE IS ord(rho), FINITE VERSUS INFINITE."  A dichotomy has two edges and
# both are tested here, on the CORRECT ratio rho = W_F^dF W_C^-dC (see m4_b), observable A_23.
#
#   EDGE 1 (finite must give the effect):  ord(rho) = 4001, large but FINITE.
#   EDGE 2 (infinite must not):            rho irrational but well approximable at denominator 4.
#
# ISOLATION LEDGER.  Held fixed across every row: carrier K1, ready state s (seed 20260816, W-07's
# own), observable A_23, dressing tree, k-range 1..4000, threshold 1e-9, code path (one function,
# `profile`, used for all rows).  Moved between rows: arg(rho)/2pi ALONE — the connection is built
# by handing that single number to the same constructor.  W_F is held at -1 in every row of C2 so
# that "the curvature leaves the value -1" is ruled OUT as the moving variable.
#
# Double precision (numpy float64).  The two precision-sensitive claims — that a stated rational
# angle gives EXACT zeros, and that a stated irrational angle NEVER does — are not decided here;
# they are decided in m4_g_exact.py by integer arithmetic.  What is decided here is the OBSERVABLE
# behaviour at K=4000 and threshold 1e-9, which is what W-07's table reports.
import numpy as np

FACE_V = {0, 1, 2}; CYC_V = {0, 3, 4}; TREE = {1: (0,), 2: (0, 1), 3: (3,), 4: (3, 4)}
def dress(s, a):
    u = np.exp(1j*np.asarray(a)); t = np.array(s, dtype=complex)
    for v, p in TREE.items():
        w = 1.0+0j
        for e in p: w *= u[e]
        t[v] = s[v]/w
    return t
rng = np.random.default_rng(20260816)
s = rng.normal(size=5) + 1j*rng.normal(size=5); s /= np.linalg.norm(s)
K = 4000; k = np.arange(1, K+1)

def connection_from(theta_rho, argWF=np.pi):
    """Build a K1 connection with arg(W_F)=argWF and rho = W_F^-1 W_C^-1 = exp(2pi i theta_rho)."""
    # rho = W_F^-1 W_C^-1  =>  arg W_C = -argWF - 2pi theta_rho
    argWC = -argWF - 2*np.pi*theta_rho
    return np.array([argWF/3]*3 + [argWC/3]*3)

def profile(a, u=2, v=3):
    WF = np.exp(1j*sum(a[:3])); WC = np.exp(1j*sum(a[3:]))
    dF = (v in FACE_V)-(u in FACE_V); dC = (v in CYC_V)-(u in CYC_V)
    t = dress(s, a); amp = abs(np.conj(t[u])*t[v])
    D = amp*np.abs(WF**(dF*k) - WC**(dC*k))
    rho = WF**dF*WC**(-dC)
    return D, amp, rho

def order_of(theta, cap=100000):
    """order of exp(2pi i theta) when theta is given as an exact rational p/q; else None."""
    from fractions import Fraction
    if isinstance(theta, Fraction): return theta.denominator
    return None

from fractions import Fraction
print("== C1  THE DICHOTOMY, BOTH EDGES.  amp is printed so nothing hides in it. ==")
print(f"  ready state = W-07 leg E's own (seed 20260816); observable A_23; K = {K}; threshold 1e-9")
print()
rows = [
    ("ord(rho) = 4     [S1 PUBLISHED]        ", Fraction(-1, 4), "finite, <= K"),
    ("ord(rho) = 4001  [FINITE, > K]         ", Fraction(1, 4001), "finite, >  K"),
    ("ord(rho) = 8000  [FINITE, > K]         ", Fraction(1, 8000), "finite, >  K"),
    ("ord(rho) = 2000  [FINITE, < K]         ", Fraction(1, 2000), "finite, <  K"),
]
print(f"  {'connection':<42} {'arg(rho)/2pi':>16} {'ord(rho)':>9} {'min D':>12} {'cells<1e-9':>11}")
for tag, th, note in rows:
    a = connection_from(float(th)); D, amp, rho = profile(a)
    print(f"  {tag:<42} {float(th):>+16.10f} {th.denominator:>9} {D.min():>12.3e} {int((D<1e-9).sum()):>11}")

# EDGE 2: irrational theta, but with a strong rational approximation at denominator 4.
lam = (np.sqrt(2)-1)*1e-13                     # irrational  =>  theta irrational  =>  ord(rho) = infinity
th_irr = 0.25 + lam
a = connection_from(th_irr); D, amp, rho = profile(a)
print(f"  {'ord(rho) = INFINITE  theta = 1/4 + (sqrt2-1)e-13':<42} {th_irr:>+16.10f} {'inf':>9} {D.min():>12.3e} {int((D<1e-9).sum()):>11}")
lam2 = (np.sqrt(2)-1)*1e-16
th_irr2 = 0.25 + lam2
a2 = connection_from(th_irr2); D2, amp2, rho2 = profile(a2)
print(f"  {'ord(rho) = INFINITE  theta = 1/4 + (sqrt2-1)e-16':<42} {th_irr2:>+16.10f} {'inf':>9} {D2.min():>12.3e} {int((D2<1e-9).sum()):>11}")
print()
print("  BOTH EDGES OF W-07's DICHOTOMY FAIL:")
print("   * FINITE order 4001 and 8000 behave EXACTLY like W-07's 'generic, infinite' rows: 0 of 4000.")
print("   * INFINITE order reproduces W-06's registered figure 1000 of 4000 with min ~1e-13.")
print("  What separates the rows is not finiteness of ord(rho).  It is  min_{k<=K} ||k*theta||,")
print("  the quality of rational approximation to theta = arg(rho)/2pi at denominators <= K.")
print()

print("== C2  THE SAME CUT WITH W_F HELD AT -1 THROUGHOUT — 'the curvature left -1' RULED OUT ==")
print("  Every row above already holds arg(W_F) = pi.  Confirming, and printing W_F per row:")
for tag, th, note in rows[:2]:
    a = connection_from(float(th)); WF = np.exp(1j*sum(a[:3])); WC = np.exp(1j*sum(a[3:]))
    D, amp, rho = profile(a)
    print(f"    {tag}  W_F = {WF:+.6f}  W_C = {WC:+.6f}  cells<1e-9 = {int((D<1e-9).sum())}")
a = connection_from(th_irr); WF = np.exp(1j*sum(a[:3])); WC = np.exp(1j*sum(a[3:])); D, amp, rho = profile(a)
print(f"    ord = INFINITE (theta=1/4+eps)               W_F = {WF:+.6f}  W_C = {WC:+.6f}  cells<1e-9 = {int((D<1e-9).sum())}")
print("  W_F = -1 in all three.  The curvature is NOT the moving variable.  Neither is <W_F,W_C>")
print("  generating a finite subgroup of U(1)^2: in the INFINITE row it generates a dense subgroup")
print("  of a circle in T^2 and the effect is present anyway.\n")

print("== C3  WHY 'FIVE FOR FIVE' COULD NOT HAVE COME OUT OTHERWISE — A ONE-LINE DENSITY COUNT ==")
print("  D_k = 2*amp*|sin(pi k theta)|, so D_k < d  <=>  ||k theta|| < arcsin(d/(2 amp))/pi.")
print("  For irrational theta, Weyl: #{k<=K : ||k theta|| < e}/K -> 2e.  Expected count = 2 K e.")
for amp_ in [0.2, 0.27, 0.3]:
    e = np.arcsin(1e-9/(2*amp_))/np.pi
    print(f"    amp = {amp_:.2f}:  e = {e:.3e}   expected #cells<1e-9 in K=4000 : {2*K*e:.3e}")
print("  A random connection returns 0 of 4000 with probability ~ 1 - 1e-6.  Drawing five and")
print("  getting 0 five times is not an isolation; it is the density estimate, restated.")
print("  BY W-07's OWN STANDARD (sec 6: 'could not have failed voids a CONTROL'), W-07's five")
print("  generic rows ARE a control and they could not have failed.  The disqualifier it applies")
print("  to W-06's sweep applies to W-07's own table, where — unlike the sweep — it is apt.\n")

print("== C4  THE THEOREM THAT KILLS THE DICHOTOMY (no computation required) ==")
print("  For any K, any n <= K and any d > 0 there is an IRRATIONAL theta with")
print("      #{k <= K : |exp(2pi i k theta) - 1| < d}  =  floor(K/n).")
print("  Proof: at theta0 = 1/n the count is floor(K/n), the qualifying terms are exactly 0 and the")
print("  non-qualifying ones are bounded below by 2|sin(pi/n)| > 0.  Both conditions are open in")
print("  theta and finite in number, so they persist on a neighbourhood of 1/n; irrationals are")
print("  dense in it.  QED.  Finiteness of ord(rho) is therefore neither necessary nor (C1 rows")
print("  4001, 8000) sufficient for the effect W-07 measures at fixed K and fixed threshold.")
print()
print("  WHAT IS TRUE, AND IT IS THE ONLY SURVIVING FORM:")
print("   (i)  D_k = 0 EXACTLY for some k <= K   <=>   ord(rho) is finite AND ord(rho) <= K.")
print("        Not 'finite versus infinite' — finite AND at most K.  Rows 4001/8000 are the")
print("        counterexamples to the version W-07 registered.")
print("   (ii) D_k < d for some k <= K  is governed by min_{k<=K} ||k theta||, a DIOPHANTINE")
print("        quantity that is continuous in theta and blind to rationality.")
