# W-08 / M4 leg G — (1) RESOLVE attained-versus-approached, which W-07 sec7 declines to resolve.
#                   (2) EXACT-ARITHMETIC checks of every precision-sensitive claim in this lane.
#                   (3) A degeneracy on S1's published point that W-07 measured and did not name.
#
# ISOLATION LEDGER.  G1 holds carrier, state, observable and threshold fixed and moves K alone.
# G2 holds carrier, state, observable and K fixed and moves the THRESHOLD alone.  G3 is exact
# arithmetic on stated closed forms — no comparison.  G4 holds the connection at S1's published
# point and moves the READY STATE alone (published p versus generic p).
import numpy as np
from fractions import Fraction

print("== G1/G2  THE RESOLUTION: recurrence DENSITY as a function of TOLERANCE ==")
print("  D_k = 2 amp |sin(pi k theta)|.  Two regimes, both computed on the same code path:")
amp = 0.271776443            # W-07 leg E's own amplitude for A_23, seed 20260816
phi = (1+5**0.5)/2
th_gen = (1/phi**2)          # badly approximable, W-07's own generic choice
print(f"  amp = {amp} (W-07 leg E's), theta_generic = 1/phi^2 = {th_gen:.12f}, theta_pub = 1/4\n")
print(f"  {'tolerance':>12} | {'ord=4: density':>16} {'count/K':>12} | {'irrational: density':>20} {'count/K':>12}")
K = 10**6
k = np.arange(1, K+1)
frac = (th_gen*k) % 1.0
Dg = 2*amp*np.abs(np.sin(np.pi*frac))
Dp = 2*amp*np.abs(np.sin(np.pi*((0.25*k) % 1.0)))
for tol in [1e-1, 1e-2, 1e-3, 1e-4, 1e-6, 1e-9]:
    pred_irr = 2*np.arcsin(min(1.0, tol/(2*amp)))/np.pi
    print(f"  {tol:>12.0e} | {0.25:>16.6f} {int((Dp<tol).sum())/K:>12.6f} | {pred_irr:>20.3e} {int((Dg<tol).sum())/K:>12.3e}")
print()
print("  ATTAINED (ord finite = n): density of {D_k < tol} is 1/n for EVERY tol > 0, including 0.")
print("  APPROACHED (irrational):   density is 2 arcsin(tol/2amp)/pi ~ tol/(pi amp), -> 0 as tol -> 0.")
print()
print("  RESOLUTION.  The distinction is REAL but W-07 states it at the wrong place.  It is not")
print("  'zero versus not-zero' (unobservable at any finite precision, and unobservable in W-07's")
print("  own float64 output, which prints 2.565e-16 / 6.729e-19 / 1.349e-16 for the same exact 0).")
print("  It is the SCALING OF RECURRENCE DENSITY WITH TOLERANCE: constant 1/n versus O(tol).")
print("  At any fixed tolerance tol and any K < ~1/tol the two are indistinguishable — which is")
print("  exactly W-07's regime (tol = 1e-9, K = 4000) and why its table splits 1000 / 0 so cleanly.")
print("  So: A DISTINCTION WITH A DIFFERENCE, but the difference is Diophantine, not categorical,")
print("  and it is measurable only by varying the tolerance — which W-07 never does.\n")

print("== G3  EXACT ARITHMETIC ON THE PRECISION-SENSITIVE CLAIMS ==")
print("  (a) S1 published: rho = -i, rho^k = 1 iff 4 | k.  Exact over Z[i], no floats:")
def gauss_pow(z, n):
    r = (1, 0)
    for _ in range(n):
        r = (r[0]*z[0]-r[1]*z[1], r[0]*z[1]+r[1]*z[0])
    return r
print("      k :", [f"{k}->{gauss_pow((0,-1),k)}" for k in range(1, 9)])
print("      exact zeros among k<=4000 :", sum(1 for k in range(1, 4001) if gauss_pow((0,-1), k) == (1, 0)))
print("      W-07 leg B reported min 2.565e-16 and leg E reported 6.729e-19 for this exact 0.")
print("      Both are float residue; neither is a measurement.  W-07 sec3's table publishes")
print("      6.729e-19 as 'min over k<=4000' where the exact value is 0.  Its own C2 says 0.\n")

print("  (b) The INFINITE-ORDER counterexample of m4_c, verified in CLOSED FORM (no float):")
print("      theta = 1/4 + lam with lam = (sqrt2 - 1)*1e-13 irrational, so ord(rho) = infinity.")
print("      For k = 4m <= 4000:  ||k theta|| = ||m + 4m lam|| = 4m lam   (since 4000 lam < 1/2),")
print("      so D_k = 2 amp sin(pi 4m lam) <= 2 amp pi 4000 lam.")
lam = (2**0.5-1)*1e-13
bound = 2*amp*np.pi*4000*lam
print(f"        upper bound on D at every multiple of 4 : {bound:.6e}  <  1e-9   -> all 1000 qualify")
print("      For k not divisible by 4:  ||k theta|| >= 1/4 - 4000 lam, so")
print(f"        D_k >= 2 amp sin(pi(1/4 - 4000 lam)) = {2*amp*np.sin(np.pi*(0.25-4000*lam)):.6f}  -> none qualify")
print("      EXACTLY 1000 of 4000, with ord(rho) = INFINITY.  This is a closed-form fact, not a")
print("      float artefact: the only inputs are 4000*lam < 1/2 and sin monotone on [0, pi/2].\n")

print("  (c) W-07 sec7's caveat number 1.4e-12 is itself forced, not measured:")
print("      for irrational theta, min_{k<=K}(1-|Z_k|) ~ C/K^2 because 1-|Z| is quadratic in the")
print("      phase and min_k ||k theta|| ~ c/K.  Check the K-scaling on W-07's own generic case:")
P0, PF, PC = 0.5, 0.0, 0.5
aF, aC = 2*np.pi*phi, 2*np.pi*phi**2
for KK in [10**4, 10**5, 10**6]:
    kk = np.arange(1, KK+1)
    z = np.abs(P0*np.exp(1j*kk*(aC-aF)) + PF*np.exp(-1j*kk*aF) + PC*np.exp(1j*kk*aC))
    d = 1-z.max()
    print(f"        K = {KK:>8}   1 - max|Z_k| = {d:.3e}   (K^2 * that) = {d*KK**2:.4f}")
print("      The product is ~constant: the '1.4e-12' is 2.5/K^2 at K=1e6, i.e. a restatement of K.\n")

print("== G4  THE DEGENERACY W-07 MEASURED AND DID NOT NAME ==")
print("  On S1's PUBLISHED connection AND S1's PUBLISHED ready state p = (1/2,0,0,1/4,1/4):")
print("    P0 = 1/2 (v0), PF = 0 (v1,v2 carry no weight), PC = 1/2 (v3,v4)")
print("    Z_k = W_C^k [ P0 conj(W_F)^k + PC ] = (-i)^k [ (1/2)(-1)^k + 1/2 ]")
for kk in range(1, 9):
    val = Fraction(1, 2)*Fraction((-1)**kk) + Fraction(1, 2)
    print(f"      k={kk}:  |Z_k| = {val}  (exact)")
print("    |Z_k| = 1 for k EVEN and Z_k = 0 EXACTLY for k ODD.")
print("    So on the published point the branch comparison does not 'recur': it alternates with")
print("    PERIOD 2 between total overlap and total annihilation, and Omega_N = 0 for all N >= 1")
print("    (S3 audit COR-D, already of record).  The operative group element there is W_F, of")
print("    ORDER 2 — not rho, of order 4.  W-07 sec4 reports 500000 of 1e6 (= K/2, period 2) in")
print("    the same table as sec3's 1000 of 4000 (= K/4, period 4) and names ONE variable for both.")
print("    Two different group elements govern the two rows.  The name is wrong in sec4 as well,")
print("    and in a second way: at the published state there is no record for recurrence to undo.")
