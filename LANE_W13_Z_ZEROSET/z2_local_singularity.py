"""LANE W-13 / Z  --  z2: THE LOCAL BEHAVIOUR OF |P| NEAR A TORUS ZERO.
Item (2) of the brief.  This is the block that decides how integrable the singularity is and
what Diophantine condition a Birkhoff sum needs.  Two independent routes to every constant."""
import sys, math
import numpy as np
from fractions import Fraction as F
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from z0_lib import (strat_exact, zeros_closed_form, jensen_mods, sublevel_measure,
                    local_alpha_beta, local_singvals, Pval, fr)

W = 96
def hdr(s):
    print("=" * W); print(s); print("=" * W)

hdr("z2  THE LOCAL SINGULARITY: |P| ~ c(theta) r AT AN ISOLATED ZERO, AND THE SUBLEVEL EXPONENT")
print("numpy", np.__version__, "; IEEE double unless a line says EXACT.\n")

K1  = fr(F(0),     F(3, 10), F(3, 10), F(2, 5))
TAN = fr(F(1, 10), F(1, 5),  F(3, 10), F(2, 5))
CUR = fr(F(3, 10), F(1, 5),  F(1, 5),  F(3, 10))
S1P = fr(F(0),     F(0),     F(1, 2),  F(1, 2))

# =======================================================================================
print("-" * W)
print("(a) THE LOCAL EXPANSION, DERIVED.")
print("""
    Put x = x0 e^{i sigma}, y = y0 e^{i tau} at a zero (x0,y0).  Then
        P = i (alpha sigma + beta tau) - (1/2)(a sigma^2 + b tau^2 + c (sigma+tau)^2) + O(r^3)
    with  a = p10 x0,  b = p01 y0,  c = p11 x0 y0, and
        alpha = x0 dP/dx = a + c = x0 (p10 + p11 y0)
        beta  = y0 dP/dy = b + c = y0 (p01 + p11 x0) = y0 B(x0).
    (sigma,tau) |-> i(alpha sigma + beta tau) is a REAL-LINEAR map R^2 -> C = R^2 with matrix
        M = [[Re alpha, Re beta],[Im alpha, Im beta]],   det M = Im(conj(alpha) beta).
    SO:
      * det M != 0  ==>  |P| = |M (sigma,tau)| (1 + O(r))  with  s_min r <= |P| <= s_max r:
        a NON-DEGENERATE (simple, conical) zero.  log|P| ~ log r.  INTEGRABLE.
        The sublevel set { |P| < eps } near this zero is an ELLIPSE of area pi eps^2/|det M|.
      * det M = 0  ==>  alpha, beta are REAL-PARALLEL; the image of the linear part is a LINE,
        and along its kernel |P| vanishes to SECOND order.  A degenerate (tangential) zero.
    THE ARITHMETIC FACT THAT MAKES THIS CLEAN: det M = 0 happens EXACTLY on the ONE-point and
    CURVE strata of z1(a).  Transversal branch crossing <=> simple zero.  Verified below.
""")

print("    ZERO-BY-ZERO LOCAL DATA (all named states of z1 with isolated zeros):")
print(f"    {'state':<9s} {'alpha':>26s} {'beta':>26s} {'s_min':>10s} {'s_max':>10s} {'det':>12s}")
from z0_lib import NAMED
for name, p in NAMED:
    zz = zeros_closed_form(p)
    if not zz:
        continue
    for (x0, y0, cs, sn2, sgn) in zz:
        al, be = local_alpha_beta(p, x0, y0)
        smx, smn, det = local_singvals(al, be)
        print(f"    {name.split()[0]:<9s} {al.real:+11.6f}{al.imag:+11.6f}i "
              f"{be.real:+11.6f}{be.imag:+11.6f}i {smn:10.6f} {smx:10.6f} {det:12.3e}")
print()
print("    ==> every TWO-point state has det != 0 (simple zeros); the TANGENT state's single")
print("        zero has det = 0.  Simplicity is DECIDED by the stratum, not assumed.\n")

# =======================================================================================
print("-" * W)
print("(b) |P| / r ACROSS SIX DECADES OF r, AT K1's REGISTERED pi.  THE PREDICTED BAND IS")
print("    [s_min, s_max]; the measured min and max over 4096 directions must sit in it and")
print("    must STOP MOVING as r falls.  A single r would prove nothing (CONVERGENCE IS NOT")
print("    A WINDOW); six decades and the trend are printed.\n")
zz = zeros_closed_form(K1)
x0, y0 = zz[0][0], zz[0][1]
al, be = local_alpha_beta(K1, x0, y0)
smx, smn, det = local_singvals(al, be)
print(f"    zero at x0 = {x0:.12f}, y0 = {y0:.12f}")
print(f"    PREDICTED  s_min = {smn:.12f}   s_max = {smx:.12f}   det = {det:.12f}")
th = np.linspace(0, 2 * np.pi, 4097)[:-1]
print(f"    {'r':>10s} {'min |P|/r':>16s} {'max |P|/r':>16s} {'min/s_min - 1':>16s} {'max/s_max - 1':>16s}")
for e in range(2, 10):
    r = 10.0 ** (-e)
    sg, ta = r * np.cos(th), r * np.sin(th)
    xs = x0 * np.exp(1j * sg); ys = y0 * np.exp(1j * ta)
    vals = np.abs(Pval(K1, xs, ys)) / r
    print(f"    {r:10.1e} {vals.min():16.12f} {vals.max():16.12f} "
          f"{vals.min()/smn - 1:16.3e} {vals.max()/smx - 1:16.3e}")
print("    TREND: both ratios converge to 1 like O(r).  |P| ~ c(theta) r with c in [s_min,s_max],")
print("    c > 0 in EVERY direction.  THE ZERO IS SIMPLE AND CONICAL.  log|P| ~ log r.\n")

# =======================================================================================
print("-" * W)
print("(c) THE DEGENERATE (TANGENTIAL) ZERO: |P| ~ r ALONG MOST DIRECTIONS BUT ~ r^2 ALONG")
print("    THE KERNEL LINE.  Same code path, same r-ladder; the ONE thing that moves is the")
print("    weight vector, hence the stratum.\n")
zz = zeros_closed_form(TAN)
x0t, y0t = zz[0][0], zz[0][1]
alt, bet = local_alpha_beta(TAN, x0t, y0t)
smxt, smnt, dett = local_singvals(alt, bet)
print(f"    TANGENT pi = {tuple(str(q) for q in TAN)}   zero at ({x0t:.6f}, {y0t:.6f})")
print(f"    alpha = {alt:.9f}   beta = {bet:.9f}   det = {dett:.3e}  (BOTH REAL -> rank 1)")
kdir = np.array([bet.real, -alt.real]); kdir = kdir / np.linalg.norm(kdir)
print(f"    kernel direction (sigma,tau) = ({kdir[0]:+.9f}, {kdir[1]:+.9f})")
print(f"    {'r':>10s} {'|P| along kernel':>20s} {'/r':>14s} {'/r^2':>14s} {'min over dirs':>16s}")
for e in range(2, 10):
    r = 10.0 ** (-e)
    xs = x0t * np.exp(1j * r * kdir[0]); ys = y0t * np.exp(1j * r * kdir[1])
    v = abs(Pval(TAN, xs, ys))
    sg, ta = r * np.cos(th), r * np.sin(th)
    allv = np.abs(Pval(TAN, x0t * np.exp(1j * sg), y0t * np.exp(1j * ta)))
    print(f"    {r:10.1e} {v:20.12e} {v/r:14.6e} {v/r**2:14.9f} {allv.min():16.6e}")
print("    TREND: |P|/r  -> 0   and   |P|/r^2 -> a POSITIVE CONSTANT along the kernel.")
print("    A SECOND-ORDER (tangential) zero.  log|P| ~ 2 log r on a line, still integrable,")
print("    but the SUBLEVEL SET is fatter -- quantified in (d).\n")

# =======================================================================================
print("-" * W)
print("(d) THE SUBLEVEL MEASURE mu(eps) = Haar{ |P| < eps }, THE OBJECT THAT ACTUALLY DECIDES")
print("    INTEGRABILITY AND THE DIOPHANTINE CONDITION.  Computed by the EXACT y-arc formula")
print("    (a 1-D quadrature of a CONTINUOUS integrand -- no 2-D grid, no noise floor).")
print("    Six decades of eps; the LOCAL slope d log mu / d log eps is printed per decade.\n")

def slope_table(p, label, predict=None):
    typ, det = strat_exact(p)
    print(f"    {label}   pi = {tuple(str(q) for q in p)}   stratum {typ}{('-'+det) if det else ''}")
    eps = [10.0 ** (-e) for e in range(1, 8)]
    mu = [sublevel_measure(p, e) for e in eps]
    print(f"      {'eps':>10s} {'mu(eps)':>16s} {'mu/eps^theta':>18s} {'local slope':>13s}")
    for i, (e, m) in enumerate(zip(eps, mu)):
        sl = float('nan')
        if i and mu[i] > 0 and mu[i - 1] > 0:
            sl = math.log(mu[i] / mu[i - 1]) / math.log(eps[i] / eps[i - 1])
        nrm = m / e ** predict if predict else float('nan')
        print(f"      {e:10.1e} {m:16.9e} {nrm:18.9e} {sl:13.6f}")
    print(f"      LAST-DECADE slope = {sl:.6f}   PREDICTED theta = {predict}")
    print()

print("    ARM 1 -- TWO simple conical zeros.  PREDICTED theta = 2, and PREDICTED constant")
print("             mu/eps^2 -> (1/(4 pi)) * SUM_zeros 1/|det M_j|.")
pred = 0.0
for (x0_, y0_, *_r) in zeros_closed_form(K1):
    a_, b_ = local_alpha_beta(K1, x0_, y0_)
    pred += 1.0 / abs(local_singvals(a_, b_)[2])
pred /= (4 * math.pi)
print(f"             PREDICTED limit of mu/eps^2 = {pred:.12f}")
slope_table(K1, "K1_REG ", predict=2)

print("    ARM 2 -- ONE tangential zero.  PREDICTED theta = 3/2.")
slope_table(TAN, "TANGENT", predict=1.5)

print("    ARM 3 -- a CURVE (stratum III).  PREDICTED theta = 1.")
slope_table(CUR, "CURVE3 ", predict=1.0)

print("    ARM 4 -- a CURVE (stratum I) = K1's OWN S1-PUBLISHED READY STATE.  PREDICTED 1.")
slope_table(S1P, "S1_PUB ", predict=1.0)

print("    ARM 5 -- EMPTY zero set (B0b's own uniform state, S4:575).  |P| >= 1/9 everywhere,")
print("             so mu(eps) = 0 for eps <= 1/9 and log|P| is CONTINUOUS AND BOUNDED.")
B0b = fr(F(4, 9), F(2, 9), F(2, 9), F(1, 9))
for e in (1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7):
    print(f"      eps = {e:8.1e}   mu = {sublevel_measure(B0b, e):.6e}")
print()

# =======================================================================================
print("-" * W)
print("(e) THE ONE-VARIABLE ISOLATION, DESIGNED FROM z1(a)'s STRATUM GEOMETRY.")
print("""
    The four types are STRATA: {EMPTY} and {TWO} are open, {ONE} is the codimension-1 wall
    between them, and the three CURVE families are codimension 2 INSIDE that wall.  So:
      * a generic one-parameter path crosses the wall transversally: EMPTY -> ONE -> TWO.
        FAMILY A below is such a path, in ONE declared scalar eta.
      * NO one-parameter path reaches a CURVE transversally.  A curve arm therefore has to
        be reached by a path that STAYS IN THE WALL.  FAMILY B below is such a path, in ONE
        declared scalar mu, with (S1-S2)(D1-D2) == 0 EXACTLY at every mu.
    Both families move exactly one scalar.  BOTH ARMS ARE DIFFED -- inputs and outputs.
""")
print("    FAMILY A  pi(eta) = (1/10, 1/5, 3/10 + eta, 2/5 - eta).   ONE scalar: eta.")
print(f"    {'eta':>10s} {'stratum':>10s} {'(S1-S2)(D1-D2)':>18s} {'theta (eps 1e-4->1e-6)':>24s}")
for eta in (F(-3, 100), F(-1, 100), F(0), F(1, 100), F(3, 100), F(1, 20)):
    p = (F(1, 10), F(1, 5), F(3, 10) + eta, F(2, 5) - eta)
    S1_, S2_ = p[0] + p[1], p[2] + p[3]
    D1_, D2_ = abs(p[0] - p[1]), abs(p[2] - p[3])
    t, d = strat_exact(p)
    m1 = sublevel_measure(p, 1e-4); m2 = sublevel_measure(p, 1e-6)
    th_ = math.log(m2 / m1) / math.log(1e-2) if m1 > 0 and m2 > 0 else float('nan')
    print(f"    {str(eta):>10s} {t:>10s} {float((S1_-S2_)*(D1_-D2_)):18.6e} {th_:24.6f}")
print()
print("    FAMILY B  pi(mu) = (1/10, 2/5, mu, 1/2 - mu).  ONE scalar: mu.  S1 = S2 = 1/2 for")
print("              EVERY mu, so the path never leaves the wall; the CURVE strata are hit")
print("              exactly at mu = 1/10 (stratum II) and mu = 2/5 (stratum III).")
print(f"    {'mu':>10s} {'stratum':>12s} {'theta (eps 1e-4->1e-6)':>24s} {'mu(1e-4)':>14s}")
for mu in (F(0), F(1, 20), F(1, 10), F(3, 20), F(1, 4), F(7, 20), F(2, 5), F(9, 20), F(1, 2)):
    p = (F(1, 10), F(2, 5), mu, F(1, 2) - mu)
    t, d = strat_exact(p)
    m1 = sublevel_measure(p, 1e-4); m2 = sublevel_measure(p, 1e-6)
    th_ = math.log(m2 / m1) / math.log(1e-2) if m1 > 0 and m2 > 0 else float('nan')
    print(f"    {str(mu):>10s} {(t+('-'+d if d else '')):>12s} {th_:24.6f} {m1:14.6e}")
print()
print("    ARMS DIFFED (the defect W-08's isolation audit calls the commonest FATAL one:")
print("    byte-identical arms reported as a confirmation).  Family A eta = -1/100 vs +1/100:")
pA = (F(1, 10), F(1, 5), F(3, 10) - F(1, 100), F(2, 5) + F(1, 100))
pB = (F(1, 10), F(1, 5), F(3, 10) + F(1, 100), F(2, 5) - F(1, 100))
print(f"      arm1 pi = {tuple(str(q) for q in pA)}   arm2 pi = {tuple(str(q) for q in pB)}")
print(f"      inputs differ in {sum(1 for a,b in zip(pA,pB) if a!=b)} of 4 coordinates")
print(f"      strata: {strat_exact(pA)[0]}  vs  {strat_exact(pB)[0]}")
print(f"      mu(1e-5): {sublevel_measure(pA,1e-5):.6e}  vs  {sublevel_measure(pB,1e-5):.6e}   "
      f"ratio {sublevel_measure(pB,1e-5)/max(sublevel_measure(pA,1e-5),1e-300):.3e}")
print()

# =======================================================================================
print("-" * W)
print("(f) QUADRATURE STABILITY, AND THE WINDOW ARTEFACT SHOWN RATHER THAN HIDDEN.")
print("    LEFT: the UNIFORM-mesh estimator at three mesh sizes.  RIGHT: the crossing-refined")
print("    estimator.  Once eps falls below the uniform mesh spacing the uniform estimator")
print("    COLLAPSES TO ZERO -- which reads as 'the singularity went away' and is exactly the")
print("    defect class COR-E convicts.  This lane uses the refined estimator everywhere.\n")
from z0_lib import sublevel_measure_uniform
print(f"    {'eps':>10s} {'uniform 2^18':>16s} {'uniform 2^20':>16s} {'uniform 2^22':>16s} {'REFINED':>16s}")
for e in (1e-3, 1e-4, 1e-5, 1e-6, 1e-7):
    row = [sublevel_measure_uniform(K1, e, n) for n in (1 << 18, 1 << 20, 1 << 22)]
    print(f"    {e:10.1e} " + " ".join(f"{v:16.9e}" for v in row)
          + f" {sublevel_measure(K1, e):16.9e}")
print()
print("    REFINED estimator, self-consistency at two fine-mesh sizes (must agree):")
for lbl, p in (("K1_REG", K1), ("CURVE3", CUR), ("TANGENT", TAN)):
    a1 = sublevel_measure(p, 1e-5, 1 << 16, 1 << 14)
    a2 = sublevel_measure(p, 1e-5, 1 << 17, 1 << 16)
    print(f"    {lbl:<8s} mu(1e-5) = {a1:.9e}  vs  {a2:.9e}   rel diff "
          f"{abs(a2-a1)/max(a2,1e-300):.2e}")
print()

print("-" * W)
print("(g) WHAT THE LOCAL BEHAVIOUR DECIDES.")
print("""
    log|P| is in L^1(T^2) in EVERY stratum -- m(P) is finite always -- but the three strata
    give three different SUBLEVEL EXPONENTS theta, and theta is the quantity that a Birkhoff
    sum feels:

        stratum          Z(P)                theta      what the orbit must avoid
        ---------------  ------------------  ---------  -----------------------------------
        EMPTY            no zeros            (none)     NOTHING.  log|P| is CONTINUOUS and
                                                        BOUNDED on T^2.
        TWO   (generic)  2 conjugate points  2          two POINTS of T^2 (codimension 2)
        ONE   (wall)     1 tangential point  3/2        one POINT, approached along a
                                                        parabola
        CURVE (cod 2)    a circle            1          a CIRCLE (codimension 1)

    The truncation error of the eps-regularisation is exactly governed by theta:
        INT (max(log|P|,log eps) - log|P|) = INT_0^eps log(eps/u) d mu(u) = C eps^theta/theta,
    so at fixed eps the CURVE case costs eps^1 where the generic case costs eps^2.

    AND THE DISCRETE SIDE, WHICH IS WHERE THE CONVERGENCE QUESTION LIVES:
      * a Kronecker orbit of length N has minimum distance ~ N^{-1/2} to a fixed POINT and
        ~ N^{-1} to a fixed CIRCLE, so the CURVE case is the one the orbit meets often;
      * an EXACT hit needs (u^k,v^k) = (x0,y0): TWO simultaneous conditions in the point
        strata -- a countable set of (u,v) -- but only ONE condition in the curve strata,
        e.g. u^k = -1 in stratum I, whose solution set is a FINITE UNION OF CIRCLES in the
        connection torus for each k, and DENSE over k.
      THAT is the sense in which the curve case is materially harder: the bad set of
      connections goes from codimension 2 to codimension 1.
""")
print("DONE z2")
