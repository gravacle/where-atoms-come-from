#!/usr/bin/env python3
"""
X_03 — THE SUDLER ANALOGY IS DRAWN AGAINST AN OBJECT THAT NEVER OCCURS IN THIS CORPUS, AND
       THE FORCING R-8 ATTRIBUTES TO CODIMENSION IS ATTRIBUTABLE TO HOMOGENEITY.

TARGET.
  R-8 (R_07_THEOREMS.txt R-8):
    "In one dimension the singularity of log|1-x| sits at the group identity and Dirichlet
     FORCES min_{k<=N}||k alpha|| <= 1/N: near-hits are unavoidable for every alpha ...
     In two dimensions the singular set ... is a set of POINTS -- codimension 2 -- and nothing
     forces a near-hit at all. ... GENERICALLY OUR PROBLEM IS EASIER THAN SUDLER'S"
  R_04_2d_arms.out.txt, closing paragraph of the diagnostic block:
    "COMPARE the ONE-dimensional diagnostic of R_03, N * min ||k a - theta*||, also O(1):
     same shape, but there Dirichlet FORCES it to be O(1) and here nothing does".

TWO DEFECTS, ONE CONSEQUENCE.

(i) THE HOMOGENEOUS SUDLER PRODUCT NEVER ARISES IN THIS CORPUS, ON ANY LINE.  Q_{m,n}(1) =
    P(1,1) = sum(pi) = 1 for EVERY (m,n) and EVERY ready state, because pi is a probability
    vector.  So the reduced one-variable polynomial never vanishes at the orbit's own
    identity: every reduced problem is INHOMOGENEOUS.  R-5 says so on the (1,1) line; R-8 then
    contrasts "our problem" with the HOMOGENEOUS product, which is not the reduced problem
    anywhere.

(ii) DIRICHLET DOES NOT FORCE THE INHOMOGENEOUS DIAGNOSTIC.  Dirichlet bounds
     min_{k<=N}||k alpha||, the HOMOGENEOUS quantity.  For a shifted target it gives nothing.
     What forces N * min_{k<=N}||k alpha - g|| = O(1) is the THREE-DISTANCE theorem plus
     BOUNDED PARTIAL QUOTIENTS -- a property of alpha that the lane's four arms happen to
     have and that R_04's sentence attributes to Dirichlet instead.

CONSEQUENCE.  The 1-D/2-D contrast R-8 draws is really a HOMOGENEOUS/INHOMOGENEOUS contrast.
Both of the lane's own objects -- the reduced 1-D problem on a resonant line and the
unreduced 2-D problem -- are inhomogeneous, and NOTHING is forced in either.  R-8's
conclusion "generically easier than Sudler's" survives only as a statement about the RATE
(N^-1/2 against N^-1, a factor 2 in log of the deepest term, both contributing O(log N / N)),
not about the existence of a forcing.

ISOLATION.  ONE VARIABLE: THE SHIFT g IN min_{k<=N}||k alpha - g||, g = 0 against g = theta*.
alpha, the N-grid, and the reducer are byte-identical between the two columns of each row.
"""
import numpy as np
from fractions import Fraction
from X_lib import PI_K1, P_eval, frac_sqrt, frac_theta_star, ExactRot, relation_lattice

PREC = 60
SQ2 = frac_sqrt(2, PREC)
TH = frac_theta_star(PREC)
THf = float(TH)
p00, p10, p01, p11 = PI_K1
DECS = [10 ** i for i in range(1, 8)]
KMAX = 10 ** 7

print("=" * 79)
print("X_03 — HOMOGENEOUS vs INHOMOGENEOUS: THE FORCING IS NOT WHERE R-8 PUTS IT")
print("=" * 79)

# ---------------------------------------------------------------- (i) Q(1) = 1 always
print("\n(i)  Q_{m,n}(1) = P(1,1) = sum(pi) = 1 FOR EVERY (m,n).  Exhaustive check.")
bad = 0; cnt = 0
for m in range(-12, 13):
    for n in range(-12, 13):
        if (m, n) == (0, 0) or np.gcd(abs(m), abs(n)) != 1:
            continue
        cnt += 1
        q1 = p00 + p10 * 1.0 + p01 * 1.0 + p11 * 1.0     # P(1^n, 1^-m)
        if abs(q1 - 1.0) > 1e-15:
            bad += 1
print("     primitive (m,n) with |m|,|n| <= 12 scanned: %d   |Q_{m,n}(1) - 1| > 1e-15 in %d of them"
      % (cnt, bad))
print("     (the lane's own count for this box is 368 -- reproduced: %s)" % ("YES" if cnt == 368 else "NO, %d" % cnt))
print("     and on the d-twisted branches of X_01, Q_r(1) = 0.3 + 0.4 zeta^r + 0.3 zeta^r:")
for d in (2, 3, 5):
    vals = [abs(p10 + (p11 + p01) * np.exp(2j * np.pi * r / d)) for r in range(d)]
    print("        d=%d  |Q_r(1)| = %s   min = %.4f  -> never 0" % (d, " ".join("%.4f" % v for v in vals), min(vals)))
print("     SO THE CLASSICAL (HOMOGENEOUS) SUDLER PRODUCT IS NOT THE REDUCED OBJECT ANYWHERE.")

# ---------------------------------------------------------------- (ii) the two diagnostics
def diag(alpha_frac, shift, decs=DECS, kmax=KMAX):
    k = np.arange(1, kmax + 1, dtype=np.int64)
    f = ExactRot(alpha_frac).frac(k)
    d = np.abs(f - shift)
    d = np.minimum(d, 1.0 - d)
    out = []
    run = np.inf
    idx = 0
    for N in decs:
        run = min(run, d[idx:N].min())
        idx = N
        out.append(N * run)
    return out

DELTA = (SQ2 - 1) / 10 ** 9           # irrational, tiny
ALPHAS = [
    ("BA-silver   a = sqrt2 - 1", SQ2 - 1),
    ("BA-golden   a = (sqrt5-1)/2", (frac_sqrt(5, PREC) - 1) / 2),
    ("BAD-alpha   a = 1/2 + (sqrt2-1)*1e-9", Fraction(1, 2) + DELTA),
]
print("\n(ii) THE TWO DIAGNOSTICS ON THE SAME alpha.  ONE VARIABLE MOVED: THE SHIFT g.")
print("     column H : N * min_{k<=N} ||k a - 0||        (HOMOGENEOUS -- Dirichlet applies)")
print("     column I : N * min_{k<=N} ||k a - theta*||   (INHOMOGENEOUS -- Dirichlet does not)")
hdr = "     %-34s" % "N" + "".join("%13d" % N for N in DECS)
for lbl, a in ALPHAS:
    H = diag(a, 0.0)
    I = diag(a, THf)
    print("\n" + hdr)
    print("     %-34s" % (lbl + "   H") + "".join("%13.4f" % v for v in H))
    print("     %-34s" % (lbl + "   I") + "".join("%13.4f" % v for v in I))
print("""
     READ.  H is <= 1 in every row and at every decade, exactly as Dirichlet forces.
     I is O(1) for the two badly-approximable rotations -- but it is NOT forced to be, and
     BAD-alpha shows it: I GROWS LINEARLY OVER ALL SEVEN DECADES while H stays under 1 on the
     SAME alpha and the SAME grid.  Only the shift moved.  BAD-alpha is irrational, and with
     beta = 1 - alpha it sits on the target lane's own singular line u v = 1 (c = f), so this
     is not an object beside the corpus: it is a connection on it (W-12 Corollary 1).""")

# ---------------------------------------------------------------- what does force it
print("\n(iii) WHAT ACTUALLY FORCES I = O(1): the THREE-DISTANCE theorem + bounded partial")
print("      quotients.  Largest gap of {k a : k <= N} times N (if this is O(1), EVERY target")
print("      is within O(1/N) and I is forced for every shift g):")
def maxgap(alpha_frac, decs, kmax):
    k = np.arange(1, kmax + 1, dtype=np.int64)
    f = ExactRot(alpha_frac).frac(k)
    out = []
    for N in decs:
        s = np.sort(f[:N])
        g = np.max(np.diff(np.concatenate(([s[-1] - 1.0], s))))
        out.append(N * g)
    return out
D2 = [10 ** i for i in range(1, 6)]
print("     %-34s" % "N" + "".join("%13d" % N for N in D2))
for lbl, a in ALPHAS:
    print("     %-34s" % lbl + "".join("%13.4f" % v for v in maxgap(a, D2, 10 ** 5)))
print("""     BA-silver and BA-golden: bounded -> the forcing is real for THEM, and it comes from
     bounded partial quotients, not from Dirichlet.  BAD-alpha: grows -> no forcing.
     BAD-alpha's partial quotients are unbounded; its rotation is still irrational and its
     orbit is still dense and equidistributed.""")

# ---------------------------------------------------------------- the Birkhoff consequence
print("\n(iv) THE CONSEQUENCE ON THE LANE'S OWN SINGULAR LINE (beta = 1 - alpha, c = f).")
print("     S_N against the (1,1) subtorus value log(0.3) = %.9f" % np.log(0.3))
k = np.arange(1, KMAX + 1, dtype=np.int64)
print("     %-34s" % "N" + "".join("%13d" % N for N in DECS))
for lbl, a in ALPHAS:
    b = 1 - a
    fa = ExactRot(a).frac(k); fb = ExactRot(b).frac(k)
    lz = np.log(np.abs(P_eval(PI_K1, np.exp(2j * np.pi * fa), np.exp(2j * np.pi * fb))))
    cs = np.cumsum(lz)
    print("     %-34s" % lbl + "".join("%13.2e" % (cs[N - 1] / N - np.log(0.3)) for N in DECS))
    del fa, fb, lz, cs
print("""     BAD-alpha has NOT reached the subtorus value inside seven decades: it sits +0.399 above
     log(0.3) at every decade from 10 to 1e7, because the approach cannot begin until
     N ~ 1/delta ~ 2.4e9, and delta is free.  Its limit IS log(0.3) -- the closure is the
     circle and equidistribution holds -- but the onset is unbounded.  So the resonant locus
     is not 'where the limit is the subtorus value' in any sense a seven-decade window can
     see either.  A big inhomogeneous approximation constant is the BENIGN direction (it makes
     the sum too large, not too small); the malign direction is R-3's ladder.  What X_03
     refutes is R-8's MECHANISM and R_04's 'Dirichlet FORCES it' sentence -- not R-8's verdict
     that the two problems differ, which survives as a statement about the RATE.""")
print("\nDONE X_03")
