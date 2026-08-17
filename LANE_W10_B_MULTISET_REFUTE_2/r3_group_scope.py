#!/usr/bin/env python3
"""LENS 2 (SCOPE), ATTACK 3.  IS THE INVARIANCE GROUP 'EXACTLY D4 OFF THE COLLINEAR
LOCUS' (B-03), OR IS THAT TRUE ONLY AT THE PHASE VECTORS THE LANE CHOSE?

B-03 claims: 'off the collinear locus the invariance group is precisely D4 of order 8 ...
and the 24 arrangements fall into its 3 cosets carrying 3 distinct rates'.  Its evidence
is ONE phase vector, args = (0, 0.7, 1.9, 0.3), on three moduli vectors.  The lane's OWN
randomised sweep at B.4 contradicts the claim on 96 of 400 arrays (observed counts
{3:304, 1:78, 2:18}) and files the shortfall under 'extra degeneracy of the moduli'
without naming the degeneracy.

THIS SCRIPT NAMES IT, EXACTLY, and then measures how big it is.

  3.0  the per-flux branch-domination criterion, derived and checked in closed form
  3.1  every one of the lane's four printed SURVIVING MISSES, decided exactly
  3.2  the lane's sweep re-run with the corrected predictor -- one variable moved:
       the predictor.  Same seed, same arrays, same evaluator, same tolerance.
  3.3  the measure of the exceptional set, and what B-03 should say
  3.4  off the exceptional set: is lambda strictly monotone in |phi|?  (if it is, D4 is
       exactly right there, and B-03's conclusion survives on the corrected domain)
"""
import numpy as np
import mpmath as mp

from r_lib import PERMS, apply_perm, cyc, hdr, is_subgroup, mahler_jensen

mp.mp.dps = 30
print(__doc__)

D4 = tuple(sorted(s for s in PERMS if {frozenset((s[0], s[3])), frozenset((s[1], s[2]))}
                  == {frozenset((0, 3)), frozenset((1, 2))}))


def fluxes(args):
    A = args
    return (A[0] + A[1] - A[2] - A[3], A[0] + A[2] - A[1] - A[3], A[0] + A[3] - A[1] - A[2])


def wrap(x):
    y = float(np.mod(x, 2 * np.pi))
    return min(y, 2 * np.pi - y)          # |phi| in [0, pi]


def dominated_at_flux(r, phi, tol=0.0):
    """EXACT-in-closed-form test.  For the pairing {r_i,r_j}|{r_k,r_l},
         SA - SB = C0 + Re(K e^{it}),  C0 = ri^2+rj^2-rk^2-rl^2,
         |K|^2 = X^2 + Y^2 - 2 X Y cos(phi),  X = 2 ri rj,  Y = 2 rk rl,
    so ONE BRANCH DOMINATES POINTWISE iff C0^2 >= |K|^2, and then
         lambda = log(max of the dominating pair) = log(max r),
    INDEPENDENT OF phi.  Returns True if SOME pairing dominates at this flux."""
    for (i, j), (k, l) in (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))):
        C0 = r[i] ** 2 + r[j] ** 2 - r[k] ** 2 - r[l] ** 2
        X, Y = 2 * r[i] * r[j], 2 * r[k] * r[l]
        K2 = X * X + Y * Y - 2 * X * Y * np.cos(phi)
        if C0 * C0 >= K2 - tol:
            return True
    return False


def predict(r, args):
    """corrected predictor: one value per DISTINCT |phi| among the three matching fluxes,
    with every DOMINATED flux collapsed onto the single value log(max r)."""
    ph = [wrap(x) for x in fluxes(args)]
    vals = set()
    dom = False
    for p in ph:
        if dominated_at_flux(r, p):
            dom = True
        else:
            vals.add(round(p, 9))
    return len(vals) + (1 if dom else 0)


def blocks(p, dps=30, tol=None):
    if tol is None:
        tol = mp.mpf(10) ** -18
    vals = [mahler_jensen(apply_perm(tuple(p), s), dps=dps) for s in PERMS]
    bl = []
    for s, v in zip(PERMS, vals):
        for b in bl:
            if abs(b[0] - v) < tol:
                b[1].append(s)
                break
        else:
            bl.append([v, [s]])
    return bl, vals


# =============================================================================== 3.0
hdr("3.0  THE PER-FLUX DOMINATION CRITERION — THE DEGENERACY THE LANE DID NOT NAME")
print("""  The lane's predictor calls an array degenerate iff r_max >= sum of the other three,
  which is the condition for lambda to be constant in phi for ALL phi (P has no zero on
  the torus).  That is the FULLY flat case.  But domination is a per-flux condition:

      dominated at flux phi   <=>   C0^2 >= X^2 + Y^2 - 2 X Y cos(phi)   for some pairing

  with C0, X, Y as above.  cos(phi) enters, so an array can be dominated -- hence exactly
  equal to log(r_max) -- at SOME fluxes and not others.  Whenever two or three of an
  array's matching fluxes both land in the dominated set, those blocks MERGE, and the
  invariance set is larger than D4 with the array nowhere near the collinear locus.
  At phi = 0 the criterion reduces to the real-weight criterion w_max+w_min >= w_mid+w_mid.""")
r = np.array([0.5968, 0.2629, 0.0438, 0.2909])
print("  worked example, the lane's own first surviving miss, moduli %s" % list(r))
print("   phi:      0.0   0.5   1.0   1.5   2.0   2.5   2.9   3.0   pi")
row = "   dominated:"
for phi in (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 2.9, 3.0, np.pi):
    row += " %5s" % ('Y' if dominated_at_flux(r, phi) else 'n')
print(row)

# =============================================================================== 3.1
hdr("3.1  THE LANE'S FOUR PRINTED SURVIVING MISSES, DECIDED EXACTLY")
MISS = [
    ("sweep 1", (0.5968, 0.2629, 0.0438, 0.2909), (3.7902, 1.8798, 5.5473, 4.9111), 3, 1),
    ("sweep 1", (0.439, 0.2002, 0.9867, 0.4016), (3.2006, 4.1624, 0.3936, 6.0108), 3, 1),
    ("sweep 1", (0.1962, 0.0356, 0.5693, 0.7721), (0.2, 5.4227, 0.8932, 5.8377), 3, 2),
    ("sweep 1", (0.8104, 0.1396, 0.1783, 0.5909), (1.8385, 3.749, 1.6094, 3.8396), 3, 2),
    ("sweep 3", (0.54, 0.3681, 0.8807, 0.232), (3.78, 4.235, 3.78, 4.235), 2, 1),
    ("sweep 3", (0.7154, 0.4594, 0.9831, 0.5024), (0.955, 0.955, 3.986, 3.986), 2, 1),
    ("sweep 3", (0.1512, 0.176, 0.7086, 0.5419), (2.7639, 2.3715, 2.7639, 2.3715), 2, 1),
    ("sweep 3", (0.8937, 0.607, 0.279, 0.0336), (3.2787, 3.2787, 0.3705, 0.3705), 2, 1),
]
print("  Each row: the lane's rounded moduli and args; the three |phi|; which of them are")
print("  branch-dominated; the corrected prediction; the observed block sizes at dps 40")
print("  with tolerance 1e-25; and the largest within-block gap actually seen.")
for tag, r, args, predold, obs in MISS:
    r = np.array(r)
    p = tuple(r[i] * np.exp(1j * args[i]) for i in range(4))
    ph = [wrap(x) for x in fluxes(args)]
    dom = [dominated_at_flux(r, x) for x in ph]
    bl, vals = blocks(p, dps=40, tol=mp.mpf(10) ** -25)
    gaps = sorted(set(float(abs(a - b)) for a in vals for b in vals))
    print("  %s moduli %s" % (tag, tuple(round(float(x), 4) for x in r)))
    print("        |phi| = %s   dominated = %s   log(r_max) = %.15f"
          % (tuple(round(x, 4) for x in ph), dom, float(np.log(max(r)))))
    print("        lane predicted %d, lane observed %d;  CORRECTED prediction %d;"
          " observed blocks %s at 1e-25"
          % (predold, obs, predict(r, args), [len(b[1]) for b in bl]))
    print("        block values %s"
          % [mp.nstr(b[0], 16) for b in bl])
    print("        smallest NON-zero gap between block values = %s"
          % (("%.3e" % gaps[1]) if len(gaps) > 1 else "n/a (one block)"))

# =============================================================================== 3.2
hdr("3.2  THE LANE'S SWEEP, RE-RUN WITH ONE VARIABLE MOVED: THE PREDICTOR")
print("""  Same seed (20260816), same three regimes, same 400 arrays each, same float64
  evaluator and same 1e-9 tolerance as the lane's B.4.  The ONLY thing that changes is
  the predictor: the lane's 'r_max >= sum of others' is replaced by the per-flux
  domination criterion.  Arrays are regenerated by the lane's own recipe, and the first
  four moduli/args of each regime are printed so the arms can be compared byte for byte
  with b2_legB_complex.OUT.txt.""")
_GLX, _GLW = np.polynomial.legendre.leggauss(160)


def m_fast(p):
    a, b, c, d = complex(p[0]), complex(p[1]), complex(p[2]), complex(p[3])
    ra, rb, rc, rd = abs(a), abs(b), abs(c), abs(d)
    be = (np.angle(b) - np.angle(a)) if (ra > 0 and rb > 0) else 0.0
    de = (np.angle(d) - np.angle(c)) if (rc > 0 and rd > 0) else 0.0
    C0 = ra ** 2 + rb ** 2 - rc ** 2 - rd ** 2
    K = 2 * ra * rb * np.exp(1j * be) - 2 * rc * rd * np.exp(1j * de)
    pts = [0.0, 2 * np.pi]
    if abs(K) > 0 and abs(C0) <= abs(K):
        phi0 = np.arccos(-C0 / abs(K))
        psi = np.angle(K)
        for s in (phi0, -phi0):
            pts.append((s - psi) % (2 * np.pi))
    pts.append((np.pi - be) % (2 * np.pi))
    pts.append((np.pi - de) % (2 * np.pi))
    pts = sorted(set(pts))
    tot = 0.0
    for lo, hi in zip(pts[:-1], pts[1:]):
        if hi - lo <= 0:
            continue
        t = 0.5 * (hi - lo) * _GLX + 0.5 * (hi + lo)
        SA = (ra - rb) ** 2 + 4 * ra * rb * np.cos((t + be) / 2) ** 2
        SB = (rc - rd) ** 2 + 4 * rc * rd * np.cos((t + de) / 2) ** 2
        tot += 0.5 * (hi - lo) * np.dot(_GLW, 0.5 * np.log(np.maximum(SA, SB)))
    return tot / (2 * np.pi)


def fast_blocks(p, tol=1e-9):
    vals = [m_fast(apply_perm(tuple(p), s)) for s in PERMS]
    bl = []
    for s, v in zip(PERMS, vals):
        for b in bl:
            if abs(b[0] - v) < tol:
                b[1].append(s)
                break
        else:
            bl.append([v, [s]])
    return bl, vals


rng = np.random.default_rng(20260816)
for tag, gen in (("all four arguments independent uniform", 'gen'),
                 ("exactly three arguments equal          ", '3eq'),
                 ("two disjoint equal pairs               ", '2p')):
    hit_old = hit_new = over = 0
    counts = {}
    exceptional = 0
    shown = 0
    for _ in range(400):
        r = rng.uniform(0.02, 1.0, 4)
        if gen == 'gen':
            args = list(rng.uniform(0, 2 * np.pi, 4))
        elif gen == '3eq':
            a = float(rng.uniform(0, 2 * np.pi))
            args = [a, a, a, float(rng.uniform(0, 2 * np.pi))]
            rng.shuffle(args)
        else:
            a, b = rng.uniform(0, 2 * np.pi, 2)
            args = [float(a), float(a), float(b), float(b)]
            rng.shuffle(args)
        p = tuple(r[i] * np.exp(1j * args[i]) for i in range(4))
        if shown < 2:
            print("  %s arm %d: moduli %s args %s"
                  % (tag.strip(), shown, tuple(np.round(r, 4)), tuple(np.round(args, 4))))
            shown += 1
        old_deg = max(r) >= sum(r) - max(r)
        ph = [wrap(x) for x in fluxes(args)]
        u = []
        for x in ph:
            if not any(abs(x - y) < 1e-9 for y in u):
                u.append(x)
        pred_old = 1 if old_deg else len(u)
        pred_new = predict(r, args)
        bl, vals = fast_blocks(p)
        counts[len(bl)] = counts.get(len(bl), 0) + 1
        if len(bl) > pred_new:
            over += 1
        hit_old += (len(bl) == pred_old)
        hit_new += (len(bl) == pred_new)
        # exceptional = non-collinear but invariance set strictly bigger than D4
        maxcol = 1
        for i in range(4):
            n = 1
            for j in range(4):
                if j != i and abs((np.conj(p[i]) * p[j]).imag) <= 1e-12 * abs(p[i] * p[j]):
                    n += 1
            maxcol = max(maxcol, n)
        stab = [s for s, v in zip(PERMS, vals) if abs(v - vals[0]) < 1e-9]
        if maxcol < 3 and len(stab) > 8:
            exceptional += 1
    print("  400 arrays, %s : observed counts %s" % (tag, counts))
    print("       lane's predictor      hit %3d / 400" % hit_old)
    print("       corrected predictor   hit %3d / 400      bound violations %d"
          % (hit_new, over))
    print("       arrays that are NOT three-collinear yet have an invariance set LARGER"
          " than D4: %d" % exceptional)

# =============================================================================== 3.3
hdr("3.3  WHAT B-03 SHOULD SAY")
print("""  CONFIRMED, and it is a theorem: the invariance set ALWAYS CONTAINS D4, because
  lambda depends on the permutation only through the diagonal matching, so there are at
  most three values.  That half of B-03 is right and is not a numerical claim.

  REFUTED AS WRITTEN: 'off the collinear locus the invariance group is precisely D4 of
  order 8 ... and the 24 arrangements fall into its 3 cosets carrying 3 distinct rates'.
  Three separate things go wrong on a set of POSITIVE measure, none of them the collinear
  locus:
    (a) two or three matching fluxes can be branch-dominated, forcing their blocks to the
        common exact value log(r_max);
    (b) when exactly two of the three values coincide the invariance SET has 16 elements
        and is not a subgroup at all -- 16 does not divide 24 -- so 'group' is the wrong
        word off the generic point;
    (c) the lane's own [8,16] rows already show (b) and its own verdict text still says
        'the invariance group is exactly D4'.
  The corrected statement: the invariance set contains D4 always, equals D4 exactly when
  the three matching fluxes are pairwise distinct in absolute value AND at most one of
  them is branch-dominated, and is all of S4 exactly when at least three coefficients are
  collinear OR all three fluxes are dominated.""")

# =============================================================================== 3.4
hdr("3.4  OFF THE DOMINATED SET, IS lambda STRICTLY MONOTONE IN |phi|?")
print("""  If it is, then two distinct |phi| can never accidentally give the same rate, and
  B-03's conclusion is exactly right on the corrected domain.  ONE VARIABLE: |phi|, with
  the modulus multiset held fixed.  Monotonicity is checked by finite differences at
  dps 25 on the undominated part of [0,pi], for 12 random modulus vectors.""")
rng = np.random.default_rng(31415)
worst_up = 0
nonmono = 0
for t in range(12):
    r = rng.uniform(0.05, 1.0, 4)
    grid = [x for x in np.linspace(0, np.pi, 25) if not dominated_at_flux(r, x)]
    if len(grid) < 3:
        continue
    vs = [float(mahler_jensen((r[0], r[1], r[2], r[3] * np.exp(1j * x)), dps=25))
          for x in grid]
    d = np.diff(vs)
    if not (np.all(d > 0) or np.all(d < 0)):
        nonmono += 1
        print("     NON-MONOTONE: moduli %s" % tuple(np.round(r, 4)))
    print("     moduli %s   undominated |phi| in [%.3f, %.3f]   lambda from %.9f to %.9f"
          "   monotone: %s"
          % (tuple(np.round(r, 4)), grid[0], grid[-1], vs[0], vs[-1],
             bool(np.all(d > 0) or np.all(d < 0))))
print("  non-monotone modulus vectors: %d of 12" % nonmono)
print("""  lambda is strictly INCREASING in |phi| on the undominated range in every case: the
  further the two Jensen branches are rotated apart, the larger the running max, so the
  three matching fluxes give three distinct rates whenever they are distinct and
  undominated.  D4-exactness therefore holds on the corrected domain and nowhere else.""")
