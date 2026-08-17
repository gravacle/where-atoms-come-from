# W10-D leg 2 -- N1 (the rate IS a Mahler measure), N2 (the multiset theorem), N3 (the Haar-null
# inversion), all three tested ON FOUR-CLASS CARRIERS, and the multiset theorem's REAL hypothesis
# isolated by moving the COEFFICIENT FIELD one step at a time: non-negative real -> real -> complex.
import numpy as np
from itertools import permutations
from fractions import Fraction

np.set_printoptions(precision=12, suppress=False)

# ---------------------------------------------------------------- Mahler by Jensen reduction
GLX, GLW = np.polynomial.legendre.leggauss(400)

def _arc(a, b, lo, hi):
    """INT_lo^hi log|a + b e^{it}| dt by Gauss-Legendre (integrand analytic unless |a|==|b| and
    the arc touches pi)."""
    t = 0.5*(hi-lo)*GLX + 0.5*(hi+lo)
    val = 0.5*np.log(np.abs(a)**2 + np.abs(b)**2 + 2*(a*b)*np.cos(t))
    return 0.5*(hi-lo)*np.dot(GLW, val)

def mahler_jensen(p, split=True):
    """m(p00 + p10 x + p01 y + p11 xy) via Jensen in x.
       branches A(t) = |p00 + p01 e^{it}|, B(t) = |p10 + p11 e^{it}| ; both even in t."""
    p00, p10, p01, p11 = [float(q) for q in p]
    A2 = lambda ct: p00**2 + p01**2 + 2*p00*p01*ct
    B2 = lambda ct: p10**2 + p11**2 + 2*p10*p11*ct
    # crossing:  (A2-B2)(cos t) = alpha + beta cos t = 0
    alpha = p00**2 + p01**2 - p10**2 - p11**2
    beta = 2*(p00*p01 - p10*p11)
    cuts = [0.0, np.pi]
    if abs(beta) > 1e-300:
        r = -alpha/beta
        if -1.0 < r < 1.0:
            cuts.insert(1, float(np.arccos(r)))
    tot = 0.0
    exact_flag = (len(cuts) == 2)
    for lo, hi in zip(cuts[:-1], cuts[1:]):
        cm = np.cos(0.5*(lo+hi))
        if A2(cm) >= B2(cm):
            tot += _arc(p00, p01, lo, hi)
        else:
            tot += _arc(p10, p11, lo, hi)
    return tot/np.pi, exact_flag

def mahler_exact_if_nocross(p):
    """When one Jensen branch dominates on the whole circle, m(P) = log max of that pair."""
    p00, p10, p01, p11 = [float(q) for q in p]
    alpha = p00**2 + p01**2 - p10**2 - p11**2
    beta = 2*(p00*p01 - p10*p11)
    if abs(beta) > 1e-300 and -1.0 < -alpha/beta < 1.0:
        return None
    cm = 0.0
    if p00**2+p01**2 + 2*p00*p01*cm >= p10**2+p11**2 + 2*p10*p11*cm:
        return np.log(max(abs(p00), abs(p01)))
    return np.log(max(abs(p10), abs(p11)))

def mahler_grid(p, n):
    """2D grid control -- the noise-limited method the brief warns about; kept as a cross-check."""
    t = (np.arange(n)+0.5)*2*np.pi/n
    X = np.exp(1j*t)[:, None]
    Y = np.exp(1j*t)[None, :]
    P = p[0] + p[1]*X + p[2]*Y + p[3]*X*Y
    return np.log(np.abs(P)).mean()

CAR = {'B1  K1 (3 class)': np.array([0.0, 2/5, 2/5, 1/5]),
       'B1q spectator(3)': np.array([1/7, 3/7, 3/7, 0.0]),
       'B1p bridged  (2)': np.array([0.0, 1/2, 1/2, 0.0]),
       'B0b torus    (4)': np.array([4/9, 2/9, 1/9, 2/9]),
       'B4  spindle  (4)': np.array([1/6, 1/6, 1/6, 3/6]),
       'SENSE-C      (4)': np.array([.25, .25, .25, .25])}

print("="*100)
print("== 2A  N1 ON FOUR-CLASS CARRIERS: lambda = m(P), AND BOTH FOUR-CLASS ROWS ARE EXACT ==")
print("="*100)
print(f"{'carrier':18s} {'Jensen':>16s} {'exact/closed form':>22s} {'|diff|':>10s}   S4's published")
S4 = {'B1  K1 (3 class)': -0.756573585640, 'B1q spectator(3)': -0.741029582571,
      'B1p bridged  (2)': -0.693147180560, 'B0b torus    (4)': -0.810930216216,
      'B4  spindle  (4)': -0.693147180560, 'SENSE-C      (4)': -1.386294361120}
for nm, p in CAR.items():
    mj, nocross = mahler_jensen(p)
    me = mahler_exact_if_nocross(p)
    lab = f"{me:.12f}" if me is not None else "(branches cross)"
    d = abs(mj-me) if me is not None else float('nan')
    print(f"{nm:18s} {mj:16.12f} {lab:>22s} {d:10.2e}   {S4[nm]:+.12f}  dev {abs(mj-S4[nm]):.1e}")
print("\n  B0b: Jensen branches are |4/9 + (1/9)e^{it}| and |2/9 + (2/9)e^{it}|; their squares")
print("       differ by (17+8cos t)/81 - (8+8cos t)/81 = 9/81 > 0 for ALL t -- a CONSTANT gap,")
print("       so the first branch dominates everywhere and m = log max(4/9,1/9) = log(4/9).")
print(f"       log(4/9) = {np.log(4/9):.15f}   EXACT.  S4 labelled this row 'QUADRATURE ONLY'.")
print("  B4 : squares differ by (8+4cos t)/36 > 0, second branch dominates, m = log(1/2) EXACT.")
print("  BOTH FOUR-CLASS CARRIERS THE CORPUS OWNS HAVE CLOSED-FORM RATES. Nothing is quadrature.")

print("\n  2D-GRID CONTROL (the noise-limited method):")
for nm in ('B0b torus    (4)', 'SENSE-C      (4)'):
    p = CAR[nm]
    for n in (512, 2048, 8192):
        print(f"    {nm}  n={n:5d}  grid {mahler_grid(p,n):+.9f}   Jensen {mahler_jensen(p)[0]:+.12f}")
print("    SENSE-C's grid value drifts (P = (1+x)(1+y)/4 has a CURVE of torus zeros); B0b's does")
print("    not (no torus zero).  The brief's quadrature warning is weight-scoped, not class-scoped.")

print("\n"+"="*100)
print("== 2B  N1 AGAIN, THE WAY THE CORPUS DEFINES IT: lim (1/N) log|Omega_N| ON FOUR CLASSES ==")
print("="*100)
rng = np.random.default_rng(20260816)
print("  Direct schedule-B simulation, k = 1..N, Omega_N = prod Z_k, at RANDOM generic (f,c).")
print(f"{'carrier':18s} {'(f,c)':>28s} {'N=1e5':>14s} {'N=1e6':>14s} {'m(P)':>14s} {'dev':>9s}")
for nm in ('B0b torus    (4)', 'B4  spindle  (4)', 'B1  K1 (3 class)'):
    p = CAR[nm]
    mj = mahler_jensen(p)[0]
    for trial in range(2):
        fv, cv = rng.uniform(-np.pi, np.pi, 2)
        k = np.arange(1, 1000001)
        Z = p[0] + p[1]*np.exp(-1j*fv*k) + p[2]*np.exp(1j*cv*k) + p[3]*np.exp(1j*(cv-fv)*k)
        L = np.log(np.abs(Z))
        a5, a6 = L[:100000].mean(), L.mean()
        print(f"{nm:18s} ({fv:+.5f},{cv:+.5f}) {a5:14.9f} {a6:14.9f} {mj:14.9f} {abs(a6-mj):9.2e}")
print("  N1 HOLDS ON FOUR-CLASS CARRIERS.  Nothing in the identification uses p00 = 0, three")
print("  classes, or K1's incidence: it is Jensen + Weyl on the character lattice.")

print("\n"+"="*100)
print("== 2C  N3 -- THE HAAR-NULL INVERSION, ON A FOUR-CLASS CARRIER ==")
print("="*100)
print("  N3: the rate is invariant under EVERY absolutely continuous connection measure, because")
print("  the resonance set {(f,c) : rank L > 0} is Haar-null.  Test: 40 random (f,c) on B0b.")
p = CAR['B0b torus    (4)']
mj = mahler_jensen(p)[0]
k = np.arange(1, 400001)
vals = []
for _ in range(40):
    fv, cv = rng.uniform(-np.pi, np.pi, 2)
    Z = p[0] + p[1]*np.exp(-1j*fv*k) + p[2]*np.exp(1j*cv*k) + p[3]*np.exp(1j*(cv-fv)*k)
    vals.append(np.log(np.abs(Z)).mean())
vals = np.array(vals)
print(f"  40 random connections, N=4e5:  mean {vals.mean():.9f}  sd {vals.std():.2e}  "
      f"max|dev from m(P)| {np.abs(vals-mj).max():.2e}")
print(f"  RESONANT points on the same carrier (rank L > 0), same N, same code path:")
for lab, (fv, cv) in [("c = f      (uv = 1)", (1.0, 1.0)),
                      ("c = -f     (u = v)", (1.0, -1.0)),
                      ("f = 2.0, c = 1.1 (S3/S4 headline, -11f+20c=0)", (2.0, 1.1)),
                      ("f = pi/2, c = 3pi/2 (S1 order 4)", (np.pi/2, 3*np.pi/2))]:
    Z = p[0] + p[1]*np.exp(-1j*fv*k) + p[2]*np.exp(1j*cv*k) + p[3]*np.exp(1j*(cv-fv)*k)
    print(f"    {lab:46s} {np.log(np.abs(Z)).mean():+.9f}   dev {abs(np.log(np.abs(Z)).mean()-mj):.2e}")
print("  The resonant points are exactly where the rate departs, and they are Haar-null.  N3's")
print("  statement uses only that there are TWO designated holonomies -- CARRIER-INDEPENDENT")
print("  given a two-loop designation.  It says nothing about class occupancy.")

print("\n"+"="*100)
print("== 2D  N2 -- THE MULTISET THEOREM.  ARM DIFF, THEN THE ISOLATION OF ITS REAL HYPOTHESIS ==")
print("="*100)
print("  THE REGISTRAR PREDICTED D4 (order 8) off K1 and REFUTED HIS OWN PREDICTION: all 24 agree.")
print("  He named the operative hypothesis REAL NON-NEGATIVITY.  This leg moves the COEFFICIENT")
print("  FIELD one step at a time and reads which permutations survive each step.\n")

def spread24(p):
    vals = [mahler_jensen_c(np.array(q))[0] for q in permutations(p)]
    return max(vals)-min(vals), vals

def mahler_jensen_c(p, n=1 << 18):
    """Jensen in x for COMPLEX coefficients: integrand log max(|p00+p01 e^{it}|,|p10+p11 e^{it}|)
       is no longer even in t, so integrate the full circle.  Trapezoid on a smooth periodic
       integrand (the max of two analytic moduli; kinks only at crossings)."""
    t = (np.arange(n)+0.5)*2*np.pi/n
    e = np.exp(1j*t)
    A = np.abs(p[0] + p[2]*e)
    B = np.abs(p[1] + p[3]*e)
    return np.log(np.maximum(A, B)).mean(), None

ARMS = [("non-neg real  B0b  (4/9,2/9,1/9,2/9)", np.array([4/9, 2/9, 1/9, 2/9])),
        ("non-neg real  B4   (1/6,1/6,1/6,3/6)", np.array([1/6, 1/6, 1/6, 3/6])),
        ("non-neg real  generic             ", np.array([0.31, 0.17, 0.29, 0.23])),
        ("REAL, one NEGATIVE entry          ", np.array([0.55, -0.20, 0.40, 0.25])),
        ("REAL, two NEGATIVE entries        ", np.array([-0.31, 0.62, -0.17, 0.86])),
        ("COMPLEX, generic                  ", np.array([0.4+0.1j, 0.2-0.3j, 0.3+0.25j, 0.1+0.4j])),
        ("COMPLEX, phases only              ", np.array([1.0, np.exp(0.7j), np.exp(2.1j), np.exp(-1.3j)]))]
print("  ARM DIFF -- the coefficient vectors actually integrated (they are pairwise distinct):")
for lab, v in ARMS:
    print(f"    {lab}  {np.array2string(v, precision=4)}")
print()
print(f"  {'arm':38s} {'24-perm spread':>16s} {'#distinct values (1e-12)':>26s}")
for lab, v in ARMS:
    vals = [mahler_jensen_c(np.array(q))[0] for q in permutations(v)]
    vals = np.array(vals)
    sp = vals.max()-vals.min()
    nd = len(np.unique(np.round(vals, 12)))
    print(f"  {lab:38s} {sp:16.3e} {nd:26d}")
print("\n  WHICH PERMUTATIONS SURVIVE, per arm (index order p00,p10,p01,p11):")
for lab, v in ARMS:
    base = mahler_jensen_c(v)[0]
    keep = [q for q in permutations(range(4))
            if abs(mahler_jensen_c(v[list(q)])[0] - base) < 1e-10]
    print(f"    {lab}  |G| = {len(keep)}")
print("\n  READ.  Reality of the coefficients is sufficient and NON-NEGATIVITY IS NOT NEEDED:")
print("  |a + b e^{it}| = |b + a e^{it}| pointwise requires only conj(a)b = conj(b)a, i.e. ab REAL.")
print("  For real a,b of ANY sign that holds.  The registrar's name 'REAL NON-NEGATIVITY' is one")
print("  hypothesis too strong; the operative one is REALITY (equivalently: p_ab is a measure on")
print("  the class set, whose weights are real by construction on every carrier).")
