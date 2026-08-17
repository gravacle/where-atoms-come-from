"""LANE W-13 / Z  --  shared library.  No figure is produced here; z1..z4 import this."""
import numpy as np
from fractions import Fraction as F
import math

# ------------------------------------------------------------------ exact stratum classifier
def strat_exact(p):
    """p = (p00,p10,p01,p11) as Fractions summing to 1, all >= 0.
    Returns (TYPE, detail).  TYPE in {'EMPTY','ONE','TWO','CURVE'}.
    THEOREM Z1 of this lane.  Pure integer/rational comparisons; no float."""
    p00, p10, p01, p11 = p
    assert sum(p) == 1 and all(q >= 0 for q in p)
    S1, S2 = p00 + p10, p01 + p11
    D1, D2 = abs(p00 - p10), abs(p01 - p11)
    # --- the three curve strata, each a pair of exact equalities
    cI   = (p00 == p10) and (p01 == p11)     # Z = {x=-1} x T   (A,B share the zero x=-1)
    cII  = (p00 == p01) and (p10 == p11)     # Z = T x {y=-1}   (P = (1+y)(p00+p10 x))
    cIII = (p00 == p11) and (p10 == p01)     # Z = graph over the x-circle
    if cI or cII or cIII:
        return 'CURVE', '+'.join(n for n, b in (('I', cI), ('II', cII), ('III', cIII)) if b)
    prod = (S1 - S2) * (D1 - D2)
    if prod > 0:
        return 'EMPTY', ''
    if prod < 0:
        return 'TWO', ''
    return 'ONE', ('S1=S2' if S1 == S2 else 'D1=D2')

def polygon_exact(p):
    """W-01's convex-hull criterion under the FREE-CHARACTER reading: max_a p_a <= 1/2."""
    return max(p) <= F(1, 2)

# ------------------------------------------------------------------ closed-form zero locations
def zeros_closed_form(p):
    """Return the list of isolated zeros (x0,y0) in C^2, computed from the closed form
    cos(s) = -C/D with C = |A(1)|^2-... .  Exact where p is rational: cos s0 is RATIONAL.
    Returns [] for EMPTY, [pt] for ONE, [pt,conj(pt)] for TWO, None for CURVE."""
    p00, p10, p01, p11 = [F(q) for q in p]
    typ, _ = strat_exact((p00, p10, p01, p11))
    if typ == 'CURVE':
        return None
    if typ == 'EMPTY':
        return []
    C = p00**2 + p10**2 - p01**2 - p11**2
    D = 2 * (p00 * p10 - p01 * p11)
    if D == 0:
        return []                       # C != 0 here (else CURVE); no crossing
    cs = F(-C, 1) / D                   # EXACT rational cos(s0)
    if cs > 1 or cs < -1:
        return []
    sn2 = 1 - cs * cs                   # exact rational sin^2
    sn = math.sqrt(float(sn2))
    out = []
    for sgn in (+1, -1):
        x0 = complex(float(cs), sgn * sn)
        A = float(p00) + float(p10) * x0
        B = float(p01) + float(p11) * x0
        y0 = -A / B if abs(B) > 0 else None
        out.append((x0, y0, cs, sn2, sgn))
        if sn2 == 0:
            break
    return out

# ------------------------------------------------------------------ Jensen-reduced evaluators
def jensen_mods(p, t):
    """|A(e^{it})|, |B(e^{it})| as float arrays."""
    p00, p10, p01, p11 = [float(q) for q in p]
    c, s = np.cos(t), np.sin(t)
    a = np.hypot(p00 + p10 * c, p10 * s)
    b = np.hypot(p01 + p11 * c, p11 * s)
    return a, b

def mahler(p, n=1 << 20):
    """m(P) by the Jensen reduction: (1/2pi) INT log max(|A|,|B|).  CONTINUOUS integrand.
    Midpoint rule on n points."""
    t = (np.arange(n) + 0.5) * (2 * np.pi / n)
    a, b = jensen_mods(p, t)
    return float(np.mean(np.log(np.maximum(a, b))))

def min_abs_P(p, n=1 << 20):
    """EXACT-in-the-limit min over T^2 of |P| = min_t ||A|-|B||.  One-dimensional; the
    minimum of a CONTINUOUS function, so a fine grid is a genuine upper bound that converges."""
    t = (np.arange(n) + 0.5) * (2 * np.pi / n)
    a, b = jensen_mods(p, t)
    return float(np.min(np.abs(a - b)))

def _arcfrac(p, t, eps):
    a, b = jensen_mods(p, t)
    den = 2.0 * a * b
    with np.errstate(divide='ignore', invalid='ignore'):
        R = np.where(den > 0, (eps ** 2 - a ** 2 - b ** 2) / np.where(den > 0, den, 1.0), -np.inf)
    R = np.clip(R, -1.0, 1.0)
    frac = 1.0 - np.arccos(R) / np.pi
    deg = den <= 0
    if np.any(deg):
        other = np.where(a > b, a, b)
        frac = np.where(deg, (other < eps).astype(float), frac)
    return frac

def _crossings(p):
    """t-values in [0,2pi) where |A(t)| = |B(t)|.  Closed form.  [] if none, None if |A|==|B|."""
    p00, p10, p01, p11 = [F(q) for q in p]
    C = p00**2 + p10**2 - p01**2 - p11**2
    D = 2 * (p00 * p10 - p01 * p11)
    if D == 0:
        return None if C == 0 else []
    z = F(-C, 1) / D
    if z > 1 or z < -1:
        return []
    t0 = math.acos(float(z))
    return [t0] if (z == 1 or z == -1) else [t0, 2 * math.pi - t0]

def sublevel_measure(p, eps, n=1 << 16, n_fine=1 << 15):
    """Haar measure of { (x,y) in T^2 : |P| < eps }, by the EXACT y-arc formula, on a mesh
    REFINED around the Jensen-branch crossings (whose locations are known in closed form).
    A uniform mesh underestimates mu once eps drops below the mesh spacing -- that is a
    window artefact, and this routine removes it rather than reporting it."""
    xs = _crossings(p)
    nodes = [np.linspace(0.0, 2 * np.pi, n, endpoint=False)]
    if xs:
        for t0 in xs:
            w = 1e-14
            while w < np.pi:
                h = abs(np.subtract(*jensen_mods(p, np.array([t0 + w]))))[0]
                if h > eps:
                    break
                w *= 2.0
            w = min(w * 4.0, np.pi)
            loc = np.linspace(-w, w, n_fine)
            nodes.append((t0 + loc) % (2 * np.pi))
    t = np.unique(np.concatenate(nodes))
    t = np.concatenate([t, [t[0] + 2 * np.pi]])
    f = _arcfrac(p, t, eps)
    return float(np.trapz(f, t) / (2 * np.pi))

def sublevel_measure_uniform(p, eps, n=1 << 20):
    """The unrefined version, kept so the window artefact can be shown rather than hidden."""
    t = (np.arange(n) + 0.5) * (2 * np.pi / n)
    a, b = jensen_mods(p, t)
    den = 2.0 * a * b
    with np.errstate(divide='ignore', invalid='ignore'):
        R = np.where(den > 0, (eps ** 2 - a ** 2 - b ** 2) / np.where(den > 0, den, 1.0), -np.inf)
    R = np.clip(R, -1.0, 1.0)
    frac = 1.0 - np.arccos(R) / np.pi          # fraction of the y-circle with |P| < eps
    deg = den <= 0                              # a or b vanishes: |P| = |the other|
    if np.any(deg):
        other = np.where(a > b, a, b)
        frac = np.where(deg, (other < eps).astype(float), frac)
    return float(np.mean(frac))

# ------------------------------------------------------------------ local data at a zero
def local_alpha_beta(p, x0, y0):
    """alpha = x dP/dx, beta = y dP/dy at the zero.  P ~ i(alpha*sigma + beta*tau)."""
    p00, p10, p01, p11 = [float(q) for q in p]
    alpha = x0 * (p10 + p11 * y0)
    beta = y0 * (p01 + p11 * x0)
    return alpha, beta

def local_singvals(alpha, beta):
    M = np.array([[alpha.real, beta.real], [alpha.imag, beta.imag]])
    sv = np.linalg.svd(M, compute_uv=False)
    return float(sv[0]), float(sv[1]), float(np.linalg.det(M))

def Pval(p, x, y):
    p00, p10, p01, p11 = [float(q) for q in p]
    return p00 + p10 * x + p01 * y + p11 * x * y

# ------------------------------------------------------------------ named weight vectors
def fr(*a):
    return tuple(F(q) for q in a)

NAMED = [
    ("K1_REG   S3/S4 sense-C on K1 (0,.3,.3,.4)  <- N1 AS REGISTERED",
     fr(F(0), F(3, 10), F(3, 10), F(2, 5))),
    ("S1_PUB   K1's published ready state (1/2,0,0,1/4,1/4)",
     fr(F(0), F(0), F(1, 2), F(1, 2))),
    ("SENSEC4  S4 sense-C, four classes (1/4,1/4,1/4,1/4)",
     fr(F(1, 4), F(1, 4), F(1, 4), F(1, 4))),
    ("B0b_U    S4:575 ring torus loops meet, uniform  {00:4,10:2,01:2,11:1}/9",
     fr(F(4, 9), F(2, 9), F(2, 9), F(1, 9))),
    ("B0a_U    S4:575 ring torus loops disjoint       {00:2,10:4,01:3}/9",
     fr(F(2, 9), F(4, 9), F(3, 9), F(0))),
    ("B4_U     S4:575 spindle, uniform                {00:1,10:1,01:1,11:3}/6",
     fr(F(1, 6), F(1, 6), F(1, 6), F(1, 2))),
    ("B1_U     S4:575 K1 (=B3,B2), uniform           {10:2,01:2,11:1}/5",
     fr(F(0), F(2, 5), F(2, 5), F(1, 5))),
    ("B1p_U    S4:575 K1-bridged, uniform            {10:3,01:3}/6",
     fr(F(0), F(1, 2), F(1, 2), F(0))),
    ("B1q_U    S4:575 K1-bridged + spectator         {00:1,10:3,01:3}/7",
     fr(F(1, 7), F(3, 7), F(3, 7), F(0))),
    ("B1s_U    S4:575 K1 subdivided                  {10:5,01:5,11:1}/11",
     fr(F(0), F(5, 11), F(5, 11), F(1, 11))),
    ("TANGENT  constructed: (1/10,1/5,3/10,2/5), D1 = D2",
     fr(F(1, 10), F(1, 5), F(3, 10), F(2, 5))),
    ("TWO4     constructed: (3/20,1/4,3/10,3/10), all four occupied",
     fr(F(3, 20), F(1, 4), F(3, 10), F(3, 10))),
    ("CURVE3   constructed: (3/10,1/5,1/5,3/10), stratum III",
     fr(F(3, 10), F(1, 5), F(1, 5), F(3, 10))),
    ("NOZERO4  constructed: (9/20,1/4,1/5,1/10), polygon holds, no zero",
     fr(F(9, 20), F(1, 4), F(1, 5), F(1, 10))),
]
