"""
W-10 SYNTHESIS LANE (register head W-09, successor registrar).  INDEPENDENT VERIFICATION PASS.

PURPOSE.  I am the synthesis and I do not want to be a pure aggregator of eight refuter
verdicts.  This script re-decides, from scratch and with my own code, the five points on which
the eight verdicts either COLLIDE with one another or carry the whole scope ruling.

CONVENTIONS (published, per the round's process rules).
  characters:  class 00 -> 1 ; 10 -> u = conj(W_F) = e^{-i f} ; 01 -> v = W_C = e^{i c} ; 11 -> uv
  P(x,y)     = p00 + p10 x + p01 y + p11 xy ;  Z_k = P(u^k, v^k) ;  lambda = m(P)
  Jensen-in-y reduction (the round's mandated evaluator, no 2D grid):
     m(P) = (1/2pi) INT_0^{2pi} log max( |p00 + p10 e^{it}| , |p01 + p11 e^{it}| ) dt
  DOUBLE PRECISION unless a line says mpmath; every precision-sensitive figure below is either
  closed-form, exact in Fractions, or recomputed at mp.dps=40.
  GRID.  Part A: deterministic 1201x1201 midpoint grid on [0,2pi)^2 (NOT random draws -- the
  three prior runs of this figure all used random draws with the same seed convention).
  SEED.  numpy default_rng(20260816) where randomness is used at all (Part D only).

PARTS
  A  firing region, three classes vs four, deterministic grid + closed form  [reproduction]
  B  S4's B4 row: is it determined by S4's own published parameters?         [A-1 vs A-2 collision]
  C  arithmetic typing of every (f,c) the corpus publishes                   [C-2's census, redone]
  D  the multiset theorem's hypothesis: three rival names, all refuted       [A-2 vs D-1 collision]
  E  the character identity, exact, at four/three/two classes                [declared CONTROL]
"""
import numpy as np
from fractions import Fraction as F
from itertools import permutations, combinations
import mpmath as mp

OUT = []
def say(s=""):
    print(s); OUT.append(s)

say("="*100)
say("W-10 SYNTHESIS LANE -- INDEPENDENT VERIFICATION.  numpy %s, mpmath %s" % (np.__version__, mp.__version__))
say("="*100)

# ----------------------------------------------------------------------------------------------
# PART A -- THE FIRING REGION.  Deterministic grid; hull test by MAXIMUM ANGULAR GAP.
# ----------------------------------------------------------------------------------------------
say("")
say("PART A -- FIRING REGION: 0 in conv{occupied characters}?   deterministic 1201x1201 midpoint grid")
say("-"*100)
say("hull test used here: for unit-modulus points, 0 in conv(S) <=> the maximum angular gap between")
say("consecutive points (sorted) is <= pi.  This is a DIFFERENT algorithm from W-09's simplex/LP and")
say("from lane D's two; it is exact up to the grid, with no linear-programming tolerance at all.")

n = 1201
tt = (np.arange(n) + 0.5) * (2*np.pi/n)
fg, cg = np.meshgrid(tt, tt, indexing='ij')

def hull_fires(angle_list):
    """angle_list: list of arrays of angles (same shape).  Returns bool array: 0 in conv?"""
    A = np.stack([np.mod(a, 2*np.pi) for a in angle_list], axis=-1)
    A = np.sort(A, axis=-1)
    gaps = np.diff(A, axis=-1)
    wrap = (A[..., 0] + 2*np.pi) - A[..., -1]
    mx = np.maximum(gaps.max(axis=-1), wrap)
    return mx <= np.pi + 1e-12

# characters: 1 -> angle 0 ; u -> -f ; v -> c ; uv -> c-f
ang_1 = np.zeros_like(fg); ang_u = -fg; ang_v = cg; ang_uv = cg - fg

three = hull_fires([ang_v, ang_u, ang_uv])          # K1: {01,10,11}
three_b1q = hull_fires([ang_1, ang_v, ang_u])       # B1q: {00,01,10}
four = hull_fires([ang_1, ang_v, ang_u, ang_uv])    # B0b / B4: all four
two = hull_fires([ang_v, ang_u])                    # B1p: {01,10}
closed = (np.cos(fg) + np.cos(cg)) <= 0

say("")
say("  occupied set          firing fraction        exact value    |dev|")
for lbl, arr, exact in [("{01,10,11}  K1/B1/B2/B3/B1s", three, 0.25),
                        ("{00,01,10}  B1q            ", three_b1q, 0.25),
                        ("{01,10}     B1p            ", two, 0.0),
                        ("{00,01,10,11} B0b/B4       ", four, 0.5)]:
    fr = arr.mean()
    say("  %s   %.9f          %.4f       %.2e" % (lbl, fr, exact, abs(fr-exact)))
say("")
say("  GRID CONVERGENCE (the boundary of the firing region is measure zero; the error is O(1/n)):")
for nn in (601, 1201, 2401, 4801):
    t2 = (np.arange(nn) + 0.5) * (2*np.pi/nn)
    f2, c2 = np.meshgrid(t2, t2, indexing='ij')
    th = hull_fires([c2, -f2, c2 - f2]).mean()
    fo = hull_fires([np.zeros_like(f2), c2, -f2, c2 - f2]).mean()
    say("     n = %4d   three-class %.9f (dev %.2e)   four-class %.9f (dev %.2e)"
        % (nn, th, abs(th - 0.25), fo, abs(fo - 0.5)))
say("")
say("  closed form cos f + cos c <= 0 agrees with the four-class hull on %d of %d grid points"
    % (int((closed == four).sum()), four.size))
say("  f -> -f flips the verdict:  three classes %d/%d      four classes %d/%d"
    % (int((hull_fires([-ang_v*0+cg, fg, cg+fg]) != three).sum()), three.size,
       int((hull_fires([ang_1, cg, fg, cg+fg]) != four).sum()), four.size))
say("")
say("  RULING: W-09's two exact values reproduce on a deterministic grid and a third hull algorithm.")
say("  NOTE, AND IT IS A VACUITY DISCLOSURE: the f -> -f column could not have come out otherwise at")
say("  four classes -- cos f + cos c is even in f, so the four-class arm is an identity, not a test.")

# ----------------------------------------------------------------------------------------------
# PART B -- IS S4's B4 ROW DETERMINED BY S4's OWN PUBLISHED PARAMETERS?
#           REFUTER A-1 says no (exhibits a second spindle).  REFUTER A-2 says the multiset is
#           FORCED.  They are answering different questions.  I build BOTH complexes.
# ----------------------------------------------------------------------------------------------
say("")
say("PART B -- S4's B4 ROW.  Two non-isomorphic 'two 2-spheres glued at two points' complexes.")
say("-"*100)

def rank_Q(M):
    """exact rank over Q by Fraction elimination"""
    M = [[F(x) for x in row] for row in M]
    rows, cols = len(M), (len(M[0]) if M else 0)
    r = 0
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if M[i][c] != 0: piv = i; break
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        pv = M[r][c]
        for i in range(rows):
            if i != r and M[i][c] != 0:
                fct = M[i][c]/pv
                M[i] = [a - fct*b for a, b in zip(M[i], M[r])]
        r += 1
        if r == rows: break
    return r

def build(verts, edges, faces):
    """edges: list of (src,tgt) vertex names.  faces: list of list of signed edge indices (1-based, sign=orientation)."""
    V, E, Fc = len(verts), len(edges), len(faces)
    vi = {v: i for i, v in enumerate(verts)}
    d1 = [[0]*E for _ in range(V)]
    for j, (s, t) in enumerate(edges):
        d1[vi[s]][j] -= 1
        d1[vi[t]][j] += 1
    d2 = [[0]*Fc for _ in range(E)]
    for k, fc in enumerate(faces):
        for se in fc:
            d2[abs(se)-1][k] += (1 if se > 0 else -1)
    return dict(verts=verts, edges=edges, faces=faces, V=V, E=E, F=Fc, d1=d1, d2=d2, vi=vi)

def topo(K):
    r1, r2 = rank_Q(K['d1']), rank_Q(K['d2'])
    b0 = K['V'] - r1
    b1 = K['E'] - r1 - r2
    b2 = K['F'] - r2
    # d1 . d2
    prod = [[sum(K['d1'][i][k]*K['d2'][k][j] for k in range(K['E'])) for j in range(K['F'])] for i in range(K['V'])]
    mx = max((abs(x) for row in prod for x in row), default=0)
    chi = K['V'] - K['E'] + K['F']
    gauge = K['V'] - b0
    inv = K['E'] - gauge
    return dict(chi=chi, b0=b0, b1=b1, b2=b2, rank_d1=r1, rank_d2=r2, d1d2=mx,
                gauge=gauge, inv=inv, curv=r2, flat=b1)

def cycle_tests(K, chain):
    """chain: dict edge_index(0-based) -> coefficient.  Returns (is_cycle, bounds)."""
    vec = [0]*K['E']
    for e, s in chain.items(): vec[e] += s
    bd = [sum(K['d1'][i][j]*vec[j] for j in range(K['E'])) for i in range(K['V'])]
    is_cycle = all(x == 0 for x in bd)
    # bounds?  vec in image(d2)?  rank test
    cols = [[K['d2'][e][k] for e in range(K['E'])] for k in range(K['F'])]
    Mt = [list(col) for col in cols]                       # F x E
    r_before = rank_Q(Mt) if K['F'] else 0
    r_after = rank_Q(Mt + [vec])
    bounds = (r_after == r_before)
    return is_cycle, bounds

def classes(K, gF_v, gC_v):
    cnt = {'00': 0, '10': 0, '01': 0, '11': 0}
    for v in K['verts']:
        a = '1' if v in gF_v else '0'
        b = '1' if v in gC_v else '0'
        cnt[a+b] += 1
    return cnt

def indep(K, ch1, ch2):
    """S4's 'independent' column: the two cycle vectors are linearly independent in the cycle
    space, i.e. (f,c) is independently assignable by a real connection."""
    v1 = [0]*K['E']; v2 = [0]*K['E']
    for e, s in ch1.items(): v1[e] += s
    for e, s in ch2.items(): v2[e] += s
    return rank_Q([v1, v2]) == 2

def lam_jensen(pi_, N=1 << 20):
    p00, p10, p01, p11 = [float(x) for x in pi_]
    t = (np.arange(N) + 0.5) * (2*np.pi/N)
    z = np.exp(1j*t)
    A = np.abs(p00 + p10*z); B = np.abs(p01 + p11*z)
    return float(np.mean(np.log(np.maximum(A, B))))

# --- B4-SQUARE: S4's own reading.  sphere A = 4-cycle p-a1-q-a2-p with two 2-cells;
#     sphere B = 4-cycle p-b1-q-b2-p with two 2-cells; glued at p and q.
vertsS = ['p', 'q', 'a1', 'a2', 'b1', 'b2']
edgesS = [('p', 'a1'), ('a1', 'q'), ('q', 'a2'), ('a2', 'p'),      # 1..4  sphere A square
          ('p', 'b1'), ('b1', 'q'), ('q', 'b2'), ('b2', 'p')]      # 5..8  sphere B square
facesS = [[1, 2, 3, 4], [-1, -2, -3, -4], [5, 6, 7, 8], [-5, -6, -7, -8]]
B4sq = build(vertsS, edgesS, facesS)

# --- B4-TRIPENT: refuter A-1's second spindle.  sphere A = triangle p-q-r with two 2-cells;
#     sphere B = pentagon p-s1-q-s2-s3-p with two 2-cells; glued at p and q.
vertsT = ['p', 'q', 'r', 's1', 's2', 's3']
edgesT = [('p', 'q'), ('q', 'r'), ('r', 'p'),                                    # 1..3 triangle
          ('p', 's1'), ('s1', 'q'), ('q', 's2'), ('s2', 's3'), ('s3', 'p')]      # 4..8 pentagon
facesT = [[1, 2, 3], [-1, -2, -3], [4, 5, 6, 7, 8], [-4, -5, -6, -7, -8]]
B4tp = build(vertsT, edgesT, facesT)

say("")
say("S4:519 publishes for B4:   V=6 E=8 F=4 chi=2 b0=1 b1=1 b2=2 gauge=5 inv=3 curv=2 flat=1,")
say("                           d1.d2=0, gF bounds, gC does NOT bound, independent")
say("S4:580 publishes for B4:   class multiset {00:1, 01:1, 10:1, 11:3}, lambda(SENSE U) = log(1/2)")
say("")
for name, K, gF, gC, gFv, gCv in [
    ("B4-SQUARE  (square|square)", B4sq, {0: 1, 1: 1, 2: 1, 3: 1}, {0: 1, 1: 1, 5: -1, 4: -1},
     {'p', 'a1', 'q', 'a2'}, {'p', 'a1', 'q', 'b1'}),
    ("B4-TRIPENT (triangle|pent)", B4tp, {0: 1, 1: 1, 2: 1}, {0: 1, 5: 1, 6: 1, 7: 1},
     {'p', 'q', 'r'}, {'p', 'q', 's2', 's3'}),
]:
    T = topo(K)
    cF, bF = cycle_tests(K, gF); cC, bC = cycle_tests(K, gC)
    ind = indep(K, gF, gC)
    cm = classes(K, gFv, gCv)
    pi_ = [F(cm['00'], K['V']), F(cm['10'], K['V']), F(cm['01'], K['V']), F(cm['11'], K['V'])]
    lam = lam_jensen(pi_)
    fact = (pi_[0]*pi_[3] == pi_[1]*pi_[2])
    say("  %s" % name)
    say("     V=%d E=%d F=%d chi=%d b0=%d b1=%d b2=%d gauge=%d inv=%d curv=%d flat=%d  max|d1.d2|=%d"
        % (K['V'], K['E'], K['F'], T['chi'], T['b0'], T['b1'], T['b2'], T['gauge'], T['inv'],
           T['curv'], T['flat'], T['d1d2']))
    say("     gF cycle=%s bounds=%s   gC cycle=%s bounds=%s   independent=%s"
        % (cF, bF, cC, bC, ind))
    say("     class multiset {00:%d, 10:%d, 01:%d, 11:%d}   pi = %s"
        % (cm['00'], cm['10'], cm['01'], cm['11'], tuple(str(x) for x in pi_)))
    say("     p00*p11 == p10*p01 (P factors)?  %s" % fact)
    say("     lambda(SENSE U), Jensen 2^20 = %.12f%s" % (lam,
        ("   [P FACTORS (p00 p11 = p10 p01); with non-negative coefficients that forces"
         " lambda = log(max p_c) = log(%s) = %.15f EXACTLY, quadrature error %.1e --"
         " the 2^20 loss is the log singularity of the torus zero the factorisation creates]"
         % (max(pi_), float(np.log(float(max(pi_)))),
            abs(lam - float(np.log(float(max(pi_))))))) if fact else
        "   [no torus zero at these weights; branch domination gives log(%s) exactly]"
        % max(pi_[0], pi_[1], pi_[2], pi_[3])))
say("")
say("  exact closed forms:  log(1/2) = %.12f    log(1/3) = %.12f" % (np.log(0.5), np.log(1/3)))
say("")
say("  RULING ON THE A-1 / A-2 COLLISION.  Both refuters are right and they answer different")
say("  questions.  A-2 asked: GIVEN the square|square complex, is gamma_C forced?  It is.")
say("  A-1 asked: is the COMPLEX forced by what S4 published?  It is not: B4-TRIPENT above matches")
say("  EVERY published column of S4's B4 row and delivers a different class multiset and a")
say("  different lambda.  S4's B4 row is UNDER-DETERMINED by its own parameters -- a COR-K defect")
say("  against S4 -- and any lane that hard-codes S4's multiset as its build target is running a")
say("  satisfiability check, not an audit.  The scope table must therefore mark the B4 row's")
say("  lambda and multiset UNDETERMINED-AT-SOURCE, not merely unaudited.")

# ----------------------------------------------------------------------------------------------
# PART C -- ARITHMETIC TYPING OF EVERY CONNECTION THE CORPUS PUBLISHES
# ----------------------------------------------------------------------------------------------
say("")
say("PART C -- ARITHMETIC TYPING OF THE CORPUS'S PUBLISHED CONNECTIONS")
say("-"*100)
say("relation lattice L = {(m,n) in Z^2 : m f + n c in 2 pi Z}.  THE TYPING RULE, which is a")
say("two-line theorem and not a measurement:  if f and c are both RATIONAL then m f + n c is")
say("rational, and 2 pi j is irrational for j != 0, so every relation has j = 0 and a nonzero")
say("(m,n) with m f + n c = 0 always exists.  EVERY RATIONAL CONNECTION IS EXACTLY RESONANT.")
say("")
conns = [
    ("S3/S4 headline           f=2.0     c=1.1",      F(2), F(11, 10)),
    ("S3 sec6 diagonal         f=2.0     c=2.0",      F(2), F(2)),
    ("S4:973 'they are generic' f=3.14159 c=1.57080", F(314159, 100000), F(157080, 100000)),
    ("lane D 'generic'         f=1.3     c=2.0",      F(13, 10), F(2)),
    ("lane D 'golden'          f=2pi/phi c=2pi/phi^2", None, None),
    ("S1 sec6 published        f=pi      c=3pi/2",    None, None),
    ("S4:603 verification      f=1.0     c=sqrt2",    None, None),
]
def prim_rel(fr_f, fr_c):
    """f = a/b, c = d/e rationals as Fractions.  Smallest (m,n) with m f + n c = 0."""
    num_f, den_f = fr_f.numerator, fr_f.denominator
    num_c, den_c = fr_c.numerator, fr_c.denominator
    # m*(num_f/den_f) = -n*(num_c/den_c) -> m*num_f*den_c = -n*num_c*den_f
    A = num_f*den_c; B = num_c*den_f
    g = np.gcd(abs(A), abs(B)) if (A or B) else 1
    m, nn = B//g, -A//g
    assert m*fr_f + nn*fr_c == 0
    return m, nn
for lbl, ff, cc in conns:
    if ff is not None:
        m, nn = prim_rel(ff, cc)
        say("  %-48s  RATIONAL  -> EXACTLY RESONANT, primitive relation (m,n) = (%d,%d)" % (lbl, m, nn))
    elif "golden" in lbl:
        phi = (1+5**0.5)/2
        say("  %-48s  1/phi + 1/phi^2 = 1 EXACTLY, so f + c = 2 pi -> EXACTLY RESONANT, (1,1)"
            % lbl)
        say("       residual |f + c - 2pi| = %.3e (float64 only; the identity is exact)"
            % abs(2*np.pi/phi + 2*np.pi/phi**2 - 2*np.pi))
    elif "S1 sec6" in lbl:
        say("  %-48s  f=pi, c=3pi/2 -> both in 2 pi Q -> FINITE ORDER (branch ratio order 4)" % lbl)
    else:
        say("  %-48s  f rational, c = sqrt2 irrational, pi transcendental -> NO RELATION: GENERIC"
            % lbl)
say("")
say("  RULING.  Of the connections the corpus publishes, EXACTLY ONE is generic in the sense the")
say("  corpus's own limit theorems require -- S4:603's f = 1.0, c = sqrt(2), the connection S4 used")
say("  to VERIFY the lambda column.  S4:973's 'the truncated decimals 3.14159 and 1.57080 ... are")
say("  GENERIC' is FALSE AT THE BYTES, and it sits inside the paragraph correcting the corpus's")
say("  first two mislabelled connections.  This is the corpus's THIRD mislabelled connection and")
say("  lane D's (1.3, 2.0) is the FOURTH.  Verified above by exact rational arithmetic.")

# ----------------------------------------------------------------------------------------------
# PART D -- THE MULTISET THEOREM'S HYPOTHESIS.  THREE RIVAL NAMES, ALL REFUTED BY ONE EXHIBIT.
# ----------------------------------------------------------------------------------------------
say("")
say("PART D -- THE MULTISET THEOREM OFF THE NON-NEGATIVE LOCUS: three names on the table, one exhibit")
say("-"*100)
say("names offered this round:")
say("  registrar (W-09 brief) : REAL NON-NEGATIVE coefficients")
say("  lane D / W10D-13       : REALITY of the coefficients ('complex drops |G| to exactly 8')")
say("  refuter A-2            : cos D1 = cos D2 = cos D3 on Delta = arg(p00 p11 / (p10 p01))")
say("  refuter A-1            : all four coefficients share one argument mod pi (collinearity)")
say("  refuter D-1            : how many of the four Jensen-adjacent pairs have conj(p_i)p_j real")
say("                           (4 or 2 -> 24, 1 -> 16, 0 -> 8)")
say("")
mp.mp.dps = 40
def m_perm(p, N=None):
    """m(P) at dps 40 by Gauss-Legendre on the Jensen-in-y reduction."""
    a, b, c, d = [mp.mpmathify(z) for z in p]
    fn = lambda t: mp.log(max(abs(a + b*mp.e**(1j*t)), abs(c + d*mp.e**(1j*t))))
    return mp.quad(fn, [0, mp.pi/2, mp.pi, 3*mp.pi/2, 2*mp.pi])/(2*mp.pi)

def spread(p):
    vals = [m_perm([p[i] for i in s]) for s in permutations(range(4))]
    return float(max(vals) - min(vals)), vals

def cosD(p):
    """the three pair-partition fluxes D = arg(prod A / prod B); a permutation carries the
    gauge-invariant phase into +-D1, +-D2, +-D3, so m depends on the permutation only through
    which partition is diagonal, and only through cos D."""
    import cmath
    a = [cmath.phase(z) for z in p]
    Ds = [a[0]+a[3]-a[1]-a[2], a[0]+a[1]-a[2]-a[3], a[0]+a[2]-a[1]-a[3]]
    return [float(np.cos(d)) for d in Ds]

def adj_real_pairs(p):
    """the four Jensen-adjacent pairs {00,10},{01,11},{00,01},{10,11}: is conj(p_i)p_j real?"""
    pairs = [(0, 1), (2, 3), (0, 2), (1, 3)]
    return sum(1 for i, j in pairs if abs((np.conj(p[i])*p[j]).imag) < 1e-13)

rng = np.random.default_rng(20260816)
say("  exhibit 1  DOMINATED, NON-COLLINEAR, ZERO adjacent-real pairs, arbitrary phases:")
ph = [0.0, 0.7, 1.9, 2.9]
p1 = [10*np.exp(1j*ph[0]), 1*np.exp(1j*ph[1]), 1*np.exp(1j*ph[2]), 1*np.exp(1j*ph[3])]
s1, _ = spread(p1)
say("     p = (10 e^{i0}, e^{i0.7}, e^{i1.9}, e^{i2.9})   adjacent-real pairs = %d" % adj_real_pairs(p1))
say("     24-permutation spread = %.6e      m = log 10 = %.12f  (|10 +- 1| >= 9 > 2 >= other branch)"
    % (s1, float(mp.log(10))))
say("     cos D1, cos D2, cos D3 = %.6f  %.6f  %.6f   -- PAIRWISE DISTINCT" % tuple(cosD(p1)))
say("  exhibit 2  same moduli, NOT dominated (all moduli comparable), same kind of phases:")
p2 = [1.0*np.exp(1j*0.0), 0.9*np.exp(1j*0.7), 1.1*np.exp(1j*1.9), 0.95*np.exp(1j*2.9)]
s2, _ = spread(p2)
say("     p = (1, 0.9 e^{i0.7}, 1.1 e^{i1.9}, 0.95 e^{i2.9})  adjacent-real pairs = %d" % adj_real_pairs(p2))
say("     24-permutation spread = %.6e     cos D = %.6f %.6f %.6f" % ((s2,) + tuple(cosD(p2))))
say("  exhibit 3  REAL NON-NEGATIVE (the corpus's actual objects: pi is a probability vector):")
p3 = [F(4, 9), F(2, 9), F(1, 9), F(2, 9)]
s3, _ = spread([float(x) for x in p3])
say("     p = B0b's own pi = (4/9,2/9,1/9,2/9)            24-permutation spread = %.6e" % s3)
say("     cos D = %.6f %.6f %.6f  -- all 1, as they are for EVERY probability vector" % tuple(cosD([float(x) for x in p3])))
say("  exhibit 4  a random NON-collinear complex draw, undominated:")
p4 = list(rng.normal(size=4) + 1j*rng.normal(size=4))
s4, _ = spread(p4)
say("     adjacent-real pairs = %d                          24-permutation spread = %.6e"
    % (adj_real_pairs(p4), s4))
say("")
say("  RULING ON THE NAMING.  Exhibit 1 has ZERO adjacent-real pairs, is not collinear, and has")
say("  cos D1, cos D2, cos D3 pairwise distinct -- and is FULLY S4-INVARIANT.  So all three of")
say("  'reality', 'collinearity/one argument mod pi', 'cos D1 = cos D2 = cos D3' and D-1's")
say("  adjacent-pair trichotomy are SUFFICIENT-ONLY at best, and D-1's trichotomy ('0 real pairs")
say("  -> exactly 8') is FALSE.  TWO INDEPENDENT MECHANISMS produce invariance: (i) coincidence of")
say("  the three matching fluxes (phase-side), (ii) branch domination, under which m collapses to")
say("  log of the largest modulus and the flux drops out entirely (modulus-side).  No name on the")
say("  table covers both.  SIX consecutive names, three of them produced THIS ROUND by lanes")
say("  commissioned to catch misnaming.  I decline to supply a seventh and mark it UNDETERMINED.")
say("  AND IT DOES NOT MATTER FOR ANY REGISTERED RESULT: pi is a probability vector, always real")
say("  and non-negative, so exhibit 3's regime is the only one the corpus ever occupies.")

# ----------------------------------------------------------------------------------------------
# PART E -- THE CHARACTER IDENTITY.  DECLARED A CONTROL THAT COULD NOT HAVE FAILED.
# ----------------------------------------------------------------------------------------------
say("")
say("PART E -- W-08's CHARACTER IDENTITY, EXACT.  DECLARED IN ADVANCE: THIS IS A CONTROL AND IT")
say("          COULD NOT HAVE FAILED.  It is an algebraic identity in the group ring, with no")
say("          class-count index anywhere in its proof:")
say("            1 - |Z|^2 = |sum_j w_j|^2 - |sum_j w_j chi_j|^2 = sum_{j<l} w_j w_l |chi_j - chi_l|^2")
say("          for ANY finite set of unit-modulus chi_j and non-negative w_j summing to 1.")
say("          I run it once, exactly, only so the number is on a page and not in a sentence.")
say("-"*100)
def check_identity(q, exps, ws, kmax=8):
    """exact in Z[x]/(x^q-1) with Fraction weights; conj is x -> x^{q-1}."""
    worst = F(0)
    for k in range(1, kmax+1):
        # Z_k = sum w_j x^{k e_j};  |Z_k|^2 = Z_k * conj(Z_k)
        lhs = [F(0)]*q
        for wj, ej in zip(ws, exps):
            for wl, el in zip(ws, exps):
                lhs[(k*ej - k*el) % q] += wj*wl
        # 1 - |Z_k|^2  vs  sum_{j<l} w_j w_l |chi_j^k - chi_l^k|^2
        rhs = [F(0)]*q
        for (wj, ej), (wl, el) in combinations(list(zip(ws, exps)), 2):
            rhs[0] += 2*wj*wl
            rhs[(k*ej - k*el) % q] -= wj*wl
            rhs[(k*el - k*ej) % q] -= wj*wl
        one_minus = [F(0)]*q
        one_minus[0] += F(1)
        for i in range(q): one_minus[i] -= lhs[i]
        res = max(abs(one_minus[i] - rhs[i]) for i in range(q))
        worst = max(worst, res)
    return worst
cases = [(12, [0, 1, 5, 7], [F(4, 9), F(2, 9), F(1, 9), F(2, 9)], "4 classes, q=12, B0b weights"),
         (7,  [0, 2, 3, 5], [F(1, 6), F(1, 6), F(1, 6), F(1, 2)], "4 classes, q=7,  B4 weights"),
         (5,  [1, 2, 4],    [F(2, 5), F(2, 5), F(1, 5)],          "3 classes, q=5,  K1 weights"),
         (37, [0, 3, 11, 29], [F(1, 4), F(1, 4), F(1, 4), F(1, 4)], "4 classes, q=37, SENSE C"),
         (6,  [0, 3],       [F(1, 2), F(1, 2)],                    "2 classes, q=6")]
for q, e, w, lbl in cases:
    say("  %-34s  worst exact residual over k<=8 : %s" % (lbl, check_identity(q, e, w)))
say("")
say("  As declared: zero, and it had to be.  It is entered as a THEOREM with a proof, not as")
say("  evidence.  Its scope is CARRIER_INDEPENDENT and CLASS-COUNT-INDEPENDENT by the proof, not")
say("  by these five rows.")

say("")
say("="*100)
say("END.  Every figure above is reproducible from this file alone.")
say("="*100)

with open(__file__.replace('.py', '.OUT.txt'), 'w') as fh:
    fh.write("\n".join(OUT) + "\n")
