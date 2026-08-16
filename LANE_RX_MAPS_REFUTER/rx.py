"""
REFUTER LANE -- "maps" axis.  Written from scratch.  numpy + mpmath only.
Conventions, seeds and incidence matrices are PUBLISHED in the output.

Objects
-------
Carrier  : vertices (names), edges (src,tgt), faces (edge-sign word), two designated
           loops gamma_F (bounds a face) and gamma_C (free), each an edge-sign word.
Fibre    : C^1 (U(1)) or C^2 (SU(2)).  charge q_v in Z on each vertex (U(1) only).
Transport: M_gamma acts on Gamma(L) = (+)_v L_v by the holonomy of gamma BASED AT v,
           and by the identity at v not on gamma.  (S3 sec 2.3; S4 sec 2.)
Z_k      : <M_F^k s, M_C^k s>.
lambda   : schedule A  lam_A = log|Z_1| ; schedule B (k_n = n) lam_B = lim (1/N) sum log|Z_n|.
"""
import numpy as np, itertools, math
from mpmath import mp, mpf, quad, log as mlog, cos as mcos, sqrt as msqrt, acos as macos, pi as mpi

mp.dps = 40
RNG_SEED = 20260816          # PUBLISHED SEED

# ----------------------------------------------------------------------------- carriers
class Carrier:
    def __init__(self, name, verts, edges, faces, gF, gC):
        self.name, self.V, self.E, self.F = name, list(verts), list(edges), list(faces)
        self.gF, self.gC = gF, gC          # list of (edge_index, +1/-1)
        self.vi = {v: i for i, v in enumerate(self.V)}
    def d1(self):
        """boundary_1 : C_1 -> C_0 , columns = edges, rows = vertices."""
        M = np.zeros((len(self.V), len(self.E)), dtype=int)
        for j, (s, t) in enumerate(self.E):
            M[self.vi[t], j] += 1
            M[self.vi[s], j] -= 1
        return M
    def d2(self):
        """boundary_2 : C_2 -> C_1 , columns = faces."""
        M = np.zeros((len(self.E), max(1, len(self.F))), dtype=int)
        for j, w in enumerate(self.F):
            for (e, sgn) in w:
                M[e, j] += sgn
        return M[:, :len(self.F)]
    def loop_vertices(self, g):
        vs = set()
        for (e, sgn) in g:
            s, t = self.E[e]
            vs.add(s); vs.add(t)
        return vs
    def classes(self):
        VF, VC = self.loop_vertices(self.gF), self.loop_vertices(self.gC)
        return {v: (int(v in VF), int(v in VC)) for v in self.V}
    def betti(self):
        V, E, F = len(self.V), len(self.E), len(self.F)
        r1 = np.linalg.matrix_rank(self.d1().astype(float)) if E else 0
        r2 = np.linalg.matrix_rank(self.d2().astype(float)) if F else 0
        b0 = V - r1
        b1 = (E - r1) - r2
        b2 = F - r2
        return b0, b1, b2, V - E + F
    def gauge_invariant_dim(self):
        """# independent gauge-invariant U(1) parameters = E - (V - b0)."""
        b0 = self.betti()[0]
        return len(self.E) - (len(self.V) - b0)

def K1():
    V = ['v0', 'v1', 'v2', 'v3', 'v4']
    E = [('v0','v1'), ('v1','v2'), ('v2','v0'),      # e1 e2 e3
         ('v0','v3'), ('v3','v4'), ('v4','v0')]      # e4 e5 e6
    Fc = [[(0,+1), (1,+1), (2,+1)]]                  # face F along e1 e2 e3
    gF = [(0,+1), (1,+1), (2,+1)]
    gC = [(3,+1), (4,+1), (5,+1)]
    return Carrier('B1 = K1', V, E, Fc, gF, gC)

def K1_subdivided():
    """B1s: every edge of K1 subdivided once.  11 vertices, 12 edges, 1 face."""
    V = ['v0','v1','v2','v3','v4','m1','m2','m3','m4','m5','m6']
    E = [('v0','m1'),('m1','v1'),          # e1a e1b   (was e1: v0->v1)
         ('v1','m2'),('m2','v2'),          # e2a e2b
         ('v2','m3'),('m3','v0'),          # e3a e3b
         ('v0','m4'),('m4','v3'),          # e4a e4b
         ('v3','m5'),('m5','v4'),          # e5a e5b
         ('v4','m6'),('m6','v0')]          # e6a e6b
    gF = [(i,+1) for i in range(0,6)]
    gC = [(i,+1) for i in range(6,12)]
    Fc = [list(gF)]
    return Carrier('B1s = K1 subdivided', V, E, Fc, gF, gC)

def K1_collapse_e4():
    """L4 : K1 / e4 .  v0 and v3 identified.  CLASS-MERGING ((1,1) with (0,1))."""
    V = ['w', 'v1', 'v2', 'v4']            # w = [v0 = v3]
    E = [('w','v1'), ('v1','v2'), ('v2','w'),   # e1 e2 e3
         ('w','v4'), ('v4','w')]                # e5 e6   (e4 collapsed)
    Fc = [[(0,+1),(1,+1),(2,+1)]]
    gF = [(0,+1),(1,+1),(2,+1)]
    gC = [(3,+1),(4,+1)]                   # image of e4.e5.e6 = (const).e5.e6
    return Carrier('L4 = K1/e4', V, E, Fc, gF, gC)

def K1_collapse_e2():
    """L2 : K1 / e2 .  v1 and v2 identified.  CLASS-COMPATIBLE (both (1,0))."""
    V = ['v0', 'w', 'v3', 'v4']            # w = [v1 = v2]
    E = [('v0','w'), ('w','v0'),           # e1 e3   (e2 collapsed)
         ('v0','v3'), ('v3','v4'), ('v4','v0')]
    Fc = [[(0,+1),(1,+1)]]
    gF = [(0,+1),(1,+1)]
    gC = [(2,+1),(3,+1),(4,+1)]
    return Carrier('L2 = K1/e2', V, E, Fc, gF, gC)

# --------------------------------------------------------------- U(1) transport, exact
def Z_k_direct(carrier, f, c, p, k, q=None):
    """Direct matrix action on Gamma(L)=C^V.  p = dict v->weight, q = dict v->charge."""
    V = carrier.V
    n = len(V)
    if q is None: q = {v: 1 for v in V}
    cls = carrier.classes()
    WF, WC = np.exp(1j*f), np.exp(1j*c)
    MF = np.diag([WF**(q[v]) if cls[v][0] else 1.0+0j for v in V])
    MC = np.diag([WC**(q[v]) if cls[v][1] else 1.0+0j for v in V])
    s = np.array([math.sqrt(p[v]) for v in V], dtype=complex)
    x = np.linalg.matrix_power(MF, k) @ s
    y = np.linalg.matrix_power(MC, k) @ s
    return np.vdot(x, y)

def class_weights(carrier, p):
    cls = carrier.classes()
    w = {(0,0):0.0, (1,0):0.0, (0,1):0.0, (1,1):0.0}
    for v in carrier.V:
        w[cls[v]] += p[v]
    return w

def Z_k_class(w, f, c, k):
    """w = dict class -> weight.  Z_k = sum_ab w_ab u^{ka} v^{kb}, u=e^{-if}, v=e^{ic}."""
    u, v = np.exp(-1j*f), np.exp(1j*c)
    return sum(w[(a,b)] * u**(k*a) * v**(k*b) for a in (0,1) for b in (0,1))

# ------------------------------------------------------- lambda_B : exact torus average
def mahler2(a1, a2, a3, a4):
    """
    m( a1 + a2 X + a3 Y + a4 XY ) for REAL non-negative a_i, by Jensen in Y:
        m = (1/2pi) int_0^{2pi} log max( |a1+a2 e^{it}| , |a3+a4 e^{it}| ) dt
    Crossing points solved in closed form; each smooth arc integrated by mpmath quad.
    """
    a1, a2, a3, a4 = [mpf(x) for x in (a1, a2, a3, a4)]
    A = a1*a1 + a2*a2 - a3*a3 - a4*a4
    B = 2*(a1*a2 - a3*a4)
    # |f|^2-|g|^2 = A + B cos t
    pts = [mpf(0), mpi]
    if B != 0:
        cc = -A/B
        if -1 < cc < 1:
            pts.append(macos(cc))
    pts = sorted(set(pts))
    def integ(t):
        f2 = a1*a1 + a2*a2 + 2*a1*a2*mcos(t)
        g2 = a3*a3 + a4*a4 + 2*a3*a4*mcos(t)
        m2 = f2 if f2 > g2 else g2
        floor = mpf(10)**(-2*mp.dps)      # integrable log singularity: floor, do not poison
        if m2 < floor: m2 = floor
        return mlog(m2)/2
    tot = mpf(0)
    for i in range(len(pts)-1):
        tot += quad(integ, [pts[i], pts[i+1]])
    return float(tot/mpi)          # (1/2pi)*2*int_0^pi  = (1/pi)*int_0^pi

def lambdaB_exact(w):
    """generic (f,c): relation lattice rank 0, lambda_B = m(p00 + p10 x + p01 y + p11 xy)."""
    return mahler2(w[(0,0)], w[(1,0)], w[(0,1)], w[(1,1)])

def lambdaB_direct(carrier, f, c, p, N, q=None):
    tot = 0.0
    for k in range(1, N+1):
        tot += math.log(abs(Z_k_direct(carrier, f, c, p, k, q)))
    return tot / N

def lambdaB_direct_class(w, f, c, N):
    tot = 0.0
    for k in range(1, N+1):
        tot += math.log(abs(Z_k_class(w, f, c, k)))
    return tot / N

# ----------------------------------------------------------------- formation criterion
def rankG(w, tol=0.0):
    S = [ab for ab, x in w.items() if x > tol]
    if len(S) <= 1: return 0, S
    # chi_(a,b) = u^a v^b ; exponent vector (a,b) in Z^2 (u,v independent generically)
    base = np.array(S[0])
    D = np.array([np.array(x) - base for x in S[1:]])
    return int(np.linalg.matrix_rank(D.astype(float))), S

def formation(w):
    r, S = rankG(w)
    return r >= 1, r, S
