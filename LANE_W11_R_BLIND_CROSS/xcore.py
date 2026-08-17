"""LANE_W11_R_BLIND_CROSS -- shared core, written for the CROSS-REFUTATION of LANE_W11_R_BLIND.

Structurally independent of BOTH prior implementations by construction:
  * the registrar builds T entrywise from (src,dst,edge-index) triples;
  * the blind lane builds T entrywise from a Carrier class holding per-step phases;
  * I build T as an explicit factorisation  T = P . D  (P = the loop's cyclic permutation
    matrix extended by the identity off the loop, D = diagonal of the per-step phases),
    and CROSS-CHECK that factorisation against an entrywise build.  The factorisation is what
    makes T^L = M a combinatorial identity rather than a numerical coincidence, so I want it
    to be the object in the code, not a remark about the code.

CONVENTIONS (published; see CONVENTIONS.txt):
  * numpy float64 (double precision) is the default everywhere.  Exact checks say so.
  * a LOOP is a directed simple cycle given as a vertex list [w_0..w_{L-1}] meaning
    w_0 -> w_1 -> ... -> w_{L-1} -> w_0, together with the L phases picked up on those steps.
  * CLASS of a vertex = (v on gamma_F ?, v on gamma_C ?), written 00/10/01/11.
  * pi = (p00, p10, p01, p11), p_ab = sum over class ab of |s_v|^2.
  * ADVANCE PAIR (a,b): the functional evaluated is Y(a,b) = <T_F^a s, T_C^b s>.
      CIRCUIT convention  = the line (a,b) = (L_F k, L_C k)   [the corpus's]
      EDGE    convention  = the line (a,b) = (n, n)           [COR-F's]
  * u = exp(i f), v = exp(i c) with f, c the loop holonomy ANGLES.  (Sign conventions for
    which character is conjugated are irrelevant to every claim tested here; |Z| is used.)
"""
import numpy as np
from fractions import Fraction as Fr

# ------------------------------------------------------------------ carrier

class Cx:
    def __init__(self, name, nV, loopF, loopC, phF, phC):
        self.name, self.nV = name, nV
        self.loopF, self.loopC = list(loopF), list(loopC)
        self.phF, self.phC = np.asarray(phF, float), np.asarray(phC, float)
        assert len(self.phF) == len(self.loopF) and len(self.phC) == len(self.loopC)
        self.LF, self.LC = len(self.loopF), len(self.loopC)
        self.f, self.c = float(self.phF.sum()), float(self.phC.sum())

    def _loop(self, w): return (self.loopF, self.phF) if w == 'F' else (self.loopC, self.phC)

    def perm(self, w):
        """P: the loop's cyclic permutation, identity off the loop.  P[w_{j+1}, w_j] = 1."""
        loop, _ = self._loop(w); L = len(loop)
        P = np.eye(self.nV, dtype=complex)
        for j in range(L): P[loop[j], loop[j]] = 0.0
        for j in range(L): P[loop[(j+1) % L], loop[j]] = 1.0
        return P

    def dphase(self, w):
        """D: diagonal, exp(i*phase of the step LEAVING each loop vertex), 1 off the loop."""
        loop, ph = self._loop(w)
        d = np.ones(self.nV, dtype=complex)
        for j, vtx in enumerate(loop): d[vtx] = np.exp(1j*ph[j])
        return np.diag(d)

    def T(self, w):
        """EDGE tick, built as P.D (COR-F's operator)."""
        return self.perm(w) @ self.dphase(w)

    def T_entrywise(self, w):
        """independent entrywise build, used only to check the P.D factorisation."""
        loop, ph = self._loop(w); L = len(loop)
        A = np.zeros((self.nV, self.nV), dtype=complex)
        onloop = set(loop)
        for x in range(self.nV):
            if x not in onloop: A[x, x] = 1.0
        for j in range(L): A[loop[(j+1) % L], loop[j]] = np.exp(1j*ph[j])
        return A

    def M(self, w):
        """CIRCUIT operator: multiply every loop vertex by the whole-circuit holonomy."""
        loop, ph = self._loop(w)
        d = np.ones(self.nV, dtype=complex)
        for vtx in loop: d[vtx] = np.exp(1j*float(ph.sum()))
        return np.diag(d)

    def S(self, w):
        """'SMEARED' diagonal L-th root of M: an L-th of the holonomy per tick, on the loop."""
        loop, ph = self._loop(w); L = len(loop)
        d = np.ones(self.nV, dtype=complex)
        for vtx in loop: d[vtx] = np.exp(1j*float(ph.sum())/L)
        return np.diag(d)

    def onF(self): 
        m = np.zeros(self.nV, bool); m[self.loopF] = True; return m
    def onC(self):
        m = np.zeros(self.nV, bool); m[self.loopC] = True; return m

    def classes(self):
        F, C = self.onF(), self.onC()
        return {(0,0): np.where(~F & ~C)[0], (1,0): np.where(F & ~C)[0],
                (0,1): np.where(~F & C)[0], (1,1): np.where(F & C)[0]}

    def pi(self, s):
        p = np.abs(s)**2; cl = self.classes()
        return np.array([p[cl[(0,0)]].sum(), p[cl[(1,0)]].sum(),
                         p[cl[(0,1)]].sum(), p[cl[(1,1)]].sum()])

ORDER = [(0,0),(1,0),(0,1),(1,1)]

# ------------------------------------------------------------------ carriers of record

def K1(f, c, splitF=(0.31,0.47,1.0), splitC=(0.23,0.91,1.0)):
    """S1 sec.1.  v0..v4.  gamma_F = v0->v1->v2->v0 (the filled triangle, e1 e2 e3);
       gamma_C = v0->v3->v4->v0 (the unfilled triangle, e4 e5 e6).  COR-F's own loop."""
    a = np.array(splitF, float); a = a/a.sum()*f
    b = np.array(splitC, float); b = b/b.sum()*c
    return Cx('K1', 5, [0,1,2], [0,3,4], a, b)

def B0b_registrar(f, c, g=0.37, splitF=(0.29,0.61,1.0), splitC=(0.44,1.0)):
    """3x3 ring torus, V=9, vertex (i,j) -> 3*j + i  (the registrar's index map).
       gamma_F = the (0,0) square face 0->1->4->3->0, L=4, BOUNDS.
       gamma_C = the j=0 horizontal row 0->1->2->0, L=3, does not bound.
       The two loops share the DIRECTED edge 0->1, traversed the SAME way by both.
       Class multiset {11:{0,1}, 10:{3,4}, 01:{2}, 00:{5,6,7,8}} = {00:4,01:1,10:2,11:2}."""
    # shared edge 0->1 carries phase g, traversed forward by BOTH loops
    a = np.array(splitF, float); a = a/a.sum()*(f - g)
    phF = np.concatenate([[g], a])                       # steps 0->1 , 1->4, 4->3, 3->0
    b = np.array(splitC, float); b = b/b.sum()*(c - g)
    phC = np.concatenate([[g], b])                       # steps 0->1 , 1->2, 2->0
    return Cx('B0b/reg', 9, [0,1,4,3], [0,1,2], phF, phC)

def B0b_blind(f, c, g=0.37, splitF=(0.29,0.61,1.0), splitC=(0.44,1.0)):
    """The blind lane's reconstruction: same square face, but gamma_C is the COLUMN 0->3->6->0.
       The two loops share the UNDIRECTED edge {0,3}, traversed in OPPOSITE directions.
       Same class multiset.  Included so that every B0b claim is run on BOTH reconstructions."""
    a = np.array(splitF, float); a = a/a.sum()*(f + g)
    phF = np.concatenate([a, [-g]])                      # 0->1, 1->4, 4->3, then 3->0 = -g
    b = np.array(splitC, float); b = b/b.sum()*(c - g)
    phC = np.concatenate([[g], b])                       # 0->3 = +g, 3->6, 6->0
    return Cx('B0b/blind', 9, [0,1,4,3], [0,3,6], phF, phC)

# ------------------------------------------------------------------ functionals

def Y(car, s, a, b, opF=None, opC=None):
    """<A_F^a s, A_C^b s>.  Defaults to the EDGE tick T."""
    AF = car.T('F') if opF is None else opF
    AC = car.T('C') if opC is None else opC
    return np.vdot(np.linalg.matrix_power(AF, a) @ s, np.linalg.matrix_power(AC, b) @ s)

def traj(car, s, N, opF, opC, adv):
    """|Z_t| for t = 1..N, per-tick operators opF/opC, advance rule adv(t) -> (a,b)."""
    out = np.empty(N)
    xF = s.copy(); xC = s.copy()
    a0 = b0 = 0
    for t in range(1, N+1):
        a, b = adv(t)
        for _ in range(a - a0): xF = opF @ xF
        for _ in range(b - b0): xC = opC @ xC
        a0, b0 = a, b
        out[t-1] = abs(np.vdot(xF, xC))
    return out

def rate_of(z):  return float(np.mean(np.log(z)))

def mahler4(c00, c10, c01, c11, ngrid=1 << 20):
    """m(c00 + c10 x + c01 y + c11 x y).  Jensen in y: m = mean_x log max(|c00+c10 x|,|c01+c11 x|).
       Independent of the blind lane's routine only in that it is written here; the reduction is
       the same standard one.  A second, Jensen-free 2-D quadrature is run in X6."""
    th = np.arange(ngrid)*(2*np.pi/ngrid); x = np.exp(1j*th)
    return float(np.mean(np.log(np.maximum(np.abs(c00 + c10*x), np.abs(c01 + c11*x)))))

def dressed_coeffs(car, s, r, opF=None, opC=None):
    """c^(r)_ab = sum_{w in class ab} conj((A_F^r s)(w)) * (A_C^r s)(w)."""
    AF = car.T('F') if opF is None else opF
    AC = car.T('C') if opC is None else opC
    tr = np.linalg.matrix_power(AF, r) @ s
    wr = np.linalg.matrix_power(AC, r) @ s
    prod = np.conj(tr)*wr; cl = car.classes()
    return np.array([prod[cl[o]].sum() for o in ORDER])

def invisible(car, states, a, b, opF=None, opC=None, tol=1e-12):
    vals = [abs(Y(car, s, a, b, opF, opC)) for s in states]
    return (max(vals) - min(vals)) < tol, max(vals) - min(vals)

def equal_pi_triple(car, base, moves, phases, rng=None):
    """Build three states with IDENTICAL pi: A = base; B = base with weight moved WITHIN
       classes; C = base with phases only.  Returns (sA,sB,sC) and asserts pi agreement."""
    w = np.array(base, float); w = w/w.sum()
    wB = w.copy()
    for (i, j, d) in moves: wB[i] += d; wB[j] -= d
    sA = np.sqrt(w).astype(complex)
    sB = np.sqrt(wB).astype(complex)
    sC = np.sqrt(w)*np.exp(1j*np.asarray(phases, float))
    assert np.allclose(car.pi(sA), car.pi(sB), atol=1e-14), "arm B moved pi -- NOT a valid arm"
    assert np.allclose(car.pi(sA), car.pi(sC), atol=1e-14), "arm C moved pi -- NOT a valid arm"
    assert np.linalg.norm(sA-sB) > 1e-3 and np.linalg.norm(sA-sC) > 1e-3, "ZERO VARIABLES MOVED"
    return sA, sB, sC


def auto_moves(car, w, frac=0.5):
    """Within-class weight moves derived FROM THE CARRIER'S OWN CLASSES, so that arm B can
       never accidentally move pi.  For every class with >=2 vertices, shift `frac` of the
       first vertex's weight onto the second."""
    mv = []
    for o in ORDER:
        idx = car.classes()[o]
        if len(idx) >= 2:
            i, j = int(idx[0]), int(idx[1])
            mv.append((j, i, frac*float(w[i])))
    return mv
