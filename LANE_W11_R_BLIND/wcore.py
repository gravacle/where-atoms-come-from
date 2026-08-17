"""
LANE_W11_R_BLIND -- shared core. Built from definitions only:
  S1_CARRIER_K1_V001.md (complex, fibres, connection, gauge)
  S3_THE_CROSSING_AUDIT_V001.md:160-215  (COR-F: the edge-tick transport T)
  REGISTER_V001.md W-01 (M_gamma), W-05 N1 (lambda = m(p00+p10x+p01y+p11xy))
No registrar code was read before this file was written.

CONVENTIONS (published):
  * numpy float64 (double precision) is the default everywhere.
  * fibre order on K1 is [v0,v1,v2,v3,v4]; on B0b it is index 3*i+j for (i,j) in Z3xZ3.
  * u := conj(W_F) = exp(-i f)      v := W_C = exp(+i c)      (W-01's characters)
  * class (a,b) := (is vertex on gamma_F ?, is vertex on gamma_C ?)
  * pi := (p00, p10, p01, p11), p_ab = sum_{v in class ab} |s_v|^2
  * P(x,y) := p00 + p10 x + p01 y + p11 x y                    (N1)
"""
import numpy as np

# ---------------------------------------------------------------- carriers

class Carrier:
    def __init__(self, name, nV, loopF, loopC, phF, phC):
        """loopF/loopC: vertex cycles [v_0,...,v_{L-1}] meaning v_0->v_1->...->v_0.
           phF/phC: the connection phase picked up traversing each of those L steps."""
        self.name, self.nV = name, nV
        self.loopF, self.loopC = list(loopF), list(loopC)
        self.phF, self.phC = np.asarray(phF, float), np.asarray(phC, float)
        assert len(self.phF) == len(self.loopF) and len(self.phC) == len(self.loopC)
        self.LF, self.LC = len(loopF), len(loopC)
        self.inF = np.zeros(nV, bool); self.inF[self.loopF] = True
        self.inC = np.zeros(nV, bool); self.inC[self.loopC] = True
        self.f = float(self.phF.sum())          # face/loop-F holonomy angle
        self.c = float(self.phC.sum())          # loop-C holonomy angle

    # ---- CIRCUIT convention (the corpus's): fibre-wise mult by whole-circuit holonomy
    def M(self, which):
        ang = self.f if which == 'F' else self.c
        onloop = self.inF if which == 'F' else self.inC
        d = np.ones(self.nV, complex)
        d[onloop] = np.exp(1j*ang)
        return np.diag(d)

    # ---- EDGE convention (COR-F's): move each fibre value one edge along the loop
    def T(self, which):
        loop = self.loopF if which == 'F' else self.loopC
        ph   = self.phF   if which == 'F' else self.phC
        L = len(loop)
        A = np.eye(self.nV, dtype=complex)
        for j in range(L):
            src, tgt = loop[j], loop[(j+1) % L]
            A[src, src] = 0.0
        for j in range(L):
            src, tgt = loop[j], loop[(j+1) % L]
            A[tgt, src] = np.exp(1j*ph[j])      # (T s)(tgt) = U_e s(src)
        return A

    def classes(self):
        """returns index arrays for classes 00,10,01,11"""
        return {(0,0): np.where(~self.inF & ~self.inC)[0],
                (1,0): np.where( self.inF & ~self.inC)[0],
                (0,1): np.where(~self.inF &  self.inC)[0],
                (1,1): np.where( self.inF &  self.inC)[0]}

    def pi(self, s):
        p = np.abs(s)**2
        cl = self.classes()
        return np.array([p[cl[(0,0)]].sum(), p[cl[(1,0)]].sum(),
                         p[cl[(0,1)]].sum(), p[cl[(1,1)]].sum()])


def K1(f, c, split_f=None, split_c=None):
    """S1's carrier. v0..v4; gamma_F = v0->v1->v2->v0 (filled), gamma_C = v0->v3->v4->v0."""
    if split_f is None: split_f = np.array([0.31, 0.47, 1.0]);  split_f = split_f/split_f.sum()*f
    if split_c is None: split_c = np.array([0.23, 0.91, 1.0]);  split_c = split_c/split_c.sum()*c
    return Carrier('K1', 5, [0,1,2], [0,3,4], split_f, split_c)


def B0b(f, c, g=0.37):
    """S4:512 'ring torus 3x3 grid, loops meet': V=9. Reconstructed from S4's own published
    class multiset {00:4, 01:1, 10:2, 11:2} and 'gF bounds True / gC bounds False'.
    vertex (i,j) -> 3i+j.  gamma_F = square face boundary 0->1->4->3->0  (L=4, bounds)
                           gamma_C = column cycle       0->3->6->0       (L=3, does not bound)
    They meet in {0,3}: 11:{0,3}=2, 10:{1,4}=2, 01:{6}=1, 00:{2,5,7,8}=4.  MATCHES S4:575."""
    # shared edge 0->3 carries phase g; gamma_F traverses it as 3->0 hence -g.
    pf = np.array([0.29, 0.61, 1.0]); pf = pf/pf.sum()*(f + g)
    phF = np.concatenate([pf, [-g]])                    # steps 0->1,1->4,4->3,3->0
    pc = np.array([0.44, 1.0]);  pc = pc/pc.sum()*(c - g)
    phC = np.concatenate([[g], pc])                     # steps 0->3,3->6,6->0
    return Carrier('B0b', 9, [0,1,4,3], [0,3,6], phF, phC)


# ---------------------------------------------------------------- functionals

def Z_circuit(car, s, K):
    """Z_k = <M_F^k s, M_C^k s>, k = 1..K, by direct matrix action (no shortcut)."""
    MF, MC = car.M('F'), car.M('C')
    a = s.copy(); b = s.copy(); out = np.empty(K, complex)
    for k in range(K):
        a = MF @ a; b = MC @ b
        out[k] = np.vdot(a, b)
    return out

def Z_edge(car, s, N):
    """Z^T_n = <T_F^n s, T_C^n s>, n = 1..N, by direct matrix action."""
    TF, TC = car.T('F'), car.T('C')
    a = s.copy(); b = s.copy(); out = np.empty(N, complex)
    for n in range(N):
        a = TF @ a; b = TC @ b
        out[n] = np.vdot(a, b)
    return out

def Z_pair(car, s, nF, nC):
    """<T_F^{nF} s, T_C^{nC} s> -- the single family both conventions sample."""
    TF, TC = car.T('F'), car.T('C')
    a = np.linalg.matrix_power(TF, nF) @ s
    b = np.linalg.matrix_power(TC, nC) @ s
    return np.vdot(a, b)

def rate(Z):
    z = np.abs(Z)
    return float(np.mean(np.log(z)))

# ---------------------------------------------------------------- Mahler measure

def mahler4(c00, c10, c01, c11, ngrid=1 << 20):
    """m(c00 + c10 x + c01 y + c11 x y) by Jensen in y:
         m = (1/2pi) \int log max(|A(x)|,|B(x)|) dtheta,  A = c00+c10 x, B = c01+c11 x.
       Valid for COMPLEX coefficients. Trapezoid on a uniform grid in theta."""
    th = np.arange(ngrid)*(2*np.pi/ngrid)
    x = np.exp(1j*th)
    A = c00 + c10*x
    B = c01 + c11*x
    mx = np.maximum(np.abs(A), np.abs(B))
    return float(np.mean(np.log(mx)))

# ---------------------------------------------------------------- gauge

def gauge(nV, rng):
    return np.exp(1j*rng.uniform(0, 2*np.pi, nV))
