# LANE W-11 R/C-CROSS — independent library.  Written from the sealed bytes; NOT imported from
# LANE_W11_R_CLOCK and NOT copied from LANE_W11_CONVENTION_TEST, so a shared implementation
# error cannot pass through.  Cross-validated against both at x1.
#   K1  : S1_CARRIER_K1_V001.md :10-25  (5 vertices, 6 edges, 2 triangles, root v0)
#   B0b : LANE_W10_A_CARRIERS/w10a_lib.py:182-213 (sealed) -> 3x3 torus grid, class multiset
#         {00:4, 01:1, 10:2, 11:2}  (S4_THE_MEASUREMENT_V001.md:575)
#   COR-F edge tick T : S3_THE_CROSSING_AUDIT_V001.md:160-209, :794
#   W-01 circuit op M : REGISTER_V001.md:31-35
# DOUBLE PRECISION IS THE DEFAULT.  Exact / closed-form checks are marked at the call site.
import numpy as np

class C:
    def __init__(self, name, nv, edges, wF, wC):
        self.name, self.nv, self.edges, self.wF, self.wC = name, nv, edges, wF, wC
        self.LF, self.LC = len(wF), len(wC)
        self.VF = tuple(sorted({u for u, _, _, _ in wF}))
        self.VC = tuple(sorted({u for u, _, _, _ in wC}))
        for w in (wF, wC):
            for i in range(len(w)):
                assert w[i][1] == w[(i+1) % len(w)][0]
            assert len({u for u, _, _, _ in w}) == len(w)
    def cls(self, v):  return (int(v in self.VF), int(v in self.VC))
    def multiset(self):
        from collections import Counter
        return dict(Counter("%d%d" % self.cls(v) for v in range(self.nv)))

def K1():
    E  = [(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)]
    wF = [(0,1,0,1),(1,2,1,1),(2,0,2,1)]
    wC = [(0,3,3,1),(3,4,4,1),(4,0,5,1)]
    return C("K1", 5, E, wF, wC)

def B0b():
    # independent re-derivation of the 3x3 grid torus: 9 vertices, 18 edges.
    V  = lambda i,j: 3*(j%3) + (i%3)
    Hi = lambda i,j: 3*(j%3) + (i%3)
    Wi = lambda i,j: 9 + 3*(j%3) + (i%3)
    E = [None]*18
    for j in range(3):
        for i in range(3):
            E[Hi(i,j)] = (V(i,j), V(i+1,j))
            E[Wi(i,j)] = (V(i,j), V(i,j+1))
    # gamma_F = boundary of the (0,0) plaquette;  gamma_C = the j=0 row (flat, wraps the torus)
    gF = [(Hi(0,0),1),(Wi(1,0),1),(Hi(0,1),-1),(Wi(0,0),-1)]
    gC = [(Hi(0,0),1),(Hi(1,0),1),(Hi(2,0),1)]
    mk = lambda g: [ (E[e][0],E[e][1],e,1) if s>0 else (E[e][1],E[e][0],e,-1) for (e,s) in g ]
    return C("B0b", 9, E, mk(gF), mk(gC))

def T_edge(K, walk, a):
    """COR-F's tick: (Ts)(dst) = U_e^{+/-1} s(src) along the loop; identity off it."""
    U = np.exp(1j*np.asarray(a, float)); T = np.zeros((K.nv,K.nv), complex)
    on = {u for u,_,_,_ in walk}
    for v in range(K.nv):
        if v not in on: T[v,v] = 1.0
    for (u,v,e,s) in walk: T[v,u] = U[e] if s>0 else np.conj(U[e])
    return T

def holo(walk, a):
    z = 1+0j
    for (_,_,e,s) in walk: z *= np.exp(1j*a[e]) if s>0 else np.exp(-1j*a[e])
    return z

def M_circ(K, vs, W):
    M = np.eye(K.nv, dtype=complex)
    for v in vs: M[v,v] = W
    return M

def D_root(K, vs, W, L, branch=0):
    """THE FIBRE-WISE L-th ROOT OF M_gamma.  D = diag(W^{1/L}) on the loop, identity off it.
       D is UNITARY, FIBRE-WISE (so it lies in the local gauge group U(1)^V, W-06's corrected N4
       mechanism), GAUGE-INVARIANT at every power, and D^L = M_gamma EXACTLY.
       It is the same KIND of object as the corpus's own M_gamma -- scalar multiplication of the
       loop's fibres -- ticking once per edge instead of once per circuit."""
    th = (np.angle(W) + 2*np.pi*branch)/L
    return M_circ(K, vs, np.exp(1j*th))

def pi_of(K, s):
    w = np.abs(s)**2; p = {(0,0):0.,(1,0):0.,(0,1):0.,(1,1):0.}
    for v in range(K.nv): p[K.cls(v)] += w[v]
    return np.array([p[(0,0)],p[(1,0)],p[(0,1)],p[(1,1)]])

def states_same_pi(K, pi, n, rng, phases=True):
    cl = {}
    for v in range(K.nv): cl.setdefault(K.cls(v), []).append(v)
    order = [(0,0),(1,0),(0,1),(1,1)]; out = []
    for _ in range(n):
        w = np.zeros(K.nv)
        for idx,c in enumerate(order):
            vs = cl.get(c, [])
            if not vs:
                assert abs(pi[idx]) < 1e-14; continue
            x = rng.dirichlet(np.ones(len(vs)))
            for v,xv in zip(vs,x): w[v] = pi[idx]*xv
        ph = rng.uniform(0,2*np.pi,K.nv) if phases else np.zeros(K.nv)
        out.append(np.sqrt(w)*np.exp(1j*ph))
    return out

def m_jensen(p, n=1<<22):
    a,b,c,d = p
    t = 2*np.pi*np.arange(n)/n; ct = np.cos(t)
    A = np.sqrt(np.maximum(a*a+b*b+2*a*b*ct, 0.0))
    B = np.sqrt(np.maximum(c*c+d*d+2*c*d*ct, 0.0))
    return float(np.log(np.maximum(np.maximum(A,B),1e-300)).mean())

def generic_conn(K, rng, f=1.0, c=2**0.5):
    """(f,c) = (1.0, sqrt2): the ONLY generic connection the corpus publishes (W-10 N-4; S4:603)."""
    a = rng.uniform(0.2, 2.0, len(K.edges))
    setC = {e for (_,_,e,_) in K.wC}; setF = {e for (_,_,e,_) in K.wF}
    eF = [e for (_,_,e,_) in K.wF if e not in setC]
    eC = [e for (_,_,e,_) in K.wC if e not in setF]
    ph = lambda w: sum(a[e] if s>0 else -a[e] for (_,_,e,s) in w)
    sF = [s for (_,_,e,s) in K.wF if e == eF[0]][0]; a[eF[0]] += sF*(f-ph(K.wF))
    sC = [s for (_,_,e,s) in K.wC if e == eC[0]][0]; a[eC[0]] += sC*(c-ph(K.wC))
    return a

def spread_over(states, AF, AC):
    v = np.array([abs(np.vdot(AF@s, AC@s)) for s in states])
    return v.max()-v.min()
