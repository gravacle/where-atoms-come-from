# W-11 RULING LANE — shared library. Built from the sealed bytes, not from any lane's code.
# Carriers: K1 from S1_CARRIER_K1_V001.md sec1 (:16-22, :24-25).
#           B0b: 3x3 ring torus, gamma_F = one face boundary (len 4), gamma_C = row j=0 (len 3);
#           class multiset checked against S4_THE_MEASUREMENT_V001.md:575 {00:4,01:1,10:2,11:2}.
# CONVENTIONS. Double precision is the default everywhere unless a leg says "Fraction".
# Seed 20260818 -- deliberately NOT the registrar's 20260817 nor any refuter's.
import numpy as np
from collections import Counter

SEED = 20260818

# ---------------------------------------------------------------- carriers
class Carrier:
    """A carrier is: n vertices, an edge list, and two loops given as DIRECTED VERTEX WALKS
    with the edge index and traversal sign at each step.  Everything else is derived."""
    def __init__(self, name, nv, edges, walkF, walkC):
        self.name, self.nv, self.edges = name, nv, edges
        self.walkF, self.walkC = walkF, walkC          # [(u,v,e,sgn)] u->v using edge e
        self.VF = {u for u, _, _, _ in walkF}
        self.VC = {u for u, _, _, _ in walkC}
        self.LF, self.LC = len(walkF), len(walkC)
        self.cls = [(int(v in self.VF), int(v in self.VC)) for v in range(nv)]
        self.CLASSES = [(0,0),(1,0),(0,1),(1,1)]       # p00,p10,p01,p11 order, as W-10/N1
        self.idx = {c: [v for v in range(nv) if self.cls[v] == c] for c in self.CLASSES}
    def multiset(self):
        return dict(Counter(f"{a}{b}" for a, b in self.cls))
    def pi_of(self, s):
        w = np.abs(s)**2
        return np.array([sum(w[v] for v in self.idx[c]) for c in self.CLASSES])

def K1():
    E = [(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)]                       # e1..e6, S1 :19-21
    wF = [(0,1,0,1),(1,2,1,1),(2,0,2,1)]                            # filled triangle
    wC = [(0,3,3,1),(3,4,4,1),(4,0,5,1)]                            # unfilled cycle
    return Carrier("K1", 5, E, wF, wC)

def B0b():
    Vx = lambda i,j: 3*(j%3)+(i%3)
    H  = lambda i,j: 3*(j%3)+(i%3)
    Wv = lambda i,j: 9+3*(j%3)+(i%3)
    E = [None]*18
    for j in range(3):
        for i in range(3):
            E[H(i,j)]  = (Vx(i,j), Vx(i+1,j))
            E[Wv(i,j)] = (Vx(i,j), Vx(i,j+1))
    def step(e, sgn):
        u, v = E[e]
        return (u,v,e,+1) if sgn > 0 else (v,u,e,-1)
    wF = [step(H(0,0),1), step(Wv(1,0),1), step(H(0,1),-1), step(Wv(0,0),-1)]   # face bdry, len 4
    wC = [step(H(0,0),1), step(H(1,0),1),  step(H(2,0),1)]                       # row j=0, len 3
    return Carrier("B0b", 9, E, wF, wC)

# ---------------------------------------------------------------- operators
def Tedge(car, walk, a):
    """COR-F's EDGE tick: move each fibre value one edge along the loop, identity off it.
    S1 :51-53 -- transport along e:u->v is z|->U_e z; reverse traversal by U_e^{-1}."""
    U = np.exp(1j*np.asarray(a, dtype=float))
    T = np.zeros((car.nv, car.nv), dtype=complex)
    on = {u for u,_,_,_ in walk}
    for v in range(car.nv):
        if v not in on: T[v,v] = 1.0
    for (u,v,e,sg) in walk:
        T[v,u] = U[e] if sg > 0 else np.conj(U[e])
    return T

def holon(walk, a):
    z = 1.0+0j
    for (_,_,e,sg) in walk:
        z *= np.exp(1j*a[e]) if sg > 0 else np.exp(-1j*a[e])
    return z

def Mcirc(car, Vs, W):
    """W-01's CIRCUIT operator: multiply s_v by the whole-circuit holonomy on the loop."""
    M = np.eye(car.nv, dtype=complex)
    for v in Vs: M[v,v] = W
    return M

def Droot(car, Vs, W, L, branch=0):
    """The FIBRE-WISE L-th root of M_gamma: diag(W^{1/L}) on the loop, identity off it.
    Requires a BRANCH of the L-th root -- a lift the carrier does not supply."""
    w = np.exp(1j*(np.angle(W) + 2*np.pi*branch)/L)
    return Mcirc(car, Vs, w)

# ---------------------------------------------------------------- observable
def Z(opF, opC, s, n):
    return np.vdot(np.linalg.matrix_power(opF,n)@s, np.linalg.matrix_power(opC,n)@s)

def Qrel(opF, opC, n):
    """The RELATIVE branch operator at tick n:  Q_n = (opF^n)^*  opC^n.   Z_n = <s, Q_n s>."""
    A = np.linalg.matrix_power(opF,n); B = np.linalg.matrix_power(opC,n)
    return A.conj().T @ B

def is_class_constant_diag(car, Q, tol=1e-9):
    """Is Q multiplication by a function of the INCIDENCE CLASS?"""
    off = np.linalg.norm(Q - np.diag(np.diag(Q)))
    if off > tol: return False, off
    d = np.diag(Q); worst = 0.0
    for c in car.CLASSES:
        vs = car.idx[c]
        if len(vs) >= 2:
            worst = max(worst, float(np.max(np.abs(d[vs] - d[vs[0]]))))
    return worst <= tol, max(off, worst)

# ---------------------------------------------------------------- same-pi state families
def same_pi_states(car, rng, base_w, m):
    """m states with EXACTLY the base class sums, differing only in the within-class split
    and in phase.  pi is held fixed BY CONSTRUCTION; nothing else is."""
    out = []
    for _ in range(m):
        w = np.zeros(car.nv)
        for c in car.CLASSES:
            vs = car.idx[c]
            if not vs: continue
            tot = sum(base_w[v] for v in vs)
            if len(vs) == 1: w[vs[0]] = tot
            else:
                x = rng.dirichlet(np.ones(len(vs)))
                for k, v in enumerate(vs): w[v] = tot*x[k]
        ph = rng.uniform(0, 2*np.pi, car.nv)
        out.append(np.sqrt(w)*np.exp(1j*ph))
    return out

def pi_spread(car, opF, opC, states, ns):
    """max over ticks in ns of the spread of |Z_n| across states of IDENTICAL pi."""
    worst = 0.0
    for n in ns:
        v = [abs(Z(opF,opC,s,n)) for s in states]
        worst = max(worst, max(v)-min(v))
    return worst

# ---------------------------------------------------------------- rates
def rate(opF, opC, s, N):
    xF = s.copy(); xC = s.copy(); tot = 0.0
    for _ in range(N):
        xF = opF@xF; xC = opC@xC
        z = abs(np.vdot(xF, xC))
        tot += np.log(z) if z > 0 else -700.0
    return tot/N

def mahler4(p, n=1<<20):
    """m(p00 + p10 x + p01 y + p11 xy) by the Jensen reduction, non-negative coefficients."""
    a,b,c,d = p
    t = 2*np.pi*np.arange(n)/n; ct = np.cos(t)
    A = np.sqrt(np.maximum(a*a+b*b+2*a*b*ct, 0)); B = np.sqrt(np.maximum(c*c+d*d+2*c*d*ct, 0))
    return float(np.log(np.maximum(A,B)+1e-300).mean())

def a_generic(car, rng, f, c):
    """A connection with holonomy exp(i f) on gamma_F and exp(i c) on gamma_C, built by
    solving for the LAST step of each walk with the correct SIGN for its traversal direction."""
    a = rng.uniform(0, 2*np.pi, len(car.edges))
    for (walk, target) in ((car.walkF, f), (car.walkC, c)):
        (_,_,elast,sglast) = walk[-1]
        a[elast] = 0.0
        rest = np.angle(holon(walk, a))
        a[elast] = ((target - rest) if sglast > 0 else (rest - target)) % (2*np.pi)
    return a
