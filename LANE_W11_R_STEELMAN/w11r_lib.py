# LANE W11-R (STEELMAN READING A) -- shared library.
# Everything here is rebuilt from the SEALED ARTIFACT BYTES, not from LANE_W11_CONVENTION_TEST's code.
#   K1 : S1_CARRIER_K1_V001.md sec1 (5 vertices, 6 edges, faces, root v0).
#   B0b: incidence identical to LANE_W10_A_CARRIERS/w10a_lib.py:182-213 (sealed; reproduced there
#        against S4:512/:539/:575 by a build lane and two refuters).  I re-derive the walk myself
#        below and CHECK the class multiset against S4:575 rather than importing the code.
# Double precision is the default.  Exact checks are done with fractions.Fraction where labelled.
import numpy as np

# ---------------------------------------------------------------- K1, from S1 sec1
K1_EDGES = [(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)]        # e1..e6, each source->target
K1_LOOP_F = [(0,1,0),(1,2,1),(2,0,2)]                   # gamma_F = e1.e2.e3 (the FILLED triangle)
K1_LOOP_C = [(0,3,3),(3,4,4),(4,0,5)]                   # gamma_C = e4.e5.e6 (the UNFILLED cycle)
K1_NV = 5

# ---------------------------------------------------------------- B0b, 3x3 torus grid
def _b0b():
    V = lambda i,j: 3*(j%3)+(i%3)
    H = lambda i,j: 3*(j%3)+(i%3)
    W = lambda i,j: 9+3*(j%3)+(i%3)
    E = [None]*18
    for j in range(3):
        for i in range(3):
            E[H(i,j)] = (V(i,j), V(i+1,j))
            E[W(i,j)] = (V(i,j), V(i,j+1))
    gF = [(H(0,0),1),(W(1,0),1),(H(0,1),-1),(W(0,0),-1)]   # face (0,0) boundary, length 4
    gC = [(H(0,0),1),(H(1,0),1),(H(2,0),1)]                # horizontal row j=0, length 3
    def walk(g):
        out=[]
        for (e,s) in g:
            u,v = E[e]
            out.append((u,v,e,+1) if s>0 else (v,u,e,-1))
        return out
    def as_loop(g):
        # convert a signed edge cycle into the (src,dst,edge,sign) walk used by Top/Mop
        return walk(g)
    return E, as_loop(gF), as_loop(gC), 9

B0B_E, B0B_LOOP_F, B0B_LOOP_C, B0B_NV = _b0b()

# ---------------------------------------------------------------- operators
def _norm_loop(loop):
    """accept (src,dst,edge) triples (K1) or (src,dst,edge,sign) quadruples (B0b)."""
    return [(t[0],t[1],t[2], t[3] if len(t)==4 else +1) for t in loop]

def loop_vertices(loop):
    return {t[0] for t in _norm_loop(loop)}

def T_edge(loop, a, NV):
    """COR-F's EDGE TICK: move each fibre value one edge along the loop, identity off it.
       (T s)(dst) = U_e^{+-1} s(src).   S3_THE_CROSSING_AUDIT_V001.md:160-209."""
    L = _norm_loop(loop); on = loop_vertices(loop)
    T = np.zeros((NV,NV), dtype=complex)
    for v in range(NV):
        if v not in on: T[v,v] = 1.0
    for (src,dst,e,sg) in L:
        T[dst,src] = np.exp(1j*a[e]) if sg>0 else np.exp(-1j*a[e])
    return T

def holonomy(loop, a):
    z = 1.0+0j
    for (src,dst,e,sg) in _norm_loop(loop):
        z *= np.exp(1j*a[e]) if sg>0 else np.exp(-1j*a[e])
    return z

def M_circuit(loop, a, NV):
    """THE CORPUS'S CONVENTION (W-01 / S3 sec2.3): multiply s_v by the WHOLE-CIRCUIT holonomy for
       every v on the loop, identity elsewhere.  Diagonal by construction."""
    W = holonomy(loop, a)
    M = np.eye(NV, dtype=complex)
    for v in loop_vertices(loop): M[v,v] = W
    return M

def D_uniform(loop, a, NV):
    """THE UNIFORM (PRINCIPAL) L-th ROOT: diag(W^{1/L} on the loop, 1 off).  D^L = M_circuit
       exactly.  Fibre-wise; equals the identity at the trivial connection."""
    L = len(_norm_loop(loop))
    W = holonomy(loop, a)
    root = np.exp(1j*np.angle(W)/L)          # principal branch
    D = np.eye(NV, dtype=complex)
    for v in loop_vertices(loop): D[v,v] = root
    return D

# ---------------------------------------------------------------- classes / pushforward
def classes(loopF, loopC, NV):
    F, C = loop_vertices(loopF), loop_vertices(loopC)
    return [(int(v in F), int(v in C)) for v in range(NV)]

def pi_of(s, loopF, loopC, NV):
    cl = classes(loopF, loopC, NV); w = np.abs(s)**2
    p = {(0,0):0.0,(1,0):0.0,(0,1):0.0,(1,1):0.0}
    for v in range(NV): p[cl[v]] += w[v]
    return np.array([p[(0,0)], p[(1,0)], p[(0,1)], p[(1,1)]])

def m_jensen(p, n=1<<20):
    """logarithmic Mahler measure of p00 + p10 x + p01 y + p11 xy, by the Jensen reduction."""
    a,b,c,d = p
    t = 2*np.pi*np.arange(n)/n; ct = np.cos(t)
    A = np.sqrt(np.maximum(a*a+b*b+2*a*b*ct, 0)); B = np.sqrt(np.maximum(c*c+d*d+2*c*d*ct, 0))
    return np.log(np.maximum(A,B)+1e-300).mean()

def Z(opF, opC, s, nF, nC):
    """the observable <branch_F, branch_C>; nF/nC are the two branches' tick counts."""
    return np.vdot(np.linalg.matrix_power(opF,nF)@s, np.linalg.matrix_power(opC,nC)@s)

def rate(opF, opC, s, N, stepF=1, stepC=1):
    """(1/N) sum_{n<=N} log|Z| accumulated stepwise, no matrix powers, no overflow."""
    AF = np.linalg.matrix_power(opF, stepF); AC = np.linalg.matrix_power(opC, stepC)
    xF = s.copy(); xC = s.copy(); tot = 0.0
    for _ in range(N):
        xF = AF@xF; xC = AC@xC
        z = abs(np.vdot(xF,xC))
        tot += np.log(z) if z > 0 else -700.0
    return tot/N

def arms_differ(*arrays):
    """DIFF YOUR ARMS.  Returns True only if every pair of arms differs in BYTES."""
    bs = [np.ascontiguousarray(x).tobytes() for x in arrays]
    return all(bs[i] != bs[j] for i in range(len(bs)) for j in range(i+1, len(bs)))
