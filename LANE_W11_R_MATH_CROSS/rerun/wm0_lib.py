# LANE_W11_R_MATH — shared library.  Double precision is the default throughout this lane;
# every precision-sensitive claim is re-checked in exact arithmetic in wm5.
#
# CONVENTIONS (published, identical to LANE_W11_CONVENTION_TEST/PUBLISHED_CONVENTIONS.txt):
#   CIRCUIT convention (corpus's, S3 sec2.3 / W-01):
#       (M_gamma s)(v) = W(gamma) s(v) for v on the loop, s(v) otherwise.  DIAGONAL.
#   EDGE convention (COR-F's, S3_THE_CROSSING_AUDIT_V001.md:160-209, :794):
#       (T_gamma s)(dst) = U_e s(src) for each directed edge (src,dst,e) of the loop walk,
#       identity off the loop.  T^L = M_gamma, L = |loop|.
#   Observable  Z = <branch_F, branch_C>  with <z,w> = sum_v conj(z_v) w_v.
import numpy as np
from math import gcd

# ---------------------------------------------------------------- carriers
def K1():
    """S1 sec1.  5 vertices, 6 edges.  gamma_F = e1e2e3 (filled), gamma_C = e4e5e6 (unfilled)."""
    NV = 5
    walkF = [(0,1,0,+1),(1,2,1,+1),(2,0,2,+1)]          # (src,dst,edge,sign)
    walkC = [(0,3,3,+1),(3,4,4,+1),(4,0,5,+1)]
    return dict(name="K1", NV=NV, NE=6, walkF=walkF, walkC=walkC)

def B0b():
    """3x3 ring torus, loops meet.  Incidence verbatim from
       LANE_W10_A_CARRIERS/w10a_lib.py:182-213 (sealed).  |gamma_F|=4, |gamma_C|=3."""
    V = lambda i,j: 3*(j%3)+(i%3)
    H = lambda i,j: 3*(j%3)+(i%3)
    Wg= lambda i,j: 9+3*(j%3)+(i%3)
    E=[None]*18
    for j in range(3):
        for i in range(3):
            E[H(i,j)]=(V(i,j),V(i+1,j)); E[Wg(i,j)]=(V(i,j),V(i,j+1))
    gF=[(H(0,0),1),(Wg(1,0),1),(H(0,1),-1),(Wg(0,0),-1)]
    gC=[(H(0,0),1),(H(1,0),1),(H(2,0),1)]
    def walk(g):
        out=[]
        for (e,s) in g:
            u,v=E[e]
            out.append((u,v,e,+1) if s>0 else (v,u,e,-1))
        return out
    return dict(name="B0b", NV=9, NE=18, walkF=walk(gF), walkC=walk(gC), E=E)

# ---------------------------------------------------------------- operators
def loop_vs(walk):  return {u for u,_,_,_ in walk}

def Top(walk, a, NV):
    """COR-F's edge tick."""
    T=np.zeros((NV,NV),dtype=complex); on=loop_vs(walk)
    for v in range(NV):
        if v not in on: T[v,v]=1.0
    for (u,v,e,sg) in walk:
        T[v,u]= np.exp(1j*a[e]) if sg>0 else np.exp(-1j*a[e])
    return T

def hol(walk, a):
    z=1.0+0j
    for (u,v,e,sg) in walk: z *= np.exp(1j*a[e]) if sg>0 else np.exp(-1j*a[e])
    return z

def Mop(walk, a, NV):
    """The corpus's whole-circuit operator."""
    W=hol(walk,a); M=np.eye(NV,dtype=complex)
    for v in loop_vs(walk): M[v,v]=W
    return M

def classes(car):
    """(is v in gamma_F?, is v in gamma_C?) -> index 0=00,1=10,2=01,3=11 in the corpus's order."""
    F,C=loop_vs(car["walkF"]),loop_vs(car["walkC"])
    idx={(0,0):0,(1,0):1,(0,1):2,(1,1):3}
    return np.array([idx[(int(v in F),int(v in C))] for v in range(car["NV"])]), F, C

def pi_of(s, cl):
    w=np.abs(s)**2; p=np.zeros(4)
    for v,c in enumerate(cl): p[c]+=w[v]
    return p

# ---------------------------------------------------------------- Mahler measure
def m_poly(c, n=1<<20):
    """log-Mahler measure of  c0 + c1 X + c2 Y + c3 XY  for COMPLEX c, by Jensen in Y:
       m = (1/2pi) int log max(|c0+c1 X|,|c2+c3 X|) d(arg X).      Exact identity, one quadrature."""
    c0,c1,c2,c3=[complex(z) for z in c]
    t=2*np.pi*np.arange(n)/n; X=np.exp(1j*t)
    A=np.abs(c0+c1*X); B=np.abs(c2+c3*X)
    return np.log(np.maximum(np.maximum(A,B),1e-300)).mean()

def m_jensen_registrar(p, n=1<<20):
    """The registrar's / N1's own reduction, for real non-negative p.  Kept verbatim for
       arithmetic cross-check against m_poly."""
    a,b,c,d=p; t=2*np.pi*np.arange(n)/n; ct=np.cos(t)
    A=np.sqrt(np.maximum(a*a+b*b+2*a*b*ct,0)); B=np.sqrt(np.maximum(c*c+d*d+2*c*d*ct,0))
    return np.log(np.maximum(np.maximum(A,B),1e-300)).mean()
