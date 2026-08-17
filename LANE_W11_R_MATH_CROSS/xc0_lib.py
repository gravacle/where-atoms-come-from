# LANE_W11_R_MATH_CROSS — shared library.  WRITTEN INDEPENDENTLY of LANE_W11_R_MATH/wm0_lib.py:
# the operators are built here from the walk lists directly, and the Mahler measure is computed
# by an EXACT CLOSED FORM (Cassaigne-Maillot, with the Bloch-Wigner dilogarithm summed from its
# own series) as well as by quadrature, so the two paths share no code.
#
# CONVENTIONS: identical to LANE_W11_CONVENTION_TEST/PUBLISHED_CONVENTIONS.txt.
#   CIRCUIT  (M s)(v) = W(gamma) s(v) on the loop, s(v) off it.  DIAGONAL.
#   EDGE     (T s)(dst) = U_e s(src) along the loop walk, identity off it.  T^L = M.
#   Z = <bF, bC>, <z,w> = sum conj(z) w.
# DOUBLE PRECISION IS THE DEFAULT.  Exact blocks are marked and use Fraction only.
import numpy as np
from math import gcd, log, pi as PI, atan2, cos, sin, sqrt

# ---------------------------------------------------------------- carriers (re-derived here)
def K1():
    return dict(name="K1", NV=5, NE=6,
                walkF=[(0,1,0,+1),(1,2,1,+1),(2,0,2,+1)],
                walkC=[(0,3,3,+1),(3,4,4,+1),(4,0,5,+1)])

def B0b():
    # 3x3 torus grid.  h(i,j) index 3j+i ; w(i,j) index 9+3j+i ; vertex (i,j) -> 3j+i.
    V=lambda i,j: 3*(j%3)+(i%3)
    E=[None]*18
    for j in range(3):
        for i in range(3):
            E[3*j+i]=(V(i,j),V(i+1,j)); E[9+3*j+i]=(V(i,j),V(i,j+1))
    gF=[(0,+1),(9+1,+1),(3,-1),(9,-1)]          # h(0,0), w(1,0), -h(0,1), -w(0,0)
    gC=[(0,+1),(1,+1),(2,+1)]                   # h(0,0) h(1,0) h(2,0)
    def walk(g):
        out=[]
        for e,s in g:
            u,v=E[e]; out.append((u,v,e,+1) if s>0 else (v,u,e,-1))
        return out
    return dict(name="B0b", NV=9, NE=18, walkF=walk(gF), walkC=walk(gC), E=E)

def SHARE2():
    """TWO 4-CYCLES SHARING ONE EDGE.  Built here for the case LANE_W11_R_MATH's Theorem U4
       declares uncovered: |class 11| >= 2 AND gcd(L_F,L_C) > 1.
       vertices 0..6 ; gamma_F = 0->1->2->3->0 ; gamma_C = 0->1->4->5->0 ; v6 is a spectator.
       edges: 0:(0,1) shared, 1:(1,2), 2:(2,3), 3:(3,0), 4:(1,4), 5:(4,5), 6:(5,0), 7:(6,0)."""
    E=[(0,1),(1,2),(2,3),(3,0),(1,4),(4,5),(5,0),(6,0)]
    walkF=[(0,1,0,+1),(1,2,1,+1),(2,3,2,+1),(3,0,3,+1)]
    walkC=[(0,1,0,+1),(1,4,4,+1),(4,5,5,+1),(5,0,6,+1)]
    return dict(name="SHARE2", NV=7, NE=8, walkF=walkF, walkC=walkC, E=E)

# ---------------------------------------------------------------- operators
def loop_vs(walk): return {u for u,_,_,_ in walk}

def Top(walk,a,NV):
    T=np.zeros((NV,NV),dtype=complex); on=loop_vs(walk)
    for v in range(NV):
        if v not in on: T[v,v]=1.0
    for (u,v,e,sg) in walk: T[v,u]=np.exp(1j*a[e]) if sg>0 else np.exp(-1j*a[e])
    return T

def hol(walk,a):
    z=1.0+0j
    for (u,v,e,sg) in walk: z*= np.exp(1j*a[e]) if sg>0 else np.exp(-1j*a[e])
    return z

def Mop(walk,a,NV):
    W=hol(walk,a); M=np.eye(NV,dtype=complex)
    for v in loop_vs(walk): M[v,v]=W
    return M

def classes(car):
    F,C=loop_vs(car["walkF"]),loop_vs(car["walkC"])
    idx={(0,0):0,(1,0):1,(0,1):2,(1,1):3}
    return np.array([idx[(int(v in F),int(v in C))] for v in range(car["NV"])]),F,C

def pi_of(s,cl):
    w=np.abs(s)**2; p=np.zeros(4)
    for v,c in enumerate(cl): p[c]+=w[v]
    return p

def same_pi_states(cl,pi,rng,k):
    NV=len(cl); out=[]
    for _ in range(k):
        w=np.zeros(NV)
        for c in range(4):
            idx=np.where(cl==c)[0]
            if len(idx)==0: continue
            q=rng.random(len(idx)); q=q/q.sum()*pi[c]; w[idx]=q
        out.append(np.sqrt(w)*np.exp(1j*rng.uniform(0,2*np.pi,NV)))
    return out

def spread(opF,opC,n,states):
    v=[abs(np.vdot(np.linalg.matrix_power(opF,n)@s,np.linalg.matrix_power(opC,n)@s)) for s in states]
    return max(v)-min(v)

# ------------------------------------------------- Mahler measure, TWO INDEPENDENT PATHS
def m_quad(c,n=1<<20):
    """Jensen in Y, uniform grid on the circle (same identity the lane under test uses,
       written here from scratch)."""
    c0,c1,c2,c3=[complex(z) for z in c]
    t=2*np.pi*np.arange(n)/n; X=np.exp(1j*t)
    return float(np.log(np.maximum(np.maximum(np.abs(c0+c1*X),np.abs(c2+c3*X)),1e-300)).mean())

from fractions import Fraction as _Fr
def _bernoulli(N=44):
    B=[_Fr(0)]*(N+1); B[0]=_Fr(1)
    from math import comb
    for m in range(1,N+1):
        s=_Fr(0)
        for k in range(m): s+=_Fr(comb(m+1,k))*B[k]
        B[m]=-s/_Fr(m+1)
    return [float(b) for b in B]
_BERN=_bernoulli(44)
_FACT=[1.0]*46
for _i in range(1,46): _FACT[_i]=_FACT[_i-1]*_i

def _Li2(z):
    """dilogarithm by the Bernoulli series in u = -log(1-z) (converges for |u| < 2 pi),
       with the inversion Li2(z) = -Li2(1/z) - pi^2/6 - (1/2) log(-z)^2 for |z| > 1.
       Written from the definition; shares no code with the lane under test."""
    z=complex(z)
    if abs(z-1.0)<1e-15: return complex(PI**2/6)
    if abs(z)<1e-18: return 0j
    if abs(z)>1.0:
        return -_Li2(1.0/z)-PI**2/6-0.5*np.log(-z)**2
    u=-np.log(1.0-z); s=0j; up=1.0+0j
    for k in range(0,44):
        up=up*u
        s+=_BERN[k]*up/_FACT[k+1]
    return s

def bloch_wigner(z):
    z=complex(z)
    if abs(z)<1e-300: return 0.0
    if abs(z)<=1.0: return float(np.imag(_Li2(z))+np.angle(1-z)*np.log(abs(z)))
    # D(1/z) = -D(z)
    return -bloch_wigner(1.0/z)

def m_CM(a,b,c):
    """EXACT Cassaigne-Maillot closed form for m(a + b x + c y), a,b,c >= 0.
       If the triangle inequality fails, m = log max(a,b,c);  else
       pi*m = alpha log a + beta log b + gamma log c + D( (a/b) e^{i gamma} )
       with alpha,beta,gamma the angles opposite the sides a,b,c."""
    a,b,c=float(a),float(b),float(c)
    if min(a,b,c)<0: raise ValueError
    if a>b+c or b>a+c or c>a+b or min(a,b,c)==0.0:
        return log(max(a,b,c))
    al=np.arccos((b*b+c*c-a*a)/(2*b*c))     # angle opposite side a
    be=np.arccos((a*a+c*c-b*b)/(2*a*c))
    ga=np.arccos((a*a+b*b-c*c)/(2*a*b))
    z=(a/b)*np.exp(1j*ga)
    return (al*log(a)+be*log(b)+ga*log(c)+bloch_wigner(z))/PI
