# LANE W-11-R-T — shared library.  Conventions in PUBLISHED_CONVENTIONS.txt.  Double precision.
import numpy as np

# ---------------- carriers ----------------
class Carrier:
    def __init__(self, name, NV, edges, loopF, loopC):
        self.name, self.NV, self.edges = name, NV, edges
        self.loopF, self.loopC = loopF, loopC          # list of (src,dst,edge_index) in cyclic order
        self.VF = [s for s,_,_ in loopF]; self.VC = [s for s,_,_ in loopC]
        self.LF, self.LC = len(loopF), len(loopC)
        self.cls = [ (int(v in self.VF), int(v in self.VC)) for v in range(NV) ]
        self.classes = {}
        for v in range(NV): self.classes.setdefault(self.cls[v], []).append(v)
    def hol(self, a, which):
        loop = self.loopF if which=='F' else self.loopC
        z = 1.0+0j
        for (_,_,e) in loop: z *= np.exp(1j*a[e])
        return z
    def M(self, a, which):
        vs = self.VF if which=='F' else self.VC
        Mm = np.eye(self.NV, dtype=complex); W = self.hol(a, which)
        for v in vs: Mm[v,v] = W
        return Mm
    def T_corf(self, a, which):
        """COR-F's edge tick: (Ts)(dst) = U_e s(src) along the loop, identity off it."""
        loop = self.loopF if which=='F' else self.loopC
        on = {s for s,_,_ in loop}
        T = np.zeros((self.NV,self.NV), dtype=complex)
        for v in range(self.NV):
            if v not in on: T[v,v] = 1.0
        for (s_,d_,e) in loop: T[d_,s_] = np.exp(1j*a[e])
        return T

# K1 — S1 sec1 verbatim
K1 = Carrier("K1", 5, [(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)],
             [(0,1,0),(1,2,1),(2,0,2)], [(0,3,3),(3,4,4),(4,0,5)])

# B0b — incidence verbatim from LANE_W10_A_CARRIERS/w10a_lib.py:182-213 as reproduced in
# LANE_W11_CONVENTION_TEST/w11_c_b0b.py; re-expressed here as directed vertex walks.
def _b0b():
    V = lambda i,j: 3*(j%3)+(i%3); H = lambda i,j: 3*(j%3)+(i%3); Wg = lambda i,j: 9+3*(j%3)+(i%3)
    E = [None]*18
    for j in range(3):
        for i in range(3):
            E[H(i,j)] = (V(i,j),V(i+1,j)); E[Wg(i,j)] = (V(i,j),V(i,j+1))
    gF = [(H(0,0),1),(Wg(1,0),1),(H(0,1),-1),(Wg(0,0),-1)]
    gC = [(H(0,0),1),(H(1,0),1),(H(2,0),1)]
    def walk(g):
        out = []
        for (e,s) in g:
            u,v = E[e]; out.append((u,v,e) if s>0 else (v,u,e))
        return out
    return E, walk(gF), walk(gC)
_E, _wF, _wC = _b0b()
class B0bCarrier(Carrier):
    """B0b traverses one edge of gamma_F backwards, so the edge phase is conjugated there."""
    def __init__(self):
        super().__init__("B0b", 9, _E, _wF, _wC)
        self.signF = [1,1,-1,-1]; self.signC = [1,1,1]
    def hol(self, a, which):
        loop, sg = (self.loopF,self.signF) if which=='F' else (self.loopC,self.signC)
        z = 1.0+0j
        for (_,_,e),s in zip(loop,sg): z *= np.exp(1j*a[e]*s)
        return z
    def T_corf(self, a, which):
        loop, sg = (self.loopF,self.signF) if which=='F' else (self.loopC,self.signC)
        on = {s_ for s_,_,_ in loop}
        T = np.zeros((self.NV,self.NV), dtype=complex)
        for v in range(self.NV):
            if v not in on: T[v,v] = 1.0
        for (s_,d_,e),sg_ in zip(loop,sg): T[d_,s_] = np.exp(1j*a[e]*sg_)
        return T
B0b = B0bCarrier()

# SYNTH-D — NOT a corpus carrier.  Two 4-cycles sharing the two vertices {0,1}, plus a spectator.
# Built ONLY to exhibit the boundary of this lane's own theorem (|J| >= 2 AND gcd(L_F,L_C) >= 2).
SYNTHD = Carrier("SYNTH-D", 7, [(0,1),(1,2),(2,3),(3,0),(1,4),(4,5),(5,0)],
                 [(0,1,0),(1,2,1),(2,3,2),(3,0,3)], [(0,1,0),(1,4,4),(4,5,5),(5,0,6)])

# ---------------- the family ----------------
def haar(n, rng):
    z = (rng.normal(size=(n,n)) + 1j*rng.normal(size=(n,n)))/np.sqrt(2)
    q,r = np.linalg.qr(z); return q*(np.diag(r)/abs(np.diag(r)))

def member(car, a, which, jvec, V):
    """T_S = V diag(rho zeta^j) V^*, identity off the loop.  ONTO the family (see leg T1.3)."""
    vs = car.VF if which=='F' else car.VC
    L  = len(vs); W = car.hol(a, which)
    rho = np.exp(1j*np.angle(W)/L); zt = np.exp(2j*np.pi/L)
    A = V@np.diag(rho*zt**np.asarray(jvec))@V.conj().T
    T = np.eye(car.NV, dtype=complex)
    for i,u in enumerate(vs):
        for k,w in enumerate(vs): T[u,w] = A[i,k]
    return T

def diag_member(car, a, which, jvec):
    return member(car, a, which, jvec, np.eye(len(car.VF if which=='F' else car.VC), dtype=complex))

def shift_member(car, a, which, phases):
    """branch (B): a unimodular weighted cyclic shift with prod n_i = W.  phases has L-1 free angles."""
    loop = car.loopF if which=='F' else car.loopC
    L = len(loop); W = car.hol(a, which)
    n = np.exp(1j*np.asarray(phases)); n = np.append(n, W/np.prod(n))
    on = {s_ for s_,_,_ in loop}
    T = np.zeros((car.NV,car.NV), dtype=complex)
    for v in range(car.NV):
        if v not in on: T[v,v] = 1.0
    for k,(s_,d_,_) in enumerate(loop): T[d_,s_] = n[k]
    return T

# ---------------- THE INVISIBILITY INSTRUMENT ----------------
def Dscore(car, TF, TC, nmax=12):
    """D = max_n ( ||offdiag Q_n||_F + max_class spread of diag Q_n ),  Q_n = (T_F^n)^* T_C^n.
       Z_n(s) depends on s only through pi FOR ALL s  <=>  D = 0.  (Proof in PUBLISHED_CONVENTIONS.)"""
    worst = 0.0; XF = np.eye(car.NV,dtype=complex); XC = np.eye(car.NV,dtype=complex)
    for n in range(1, nmax+1):
        XF = TF@XF; XC = TC@XC
        Q = XF.conj().T@XC
        off = np.linalg.norm(Q - np.diag(np.diag(Q)))
        d = np.diag(Q); spread = 0.0
        for c,vs in car.classes.items():
            if len(vs) > 1:
                spread = max(spread, max(abs(d[u]-d[w]) for u in vs for w in vs))
        worst = max(worst, off + spread)
    return worst

def Dper_n(car, TF, TC, nmax=12):
    out = []; XF = np.eye(car.NV,dtype=complex); XC = np.eye(car.NV,dtype=complex)
    for n in range(1, nmax+1):
        XF = TF@XF; XC = TC@XC
        Q = XF.conj().T@XC
        off = np.linalg.norm(Q - np.diag(np.diag(Q)))
        d = np.diag(Q); spread = 0.0
        for c,vs in car.classes.items():
            if len(vs) > 1:
                spread = max(spread, max(abs(d[u]-d[w]) for u in vs for w in vs))
        out.append(off+spread)
    return out

def pi_of(car, s):
    w = np.abs(s)**2; p = {}
    for v in range(car.NV): p[car.cls[v]] = p.get(car.cls[v],0.0)+w[v]
    return np.array([p.get((0,0),0.0), p.get((1,0),0.0), p.get((0,1),0.0), p.get((1,1),0.0)])

def a_with_holonomies(car, f, c, rng):
    """draw edge phases uniformly, then fix ONE edge of each loop so the holonomies are exactly
       (f,c).  Gauge-inequivalent connections with the same (f,c) do not exist on these carriers."""
    a = rng.uniform(0,2*np.pi,len(car.edges))
    lastF = car.loopF[-1][2]; a[lastF] = 0.0
    a[lastF] = (f - np.angle(car.hol(a,'F'))) % (2*np.pi)
    lastC = car.loopC[-1][2]
    if lastC != lastF:
        a[lastC] = 0.0; a[lastC] = (c - np.angle(car.hol(a,'C'))) % (2*np.pi)
    return a
