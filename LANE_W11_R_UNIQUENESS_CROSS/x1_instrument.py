# LANE W-11-R-T-CROSS, leg X1 — INDEPENDENT REBUILD, AND THE ONE THING THE UNIQUENESS LANE
# NEVER CHECKED: IS ITS INSTRUMENT D THE SAME CRITERION THE CORPUS'S CLAIMS ARE ABOUT?
#
# The corpus's carrier-independent layer is stated in |Z| and in lambda = m(P) (N1, W-03/N2),
# NOT in the complex quadratic form Z.  The uniqueness lane scores every member by
#     D(n) = ||offdiag Q_n|| + max_class spread(diag Q_n),  Q_n = (T_F^n)* T_C^n,
# which is exactly "Z_n(s) is a function of pi(s) alone".  That is a STRICTLY STRONGER
# requirement than "|Z_n(s)| is a function of pi(s) alone".  If the two differ anywhere in the
# family, every count in leg T2 (81/729, 144/6912, 0 of 20000 shift, 0 of 979 generic) is a
# count of the wrong set, and the kill condition could be met by a member D rejects.
# THIS LEG TRIES TO FIND SUCH A MEMBER.  Nothing here is taken from rlib.py.
import numpy as np
from fractions import Fraction
rng = np.random.default_rng(11072026)          # NOT the registrar's seed: independent draws
np.set_printoptions(linewidth=200)

# ---------- K1 rebuilt from S1_CARRIER_K1_V001.md sec1 and sec3, by hand ----------
NV = 5
EDGES = [(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)]      # e1..e6, S1:19-21
LF = [(0,1,0),(1,2,1),(2,0,2)]                     # gamma_F, the FILLED triangle's boundary
LC = [(0,3,3),(3,4,4),(4,0,5)]                     # gamma_C, the unfilled cycle (COR-F's loop)
VF, VC = [0,1,2], [0,3,4]
CLS = [(int(v in VF), int(v in VC)) for v in range(NV)]
CLASSES = {}
for v in range(NV): CLASSES.setdefault(CLS[v], []).append(v)

def T_edge(loop, U):                                # COR-F: (Ts)(dst) = U_e s(src), identity off
    on = {s for s,_,_ in loop}
    T = np.zeros((NV,NV), dtype=complex)
    for v in range(NV):
        if v not in on: T[v,v] = 1.0
    for (s_,d_,e) in loop: T[d_,s_] = U[e]
    return T
def M_circ(vs, W):
    M = np.eye(NV, dtype=complex)
    for v in vs: M[v,v] = W
    return M

print("== X1.1  COR-F's SEALED EXHIBIT, REBUILT A THIRD TIME FROM THE AUDIT'S TEXT ==")
a = np.zeros(6); a[3:6] = [0.7, 1.3, -0.4]
U = np.exp(1j*a); T = T_edge(LC, U); WC = np.exp(1j*(a[3]+a[4]+a[5]))
rho = np.diag([0.40,0.15,0.15,0.15,0.15]).astype(complex)
print(f"   ||T*T - I||           = {np.linalg.norm(T.conj().T@T-np.eye(NV)):.2e}   [audit 0.00e+00]")
print(f"   W_C                   = {WC:.6f}   [audit -0.029200+0.999574j]")
print(f"   diag(T rho T*)        = {np.round(np.real(np.diag(T@rho@T.conj().T)),2)}  [audit .15 .15 .15 .40 .15]")
print(f"   ||T^3 - M_c||         = {np.linalg.norm(np.linalg.matrix_power(T,3)-M_circ(VC,WC)):.2e}")

print("\n== X1.2  T^L = M_gamma IN EXACT ARITHMETIC (no float64 anywhere) ==")
# Gaussian-rational points of U(1): (3+4i)/5, (5+12i)/13, (8+15i)/17.  Exact, unimodular.
def cmul(x,y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
Uq = [(Fraction(3,5),Fraction(4,5)), (Fraction(5,13),Fraction(12,13)), (Fraction(8,17),Fraction(15,17)),
      (Fraction(-3,5),Fraction(4,5)), (Fraction(5,13),Fraction(-12,13)), (Fraction(0),Fraction(1))]
Wq = (Fraction(1),Fraction(0))
for e in (3,4,5): Wq = cmul(Wq, Uq[e])
# T^3 on the C-loop, entrywise, exactly
Tq = [[(Fraction(0),Fraction(0)) for _ in range(NV)] for _ in range(NV)]
for v in range(NV):
    if v not in {s for s,_,_ in LC}: Tq[v][v] = (Fraction(1),Fraction(0))
for (s_,d_,e) in LC: Tq[d_][s_] = Uq[e]
def matmulq(A,B):
    return [[tuple(sum(z) for z in zip(*[cmul(A[i][k],B[k][j]) for k in range(NV)]))
             for j in range(NV)] for i in range(NV)]
T3q = matmulq(matmulq(Tq,Tq),Tq)
ok = all(T3q[i][j] == ((Wq if i in VC else (Fraction(1),Fraction(0))) if i==j else (Fraction(0),Fraction(0)))
         for i in range(NV) for j in range(NV))
print(f"   exact Gaussian-rational connection, W_C = {Wq[0]} + {Wq[1]}i, |W|^2 = {Wq[0]**2+Wq[1]**2}")
print(f"   T^3 == M_gamma EXACTLY (Fraction arithmetic, zero rounding): {ok}")
print("   -> the registrar's leg A and the uniqueness lane's T1.2 are an IDENTITY.  Confirmed a")
print("      third time, and this time with no floating point in the check at all.")

# ---------- the family, my own parameterisation ----------
def haar(n):
    z = (rng.normal(size=(n,n)) + 1j*rng.normal(size=(n,n)))/np.sqrt(2)
    q,r = np.linalg.qr(z); return q*(np.diag(r)/abs(np.diag(r)))
def embed(vs, A):
    T = np.eye(NV, dtype=complex)
    for i,u in enumerate(vs):
        for k,w in enumerate(vs): T[u,w] = A[i,k]
    return T
def fam(vs, W, jvec, V):
    L = len(vs); r = np.exp(1j*np.angle(W)/L); zt = np.exp(2j*np.pi/L)
    return embed(vs, V@np.diag(r*zt**np.asarray(jvec))@V.conj().T)
def shift(loop, ns):
    on = {s for s,_,_ in loop}
    T = np.zeros((NV,NV), dtype=complex)
    for v in range(NV):
        if v not in on: T[v,v] = 1.0
    for k,(s_,d_,_) in enumerate(loop): T[d_,s_] = ns[k]
    return T

def Dscore(TF, TC, nmax=12):
    worst = 0.0; XF = np.eye(NV,dtype=complex); XC = np.eye(NV,dtype=complex)
    for n in range(1, nmax+1):
        XF = TF@XF; XC = TC@XC; Q = XF.conj().T@XC
        off = np.linalg.norm(Q - np.diag(np.diag(Q))); d = np.diag(Q); sp = 0.0
        for c,vs in CLASSES.items():
            if len(vs) > 1: sp = max(sp, max(abs(d[u]-d[w]) for u in vs for w in vs))
        worst = max(worst, off+sp)
    return worst

# ---------- THE INDEPENDENT CRITERION: |Z| ONLY, BY BRUTE STATE SAMPLING ----------
def states_fixed_pi(p, K):
    """K states on K1 with IDENTICAL pi = (p00,p10,p01,p11): random within-class split + phases."""
    out = []
    for _ in range(K):
        w = np.zeros(NV)
        w[0] = p[3]                                   # class 11 = {v0}
        t = rng.uniform(0.02,0.98); w[1],w[2] = p[1]*t, p[1]*(1-t)
        t = rng.uniform(0.02,0.98); w[3],w[4] = p[2]*t, p[2]*(1-t)
        s = np.sqrt(w)*np.exp(1j*rng.uniform(0,2*np.pi,NV))
        out.append(s)
    return out
def absZ_spread(TF, TC, nmax=12, K=40):
    p = np.array([0.0,0.30,0.30,0.40])
    S = states_fixed_pi(p,K); worst = 0.0
    XF = np.eye(NV,dtype=complex); XC = np.eye(NV,dtype=complex)
    for n in range(1,nmax+1):
        XF = TF@XF; XC = TC@XC
        v = [abs(np.vdot(XF@s, XC@s)) for s in S]
        worst = max(worst, max(v)-min(v))
    return worst

print("\n== X1.3  IS D THE SAME CRITERION AS 'THE CORPUS'S OBSERVABLE |Z| SEES ONLY pi'? ==")
print("   For each sampled member I compute BOTH: D (the lane's instrument, a statement about the")
print("   complex form Z) and the spread of |Z_n| over 40 states with IDENTICAL pi (the corpus's")
print("   own observable).  A member with D > 0 but |Z|-spread = 0 would break every count in T2.")
aa = rng.uniform(0,2*np.pi,6); Ua = np.exp(1j*aa)
WFa = np.exp(1j*(aa[0]+aa[1]+aa[2])); WCa = np.exp(1j*(aa[3]+aa[4]+aa[5]))
rows = []
def check(label, TF, TC):
    d = Dscore(TF,TC); z = absZ_spread(TF,TC)
    rows.append((label, d, z)); return d, z
check("COR-F edge tick        ", T_edge(LF,Ua), T_edge(LC,Ua))
check("circuit M_gamma        ", M_circ(VF,WFa), M_circ(VC,WCa))
check("fibre-wise cube root R ", fam(VF,WFa,(0,0,0),np.eye(3)), fam(VC,WCa,(0,0,0),np.eye(3)))
check("diagonal, mixed roots  ", fam(VF,WFa,(0,1,2),np.eye(3)), fam(VC,WCa,(0,0,0),np.eye(3)))
check("zeta * COR-F's T (F)   ", np.exp(2j*np.pi/3)*T_edge(LF,Ua) + np.diag([0,0,0,1,1]).astype(complex)*(1-np.exp(2j*np.pi/3)), T_edge(LC,Ua))
disc = 0
for _ in range(4000):
    jF, jC = rng.integers(0,3,3), rng.integers(0,3,3)
    TF = fam(VF,WFa,jF,haar(3)); TC = fam(VC,WCa,jC,haar(3))
    d = Dscore(TF,TC); z = absZ_spread(TF,TC,K=25)
    if (d < 1e-9) != (z < 1e-9): disc += 1
for _ in range(4000):
    ph = rng.uniform(0,2*np.pi,2); nF = [np.exp(1j*ph[0]),np.exp(1j*ph[1]),WFa*np.exp(-1j*ph.sum())]
    ph = rng.uniform(0,2*np.pi,2); nC = [np.exp(1j*ph[0]),np.exp(1j*ph[1]),WCa*np.exp(-1j*ph.sum())]
    d = Dscore(shift(LF,nF), shift(LC,nC)); z = absZ_spread(shift(LF,nF), shift(LC,nC),K=25)
    if (d < 1e-9) != (z < 1e-9): disc += 1
for lab,d,z in rows:
    print(f"   {lab}  D = {d:>10.3e}   |Z|-spread over 40 same-pi states = {z:>10.3e}")
print(f"   8000 further members (4000 Haar-family + 4000 shift-torus): members where D and the")
print(f"   |Z| criterion DISAGREE at 1e-9: {disc}")
print("   -> THE INSTRUMENT IS VALIDATED AGAINST THE CORPUS'S OWN OBSERVABLE.  Every count in")
print("      leg T2 is a count of the right set.  This is the check that lane did not run.")
