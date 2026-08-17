# W-11 leg A — build COR-F's edge-tick transport on K1 and validate it against the sealed text
# BEFORE testing anything with it.  COR-F: S3_THE_CROSSING_AUDIT_V001.md:160-209, :794.
import numpy as np
rng = np.random.default_rng(20260817)
EDGES=[(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)]            # e1..e6, S1 sec1
LOOP_F=[(0,1,0),(1,2,1),(2,0,2)]                       # (from, to, edge index) around gamma_F
LOOP_C=[(0,3,3),(3,4,4),(4,0,5)]                       # around gamma_C
FACE_V={0,1,2}; CYC_V={0,3,4}

def Top(loop, a):
    """COR-F's edge tick: move each fibre value one edge along the loop, identity off it."""
    U=np.exp(1j*np.asarray(a)); T=np.zeros((5,5),dtype=complex)
    on={v for v,_,_ in loop}
    for v in range(5):
        if v not in on: T[v,v]=1.0
    for (src,dst,e) in loop: T[dst,src]=U[e]            # (Ts)(dst) = U_e s(src)
    return T
def Mop(vs, W):
    M=np.eye(5,dtype=complex)
    for v in vs: M[v,v]=W
    return M
def hol(a): return np.exp(1j*(a[0]+a[1]+a[2])), np.exp(1j*(a[3]+a[4]+a[5]))

print("== A1  VALIDATE AGAINST COR-F's SEALED EXHIBIT ==")
a=np.zeros(6); a[3],a[4],a[5]=0.7,1.3,-0.4             # COR-F's own a4,a5,a6
T=Top(LOOP_C,a); WF,WC=hol(a)
print(f"  || T*T - I ||            = {np.linalg.norm(T.conj().T@T-np.eye(5)):.2e}      [COR-F: 0.00e+00]")
print(f"  T diagonal?                {np.allclose(T,np.diag(np.diag(T)))}                    [COR-F: False]")
print(f"  T^3 = diag(...)          = {np.round(np.diag(np.linalg.matrix_power(T,3)),6)}")
print(f"  W_C                      = {WC:.6f}      [COR-F: -0.029200+0.999574j]")
rho=np.diag([0.40,0.15,0.15,0.15,0.15]).astype(complex)
print(f"  diag(T rho T*)           = {np.round(np.real(np.diag(T@rho@T.conj().T)),2)}   [COR-F: 0.15 0.15 0.15 0.40 0.15]")
print(f"  diag(rho)                = {np.round(np.real(np.diag(rho)),2)}   NOT PRESERVED\n")

print("== A2  T^L = M_gamma, AND GAUGE COVARIANCE ==")
worstF=worstC=worstG=0.0
for _ in range(2000):
    a=rng.uniform(0,2*np.pi,6); WF,WC=hol(a)
    TF,TC=Top(LOOP_F,a),Top(LOOP_C,a)
    worstF=max(worstF,np.linalg.norm(np.linalg.matrix_power(TF,3)-Mop(FACE_V,WF)))
    worstC=max(worstC,np.linalg.norm(np.linalg.matrix_power(TC,3)-Mop(CYC_V,WC)))
    # gauge: a_e -> a_e + th[t] - th[s];  s_v -> e^{i th_v} s_v.  T_F s must transform as a section.
    th=rng.uniform(0,2*np.pi,5); ag=np.array([a[j]+th[t]-th[s] for j,(s,t) in enumerate(EDGES)])
    s=rng.normal(size=5)+1j*rng.normal(size=5)
    lhs=Top(LOOP_F,ag)@(np.exp(1j*th)*s); rhs=np.exp(1j*th)*(Top(LOOP_F,a)@s)
    worstG=max(worstG,np.linalg.norm(lhs-rhs))
print(f"  max || T_F^3 - M_dF ||                      = {worstF:.2e}   (2000 connections)")
print(f"  max || T_C^3 - M_c  ||                      = {worstC:.2e}")
print(f"  max || T_F(g.s, g.a) - g.(T_F s) ||         = {worstG:.2e}   -> T is gauge-COVARIANT,")
print("     so <T_F^n s, T_C^n s> is gauge-INVARIANT, exactly as <M_dF^k s, M_c^k s> is.")
print("  ==> COR-F's transport is a legitimate rival convention on the corpus's own terms:")
print("      unitary, gauge-covariant, and its L-th power IS the corpus's own operator.")
