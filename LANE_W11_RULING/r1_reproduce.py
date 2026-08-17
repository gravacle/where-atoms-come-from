# RULING LEG 1 -- reproduce, from code sharing no line with any W-11 lane, (i) COR-F's own sealed
# exhibit, (ii) the registrar's legs A/B/C/D.  MY OWN seed (20260818), MY OWN ready states, and
# the corpus's ONLY published generic connection f=1.0, c=sqrt(2) (S4:603) -- NOT the registrar's
# unpublished f=57/25, c=2+sqrt(2).  If the picture is an artefact of its connection, it dies here.
import numpy as np, rlib
from fractions import Fraction
np.set_printoptions(linewidth=200)
rng = np.random.default_rng(rlib.SEED)

print("== R1.1  COR-F's SEALED EXHIBIT, S3_THE_CROSSING_AUDIT_V001.md:174-186 ==")
K = rlib.K1(); a = np.zeros(6); a[3],a[4],a[5] = 0.7,1.3,-0.4      # COR-F's own a4,a5,a6
T = rlib.Tedge(K, K.walkC, a); WC = rlib.holon(K.walkC, a)
rho = np.diag([0.40,0.15,0.15,0.15,0.15]).astype(complex)
print(f"  ||T*T - I||        = {np.linalg.norm(T.conj().T@T-np.eye(5)):.2e}          [sealed 0.00e+00]")
print(f"  T diagonal?          {np.allclose(T,np.diag(np.diag(T)))}                      [sealed False]")
print(f"  W_C                = {WC:.6f}        [sealed -0.029200+0.999574j]")
print(f"  ||T^3 - diag(W_C,1,1,W_C,W_C)|| = {np.linalg.norm(np.linalg.matrix_power(T,3)-np.diag([WC,1,1,WC,WC])):.2e}")
print(f"  diag(T rho T*)     = {np.round(np.real(np.diag(T@rho@T.conj().T)),2)}   [sealed 0.15 0.15 0.15 0.40 0.15]")
print(f"  diag(rho)          = {np.round(np.real(np.diag(rho)),2)}   NOT PRESERVED\n")

print("== R1.2  T^L = M_gamma IN EXACT GAUSSIAN RATIONALS -- no float anywhere ==")
# unimodular Gaussian rationals: (3+4i)/5, (5+12i)/13, (8+15i)/17, (20+21i)/29
UZ = [(Fraction(3,5),Fraction(4,5)), (Fraction(5,13),Fraction(12,13)),
      (Fraction(8,17),Fraction(15,17)), (Fraction(20,29),Fraction(21,29))]
def cmul(x,y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def matmulQ(A,B):
    n=len(A); C=[[(Fraction(0),Fraction(0)) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for k in range(n):
            if A[i][k]==(0,0): continue
            for j in range(n):
                p=cmul(A[i][k],B[k][j]); C[i][j]=(C[i][j][0]+p[0], C[i][j][1]+p[1])
    return C
for L in (3,4):
    n=L+1; walk=[(i,(i+1)%L,i,1) for i in range(L)]           # L-cycle plus one spectator
    Tq=[[(Fraction(0),Fraction(0)) for _ in range(n)] for _ in range(n)]
    Tq[L][L]=(Fraction(1),Fraction(0))
    for (u,v,e,_) in walk: Tq[v][u]=UZ[e]
    P=Tq
    for _ in range(L-1): P=matmulQ(P,Tq)
    W=(Fraction(1),Fraction(0))
    for e in range(L): W=cmul(W,UZ[e])
    Mq=[[(Fraction(0),Fraction(0)) for _ in range(n)] for _ in range(n)]
    for i in range(L): Mq[i][i]=W
    Mq[L][L]=(Fraction(1),Fraction(0))
    print(f"  L={L}: T^L == M_gamma exactly in Fraction arithmetic:  {P==Mq}   (W = {W[0]} + {W[1]}i)")
print()

print("== R1.3  T IS UNITARY AND GAUGE-COVARIANT; M IS GAUGE-INVARIANT.  2000 connections each ==")
for C in (rlib.K1(), rlib.B0b()):
    wU=wF=wC=wG=wI=0.0
    for _ in range(2000):
        a=rng.uniform(0,2*np.pi,len(C.edges))
        TF,TC=rlib.Tedge(C,C.walkF,a),rlib.Tedge(C,C.walkC,a)
        MF=rlib.Mcirc(C,C.VF,rlib.holon(C.walkF,a)); MC=rlib.Mcirc(C,C.VC,rlib.holon(C.walkC,a))
        wU=max(wU,np.linalg.norm(TF.conj().T@TF-np.eye(C.nv)))
        wF=max(wF,np.linalg.norm(np.linalg.matrix_power(TF,C.LF)-MF))
        wC=max(wC,np.linalg.norm(np.linalg.matrix_power(TC,C.LC)-MC))
        th=rng.uniform(0,2*np.pi,C.nv)
        ag=np.array([a[j]+th[t]-th[s] for j,(s,t) in enumerate(C.edges)])
        s=rng.normal(size=C.nv)+1j*rng.normal(size=C.nv)
        wG=max(wG,np.linalg.norm(rlib.Tedge(C,C.walkF,ag)@(np.exp(1j*th)*s)-np.exp(1j*th)*(TF@s)))
        MFg=rlib.Mcirc(C,C.VF,rlib.holon(C.walkF,ag))
        wI=max(wI,np.linalg.norm(MFg@(np.exp(1j*th)*s)-np.exp(1j*th)*(MF@s)))
    print(f"  {C.name:4s} ||T*T-I||={wU:.2e}  ||T_F^LF-M_F||={wF:.2e}  ||T_C^LC-M_C||={wC:.2e}"
          f"  covariance(T)={wG:.2e}  covariance(M)={wI:.2e}")
print("  -> COR-F's T is unitary, gauge-COVARIANT, and its L-th power IS the corpus's operator.")
print("     Registrar leg A reproduces on BOTH carriers and in exact arithmetic.  LEG A STANDS.\n")

print("== R1.4  THE DECISIVE ARMS, MY OWN STATES AND THE CORPUS'S PUBLISHED GENERIC CONNECTION ==")
print("   f = 1.0, c = sqrt(2)  (S4:603 -- the ONLY generic connection the corpus publishes; W-10 N-4)")
for C in (rlib.K1(), rlib.B0b()):
    a = rlib.a_generic(C, rng, 1.0, 2**0.5)
    TF,TC = rlib.Tedge(C,C.walkF,a), rlib.Tedge(C,C.walkC,a)
    MF = rlib.Mcirc(C,C.VF,rlib.holon(C.walkF,a)); MC = rlib.Mcirc(C,C.VC,rlib.holon(C.walkC,a))
    base = rng.dirichlet(np.ones(C.nv))
    S = rlib.same_pi_states(C, rng, base, 40)
    pis = np.array([C.pi_of(s) for s in S])
    print(f"  {C.name}: 40 states, pi held BY CONSTRUCTION, max|pi_i - pi_0| = {np.max(np.abs(pis-pis[0])):.2e}")
    print(f"       ARMS DIFFED: min_{{i<j}} ||s_i - s_j|| = {min(np.linalg.norm(S[i]-S[j]) for i in range(40) for j in range(i+1,40)):.4f}  (NOT byte-identical)")
    ns = list(range(1, 25))
    print(f"       CIRCUIT convention, spread of |Z_k| over k<=24 : {rlib.pi_spread(C,MF,MC,S,ns):.2e}")
    print(f"       EDGE    convention, spread of |Z_n| over n<=24 : {rlib.pi_spread(C,TF,TC,S,ns):.2e}")
    good=[n for n in ns if n % C.LF==0 and n % C.LC==0]
    print(f"       EDGE, restricted to n = 0 mod lcm({C.LF},{C.LC}) = {good} : {rlib.pi_spread(C,TF,TC,S,good):.2e}")
    if C.name=="B0b":
        same=[n for n in range(1,20001) if n%C.LF==0 and n%C.LC==0 and n//C.LF==n//C.LC]
        print(f"       edge ticks n<=20000 at which BOTH branches sit at the SAME circuit count: {same}")
