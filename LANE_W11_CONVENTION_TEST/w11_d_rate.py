# W-11 leg D — N1 is the result proposed for publication. Does the RATE survive the convention?
import numpy as np
rng=np.random.default_rng(20260817)
LOOP_F=[(0,1,0),(1,2,1),(2,0,2)]; LOOP_C=[(0,3,3),(3,4,4),(4,0,5)]
FACE_V={0,1,2}; CYC_V={0,3,4}
def Top(loop,a):
    U=np.exp(1j*np.asarray(a)); T=np.zeros((5,5),dtype=complex); on={v for v,_,_ in loop}
    for v in range(5):
        if v not in on: T[v,v]=1.0
    for (s_,d_,e) in loop: T[d_,s_]=U[e]
    return T
def Mop(vs,W):
    M=np.eye(5,dtype=complex)
    for v in vs: M[v,v]=W
    return M
def hol(a): return np.exp(1j*(a[0]+a[1]+a[2])), np.exp(1j*(a[3]+a[4]+a[5]))
def m_jensen(p,n=1<<20):
    a,b,c,d=p; t=2*np.pi*np.arange(n)/n; ct=np.cos(t)
    A=np.sqrt(np.maximum(a*a+b*b+2*a*b*ct,0)); B=np.sqrt(np.maximum(c*c+d*d+2*c*d*ct,0))
    return np.log(np.maximum(A,B)+1e-300).mean()

a=np.array([1.0,0.37,0.91,2**0.5,0.23,1.77])
WF,WC=hol(a); TF,TC=Top(LOOP_F,a),Top(LOOP_C,a); MF,MC=Mop(FACE_V,WF),Mop(CYC_V,WC)
sA=np.sqrt(np.array([0.40,0.15,0.15,0.15,0.15]))+0j
sB=np.sqrt(np.array([0.40,0.30,0.00,0.05,0.25]))+0j
sC=sA*np.exp(1j*np.array([0.0,1.3,-0.7,2.2,0.4]))
PI=np.array([0.0,0.30,0.30,0.40])       # p00,p10,p01,p11 -- identical for A, B, C

def rate(s,op_F,op_C,N):
    """(1/N) sum log|Z| accumulated by repeated application -- no matrix powers, no overflow"""
    xF=s.copy(); xC=s.copy(); tot=0.0
    for _ in range(N):
        xF=op_F@xF; xC=op_C@xC
        z=abs(np.vdot(xF,xC))
        tot+=np.log(z) if z>0 else -700.0
    return tot/N

print("== D1  THE RATE UNDER EACH CONVENTION, SAME pi, SAME CONNECTION ==")
print(f"  N1's registered value      lambda = m(P) = {m_jensen(PI):.12f}   for pi = {PI}")
print()
print(f"  {'':<28}{'state A':>16}{'state B':>16}{'state C':>16}{'spread':>11}")
for N in (2000, 20000, 200000):
    rc=[rate(s,MF,MC,N) for s in (sA,sB,sC)]
    re_=[rate(s,TF,TC,N) for s in (sA,sB,sC)]
    print(f"  CIRCUIT  N={N:<7} per circuit {rc[0]:>15.9f}{rc[1]:>16.9f}{rc[2]:>16.9f}{max(rc)-min(rc):>11.1e}")
    print(f"  EDGE     N={N:<7} per tick    {re_[0]:>15.9f}{re_[1]:>16.9f}{re_[2]:>16.9f}{max(re_)-min(re_):>11.1e}")
    print(f"  EDGE, rescaled x3 (per circuit, |gamma|=3)  {3*re_[0]:>13.9f}{3*re_[1]:>16.9f}{3*re_[2]:>16.9f}")
    print()
print("  -> the CIRCUIT rate is m(P) and is blind to which state it is, as registered.")
print("     the EDGE rate is NOT m(P)/3 and is NOT the same number for the three states.")
