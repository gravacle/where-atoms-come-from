# W-11 leg C — the same test on B0b, where the two loops have DIFFERENT LENGTHS (4 and 3).
# B0b's incidence is taken verbatim from LANE_W10_A_CARRIERS/w10a_lib.py:182-213 (sealed, and
# reproduced there against S4:512/:539/:575 by a build lane and two refuters).
import numpy as np
rng=np.random.default_rng(20260817)
def V(i,j): return 3*(j%3)+(i%3)
def H(i,j): return 3*(j%3)+(i%3)
def Wg(i,j): return 9+3*(j%3)+(i%3)
E=[None]*18
for j in range(3):
    for i in range(3):
        E[H(i,j)]=(V(i,j),V(i+1,j)); E[Wg(i,j)]=(V(i,j),V(i,j+1))
gF=[(H(0,0),1),(Wg(1,0),1),(H(0,1),-1),(Wg(0,0),-1)]     # face (0,0) boundary, length 4
gC=[(H(0,0),1),(H(1,0),1),(H(2,0),1)]                    # horizontal row j=0, length 3
NV=9
def walk(g):
    """directed vertex walk and per-step transport, following the signed edge cycle"""
    steps=[]; 
    for (e,s) in g:
        u,v=E[e]
        steps.append((u,v,e,+1) if s>0 else (v,u,e,-1))
    return steps
def loop_vs(g): return {u for u,_,_,_ in walk(g)}
WF_V, WC_V = loop_vs(gF), loop_vs(gC)
print(f"  gamma_F length {len(gF)}, vertices {sorted(WF_V)}")
print(f"  gamma_C length {len(gC)}, vertices {sorted(WC_V)}")
cls=lambda v:(int(v in WF_V),int(v in WC_V))
from collections import Counter
print(f"  class multiset {dict(Counter(''.join(map(str,cls(v))) for v in range(NV)))}   [S4:575 {{00:4,01:1,10:2,11:2}}]\n")

def Top(g,a):
    U=np.exp(1j*np.asarray(a)); T=np.zeros((NV,NV),dtype=complex); on=loop_vs(g)
    for v in range(NV):
        if v not in on: T[v,v]=1.0
    for (u,v,e,sgn) in walk(g): T[v,u]= U[e] if sgn>0 else np.conj(U[e])
    return T
def Mop(vs,W):
    M=np.eye(NV,dtype=complex)
    for v in vs: M[v,v]=W
    return M
def hol(g,a):
    z=1.0+0j
    for (u,v,e,sgn) in walk(g): z*= np.exp(1j*a[e]) if sgn>0 else np.exp(-1j*a[e])
    return z

a=rng.uniform(0,2*np.pi,18)
TF,TC=Top(gF,a),Top(gC,a); MF,MC=Mop(WF_V,hol(gF,a)),Mop(WC_V,hol(gC,a))
print("== C1  T^L = M_gamma ON B0b, WITH THE TWO LENGTHS DIFFERENT ==")
print(f"  || T_F^4 - M_dF || = {np.linalg.norm(np.linalg.matrix_power(TF,4)-MF):.2e}   (|gamma_F| = 4)")
print(f"  || T_C^3 - M_c  || = {np.linalg.norm(np.linalg.matrix_power(TC,3)-MC):.2e}   (|gamma_C| = 3)\n")

print("== C2  IS THE CIRCUIT CONVENTION A SUBSEQUENCE OF THE EDGE CONVENTION HERE? ==")
print("  It is on K1 only because both loops have length 3. Here lcm(4,3) = 12, and at n = 12")
print("  T_F^12 = M_dF^3 while T_C^12 = M_c^4 -- the two BRANCHES are at DIFFERENT circuit counts.")
print(f"  || T_F^12 - M_dF^3 || = {np.linalg.norm(np.linalg.matrix_power(TF,12)-np.linalg.matrix_power(MF,3)):.2e}")
print(f"  || T_C^12 - M_c^4  || = {np.linalg.norm(np.linalg.matrix_power(TC,12)-np.linalg.matrix_power(MC,4)):.2e}")
same=[n for n in range(1,2001) if n%4==0 and n%3==0 and (n//4)==(n//3)]
print(f"  edge ticks n <= 2000 at which BOTH branches sit at the SAME circuit count: {same}")
print("  -> NONE. On B0b the circuit convention is not a subsequence of the edge convention.\n")

print("== C3  THE DECISIVE TEST ON B0b.  SAME pi, ONE VARIABLE: THE CONVENTION ==")
def pi_of(s):
    w=np.abs(s)**2; p={ '00':0.,'10':0.,'01':0.,'11':0.}
    for v in range(NV): p[f"{cls(v)[0]}{cls(v)[1]}"]+=w[v]
    return np.array([p['00'],p['10'],p['01'],p['11']])
wA=np.array([.10,.12,.09,.14,.11,.11,.11,.11,.11]); wA/=wA.sum()
# same class sums, weight shuffled WITHIN each class
wB=wA.copy()
wB[0],wB[1]=wA[0]+wA[1],0.0                 # class 11 = {0,1}
wB[3],wB[4]=0.0,wA[3]+wA[4]                 # class 10 = {3,4}
wB[5],wB[8]=wA[5]+wA[8],0.0                 # class 00 = {5,6,7,8}
sA,sB=np.sqrt(wA)+0j, np.sqrt(wB)+0j
sC=sA*np.exp(1j*rng.uniform(0,2*np.pi,NV))  # same moduli, different phases
print(f"  pi(A) = {np.round(pi_of(sA),12)}")
print(f"  pi(B) = {np.round(pi_of(sB),12)}")
print(f"  pi(C) = {np.round(pi_of(sC),12)}")
assert np.allclose(pi_of(sA),pi_of(sB)) and np.allclose(pi_of(sA),pi_of(sC))
print()
def Zc(s,k): return np.vdot(np.linalg.matrix_power(MF,k)@s, np.linalg.matrix_power(MC,k)@s)
def Ze(s,n): return np.vdot(np.linalg.matrix_power(TF,n)@s, np.linalg.matrix_power(TC,n)@s)
print(f"  CIRCUIT convention (M):  {'k':>2} {'|Z(A)|':>15} {'|Z(B)|':>15} {'|Z(C)|':>15} {'spread':>10}")
for k in range(1,6):
    v=[abs(Zc(s,k)) for s in (sA,sB,sC)]
    print(f"  {'':<25}{k:>2} {v[0]:>15.12f} {v[1]:>15.12f} {v[2]:>15.12f} {max(v)-min(v):>10.1e}")
print(f"  EDGE convention (T):     {'n':>2} {'|Z(A)|':>15} {'|Z(B)|':>15} {'|Z(C)|':>15} {'spread':>10}")
for n in range(1,13):
    v=[abs(Ze(s,n)) for s in (sA,sB,sC)]
    print(f"  {'':<25}{n:>2} {v[0]:>15.12f} {v[1]:>15.12f} {v[2]:>15.12f} {max(v)-min(v):>10.1e}")
