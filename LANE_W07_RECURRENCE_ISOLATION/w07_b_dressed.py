# W-07 leg B — rebuild W-06's dressed gauge-invariant observable and re-run ITS recurrence test.
# ISOLATION LEDGER: same carrier, same ready state, same observable, same k-range, same code path.
#   THE ONE VARIABLE THAT MOVES IS THE ARITHMETIC OF THE CONNECTION. Nothing else is touched.
import numpy as np
np.set_printoptions(precision=6, suppress=True)
rng = np.random.default_rng(20260816)

EDGES=[(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)]
FACE_V={0,1,2}; CYC_V={0,3,4}
TREE={1:(0,),2:(0,1),3:(3,),4:(3,4)}          # tree paths v0->v as edge-index lists (e1,e2,e4,e5)

def U(a): return np.exp(1j*np.asarray(a))
def dress(s,a):                                # t_v = W(tree path v0->v)^-1 s_v
    u=U(a); t=np.array(s,dtype=complex)
    for v,path in TREE.items():
        w=1.0+0j
        for e in path: w*=u[e]
        t[v]=s[v]/w
    return t
def gauge(s,a,th):                             # FULL action: S1:63 on the connection AND on the section
    a2=np.array(a,dtype=float)
    for j,(p,q) in enumerate(EDGES): a2[j]=a[j]+th[q]-th[p]
    return np.exp(1j*np.asarray(th))*s, a2
def A(s,a,u,v): t=dress(s,a); return np.conj(t[u])*t[v]      # the dressed observable

print("== B1  GAUGE INVARIANCE OF THE DRESSED OBSERVABLE, UNDER THE FULL ACTION ==")
a0=np.array([np.pi/3]*3+[np.pi/2]*3); s0=rng.normal(size=5)+1j*rng.normal(size=5); s0/=np.linalg.norm(s0)
worst=0.0
for _ in range(2000):
    th=rng.uniform(0,2*np.pi,5); sg,ag=gauge(s0,a0,th)
    for (u,v) in [(2,3),(1,4),(0,3),(2,4)]:
        worst=max(worst,abs(A(sg,ag,u,v)-A(s0,a0,u,v)))
print(f"  max |A(gauge.s, gauge.a) - A(s,a)| over 2000 gauge transforms, 4 observables = {worst:.3e}")
print("  -> gauge-invariant. W-06's structural claim REPRODUCED, independently, in this lineage.\n")

def sep_profile(a, s, u, v, K=4000):
    """|A[M_dF^k s] - A[M_c^k s]| for k=1..K.  M_dF multiplies FACE_V by W_F; M_c multiplies CYC_V by W_C."""
    WF=np.exp(1j*(a[0]+a[1]+a[2])); WC=np.exp(1j*(a[3]+a[4]+a[5]))
    dF=(v in FACE_V)-(u in FACE_V); dC=(v in CYC_V)-(u in CYC_V)
    amp=abs(A(s,a,u,v)); k=np.arange(1,K+1)
    return amp*np.abs(WF**(dF*k)-WC**(dC*k)), WF, WC, amp

def report(tag, a, s, u=2, v=3, K=4000):
    D,WF,WC,amp=sep_profile(a,s,u,v,K)
    ratio=np.conj(WF)/WC
    print(f"  {tag}")
    print(f"    W_F={WF:+.6f}  W_C={WC:+.6f}   conj(W_F)/W_C = {ratio:+.6f}  arg/2pi = {np.angle(ratio)/(2*np.pi):.9f}")
    print(f"    separation at k=1 : {D[0]:.12f}      max over k<=4000 : {D.max():.12f}")
    print(f"    min over k<=4000  : {D.min():.3e}     cells with D < 1e-9 : {int((D<1e-9).sum())} of {K}")
    return D

print("== B2  THE RECURRENCE TEST, RE-RUN.  ONE VARIABLE: THE CONNECTION'S ARITHMETIC ==")
s=rng.normal(size=5)+1j*rng.normal(size=5); s/=np.linalg.norm(s)   # generic ready state, held FIXED below
print("  ready state held fixed:  |s| =", np.round(np.abs(s),6), "\n")

a_pub=np.array([np.pi/3]*3+[np.pi/2]*3)                            # S1 sec6, PUBLISHED
report("(i)  S1's PUBLISHED connection  a=(pi/3 x3, pi/2 x3)", a_pub, s)
print()
phi=(1+5**0.5)/2                                                    # badly approximable: worst case for near-returns
a_gen=np.array([2*np.pi*phi/3]*3+[2*np.pi*(phi**2)/3]*3)
report("(ii) a GENERIC connection, same carrier, same state, same observable", a_gen, s)
print()
a_rnd=rng.uniform(0,2*np.pi,6)
report("(iii) a RANDOM connection (seed 20260816)", a_rnd, s)
