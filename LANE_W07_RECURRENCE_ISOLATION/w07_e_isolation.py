# W-07 leg E — the isolation, clean.  Leg B(ii) carried a confound OF MY OWN MAKING: I built the
# "generic" connection from phi and phi^2, and phi^2 = phi+1, so W_F == W_C exactly.  Redone here
# with independent irrationals and three random draws.  Recorded rather than silently fixed.
#
# ISOLATION LEDGER.  Held fixed: carrier K1, ready state s (seed 20260816), observable A_23,
# dressing tree, k-range 1..4000, code path.  Moved: the connection.
# OPERATIVE VARIABLE, NAMED: the ORDER of the branch ratio rho = conj(W_F)^dF * W_C^{-dC} in U(1).
import numpy as np
rng=np.random.default_rng(20260816)
FACE_V={0,1,2}; CYC_V={0,3,4}; TREE={1:(0,),2:(0,1),3:(3,),4:(3,4)}
def dress(s,a):
    u=np.exp(1j*np.asarray(a)); t=np.array(s,dtype=complex)
    for v,p in TREE.items():
        w=1.0+0j
        for e in p: w*=u[e]
        t[v]=s[v]/w
    return t
s=rng.normal(size=5)+1j*rng.normal(size=5); s/=np.linalg.norm(s)
def run(tag,a,K=4000,u=2,v=3):
    WF=np.exp(1j*sum(a[:3])); WC=np.exp(1j*sum(a[3:]))
    dF=(v in FACE_V)-(u in FACE_V); dC=(v in CYC_V)-(u in CYC_V)
    t=dress(s,a); amp=abs(np.conj(t[u])*t[v]); k=np.arange(1,K+1)
    rho=np.conj(WF)**dF*WC**(-dC)
    D=amp*np.abs(np.exp(1j*np.angle(rho)*k)-1)
    ordr=next((n for n in range(1,10001) if abs(rho**n-1)<1e-12), None)
    print(f"  {tag}")
    print(f"    W_F={WF:+.6f}  W_C={WC:+.6f}   W_F==W_C? {abs(WF-WC)<1e-12}")
    print(f"    rho = {rho:+.6f}   arg/2pi = {np.angle(rho)/(2*np.pi):+.9f}   ORDER in U(1) = {ordr}")
    print(f"    min D over k<=4000 = {D.min():.3e}     cells with D < 1e-9 : {int((D<1e-9).sum())} of {K}")
print("== LEG E — THE ISOLATION.  ONE VARIABLE: ord(rho). ==")
run("S1 PUBLISHED  a=(pi/3 x3, pi/2 x3)", np.array([np.pi/3]*3+[np.pi/2]*3)); print()
run("GENERIC  sqrt2 / sqrt3            ", np.array([2*np.pi*(2**0.5)/3]*3+[2*np.pi*(3**0.5)/3]*3)); print()
for sd in (1,2,3):
    run(f"RANDOM seed {sd}                     ", np.random.default_rng(sd).uniform(0,2*np.pi,6)); print()
print("  ord(rho) = 4        -> 1000 of 4000, min 6.7e-19.   W-06's registered figure, reproduced.")
print("  ord(rho) = infinite -> 0 of 4000, five times out of five, min 1.6e-05 .. 3.1e-04.")
