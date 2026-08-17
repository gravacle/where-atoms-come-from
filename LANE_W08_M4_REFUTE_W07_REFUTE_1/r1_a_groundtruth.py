# W-08 / M4-REFUTE-1  leg A — GROUND TRUTH BY BRUTE FORCE, NO ALGEBRA.
# LENS: STEELMAN W-07.  Before crediting M4's central algebraic finding (M4-2, "W-07 leg E computes
# the wrong group element"), verify the governing element by literally applying the branch operators
# and measuring, with no closed form used anywhere.
#
# ISOLATION LEDGER.  Held fixed: carrier K1, dressing tree {e1,e2,e4,e5} rooted v0, ready state
# (numpy default_rng(20260816), W-07 leg E's own), observable A_23, k-range 1..4000.
# MOVED: nothing.  This leg draws no comparison between conditions; it compares three FORMULAS
# (brute force / W-07 leg B / W-07 leg E) against the same single computation.
# PRECISION: numpy float64 (double).  Stated.  The one exact-zero claim is deferred to leg D.
import numpy as np
np.set_printoptions(precision=6, suppress=True)

FACE_V={0,1,2}; CYC_V={0,3,4}; TREE={1:(0,),2:(0,1),3:(3,),4:(3,4)}
def dress(s,a):
    u=np.exp(1j*np.asarray(a)); t=np.array(s,dtype=complex)
    for v,p in TREE.items():
        w=1.0+0j
        for e in p: w*=u[e]
        t[v]=s[v]/w
    return t
def A(s,a,u,v):
    t=dress(s,a); return np.conj(t[u])*t[v]

rng=np.random.default_rng(20260816)
s=rng.normal(size=5)+1j*rng.normal(size=5); s/=np.linalg.norm(s)
K=4000; u,v=2,3
dF=(v in FACE_V)-(u in FACE_V); dC=(v in CYC_V)-(u in CYC_V)
print(f"== A0  brute force: apply M_dF^k and M_c^k literally, dress, subtract.  (u,v)=({u},{v}), dF={dF}, dC={dC}")

def brute(a,K=K):
    WF=np.exp(1j*sum(a[:3])); WC=np.exp(1j*sum(a[3:]))
    out=np.empty(K)
    x=np.array(s,dtype=complex); y=np.array(s,dtype=complex)
    for kk in range(1,K+1):
        for w_ in range(5):
            if w_ in FACE_V: x[w_]*=WF
            if w_ in CYC_V:  y[w_]*=WC
        out[kk-1]=abs(A(x,a,u,v)-A(y,a,u,v))
    return out

def legB(a,K=K):                      # w07_b_dressed.sep_profile, verbatim algebra
    WF=np.exp(1j*sum(a[:3])); WC=np.exp(1j*sum(a[3:]))
    amp=abs(A(s,a,u,v)); k=np.arange(1,K+1)
    return amp*np.abs(WF**(dF*k)-WC**(dC*k))

def legE(a,K=K):                      # w07_e_isolation.run, verbatim
    WF=np.exp(1j*sum(a[:3])); WC=np.exp(1j*sum(a[3:]))
    amp=abs(A(s,a,u,v)); k=np.arange(1,K+1)
    rho=np.conj(WF)**dF*WC**(-dC)
    return amp*np.abs(np.exp(1j*np.angle(rho)*k)-1)

def rho_true(a):
    WF=np.exp(1j*sum(a[:3])); WC=np.exp(1j*sum(a[3:])); return WF**dF*WC**(-dC)
def rho_E(a):
    WF=np.exp(1j*sum(a[:3])); WC=np.exp(1j*sum(a[3:])); return np.conj(WF)**dF*WC**(-dC)

cases=[("S1 PUBLISHED (pi/3 x3, pi/2 x3)", np.array([np.pi/3]*3+[np.pi/2]*3)),
       ("GENERIC sqrt2/sqrt3            ", np.array([2*np.pi*(2**0.5)/3]*3+[2*np.pi*(3**0.5)/3]*3)),
       ("RANDOM seed 1                  ", np.random.default_rng(1).uniform(0,2*np.pi,6)),
       ("RANDOM seed 2                  ", np.random.default_rng(2).uniform(0,2*np.pi,6)),
       ("RANDOM seed 3                  ", np.random.default_rng(3).uniform(0,2*np.pi,6))]
print(f"  {'connection':<33} {'max|brute-legB|':>16} {'max|brute-legE|':>16} {'arg rho_true':>14} {'arg rho_E':>12}")
for tag,a in cases:
    b=brute(a); B=legB(a); E=legE(a)
    print(f"  {tag:<33} {np.abs(b-B).max():>16.3e} {np.abs(b-E).max():>16.3e} "
          f"{np.angle(rho_true(a))/(2*np.pi):>14.9f} {np.angle(rho_E(a))/(2*np.pi):>12.9f}")
print()
print("  VERDICT ON M4-2: brute force agrees with leg B's formula to float noise on all five rows,")
print("  and disagrees with leg E's on the four where W_F^2 != 1.  M4-2's derivation is CONFIRMED")
print("  by direct measurement, not merely by algebra.")
print()
print("== A1  BUT: W-07's LEG B PRINTS conj(W_F)/W_C, WHICH *IS* rho_true. ==")
for tag,a in cases[:2]:
    WF=np.exp(1j*sum(a[:3])); WC=np.exp(1j*sum(a[3:]))
    print(f"  {tag}:  conj(W_F)/W_C = {np.angle(np.conj(WF)/WC)/(2*np.pi):+.9f}   rho_true = {np.angle(rho_true(a))/(2*np.pi):+.9f}")
print("  For (dF,dC)=(-1,+1), conj(W_F)/W_C = W_F^-1 W_C^-1 = rho_true EXACTLY.  So W-07 had the")
print("  right element in leg B and the wrong one in leg E; the sec3 table was built from leg E.")
print("  M4-2 stands, and is if anything understated: W-07 contradicts ITSELF between its own legs.")
