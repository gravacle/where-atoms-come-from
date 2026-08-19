"""PF-5 / O-9: can an ORDINARY LOCAL environment form a record?

Pre-registered (commit above):
  1. a SINGLE local term      lam * Z_e (x) X_0            -> chi = 0 at every lam   [CONTROL]
  2. a SUM of local terms     lam * sum_l Z_l (x) X_{l%3}  -> chi > 0, slope ~ 2d
  3. a weight-d coupling      lam * Zbar (x) sum_j X_j     -> chi > 0, slope ~ 2
Carrier: toric code 2x2 (dim 256) + 3-qubit bath, evolved unitarily from a PRODUCT state."""
import sys, numpy as np
def say(*a): print(*a); sys.stdout.flush()
exec(open('/Users/bgm/MB Work/where-atoms-come-from/LANE_F7_OCCUPANCY/f7_davies.py').read().split('say("="*104); say("0.')[0])
nS=2**L; nq=3; nB=2**nq; beta=2.0
I2b=np.eye(2); Xb=np.array([[0,1],[1,0]],dtype=complex); Zb=np.array([[1,0],[0,-1]],dtype=complex)
def bop(j,P):
    M=np.array([[1]],dtype=complex)
    for k in range(nq): M=np.kron(M,P if k==j else I2b)
    return M
HB=sum(w*bop(j,Zb) for j,w in enumerate([1.0,1.4,0.7]))
E0,V0=np.linalg.eigh(H0); gs=int(np.sum(np.abs(E0-E0[0])<1e-9)); Pg=V0[:,:gs]@V0[:,:gs].conj().T
def vN(r):
    e=np.linalg.eigvalsh(r); e=e[e>1e-13]; return float(-(e*np.log2(e)).sum())
def chi_at(HINT,lam,t=4.0):
    Ht=np.kron(H0,np.eye(nB))+np.kron(np.eye(nS),HB)+lam*HINT
    w,U=np.linalg.eigh(Ht)
    wB,VB=np.linalg.eigh(HB); pB=np.exp(-beta*wB); pB/=pB.sum()
    r0=np.kron(Pg/gs,(VB*pB)@VB.conj().T)
    Uc=U.conj().T@r0@U; ph=np.exp(-1j*w*t)
    r=U@(ph[:,None]*Uc*ph.conj()[None,:])@U.conj().T
    out=[]
    for s in (+1,-1):
        P=np.kron((np.eye(nS)+s*Zbar)/2,np.eye(nB)); blk=P@r@P; p=np.real(np.trace(blk))
        if p<1e-12: continue
        out.append((p,(blk/p).reshape(nS,nB,nS,nB).trace(axis1=0,axis2=2)))
    if len(out)<2: return 0.0
    av=sum(p*rb for p,rb in out)
    return max(vN(av)-sum(p*vN(rb) for p,rb in out),0.0)

SINGLE = np.kron(op({ind[('h',0,0)]:Z},L), bop(0,Xb))
SUMLOC = sum(np.kron(op({l:Z},L), bop(l%nq,Xb)) for l in range(L))
WEIGHTD= np.kron(Zbar, sum(bop(j,Xb) for j in range(nq)))
say("="*100); say("PF-5  CAN AN ORDINARY LOCAL ENVIRONMENT FORM A RECORD?"); say("="*100)
say(f"  carrier dim {nS}, ground degeneracy {gs}, d = 2; bath {nq} qubits; t = 4.0")
lams=[0.02,0.05,0.1,0.2,0.4]
say(f"\n  {'lambda':>8}{'1. SINGLE local term':>24}{'2. SUM of local terms':>25}{'3. weight-d coupling':>24}")
rows=[]
for lam in lams:
    a,b,c = chi_at(SINGLE,lam), chi_at(SUMLOC,lam), chi_at(WEIGHTD,lam)
    rows.append((lam,a,b,c))
    say(f"  {lam:>8.3f}{a:>24.3e}{b:>25.8f}{c:>24.8f}")
say("")
import numpy as _np
lg=_np.log([r[0] for r in rows])
for i,(nm,pred) in enumerate([("1. SINGLE local term","exactly 0"),("2. SUM of local terms","~2d = 4"),("3. weight-d coupling","~2")],start=1):
    v=_np.array([r[i] for r in rows])
    if v.max()<1e-10: say(f"  {nm:<24} max = {v.max():.3e}  -> EXACTLY ZERO at every lambda   (predicted {pred})")
    else:
        sl=_np.polyfit(lg[:3],_np.log(v[:3]),1)[0]
        say(f"  {nm:<24} log-log slope in lambda = {sl:.4f}   (predicted {pred})")
say("")
say(f"  NOISE FLOOR (lam = 0): {chi_at(SUMLOC,0.0):.3e}")
