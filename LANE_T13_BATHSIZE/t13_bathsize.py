"""T-13 / O-15: chi(t) recurred on a 3-qubit bath (0.975 at t=1, 0.787 at t=2, 0.908 at t=4).
A macroscopic environment would not recur -- but that was ASSERTED, never measured.

The dimension constraint is real: the toric 2x2 is dim 256, so a 7-qubit bath would need
eigh on 32768x32768. Shrink the SYSTEM instead. [[5,1,3]] is dim 32 with verified logicals
Zbar = ZZZZZ and Xbar = XXXXX, so a 7-qubit bath gives 32 x 128 = 4096 -- one eigendecomposition,
then every time point is a phase multiplication.

PREDICTION, before the run: chi(t) recurs on small baths and approaches a MONOTONE PLATEAU as
the bath grows; fragment redundancy evens out."""
import sys, numpy as np
sys.path.insert(0,'/Users/bgm/MB Work/where-atoms-come-from/model')
from record_model import RecordModel, Environment
def say(*a): print(*a); sys.stdout.flush()
I2=np.eye(2); X=np.array([[0,1],[1,0]],dtype=complex); Z=np.array([[1,0],[0,-1]],dtype=complex)
def pauli(s):
    M=np.array([[1]],dtype=complex)
    for c in s: M=np.kron(M,{'I':I2,'X':X,'Z':Z,'Y':1j*(X@Z)}[c])
    return M
S=[pauli(g) for g in ['XZZXI','IXZZX','XIXZZ','ZXIXZ']]
H0=-sum(S); Zbar=pauli('ZZZZZ'); Xbar=pauli('XXXXX'); nS=32
m=RecordModel(H0,[])
say("="*100); say("T-13 / O-15   DOES chi(t) STOP RECURRING AS THE BATH GROWS?"); say("="*100)
say(f"  system [[5,1,3]] dim {nS}   ||[Zbar,H0]|| = {np.linalg.norm(Zbar@H0-H0@Zbar):.1e}"
    f"   ||{{Zbar,Xbar}}|| = {np.linalg.norm(Zbar@Xbar+Xbar@Zbar):.1e}")
ch=m.channel(Zbar,Zbar)
say(f"  channel(Zbar, Zbar): component {ch['component']:.3f}, opens = {ch['opens_channel']}")

def chi_series(nq, ts, lam=0.8, beta=2.0, seed=1):
    rng=np.random.default_rng(seed)
    en=tuple(0.6+1.2*rng.random(nq))
    env=Environment(nq=nq, energies=en, beta=beta)
    Pg,k=m.ground_space(); rho0=Pg/k
    nB=env.dim
    Ht=np.kron(H0,np.eye(nB))+np.kron(np.eye(nS),env.HB)+lam*np.kron(Zbar,env.probe)
    w,U=np.linalg.eigh(Ht)
    r0=np.kron(rho0,env.thermal()); Uc=U.conj().T@r0@U
    out=[]
    for t in ts:
        ph=np.exp(-1j*w*t)
        r=U@(ph[:,None]*Uc*ph.conj()[None,:])@U.conj().T
        out.append(env.holevo(r,Zbar,nS))
    return env,w,U,Uc,out

ts=[0.0,0.5,1.0,1.5,2.0,3.0,4.0,6.0,8.0,12.0,16.0,24.0]
say("")
say(f"  {'t':>6}" + "".join(f"{'nq='+str(q):>12}" for q in (3,5,7)))
series={}
for q in (3,5,7):
    _,_,_,_,vals = chi_series(q, ts)
    series[q]=vals
for i,t in enumerate(ts):
    say(f"  {t:>6.1f}" + "".join(f"{series[q][i]:>12.6f}" for q in (3,5,7)))
say("")
say(f"  {'bath':<10}{'max chi':>10}{'final chi':>11}{'# decreases':>13}{'largest dip':>13}")
for q in (3,5,7):
    v=np.array(series[q]); d=np.diff(v)
    say(f"  nq = {q:<5}{v.max():>10.6f}{v[-1]:>11.6f}{int((d<-1e-6).sum()):>13}{(-d.min() if (d<0).any() else 0):>13.6f}")
say("")
say("  FRAGMENT REDUNDANCY at t = 8.0")
for q in (3,7):
    rng=np.random.default_rng(1); en=tuple(0.6+1.2*rng.random(q))
    env=Environment(nq=q,energies=en,beta=2.0)
    fr=[m.formation(Zbar,Zbar,env,lam=0.8,t=8.0,fragment=[j]) for j in range(q)]
    fr=np.array(fr)
    say(f"    nq = {q}:  per-fragment chi = [{', '.join(f'{x:.4f}' for x in fr)}]")
    say(f"             mean {fr.mean():.4f}   spread {fr.max()-fr.min():.4f}   "
        f"relative spread {(fr.max()-fr.min())/max(fr.mean(),1e-12):.3f}")
