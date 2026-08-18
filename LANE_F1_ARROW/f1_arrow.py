"""F-1 / T-III.6: is there an arrow, and is it OURS or imported decoherence?

PRE-REGISTERED (commit cf2cf8f):
  1. chi(Zbar:B) > 0 for a weight-d coupling, matching a closed form from Z_B(+-1)
  2. chi(Zbar:B) = 0 for a weight-1 coupling
  3. local unitaries leave chi EXACTLY invariant
  4. DECIDING: sweeping observables by weight, chi turns on at weight = d and is zero below.
     If low-weight observables also carry chi, the arrow is AMBIENT DECOHERENCE and NOT ours.
Carrier: toric code 2x2 (dim 256) + 4-level bath, mean-force Gibbs state."""
import sys, itertools, numpy as np
def say(*a): print(*a); sys.stdout.flush()
exec(open('/Users/bgm/MB Work/where-atoms-come-from/LANE_F7_OCCUPANCY/f7_davies.py').read().split('say("="*104); say("0.')[0])
nS=2**L; E0,V0=np.linalg.eigh(H0); gs=int(np.sum(np.abs(E0-E0[0])<1e-9)); Pg=V0[:,:gs]@V0[:,:gs].conj().T
nB=4; EB=np.array([0.0,0.7,1.3,2.1]); bb=np.array([1.0,0.3,-0.2,-0.9]); HB=np.diag(EB); Bop=np.diag(bb)
beta=2.0

def rho_SB(A,lam):
    Ht=np.kron(H0,np.eye(nB))+np.kron(np.eye(nS),HB)+lam*np.kron(A,Bop)
    w,U=np.linalg.eigh(Ht); w=w-w.min(); M=(U*np.exp(-beta*w))@U.conj().T
    return M/np.trace(M)

def vN(r):
    e=np.linalg.eigvalsh(r); e=e[e>1e-13]; return float(-(e*np.log2(e)).sum())

def chi(r,O):
    """Holevo info the bath holds about the +-1 observable O (O^2 = I)."""
    I2=np.eye(nS); out=[]
    for s in (+1,-1):
        P=(I2+s*O)/2
        blk=np.kron(P,np.eye(nB))@r@np.kron(P,np.eye(nB))
        p=np.real(np.trace(blk))
        if p<1e-12: out.append((0.0,np.zeros((nB,nB)))); continue
        rB=(blk/p).reshape(nS,nB,nS,nB).trace(axis1=0,axis2=2)
        out.append((p,rB))
    rBav=sum(p*rB for p,rB in out)
    return max(vN(rBav)-sum(p*vN(rB) for p,rB in out),0.0)

say("="*100); say("F-1  THE ARROW: DOES THE BATH HOLD A COPY, AND AT WHAT THRESHOLD?"); say("="*100)
say(f"  SELF-CHECK ||[Zbar,H_S]|| = {np.linalg.norm(Zbar@H0-H0@Zbar):.1e}   ground degeneracy {gs}   d = 2")
lam=0.8
say("")
say("1-2.  HOLEVO INFORMATION THE BATH HOLDS ABOUT THE RECORD")
say(f"  {'coupling A':<22}{'weight':>8}{'||[A,Xbar]||':>14}{'chi(Zbar:B) bits':>19}{'closed form':>14}")
Ze=op({ind[('h',0,0)]:Z},L)
for nm,A,wt in (("Zbar  (logical)",Zbar,2),("Zbar2 (logical)",Zbar2,2),("Z_e   (single site)",Ze,1),("identity (no coupling)",np.eye(nS),0)):
    r=rho_SB(A,lam); c=chi(r,Zbar)
    cf=""
    if nm.startswith("Zbar "):
        pB=np.exp(-beta*(EB+lam*bb)); pM=np.exp(-beta*(EB-lam*bb))
        Zp,Zm=pB.sum(),pM.sum(); pp,pm=Zp/(Zp+Zm),Zm/(Zp+Zm)
        av=pp*(pB/Zp)+pm*(pM/Zm)
        H=lambda v:float(-(v[v>1e-13]*np.log2(v[v>1e-13])).sum())
        cf=f"{H(av)-pp*H(pB/Zp)-pm*H(pM/Zm):.8f}"
    say(f"  {nm:<22}{wt:>8}{np.linalg.norm(A@Xbar-Xbar@A):>14.3f}{c:>19.8f}{cf:>14}")

say("")
say("3.  DO LOCAL (SYSTEM-ONLY) UNITARIES CHANGE chi?   [if not, the system cannot unmake the copy]")
r=rho_SB(Zbar,lam); c0=chi(r,Zbar); rng=np.random.default_rng(5); worst=0.0
for _ in range(12):
    Mx=rng.normal(size=(nS,nS))+1j*rng.normal(size=(nS,nS))
    Q,_=np.linalg.qr(Mx); Us=np.kron(Q,np.eye(nB))
    rp=Us@r@Us.conj().T
    worst=max(worst,abs(chi(rp,Q.conj().T@Zbar@Q)-c0))
say(f"     chi before = {c0:.8f}   max |change| over 12 random system-only unitaries = {worst:.3e}")
say(f"     -> {'INVARIANT: no system-only operation can erase the bath copy' if worst<1e-8 else 'CHANGES -- prediction 3 falsified'}")

say("")
say("4.  THE DECIDING TEST -- does chi carry the RECORD'S threshold, or is it ambient decoherence?")
say("     sweeping ALL observables by weight, in the SAME state (coupling A = Zbar, lambda = 0.8)")
P1={'X':X,'Y':1j*(X@Z),'Z':Z}
say(f"  {'observable weight':>19}{'# swept':>10}{'max chi(O:B) bits':>21}{'argmax':>14}")
r=rho_SB(Zbar,lam)
for w in (1,2):
    mx=0.0; arg=None; n=0
    for sites in itertools.combinations(range(L),w):
        for lets in itertools.product('XYZ',repeat=w):
            O=op({s:P1[c] for s,c in zip(sites,lets)},L); n+=1
            v=chi(r,O)
            if v>mx: mx=v; arg=''.join(lets)+str(list(sites))
    say(f"  {w:>19}{n:>10}{mx:>21.8f}{str(arg):>14}")
say(f"  {'logical (weight 2)':>19}{2:>10}{max(chi(r,Zbar),chi(r,Zbar2)):>21.8f}{'Zbar/Zbar2':>14}")
