"""O-18: why should the world's noise be shaped so a record exists?
Candidate answer: it is EM's own constraint. An environment coupling to a gauge field must
couple GAUGE-INVARIANTLY -- commute with every Gauss operator A_v. Pre-registered above."""
import sys, itertools, numpy as np
def say(*a): print(*a); sys.stdout.flush()
exec(open('/Users/bgm/MB Work/where-atoms-come-from/LANE_F7_OCCUPANCY/f7_davies.py').read().split('say("="*104); say("0.')[0])
N=2**L
AV=[op({l:X for l in s},L) for s in STAR]        # Gauss operators
BP=[op({l:Z for l in p},L) for p in PLAQ]        # plaquettes
E0,V0=np.linalg.eigh(H0); gs=int(np.sum(np.abs(E0-E0[0])<1e-9)); Pg=V0[:,:gs]@V0[:,:gs].conj().T
say("="*98); say("O-18  DOES GAUGE INVARIANCE SUPPLY THE STRUCTURE THE NOISE NEEDS?"); say("="*98)
say(f"  toric 2x2, dim {N}, ground degeneracy {gs}")
say(f"  SELF-CHECK ||[Zbar,H0]|| = {np.linalg.norm(Zbar@H0-H0@Zbar):.1e}  "
    f"||{{Zbar,Xbar}}|| = {np.linalg.norm(Zbar@Xbar+Xbar@Zbar):.1e}")

P1={'X':X,'Y':1j*(X@Z),'Z':Z}
def gauge_inv(A): return max(np.linalg.norm(A@G-G@A) for G in AV) < 1e-9
def trivial_on_code(A):
    M=Pg@A@Pg; c=np.trace(M)/gs
    return np.linalg.norm(M-c*Pg) < 1e-9

say("\n1.  WHICH LOCAL OPERATORS ARE GAUGE-INVARIANT, AND DO THEY ACT ON THE CODE SPACE?")
say(f"  {'weight':>7}{'# Paulis':>10}{'gauge-invariant':>18}{'of those, TRIVIAL on code':>28}")
for w in (1,2,3,4):
    tot=gi=gitriv=0
    for sites in itertools.combinations(range(L),w):
        for lets in itertools.product('XYZ',repeat=w):
            A=op({s:P1[c] for s,c in zip(sites,lets)},L); tot+=1
            if gauge_inv(A):
                gi+=1
                if trivial_on_code(A): gitriv+=1
    say(f"  {w:>7}{tot:>10}{gi:>18}{f'{gitriv} of {gi}':>28}")

say("\n2.  IS PF-5's COUPLING GAUGE-INVARIANT?   (its terms are single-site Z_l)")
Ze=op({ind[('h',0,0)]:Z},L)
say(f"    single Z_l:  max ||[Z_l, A_v]|| = {max(np.linalg.norm(Ze@G-G@Ze) for G in AV):.4f}   "
    f"{'gauge-invariant' if gauge_inv(Ze) else 'NOT gauge-invariant'}")
say(f"    a plaquette: max ||[B_p, A_v]|| = {max(np.linalg.norm(BP[0]@G-G@BP[0]) for G in AV):.4f}   "
    f"{'gauge-invariant' if gauge_inv(BP[0]) else 'NOT gauge-invariant'}")
say(f"    Zbar:        max ||[Zbar,A_v]|| = {max(np.linalg.norm(Zbar@G-G@Zbar) for G in AV):.4f}   "
    f"{'gauge-invariant' if gauge_inv(Zbar) else 'NOT gauge-invariant'}")

say("\n3.  CAN GAUGE-INVARIANT LOCAL NOISE FORM A RECORD?   chi(Zbar : bath), 3-qubit bath, t = 4")
nq=3; nB=2**nq; beta=2.0
I2b=np.eye(2); Xb2=np.array([[0,1],[1,0]],dtype=complex); Zb2=np.array([[1,0],[0,-1]],dtype=complex)
def bop(j,P):
    M=np.array([[1]],dtype=complex)
    for k in range(nq): M=np.kron(M,P if k==j else I2b)
    return M
HB=sum(w*bop(j,Zb2) for j,w in enumerate([1.0,1.4,0.7]))
def vN(r):
    e=np.linalg.eigvalsh(r); e=e[e>1e-13]; return float(-(e*np.log2(e)).sum())
def chi_of(HINT,lam,t=4.0):
    Ht=np.kron(H0,np.eye(nB))+np.kron(np.eye(N),HB)+lam*HINT
    w,U=np.linalg.eigh(Ht); wB,VB=np.linalg.eigh(HB); pB=np.exp(-beta*wB); pB/=pB.sum()
    r0=np.kron(Pg/gs,(VB*pB)@VB.conj().T); Uc=U.conj().T@r0@U; ph=np.exp(-1j*w*t)
    r=U@(ph[:,None]*Uc*ph.conj()[None,:])@U.conj().T
    out=[]
    for s in (+1,-1):
        P=np.kron((np.eye(N)+s*Zbar)/2,np.eye(nB)); blk=P@r@P; p=np.real(np.trace(blk))
        if p<1e-12: continue
        out.append((p,(blk/p).reshape(N,nB,N,nB).trace(axis1=0,axis2=2)))
    if len(out)<2: return 0.0
    av=sum(p*rb for p,rb in out)
    return max(vN(av)-sum(p*vN(rb) for p,rb in out),0.0)
COUPLINGS=[("gauge-invariant local (plaquettes)", sum(np.kron(BP[i],bop(i%nq,Xb2)) for i in range(len(BP)))),
           ("NOT gauge-invariant (PF-5's sum Z_l)", sum(np.kron(op({l:Z},L),bop(l%nq,Xb2)) for l in range(L))),
           ("gauge-invariant NON-local (Zbar)",     np.kron(Zbar, sum(bop(j,Xb2) for j in range(nq))))]
say(f"  {'coupling':<40}{'gauge-inv?':>12}{'local?':>9}{'chi bits at lam=0.8':>22}")
for nm,HINT in COUPLINGS:
    gi = 'yes' if 'NOT' not in nm else 'no'
    loc= 'no' if 'NON-local' in nm else 'yes'
    say(f"  {nm:<40}{gi:>12}{loc:>9}{chi_of(HINT,0.8):>22.8f}")
say(f"\n  NOISE FLOOR (lam = 0): {chi_of(COUPLINGS[0][1],0.0):.3e}")
