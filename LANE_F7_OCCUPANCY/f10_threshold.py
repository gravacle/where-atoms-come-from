"""F-10: what is the MINIMUM WEIGHT of a bath coupling that can FORM a record?

PRE-REGISTERED PREDICTION (REGISTER_V001, before this run):
  by Knill-Laflamme, wt(C) < d implies P_g C P_g ~ P_g, so C cannot distinguish record
  states.  Hence NO coupling of weight < d forms a record, and the FORMATION threshold
  equals the DESTRUCTION threshold: both are d.  At d=2: weight 1 gives exactly 0.

Non-circular by construction: we do not hand the bath a logical operator and observe
success; we sweep ALL Pauli couplings by weight and locate the threshold.
LIMITATION stated in advance: only d=2 is reachable (dim 256). Scaling is NOT tested."""
import sys, itertools, numpy as np
def say(*a): print(*a); sys.stdout.flush()
exec(open('f7_davies.py').read().split('say("="*104); say("0.')[0])   # carrier + computed logicals

E,V=np.linalg.eigh(H0); gs=int(np.sum(np.abs(E-E[0])<1e-9))
dE=E[None,:]-E[:,None]; ws=np.unique(np.round(dE,6))
Pge=np.zeros((256,256),dtype=complex); Pge[:gs,:gs]=np.eye(gs)      # code space in eigenbasis
Zb=V.conj().T@Zbar@V; Zb2=V.conj().T@Zbar2@V; Xb=V.conj().T@Xbar@V
rho0=Pge/gs
beta=2.0
def kms(w,b=beta,g=1.0): return g*np.exp(b*w/2)/(2*np.cosh(b*w/2))

def rate(C):
    """d<Zbar>/dt at the maximally mixed code state, for a bath coupling C."""
    Cp=V.conj().T@C@V; d=np.zeros_like(rho0)
    for w in ws:
        g=kms(w)
        if g<1e-14: continue
        M=np.where(np.abs(dE-w)<1e-8,Cp,0.0)
        if np.abs(M).max()<1e-12: continue
        J=np.sqrt(g)*M; Jd=J.conj().T
        d+=J@rho0@Jd-0.5*(Jd@J@rho0+rho0@Jd@J)
    return (np.real(np.trace(d@Zb)), np.real(np.trace(d@Zb2)), np.real(np.trace(d@Xb)))

P1={'X':X,'Y':1j*(X@Z),'Z':Z}
say("="*100); say("F-10  MINIMUM COUPLING WEIGHT THAT FORMS A RECORD    (H = H0, exactly degenerate)")
say(f"      carrier: 2x2 torus, dim {2**L}, ground degeneracy {gs}, code distance d = 2")
say(f"      logical weights: Zbar={LOGINFO[0]}, Zbar2={LOGINFO[1]}, Xbar={LOGINFO[2]}"); say("="*100)

say(""); say("  SELF-CHECK: does the machinery register a nonzero rate for ANY coupling?")
say("  (a zero everywhere would be a dead instrument, not a result)")

best={}
for w in (1,2):
    mx=0.0; arg=None; nz=0; tot=0
    for sites in itertools.combinations(range(L),w):
        for letters in itertools.product('XYZ',repeat=w):
            C=op({s:P1[c] for s,c in zip(sites,letters)},L)
            r=rate(C); m=max(abs(x) for x in r); tot+=1
            if m>1e-10: nz+=1
            if m>mx: mx=m; arg=(sites,''.join(letters),r)
    best[w]=(mx,arg,nz,tot)
    say("")
    say(f"  WEIGHT {w}:  {tot} couplings swept   nonzero rate in {nz} of them")
    say(f"     max |d<logical>/dt| = {mx:.3e}"
        + (f"   at {arg[1]} on sites {arg[0]}" if arg else ""))
    if arg: say(f"     (d<Zbar>/dt, d<Zbar2>/dt, d<Xbar>/dt) = "
                f"({arg[2][0]:+.3e}, {arg[2][1]:+.3e}, {arg[2][2]:+.3e})")

say(""); say("-"*100)
say(f"  weight 1 : max rate {best[1][0]:.3e}   -> {'ZERO -- cannot form a record' if best[1][0]<1e-10 else 'NONZERO -- PREDICTION FALSIFIED'}")
say(f"  weight 2 : max rate {best[2][0]:.3e}   -> {'NONZERO -- can form a record' if best[2][0]>1e-10 else 'ZERO'}")
say("")
thr = 1 if best[1][0]>1e-10 else (2 if best[2][0]>1e-10 else None)
say(f"  MEASURED FORMATION THRESHOLD = {thr}      CODE DISTANCE d = 2")
say(f"  PREDICTION (threshold == d): {'CONFIRMED' if thr==2 else 'FALSIFIED'}")

say(""); say("="*100); say("  CROSS-CHECK -- Knill-Laflamme directly: is P_g C P_g proportional to P_g?"); say("="*100)
say(f"  {'weight':>8}{'max ||P C P - (tr/4) P||':>28}{'KL satisfied?':>18}")
for w in (1,2):
    mx=0.0
    for sites in itertools.combinations(range(L),w):
        for letters in itertools.product('XYZ',repeat=w):
            C=V.conj().T@op({s:P1[c] for s,c in zip(sites,letters)},L)@V
            M=Pge@C@Pge; c=np.trace(M)/gs
            mx=max(mx,np.linalg.norm(M-c*Pge))
    say(f"  {w:>8}{mx:>28.3e}{('YES -- acts as a scalar' if mx<1e-10 else 'NO -- distinguishes code states'):>18}")
